import numpy as np

def particle_filter_analysis(data, data_cov, param, ensemble, observation, **kwargs):
    """
    Particle Filter analysis update with optional hybrid EnKF update
    
    Args:
        data: (meas_size, 1) or (meas_size,) - actual observations
        data_cov: (meas_size, meas_size) - observation error covariance
        param: (par_size, ens_size) - parameter ensemble (MATCHES YOUR ORIGINAL FORMAT)
        ensemble: (sim_size, ens_size) - state ensemble (e.g., saturation)
        observation: (ens_size, meas_size) or (meas_size, ens_size) - predicted observations
        **kwargs:
            use_enkf_update: bool - use weighted EnKF update after resampling
            resample_threshold: float - resample when N_eff < threshold * ens_size
            jitter_std_param: float - jittering for parameters after resampling
            jitter_std_state: float - jittering for states after resampling
            alpha: float - covariance inflation (for hybrid mode)
    
    Returns:
        dict with:
            - Analysis: (sim_size, ens_size) - updated state ensemble
            - Analysisparam: (ens_size, par_size) - updated parameter ensemble
            - weights: (ens_size,) - particle weights
            - n_eff: float - effective sample size
            - resampled: bool - whether resampling occurred
    """
    
    # Default parameters
    use_enkf_update = kwargs.get('use_enkf_update', False)
    resample_threshold = kwargs.get('resample_threshold', 0.5)
    jitter_std_param = kwargs.get('jitter_std_param', 0.01)
    jitter_std_state = kwargs.get('jitter_std_state', 0.005)
    use_log_K = kwargs.get('use_log_K',False)
    param_bounds = kwargs.get('param_bounds',False)
    alpha = kwargs.get('alpha', 1.0)
    
    # Collect data sizes
    ens_size = ensemble.shape[1]
    sim_size = ensemble.shape[0]
    par_size = param.shape[0]  # param is (par_size, ens_size)
    
    # Ensure observation has correct shape: (meas_size, ens_size)
    if observation.shape[0] == ens_size:
        observation = observation.T
    
    # Ensure data has shape (meas_size,) or (meas_size, 1)
    data = np.atleast_1d(data.flatten())
    meas_size = data.shape[0]
    
    # Extract observation error standard deviations
    if data_cov.ndim == 2:
        obs_std = np.sqrt(np.diag(data_cov))
    else:
        obs_std = np.sqrt(data_cov)
    
    print(f"Particle Filter Analysis: ens_size={ens_size}, meas_size={meas_size}, par_size={par_size}")
    print(f"  param shape: {param.shape}, ensemble shape: {ensemble.shape}")
    
    # =========================================================================
    # STEP 1: Compute particle weights using exact likelihood
    # =========================================================================
    
    log_weights = np.zeros(ens_size)
    
    for i in range(ens_size):
        # Innovation: difference between actual and predicted observation
        innovation = data - observation[:, i]
        
        # Log-likelihood (assuming Gaussian observation errors)
        # log p(y|x) = -0.5 * sum((y - h(x))^2 / sigma^2)
        log_likelihood = -0.5 * np.sum((innovation / obs_std) ** 2)
        log_weights[i] = log_likelihood
    
    # Normalize weights (subtract max for numerical stability)
    log_weights -= log_weights.max()
    weights = np.exp(log_weights)
    weights /= weights.sum()
    
    # Compute effective sample size
    n_eff = 1.0 / np.sum(weights ** 2)
    print(f"  Effective sample size: N_eff = {n_eff:.1f} / {ens_size}")
    print(f"  Weight range: [{weights.min():.3e}, {weights.max():.3e}]")
    print(f"  Max/Min weight ratio: {weights.max()/weights.min():.3e}")
    
    # Diagnose issues
    diag = diagnose_particle_filter_issues(param, weights, observation, data)
    for warning in diag['warnings']:
        print(f"  ⚠️  {warning}")
    
    # =========================================================================
    # STEP 2: Resample if needed (when particle degeneracy occurs)
    # =========================================================================
    
    resampled = False
    if n_eff < resample_threshold * ens_size:
        print(f"  Resampling (N_eff < {resample_threshold * ens_size:.1f})...")
        resampled = True
        
        # Systematic resampling
        indices = systematic_resample(weights, ens_size)
        
        # Resample ensemble and parameters (both are indexed by columns)
        ensemble = ensemble[:, indices]  # (sim_size, ens_size)
        param = param[:, indices]  # (par_size, ens_size)
        observation = observation[:, indices]  # (meas_size, ens_size)
        
        # Reset weights to uniform
        weights = np.ones(ens_size) / ens_size
    
    # =========================================================================
    # STEP 3: Optional hybrid weighted EnKF update
    # =========================================================================
    
    if use_enkf_update:
        print("  Applying weighted EnKF update...")
        ensemble, param = weighted_enkf_update(
            data, data_cov, param, ensemble, observation, weights, alpha
        )
    
    # =========================================================================
    # STEP 4: Add jitter to prevent particle collapse
    # =========================================================================
    
    if resampled or jitter_std_param > 0:  # Can apply jitter even without resampling
        print("  Adding jitter to maintain diversity...")
        
        # Jitter for states (ensemble)
        if jitter_std_state > 0:
            ensemble += np.random.randn(*ensemble.shape) * jitter_std_state
        
        # Jitter for parameters (param is par_size x ens_size)
        if jitter_std_param > 0:
            if use_log_K:
                # Apply multiplicative jitter in log-space (better for K)
                # This keeps K positive and maintains scale
                # Assumes param has [phi_cells, K_cells] stacked
                n_cells = par_size // 2
                
                # Additive jitter for porosity (first half)
                param[:n_cells, :] += np.random.randn(n_cells, ens_size) * jitter_std_param
                
                # Multiplicative jitter for K (second half)
                param[n_cells:, :] *= np.exp(np.random.randn(n_cells, ens_size) * jitter_std_param)
            else:
                # Additive jitter for all parameters
                param += np.random.randn(*param.shape) * jitter_std_param
        
        print(f"  param range AFTER jitter: [{param.min():.3e}, {param.max():.3e}]")
    
    # =========================================================================
    # STEP 5: Enforce parameter bounds
    # =========================================================================
    
    # if param_bounds is not None:
    #     print("  Enforcing parameter bounds...")
        
    #     # Assuming param has [phi_cells, K_cells] stacked
    #     # You need to specify how your param is organized
    #     if 'phi' in param_bounds and 'K' in param_bounds:
    #         n_cells = par_size // 2
    #         phi_min, phi_max = param_bounds['phi']
    #         K_min, K_max = param_bounds['K']
            
    #         # Clip porosity (first half)
    #         param[:n_cells, :] = np.clip(param[:n_cells, :], phi_min, phi_max)
            
    #         # Clip permeability (second half)
    #         param[n_cells:, :] = np.clip(param[n_cells:, :], K_min, K_max)
            
    #         print(f"  param range AFTER bounds: [{param.min():.3e}, {param.max():.3e}]")
    #     elif 'all' in param_bounds:
    #         # Simple bounds for all parameters
    #         param_min, param_max = param_bounds['all']
    #         param = np.clip(param, param_min, param_max)
    #         print(f"  param range AFTER bounds: [{param.min():.3e}, {param.max():.3e}]")
    
    # =========================================================================
    # Return results
    # =========================================================================
    
    Analysis = ensemble
    Analysisparam = param
    
    result = {
        'Analysis': Analysis,
        'Analysisparam': Analysisparam,
        'weights': weights,
        'n_eff': n_eff,
        'resampled': resampled,
        'observation': observation  # Updated observations after resampling
    }
    
    return result


def systematic_resample(weights, n_particles):
    """
    Systematic resampling algorithm
    
    Args:
        weights: (n_particles,) normalized weights
        n_particles: number of particles to resample
    
    Returns:
        indices: (n_particles,) resampled indices
    """
    positions = (np.arange(n_particles) + np.random.rand()) / n_particles
    cumsum = np.cumsum(weights)
    indices = np.searchsorted(cumsum, positions)
    return indices


def diagnose_particle_filter_issues(param, weights, observation, data):
    """
    Diagnose common issues causing parameter explosion in PF
    
    Args:
        param: (par_size, ens_size) parameter ensemble
        weights: (ens_size,) particle weights
        observation: (meas_size, ens_size) predicted observations
        data: (meas_size,) actual observations
    
    Returns:
        dict of diagnostic information
    """
    diagnostics = {}
    
    # 1. Weight concentration
    diagnostics['weight_concentration'] = weights.max() / weights.min()
    diagnostics['n_eff'] = 1.0 / np.sum(weights**2)
    
    # 2. Parameter variance
    diagnostics['param_std'] = param.std(axis=1).mean()
    diagnostics['param_range'] = param.max() - param.min()
    
    # 3. Innovation magnitude
    innovation = data[:, np.newaxis] - observation
    diagnostics['innovation_norm'] = np.linalg.norm(innovation, axis=0).mean()
    
    # 4. Check for outliers in parameters
    param_median = np.median(param, axis=1, keepdims=True)
    param_mad = np.median(np.abs(param - param_median), axis=1, keepdims=True)
    outlier_threshold = 5  # MAD units
    n_outliers = np.sum(np.abs(param - param_median) > outlier_threshold * param_mad)
    diagnostics['n_outliers'] = n_outliers
    
    # 5. Warnings
    warnings = []
    if diagnostics['weight_concentration'] > 1e6:
        warnings.append("CRITICAL: Extreme weight concentration - almost all weight on one particle")
    if diagnostics['n_eff'] < 10:
        warnings.append("WARNING: Very low N_eff - particle degeneracy")
    if diagnostics['param_range'] > 100 * diagnostics['param_std']:
        warnings.append("WARNING: Extreme parameter spread - likely has outliers")
    if n_outliers > param.shape[1] * 0.1:
        warnings.append(f"WARNING: {n_outliers} parameter outliers detected")
    
    diagnostics['warnings'] = warnings
    
    return diagnostics


def weighted_enkf_update(data, data_cov, param, ensemble, observation, 
                         weights, alpha=1.0, gain_factor=1.0):
    """
    Weighted EnKF update using particle weights
    
    This combines particle filter weighting with EnKF spatial correlation structure
    
    Args:
        data: (meas_size,) observations
        data_cov: (meas_size, meas_size) observation covariance
        param: (par_size, ens_size) parameter ensemble
        ensemble: (sim_size, ens_size) state ensemble
        observation: (meas_size, ens_size) predicted observations
        weights: (ens_size,) particle weights
        alpha: covariance inflation factor
        gain_factor: reduce Kalman gain (0-1) to prevent overshooting
    
    Returns:
        ensemble_updated: (sim_size, ens_size)
        param_updated: (par_size, ens_size)
    """
    
    ens_size = ensemble.shape[1]
    sim_size = ensemble.shape[0]
    meas_size = data.shape[0]
    par_size = param.shape[0]
    
    print(f"    Weighted EnKF: gain_factor={gain_factor}, alpha={alpha}")
    
    # Combine ensemble and parameters (both are already column-wise)
    A = np.vstack([ensemble, param])  # (sim_size + par_size, ens_size)
    
    # Compute weighted mean
    Amean = np.sum(weights[np.newaxis, :] * A, axis=1, keepdims=True)
    
    # Compute weighted observation mean
    MeasAvg = np.sum(weights[np.newaxis, :] * observation, axis=1, keepdims=True)
    
    # Inflate ensemble (ONLY for computing covariances, not the actual values)
    A_inflated = np.sqrt(alpha) * (A - Amean) + Amean
    observation_inflated = np.sqrt(alpha) * (observation - MeasAvg) + MeasAvg
    
    # Compute weighted anomalies
    dA = A_inflated - Amean  # (sim_size + par_size, ens_size)
    
    # Data perturbation
    # Add noise to observations for each ensemble member
    data_pert = data[:, np.newaxis] + np.random.randn(meas_size, ens_size) * np.sqrt(np.diag(data_cov))[:, np.newaxis]
    dD = data_pert - observation_inflated  # (meas_size, ens_size)
    
    # Observation anomalies
    S = observation_inflated - MeasAvg  # (meas_size, ens_size)
    
    # Weight the anomalies
    W_sqrt = np.sqrt(weights)[np.newaxis, :]  # (1, ens_size)
    S_weighted = S * W_sqrt  # (meas_size, ens_size)
    dA_weighted = dA * W_sqrt  # (sim_size + par_size, ens_size)
    
    # Compute weighted covariance
    COV = S_weighted @ S_weighted.T + data_cov
    
    # Add regularization to avoid ill-conditioning
    COV += 1e-6 * np.eye(meas_size)
    
    # Solve for Kalman gain times innovation
    B = np.linalg.solve(COV, dD)  # (meas_size, ens_size)
    
    # Cross-covariance
    dAS = dA_weighted @ S_weighted.T  # (sim_size + par_size, meas_size)
    
    # Analysis update with reduced gain
    Analysis = A + gain_factor * dAS @ B  # (sim_size + par_size, ens_size)
    
    # Extract state and parameters
    ensemble_updated = Analysis[:sim_size, :]  # (sim_size, ens_size)
    param_updated = Analysis[sim_size:, :]  # (par_size, ens_size)
    
    # Check for extreme updates
    param_change = np.abs(param_updated - param)
    max_change = param_change.max()
    mean_change = param_change.mean()
    print(f"    Max param change: {max_change:.3e}, Mean change: {mean_change:.3e}")
    
    if max_change > 10 * param.std():
        print(f"    WARNING: Large parameter update detected! Consider reducing gain_factor or alpha")
    
    return ensemble_updated, param_updated


# ============================================================================
# SIMPLE WRAPPER - Easy to use version
# ============================================================================

def particle_filter_simple(data, data_cov, param, ensemble, observation, 
                          ens_size=None, reduce_gain=True):
    """
    Simplified particle filter interface - just works without tuning
    
    Args:
        data: observations
        data_cov: observation error covariance
        param: (par_size, ens_size) - parameter ensemble
        ensemble: (sim_size, ens_size) - state ensemble
        observation: predicted observations
        ens_size: ensemble size (auto-detected if None)
        reduce_gain: if True, use conservative settings to avoid K explosion
    
    Returns:
        Analysis, Analysisparam (same as your original EnKF function)
    """
    
    if ens_size is None:
        ens_size = ensemble.shape[1]
    
    par_size = param.shape[0]
    n_cells = par_size // 2  # Assuming [phi, K] stacked
    
    if reduce_gain:
        # Conservative settings to prevent K explosion
        kwargs = {
            'use_enkf_update': True,
            'enkf_gain_factor': 0.3,  # Reduced gain
            'use_log_K': True,
            'jitter_std_param': 0.02,
            'jitter_std_state': 0.005,
            'resample_threshold': 0.5,
            'alpha': 1.05,
            'param_bounds': {
                'phi': (0.01, 0.99),
                'K': (1e-6, 1e6)
            }
        }
    else:
        # Standard settings
        kwargs = {
            'use_enkf_update': False,
            'use_log_K': True,
            'jitter_std_param': 0.05,
            'jitter_std_state': 0.005,
            'resample_threshold': 0.5,
            'param_bounds': {
                'phi': (0.01, 0.99),
                'K': (1e-6, 1e6)
            }
        }
    
    result = particle_filter_analysis(
        data, data_cov, param, ensemble, observation, **kwargs
    )
    
    return result['Analysis'], result['Analysisparam']


# ============================================================================
# Example usage matching your original interface
# ============================================================================

# def example_usage():
#     """
#     Example demonstrating how to use the particle filter with your data structure
#     AND how to avoid K explosion issues
#     """
    
#     # Setup
#     ens_size = 100
#     sim_size = 50  # Number of grid cells for saturation
#     meas_size = 10  # Number of ERT observations
#     par_size = 2 * sim_size  # Porosity and permeability for each cell
    
#     # Create synthetic data (matching your original format)
#     ensemble = np.random.rand(sim_size, ens_size) * 0.5 + 0.3  # Saturation: (sim_size, ens_size)
    
#     # IMPORTANT: Initialize K in log-space for better stability
#     log_K_mean = 0.0
#     log_K_std = 1.0
#     param = np.zeros((par_size, ens_size))
#     param[:sim_size, :] = np.random.rand(sim_size, ens_size) * 0.3 + 0.1  # Porosity [0.1, 0.4]
#     param[sim_size:, :] = np.exp(np.random.randn(sim_size, ens_size) * log_K_std + log_K_mean)  # K log-normal
    
#     # Predicted observations (from Archie's law)
#     observation = np.random.rand(ens_size, meas_size) * 50 + 10  # Resistivity
    
#     # Actual observations
#     data = np.random.rand(meas_size) * 50 + 10
    
#     # Observation error covariance
#     data_cov = np.diag(np.ones(meas_size) * 2.0)  # 2 Ohm-m standard deviation
    
#     # ========================================================================
#     # TEST 1: Pure PF (most likely to have K explosion)
#     # ========================================================================
#     print("\n" + "="*70)
#     print("TEST 1: PURE PARTICLE FILTER (Watch for K explosion)")
#     print("="*70)
    
#     param_test1 = param.copy()
#     ensemble_test1 = ensemble.copy()
    
#     result_pf = particle_filter_analysis(
#         data, data_cov, param_test1, ensemble_test1, observation,
#         use_enkf_update=False,
#         resample_threshold=0.5,
#         jitter_std_param=0.05,  # Multiplicative in log space
#         jitter_std_state=0.005,
#         use_log_K=True,  # KEY: Use multiplicative jitter for K
#         param_bounds={'phi': (0.01, 0.99), 'K': (1e-6, 1e6)}
#     )
    
#     print(f"✓ Final K range: [{result_pf['Analysisparam'][sim_size:, :].min():.3e}, "
#           f"{result_pf['Analysisparam'][sim_size:, :].max():.3e}]")
    
#     # ========================================================================
#     # TEST 2: Hybrid with REDUCED GAIN (recommended for your problem)
#     # ========================================================================
#     print("\n" + "="*70)
#     print("TEST 2: HYBRID PF-EnKF with REDUCED GAIN (Recommended)")
#     print("="*70)
    
#     param_test2 = param.copy()
#     ensemble_test2 = ensemble.copy()
    
#     result_hybrid = particle_filter_analysis(
#         data, data_cov, param_test2, ensemble_test2, observation,
#         use_enkf_update=True,
#         resample_threshold=0.5,
#         jitter_std_param=0.02,  # Smaller jitter when using EnKF
#         jitter_std_state=0.005,
#         use_log_K=True,
#         alpha=1.05,  # Slight inflation
#         enkf_gain_factor=0.3,  # CRITICAL: Reduce gain to prevent overshooting
#         param_bounds={'phi': (0.01, 0.99), 'K': (1e-6, 1e6)}
#     )
    
#     print(f"✓ Final K range: [{result_hybrid['Analysisparam'][sim_size:, :].min():.3e}, "
#           f"{result_hybrid['Analysisparam'][sim_size:, :].max():.3e}]")
    
#     # ========================================================================
#     # TEST 3: Conservative approach for first few assimilations
#     # ========================================================================
#     print("\n" + "="*70)
#     print("TEST 3: CONSERVATIVE (For early assimilation cycles)")
#     print("="*70)
    
#     param_test3 = param.copy()
#     ensemble_test3 = ensemble.copy()
    
#     result_conservative = particle_filter_analysis(
#         data, data_cov, param_test3, ensemble_test3, observation,
#         use_enkf_update=True,
#         resample_threshold=0.3,  # Resample less often
#         jitter_std_param=0.01,  # Very small jitter
#         jitter_std_state=0.002,
#         use_log_K=True,
#         alpha=1.0,  # No inflation
#         enkf_gain_factor=0.1,  # Very conservative gain
#         param_bounds={'phi': (0.05, 0.45), 'K': (1e-3, 1e3)}  # Tighter bounds
#     )
    
#     print(f"✓ Final K range: [{result_conservative['Analysisparam'][sim_size:, :].min():.3e}, "
#           f"{result_conservative['Analysisparam'][sim_size:, :].max():.3e}]")
    
#     # ========================================================================
#     # Print recommendations
#     # ========================================================================
#     print("\n" + "="*70)
#     print("RECOMMENDATIONS TO AVOID K EXPLOSION:")
#     print("="*70)
#     print("1. Use use_log_K=True (multiplicative jitter in log-space)")
#     print("2. Set enkf_gain_factor=0.1-0.5 (reduce Kalman gain)")
#     print("3. Use tight param_bounds, especially for K")
#     print("4. Start with small jitter_std_param (0.01-0.05)")
#     print("5. Increase gain_factor gradually as filter stabilizes")
#     print("6. Monitor diagnostics: if N_eff < 10, you have problems")
#     print("="*70)
    
#     return result_hybrid


# # ============================================================================
# # Wrapper to maintain compatibility with your original function signature
# # ============================================================================

# def enkf_analysis_inflation_pf_version(data, data_cov, param, ensemble, 
#                                        observation, **kwargs):
#     """
#     Drop-in replacement for your enkf_analysis_inflation function
    
#     This provides the same return structure but uses particle filtering
    
#     Args:
#         param: (par_size, ens_size) - matches your original format
#         ensemble: (sim_size, ens_size)
    
#     Returns list matching original: 
#     [A, Amean, dA, dD, MeasAvg, S, COV, B, dAS, Analysis, Analysisparam]
#     """
    
#     # Run particle filter
#     result = particle_filter_analysis(data, data_cov, param, ensemble, 
#                                      observation, **kwargs)
    
#     # For compatibility, we'll create dummy values for EnKF-specific outputs
#     # that don't have direct equivalents in particle filtering
#     ens_size = ensemble.shape[1]
#     sim_size = ensemble.shape[0]
#     par_size = param.shape[0]
    
#     # Combine for output format (both are column-wise)
#     A = np.vstack([ensemble, param])  # (sim_size + par_size, ens_size)
#     Amean = np.mean(A, axis=1, keepdims=True)
#     dA = A - Amean
    
#     # Dummy values (not used in PF but returned for interface compatibility)
#     dD = np.zeros((data.shape[0], ens_size))
#     MeasAvg = np.mean(observation, axis=1, keepdims=True) if observation.shape[0] < observation.shape[1] else np.mean(observation, axis=0, keepdims=True)
#     S = np.zeros_like(observation)
#     COV = data_cov
#     B = np.zeros_like(dD)
#     dAS = np.zeros((sim_size + par_size, data.shape[0]))
    
#     Analysis = result['Analysis']  # (sim_size, ens_size)
#     Analysisparam = result['Analysisparam']  # (par_size, ens_size)
    
#     return [A, Amean, dA, dD, MeasAvg, S, COV, B, dAS, Analysis, Analysisparam]


# if __name__ == "__main__":
#     example_usage()