import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from dataclasses import dataclass
from typing import Tuple, Optional

@dataclass
class ArchieParams:
    """Archie's law parameters"""
    a: float = 1.0  # Tortuosity factor
    m: float = 2.0  # Cementation exponent
    n: float = 2.0  # Saturation exponent
    rho_w: float = 10.0  # Water resistivity (Ohm-m)

class DynamicModel:
    """Simple 1D flow model for saturation evolution"""
    
    def __init__(self, n_cells: int, dt: float = 0.1):
        self.n_cells = n_cells
        self.dt = dt
        
    def step(self, Sw: np.ndarray, phi: np.ndarray, K: np.ndarray, 
             injection_rate: float = 0.1) -> np.ndarray:
        """
        Simple advection-diffusion for saturation
        Sw: (n_cells,) saturation
        phi: (n_cells,) porosity
        K: (n_cells,) permeability
        """
        Sw_new = Sw.copy()
        
        # Simple left-to-right flow with diffusion
        for i in range(1, self.n_cells - 1):
            # Advection term (proportional to K)
            velocity = K[i] / phi[i]
            advection = -velocity * (Sw[i] - Sw[i-1])
            
            # Diffusion term
            diffusion = 0.01 * (Sw[i+1] - 2*Sw[i] + Sw[i-1])
            
            Sw_new[i] = Sw[i] + self.dt * (advection + diffusion)
        
        # Injection at left boundary
        Sw_new[0] += self.dt * injection_rate
        
        # Bound saturation [0, 1]
        Sw_new = np.clip(Sw_new, 0.01, 0.99)
        
        return Sw_new

def archie_forward(Sw: np.ndarray, phi: np.ndarray, 
                   params: ArchieParams) -> np.ndarray:
    """
    Archie's law: rho = a * phi^(-m) * Sw^(-n) * rho_w
    
    Args:
        Sw: (n_cells,) saturation
        phi: (n_cells,) porosity
        params: Archie parameters
    
    Returns:
        rho: (n_cells,) resistivity
    """
    return params.a * phi**(-params.m) * Sw**(-params.n) * params.rho_w

def ert_observation_operator(rho_cells: np.ndarray, 
                            obs_locations: np.ndarray,
                            cell_centers: np.ndarray) -> np.ndarray:
    """
    Simulate ERT observations as weighted averages of cell resistivities
    
    Args:
        rho_cells: (n_cells,) cell resistivities
        obs_locations: (n_obs,) observation locations
        cell_centers: (n_cells,) cell center locations
    
    Returns:
        rho_obs: (n_obs,) observed resistivities
    """
    n_obs = len(obs_locations)
    rho_obs = np.zeros(n_obs)
    
    for i, obs_loc in enumerate(obs_locations):
        # Gaussian sensitivity kernel
        distances = np.abs(cell_centers - obs_loc)
        weights = np.exp(-distances**2 / (2 * 0.5**2))  # sigma = 0.5
        weights /= weights.sum()
        
        rho_obs[i] = np.sum(weights * rho_cells)
    
    return rho_obs

def gaspari_cohn(r: np.ndarray, L: float) -> np.ndarray:
    """
    Gaspari-Cohn localization function
    
    Args:
        r: distance / localization radius
        L: localization radius
    
    Returns:
        rho: correlation values in [0, 1]
    """
    r_norm = r / L
    rho = np.zeros_like(r_norm)
    
    mask1 = (r_norm >= 0) & (r_norm < 1)
    mask2 = (r_norm >= 1) & (r_norm < 2)
    
    r1 = r_norm[mask1]
    rho[mask1] = 1 - 5/3*r1**2 + 5/8*r1**3 + 1/2*r1**4 - 1/4*r1**5
    
    r2 = r_norm[mask2]
    rho[mask2] = 4 - 5*r2 + 5/3*r2**2 + 5/8*r2**3 - 1/2*r2**4 + 1/12*r2**5 - 2/(3*r2)
    
    return rho

def logit_transform(x: np.ndarray) -> np.ndarray:
    """Transform bounded [0,1] to unbounded"""
    x_safe = np.clip(x, 1e-6, 1-1e-6)
    return np.log(x_safe / (1 - x_safe))

def inverse_logit(xi: np.ndarray) -> np.ndarray:
    """Transform unbounded to bounded [0,1]"""
    return 1 / (1 + np.exp(-xi))

class EnKF:
    """Ensemble Kalman Filter with localization and parameter transformation"""
    
    def __init__(self, n_ens: int, n_cells: int, obs_error_std: float,
                 localization_radius: float, use_transformation: bool = True):
        self.n_ens = n_ens
        self.n_cells = n_cells
        self.obs_error_std = obs_error_std
        self.R = obs_error_std**2 * np.eye(1)  # Will be adjusted for n_obs
        self.loc_radius = localization_radius
        self.use_transformation = use_transformation
        
    def update(self, ensemble_Sw: np.ndarray, ensemble_phi: np.ndarray,
               ensemble_K: np.ndarray, observations: np.ndarray,
               obs_locations: np.ndarray, cell_centers: np.ndarray,
               archie_params: ArchieParams) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        EnKF analysis update
        
        Args:
            ensemble_Sw: (n_ens, n_cells) saturation ensemble
            ensemble_phi: (n_ens, n_cells) porosity ensemble
            ensemble_K: (n_ens, n_cells) permeability ensemble
            observations: (n_obs,) ERT observations
            obs_locations: (n_obs,) observation locations
            cell_centers: (n_cells,) cell centers
            archie_params: Archie parameters
        
        Returns:
            Updated ensembles (Sw, phi, K)
        """
        n_obs = len(observations)
        self.R = self.obs_error_std**2 * np.eye(n_obs)
        
        # Transform parameters to unbounded space
        if self.use_transformation:
            ensemble_phi_trans = logit_transform(ensemble_phi)
            ensemble_K_trans = np.log(np.maximum(ensemble_K, 1e-6))
        else:
            ensemble_phi_trans = ensemble_phi.copy()
            ensemble_K_trans = ensemble_K.copy()
        
        # Predict observations for each ensemble member
        predicted_obs = np.zeros((self.n_ens, n_obs))
        for i in range(self.n_ens):
            rho_cells = archie_forward(ensemble_Sw[i], ensemble_phi[i], archie_params)
            predicted_obs[i] = ert_observation_operator(rho_cells, obs_locations, cell_centers)
        
        # Perturb observations
        obs_pert = observations[np.newaxis, :] + np.random.randn(self.n_ens, n_obs) * self.obs_error_std
        
        # Compute ensemble means
        mean_Sw = ensemble_Sw.mean(axis=0)
        mean_phi_trans = ensemble_phi_trans.mean(axis=0)
        mean_K_trans = ensemble_K_trans.mean(axis=0)
        mean_pred_obs = predicted_obs.mean(axis=0)
        
        # Compute anomalies
        A_Sw = ensemble_Sw - mean_Sw[np.newaxis, :]  # (n_ens, n_cells)
        A_phi = ensemble_phi_trans - mean_phi_trans[np.newaxis, :]  # (n_ens, n_cells)
        A_K = ensemble_K_trans - mean_K_trans[np.newaxis, :]  # (n_ens, n_cells)
        A_obs = predicted_obs - mean_pred_obs[np.newaxis, :]  # (n_ens, n_obs)
        
        # Stack state-parameter anomalies
        A_x = np.hstack([A_Sw, A_phi, A_K])  # (n_ens, 3*n_cells)
        
        # Compute covariances
        C_xy = (A_x.T @ A_obs) / (self.n_ens - 1)  # (3*n_cells, n_obs)
        C_yy = (A_obs.T @ A_obs) / (self.n_ens - 1)  # (n_obs, n_obs)
        
        # Apply localization
        if self.loc_radius > 0:
            loc_matrix = self._compute_localization_matrix(cell_centers, obs_locations)
            C_xy = C_xy * loc_matrix  # Element-wise multiplication (3*n_cells, n_obs)
        
        # Kalman gain
        K_gain = C_xy @ np.linalg.inv(C_yy + self.R)  # (3*n_cells, n_obs)
        
        # Update each ensemble member
        ensemble_Sw_new = ensemble_Sw.copy()
        ensemble_phi_new = ensemble_phi.copy()
        ensemble_K_new = ensemble_K.copy()
        
        for i in range(self.n_ens):
            innovation = obs_pert[i] - predicted_obs[i]
            
            # Stack state-parameters
            x_i = np.hstack([ensemble_Sw[i], ensemble_phi_trans[i], ensemble_K_trans[i]])
            
            # Update
            x_i_new = x_i + K_gain @ innovation
            
            # Extract updated components
            ensemble_Sw_new[i] = x_i_new[:self.n_cells]
            
            if self.use_transformation:
                ensemble_phi_new[i] = inverse_logit(x_i_new[self.n_cells:2*self.n_cells])
                ensemble_K_new[i] = np.exp(x_i_new[2*self.n_cells:])
            else:
                ensemble_phi_new[i] = x_i_new[self.n_cells:2*self.n_cells]
                ensemble_K_new[i] = x_i_new[2*self.n_cells:]
                # Apply hard bounds if not using transformation
                ensemble_phi_new[i] = np.clip(ensemble_phi_new[i], 0.01, 0.99)
                ensemble_K_new[i] = np.maximum(ensemble_K_new[i], 0.01)
            
            # Bound saturation
            ensemble_Sw_new[i] = np.clip(ensemble_Sw_new[i], 0.01, 0.99)
        
        return ensemble_Sw_new, ensemble_phi_new, ensemble_K_new
    
    def _compute_localization_matrix(self, cell_centers: np.ndarray,
                                    obs_locations: np.ndarray) -> np.ndarray:
        """
        Compute Gaspari-Cohn localization matrix
        
        Returns:
            loc_matrix: (3*n_cells, n_obs) localization weights
        """
        # Distance matrix between all cells and observations
        distances = cdist(cell_centers.reshape(-1, 1), 
                         obs_locations.reshape(-1, 1))  # (n_cells, n_obs)
        
        # Apply Gaspari-Cohn
        loc_weights = gaspari_cohn(distances, self.loc_radius)  # (n_cells, n_obs)
        
        # Repeat for Sw, phi, K (same spatial localization for all)
        loc_matrix = np.vstack([loc_weights, loc_weights, loc_weights])  # (3*n_cells, n_obs)
        
        return loc_matrix

class ParticleFilter:
    """Particle Filter for nonlinear Archie problem"""
    
    def __init__(self, n_particles: int, n_cells: int, obs_error_std: float,
                 localization_radius: float = 0.0, use_enkf_update: bool = False):
        self.n_particles = n_particles
        self.n_cells = n_cells
        self.obs_error_std = obs_error_std
        self.loc_radius = localization_radius
        self.use_enkf_update = use_enkf_update
        self.weights = np.ones(n_particles) / n_particles
        
    def update(self, ensemble_Sw: np.ndarray, ensemble_phi: np.ndarray,
               ensemble_K: np.ndarray, observations: np.ndarray,
               obs_locations: np.ndarray, cell_centers: np.ndarray,
               archie_params: ArchieParams) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Particle filter update with optional EnKF-style analysis
        
        Returns:
            Updated ensembles (Sw, phi, K)
        """
        n_obs = len(observations)
        
        # Predict observations for each particle
        predicted_obs = np.zeros((self.n_particles, n_obs))
        for i in range(self.n_particles):
            rho_cells = archie_forward(ensemble_Sw[i], ensemble_phi[i], archie_params)
            predicted_obs[i] = ert_observation_operator(rho_cells, obs_locations, cell_centers)
        
        # Compute weights using exact likelihood
        log_weights = np.zeros(self.n_particles)
        for i in range(self.n_particles):
            innovation = observations - predicted_obs[i]
            log_weights[i] = -0.5 * np.sum((innovation / self.obs_error_std)**2)
        
        # Normalize weights (subtract max for numerical stability)
        log_weights -= log_weights.max()
        self.weights = np.exp(log_weights)
        self.weights /= self.weights.sum()
        
        # Compute effective sample size
        n_eff = 1.0 / np.sum(self.weights**2)
        
        print(f"  N_eff = {n_eff:.1f} / {self.n_particles}")
        
        # Resample if needed
        if n_eff < self.n_particles / 2:
            print("  Resampling...")
            indices = self._systematic_resample()
            ensemble_Sw = ensemble_Sw[indices]
            ensemble_phi = ensemble_phi[indices]
            ensemble_K = ensemble_K[indices]
            self.weights = np.ones(self.n_particles) / self.n_particles
        
        # Optional: Apply weighted EnKF update
        if self.use_enkf_update:
            ensemble_Sw, ensemble_phi, ensemble_K = self._weighted_enkf_update(
                ensemble_Sw, ensemble_phi, ensemble_K,
                observations, obs_locations, cell_centers, archie_params
            )
        
        # Add jitter to prevent particle collapse
        ensemble_phi += np.random.randn(*ensemble_phi.shape) * 0.005
        ensemble_K *= np.exp(np.random.randn(*ensemble_K.shape) * 0.02)
        
        # Enforce bounds
        ensemble_phi = np.clip(ensemble_phi, 0.01, 0.99)
        ensemble_K = np.maximum(ensemble_K, 0.01)
        
        return ensemble_Sw, ensemble_phi, ensemble_K
    
    def _systematic_resample(self) -> np.ndarray:
        """Systematic resampling"""
        positions = (np.arange(self.n_particles) + np.random.rand()) / self.n_particles
        cumsum = np.cumsum(self.weights)
        indices = np.searchsorted(cumsum, positions)
        return indices
    
    def _weighted_enkf_update(self, ensemble_Sw, ensemble_phi, ensemble_K,
                             observations, obs_locations, cell_centers,
                             archie_params):
        """EnKF update using weighted covariances"""
        n_obs = len(observations)
        
        # Predict observations
        predicted_obs = np.zeros((self.n_particles, n_obs))
        for i in range(self.n_particles):
            rho_cells = archie_forward(ensemble_Sw[i], ensemble_phi[i], archie_params)
            predicted_obs[i] = ert_observation_operator(rho_cells, obs_locations, cell_centers)
        
        # Weighted means
        mean_Sw = np.sum(self.weights[:, np.newaxis] * ensemble_Sw, axis=0)
        mean_phi = np.sum(self.weights[:, np.newaxis] * ensemble_phi, axis=0)
        mean_K = np.sum(self.weights[:, np.newaxis] * ensemble_K, axis=0)
        mean_obs = np.sum(self.weights[:, np.newaxis] * predicted_obs, axis=0)
        
        # Weighted anomalies (centered using weighted mean)
        A_Sw = ensemble_Sw - mean_Sw[np.newaxis, :]  # (n_particles, n_cells)
        A_phi = ensemble_phi - mean_phi[np.newaxis, :]  # (n_particles, n_cells)
        A_K = ensemble_K - mean_K[np.newaxis, :]  # (n_particles, n_cells)
        A_obs = predicted_obs - mean_obs[np.newaxis, :]  # (n_particles, n_obs)
        
        # Stack into single matrix
        A_x = np.hstack([A_Sw, A_phi, A_K])  # (n_particles, 3*n_cells)
        
        # Weighted covariances
        W_sqrt = np.sqrt(self.weights)[:, np.newaxis]  # (n_particles, 1)
        A_x_weighted = A_x * W_sqrt  # (n_particles, 3*n_cells)
        A_obs_weighted = A_obs * W_sqrt  # (n_particles, n_obs)
        
        C_xy = A_x_weighted.T @ A_obs_weighted  # (3*n_cells, n_obs)
        C_yy = A_obs_weighted.T @ A_obs_weighted  # (n_obs, n_obs)
        
        # Kalman gain with inflation
        R = (self.obs_error_std**2) * np.eye(n_obs)
        K_gain = C_xy @ np.linalg.inv(C_yy + R)
        
        # Update
        for i in range(self.n_particles):
            innovation = observations - predicted_obs[i]
            x_i = np.hstack([ensemble_Sw[i], ensemble_phi[i], ensemble_K[i]])
            x_i_new = x_i + 0.5 * K_gain @ innovation  # Reduced gain for stability
            
            ensemble_Sw[i] = np.clip(x_i_new[:self.n_cells], 0.01, 0.99)
            ensemble_phi[i] = np.clip(x_i_new[self.n_cells:2*self.n_cells], 0.01, 0.99)
            ensemble_K[i] = np.maximum(x_i_new[2*self.n_cells:], 0.01)
        
        return ensemble_Sw, ensemble_phi, ensemble_K

def run_experiment(method: str = 'enkf', n_ensemble: int = 100,
                   use_transformation: bool = True, 
                   use_localization: bool = True):
    """
    Run data assimilation experiment
    
    Args:
        method: 'enkf', 'enkf_notrans', 'pf', or 'pf_enkf'
        n_ensemble: number of ensemble members/particles
        use_transformation: use logit/log transforms (EnKF only)
        use_localization: use covariance localization
    """
    # Setup
    np.random.seed(42)
    n_cells = 20
    n_steps = 10
    n_obs = 5
    dt = 0.1
    
    # True parameters
    phi_true = 0.2 + 0.1 * np.sin(np.linspace(0, 2*np.pi, n_cells))
    K_true = np.exp(0.5 * np.random.randn(n_cells))
    Sw_true = 0.3 * np.ones(n_cells)
    
    # Observation setup
    cell_centers = np.linspace(0, 1, n_cells)
    obs_locations = np.linspace(0.1, 0.9, n_obs)
    obs_error_std = 0.5
    
    # Archie parameters
    archie_params = ArchieParams()
    
    # Dynamic model
    model = DynamicModel(n_cells, dt)
    
    # Initialize ensemble
    ensemble_Sw = 0.3 + 0.05 * np.random.randn(n_ensemble, n_cells)
    ensemble_phi = 0.25 + 0.05 * np.random.randn(n_ensemble, n_cells)
    ensemble_K = np.exp(0.5 * np.random.randn(n_ensemble, n_cells))
    
    ensemble_Sw = np.clip(ensemble_Sw, 0.01, 0.99)
    ensemble_phi = np.clip(ensemble_phi, 0.05, 0.45)
    ensemble_K = np.maximum(ensemble_K, 0.01)
    
    # Initialize filter
    loc_radius = 0.3 if use_localization else 0.0
    
    if method == 'enkf':
        filter_obj = EnKF(n_ensemble, n_cells, obs_error_std, loc_radius, 
                         use_transformation=True)
    elif method == 'enkf_notrans':
        filter_obj = EnKF(n_ensemble, n_cells, obs_error_std, loc_radius,
                         use_transformation=False)
    elif method == 'pf':
        filter_obj = ParticleFilter(n_ensemble, n_cells, obs_error_std, 
                                    loc_radius, use_enkf_update=False)
    elif method == 'pf_enkf':
        filter_obj = ParticleFilter(n_ensemble, n_cells, obs_error_std,
                                    loc_radius, use_enkf_update=True)
    
    # Storage
    phi_mean_history = []
    phi_std_history = []
    K_mean_history = []
    rmse_history = []
    
    # Main DA loop
    print(f"\nRunning {method.upper()} with N={n_ensemble}...")
    for step in range(n_steps):
        print(f"Step {step+1}/{n_steps}")
        
        # Forecast: propagate ensemble
        for i in range(n_ensemble):
            ensemble_Sw[i] = model.step(ensemble_Sw[i], ensemble_phi[i], 
                                       ensemble_K[i], injection_rate=0.1)
        
        # Propagate truth
        Sw_true = model.step(Sw_true, phi_true, K_true, injection_rate=0.1)
        
        # Generate synthetic observations
        rho_true = archie_forward(Sw_true, phi_true, archie_params)
        obs_true = ert_observation_operator(rho_true, obs_locations, cell_centers)
        observations = obs_true + obs_error_std * np.random.randn(n_obs)
        
        # Analysis: update ensemble
        ensemble_Sw, ensemble_phi, ensemble_K = filter_obj.update(
            ensemble_Sw, ensemble_phi, ensemble_K,
            observations, obs_locations, cell_centers, archie_params
        )
        
        # Compute statistics
        phi_mean = ensemble_phi.mean(axis=0)
        phi_std = ensemble_phi.std(axis=0)
        K_mean = ensemble_K.mean(axis=0)
        
        rmse_phi = np.sqrt(np.mean((phi_mean - phi_true)**2))
        
        phi_mean_history.append(phi_mean)
        phi_std_history.append(phi_std)
        K_mean_history.append(K_mean)
        rmse_history.append(rmse_phi)
        
        print(f"  RMSE(phi) = {rmse_phi:.4f}")
        print(f"  phi range: [{ensemble_phi.min():.3f}, {ensemble_phi.max():.3f}]")
        print(f"  K range: [{ensemble_K.min():.3f}, {ensemble_K.max():.3f}]")
    
    # Plotting
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Porosity evolution
    ax = axes[0, 0]
    for i in range(min(10, n_ensemble)):
        ax.plot(cell_centers, ensemble_phi[i], 'b-', alpha=0.3, linewidth=0.5)
    ax.plot(cell_centers, phi_mean_history[-1], 'r-', linewidth=2, label='Mean')
    ax.plot(cell_centers, phi_true, 'k--', linewidth=2, label='True')
    ax.fill_between(cell_centers, 
                    phi_mean_history[-1] - phi_std_history[-1],
                    phi_mean_history[-1] + phi_std_history[-1],
                    alpha=0.3, color='red', label='±1 std')
    ax.set_xlabel('Spatial Location')
    ax.set_ylabel('Porosity')
    ax.set_title(f'Final Porosity - {method.upper()}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1])
    
    # Permeability evolution
    ax = axes[0, 1]
    for i in range(min(10, n_ensemble)):
        ax.plot(cell_centers, ensemble_K[i], 'b-', alpha=0.3, linewidth=0.5)
    ax.plot(cell_centers, K_mean_history[-1], 'r-', linewidth=2, label='Mean')
    ax.plot(cell_centers, K_true, 'k--', linewidth=2, label='True')
    ax.set_xlabel('Spatial Location')
    ax.set_ylabel('Permeability')
    ax.set_title(f'Final Permeability - {method.upper()}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    
    # RMSE evolution
    ax = axes[1, 0]
    ax.plot(range(1, n_steps+1), rmse_history, 'b-o', linewidth=2)
    ax.set_xlabel('Assimilation Step')
    ax.set_ylabel('RMSE (Porosity)')
    ax.set_title('Porosity RMSE Evolution')
    ax.grid(True, alpha=0.3)
    
    # Saturation field
    ax = axes[1, 1]
    for i in range(min(10, n_ensemble)):
        ax.plot(cell_centers, ensemble_Sw[i], 'b-', alpha=0.3, linewidth=0.5)
    ax.plot(cell_centers, ensemble_Sw.mean(axis=0), 'r-', linewidth=2, label='Mean')
    ax.plot(cell_centers, Sw_true, 'k--', linewidth=2, label='True')
    ax.set_xlabel('Spatial Location')
    ax.set_ylabel('Saturation')
    ax.set_title(f'Final Saturation - {method.upper()}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1])
    
    plt.tight_layout()
    plt.savefig(f'results_{method}.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    return rmse_history[-1]

# Run comparisons
if __name__ == "__main__":
    print("="*70)
    print("ERT-ARCHIE DATA ASSIMILATION COMPARISON")
    print("="*70)
    
    results = {}
    
    # 1. EnKF without transformation (will show the problem)
    print("\n1. EnKF WITHOUT transformation (expect issues):")
    results['enkf_notrans'] = run_experiment('enkf_notrans', n_ensemble=100, 
                                            use_transformation=False)
    
    # 2. EnKF with transformation (should work better)
    print("\n2. EnKF WITH transformation (logit/log):")
    results['enkf_trans'] = run_experiment('enkf', n_ensemble=100,
                                          use_transformation=True)
    
    # 3. Particle filter (handles nonlinearity naturally)
    print("\n3. Particle Filter (exact likelihood):")
    results['pf'] = run_experiment('pf', n_ensemble=100)
    
    # 4. Hybrid: PF weights + EnKF update
    print("\n4. Hybrid PF-EnKF:")
    results['pf_enkf'] = run_experiment('pf_enkf', n_ensemble=100)
    
    # Summary
    print("\n" + "="*70)
    print("FINAL RESULTS (Porosity RMSE):")
    print("="*70)
    for method, rmse in results.items():
        print(f"{method:20s}: {rmse:.4f}")
    print("="*70)