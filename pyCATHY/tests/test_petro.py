import numpy as np
import matplotlib.pyplot as plt
from pyCATHY import petro 

# VGP_predict_CATHY = petro.predict_unsat_soil_hydr_param(data=[[80,10,10]])


# def predict_unsat_soil_hydr_param(data=[[20, 20, 60]]):
#     """
#     Predict Van Genuchten hydraulic parameters using Rosetta pedotransfer function
#     and convert them to CATHY-compatible parameters.
#     """
#     from rosetta import rosetta, SoilData
    
#     soildata = SoilData.from_array(data)
#     mean, stdev, codes = rosetta(3, soildata)
    
#     theta_r = float(mean[:, 0])
#     theta_s = float(mean[:, 1])
#     alpha = 10 ** float(mean[:, 2])  # 1/cm
#     n = 10 ** float(mean[:, 3])
#     ksat = 10 ** float(mean[:, 4])  # cm/day
    
#     VGP_predict_CATHY = {}
#     VGP_predict_CATHY["POROS"] = theta_s
#     VGP_predict_CATHY["VGNCELL"] = n
#     VGP_predict_CATHY["VGM"] = 1 - 1 / n
#     VGP_predict_CATHY["VGRMCCELL"] = theta_r
#     VGP_predict_CATHY["VGPSATCELL"] = -1 / alpha / 100.0  # cm to m
#     VGP_predict_CATHY["PERMX"] = ksat * (1e-2 / 86400)
#     VGP_predict_CATHY["PERMY"] = ksat * (1e-2 / 86400)
#     VGP_predict_CATHY["PERMZ"] = ksat * (1e-2 / 86400)
#     VGP_predict_CATHY["ELSTOR"] = 1e-5
#     VGP_predict_CATHY["alpha"] = alpha
    
#     return VGP_predict_CATHY


def van_genuchten_theta(psi, theta_r, theta_s, alpha, n):
    """
    Van Genuchten water retention curve: theta(psi)
    
    Parameters
    ----------
    psi : array
        Pressure head [cm] (negative for unsaturated)
    theta_r : float
        Residual water content [-]
    theta_s : float
        Saturated water content [-]
    alpha : float
        Van Genuchten parameter [1/cm]
    n : float
        Van Genuchten parameter [-]
    
    Returns
    -------
    theta : array
        Volumetric water content [-]
    """
    m = 1 - 1/n
    Se = 1 / (1 + (alpha * np.abs(psi))**n)**m
    theta = theta_r + (theta_s - theta_r) * Se
    return theta


# def plot_multiple_soil_types(soil_dict):
#     """
#     Plot water retention curves for multiple soil types in subplots
    
#     Parameters
#     ----------
#     soil_dict : dict
#         Dictionary with soil names as keys and [sand, silt, clay] as values
#         Example: {'Sand': [85, 10, 5], 'Loam': [40, 40, 20]}
#     """
#     n_soils = len(soil_dict)
    
#     # Create subplots (2 rows x n_soils columns)
#     fig, axes = plt.subplots(2, n_soils, figsize=(6*n_soils, 10))
    
#     # Handle single soil case
#     if n_soils == 1:
#         axes = axes.reshape(2, 1)
    
#     # Pressure head range
#     psi = np.logspace(0, 5, 200)  # 1 to 10^5 cm
    
#     for idx, (name, texture) in enumerate(soil_dict.items()):
#         # Predict parameters
#         VGP = petro.predict_unsat_soil_hydr_param([texture])
        
#         theta_r = VGP["VGRMCCELL"]
#         theta_s = VGP["POROS"]
#         # alpha = VGP["VGPSATCELL"]
#         alpha = -1.0 / (100.0 * VGP["VGPSATCELL"])  # cm⁻¹
#         n = VGP["VGNCELL"]
        
#         # Calculate water content
#         theta = van_genuchten_theta(-psi, theta_r, theta_s, alpha, n)
        
#         # Top row: Semi-log plot
#         axes[0, idx].semilogx(psi, theta, 'b-', linewidth=2.5)
#         axes[0, idx].axhline(y=theta_s, color='r', linestyle='--', 
#                             linewidth=1.5, alpha=0.7, label=f'θs={theta_s:.3f}')
#         axes[0, idx].axhline(y=theta_r, color='g', linestyle='--', 
#                             linewidth=1.5, alpha=0.7, label=f'θr={theta_r:.3f}')
#         axes[0, idx].set_xlabel('|Pressure Head| (cm)', fontsize=11)
#         axes[0, idx].set_ylabel('Water Content θ (-)', fontsize=11)
#         axes[0, idx].set_title(f'{name}\n[S/Si/C: {texture[0]}/{texture[1]}/{texture[2]}%]', 
#                               fontsize=12, fontweight='bold')
#         axes[0, idx].grid(True, alpha=0.3)
#         axes[0, idx].legend(loc='best', fontsize=9)
#         axes[0, idx].set_xlim([1, 1e5])
#         axes[0, idx].set_ylim([0, max(theta_s * 1.1, 0.6)])
        
#         # Bottom row: Linear-log plot
#         axes[1, idx].plot(np.log10(psi), theta, 'b-', linewidth=2.5)
#         axes[1, idx].axhline(y=theta_s, color='r', linestyle='--', 
#                             linewidth=1.5, alpha=0.7)
#         axes[1, idx].axhline(y=theta_r, color='g', linestyle='--', 
#                             linewidth=1.5, alpha=0.7)
#         axes[1, idx].set_xlabel('log₁₀|Pressure Head| (log₁₀ cm)', fontsize=11)
#         axes[1, idx].set_ylabel('Water Content θ (-)', fontsize=11)
#         axes[1, idx].grid(True, alpha=0.3)
#         axes[1, idx].set_ylim([0, max(theta_s * 1.1, 0.6)])
        
#         # Add parameter text box
#         param_text = f'n = {n:.3f}\nm = {1-1/n:.3f}\nα = {alpha:.4f} cm⁻¹'
#         axes[1, idx].text(0.05, 0.95, param_text, transform=axes[1, idx].transAxes,
#                          fontsize=10, verticalalignment='top',
#                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.6))
    
#     plt.tight_layout()
#     return fig

def plot_comparison_single_panel_swapped(soil_dict):
    """
    Plot all soil types on the same axes for comparison
    with water content θ on the x-axis and pressure head on the y-axis.

    Parameters
    ----------
    soil_dict : dict
        Dictionary with soil names as keys and [sand, silt, clay] as values
    """
    import matplotlib.pyplot as plt
    import numpy as np

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Pressure head in meters (absolute value)
    psi = np.logspace(0, 6, 200) / 100.0  # cm → m
    colors = plt.cm.tab10(np.linspace(0, 1, len(soil_dict)))
    
    for idx, (name, texture) in enumerate(soil_dict.items()):
        VGP = petro.predict_unsat_soil_hydr_param([texture])
        
        theta_r = VGP["VGRMCCELL"]
        theta_s = VGP["POROS"]
        # convert VGPSATCELL (m) to alpha in m^-1
        alpha = -1.0 / VGP["VGPSATCELL"]  
        n = VGP["VGNCELL"]
        
        # Van Genuchten water retention
        theta = van_genuchten_theta(-psi, theta_r, theta_s, alpha, n)
        
        label = f'{name} (S/Si/C: {texture[0]}/{texture[1]}/{texture[2]}%)'
        
        # Swapped axes: θ on x, |psi| on y
        ax1.semilogy(theta, np.abs(psi), linewidth=2.5, color=colors[idx], label=label)
        ax2.plot(theta, np.log10(np.abs(psi)), linewidth=2.5, color=colors[idx], label=label)
    
    # Linear-log plot
    ax1.set_xlabel('Water Content θ (-)', fontsize=12)
    ax1.set_ylabel('|Pressure Head| (m)', fontsize=12)
    ax1.set_title('Water Retention Curves - θ vs Pressure Head', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='best', fontsize=9)
    
    # Linear-log plot
    ax2.set_xlabel('Water Content θ (-)', fontsize=12)
    ax2.set_ylabel('log10 |Pressure Head| (m)', fontsize=12)
    ax2.set_title('Water Retention Curves - θ vs log10(Pressure Head)', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best', fontsize=9)
    
    plt.tight_layout()
    return fig



def plot_comparison_single_panel(soil_dict):
    """
    Plot all soil types on the same axes for comparison
    
    Parameters
    ----------
    soil_dict : dict
        Dictionary with soil names as keys and [sand, silt, clay] as values
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    psi = np.logspace(0, 5, 200)
    colors = plt.cm.tab10(np.linspace(0, 1, len(soil_dict)))
    
    for idx, (name, texture) in enumerate(soil_dict.items()):
        VGP = petro.predict_unsat_soil_hydr_param([texture])
        
        theta_r = VGP["VGRMCCELL"]
        theta_s = VGP["POROS"]
        alpha = -1.0 / (100.0 * VGP["VGPSATCELL"])  # cm⁻¹
        n = VGP["VGNCELL"]
        
        theta = van_genuchten_theta(-psi, theta_r, theta_s, alpha, n)
        
        label = f'{name} (S/Si/C: {texture[0]}/{texture[1]}/{texture[2]}%)'
        
        ax1.semilogx(psi, theta, linewidth=2.5, color=colors[idx], label=label)
        ax2.plot(np.log10(psi), theta, linewidth=2.5, color=colors[idx], label=label)
    
    ax1.set_xlabel('|Pressure Head| (cm)', fontsize=12)
    ax1.set_ylabel('Water Content θ (-)', fontsize=12)
    ax1.set_title('Water Retention Curves - Comparison', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='best', fontsize=9)
    ax1.set_xlim([1, 1e5])
    
    ax2.set_xlabel('log₁₀|Pressure Head| (log₁₀ cm)', fontsize=12)
    ax2.set_ylabel('Water Content θ (-)', fontsize=12)
    ax2.set_title('Water Retention Curves - Comparison (Linear-Log)', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best', fontsize=9)
    
    plt.tight_layout()
    return fig


# Example usage
if __name__ == "__main__":
    
    # Define multiple soil types
    soil_types = {
        'Sand': [85, 10, 5],
        'Loam': [40, 40, 20],
        'Clay': [20, 20, 60],
        'Silty Clay': [10, 55, 35]
    }
    
    # Plot 1: Individual subplots for each soil type
    # print("Creating individual subplots for each soil type...")
    # fig1 = plot_multiple_soil_types(soil_types)
    # plt.show()
    
    # Plot 2: All soils on same axes for comparison
    print("Creating comparison plot...")
    fig2 = plot_comparison_single_panel_swapped(soil_types)
    plt.show()
    
    # # Example: Single soil analysis
    # print("\nCATHY Hydraulic Parameters for Clay soil:")
    # print("-" * 50)
    # VGP = petro.predict_unsat_soil_hydr_param([[20, 20, 60]])
    # for key, value in VGP.items():
    #     if key != 'alpha':
    #         print(f"{key:15s} = {value:.6e}")
    
    # # Single soil subplot
    # single_soil = {'Clay': [20, 20, 60]}
    # fig3 = plot_multiple_soil_types(single_soil)
    # plt.show()