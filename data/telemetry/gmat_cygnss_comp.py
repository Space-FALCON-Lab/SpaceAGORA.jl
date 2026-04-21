"""
GMAT-CYGNSS 48hr Comparison Plotter

Compares actual CYGNSS telemetry from 48hr_PVT.csv with GMAT simulation 
results from Sim_CYGNSS_Comparison.csv.

Plots orbital elements, position errors, and velocity errors.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import os

# Configuration
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CYGNSS_DATA_FILE = os.path.join(SCRIPT_DIR, "CYGNSS/cygnss_data_48hr.feather")
GMAT_SIM_FILE = os.path.join(SCRIPT_DIR, "GMAT_Examples/Sim_CYGNSS_Comparison.csv")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "plots")

def load_cygnss_data():
    """Load CYGNSS 48hr telemetry data from feather format."""
    try:
        # Try pandas native feather reader first (uses pyarrow backend)
        df = pd.read_feather(CYGNSS_DATA_FILE)
    except Exception as e:
        raise IOError(f"Failed to load feather file {CYGNSS_DATA_FILE}: {e}")
    
    # Extract position and velocity columns
    # CYGNSS feather schema: time, pos_ii_{1,2,3}, vel_ii_{1,2,3}
    # Positions are in METERS, velocities in m/s
    cygnss_data = {
        'time_offset': df['time'].values,
        'pos_x': df['pos_ii_1'].values / 1000,  # Convert meters to km
        'pos_y': df['pos_ii_2'].values / 1000,
        'pos_z': df['pos_ii_3'].values / 1000,
        'vel_x': df['vel_ii_1'].values / 1000,  # Convert m/s to km/s
        'vel_y': df['vel_ii_2'].values / 1000,
        'vel_z': df['vel_ii_3'].values / 1000,
    }
    
    return cygnss_data

def load_gmat_simulation():
    """Load GMAT simulation results."""
    df = pd.read_csv(GMAT_SIM_FILE)
    
    # GMAT uses: DefaultSC.EarthMJ2000Eq coordinates (in km for position, km/s for velocity)
    gmat_data = {
        'time': df['DefaultSC.ElapsedSecs'].values,
        'pos_x': df['DefaultSC.EarthMJ2000Eq.X'].values,
        'pos_y': df['DefaultSC.EarthMJ2000Eq.Y'].values,
        'pos_z': df['DefaultSC.EarthMJ2000Eq.Z'].values,
        'vel_x': df['DefaultSC.EarthMJ2000Eq.VX'].values,
        'vel_y': df['DefaultSC.EarthMJ2000Eq.VY'].values,
        'vel_z': df['DefaultSC.EarthMJ2000Eq.VZ'].values,
        'sma': df['DefaultSC.Earth.SMA'].values,
        'ecc': df['DefaultSC.Earth.ECC'].values,
        'inc': df['DefaultSC.EarthMJ2000Eq.INC'].values,  # EarthMJ2000Eq inclination
    }
    
    return gmat_data

def interpolate_to_cygnss_times(gmat_data, cygnss_times):
    """Interpolate GMAT data to CYGNSS measurement times."""
    gmat_interp = {}
    
    for key in gmat_data:
        if key != 'time':
            gmat_interp[key] = np.interp(cygnss_times, gmat_data['time'], gmat_data[key])
    
    return gmat_interp

def compute_orbital_elements(pos, vel, mu=398600.4418):
    """
    Compute orbital elements from position and velocity vectors.
    
    Parameters:
    - pos: [x, y, z] position in km
    - vel: [vx, vy, vz] velocity in km/s
    - mu: GM of Earth in km^3/s^2
    """
    r = np.linalg.norm(pos)
    v = np.linalg.norm(vel)
    
    # Specific orbital energy
    energy = 0.5 * v**2 - mu / r
    
    # Semi-major axis
    a = -mu / (2 * energy)
    
    # Angular momentum
    h = np.cross(pos, vel)
    h_mag = np.linalg.norm(h)
    
    # Eccentricity
    node = np.cross([0, 0, 1], h)
    e_vec = (np.cross(vel, h) / mu) - (pos / r)
    e = np.linalg.norm(e_vec)
    
    # Inclination (radians)
    inc = np.arccos(h[2] / h_mag)
    
    return a, e, inc

def create_output_dir():
    """Create output directory for plots."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)

def plot_position_comparison(cygnss_data, gmat_interp, time_hours):
    """Plot position comparison between CYGNSS and GMAT."""
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    
    components = [('pos_x', 'X'), ('pos_y', 'Y'), ('pos_z', 'Z')]
    
    for idx, (col, label) in enumerate(components):
        ax = axes[idx]
        
        pos_cygnss = cygnss_data[col]
        pos_gmat = gmat_interp[col]
        diff = pos_cygnss - pos_gmat
        
        ax.plot(time_hours, pos_cygnss, 'b.-', label='CYGNSS Telemetry', linewidth=1, markersize=4)
        ax.plot(time_hours, pos_gmat, 'r--', label='GMAT Simulation', linewidth=1.5)
        ax.set_ylabel(f'Position {label} (km)', fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best')
        
        rmse = np.sqrt(np.mean(diff**2))
        max_err = np.max(np.abs(diff))
        ax.text(0.02, 0.95, f'RMSE={rmse:.2f} km, Max={max_err:.2f} km',
                transform=ax.transAxes, fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    axes[-1].set_xlabel('Time (hours)', fontsize=11)
    plt.suptitle('CYGNSS vs GMAT: Position Comparison (48hr)', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'position_comparison.png'), dpi=150)
    print("✓ Saved: position_comparison.png")
    plt.close()

def plot_velocity_comparison(cygnss_data, gmat_interp, time_hours):
    """Plot velocity comparison between CYGNSS and GMAT."""
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    
    components = [('vel_x', 'X'), ('vel_y', 'Y'), ('vel_z', 'Z')]
    
    for idx, (col, label) in enumerate(components):
        ax = axes[idx]
        
        vel_cygnss = cygnss_data[col]
        vel_gmat = gmat_interp[col]
        diff = vel_cygnss - vel_gmat
        
        ax.plot(time_hours, vel_cygnss, 'b.-', label='CYGNSS Telemetry', linewidth=1, markersize=4)
        ax.plot(time_hours, vel_gmat, 'r--', label='GMAT Simulation', linewidth=1.5)
        ax.set_ylabel(f'Velocity {label} (km/s)', fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best')
        
        rmse = np.sqrt(np.mean(diff**2)) * 1000  # Convert to m/s
        max_err = np.max(np.abs(diff)) * 1000
        ax.text(0.02, 0.95, f'RMSE={rmse:.2f} m/s, Max={max_err:.2f} m/s',
                transform=ax.transAxes, fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    axes[-1].set_xlabel('Time (hours)', fontsize=11)
    plt.suptitle('CYGNSS vs GMAT: Velocity Comparison (48hr)', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'velocity_comparison.png'), dpi=150)
    print("✓ Saved: velocity_comparison.png")
    plt.close()

def plot_position_error(cygnss_data, gmat_interp, time_hours):
    """Plot absolute position error."""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    pos_cygnss = np.array([cygnss_data['pos_x'], cygnss_data['pos_y'], cygnss_data['pos_z']])
    pos_gmat = np.array([gmat_interp['pos_x'], gmat_interp['pos_y'], gmat_interp['pos_z']])
    
    # Position error magnitude
    err_mag = np.linalg.norm(pos_cygnss - pos_gmat, axis=0)
    
    ax.plot(time_hours, err_mag, 'r-', linewidth=2, label='Position Error')
    ax.fill_between(time_hours, 0, err_mag, alpha=0.3, color='red')
    
    # Statistics
    rmse = np.sqrt(np.mean(err_mag**2))
    max_err = np.max(err_mag)
    mean_err = np.mean(err_mag)
    
    ax.set_ylabel('Position Error (km)', fontsize=12)
    ax.set_xlabel('Time (hours)', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11)
    
    stats_text = f'RMSE={rmse:.2f} km\nMax={max_err:.2f} km\nMean={mean_err:.2f} km'
    ax.text(0.02, 0.95, stats_text, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
    
    plt.title('CYGNSS vs GMAT: Position Error Magnitude (48hr)', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'position_error.png'), dpi=150)
    print("✓ Saved: position_error.png")
    plt.close()

def plot_orbital_elements_comparison(cygnss_data, gmat_data, gmat_interp, time_hours):
    """Plot orbital elements comparison."""
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    
    # Compute orbital elements for CYGNSS data
    cygnss_orbital = {
        'sma': [],
        'ecc': [],
        'inc': [],
    }
    
    for i in range(len(cygnss_data['pos_x'])):
        pos = np.array([cygnss_data['pos_x'][i], cygnss_data['pos_y'][i], cygnss_data['pos_z'][i]])
        vel = np.array([cygnss_data['vel_x'][i], cygnss_data['vel_y'][i], cygnss_data['vel_z'][i]])
        
        try:
            a, e, inc = compute_orbital_elements(pos, vel)
            cygnss_orbital['sma'].append(a)
            cygnss_orbital['ecc'].append(e)
            cygnss_orbital['inc'].append(np.degrees(inc))
        except:
            # If computation fails, use NaN
            cygnss_orbital['sma'].append(np.nan)
            cygnss_orbital['ecc'].append(np.nan)
            cygnss_orbital['inc'].append(np.nan)
    
    cygnss_orbital['sma'] = np.array(cygnss_orbital['sma'])
    cygnss_orbital['ecc'] = np.array(cygnss_orbital['ecc'])
    cygnss_orbital['inc'] = np.array(cygnss_orbital['inc'])
    
    # Interpolate GMAT orbital elements to CYGNSS times
    gmat_sma_interp = np.interp(cygnss_data['time_offset'], gmat_data['time'], gmat_data['sma'])
    gmat_ecc_interp = np.interp(cygnss_data['time_offset'], gmat_data['time'], gmat_data['ecc'])
    gmat_inc_interp = np.interp(cygnss_data['time_offset'], gmat_data['time'], gmat_data['inc'])
    
    # SMA comparison
    ax = axes[0]
    ax.plot(time_hours, cygnss_orbital['sma'], 'b.-', label='CYGNSS-derived', linewidth=1, markersize=4)
    ax.plot(time_hours, gmat_sma_interp, 'r--', label='GMAT Simulation', linewidth=1.5)
    ax.set_ylabel('Semi-Major Axis (km)', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best')
    
    sma_rmse = np.sqrt(np.nanmean((cygnss_orbital['sma'] - gmat_sma_interp)**2))
    ax.text(0.02, 0.95, f'RMSE={sma_rmse:.2f} km',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Eccentricity comparison
    ax = axes[1]
    ax.plot(time_hours, cygnss_orbital['ecc'], 'b.-', label='CYGNSS-derived', linewidth=1, markersize=4)
    ax.plot(time_hours, gmat_ecc_interp, 'r--', label='GMAT Simulation', linewidth=1.5)
    ax.set_ylabel('Eccentricity', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best')
    
    ecc_rmse = np.sqrt(np.nanmean((cygnss_orbital['ecc'] - gmat_ecc_interp)**2))
    ax.text(0.02, 0.95, f'RMSE={ecc_rmse:.1e}',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Inclination comparison
    ax = axes[2]
    ax.plot(time_hours, cygnss_orbital['inc'], 'b.-', label='CYGNSS-derived', linewidth=1, markersize=4)
    ax.plot(time_hours, gmat_inc_interp, 'r--', label='GMAT Simulation', linewidth=1.5)
    ax.set_ylabel('Inclination (degrees)', fontsize=11)
    ax.set_xlabel('Time (hours)', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best')
    
    inc_rmse = np.sqrt(np.nanmean((cygnss_orbital['inc'] - gmat_inc_interp)**2))
    ax.text(0.02, 0.95, f'RMSE={inc_rmse:.4f}°',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.suptitle('CYGNSS vs GMAT: Orbital Elements Comparison (48hr)', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'orbital_elements_comparison.png'), dpi=150)
    print("✓ Saved: orbital_elements_comparison.png")
    plt.close()

def print_summary_statistics(cygnss_data, gmat_interp, time_hours):
    """Print summary statistics."""
    print("\n" + "="*70)
    print("CYGNSS vs GMAT 48hr Comparison Summary")
    print("="*70)
    
    # Position errors
    pos_cygnss = np.array([cygnss_data['pos_x'], cygnss_data['pos_y'], cygnss_data['pos_z']])
    pos_gmat = np.array([gmat_interp['pos_x'], gmat_interp['pos_y'], gmat_interp['pos_z']])
    pos_err_mag = np.linalg.norm(pos_cygnss - pos_gmat, axis=0)
    
    print(f"\n{'Position Error (km)':.<40}")
    print(f"  {'RMSE':.<35} {np.sqrt(np.mean(pos_err_mag**2)):.2f}")
    print(f"  {'Mean':.<35} {np.mean(pos_err_mag):.2f}")
    print(f"  {'Max':.<35} {np.max(pos_err_mag):.2f}")
    print(f"  {'Min':.<35} {np.min(pos_err_mag):.2f}")
    
    # Velocity errors
    vel_cygnss = np.array([cygnss_data['vel_x'], cygnss_data['vel_y'], cygnss_data['vel_z']])
    vel_gmat = np.array([gmat_interp['vel_x'], gmat_interp['vel_y'], gmat_interp['vel_z']])
    vel_err_mag = np.linalg.norm(vel_cygnss - vel_gmat, axis=0) * 1000  # m/s
    
    print(f"\n{'Velocity Error (m/s)':.<40}")
    print(f"  {'RMSE':.<35} {np.sqrt(np.mean(vel_err_mag**2)):.2f}")
    print(f"  {'Mean':.<35} {np.mean(vel_err_mag):.2f}")
    print(f"  {'Max':.<35} {np.max(vel_err_mag):.2f}")
    print(f"  {'Min':.<35} {np.min(vel_err_mag):.2f}")
    
    # Data quality
    print(f"\n{'Data Quality':.<40}")
    print(f"  {'CYGNSS points':.<35} {len(cygnss_data['pos_x'])}")
    print(f"  {'GMAT points (interpolated)':.<35} {len(gmat_interp['pos_x'])}")
    print(f"  {'Mission duration':.<35} {time_hours[-1]:.1f} hours")
    
    print("="*70)

def main():
    """Main execution function."""
    print("Loading data...")
    cygnss_data = load_cygnss_data()
    gmat_data = load_gmat_simulation()
    
    print(f"  CYGNSS: {len(cygnss_data['pos_x'])} measurements")
    print(f"  GMAT:    {len(gmat_data['pos_x'])} snapshots")
    
    # Convert CYGNSS time offset to hours
    time_hours = cygnss_data['time_offset'] / 3600.0
    
    # Interpolate GMAT to CYGNSS times
    print("Interpolating GMAT data to CYGNSS measurement times...")
    gmat_interp = interpolate_to_cygnss_times(gmat_data, cygnss_data['time_offset'])
    
    # Create output directory
    create_output_dir()
    
    # Generate plots
    print("\nGenerating plots...")
    plot_position_comparison(cygnss_data, gmat_interp, time_hours)
    plot_velocity_comparison(cygnss_data, gmat_interp, time_hours)
    plot_position_error(cygnss_data, gmat_interp, time_hours)
    plot_orbital_elements_comparison(cygnss_data, gmat_data, gmat_interp, time_hours)
    
    # Print statistics
    print_summary_statistics(cygnss_data, gmat_interp, time_hours)
    
    print(f"\n✓ All plots saved to '{OUTPUT_DIR}/' directory")

if __name__ == "__main__":
    main()
