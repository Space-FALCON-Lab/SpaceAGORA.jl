# GMAT-CYGNSS 48hr Comparison Script

This script compares CYGNSS actual telemetry data with GMAT simulation results for a 48-hour mission profile.

## Usage

```bash
python gmat_cygnss_comp.py
```

The script will:
1. Load CYGNSS actual telemetry from `CYGNSS/cygnss_data_48hr.feather`
2. Load GMAT simulation results from `GMAT_Examples/Sim_CYGNSS_Comparison.csv`
3. Interpolate GMAT data to match CYGNSS measurement timestamps
4. Generate comparison plots and print summary statistics

## Generated Plots

All plots are saved to the `plots/` directory:

- **position_comparison.png**: X, Y, Z position comparison with RMSE statistics
- **velocity_comparison.png**: X, Y, Z velocity comparison with RMSE statistics  
- **position_error.png**: Magnitude of position error over time
- **orbital_elements_comparison.png**: Semi-major axis, eccentricity, and inclination comparison

## Output Statistics

The script prints a summary table containing:
- Position error metrics (RMSE, Mean, Max, Min) in km
- Velocity error metrics in m/s
- Data quality information (number of points, mission duration)

## Dependencies

- pandas
- numpy
- matplotlib
- pyarrow (for feather file support)

Install with: `pip install pyarrow`

## Data Files

**CYGNSS cygnss_data_48hr.feather columns:**
- `time`: Time offset in seconds
- `pos_ii_1`, `pos_ii_2`, `pos_ii_3`: Position in meters (ECI frame)
- `vel_ii_1`, `vel_ii_2`, `vel_ii_3`: Velocity in m/s (ECI frame)

**GMAT Sim_CYGNSS_Comparison.csv columns:**
- DefaultSC.ElapsedSecs: Elapsed time in seconds
- DefaultSC.EarthMJ2000Eq.{X,Y,Z}: Position in km (ECI)
- DefaultSC.EarthMJ2000Eq.{VX,VY,VZ}: Velocity in km/s (ECI)
- DefaultSC.Earth.SMA: Semi-major axis in km
- DefaultSC.Earth.ECC: Eccentricity
- DefaultSC.EarthMJ2000Eq.INC: Inclination in degrees

## Performance Notes

With ~172,797 CYGNSS data points interpolated against 17,281 GMAT snapshots:
- Position RMSE typically ~11 km (depends on model fidelity)
- Velocity RMSE typically ~12 m/s
- Execution time: ~20-30 seconds on typical hardware

