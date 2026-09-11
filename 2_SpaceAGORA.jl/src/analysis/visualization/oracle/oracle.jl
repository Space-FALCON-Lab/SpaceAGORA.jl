module OracleAnalysis

# Analysis/diagnostic functions exported for use by example runners and user code.
export _FlatSol, _make_flat_sol, _make_flat_sol_from_feather
export _orbit_count_from_flat_sol
export orbital_energy, total_momentum, angular_momentum
export evaluate_laser_exchanges, laser_exchanges_time_series
export _write_csv!, _write_feather!, _case_id, _print_summary, _scenario_stem
export _target_rv_at, _build_timeseries_dataframe
export elements_time_series, print_initial_final_elements
export r_v_a_RTN_time_series, a_gravity_laser_RTN_time_series
export delta_v_RTN_time_series, laser_force_RTN_time_series, link_status_time_series
# plot functions
export plot_orbits, plot_orbit_energy, plot_orbit_energy_individual_satellites
export plot_orbit_energy_total, plot_orbital_energy_change
export plot_angmom_two_axes, plot_momentum_two_axes
export plot_laser_dE_time_series, plot_laser_dP_time_series
export report_and_plot_r_RTN, report_and_plot_v_RTN, report_and_plot_a_RTN
export report_and_plot_dv_RTN, report_and_plot_F_RTN, report_and_plot_delta_P_RTN
export report_and_plot_OE, report_and_plot_OE_diff, report_and_plot_rp_ra

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

using LinearAlgebra
using StaticArrays
using DataFrames
using Arrow
using CSV
using Printf
using Plots
import Plots: plot, plot!, scatter!, savefig, twinx

using ..SimulationModel
using ..SimulationEngine

# _validate_options and _with are private helpers in SimulationModel (oracle_types.jl).
import ..SimulationModel: _validate_options, _with

include(joinpath(@__DIR__, "oracle_diagnostics.jl"))
include(joinpath(@__DIR__, "oracle_oe_analysis.jl"))
include(joinpath(@__DIR__, "oracle_plots.jl"))

end # module OracleAnalysis

