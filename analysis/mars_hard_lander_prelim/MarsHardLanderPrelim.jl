module MarsHardLanderPrelim

using CSV
using DataFrames
using GRAMSuite
using LinearAlgebra
using OrdinaryDiffEq
using Plots
using Printf
using SpaceAGORA
using StaticArrays
using Statistics

const SM = SpaceAGORA.SimulationModel
const AERO = SM.DynamicEffectors.AerodynamicEffectors
const ENV = SM.EnvironmentModels

include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "atmosphere.jl"))
include(joinpath(@__DIR__, "aero_states.jl"))
include(joinpath(@__DIR__, "entry_3dof.jl"))
include(joinpath(@__DIR__, "entry_3d_point_mass.jl"))
include(joinpath(@__DIR__, "sweeps.jl"))
include(joinpath(@__DIR__, "guidance.jl"))
include(joinpath(@__DIR__, "crossrange_proxy.jl"))
include(joinpath(@__DIR__, "alpha_sensitivity.jl"))
include(joinpath(@__DIR__, "crossrange_sensitivity.jl"))
include(joinpath(@__DIR__, "monte_carlo.jl"))
include(joinpath(@__DIR__, "trade_metrics.jl"))
include(joinpath(@__DIR__, "shield_trim_state.jl"))
include(joinpath(@__DIR__, "shield_rim_flap_study.jl"))
include(joinpath(@__DIR__, "shield_body_flap_study.jl"))
include(joinpath(@__DIR__, "shield_shoulder_strake_study.jl"))
include(joinpath(@__DIR__, "plots.jl"))
include(joinpath(@__DIR__, "workflow.jl"))

export VehicleGeometryConfig,
    PrelimConfig,
    PrelimResults,
    default_config,
    shield_published_surrogate_geometry,
    build_atmosphere_adapter,
    calibrate_aero_case,
    derive_fixed_mass_beta_targets,
    derive_fixed_mass_beta_targets_from_deployed_geometry,
    deployed_drag_skirt_equivalent_area,
    run_shield_state3_trim_study,
    run_shield_rim_flap_study,
    run_shield_body_flap_study,
    run_shield_body_flap_stability_study,
    run_shield_shoulder_strake_stability_study,
    run_shield_shoulder_strake_trajectory_study,
    solve_entry_trajectory,
    find_switch_altitude_for_target_range,
    run_prelim,
    main

end
