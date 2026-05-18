#!/usr/bin/env julia

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const PROPOSAL_OUTPUT_ROOT = normpath(
    "/Users/josephine/Research/Space-FALCON/_Proposals/[FY27]SURP/AERO-08_Preliminary_Results",
)
include(joinpath(REPO_ROOT, "analysis", "mars_hard_lander_prelim", "MarsHardLanderPrelim.jl"))

using .MarsHardLanderPrelim

cfg = default_config(
    output_root=PROPOSAL_OUTPUT_ROOT,
    atmosphere_mode=:gram,
    beta_high_targets=[150.0],
    beta_ratios=[4.0],
    h_jettison_grid_m=[20e3],
    h_switch_grid_m=collect(120e3:-10e3:0.0),
    representative_beta_high=150.0,
    representative_beta_ratio=4.0,
    representative_h_jettison_m=20e3,
    representative_switches_m=[120e3, 80e3, 40e3, 20e3, 0.0],
    monte_carlo_samples=12,
    generate_plots=true,
)

run_prelim(cfg)
