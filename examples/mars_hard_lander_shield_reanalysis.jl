include(joinpath(dirname(@__DIR__), "analysis", "mars_hard_lander_prelim", "MarsHardLanderPrelim.jl"))

using .MarsHardLanderPrelim

const SHIELD_OUTPUT_ROOT = "/Users/josephine/Research/Space-FALCON/_Proposals/[FY27]SURP/SHIELD_Reanalysis_Preliminary_Results"
const SUBSONIC_SWITCH_STEP_M = 500.0
const REPRESENTATIVE_SWITCH_STEP_M = 1_000.0

function _floor_to_step(x_m::Real, step_m::Real)
    return step_m * floor(Float64(x_m) / step_m)
end

function _subsonic_switch_band(config, adapter, target_beta_high, target_beta_ratio)
    aero_case = MarsHardLanderPrelim.calibrate_aero_case(config, adapter, target_beta_high, target_beta_ratio)
    body_policy = MarsHardLanderPrelim.ControlPolicy("body_only", :body_only, NaN, NaN)
    body_result = MarsHardLanderPrelim.solve_entry_trajectory(
        config,
        adapter,
        aero_case,
        body_policy;
        save_trajectory=true,
    )
    subsonic_idx = findfirst(x -> x <= 1.0, body_result.trajectory.mach)
    subsonic_idx === nothing && error("Published-size SHIELD surrogate never becomes subsonic before impact.")

    subsonic_onset_m = body_result.trajectory.altitude_km[subsonic_idx] * 1e3
    switch_ceiling_m = max(0.0, _floor_to_step(subsonic_onset_m, SUBSONIC_SWITCH_STEP_M))
    h_switch_grid_m = collect(switch_ceiling_m:-SUBSONIC_SWITCH_STEP_M:0.0)
    if isempty(h_switch_grid_m) || h_switch_grid_m[end] != 0.0
        push!(h_switch_grid_m, 0.0)
    end

    representative_switches_m = collect(switch_ceiling_m:-REPRESENTATIVE_SWITCH_STEP_M:0.0)
    if isempty(representative_switches_m) || representative_switches_m[end] != 0.0
        push!(representative_switches_m, 0.0)
    end

    return (
        subsonic_onset_m=subsonic_onset_m,
        switch_ceiling_m=switch_ceiling_m,
        h_switch_grid_m=h_switch_grid_m,
        representative_switches_m=representative_switches_m,
    )
end

geom = MarsHardLanderPrelim.shield_published_surrogate_geometry()

probe_cfg = MarsHardLanderPrelim.default_config(
    output_root=SHIELD_OUTPUT_ROOT,
    atmosphere_mode=:gram,
    run_secondary_sweeps=false,
    include_reverse_bangbang_sweep=false,
    enable_crossrange_proxy=true,
    enable_target_range_guidance=true,
    enable_monte_carlo=true,
    monte_carlo_samples=20,
    generate_plots=true,
    h_switch_grid_m=collect(40e3:-2e3:0.0),
    h_jettison_grid_m=[0.0],
    representative_switches_m=[40e3, 32e3, 24e3, 16e3, 8e3, 0.0],
    geometry=geom,
)

adapter = MarsHardLanderPrelim.build_atmosphere_adapter(probe_cfg)
targets = MarsHardLanderPrelim.derive_fixed_mass_beta_targets_from_deployed_geometry(
    probe_cfg,
    adapter;
    entry_mass_kg=120.0,
)
switch_band = _subsonic_switch_band(probe_cfg, adapter, targets.beta_high_kg_m2, targets.beta_ratio)

cfg = MarsHardLanderPrelim.default_config(
    output_root=SHIELD_OUTPUT_ROOT,
    atmosphere_mode=:gram,
    run_secondary_sweeps=false,
    include_reverse_bangbang_sweep=false,
    enable_crossrange_proxy=true,
    enable_target_range_guidance=true,
    enable_monte_carlo=true,
    monte_carlo_samples=20,
    generate_plots=true,
    h_switch_grid_m=switch_band.h_switch_grid_m,
    h_jettison_grid_m=[0.0],
    beta_high_targets=[targets.beta_high_kg_m2],
    beta_ratios=[targets.beta_ratio],
    representative_beta_high=targets.beta_high_kg_m2,
    representative_beta_ratio=targets.beta_ratio,
    representative_h_jettison_m=0.0,
    representative_switches_m=switch_band.representative_switches_m,
    geometry=geom,
)

MarsHardLanderPrelim.run_prelim(cfg)
