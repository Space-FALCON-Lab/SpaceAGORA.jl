include(joinpath(dirname(@__DIR__), "analysis", "mars_hard_lander_prelim", "MarsHardLanderPrelim.jl"))

using .MarsHardLanderPrelim

const SHIELD_BODY_FLAP_OUTPUT_ROOT = "/Users/josephine/Research/Space-FALCON/_Proposals/[FY27]SURP/SHIELD_BodyFlap_Study"

geom = MarsHardLanderPrelim.shield_published_surrogate_geometry()

cfg = MarsHardLanderPrelim.default_config(
    output_root=SHIELD_BODY_FLAP_OUTPUT_ROOT,
    atmosphere_mode=:gram,
    run_secondary_sweeps=false,
    include_reverse_bangbang_sweep=false,
    enable_crossrange_proxy=false,
    enable_target_range_guidance=false,
    enable_monte_carlo=false,
    generate_plots=false,
    h_jettison_grid_m=[0.0],
    geometry=geom,
)

adapter = MarsHardLanderPrelim.build_atmosphere_adapter(cfg)

MarsHardLanderPrelim.run_shield_body_flap_study(
    cfg,
    adapter;
    entry_mass_kg=120.0,
    flap_area_fractions_of_stowed_ref=[0.01, 0.02, 0.05],
    deflection_grid_deg=[5.0, 10.0, 15.0],
)
