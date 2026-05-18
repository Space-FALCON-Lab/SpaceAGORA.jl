include(joinpath(dirname(@__DIR__), "analysis", "mars_hard_lander_prelim", "MarsHardLanderPrelim.jl"))

using .MarsHardLanderPrelim

const SHIELD_SHOULDER_STRAKE_TRAJECTORY_OUTPUT_ROOT = "/Users/josephine/Research/Space-FALCON/_Proposals/[FY27]SURP/SHIELD_ShoulderStrake_Trajectory_Study"

geom = MarsHardLanderPrelim.shield_published_surrogate_geometry()

cfg = MarsHardLanderPrelim.default_config(
    output_root=SHIELD_SHOULDER_STRAKE_TRAJECTORY_OUTPUT_ROOT,
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

MarsHardLanderPrelim.run_shield_shoulder_strake_trajectory_study(
    cfg,
    adapter;
    entry_mass_kg=120.0,
    cg_axial_fraction_of_body_length=0.60,
    static_margin_fraction_of_diameter=0.10,
    strake_area_fraction_of_stowed_ref=0.05,
    control_deflection_deg=15.0,
    trim_check_angle_deg=5.0,
    cl_alpha_per_rad=3.5,
    cd0=0.05,
    induced_drag_factor=0.15,
    cl_max=1.2,
    pitch_command_grid=collect(-1.0:0.1:1.0),
    yaw_command_grid=[-1.0, 1.0],
)
