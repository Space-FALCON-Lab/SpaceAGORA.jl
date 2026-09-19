# examples/AGORA_Odyssey_Energy_Depletion_Control_Test.jl configuration, propagated for two orbits.
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
using SpaceAGORA, StaticArrays
import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft
const SM = SpaceAGORA.SimulationModel
planet = make_no_gram_planet(:mars)
ra = 28_559.615e3; rp = planet.Rp_e + 77e3
period = 2pi * sqrt(((ra + rp) / 2)^3 / planet.μ)
sc = make_three_body_spacecraft(bus_dims=(2.2, 2.6, 1.7), panel_dims=(0.01, 3.76 / 2.0, 1.93), bus_mass=391.0, panel_mass_each=10.0,
    panel_offset_y=2.6 / 2.0 + 3.76 / 4.0, ic=SM.InitialCondition(ra=ra, rp=rp, i=93.6, ω=109.7454, Ω=28.1517, ν=180.0),
    reflection_coefficient=0.9, prop_mass=50.0, id=101)
cfg = SM.AerobrakingEnergyDepletionConfig(guidance_modes=(:targeting, :max_energy_depletion), max_energy_submodes=(:heat_rate, :heat_load),
    heat_load_switch_solver=:tpbvp_integration, controlled_panel_links=(2, 3), target_apoapsis_radius_m=26_750e3,
    max_alpha_rad=pi / 2, min_alpha_rad=1e-4, heat_rate_limit_w_cm2=0.15, heat_load_limit_j_cm2=30.0, structural_load_limit_pa=0.5)
st = SM.AerobrakingEnergyDepletionState(num_sats=1)
for link in sc.links; link.α = cfg.min_alpha_rad; end
include(joinpath(@__DIR__, "demo_options.jl"))
options = viewer_demo_options("odyssey_two_orbits", 2 * period)
outdir = options.output_dir
base = make_example_config(planet=planet, spacecraft=sc, mission_time=options.duration_s,
    initial_time=SM.InitialTime(year=2001, month=11, day=6, hour=19, minute=0, second=0.0),
    dynamic_effectors=(SM.InverseSquaredGravityModel(), SM.AerodynamicCoefficientfM()), density_model=SM.ExponentialAtmosphereModel(planet),
    ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=false, keplerian=true, EI_km=160.0, verbose=false,
    results=true, results_directory=outdir)
# Two periods require over a million 0.1-second control stops. Allow those
# accepted steps plus a margin for subdivisions/rejections, without weakening tolerances.
solver = SM.SolverConfig(maxiters=2 * ceil(Int, options.duration_s / 0.1) + 1000)
args = SM.SimConfig._with_configuration(base; solver_config=solver,
    mission_configuration=SM.MissionConfiguration(mission_type=base.mission_configuration.mission_type, keplerian=true, number_of_orbits=1,
        mission_time=options.duration_s, orientation_sim=false, num_steps_to_save=1000, data_rate=2.0),
    guidance_model=SM.GuidanceModel(guidance_effectors=(SM.AerobrakingEnergyDepletionGuidanceModel(cfg, st),), guidance_rates=[3.0]),
    control_model=SM.ControlModel(control_effectors=(SM.AerobrakingEnergyDepletionControlModel(cfg, st),), control_rates=[0.1]))
println("period_s=", period, " mission_s=", options.duration_s)
@time run_simulation(args; visualization=true)
html = export_visualization(joinpath(outdir, "simulation_results"); max_frames=6000, trail_orbits=1, title="AGORA Odyssey · energy-depletion control")
println("html: ", html, " ", filesize(html))
using Arrow, DataFrames
df = DataFrame(Arrow.Table(joinpath(outdir, "simulation_results.feather")))
println("rows=", nrow(df), " alt_min=", minimum(df.sc1_altitude), " periapsis passes (alt<160km segments)=", count(i -> df.sc1_altitude[i] < 160e3 && df.sc1_altitude[i-1] >= 160e3, 2:nrow(df)))
println("sampled maximum radius km: ", maximum(sqrt.(df.sc1_pos_1.^2 .+ df.sc1_pos_2.^2 .+ df.sc1_pos_3.^2)) / 1e3)
println("panel qy changes=", count(!=(0.0), diff(df[!, "sc1_link_pose_5"])))
