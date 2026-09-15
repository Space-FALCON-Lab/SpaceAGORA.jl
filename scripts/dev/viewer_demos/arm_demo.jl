# Robot arm on an orbiting bus: the Robot_Arm_Planner_Cloth_Demo single-arm case
# (q_start -> target through the cloth planner) coupled into run_simulation.
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DEMO_OUT_ROOT = get(ENV, "SPACEAGORA_VIEWER_DEMO_OUT", joinpath(REPO_ROOT, "output", "viewer_demos"))
using SpaceAGORA, StaticArrays, Arrow, DataFrames, JSON
import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft
const SM = SpaceAGORA.SimulationModel
planet = make_no_gram_planet(:earth)
sc = make_three_body_spacecraft(bus_dims=(1.2, 1.2, 1.5), panel_dims=(0.01, 1.4, 0.8), bus_mass=180.0, panel_mass_each=6.0, panel_offset_y=1.3,
    ic=SM.InitialCondition(ra=planet.Rp_e + 410e3, rp=planet.Rp_e + 400e3, i=51.6, ω=0.0, Ω=30.0, ν=0.0), prop_mass=0.0, id=1)
arm = SM.default_cloth_arm_model(link_lengths_m=(1.2, 1.0, 0.7), link_radii_m=(0.07, 0.06, 0.05), link_masses_kg=(8.0, 5.0, 2.5), mount_offset_body=(0.6, 0.0, 0.75))
base = SM.ClothArmBasePose(SVector{3, Float64}(0.0, 0.0, 0.0), SVector{4, Float64}(0.0, 0.0, 0.0, 1.0))
# The planner targets an end-effector position; the demo's goal joint angles give it through FK.
target = SM.cloth_fk(arm, base, [-0.12, -0.08, 0.06]).end_effector_position
plan = SM.plan_robot_arm_motion(arm, base, [0.08, 0.95, -0.85], target; config=SM.RobotArmPlannerConfig(dt_s=0.1, duration_s=40.0))
eff = SM.RobotArmControlEffector(plan=plan, spacecraft_idx=1, controller=SM.init_robot_arm_joint_mpc(plan; dt_s=0.1, horizon=8), control_dt_s=0.1)
outdir = joinpath(DEMO_OUT_ROOT, "arm_demo"); rm(outdir; force=true, recursive=true)
basecfg = make_example_config(planet=planet, spacecraft=sc, mission_time=60.0,
    initial_time=SM.InitialTime(year=2024, month=3, day=1, hour=12, minute=0, second=0.0),
    dynamic_effectors=(SM.InverseSquaredJ2GravityModel(),), density_model=SM.NoAtmosphereModel(), ephemerides_model=SM.SimpleEphemeridesModel(),
    orientation_sim=true, keplerian=true, EI_km=120.0, verbose=false, results=true, results_directory=outdir)
args = SM.SimulationConfiguration(file_paths=basecfg.file_paths, simulation_settings=basecfg.simulation_settings,
    mission_configuration=SM.MissionConfiguration(mission_type=basecfg.mission_configuration.mission_type, keplerian=true, number_of_orbits=1,
        mission_time=60.0, orientation_sim=true, num_steps_to_save=1000, data_rate=0.1),
    environment_model=basecfg.environment_model, dynamics_model=basecfg.dynamics_model, guidance_model=basecfg.guidance_model,
    navigation_model=basecfg.navigation_model, control_model=SM.ControlModel(control_effectors=(eff,), control_rates=[0.1]),
    initial_time=basecfg.initial_time, integration_tolerances=basecfg.integration_tolerances)
@time run_simulation(args; visualization=true)
df = DataFrame(Arrow.Table(joinpath(outdir, "simulation_results.feather")))
arm_cols = filter(n -> startswith(n, "sc1_arm_pose"), names(df))
println("rows=", nrow(df), " arm cols=", length(arm_cols), " link cols=", count(startswith("sc1_link_pose"), names(df)))
println("arm link1 rel pos first/last: ", collect(df[1, arm_cols[1:3]]), " ", collect(df[end, arm_cols[1:3]]))
println("arm link3 rel pos first/last: ", collect(df[1, arm_cols[15:17]]), " ", collect(df[end, arm_cols[15:17]]))
scene = read_visualization_scene(joinpath(outdir, "simulation_results_scene.json"))
println("arm geometry: ", scene.spacecraft[1].arm === nothing ? "none" : (length(scene.spacecraft[1].arm.links), scene.spacecraft[1].arm.reach_m), " bounding=", scene.spacecraft[1].bounding_radius_m)
html = export_visualization(joinpath(outdir, "simulation_results"); max_frames=1200, trail_orbits=0.02, title="AGORA robot arm · planner motion on an orbiting bus", texture_resolution="4k")
println("html: ", html, " ", filesize(html))
s = read(html, String); i = findfirst("window.SPACEAGORA_VIEWER = ", s); j = findnext(";\n</script>", s, last(i))
payload = JSON.parse(s[last(i)+1:first(j)-1]); ap = payload["frames"]["arm_pose"]
println("payload arm_pose: ", ap === nothing ? "none" : (ap["counts"], ap["total"]))
