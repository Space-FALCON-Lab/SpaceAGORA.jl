using Test
using SpaceAGORA
using StaticArrays
using LinearAlgebra
using Arrow
using DataFrames
using JSON
using Base64

import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft

const SM = SpaceAGORA.SimulationModel
const SV = SM.SceneVisualization

# The Robot_Arm_Planner_Cloth_Demo single-arm case (q_start -> target through
# the cloth planner), mounted on an orbiting bus and coupled into the engine
# through RobotArmControlEffector, which is what puts `arm_r`/`arm_q` in the state.
function _arm_plan(; duration_s::Float64=10.0)
    arm = SM.default_cloth_arm_model(link_lengths_m=(0.9, 0.8, 0.6), link_radii_m=(0.06, 0.05, 0.04), link_masses_kg=(6.0, 4.0, 2.0), mount_offset_body=(0.5, 0.0, 0.6))
    base = SM.ClothArmBasePose(SVector{3, Float64}(0.0, 0.0, 0.0), SVector{4, Float64}(0.0, 0.0, 0.0, 1.0))
    # The planner targets an end-effector position; the demo's goal joint angles give it through FK.
    target = SM.cloth_fk(arm, base, [-0.12, -0.08, 0.06]).end_effector_position
    return SM.plan_robot_arm_motion(arm, base, [0.08, 0.95, -0.85], target; config=SM.RobotArmPlannerConfig(dt_s=0.1, duration_s=duration_s))
end

function _arm_config(; results_directory::String, mission_time::Float64=12.0, with_arm::Bool=true)
    planet = make_no_gram_planet(:earth)
    sc = make_three_body_spacecraft(
        bus_dims=(1.0, 1.0, 1.2), panel_dims=(0.01, 1.0, 0.6), bus_mass=120.0, panel_mass_each=4.0, panel_offset_y=1.0,
        ic=SM.InitialCondition(ra=planet.Rp_e + 410e3, rp=planet.Rp_e + 400e3, i=51.6, ω=0.0, Ω=30.0, ν=0.0), prop_mass=0.0, id=1
    )
    base = make_example_config(
        planet=planet, spacecraft=sc, mission_time=mission_time,
        initial_time=SM.InitialTime(year=2024, month=3, day=1, hour=12, minute=0, second=0.0),
        dynamic_effectors=(SM.InverseSquaredJ2GravityModel(),), density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=true, keplerian=true, EI_km=120.0,
        verbose=false, results=true, results_directory=results_directory
    )
    effectors = if with_arm
        plan = _arm_plan(duration_s=min(10.0, mission_time))
        (SM.RobotArmControlEffector(plan=plan, spacecraft_idx=1, controller=SM.init_robot_arm_joint_mpc(plan; dt_s=0.1, horizon=6), control_dt_s=0.1),)
    else
        ()
    end
    return SM.SimulationConfiguration(
        file_paths=base.file_paths, simulation_settings=base.simulation_settings,
        mission_configuration=SM.MissionConfiguration(
            mission_type=base.mission_configuration.mission_type, keplerian=true, number_of_orbits=1,
            mission_time=mission_time, orientation_sim=true, num_steps_to_save=1000, data_rate=0.2
        ),
        environment_model=base.environment_model, dynamics_model=base.dynamics_model, guidance_model=base.guidance_model,
        navigation_model=base.navigation_model,
        control_model=SM.ControlModel(control_effectors=effectors, control_rates=with_arm ? [0.1] : Float64[]),
        initial_time=base.initial_time, integration_tolerances=base.integration_tolerances
    )
end

_decode_f32(b64) = collect(reinterpret(Float32, base64decode(b64)))

@testset "Robot arm in the viewer" begin
    @testset "geometry and plan lookup" begin
        plan = _arm_plan()
        geom = arm_geometry(plan)
        @test length(geom.links) == 3
        @test geom.links[1].name == "link_1"
        @test geom.links[1].vector_m == SVector(0.9, 0.0, 0.0)
        @test geom.links[1].com_offset_m == SVector(0.45, 0.0, 0.0)
        @test geom.links[2].radius_m == 0.05 && geom.links[3].mass_kg == 2.0
        @test geom.mount_offset_body_m == SVector(0.5, 0.0, 0.6)
        @test geom.reach_m ≈ 0.9 + 0.8 + 0.6

        args = _arm_config(results_directory=mktempdir())
        @test SV.robot_arm_plan_for(args, 1) === plan || SV.robot_arm_plan_for(args, 1) isa SM.RobotArmPlan
        @test SV.robot_arm_plan_for(args, 2) === nothing
        @test SV.robot_arm_plan_for(_arm_config(results_directory=mktempdir(), with_arm=false), 1) === nothing

        scene = build_visualization_scene(with_visualization_scene(args, true); rotation_max_samples=4)
        g = scene.spacecraft[1]
        @test g.arm !== nothing && length(g.arm.links) == 3
        # The bounding radius grows to the arm's reach so the close-up view switches in early enough.
        @test g.bounding_radius_m ≈ norm(SVector(0.5, 0.0, 0.6)) + 2.3
        d = SV.scene_dict(scene)
        @test d["spacecraft"][1]["arm"]["reach_m"] ≈ 2.3
        @test d["results"]["arm_pose"]["stride"] == 7
        path = write_visualization_scene(joinpath(mktempdir(), "scene.json"), scene)
        @test read_visualization_scene(path) == scene
        plain = spacecraft_geometry(args.dynamics_model.spacecraft[1])
        @test plain.arm === nothing
        @test SV.scene_dict(build_visualization_scene(_arm_config(results_directory=mktempdir(), with_arm=false)))["spacecraft"][1]["arm"] === nothing
    end

    @testset "arm_pose vector from a state view" begin
        view = (arm_r=[1.0 2.0; 3.0 4.0; 5.0 6.0], arm_q=[0.0 0.0; 0.0 0.0; 0.0 1.0; 1.0 0.0])
        v = SV.arm_pose_vector(view, SVector(0.5, 0.5, 0.5))
        @test v == [0.5, 2.5, 4.5, 0.0, 0.0, 0.0, 1.0, 1.5, 3.5, 5.5, 0.0, 0.0, 1.0, 0.0]
        @test isempty(SV.arm_pose_vector((pos=zeros(3),), zeros(3)))
    end

    @testset "coupled run records the arm and the page carries it" begin
        dir = mktempdir()
        args = with_visualization_scene(_arm_config(results_directory=dir), true)
        names_on = Symbol[f.name for f in SM.default_save_fields(args)]
        @test :arm_pose in names_on && :link_pose in names_on && :density ∉ names_on
        run_simulation(args)
        df = DataFrame(Arrow.Table(joinpath(dir, "simulation_results.feather")))
        cols = ["sc1_arm_pose_$(k)" for k in 1:21]
        @test all(c -> c in names(df), cols)
        @test "sc1_arm_pose_22" ∉ names(df)
        # Relative positions stay within the arm's reach of the bus, and the arm moves.
        rel1 = [norm([df[r, "sc1_arm_pose_1"], df[r, "sc1_arm_pose_2"], df[r, "sc1_arm_pose_3"]]) for r in 1:nrow(df)]
        @test all(rel1 .< 3.0)
        tip_start = [df[1, "sc1_arm_pose_15"], df[1, "sc1_arm_pose_16"], df[1, "sc1_arm_pose_17"]]
        tip_end = [df[end, "sc1_arm_pose_15"], df[end, "sc1_arm_pose_16"], df[end, "sc1_arm_pose_17"]]
        @test norm(tip_end - tip_start) > 0.1
        for r in (1, nrow(df)), link in 0:2
            q = [df[r, "sc1_arm_pose_$(7 * link + c)"] for c in 4:7]
            @test norm(q) ≈ 1.0 atol=1e-6
        end

        page = export_visualization(joinpath(dir, "simulation_results"); textures=false, max_frames=40)
        html = read(page, String)
        start = findfirst("window.SPACEAGORA_VIEWER = ", html)
        stop = findnext(";\n</script>", html, last(start))
        payload = JSON.parse(html[last(start)+1:first(stop)-1])
        ap = payload["frames"]["arm_pose"]
        @test ap !== nothing
        @test ap["counts"] == [3] && ap["offsets"] == [0] && ap["total"] == 21
        data = _decode_f32(ap["data"])
        @test length(data) == 21 * payload["frames"]["count"]
        @test payload["scene"]["spacecraft"][1]["arm"]["links"][3]["radius_m"] == 0.04

        # A run without an arm has no arm block and no arm columns.
        dir2 = mktempdir()
        run_simulation(with_visualization_scene(_arm_config(results_directory=dir2, with_arm=false), true))
        df2 = DataFrame(Arrow.Table(joinpath(dir2, "simulation_results.feather")))
        @test !any(startswith("sc1_arm_pose"), names(df2))
    end
end
