using Test
using SpaceAGORA
using StaticArrays
using LinearAlgebra
using Arrow
using DataFrames
using JSON

import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft

const SM = SpaceAGORA.SimulationModel
const SV = SM.SceneVisualization

# Bus with two panels, the same shape the quickstart uses, but built directly
# so the test controls every pose.
function _three_link_model(; panel_offset_y::Float64=1.5)
    bus = SM.Link(root=true, m=100.0, dims=MVector{3, Float64}(1.0, 1.2, 1.4))
    left = SM.Link(root=false, m=5.0, dims=MVector{3, Float64}(0.02, 2.0, 1.0), r=MVector{3, Float64}(0.0, -panel_offset_y, 0.0))
    right = SM.Link(root=false, m=5.0, dims=MVector{3, Float64}(0.02, 2.0, 1.0), r=MVector{3, Float64}(0.0, panel_offset_y, 0.0))
    joint_l = SM.Joint(bus, SVector{3, Float64}(0.0, -0.6, 0.0), left, SVector{3, Float64}(0.0, 1.0, 0.0))
    joint_r = SM.Joint(bus, SVector{3, Float64}(0.0, 0.6, 0.0), right, SVector{3, Float64}(0.0, -1.0, 0.0))
    return SM.SpacecraftModel(links=[bus, left, right], joints=[joint_l, joint_r], root=bus, id=7)
end

_with_scene_flag(args::SM.SimulationConfiguration, flag::Bool) = with_visualization_scene(args, flag)

function _short_config(; results_directory::String, orientation_sim::Bool=false, flag::Bool=true)
    planet = make_no_gram_planet(:earth)
    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.05, 2.05, 2.8),
        panel_dims=(0.01, 5.7 / 2.0, 1.0),
        bus_mass=620.0,
        panel_mass_each=10.0,
        panel_offset_y=2.05 / 2.0 + 5.7 / 4.0,
        ic=SM.InitialCondition(ra=8_000.0e3, rp=7_000.0e3, i=45.0, ω=10.0, Ω=20.0, ν=0.0),
        prop_mass=0.0,
        id=1
    )
    args = make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=600.0,
        initial_time=SM.InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0),
        dynamic_effectors=(SM.InverseSquaredJ2GravityModel(),),
        density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(),
        orientation_sim=orientation_sim,
        keplerian=true,
        EI_km=120.0,
        verbose=false,
        results=true,
        results_directory=results_directory
    )
    return _with_scene_flag(args, flag)
end

struct ViewerReportedThruster <: SM.AbstractControlEffectorModel end
SM.ControlHooks.control_thruster_levels(::ViewerReportedThruster, i::Int) = i == 1 ? [0.25] : nothing

@testset "SceneVisualization" begin
    @testset "configuration preservation and opt-in fields" begin
        args = _short_config(results_directory=mktempdir(), flag=false)
        solver = SM.SolverConfig(solver_mode=:tsit5, maxiters=12345)
        paths = SM.FilePaths(results="viewer-field-preservation")
        args = SM.SimConfig._with_configuration(args; solver_config=solver, file_paths=paths)
        enabled = with_visualization_scene(args, true)
        @test !args.simulation_settings.save_visualization_scene
        @test enabled.simulation_settings.save_visualization_scene
        @test with_visualization_scene(enabled, true) === enabled
        for name in fieldnames(typeof(args))
            name === :simulation_settings && continue
            @test getfield(enabled, name) === getfield(args, name)
        end
        for name in fieldnames(typeof(args.simulation_settings))
            name === :save_visualization_scene && continue
            @test isequal(getfield(enabled.simulation_settings, name), getfield(args.simulation_settings, name))
        end
        moved = with_results_directory(enabled, "viewer-isolated-results")
        @test moved.solver_config === solver
        @test moved.file_paths === paths
        @test moved.simulation_settings.results_directory == "viewer-isolated-results"
        @test moved.simulation_settings.save_visualization_scene
        @test moved.environment_model === args.environment_model

        push!(args.dynamics_model.spacecraft[1].root.thrusters,
            SM.Thruster(max_thrust=2.0, location=MVector{3, Float64}(0.0, 0.0, 0.0),
                direction=MVector{3, Float64}(1.0, 0.0, 0.0)))
        reporting = SM.SimConfig._with_configuration(args;
            control_model=SM.ControlModel(control_effectors=(ViewerReportedThruster(),), control_rates=[1.0]))
        @test SM.ControlHooks.control_thruster_levels(nothing, 1) === nothing
        @test isempty(SM.SimulationCallbacks.visualization_save_fields(reporting))
        off_fields = SM.SimulationCallbacks.default_save_fields(reporting)
        @test :thruster_level ∉ getproperty.(off_fields, :name)
        @test :sun_dir ∉ getproperty.(off_fields, :name)
        reporting_on = with_visualization_scene(reporting, true)
        extras = SM.SimulationCallbacks.visualization_save_fields(reporting_on)
        @test :thruster_level in getproperty.(extras, :name)
        @test :sun_dir in getproperty.(extras, :name)
        @test SM.SimulationCallbacks.thruster_level_counts(reporting_on) == [1]
        integrator = (p=(args=reporting_on,),)
        @test SM.SimulationCallbacks._save_thruster_levels(1, [1], nothing, 0.0, integrator) == [[0.25]]
    end
    @testset "spacecraft geometry from link tree" begin
        model = _three_link_model()
        geometry = spacecraft_geometry(model)
        @test geometry.id == 7
        @test geometry.name == "sc7"
        @test length(geometry.links) == 3
        @test geometry.links[1].root
        @test geometry.links[1].r_m == SVector(0.0, 0.0, 0.0)
        @test geometry.links[1].q == SVector(0.0, 0.0, 0.0, 1.0)
        @test geometry.links[1].dims_m == SVector(1.0, 1.2, 1.4)
        @test !geometry.links[2].root
        @test geometry.links[2].r_m == SVector(0.0, -1.5, 0.0)
        @test geometry.links[3].r_m == SVector(0.0, 1.5, 0.0)
        @test geometry.links[2].mass_kg == 5.0
        # |r| + half the box diagonal of a panel: 1.5 + 0.5*sqrt(0.02^2 + 4 + 1)
        @test geometry.bounding_radius_m ≈ 1.5 + 0.5 * sqrt(0.02^2 + 2.0^2 + 1.0^2)
        @test geometry.stl_path === nothing
        @test length(geometry.joints) == 2
        @test geometry.joints[1].link1 == 1 && geometry.joints[1].link2 == 2
        @test geometry.joints[2].link1 == 1 && geometry.joints[2].link2 == 3
        @test geometry.joints[1].p2_m == SVector(0.0, 1.0, 0.0)
        @test isempty(geometry.thrusters)
        @test isempty(geometry.facets)

        named = spacecraft_geometry(model; name="probe", stl_path="cad/probe.stl")
        @test named.name == "probe"
        @test named.stl_path == "cad/probe.stl"

        single = SM.SpacecraftModel(root=SM.Link(root=true, m=3.0, dims=MVector{3, Float64}(0.5, 0.5, 0.1)), id=2)
        single_geometry = spacecraft_geometry(single)
        @test length(single_geometry.links) == 1
        @test single_geometry.bounding_radius_m ≈ 0.5 * sqrt(0.5^2 + 0.5^2 + 0.1^2)
    end

    @testset "thruster and facet glyphs carry their link index" begin
        model = _three_link_model()
        push!(model.root.thrusters, SM.Thruster(max_thrust=2.0, location=MVector{3, Float64}(0.5, 0.0, 0.0), direction=MVector{3, Float64}(-1.0, 0.0, 0.0)))
        push!(model.links[3].SRP_facets, SM.Facet(area=2.0, normal_vector=MVector{3, Float64}(1.0, 0.0, 0.0), name="panel_r"))
        geometry = spacecraft_geometry(model)
        @test length(geometry.thrusters) == 1
        @test geometry.thrusters[1].link == 1
        @test geometry.thrusters[1].max_thrust_n == 2.0
        @test geometry.thrusters[1].direction == SVector(-1.0, 0.0, 0.0)
        @test length(geometry.facets) == 1
        @test geometry.facets[1].link == 3
        @test geometry.facets[1].name == "panel_r"
        @test geometry.facets[1].area_m2 == 2.0
    end

    @testset "link pose snapshot follows rotate_link" begin
        model = _three_link_model()
        @test SV.link_pose_link_indices(model) == [2, 3]
        pose = SV.link_pose_vector(model)
        @test length(pose) == 14
        @test pose[1:3] == [0.0, -1.5, 0.0]
        @test pose[4:7] == [0.0, 0.0, 0.0, 1.0]
        @test pose[8:10] == [0.0, 1.5, 0.0]

        SM.rotate_link(model.links[2], SVector{3, Float64}(0.0, 1.0, 0.0), pi / 2)
        rotated = SV.link_pose_vector(model)
        @test rotated[1:3] == pose[1:3]
        @test rotated[4:7] ≈ [0.0, sin(pi / 4), 0.0, cos(pi / 4)]
        @test rotated[8:14] == pose[8:14]

        single = SM.SpacecraftModel(root=SM.Link(root=true), id=1)
        @test isempty(SV.link_pose_link_indices(single))
        @test isempty(SV.link_pose_vector(single))

        # The save-field getter reads the live Link objects through integrator.p.args.
        args = _short_config(results_directory=mktempdir())
        fake_integrator = (p=(args=args,),)
        poses = SM.SimulationCallbacks._save_link_poses(1, nothing, 0.0, fake_integrator)
        @test length(poses) == 1
        @test poses[1] == SV.link_pose_vector(args.dynamics_model.spacecraft[1])
        @test length(poses[1]) == 14
    end

    @testset "link_pose save field is opt-in" begin
        dir = mktempdir()
        off = _short_config(results_directory=dir, flag=false)
        on = _short_config(results_directory=dir, flag=true)
        names_off = Symbol[f.name for f in SM.default_save_fields(off)]
        names_on = Symbol[f.name for f in SM.default_save_fields(on)]
        @test :link_pose ∉ names_off
        @test :link_pose ∈ names_on
        @test filter(n -> n ∉ (:sun_dir, :link_pose), names_on) == names_off
        field = last(SM.default_save_fields(on))
        @test field.per_satellite
        @test field.column_prefix == "link_pose"

        # A single-body spacecraft has nothing to record even with the flag on.
        single = SM.SpacecraftModel(root=SM.Link(root=true, m=3.0), initial_condition=on.dynamics_model.spacecraft[1].initial_condition, id=1)
        single_args = SM.SimConfig._with_configuration(on;
        dynamics_model=SM.DynamicsModel([single], on.dynamics_model.dynamic_effectors),
    )
        @test :link_pose ∉ Symbol[f.name for f in SM.default_save_fields(single_args)]
    end

    @testset "rotation sample times" begin
        @test SV.rotation_sample_times(0.0, 10.0) == [0.0]
        t = SV.rotation_sample_times(100.0, 10.0)
        @test t[1] == 0.0 && t[end] == 100.0
        @test length(t) == 11
        coarse = SV.rotation_sample_times(10_000.0, 1.0; max_samples=101)
        @test length(coarse) == 101
        @test coarse[end] == 10_000.0
        @test all(diff(coarse) .≈ 100.0)
        odd = SV.rotation_sample_times(95.0, 10.0)
        @test odd[end] == 95.0 && odd[end-1] == 90.0
        @test_throws ArgumentError SV.rotation_sample_times(100.0, 10.0; max_samples=1)
        # A zero data rate falls back to the max_samples cadence instead of failing.
        @test length(SV.rotation_sample_times(100.0, 0.0; max_samples=11)) == 11
    end

    @testset "planet rotation table matches planet_frame_lpi" begin
        model = SM.SimpleEphemeridesModel()
        et_start = 4.5e8
        times = [0.0, 3600.0, 7200.0, 86_400.0 * 3]
        for planet in (make_no_gram_planet(:earth), make_no_gram_planet(:mars))
            table = planet_rotation_table(planet, model, et_start, times)
            @test length(table) == length(times)
            for (k, t) in enumerate(times)
                q = table[k]
                @test norm(q) ≈ 1.0
                dcm = SM.planet_frame_lpi(planet, et_start + t, model)
                @test SM.rot(q) ≈ dcm atol=1e-12
                # The table rotates an inertial position into the body-fixed frame
                # exactly as the runtime latitude/longitude columns do.
                r_i = SVector{3, Float64}(7.0e6, 1.0e6, -2.0e6)
                @test SM.rot(q) * r_i ≈ dcm * r_i
            end
            # sign continuity: neighbours are never more than 90 deg apart in quaternion space
            @test all(dot(table[k], table[k+1]) >= 0.0 for k in 1:length(table)-1)
        end

        spec = planet_spec(make_no_gram_planet(:mars), model, et_start, times)
        @test spec.name == "Mars"
        @test spec.texture == "mars"
        @test spec.equatorial_radius_m == 3.396190e6
        @test spec.polar_radius_m == 3.3762e6
        @test spec.spin_rad_s[3] ≈ 7.08823596e-5
        @test spec.inertial_frame == "J2000"
        @test spec.rotation_t_s == times
        @test length(spec.rotation_q_pi) == length(times)
    end

    @testset "velocity aligned attitude" begin
        r = SVector{3, Float64}(7.0e6, 0.0, 0.0)
        v = SVector{3, Float64}(0.0, 7.5e3, 1.0e3)
        q = velocity_aligned_quaternion(r, v)
        @test norm(q) ≈ 1.0
        R = SM.rot(q)
        @test R * R' ≈ I atol=1e-12
        @test det(R) ≈ 1.0
        @test R * (v / norm(v)) ≈ SVector(1.0, 0.0, 0.0) atol=1e-12
        nadir_b = R * (-r / norm(r))
        @test nadir_b[1] ≈ 0.0 atol=1e-12
        @test nadir_b[2] ≈ 0.0 atol=1e-12
        @test nadir_b[3] > 0.0
        # Circular-orbit velocity is already orthogonal to r: body z is exactly nadir.
        v_circ = SVector{3, Float64}(0.0, 7.5e3, 0.0)
        R_circ = SM.rot(velocity_aligned_quaternion(r, v_circ))
        @test R_circ * (-r / norm(r)) ≈ SVector(0.0, 0.0, 1.0) atol=1e-12
        # Radial flight still yields a proper rotation.
        R_radial = SM.rot(velocity_aligned_quaternion(r, SVector{3, Float64}(1.0e3, 0.0, 0.0)))
        @test R_radial * R_radial' ≈ I atol=1e-12
        @test R_radial * SVector(1.0, 0.0, 0.0) ≈ SVector(1.0, 0.0, 0.0) atol=1e-12
        @test_throws ArgumentError velocity_aligned_quaternion(r, zeros(3))
    end

    @testset "frame budget arithmetic" begin
        b = visualization_frame_budget(10, 1)
        @test b.frames == 10 && b.stride == 1 && b.bytes == 120
        b = visualization_frame_budget(10_000, 1; max_frames=2000)
        @test b.frames == 2000 && b.stride == 5
        # 32768 sats * 12 B = 393216 B/frame; 150 MB allows 381 frames
        b = visualization_frame_budget(100_000, 32_768)
        @test b.stride == cld(100_000, 381)
        @test b.frames == cld(100_000, b.stride)
        @test b.bytes <= 150.0e6
        @test visualization_frame_budget(0, 4).frames == 0
        @test visualization_frame_budget(1, 4).frames == 1
        @test visualization_frame_budget(2, 4; data_budget_mb=1e-9).frames == 2
        @test_throws ArgumentError visualization_frame_budget(10, 0)
        @test_throws ArgumentError visualization_frame_budget(10, 1; max_frames=1)
    end

    @testset "scene build, JSON round trip, sidecar path" begin
        dir = mktempdir()
        args = _short_config(results_directory=dir)
        scene = build_visualization_scene(args; rotation_max_samples=16)
        @test scene.schema == 1
        @test scene.epoch_utc == "2014-05-27T05:00:00.000Z"
        @test scene.epoch_et_start_s ≈ SM.ephemerides_time_seconds(args.initial_time, args.environment_model.ephemerides_model)
        @test scene.planet.name == "Earth"
        @test scene.planet.rotation_t_s[1] == 0.0
        @test scene.planet.rotation_t_s[end] == 600.0
        @test length(scene.planet.rotation_t_s) <= 16
        @test length(scene.spacecraft) == 1
        @test length(scene.spacecraft[1].links) == 3
        @test scene.orientation_sim == false
        @test scene.results_feather == "simulation_results.feather"
        @test scene.link_pose_field == "link_pose"
        @test scene.link_pose_stride == 7

        @test visualization_scene_path(args) == joinpath(dir, "simulation_results_scene.json")

        path = joinpath(dir, "nested", "scene.json")
        @test write_visualization_scene(path, scene) == path
        @test isfile(path)
        raw = JSON.parsefile(path)
        @test raw["schema"] == 1
        @test raw["results"]["link_pose"]["layout"] == ["rx", "ry", "rz", "qx", "qy", "qz", "qw"]
        @test raw["planet"]["rotation"]["q_pi"][1] isa Vector
        @test raw["spacecraft"][1]["stl_path"] === nothing
        @test read_visualization_scene(path) == scene

        bad = copy(raw)
        bad["schema"] = 99
        bad_path = joinpath(dir, "bad.json")
        open(io -> JSON.print(io, bad), bad_path, "w")
        @test_throws ArgumentError read_visualization_scene(bad_path)

        # The engine-facing writer honors the flag.
        off = _short_config(results_directory=dir, flag=false)
        @test SV.write_visualization_scene!(off) === nothing
        @test !isfile(visualization_scene_path(off))
        @test SV.write_visualization_scene!(args) == visualization_scene_path(args)
        @test read_visualization_scene(visualization_scene_path(args)) == build_visualization_scene(args)
    end

    @testset "end-to-end: flagged run writes sidecar and link_pose columns" begin
        dir = mktempdir()
        args = _short_config(results_directory=dir, flag=true)
        run_simulation(args; visualization=false)
        @test isfile(joinpath(dir, "simulation_results.csv"))
        feather = joinpath(dir, "simulation_results.feather")
        @test isfile(feather)
        sidecar = joinpath(dir, "simulation_results_scene.json")
        @test isfile(sidecar)
        scene = read_visualization_scene(sidecar)
        @test scene.planet.name == "Earth"
        @test scene.results_feather == basename(feather)
        df = DataFrame(Arrow.Table(feather))
        for k in 1:14
            @test "sc1_link_pose_$(k)" in names(df)
        end
        @test "sc1_link_pose_15" ∉ names(df)
        @test nrow(df) > 1
        # Panels are static in this run: every row carries the configured pose.
        expected = SV.link_pose_vector(args.dynamics_model.spacecraft[1])
        @test all(all(df[!, "sc1_link_pose_$(k)"] .== expected[k]) for k in 1:14)
        @test "sc1_pos_1" in names(df)

        # Default run: no sidecar, no link_pose columns.
        dir_off = mktempdir()
        run_simulation(_short_config(results_directory=dir_off, flag=false); visualization=false)
        @test !isfile(joinpath(dir_off, "simulation_results_scene.json"))
        df_off = DataFrame(Arrow.Table(joinpath(dir_off, "simulation_results.feather")))
        @test !any(startswith("sc1_link_pose"), names(df_off))
        @test !any(startswith("sun_dir"), names(df_off))
        @test !any(occursin("thruster_level"), names(df_off))
        @test "sun_dir_1" in names(df)
        # Saving viewer data must not alter any existing result column.
        @test all(isequal(df[!, name], df_off[!, name]) for name in names(df_off))
    end
end
