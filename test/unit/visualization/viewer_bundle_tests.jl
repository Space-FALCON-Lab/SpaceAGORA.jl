using Test
using SpaceAGORA
using StaticArrays
using DataFrames
using Arrow
using JSON
using Base64
using Random

import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft

const SM = SpaceAGORA.SimulationModel
const SV = SM.SceneVisualization
const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))

function _decode_f32(b64::AbstractString)
    bytes = base64decode(b64)
    return collect(reinterpret(Float32, bytes))
end

function _decode_f64(b64::AbstractString)
    bytes = base64decode(b64)
    return collect(reinterpret(Float64, bytes))
end

# One-triangle binary STL, enough to exercise the override plumbing.
function _write_tiny_stl(path::AbstractString)
    open(path, "w") do io
        write(io, zeros(UInt8, 80))
        write(io, UInt32(1))
        for v in (0.0f0, 0.0f0, 1.0f0,  0.0f0, 0.0f0, 0.0f0,  1.0f0, 0.0f0, 0.0f0,  0.0f0, 1.0f0, 0.0f0)
            write(io, v)
        end
        write(io, UInt16(0))
    end
    return path
end

# JPEG dimensions from the first SOF marker; enough to assert the 2:1 aspect.
function _jpeg_size(path::AbstractString)
    bytes = read(path)
    (bytes[1] == 0xFF && bytes[2] == 0xD8) || return nothing
    i = 3
    while i + 9 <= length(bytes)
        bytes[i] == 0xFF || return nothing
        marker = bytes[i + 1]
        if marker in (0xC0, 0xC1, 0xC2)
            h = Int(bytes[i + 5]) << 8 | Int(bytes[i + 6])
            w = Int(bytes[i + 7]) << 8 | Int(bytes[i + 8])
            return (w, h)
        end
        seg = Int(bytes[i + 2]) << 8 | Int(bytes[i + 3])
        i += 2 + seg
    end
    return nothing
end

function _viewer_config(; results_directory::String, mission_time::Float64=600.0, orientation_sim::Bool=false)
    planet = make_no_gram_planet(:mars)
    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.0, 2.0, 2.5),
        panel_dims=(0.01, 2.5, 1.0),
        bus_mass=500.0,
        panel_mass_each=10.0,
        panel_offset_y=2.3,
        ic=SM.InitialCondition(ra=4_500.0e3, rp=3_800.0e3, i=30.0, ω=0.0, Ω=45.0, ν=0.0),
        prop_mass=0.0,
        id=1
    )
    args = make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=mission_time,
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
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
    return args
end

# A synthetic results table with the same column layout the run writer produces.
function _synthetic_results(scene::SV.VisualizationScene; n_rows::Int=25, with_vel::Bool=true, with_q::Bool=false, with_link_pose::Bool=true)
    S = length(scene.spacecraft)
    df = DataFrame(time=collect(0.0:10.0:(10.0 * (n_rows - 1))))
    for i in 1:S
        for c in 1:3
            df[!, "sc$(i)_pos_$(c)"] = [1.0e6 * (i + c) + 1000.0 * r for r in 1:n_rows]
            with_vel && (df[!, "sc$(i)_vel_$(c)"] = [7000.0 * c + r for r in 1:n_rows])
        end
        if with_q
            for c in 1:4
                df[!, "sc$(i)_q_$(c)"] = fill(c == 4 ? 1.0 : 0.0, n_rows)
            end
        end
        df[!, "sc$(i)_mass"] = [500.0 - r for r in 1:n_rows]
        n_links = length(scene.spacecraft[i].links) - 1
        if with_link_pose && n_links > 0
            for k in 1:(7 * n_links)
                df[!, "sc$(i)_link_pose_$(k)"] = [0.1 * k + 0.001 * r for r in 1:n_rows]
            end
        end
    end
    return df
end

@testset "ViewerBundle" begin
    @testset "kept rows" begin
        @test SV.kept_row_indices(0, 1) == Int[]
        @test SV.kept_row_indices(1, 3) == [1]
        @test SV.kept_row_indices(10, 1) == collect(1:10)
        @test SV.kept_row_indices(10, 3) == [1, 4, 7, 10]
        @test SV.kept_row_indices(11, 3) == [1, 4, 7, 10, 11]
    end

    @testset "texture manifest covers the five launch bodies" begin
        manifest = SV.texture_manifest()
        for body in ("earth", "mars", "venus", "titan", "moon")
            @test haskey(manifest, body)
            for (tier, entry) in manifest[body]
                @test isfile(entry["path"])
                size = _jpeg_size(entry["path"])
                @test size !== nothing
                @test size[1] == 2 * size[2]
                @test size[1] == entry["width"] && size[2] == entry["height"]
                @test size[1] == SV._tier_pixels(tier)
                @test -180.0 <= entry["lon_left_deg"] <= 180.0
                @test occursin("public-domain", entry["license"])
            end
            @test haskey(manifest[body], "4k")
            @test filesize(manifest[body]["4k"]["path"] ) < 4_000_000
        end
        @test manifest["titan"]["4k"]["lon_left_deg"] == 0.0
        @test manifest["earth"]["4k"]["lon_left_deg"] == -180.0
        @test SV.texture_entry("earth"; resolution="4k")["resolution"] == "4k"
        @test SV.texture_entry("earth"; resolution=:best)["resolution"] == "4k"
        @test SV.texture_entry("earth"; resolution="16k")["resolution"] == "4k"   # largest at or below
        @test SV.texture_entry("venus"; resolution="8k")["resolution"] == "4k"     # only tier
        @test SV.texture_entry("titan"; resolution="2k")["resolution"] == "4k"     # smallest available
        @test_throws ArgumentError SV.texture_entry("earth"; resolution="huge")

        payload = SV.texture_payload("Mars"; resolution="4k")
        @test payload !== nothing
        @test startswith(payload["url"], "data:image/jpeg;base64,")
        @test payload["lon_left_deg"] == -180.0
        @test payload["resolution"] == "4k" && payload["width"] == 4096 && payload["height"] == 2048
        @test SV.texture_payload("Mars")["resolution"] == "4k"
        @test SV.texture_payload("pluto") === nothing
        @test isempty(SV.texture_manifest(mktempdir()))
    end

    @testset "optional texture tiers and missing files" begin
        dir = mktempdir()
        # The selector only reads metadata and file presence; these small
        # sentinels exercise optional tiers without shipping large textures.
        write(joinpath(dir, "4k.jpg"), UInt8[0xff, 0xd8, 0xff, 0xd9])
        write(joinpath(dir, "8k.jpg"), UInt8[0xff, 0xd8, 0xff, 0xd9])
        write(joinpath(dir, "manifest.toml"), """
        [[texture]]
        body = "test"
        file = "4k.jpg"
        resolution = "4k"
        [[texture]]
        body = "test"
        file = "8k.jpg"
        resolution = "8k"
        [[texture]]
        body = "test"
        file = "missing16k.jpg"
        resolution = "16k"
        """)
        @test SV.texture_entry("test"; dir=dir, resolution=:best)["resolution"] == "8k"
        @test SV.texture_entry("test"; dir=dir, resolution="16k")["resolution"] == "8k"
        @test SV.texture_entry("test"; dir=dir, resolution="4k")["resolution"] == "4k"
        @test SV.texture_entry("test"; dir=dir, resolution="2k")["resolution"] == "4k"
        rm(joinpath(dir, "8k.jpg"))
        @test SV.texture_entry("test"; dir=dir, resolution=:best)["resolution"] == "4k"
        rm(joinpath(dir, "4k.jpg"))
        @test SV.texture_entry("test"; dir=dir, resolution=:best) === nothing
    end

    @testset "frame payload encodes and decimates" begin
        dir = mktempdir()
        args = _viewer_config(results_directory=dir)
        scene = build_visualization_scene(args; rotation_max_samples=8)
        df = _synthetic_results(scene; n_rows=25)

        frames = SV.build_viewer_frames(df, scene)
        @test frames["count"] == 25 && frames["sats"] == 1 && frames["stride_rows"] == 1
        @test frames["source_rows"] == 25
        @test frames["t_dtype"] == "f64"
        t = _decode_f64(frames["t_s"])
        @test length(t) == 25 && t[end] == 240.0
        # A single spacecraft gets Float64 positions (exact to the meter at 7000 km).
        @test frames["pos_dtype"] == "f64"
        pos = _decode_f64(frames["pos_km"])
        @test length(pos) == 25 * 3
        @test pos[1] == (1.0e6 * 2 + 1000.0) / 1000.0
        @test pos[3 * 24 + 3] == (1.0e6 * 4 + 25_000.0) / 1000.0
        mass = _decode_f32(frames["mass_kg"])
        @test length(mass) == 25 && mass[1] == 499.0f0 && mass[end] == 475.0f0
        vel = _decode_f32(frames["vel_kms"])
        @test vel[1] ≈ Float32(7001.0 / 1000.0)
        @test frames["q"] === nothing
        lp = frames["link_pose"]
        @test lp["stride"] == 7 && lp["counts"] == [2] && lp["offsets"] == [0] && lp["total"] == 14
        lpd = _decode_f32(lp["data"])
        @test length(lpd) == 25 * 14
        @test lpd[14 + 3] ≈ Float32(0.3 + 0.002)

        # Decimation keeps the first and last rows.
        small = SV.build_viewer_frames(df, scene; max_frames=6)
        @test small["stride_rows"] == 5
        # rows 1, 6, 11, 16, 21 plus the final row 25
        @test small["count"] == 6
        ts = _decode_f64(small["t_s"])
        @test ts[1] == 0.0 && ts[end] == 240.0
        # A tiny byte budget wins over max_frames.
        tiny = SV.build_viewer_frames(df, scene; max_frames=2000, data_budget_mb=1e-4)
        @test tiny["count"] < 25 && tiny["count"] >= 2

        # Quaternions are included when every spacecraft has them; link poses
        # drop out when their columns are missing.
        df_q = _synthetic_results(scene; with_q=true, with_link_pose=false)
        fq = SV.build_viewer_frames(df_q, scene)
        @test length(_decode_f32(fq["q"])) == 25 * 4
        @test fq["link_pose"] === nothing
        @test_throws ArgumentError SV.build_viewer_frames(DataFrame(time=[0.0]), scene)

        # Above the Float64 cutoff the positions fall back to Float32.
        many = SM.SpacecraftModel[SM.SpacecraftModel(root=SM.Link(root=true, m=3.0), initial_condition=args.dynamics_model.spacecraft[1].initial_condition, id=i) for i in 1:(SV.FLOAT64_POSITION_MAX_SPACECRAFT + 1)]
        many_args = SM.SimConfig._with_configuration(args;
        dynamics_model=SM.DynamicsModel(many, args.dynamics_model.dynamic_effectors),
    )
        many_scene = build_visualization_scene(many_args; rotation_max_samples=4)
        many_df = _synthetic_results(many_scene; n_rows=4, with_vel=false, with_link_pose=false)
        many_frames = SV.build_viewer_frames(many_df, many_scene)
        @test many_frames["pos_dtype"] == "f32"
        @test length(_decode_f32(many_frames["pos_km"])) == 4 * length(many) * 3
        @test many_frames["vel_kms"] === nothing
    end

    @testset "plume payload" begin
        dir = mktempdir()
        args = _viewer_config(results_directory=dir)
        scene = build_visualization_scene(args; rotation_max_samples=8)
        df = _synthetic_results(scene; n_rows=25)

        # Without the sc{i}_plume_* columns the block is absent and the rest of
        # the payload is unchanged.
        @test SV.build_viewer_frames(df, scene)["plume"] === nothing

        for (k, f) in enumerate(SV.PLUME_FRAME_FIELDS)
            df[!, "sc1_plume_$(f)"] = [10.0 * k + r for r in 1:25]
        end
        frames = SV.build_viewer_frames(df, scene)
        plume = frames["plume"]
        @test plume !== nothing
        @test sort(collect(keys(plume))) == sort(collect(String.(SV.PLUME_FRAME_FIELDS)))
        for (k, f) in enumerate(SV.PLUME_FRAME_FIELDS)
            values = _decode_f32(plume[f])
            @test length(values) == 25                      # one spacecraft, frame-major
            @test values[1] ≈ Float32(10.0 * k + 1)
            @test values[end] ≈ Float32(10.0 * k + 25)
        end
        # Decimation applies to the block like every other frame array.
        small = SV.build_viewer_frames(df, scene; max_frames=6)
        @test length(_decode_f32(small["plume"]["height_m"])) == small["count"]

        # One spacecraft short of the columns drops the whole block.
        partial = select(df, Not("sc1_plume_eroded_kg"))
        @test SV.build_viewer_frames(partial, scene)["plume"] === nothing
    end

    @testset "thruster levels" begin
        dir = mktempdir()
        args = _viewer_config(results_directory=dir)
        model = args.dynamics_model.spacecraft[1]
        # straight onto the link, the way Link(thrusters=...) carries them
        push!(model.root.thrusters, SM.Thruster(max_thrust=400.0, location=MVector{3, Float64}(0.0, 0.0, 1.0), direction=MVector{3, Float64}(0.0, 0.0, 1.0)))
        push!(model.root.thrusters, SM.Thruster(max_thrust=20.0, location=MVector{3, Float64}(1.0, 0.0, 0.0), direction=MVector{3, Float64}(1.0, 0.0, 0.0)))
        scene = build_visualization_scene(args; rotation_max_samples=4)
        @test length(scene.spacecraft[1].thrusters) == 2
        df = _synthetic_results(scene; n_rows=25)

        # No `sc{i}_thruster_level_{k}` columns: the block is absent and the viewer runs without it.
        bare = SV.build_viewer_frames(df, scene)
        @test bare["thruster_level"] === nothing
        @test bare["thruster_counts"] === nothing

        df[!, "sc1_thruster_level_1"] = [r / 25 for r in 1:25]
        df[!, "sc1_thruster_level_2"] = [1.0 - r / 25 for r in 1:25]
        frames = SV.build_viewer_frames(df, scene)
        @test frames["thruster_counts"] == [2]
        levels = _decode_f32(frames["thruster_level"])
        # frame-major, then spacecraft, then that spacecraft's thrusters in scene order
        @test length(levels) == 25 * 2
        @test levels[1] ≈ Float32(1 / 25)
        @test levels[2] ≈ Float32(1 - 1 / 25)
        @test levels[end - 1] ≈ 1.0f0
        @test levels[end] ≈ 0.0f0

        # Values outside 0..1 are clipped on the way into the payload.
        df[!, "sc1_thruster_level_1"] = fill(1.7, 25)
        df[!, "sc1_thruster_level_2"] = fill(-0.3, 25)
        clipped = _decode_f32(SV.build_viewer_frames(df, scene)["thruster_level"])
        @test all(==(1.0f0), clipped[1:2:end])
        @test all(==(0.0f0), clipped[2:2:end])

        # Decimation keeps the block aligned with the kept frames.
        small = SV.build_viewer_frames(df, scene; max_frames=6)
        @test length(_decode_f32(small["thruster_level"])) == small["count"] * 2

        # A spacecraft the scene draws thrusters for but the run saved no levels for stays at zero count.
        partial = select(df, Not(["sc1_thruster_level_2"]))
        @test SV.build_viewer_frames(partial, scene)["thruster_level"] === nothing
    end

    @testset "sun direction payload" begin
        dir = mktempdir()
        args = _viewer_config(results_directory=dir)
        scene = build_visualization_scene(args; rotation_max_samples=8)
        df = _synthetic_results(scene; n_rows=25)

        # Without the columns the block is simply absent and the viewer falls
        # back to its fixed light.
        @test SV.build_viewer_frames(df, scene)["sun_dir"] === nothing

        # `sun_dir_1..3` is one direction per row, shared by the whole scene.
        for c in 1:3
            df[!, "sun_dir_$(c)"] = [c == 1 ? 1.0 : 0.0 for _ in 1:25]
        end
        df[!, "sun_dir_1"] = [r <= 13 ? 1.0 : 0.0 for r in 1:25]
        df[!, "sun_dir_2"] = [r <= 13 ? 0.0 : 1.0 for r in 1:25]
        frames = SV.build_viewer_frames(df, scene)
        sun = _decode_f32(frames["sun_dir"])
        @test length(sun) == 25 * 3
        @test sun[1] == 1.0f0 && sun[2] == 0.0f0 && sun[3] == 0.0f0
        @test sun[3 * 24 + 1] == 0.0f0 && sun[3 * 24 + 2] == 1.0f0

        # Decimation keeps the block aligned with the rows it kept.
        small = SV.build_viewer_frames(df, scene; max_frames=6)
        small_sun = _decode_f32(small["sun_dir"])
        @test length(small_sun) == 3 * small["count"]
        @test small_sun[1] == 1.0f0
        @test small_sun[3 * (small["count"] - 1) + 2] == 1.0f0
    end

    @testset "sun direction save field" begin
        # The simple ephemerides model at Earth resolves the Sun, so the run
        # writes the columns; at Mars it cannot and they are left out.
        base = _viewer_config(results_directory=mktempdir())
        env = base.environment_model
        earth_env = SM.EnvironmentModel(make_no_gram_planet(:earth), env.EI, env.density_model,
            SM.SimpleEphemeridesModel(), env.topography, env.topo_degree, env.topo_order, env.wind, env.thermal_model)
        earth_args = SM.SimConfig._with_configuration(base;
        environment_model=earth_env,
    )
        names_of(a) = Symbol[f.name for f in SM.SimulationCallbacks.default_save_fields(a)]
        @test :sun_dir in names_of(with_visualization_scene(earth_args, true))
        @test :sun_dir ∉ names_of(with_visualization_scene(earth_args, false))

        mars_args = _viewer_config(results_directory=mktempdir())
        @test mars_args.environment_model.planet.name == "Mars"
        @test :sun_dir ∉ names_of(mars_args)
    end

    @testset "model formats" begin
        @test SV.model_format("bus.stl") == ("stl", "model/stl")
        @test SV.model_format("BUS.OBJ") == ("obj", "model/obj")
        @test SV.model_format("iss.glb") == ("glb", "model/gltf-binary")
        @test SV.model_format("a.gltf") == ("gltf", "model/gltf+json")
        @test_throws ArgumentError SV.model_format("mesh.fbx")

        dir = mktempdir()
        args = _viewer_config(results_directory=dir)
        scene = build_visualization_scene(args; rotation_max_samples=4)
        obj = joinpath(dir, "bus.obj")
        write(obj, "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n")
        fake_glb = joinpath(dir, "bus.glb")
        write(fake_glb, vcat(Vector{UInt8}("glTF"), zeros(UInt8, 20)))
        models = SV.model_payloads(scene; models=Dict(1 => obj), model_scale=Dict(1 => 0.5), model_rotation_deg=Dict(1 => (-90.0, 0.0, -90.0)))
        @test models["1"]["format"] == "obj"
        @test startswith(models["1"]["url"], "data:model/obj;base64,")
        @test models["1"]["scale"] == 0.5
        @test models["1"]["rotation_deg"] == [-90.0, 0.0, -90.0]
        glb = SV.model_payloads(scene; models=Dict(1 => fake_glb), model_scale=3.0, model_center=false)   # not a real glTF: no centeing
        @test glb["1"]["format"] == "glb" && glb["1"]["scale"] == 3.0 && glb["1"]["rotation_deg"] == [0.0, 0.0, 0.0]
        @test startswith(glb["1"]["url"], "data:model/gltf-binary;base64,")
        @test_throws ArgumentError SV.model_payloads(scene; models=Dict(1 => obj), model_rotation_deg=Dict(1 => (1.0, 2.0)))
        # A glTF that needs a Draco decoder is refused with guidance, not embedded to fail silently.
        draco_json = "{\"asset\":{\"version\":\"2.0\"},\"extensionsRequired\":[\"KHR_draco_mesh_compression\"]}"
        draco_gltf = joinpath(dir, "packed.gltf")
        write(draco_gltf, draco_json)
        @test SV.gltf_required_extensions(draco_gltf) == ["KHR_draco_mesh_compression"]
        err = try; SV.model_payloads(scene; models=Dict(1 => draco_gltf)); nothing; catch e; e; end
        @test err isa ArgumentError && occursin("gltf-transform", err.msg)
        draco_glb = joinpath(dir, "packed.glb")
        open(draco_glb, "w") do io
            payload = Vector{UInt8}(draco_json)
            while length(payload) % 4 != 0; push!(payload, UInt8(' ')); end
            write(io, "glTF"); write(io, UInt32(2)); write(io, UInt32(12 + 8 + length(payload)))
            write(io, UInt32(length(payload))); write(io, "JSON"); write(io, payload)
        end
        @test SV.gltf_required_extensions(draco_glb) == ["KHR_draco_mesh_compression"]
        @test_throws ArgumentError SV.model_payloads(scene; models=Dict(1 => draco_glb))
        @test SV.gltf_required_extensions(fake_glb) == String[]
        @test_throws ArgumentError SV.model_payloads(scene; models=Dict(1 => joinpath(dir, "x.fbx")))

        # The shipped NASA ISS model embeds as GLB.
        iss = joinpath(REPO, "data", "models", "iss_nasa_3d_resources_b.glb")
        @test isfile(iss)
        @test read(iss, 4) == Vector{UInt8}("glTF")
        @test SV.gltf_required_extensions(iss) == String[]   # shipped decompressed; the NASA original needs Draco
    end

    @testset "model geometry readers and point clouds" begin
        dir = mktempdir()
        obj = joinpath(dir, "quad.obj")
        write(obj, "v 0 0 0\nv 2 0 0\nv 2 1 0\nv 0 1 0\nf 1 2 3 4\n")
        tris = load_model_triangles(obj)
        @test size(tris) == (3, 6)          # one quad, fan-triangulated
        lo, hi = model_bounding_box(obj)
        @test lo == SVector(0.0, 0.0, 0.0) && hi == SVector(2.0, 1.0, 0.0)
        @test SV.model_bounding_box_center(obj) == SVector(1.0, 0.5, 0.0)
        # scale then XYZ Euler rotation: (-90, 0, -90) sends model x -> body z, y -> x, z -> y
        rotated = load_model_triangles(obj; scale=2.0, rotation_deg=(-90, 0, -90))
        @test maximum(rotated[3, :]) ≈ 4.0 atol=1e-12
        @test maximum(rotated[1, :]) ≈ 2.0 atol=1e-12
        @test all(abs.(rotated[2, :]) .< 1e-12)

        stl = _write_tiny_stl(joinpath(dir, "tri.stl"))
        @test size(load_model_triangles(stl)) == (3, 3)
        @test model_bounding_box(stl)[2] == SVector(1.0, 1.0, 0.0)

        iss = joinpath(REPO, "data", "models", "iss_nasa_3d_resources_b.glb")
        iss_tris = load_model_triangles(iss)
        @test size(iss_tris, 2) > 100_000 && size(iss_tris, 2) % 3 == 0
        lo, hi = model_bounding_box(iss)
        @test lo ≈ SVector(-11.13, -4.17, -22.86) atol=0.01
        @test hi ≈ SVector(3.29, 41.26, 22.68) atol=0.01
        cloud = sample_model_pointcloud(iss; n_points=2000, rng=MersenneTwister(3), scale=2.4, rotation_deg=(-90, 0, -90))
        @test size(cloud) == (3, 2000)
        # centerd and rotated: the ~46-unit truss now spans about 109 m along body y, centerd on zero
        @test maximum(cloud[2, :]) - minimum(cloud[2, :]) > 100.0
        @test abs(maximum(cloud[2, :]) + minimum(cloud[2, :])) < 6.0
        @test maximum(abs.(cloud[3, :])) < 20.0
        raw = sample_model_pointcloud(iss; n_points=500, rng=MersenneTwister(3), center=false)
        @test minimum(raw[2, :]) > -5.0 && maximum(raw[2, :]) > 30.0   # model units, uncentred

        # A glTF whose buffers live in an external file is refused.
        write(joinpath(dir, "ext.gltf"), "{\"asset\":{\"version\":\"2.0\"},\"buffers\":[{\"uri\":\"data.bin\",\"byteLength\":4}]}")
        @test_throws ArgumentError load_model_triangles(joinpath(dir, "ext.gltf"))

        # The payload carries the center the viewer subtracts, unless centeing is off.
        args = _viewer_config(results_directory=dir)
        scene = build_visualization_scene(args; rotation_max_samples=4)
        m = SV.model_payloads(scene; models=Dict(1 => iss), model_scale=2.4)
        @test m["1"]["center"] ≈ [-3.92, 18.55, -0.09] atol=0.01
        @test SV.model_payloads(scene; models=Dict(1 => iss), model_center=false)["1"]["center"] == [0.0, 0.0, 0.0]
        @test SV.model_payloads(scene; models=Dict(1 => iss), model_center=Dict(1 => false))["1"]["center"] == [0.0, 0.0, 0.0]
    end

    @testset "reference paths" begin
        pts = [0.0 100.0 200.0; 0.0 10.0 0.0; 0.0 0.0 5.0]
        out = SV.path_payloads([(name="plan", points_m=pts, frame=:rtn, target=2, color="#ff0000", dashed=false)])
        @test length(out) == 1 && out[1]["name"] == "plan" && out[1]["frame"] == "rtn" && out[1]["target"] == 2
        @test out[1]["count"] == 3 && !out[1]["dashed"] && out[1]["color"] == "#ff0000"
        @test _decode_f32(out[1]["points_km"]) ≈ Float32.(vec(pts) ./ 1000.0)
        d = SV.path_payloads([Dict("points_m" => pts)])
        @test d[1]["frame"] == "inertial" && d[1]["dashed"] && d[1]["name"] == "path 1"
        @test_throws ArgumentError SV.path_payloads([(points_m=pts, frame=:lvlh)])
        @test_throws ArgumentError SV.path_payloads([(points_m=pts[1:2, :],)])
        @test SV.path_payloads(()) == Dict{String, Any}[]
    end

    @testset "reference ghosts" begin
        dir = mktempdir()
        args = with_visualization_scene(_viewer_config(results_directory=dir), true)
        scene = build_visualization_scene(args; rotation_max_samples=4)
        t = [0.0, 10.0, 20.0, 30.0]
        pos = [3.9e6 3.91e6 3.92e6 3.93e6; 0.0 1.0e4 2.0e4 3.0e4; 100.0 200.0 300.0 400.0]
        vel = [0.0 0.0 0.0 0.0; 3.4e3 3.4e3 3.4e3 3.4e3; 0.0 0.0 0.0 0.0]
        q = repeat([0.0, 0.0, 0.0, 1.0], 1, 4)
        out = SV.reference_payloads([(name="SPICE", t_s=t, pos_m=pos, vel_mps=vel, q=q, target=1, color="#123456", opacity=0.3, trail=false)], scene)
        @test length(out) == 1
        r = out[1]
        @test r["name"] == "SPICE" && r["target"] == 1 && r["count"] == 4 && r["color"] == "#123456" && r["opacity"] == 0.3 && !r["trail"]
        @test _decode_f64(r["t_s"]) == t
        @test _decode_f64(r["pos_km"]) ≈ vec(pos) ./ 1000.0
        @test _decode_f32(r["vel_kms"]) ≈ Float32.(vec(vel) ./ 1000.0)
        @test _decode_f32(r["q"]) == Float32.(vec(q))
        # Defaults, dictionary spelling, and optional blocks left out.
        d = SV.reference_payloads([Dict("t_s" => t, "pos_m" => pos)], scene)
        @test d[1]["name"] == "reference 1" && d[1]["target"] == 1 && d[1]["vel_kms"] === nothing && d[1]["q"] === nothing
        @test d[1]["opacity"] == 0.45 && d[1]["trail"] && d[1]["color"] == "#ff8c69"
        @test SV.reference_payloads((), scene) == Dict{String, Any}[]
        # Validation: shapes must agree with t_s, times ordered, target in range.
        @test_throws ArgumentError SV.reference_payloads([(pos_m=pos,)], scene)
        @test_throws ArgumentError SV.reference_payloads([(t_s=t,)], scene)
        @test_throws ArgumentError SV.reference_payloads([(t_s=t, pos_m=pos[:, 1:3])], scene)
        @test_throws ArgumentError SV.reference_payloads([(t_s=reverse(t), pos_m=pos)], scene)
        @test_throws ArgumentError SV.reference_payloads([(t_s=t, pos_m=pos, vel_mps=vel[:, 1:2])], scene)
        @test_throws ArgumentError SV.reference_payloads([(t_s=t, pos_m=pos, q=q[1:3, :])], scene)
        @test_throws ArgumentError SV.reference_payloads([(t_s=t, pos_m=pos, target=2)], scene)
        @test_throws ArgumentError SV.reference_payloads([(t_s=t, pos_m=pos, opacity=1.5)], scene)
        @test_throws ArgumentError SV.reference_payloads([(t_s=Float64[], pos_m=zeros(3, 0))], scene)
        # The payload carries the block and the page ships the module that draws it.
        df = _synthetic_results(scene)
        payload = SV.viewer_payload(scene, df; include_textures=false, references=[(t_s=t, pos_m=pos)])
        @test length(payload["references"]) == 1
        @test haskey(SV.viewer_import_map()["imports"], "viewer/references.js")
        @test haskey(SV.viewer_import_map()["imports"], "viewer/plots.js")
        @test haskey(SV.viewer_import_map()["imports"], "viewer/terrain.js")
    end

    @testset "STL overrides" begin
        dir = mktempdir()
        args = _viewer_config(results_directory=dir)
        stl = _write_tiny_stl(joinpath(dir, "bus.stl"))
        scene = build_visualization_scene(args; rotation_max_samples=4, stl_paths=Dict(1 => stl))
        @test scene.spacecraft[1].stl_path == stl
        models = SV.model_payloads(scene)
        @test haskey(models, "1")
        @test startswith(models["1"]["url"], "data:model/stl;base64,")
        @test models["1"]["scale"] == 1.0
        @test models["1"]["source"] == "bus.stl"
        @test length(base64decode(split(models["1"]["url"], ",")[2])) == 84 + 50

        plain = build_visualization_scene(args; rotation_max_samples=4)
        @test isempty(SV.model_payloads(plain))
        override = SV.model_payloads(plain; stl=Dict(1 => stl), stl_scale=0.001)
        @test haskey(override, "1") && override["1"]["scale"] == 0.001
        @test_throws ArgumentError SV.model_payloads(plain; stl=Dict(1 => joinpath(dir, "missing.stl")))

        df = _synthetic_results(scene; n_rows=4)
        payload = SV.viewer_payload(scene, df; include_textures=false)
        @test haskey(payload["models"], "1")
        # The sidecar keeps the path; the page keeps the bytes.
        raw = SV.scene_dict(scene)
        @test raw["spacecraft"][1]["stl_path"] == stl
    end

    @testset "import map and HTML assembly" begin
        imports = SV.viewer_import_map()["imports"]
        for key in ("three", "three/addons/controls/OrbitControls.js", "three/addons/loaders/STLLoader.js", "three/addons/loaders/OBJLoader.js", "three/addons/loaders/GLTFLoader.js", "three/addons/utils/BufferGeometryUtils.js", "viewer/main.js", "viewer/globe.js", "viewer/data.js", "viewer/spacecraft.js", "viewer/lod.js", "viewer/timeline.js", "viewer/ui.js")
            @test haskey(imports, key)
            @test startswith(imports[key], "data:text/javascript;base64,")
        end
        @test_throws ArgumentError SV.viewer_import_map(mktempdir())

        dir = mktempdir()
        args = _viewer_config(results_directory=dir)
        scene = build_visualization_scene(args; rotation_max_samples=8)
        df = _synthetic_results(scene)
        payload = SV.viewer_payload(scene, df; texture_resolution="4k", options=Dict("frame" => "planet_fixed", "title" => "a</script>b"))
        @test haskey(payload["textures"], "mars")
        @test payload["textures"]["mars"]["resolution"] == "4k"
        html = SV.render_viewer_html(payload; title="T<1>")
        @test occursin("<title>T&lt;1&gt;</title>", html)
        @test occursin("window.SPACEAGORA_VIEWER = {", html)
        @test occursin("\"planet_fixed\"", html)
        # The only closing script tags are the template's own.
        @test count("</script>", html) == 3
        @test occursin("a<\\/script>b", html)
        @test occursin("<script type=\"importmap\">{", html)
        @test !occursin("__PAYLOAD__", html) && !occursin("__IMPORTMAP__", html) && !occursin("__TITLE__", html)

        no_tex = SV.viewer_payload(scene, df; include_textures=false)
        @test isempty(no_tex["textures"])
        @test_throws ArgumentError SV.render_viewer_html(payload; viewer_dir=mktempdir())
    end

    @testset "end-to-end export from a flagged run" begin
        dir = mktempdir()
        args = with_visualization_scene(_viewer_config(results_directory=dir), true)
        @test args.simulation_settings.save_visualization_scene
        @test with_visualization_scene(args, true) === args
        run_simulation(args)
        prefix = joinpath(dir, "simulation_results")
        out = export_visualization(prefix; max_frames=50, title="unit run", texture_resolution="4k")
        @test out == prefix * "_viewer.html"
        @test isfile(out)
        html = read(out, String)
        @test occursin("<title>unit run</title>", html)
        @test occursin("data:image/jpeg;base64,", html)
        @test filesize(out) > 1_000_000
        # Payload round trip through the page text.
        start = findfirst("window.SPACEAGORA_VIEWER = ", html)
        stop = findnext(";\n</script>", html, last(start))
        payload = JSON.parse(html[last(start)+1:first(stop)-1])
        @test payload["scene"]["planet"]["name"] == "Mars"
        @test payload["frames"]["count"] <= 50
        @test payload["frames"]["sats"] == 1
        @test payload["frames"]["link_pose"]["counts"] == [2]
        @test payload["options"]["frame"] == "inertial"

        iss = joinpath(REPO, "data", "models", "iss_nasa_3d_resources_b.glb")
        iss_page = export_visualization(args; out=joinpath(dir, "iss.html"), textures=false, models=Dict(1 => iss), model_scale=2.4, model_rotation_deg=Dict(1 => (-90, 0, -90)))
        iss_html = read(iss_page, String)
        @test occursin("data:model/gltf-binary;base64,", iss_html)
        @test occursin("\"format\":\"glb\"", iss_html)
        @test occursin("iss_nasa_3d_resources_b.glb", iss_html)
        @test occursin("GLTFLoader", iss_html)
        cli_out = IOBuffer()
        @test run_cli(["visualize", "--run=$(dir)", "--out=$(joinpath(dir, "cli_iss.html"))", "--no-textures", "--model=1=$(iss)", "--model-scale=2.4"]; io=cli_out) == 0
        @test occursin("\"scale\":2.4", read(joinpath(dir, "cli_iss.html"), String))
        with_path = export_visualization(args; out=joinpath(dir, "paths.html"), textures=false, paths=[(name="ref", points_m=[0.0 10.0; 0.0 0.0; 0.0 1.0], frame=:rtn, target=1)])
        @test occursin("\"paths\":[{", read(with_path, String)) || occursin("\"paths\":[", read(with_path, String))
        ghost_times = [0.0, 100.0, 200.0]
        ghost = export_visualization(args; out=joinpath(dir, "ghost.html"), textures=false,
            references=[(name="ghost", t_s=ghost_times, pos_m=fill(4.0e6, 3, 3), target=1)])
        ghost_html = read(ghost, String)
        @test occursin("\"references\":[{", ghost_html)
        @test occursin("\"name\":\"ghost\"", ghost_html)
        @test occursin("viewer/references.js", ghost_html)
        @test occursin("viewer/video.js", ghost_html) && occursin("mp4-muxer", ghost_html)
        @test occursin("viewer/plots.js", ghost_html)
        @test_throws ArgumentError run_cli(["visualize", "--run=$(dir)", "--model=1"]; io=devnull)

        stl = _write_tiny_stl(joinpath(dir, "bus.stl"))
        custom = export_visualization(args; out=joinpath(dir, "custom.html"), textures=false, frame=:planet_fixed, trail_orbits=5, speed=10.0, stl=Dict(1 => stl), stl_scale=0.01)
        @test isfile(custom)
        custom_html = read(custom, String)
        @test !occursin("data:image/jpeg", custom_html)
        @test occursin("data:model/stl;base64,", custom_html)
        @test occursin("\"trail_orbits\":5.0", custom_html)
        @test_throws ArgumentError export_visualization(prefix; frame=:sideways)

        dev = write_viewer_dev_payload(prefix, joinpath(dir, "dev_data.js"); textures=false)
        @test startswith(read(dev, String), "window.SPACEAGORA_VIEWER = {")

        # CLI: `visualize --run=<dir>` builds the page in-process; `run --visualize` sets the env switch.
        cli_out = IOBuffer()
        @test run_cli(["visualize", "--run=$(dir)", "--out=$(joinpath(dir, "cli.html"))", "--max-frames=20", "--frame=planet_fixed", "--texture=4k", "--no-textures", "--trail-orbits=2"]; io=cli_out) == 0
        @test occursin("viewer=$(joinpath(dir, "cli.html"))", String(take!(cli_out)))
        @test isfile(joinpath(dir, "cli.html"))
        @test_throws ArgumentError run_cli(["visualize"]; io=devnull)
        @test_throws ArgumentError run_cli(["visualize", "--run=$(dir)", "--bogus"]; io=devnull)
        print_out = IOBuffer()
        @test run_cli(["run", "--example=AGORA_Basic_Quickstart.jl", "--visualize", "--print-only"]; io=print_out) == 0
        @test occursin("SPACEAGORA_VISUALIZATION=1", replace(String(take!(print_out)), " => " => "="))

        # A run without the flag has no sidecar to export.
        dir_off = mktempdir()
        run_simulation(_viewer_config(results_directory=dir_off))
        @test_throws ArgumentError export_visualization(joinpath(dir_off, "simulation_results"))

        # visualization=true builds the page in one go, leaving the caller's args untouched.
        dir_kw = mktempdir()
        plain = _viewer_config(results_directory=dir_kw)
        run_simulation(plain; visualization=true)
        @test !plain.simulation_settings.save_visualization_scene
        @test isfile(joinpath(dir_kw, "simulation_results_scene.json"))
        @test isfile(joinpath(dir_kw, "simulation_results_viewer.html"))
    end
end

@testset "ground tracks" begin
    # The option reaches the page, and defaults to off so existing pages are
    # unchanged.
    @test SV._viewer_options(; trail_s=nothing, trail_orbits=nothing, frame=:inertial,
        speed=nothing, title=nothing)["ground_tracks"] === false
    @test SV._viewer_options(; trail_s=nothing, trail_orbits=nothing, frame=:inertial,
        speed=nothing, title=nothing, ground_tracks=true)["ground_tracks"] === true

    dir = mktempdir()
    args = with_visualization_scene(_viewer_config(results_directory=dir), true)
    run_simulation(args)
    prefix = joinpath(dir, "simulation_results")
    on = export_visualization(prefix; out=joinpath(dir, "tracks.html"), textures=false,
        max_frames=20, ground_tracks=true)
    @test occursin("\"ground_tracks\":true", read(on, String))
    off = export_visualization(prefix; out=joinpath(dir, "no_tracks.html"), textures=false, max_frames=20)
    @test occursin("\"ground_tracks\":false", read(off, String))

    # The module is registered in every list that has to carry it, in the same
    # relative position: the bundler, the development harness, the standalone
    # builder. Each shipped entrypoint must carry the same module.
    @test "groundtrack.js" in SV.VIEWER_MODULES
    @test isfile(joinpath(SV.VIEWER_DIR, "src", "groundtrack.js"))
    for (path, needle) in (
        (joinpath(REPO, "viewer", "dev.html"), "\"viewer/groundtrack.js\""),
        (joinpath(REPO, "viewer", "build_standalone.py"), "\"groundtrack.js\""),
    )
        @test occursin(needle, read(path, String))
    end
end
