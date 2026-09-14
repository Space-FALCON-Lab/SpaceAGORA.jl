using Test
using SpaceAGORA
using StaticArrays
using DataFrames
using Arrow
using JSON
using Base64

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
        for body in ("earth", "mars", "moon")
            @test haskey(manifest[body], "8k")
            @test filesize(manifest[body]["8k"]["path"]) < 12_000_000
        end

        @test SV.texture_entry("earth"; resolution="4k")["resolution"] == "4k"
        @test SV.texture_entry("earth"; resolution=:best)["resolution"] == "8k"
        @test SV.texture_entry("earth"; resolution="16k")["resolution"] == "8k"   # largest at or below
        @test SV.texture_entry("venus"; resolution="8k")["resolution"] == "4k"     # only tier
        @test SV.texture_entry("titan"; resolution="2k")["resolution"] == "4k"     # smallest available
        @test_throws ArgumentError SV.texture_entry("earth"; resolution="huge")

        payload = SV.texture_payload("Mars"; resolution="4k")
        @test payload !== nothing
        @test startswith(payload["url"], "data:image/jpeg;base64,")
        @test payload["lon_left_deg"] == -180.0
        @test payload["resolution"] == "4k" && payload["width"] == 4096 && payload["height"] == 2048
        @test SV.texture_payload("Mars")["resolution"] == "8k"
        @test SV.texture_payload("pluto") === nothing
        @test isempty(SV.texture_manifest(mktempdir()))
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
        # A single spacecraft gets Float64 positions (exact to the metre at 7000 km).
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
        many_args = SM.SimulationConfiguration(
            file_paths=args.file_paths, simulation_settings=args.simulation_settings, mission_configuration=args.mission_configuration,
            environment_model=args.environment_model, dynamics_model=SM.DynamicsModel(many, args.dynamics_model.dynamic_effectors),
            guidance_model=args.guidance_model, navigation_model=args.navigation_model, control_model=args.control_model,
            initial_time=args.initial_time, integration_tolerances=args.integration_tolerances, solver_config=args.solver_config
        )
        many_scene = build_visualization_scene(many_args; rotation_max_samples=4)
        many_df = _synthetic_results(many_scene; n_rows=4, with_vel=false, with_link_pose=false)
        many_frames = SV.build_viewer_frames(many_df, many_scene)
        @test many_frames["pos_dtype"] == "f32"
        @test length(_decode_f32(many_frames["pos_km"])) == 4 * length(many) * 3
        @test many_frames["vel_kms"] === nothing
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
        glb = SV.model_payloads(scene; models=Dict(1 => fake_glb), model_scale=3.0)
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
