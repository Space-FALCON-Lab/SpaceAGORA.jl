module LongStreamingTests
# A long-layout run with checkpoints streams its results: each checkpoint
# appends its rows to the results streams and the CSV, and the run keeps only
# its last row in memory. The streamed files hold the rows an unstreamed run
# writes, also when thinned, after a crash and its resume, and across layouts.
using Test
using SpaceAGORA
using SpaceAGORA.SimulationModel
using DataFrames
using Arrow
using CSV
using TOML
using JSON
using Logging
using OrdinaryDiffEq: DiscreteCallback
import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const PLANET = make_no_gram_planet(:earth)

quiet(f) = with_logger(f, ConsoleLogger(stderr, Logging.Error))

# Three spacecraft, 500 s, a checkpoint every 100 s, written to `dir` in `layout`.
function config(dir; layout=:long, checkpoints=true, resume=false, mission_time=500.0)
    spacecraft(index) = make_three_body_spacecraft(bus_dims=(1.0, 1.0, 1.0), panel_dims=(0.01, 1.0, 0.5),
        bus_mass=200.0, panel_mass_each=5.0, panel_offset_y=1.0,
        ic=InitialCondition(ra=7.0e6, rp=7.0e6, i=53.0, ω=0.0, Ω=0.0, ν=30.0 * index), prop_mass=0.0, id=index)
    base = make_example_config(planet=PLANET, spacecraft=spacecraft(1), mission_time=mission_time,
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(SM.InverseSquaredJ2GravityModel(),), density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=false, keplerian=true, EI_km=120.0,
        verbose=false, results=true, results_directory=dir, solver_config=SM.SolverConfig(solver_mode=:tsit5))
    settings = SM.SimConfig._with_fields(base.simulation_settings; results_layout=layout,
        checkpoint_enabled=checkpoints, checkpoint_interval_s=100.0, resume_from_checkpoint=resume)
    return SM.SimConfig._with_configuration(base; simulation_settings=settings,
        dynamics_model=SM.DynamicsModel([spacecraft(k) for k in 1:3], base.dynamics_model.dynamic_effectors))
end

# The default fields plus a run-wide field, so that both tables are streamed.
fields(args) = vcat(SM.default_save_fields(args), [SM.SaveField(:clock, (u, t, integrator) -> t^2)])

run(args; kwargs...) = quiet(() -> run_simulation(args; return_results=true, save_fields=fields(args), kwargs...))

csvs(dir) = (read(joinpath(dir, "simulation_results_spacecraft.csv")), read(joinpath(dir, "simulation_results_run.csv")))
feather(dir, part) = DataFrame(Arrow.Table(joinpath(dir, "simulation_results_$(part).feather")))
entry(args) = SE._load_checkpoint(args).runtime_state.results_rows

# The same run unstreamed: written without the bundle, which a streamed run needs.
unstreamed(args) = withenv(() -> run(args), "SPACEAGORA_SAVE_BUNDLE" => "0")

# Rows of the spacecraft stream beside the checkpoint while the run goes on.
struct Watch
    at::Vector{Float64}
    rows::Vector{Int}
end
function watcher(watch::Watch, path)
    return DiscreteCallback((u, t, integrator) -> any(s -> t > s && !(s in watch.at), (150.0, 250.0, 350.0)),
        integrator -> begin
            push!(watch.at, maximum(s for s in (150.0, 250.0, 350.0) if integrator.t > s))
            push!(watch.rows, isfile(path) ? length(Arrow.Table(path).time) : 0)
        end; save_positions=(false, false))
end

@testset "streaming keeps absent fields distinct from missing measurements" begin
    for short in (1, 3)
        with_scene(args) = SM.SimConfig._with_configuration(args;
            simulation_settings=SM.SimConfig._with_fields(args.simulation_settings; save_visualization_scene=true))
        gap = SM.SaveField(:gap, (u, t, integrator) -> [i == short ? missing : t for i in 1:3];
            per_satellite=true, column_prefix="gap")
        ragged = SM.SaveField(:ragged, (u, t, integrator) -> [i == short ? Float64[] : [t] for i in 1:3];
            per_satellite=true, column_prefix="ragged")
        run_channels(args; kwargs...) = quiet(() -> run_simulation(args; return_results=true,
            save_fields=vcat(fields(args), [gap, ragged]), kwargs...))
        wide_dir, straight_dir, resumed_dir = mktempdir(), mktempdir(), mktempdir()
        wide = run_channels(with_scene(config(wide_dir; layout=:wide)))
        straight = run_channels(with_scene(config(straight_dir)))
        crash = DiscreteCallback((u, t, integrator) -> t > 350.0,
            integrator -> error("metadata resume interruption"); save_positions=(false, false))
        @test_throws Exception run_channels(with_scene(config(resumed_dir)); extra_callbacks=(crash,))
        resumed = run_channels(with_scene(config(resumed_dir; resume=true)))
        again = run_channels(with_scene(config(resumed_dir; resume=true)))
        for result in (straight, resumed, again)
            restored = SM.SceneVisualization._wide_view(result.table, run_table(result), 3)
            same_columns = Set(names(restored)) == Set(names(wide.table))
            @test same_columns
            @test same_columns && all(name -> isequal(restored[!, name], wide.table[!, name]), names(wide.table))
        end
        channels = [(column="gap",), (column="ragged_1",)]
        payload(dir) = JSON.parsefile(export_viewer_payload(joinpath(dir, "simulation_results");
            textures=false, channels=channels))
        expected = payload(wide_dir)
        for dir in (straight_dir, resumed_dir)
            actual = payload(dir)
            @test [channel["name"] for channel in get(actual["frames"], "channels", [])] == ["gap"]
            @test isequal(actual["frames"], expected["frames"])
        end
    end
end

@testset "streaming refuses changed column ownership before appending" begin
    io = SM.IOOutputs
    dir = mktempdir()
    path, csv = joinpath(dir, "test.arrows"), joinpath(dir, "test.csv")
    table = DataFrame(spacecraft=Int32[1, 2], time=[0.0, 0.0], value=[missing, 1.0])
    io._set_long_column_presence!(table, Dict("value" => [2]))
    io._append_long_stream!(path, csv, table, "ownership-test")
    before = (read(path), read(csv))
    io._set_long_column_presence!(table, Dict("value" => [1]))
    @test_throws r"column ownership changed" io._append_long_stream!(path, csv, table, "ownership-test")
    @test (read(path), read(csv)) == before
end

@testset "streamed long results" begin
    @testset "the streamed files hold the unstreamed run's rows" begin
        dir = mktempdir()
        args = config(dir)
        ck = SM.IOConfig._checkpoint_paths(args)
        watch = Watch(Float64[], Int[])
        streamed = run(args; extra_callbacks=(watcher(watch, ck.stream_spacecraft),))
        plain_dir = mktempdir()
        plain = unstreamed(config(plain_dir))
        # The same tables, in memory and on disk, and the same CSV bytes.
        @test isequal(streamed.table, plain.table) && isequal(run_table(streamed), run_table(plain))
        @test isequal(feather(dir, "spacecraft"), plain.table) && isequal(feather(dir, "run"), run_table(plain))
        @test csvs(dir) == csvs(plain_dir)
        # The values are those of the wide layout, rearranged.
        wide = run(config(mktempdir(); layout=:wide))
        @test isequal(streamed.table, long_table(wide))
        @test isequal(run_table(streamed).clock, wide.table.clock)
        # The rows were written as the run went: before the checkpoint at 200 s,
        # those up to 90 s; before the one at 300 s, those up to 190 s.
        @test watch.at == [150.0, 250.0, 350.0] && watch.rows == 3 .* [10, 20, 30]
        # The checkpoint keeps only the last row beside it; the streams hold the rest.
        last_entry = entry(args)
        @test last_entry.version == 2 && last_entry.rows == 1
        @test last_entry.stream.written == nrow(run_table(plain)) - 1 && last_entry.stream.released == nrow(run_table(plain)) - 1
        kept = SM.IOSerialization._read_checkpoint_rows(ck.rows, last_entry.bytes)
        @test kept.times == [500.0]
        # The returned tables read the feathers in place.
        @test !(streamed.table.pos_1 isa Vector) && streamed.table.pos_1 == plain.table.pos_1
        @test ck.stream_spacecraft in streamed.files && ck.stream_run in streamed.files
        # The manifest is the long layout's, written when the run ended.
        manifest = TOML.parsefile(joinpath(dir, "simulation_results.manifest.toml"))
        @test manifest["layout"] == "long" && manifest["steps"] == nrow(run_table(plain)) && !haskey(manifest, "streaming")
        @test manifest["files"]["spacecraft_feather"]["sha256"] ==
            SM.IOSerialization._sha256_hex(joinpath(dir, "simulation_results_spacecraft.feather"))
    end

    @testset "thinned rows" begin
        for (name, value) in (("SPACEAGORA_RESULTS_THIN_STRIDE", "3"), ("SPACEAGORA_RESULTS_FINAL_ONLY", "1"))
            withenv(name => value) do
                dir, plain_dir = mktempdir(), mktempdir()
                streamed = run(config(dir))
                plain = unstreamed(config(plain_dir))
                @test isequal(streamed.table, plain.table) && isequal(run_table(streamed), run_table(plain))
                @test csvs(dir) == csvs(plain_dir)
                @test last(run_table(streamed).time) == 500.0
                name == "SPACEAGORA_RESULTS_FINAL_ONLY" && @test run_table(streamed).time == [500.0]
                name == "SPACEAGORA_RESULTS_THIN_STRIDE" && @test run_table(streamed).time == [0.0:30.0:480.0; 500.0]
            end
        end
    end

    @testset "a crash after the third checkpoint, and the resume" begin
        straight_dir = mktempdir()
        straight = run(config(straight_dir))
        dir = mktempdir()
        crash = DiscreteCallback((u, t, integrator) -> t > 350.0, integrator -> error("crashed"); save_positions=(false, false))
        @test_throws Exception run(config(dir); extra_callbacks=(crash,))
        # The failed run wrote its rows so far, past the checkpoint at 300 s.
        @test last(feather(dir, "run").time) == 300.0
        @test entry(config(dir)).stream.written == 30                     # the checkpoint's rows end before 300 s
        resumed = run(config(dir; resume=true))
        @test isequal(resumed.table, straight.table) && isequal(run_table(resumed), run_table(straight))
        @test isequal(feather(dir, "spacecraft"), straight.table)
        @test csvs(dir) == csvs(straight_dir)
        manifest = TOML.parsefile(joinpath(dir, "simulation_results.manifest.toml"))
        @test manifest["resumed_from_s"] == 300.0 && manifest["first_time_s"] == 0.0
        # A resume of the finished run solves nothing, returns the same tables and
        # leaves the files as they were.
        again = run(config(dir; resume=true))
        @test isequal(again.table, straight.table) && isequal(run_table(again), run_table(straight))
        @test csvs(dir) == csvs(straight_dir) && isequal(feather(dir, "spacecraft"), straight.table)
    end

    @testset "resumes across layouts" begin
        # A wide checkpoint, resumed in the long layout: the rows it kept are
        # streamed with the rest.
        straight = run(config(mktempdir()))
        dir = mktempdir()
        first_leg = run(config(dir; layout=:wide, mission_time=300.0))
        resumed = run(config(dir; resume=true))
        @test isequal(resumed.table, straight.table)
        # A streamed checkpoint, resumed in the wide layout: its rows are in the
        # streams, so the resume warns and writes beside them.
        dir = mktempdir()
        run(config(dir; mission_time=300.0))
        @test_logs (:warn, r"does not stream its results") match_mode=:any run_simulation(
            config(dir; layout=:wide, resume=true); save_fields=fields(config(dir)))
        @test isdir(joinpath(dir, "resumed_at_300.0s"))
        # The streams it leaves beside the checkpoint are kept, under other names.
        ck = SM.IOConfig._checkpoint_paths(config(dir))
        @test !isfile(ck.stream_spacecraft)
        @test isfile(replace(ck.stream_spacecraft, ".arrows" => ".before_resume_at_300.0s.arrows"))
        # A stream another run wrote is not taken, and not cut: it is kept whole.
        dir = mktempdir()
        run(config(dir; mission_time=300.0))
        other = mktempdir()
        run(config(other; mission_time=300.0))
        ck, other_ck = SM.IOConfig._checkpoint_paths(config(dir)), SM.IOConfig._checkpoint_paths(config(other))
        cp(other_ck.stream_spacecraft, ck.stream_spacecraft; force=true)
        @test_logs (:warn, r"written by another run") match_mode=:any run_simulation(
            config(dir; resume=true); save_fields=fields(config(dir)))
        @test read(replace(ck.stream_spacecraft, ".arrows" => ".before_resume_at_300.0s.arrows")) ==
            read(other_ck.stream_spacecraft)
    end

    @testset "a resume that cannot take the rows keeps the earlier streams" begin
        # Without CSVs, the streams and the rows file beside the checkpoint hold
        # a run's rows until it finishes; a resume with other save fields keeps them.
        dir = mktempdir()
        args = config(dir; mission_time=300.0)
        args = SM.SimConfig._with_configuration(args;
            simulation_settings=SM.SimConfig._with_fields(args.simulation_settings; save_csv=false))
        run(args)
        ck = SM.IOConfig._checkpoint_paths(args)
        earlier = Dict(path => read(path) for path in (ck.rows, ck.stream_spacecraft, ck.stream_run))
        resume = SM.SimConfig._with_configuration(config(dir; resume=true);
            simulation_settings=SM.SimConfig._with_fields(config(dir; resume=true).simulation_settings; save_csv=false))
        resumed = @test_logs (:warn, r"kept as .*before_resume_at_300\.0s") match_mode=:any run_simulation(resume;
            return_results=true, save_fields=SM.default_save_fields(resume))
        for (path, bytes) in earlier
            stem, extension = splitext(path)
            @test read(stem * ".before_resume_at_300.0s" * extension) == bytes
        end
        kept = DataFrame(Arrow.Table(replace(ck.stream_spacecraft, ".arrows" => ".before_resume_at_300.0s.arrows")))
        @test first(kept.time) == 0.0 && last(kept.time) == 300.0
        @test first(resumed.table.time) == 300.0                          # the resumed run's own rows
    end

    @testset "a column whose type changes between checkpoints" begin
        # Missing before 150 s and a number after it: the stream's column was
        # written as Missing, and the run stops naming it.
        late = SM.SaveField(:late, (u, t, integrator) -> t < 150.0 ? missing : t)
        args = config(mktempdir())
        @test_throws r"column \"late\" was written as Missing" quiet(() -> run_simulation(args;
            save_fields=vcat(SM.default_save_fields(args), [late])))
        # Missing only in the first segment: later segments' numbers join the
        # stream's nullable column, and the files are those of the unstreamed run.
        early = SM.SaveField(:early, (u, t, integrator) -> t < 50.0 ? missing : t)
        dir, plain_dir = mktempdir(), mktempdir()
        with_early(args) = vcat(SM.default_save_fields(args), [early])
        streamed = quiet(() -> run_simulation(config(dir); return_results=true, save_fields=with_early(config(dir))))
        plain = withenv("SPACEAGORA_SAVE_BUNDLE" => "0") do
            quiet(() -> run_simulation(config(plain_dir); return_results=true, save_fields=with_early(config(plain_dir))))
        end
        @test isequal(run_table(streamed), run_table(plain)) && eltype(run_table(streamed).early) == Union{Missing, Float64}
        @test csvs(dir) == csvs(plain_dir)
        # A resume takes the streams' column types back with the checkpoint, so
        # its numbers join the nullable column too.
        crash = DiscreteCallback((u, t, integrator) -> t > 350.0, integrator -> error("crashed"); save_positions=(false, false))
        crashed_dir = mktempdir()
        @test_throws Exception quiet(() -> run_simulation(config(crashed_dir);
            save_fields=with_early(config(crashed_dir)), extra_callbacks=(crash,)))
        resumed = quiet(() -> run_simulation(config(crashed_dir; resume=true); return_results=true,
            save_fields=with_early(config(crashed_dir))))
        @test isequal(run_table(resumed), run_table(plain)) && csvs(crashed_dir) == csvs(plain_dir)
    end

    @testset "a run killed while writing its files is finished by its resume" begin
        dir, straight_dir = mktempdir(), mktempdir()
        straight = run(config(straight_dir))
        run(config(dir))
        files_before = csvs(dir)
        # As a kill during the final write leaves it: the streams hold the last
        # row past the checkpoint, the manifest still says streaming, no Feather.
        rm.(joinpath.(dir, ["simulation_results_spacecraft.feather", "simulation_results_run.feather"]))
        manifest_path = joinpath(dir, "simulation_results.manifest.toml")
        manifest = TOML.parsefile(manifest_path)
        manifest["streaming"] = true
        open(io -> TOML.print(io, manifest), manifest_path, "w")
        again = run(config(dir; resume=true))
        @test isequal(again.table, straight.table) && isequal(run_table(again), run_table(straight))
        @test isequal(feather(dir, "spacecraft"), straight.table) && csvs(dir) == files_before
        @test !haskey(TOML.parsefile(manifest_path), "streaming")
    end

    @testset "a resume that solves nothing changes nothing it should not" begin
        with_csv(args, save_csv) = SM.SimConfig._with_configuration(args;
            simulation_settings=SM.SimConfig._with_fields(args.simulation_settings; save_csv=save_csv))
        streams(ck) = (read(ck.stream_spacecraft), read(ck.stream_run))
        feathers(dir) = read.(joinpath.(dir, ["simulation_results_spacecraft.feather", "simulation_results_run.feather"]))
        # As a kill during the final write leaves a run: no Feather files, the manifest still streaming.
        function kill_final_write!(dir)
            rm.(joinpath.(dir, ["simulation_results_spacecraft.feather", "simulation_results_run.feather"]))
            manifest_path = joinpath(dir, "simulation_results.manifest.toml")
            manifest = TOML.parsefile(manifest_path)
            manifest["streaming"] = true
            open(io -> TOML.print(io, manifest), manifest_path, "w")
        end
        straight = run(config(mktempdir()))

        # A finished run resumed again, without returning results: nothing changes,
        # and a later resume that returns results returns the whole run.
        dir = mktempdir()
        run(config(dir))
        ck = SM.IOConfig._checkpoint_paths(config(dir))
        before = (streams(ck), feathers(dir), csvs(dir))
        quiet(() -> run_simulation(config(dir; resume=true); save_fields=fields(config(dir))))
        @test (streams(ck), feathers(dir), csvs(dir)) == before
        @test !isdir(joinpath(dir, "resumed_at_500.0s"))
        again = run(config(dir; resume=true))
        @test isequal(again.table, straight.table) && isequal(run_table(again), run_table(straight))

        # A run killed while writing its files, resumed without returning results,
        # with and without CSVs: the resume writes its files from the streams.
        for save_csv in (true, false)
            dir, straight_dir = mktempdir(), mktempdir()
            run(with_csv(config(straight_dir), save_csv))
            run(with_csv(config(dir), save_csv))
            kill_final_write!(dir)
            quiet(() -> run_simulation(with_csv(config(dir; resume=true), save_csv); save_fields=fields(config(dir))))
            @test isequal(feather(dir, "spacecraft"), straight.table) && isequal(feather(dir, "run"), run_table(straight))
            @test !haskey(TOML.parsefile(joinpath(dir, "simulation_results.manifest.toml")), "streaming")
            @test !isdir(joinpath(dir, "resumed_at_500.0s"))
            save_csv && @test csvs(dir) == csvs(straight_dir)
        end

        # The same kill, resumed with other save fields: the resume cannot take the
        # rows, writes nothing, and leaves the streams for a resume that can.
        dir = mktempdir()
        run(config(dir))
        kill_final_write!(dir)
        ck = SM.IOConfig._checkpoint_paths(config(dir))
        kept = streams(ck)
        other = @test_logs (:warn, r"solves nothing, so it returns an empty table and writes no files") match_mode=:any run_simulation(
            config(dir; resume=true); return_results=true, save_fields=SM.default_save_fields(config(dir)))
        @test nrow(other.table) == 0
        @test streams(ck) == kept && !isdir(joinpath(dir, "resumed_at_500.0s"))
        finished = run(config(dir; resume=true))
        @test isequal(finished.table, straight.table) && isequal(feather(dir, "spacecraft"), straight.table)
    end

    @testset "a resume of a run that wrote beside the earlier results" begin
        # A crashes after the checkpoint at 300 s; B resumes it with other save
        # fields, so it writes to resumed_at_300.0s, and crashes after 400 s; C
        # resumes B. C continues B's files, as B would have written them uninterrupted,
        # and leaves A's as they are.
        crash(at) = DiscreteCallback((u, t, integrator) -> t > at, integrator -> error("crashed"); save_positions=(false, false))
        other_fields(args) = vcat(SM.default_save_fields(args), [SM.SaveField(:tock, (u, t, integrator) -> 2t)])
        dir, plain_dir = mktempdir(), mktempdir()
        for d in (dir, plain_dir)
            @test_throws Exception run(config(d); extra_callbacks=(crash(350.0),))
        end
        first_leg = csvs(dir)
        @test_throws Exception quiet(() -> run_simulation(config(dir; resume=true); save_fields=other_fields(config(dir)),
            extra_callbacks=(crash(450.0),)))
        resumed = quiet(() -> run_simulation(config(dir; resume=true); return_results=true, save_fields=other_fields(config(dir))))
        uninterrupted = quiet(() -> run_simulation(config(plain_dir; resume=true); return_results=true,
            save_fields=other_fields(config(plain_dir))))
        @test csvs(dir) == first_leg
        @test isequal(resumed.table, uninterrupted.table) && isequal(run_table(resumed), run_table(uninterrupted))
        @test csvs(joinpath(dir, "resumed_at_300.0s")) == csvs(joinpath(plain_dir, "resumed_at_300.0s"))
        @test isequal(feather(joinpath(dir, "resumed_at_300.0s"), "spacecraft"), uninterrupted.table)
    end

    @testset "a fallback interrupted before its first checkpoint preserves the old run" begin
        snapshot(dir) = Dict(joinpath(root, file) => (read(joinpath(root, file)), stat(joinpath(root, file)).mtime)
            for (root, dirs, files) in walkdir(dir) for file in files)
        crash = DiscreteCallback((u, t, integrator) -> t > 350.0,
            integrator -> error("interrupted fallback"); save_positions=(false, false))
        reference = run(config(mktempdir()))
        for save_csv in (false, true), recording in (false, true)
            dir = mktempdir()
            original = config(dir; mission_time=300.0)
            original = SM.SimConfig._with_configuration(original;
                simulation_settings=SM.SimConfig._with_fields(original.simulation_settings; save_csv))
            run(original)
            checkpoint_dir = dirname(SM.IOConfig._checkpoint_paths(original).rows)
            before = snapshot(checkpoint_dir)
            fallback = config(dir; resume=true)
            fallback = SM.SimConfig._with_configuration(fallback;
                simulation_settings=SM.SimConfig._with_fields(fallback.simulation_settings;
                    save_csv, results=recording))
            @test_throws ErrorException quiet(() -> run_simulation(fallback;
                return_results=recording, save_fields=SM.default_save_fields(fallback), extra_callbacks=(crash,)))
            @test snapshot(checkpoint_dir) == before
            recovered = run(config(dir; resume=true))
            @test isequal(recovered.table, reference.table)
            @test isequal(run_table(recovered), run_table(reference))
        end
    end

    @testset "what does not stream" begin
        # Without checkpoints, or in the wide layout, nothing streams.
        dir = mktempdir()
        run(config(dir; checkpoints=false))
        @test !isfile(SM.IOConfig._checkpoint_paths(config(dir)).stream_spacecraft)
        dir = mktempdir()
        run(config(dir; layout=:wide))
        @test !isfile(SM.IOConfig._checkpoint_paths(config(dir)).stream_spacecraft)
        @test entry(config(dir; layout=:wide)).version == 1
    end
end

end # module LongStreamingTests
