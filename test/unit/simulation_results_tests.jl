module SimulationResultsTests

using Test
using CSV
using DataFrames
using Logging
using Serialization
using StaticArrays
using TOML
using SpaceAGORA
using SpaceAGORA.SimulationModel
using SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft

# Minimal no-GRAM, no-SPICE run: J2 gravity, no atmosphere, one hour.
_test_spacecraft(id::Int; ν::Float64=0.0) = make_three_body_spacecraft(
    bus_dims=(2.05, 2.05, 2.8), panel_dims=(0.01, 2.85, 1.0), bus_mass=620.0,
    panel_mass_each=10.0, panel_offset_y=2.05 / 2.0 + 5.7 / 4.0,
    ic=InitialCondition(ra=7.0e6, rp=7.0e6, i=45.0, ω=0.0, Ω=0.0, ν=ν),
    prop_mass=200.0, id=id)
function _results_test_config(; results::Bool, dir::String, mission_time::Float64=3600.0)
    planet = make_no_gram_planet(:earth)
    sc = _test_spacecraft(1)
    return make_example_config(
        planet=planet, spacecraft=sc, mission_time=mission_time,
        initial_time=InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0),
        dynamic_effectors=(InverseSquaredJ2GravityModel(),),
        density_model=NoAtmosphereModel(), ephemerides_model=SimpleEphemeridesModel(),
        orientation_sim=false, keplerian=true, EI_km=120.0, verbose=false,
        results=results, results_directory=dir)
end

# Same run with checkpointing on (1200 s segments over the 3600 s mission:
# three segments), via SimulationSettings fields.
function _checkpointed_config(; results::Bool, dir::String)
    a = _results_test_config(results=results, dir=dir)
    s = a.simulation_settings
    names = fieldnames(typeof(s))
    vals = NamedTuple{names}(map(n -> getfield(s, n), names))
    s2 = SpaceAGORA.SimulationModel.SimulationSettings(;
        merge(vals, (checkpoint_enabled=true, checkpoint_interval_s=1200.0))...)
    return SpaceAGORA.SimulationModel.SimConfig._with_configuration(a; simulation_settings=s2)
end

@testset "SimulationResults" begin
    csv_dir = mktempdir()
    run_simulation(_results_test_config(results=true, dir=csv_dir))
    csv_path = joinpath(csv_dir, "simulation_results.csv")
    @test isfile(csv_path)
    csv = CSV.read(csv_path, DataFrame)

    mem_dir = mktempdir()
    r = run_simulation(_results_test_config(results=false, dir=mem_dir); return_results=true)
    @test r isa SimulationResults
    @test isempty(readdir(mem_dir))                  # nothing written
    @test isempty(r.files)
    @test r.solution === nothing
    @test nrow(r.table) > 0
    @test names(r.table) == names(csv)
    @test nrow(r.table) == nrow(csv)
    for c in names(csv)
        eltype(csv[!, c]) <: Real || continue
        @test all(isapprox.(Float64.(r.table[!, c]), Float64.(csv[!, c]); rtol=1e-12, atol=1e-12))
    end

    # With results on, the returned table is the written one and the files are listed.
    csv_dir2 = mktempdir()
    r2 = run_simulation(_results_test_config(results=true, dir=csv_dir2); return_results=true)
    @test joinpath(csv_dir2, "simulation_results.csv") in r2.files
    @test names(r2.table) == names(csv)

    # Stale bundle from an earlier run in the same directory is not listed.
    r_stale = run_simulation(_results_test_config(results=false, dir=csv_dir); return_results=true)
    @test isfile(joinpath(csv_dir, "simulation_results.feather"))
    @test isempty(r_stale.files)

    # Checkpointed solve with results=false fills the table from saved_values.
    ck_a = run_simulation(_checkpointed_config(results=true, dir=mktempdir()); return_results=true)
    ck_b_args = _checkpointed_config(results=false, dir=mktempdir())
    ck_b = run_simulation(ck_b_args; return_results=true)
    @test nrow(ck_b.table) > 0
    @test nrow(ck_b.table) == nrow(ck_a.table)
    @test names(ck_b.table) == names(ck_a.table)

    # A completed checkpoint is input to a resume, not output written by it. The
    # run wrote the checkpoint, its manifest and the rows kept beside it.
    @test length(ck_b.files) == 3
    @test any(endswith("simulation_checkpoint.results.bin"), ck_b.files)
    checkpoint_before = Dict(f => (stat(f).mtime, read(f)) for f in ck_b.files)
    settings = ck_b_args.simulation_settings
    fields = fieldnames(typeof(settings))
    values = NamedTuple{fields}(map(n -> getfield(settings, n), fields))
    resume_settings = SimulationSettings(; merge(values, (resume_from_checkpoint=true,))...)
    resume_args = SpaceAGORA.SimulationModel.SimConfig._with_configuration(
        ck_b_args; simulation_settings=resume_settings)
    resumed = run_simulation(resume_args; return_results=true)
    @test isempty(resumed.files)
    @test all(checkpoint_before[f] == (stat(f).mtime, read(f)) for f in ck_b.files)

    # Default return value is unchanged.
    @test run_simulation(_results_test_config(results=false, dir=mktempdir())) === nothing
    sol = run_simulation(_results_test_config(results=false, dir=mktempdir()); return_solution=true)
    @test length(sol.t) > 1

    # return_solution + return_results puts the solution on the struct.
    r3 = run_simulation(_results_test_config(results=false, dir=mktempdir());
                        return_results=true, return_solution=true)
    @test r3.solution !== nothing
    @test length(r3.solution.t) > 1
    @test_throws ArgumentError run_simulation(_results_test_config(results=false, dir=mktempdir());
                                              return_results=true, return_solver_metadata=true)
end

# Failing runs: a callback that throws its own exception, or one that ends the
# solve without success, part way through the mission. What the caller gets and
# what the run leaves behind, for a run in one solve and a checkpointed one
# (1200 s segments).
const DEB = SpaceAGORA.SimulationModel.SimulationCallbacks.DiffEqBase
struct Boom <: Exception
    t::Float64
end
@noinline boom!(integrator) = throw(Boom(integrator.t))
throws_after(t_s) = DEB.DiscreteCallback((u, t, integrator) -> t > t_s, boom!; save_positions=(false, false))
fails_after(t_s) = DEB.DiscreteCallback((u, t, integrator) -> t > t_s,
    integrator -> DEB.terminate!(integrator, DEB.SciMLBase.ReturnCode.Failure); save_positions=(false, false))
function failed_run(config, callback)
    dir = mktempdir()
    args = config(results=true, dir=dir)
    caught = try
        run_simulation(args; extra_callbacks=(callback,))
        nothing
    catch e
        (exception=e, backtrace=catch_backtrace())
    end
    csv = joinpath(dir, "simulation_results.csv")
    times = isfile(csv) ? CSV.read(csv, DataFrame).time : Float64[]
    checkpoint = SpaceAGORA.SimulationEngine._load_checkpoint(args)
    return (; caught, times, checkpoint_t=checkpoint === nothing ? nothing : checkpoint.t)
end

@testset "A failing run reports its failure and leaves the results so far" begin
    # The callback's own exception reaches the caller unchanged, with a
    # backtrace through the callback that threw it, and the results so far are
    # written first.
    for (config, checkpointed) in ((_results_test_config, false), (_checkpointed_config, true))
        run = failed_run(config, throws_after(1500.0))
        @test run.caught.exception isa Boom && run.caught.exception.t > 1500.0
        @test any(frame -> frame.func === :boom!, stacktrace(run.caught.backtrace))
        @test !isempty(run.times) && issorted(run.times)
        if checkpointed
            # The rows of the segments it finished, and the checkpoint written
            # after them: the failing segment's rows are not kept.
            @test maximum(run.times) == 1200.0 && run.checkpoint_t == 1200.0
        else
            # The rows saved up to the failure.
            @test 1400.0 < maximum(run.times) <= run.caught.exception.t && run.checkpoint_t === nothing
        end
    end
    # A solve that ends without success raises an ErrorException naming its
    # retcode, after writing the results so far: a checkpointed run's include
    # the failing segment's rows, and its checkpoint stays the last one written.
    for (config, message, checkpoint_t, segment_end) in
            ((_results_test_config, "Solve failed with retcode=Failure.", nothing, 3600.0),
             (_checkpointed_config, "Checkpointed solve failed with retcode=Failure.", 1200.0, 2400.0))
        run = failed_run(config, fails_after(1500.0))
        @test run.caught.exception isa ErrorException && run.caught.exception.msg == message
        # One solve's rows reach the last save time before the failure; a
        # checkpointed run's also hold the failing step's final state.
        @test !isempty(run.times) && issorted(run.times) && 1500.0 <= maximum(run.times) < segment_end
        @test run.checkpoint_t == checkpoint_t
    end
    # A run that succeeds returns its outcome and one trace entry per solve.
    meta = run_simulation(_results_test_config(results=false, dir=mktempdir()); return_solver_metadata=true)
    @test meta.retcode == "Success" && meta.solution === nothing && length(meta.solver_trace) == 1
end

# Solver metadata of a checkpointed run (1200 s segments over 3600 s), and of runs
# resumed from its final checkpoint: to the same end, to a shorter mission and to a
# longer one.
function _with_settings(args; kwargs...)
    s = args.simulation_settings
    names = fieldnames(typeof(s))
    vals = NamedTuple{names}(map(n -> getfield(s, n), names))
    return SpaceAGORA.SimulationModel.SimConfig._with_configuration(args;
        simulation_settings=SpaceAGORA.SimulationModel.SimulationSettings(; merge(vals, kwargs)...))
end
function _with_mission_time(args, mission_time)
    m = args.mission_configuration
    return SpaceAGORA.SimulationModel.SimConfig._with_configuration(args;
        mission_configuration=SpaceAGORA.SimulationModel.MissionConfiguration(m.mission_type, m.keplerian,
            m.number_of_orbits, mission_time, m.orientation_sim, m.num_steps_to_save, m.data_rate))
end

@testset "Solver metadata of a checkpointed run" begin
    dir = mktempdir()
    args = _checkpointed_config(results=false, dir=dir)
    meta = run_simulation(args; return_solver_metadata=true)
    @test meta.retcode == "Success" && meta.solution === nothing
    @test length(meta.solver_trace) == 3                     # one entry per segment solved
    # Stopped by a termination event in its second segment: two segments solved.
    stops = DEB.DiscreteCallback((u, t, integrator) -> t > 1500.0, DEB.terminate!; save_positions=(false, false))
    stopped = run_simulation(_checkpointed_config(results=false, dir=mktempdir());
        extra_callbacks=(stops,), return_solver_metadata=true)
    @test stopped.retcode == "Terminated" && length(stopped.solver_trace) == 2
    # With return_solution=true, the final segment's solution.
    final = run_simulation(_checkpointed_config(results=false, dir=mktempdir());
        return_solver_metadata=true, return_solution=true)
    @test final.retcode == "Success" && final.solution.t[1] == 2400.0 && final.solution.t[end] == 3600.0
    # Resumed from the checkpoint at or past the mission end: nothing to solve.
    resumed = _with_settings(args; resume_from_checkpoint=true)
    for mission in (resumed, _with_mission_time(resumed, 2400.0))
        none = run_simulation(mission; return_solver_metadata=true, return_solution=true)
        @test none.retcode === nothing && none.solution === nothing && isempty(none.solver_trace)
    end
    # Resumed to a longer mission: one more 1200 s segment.
    longer = run_simulation(_with_mission_time(resumed, 4800.0); return_solver_metadata=true)
    @test longer.retcode == "Success" && length(longer.solver_trace) == 1
end

# RPO command log reachable through the configuration that ran.
module RPOResultsExample
using SpaceAGORA
using SpaceAGORA.SimulationModel
end
const _RPO_SPICE = joinpath(@__DIR__, "..", "..", "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
if !isdir(_RPO_SPICE)
    @info "SPICE kernels absent; skipping the RPO command-log results test" _RPO_SPICE
else
    # The RPO example needs the HYPR package; the unit runner has loaded it already.
    isdefined(Main, :HYPR) || Base.include(Main, joinpath(@__DIR__, "..", "helpers", "load_hypr.jl"))
    Base.include(RPOResultsExample, normpath(joinpath(@__DIR__, "..", "..", "examples", "Earth_RPO_CubeSat_MPC.jl")))
    @testset "SimulationResults RPO command log" begin
        dt = Base.invokelatest() do
            RPOResultsExample.build_rpo_cubesat_mpc_demo(;
            pso_n_particles=8, pso_n_iters=2, n_station_points=256, verbose=false,
            results_directory=mktempdir(), mission_time=10.0, data_rate_s=1.0,
            start_rtn=SVector(-8.0, -4.0, 2.0), goal_rtn=SVector(-5.5, -2.5, 1.0),
            record_control_commands=true)
        end
        # The planner plans in the run; the stop callback ends it 20 s after the reference.
        r = run_simulation(dt.args; return_results=true, extra_callbacks=(dt.stop_callback,))   # default isolate_state=true
        ran = only(c for c in r.configuration.control_model.control_effectors
                   if c isa SpaceAGORA.SimulationModel.RPOMPCControlModel)
        @test ran !== dt.control
        @test ran.command_log !== nothing
        @test !isempty(ran.command_log.t_s)
        @test isempty(dt.control.command_log.t_s)           # caller's object untouched
    end
end

# A resumed run's results hold the whole run: the rows kept beside each checkpoint
# come back first. Against an uninterrupted checkpointed run, for a run extended
# from 2400 s to 3600 s (A), a run that crashed after its 2400 s checkpoint (B),
# and a second resume from the finished checkpoint (C).
function _resumable_config(dir; results::Bool, mission_time::Float64=3600.0, resume::Bool=false)
    a = _results_test_config(results=results, dir=dir, mission_time=mission_time)
    s = a.simulation_settings
    names = fieldnames(typeof(s))
    vals = NamedTuple{names}(map(n -> getfield(s, n), names))
    s2 = SimulationSettings(; merge(vals, (checkpoint_enabled=true, checkpoint_interval_s=1200.0,
        resume_from_checkpoint=resume))...)
    return SpaceAGORA.SimulationModel.SimConfig._with_configuration(a; simulation_settings=s2)
end
_resumed(dir; results, kw...) = run_simulation(_resumable_config(dir; results, resume=true); return_results=true, kw...)
_quietly(f) = with_logger(f, ConsoleLogger(stderr, Logging.Error))
_files(dir) = Dict(f => (mtime(f), read(f)) for f in readdir(dir; join=true) if isfile(f))
_manifest(dir) = TOML.parsefile(joinpath(dir, "simulation_results.manifest.toml"))
_csv(dir) = read(joinpath(dir, "simulation_results.csv"))
# A checkpoint as SpaceAGORA 0.5.1 wrote it: no record of rows beside it.
function _drop_rows_entry!(dir)
    data = joinpath(dir, "checkpoints", "simulation_checkpoint.bin")
    payload = deserialize(data)
    runtime = payload.runtime_state
    kept = NamedTuple{Tuple(k for k in keys(runtime) if k !== :results_rows)}(runtime)
    serialize(data, merge(payload, (runtime_state=kept,)))
end
# A checkpoint whose record names another rows stream than the one beside it.
function _set_rows_token!(dir, token)
    data = joinpath(dir, "checkpoints", "simulation_checkpoint.bin")
    payload = deserialize(data)
    runtime = payload.runtime_state
    entry = merge(runtime.results_rows, (token=token,))
    serialize(data, merge(payload, (runtime_state=merge(runtime, (results_rows=entry,)),)))
end

@testset "A resumed run's results hold the whole run" begin
    for results in (true, false)
        reference_dir = mktempdir()
        reference = run_simulation(_resumable_config(reference_dir; results); return_results=true)
        @test nrow(reference.table) == 361 && reference.table.time[1] == 0.0
        whole(res) = isequal(res.table, reference.table)
        same_files(dir) = !results || (_csv(dir) == _csv(reference_dir) &&
            _manifest(dir)["steps"] == 361 && _manifest(dir)["first_time_s"] == 0.0)

        # A: extend the run.
        dir_a = mktempdir()
        run_simulation(_resumable_config(dir_a; results, mission_time=2400.0); return_results=true)
        a = _resumed(dir_a; results)
        @test whole(a) && same_files(dir_a)
        @test !results || _manifest(dir_a)["resumed_from_s"] == 2400.0
        @test !haskey(results ? _manifest(reference_dir) : Dict(), "resumed_from_s")

        # B: a crash after the 2400 s checkpoint, then a resume.
        dir_b = mktempdir()
        @test_throws Boom _quietly(() -> run_simulation(_resumable_config(dir_b; results);
            extra_callbacks=(throws_after(2500.0),), return_results=true))
        b = _resumed(dir_b; results)
        @test whole(b) && same_files(dir_b)

        # C: a second resume from the finished checkpoint solves nothing, writes
        # nothing, and returns the whole run.
        before = _files(dir_a)
        c = _resumed(dir_a; results)
        @test whole(c) && isempty(c.files) && _files(dir_a) == before
        # Likewise to a mission that ends before the checkpoint: the rows up to it.
        past = run_simulation(_resumable_config(dir_a; results, mission_time=2400.0, resume=true); return_results=true)
        @test whole(past) && isempty(past.files) && _files(dir_a) == before

        # Rows a crash wrote after the last checkpoint are cut away on resume.
        dir_d = mktempdir()
        run_simulation(_resumable_config(dir_d; results, mission_time=2400.0); return_results=true)
        rows_path = joinpath(dir_d, "checkpoints", "simulation_checkpoint.results.bin")
        SpaceAGORA.SimulationModel.IOSerialization._append_checkpoint_rows!(rows_path, [2410.0],
            [SpaceAGORA.SimulationModel.SaveData()])
        @test whole(_resumed(dir_d; results))
    end

    # A checkpoint without rows (SpaceAGORA 0.5.1 or earlier) still resumes, with a
    # warning; its results start at the resume and go beside the earlier ones,
    # which stay as they were.
    dir = mktempdir()
    run_simulation(_resumable_config(dir; results=true, mission_time=2400.0))
    _drop_rows_entry!(dir)
    before = _files(dir)
    old = @test_logs (:warn, r"cannot include the rows before its checkpoint at t = 2400\.0 s") match_mode=:any _resumed(dir; results=true)
    @test _files(dir) == before
    folder = joinpath(dir, "resumed_at_2400.0s")
    @test old.table.time[1] == 2400.0 && nrow(old.table) == 121
    @test isfile(joinpath(folder, "simulation_results.csv")) && joinpath(folder, "simulation_results.csv") in old.files
    @test _manifest(folder)["first_time_s"] == 2400.0 && _manifest(folder)["resumed_from_s"] == 2400.0
    @test old.configuration.simulation_settings.results_directory == folder

    # A rows stream another run wrote is not mixed in: the stream and the
    # checkpoint carry the same token.
    dir = mktempdir()
    run_simulation(_resumable_config(dir; results=true, mission_time=2400.0))
    _set_rows_token!(dir, "another-run")
    before = _files(dir)
    stray = @test_logs (:warn, r"written by another run") match_mode=:any _resumed(dir; results=true)
    @test _files(dir) == before && stray.table.time[1] == 2400.0

    # A resume that solves nothing and returns no table reads no rows: a damaged
    # rows file draws no warning, and no file is written.
    dir = mktempdir()
    run_simulation(_resumable_config(dir; results=true, mission_time=2400.0))
    write(joinpath(dir, "checkpoints", "simulation_checkpoint.results.bin"), "damaged")
    before = _files(dir)
    @test_logs min_level=Logging.Warn run_simulation(_resumable_config(dir; results=true, mission_time=2400.0, resume=true))
    @test _files(dir) == before
    # Returning the table, it warns that the rows cannot be read, and the table is empty.
    empty_table = @test_logs (:warn, r"solves nothing, so it returns an empty table") match_mode=:any run_simulation(
        _resumable_config(dir; results=true, mission_time=2400.0, resume=true); return_results=true)
    @test nrow(empty_table.table) == 0 && _files(dir) == before

    # Rows saved with other save fields are not mixed in either.
    dir = mktempdir()
    run_simulation(_resumable_config(dir; results=true, mission_time=2400.0))
    before = _files(dir)
    fields = SpaceAGORA.SimulationModel.default_save_fields(_resumable_config(dir; results=true))[1:2]
    other = @test_logs (:warn, r"other save fields") match_mode=:any _resumed(dir; results=true, save_fields=fields)
    @test _files(dir) == before && other.table.time[1] == 2400.0

    # A resume that cannot include the earlier rows keeps the rows file beside the
    # checkpoint, which after a killed run is their only copy, and starts its own.
    dir = mktempdir()
    run_simulation(_resumable_config(dir; results=true, mission_time=2400.0))
    rows_path = joinpath(dir, "checkpoints", "simulation_checkpoint.results.bin")
    kept = joinpath(dir, "checkpoints", "simulation_checkpoint.results.before_resume_at_2400.0s.bin")
    earlier = read(rows_path)
    other = @test_logs (:warn, r"rows file beside the checkpoint is kept as .*before_resume_at_2400\.0s\.bin") match_mode=:any _resumed(
        dir; results=true, save_fields=fields)
    @test read(kept) == earlier && other.table.time[1] == 2400.0
    @test isfile(rows_path) && read(rows_path) != earlier
    kept_rows = SpaceAGORA.SimulationModel.IOSerialization._read_checkpoint_rows(kept, filesize(kept))
    @test first(kept_rows.times) == 0.0 && last(kept_rows.times) == 2400.0
    # A second such resume keeps the first one's rows under the next name.
    run_simulation(_resumable_config(dir; results=true, mission_time=2400.0))
    again = read(rows_path)
    @test_logs (:warn, r"before_resume_at_2400\.0s_2\.bin") match_mode=:any _resumed(dir; results=true, save_fields=fields)
    @test read(kept) == earlier && read(replace(kept, ".bin" => "_2.bin")) == again
    # A resume that saves no results keeps the rows file rather than removing it.
    dir = mktempdir()
    run_simulation(_resumable_config(dir; results=true, mission_time=2400.0))
    rows_path = joinpath(dir, "checkpoints", "simulation_checkpoint.results.bin")
    earlier = read(rows_path)
    run_simulation(_resumable_config(dir; results=false, resume=true))
    @test !isfile(rows_path)
    @test read(joinpath(dir, "checkpoints", "simulation_checkpoint.results.before_resume_at_2400.0s.bin")) == earlier
    # A resume that continues the rows file renames nothing.
    dir = mktempdir()
    run_simulation(_resumable_config(dir; results=true, mission_time=2400.0))
    _resumed(dir; results=true)
    @test readdir(joinpath(dir, "checkpoints")) ==
        ["simulation_checkpoint.bin", "simulation_checkpoint.manifest.toml", "simulation_checkpoint.results.bin"]

    # A fallback resume can itself fail before its first checkpoint. Until it
    # reaches that checkpoint, the old checkpoint must still find its own rows.
    # Exercise both paths that would otherwise replace/remove the old stream.
    for recording in (true, false)
        dir = mktempdir()
        run_simulation(_resumable_config(dir; results=true, mission_time=2400.0))
        checkpoint_dir = joinpath(dir, "checkpoints")
        before_checkpoint = _files(checkpoint_dir)
        @test_throws Boom _quietly(() -> run_simulation(
            _resumable_config(dir; results=recording, resume=true);
            return_results=recording, save_fields=fields, extra_callbacks=(throws_after(2500.0),)))
        @test _files(checkpoint_dir) == before_checkpoint
        recovered = _resumed(dir; results=true)
        @test recovered.table.time[1] == 0.0 && nrow(recovered.table) == 361
        @test readdir(checkpoint_dir) == sort!(collect(basename(f) for f in keys(before_checkpoint)))
    end

    # A no-op fallback does not rename a stream even when it returns a table.
    dir = mktempdir()
    run_simulation(_resumable_config(dir; results=true, mission_time=2400.0))
    checkpoint_dir = joinpath(dir, "checkpoints")
    before_checkpoint = _files(checkpoint_dir)
    _quietly(() -> run_simulation(
        _resumable_config(dir; results=true, mission_time=2400.0, resume=true);
        return_results=true, save_fields=fields))
    @test _files(checkpoint_dir) == before_checkpoint

    # Ensemble members keep their rows beside their own checkpoints.
    function pair(dir; mission_time, resume=false)
        a = _resumable_config(dir; results=true, mission_time, resume)
        second = _test_spacecraft(2; ν=90.0)
        return SpaceAGORA.SimulationModel.SimConfig._with_configuration(a;
            dynamics_model=DynamicsModel([a.dynamics_model.spacecraft[1], second], a.dynamics_model.dynamic_effectors))
    end
    reference_dir = mktempdir()
    run_constellation_ensemble(pair(reference_dir; mission_time=3600.0))
    dir = mktempdir()
    run_constellation_ensemble(pair(dir; mission_time=2400.0))
    run_constellation_ensemble(pair(dir; mission_time=3600.0, resume=true))
    for member in ("sat_1_id_1", "sat_2_id_2")
        @test _csv(joinpath(dir, member)) == _csv(joinpath(reference_dir, member))
    end
end

end # module SimulationResultsTests
