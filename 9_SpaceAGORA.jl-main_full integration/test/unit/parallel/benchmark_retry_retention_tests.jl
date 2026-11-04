module BenchmarkRetryRetentionTests
using Test
using CSV
using DataFrames
using Dates
using Statistics

const RUNTIME = normpath(joinpath(@__DIR__, "..", "..", "..", "benchmarks", "studies", "performance_runtime_analysis"))
const REPORT_FIXTURE = joinpath(@__DIR__, "benchmark_report_tests.jl")

# Select definitions only, including helpers inside an existing fixture module.
# No other fixture testsets, native imports or benchmark entrypoints execute.
function load_definition(file::String, wanted::Symbol; required=true)
    nodes = Meta.parseall(read(file, String)).args
    for expr in nodes
        if expr isa Expr && expr.head == :module
            nodes = expr.args[end].args
            break
        end
    end
    for expr in nodes
        expr isa Expr || continue
        node = expr.head == :macrocall ? expr.args[end] : expr
        node isa Expr || continue
        name = if node.head == :struct
            node.args[2]
        elseif node.head in (:function, :(=))
            signature = node.args[1]
            signature isa Expr || continue
            signature = signature.head == :(::) ? signature.args[1] : signature
            signature.head == :call ? signature.args[1] : nothing
        else
            nothing
        end
        if name == wanted
            Core.eval(@__MODULE__, expr)
            return true
        end
    end
    required && error("Required retry fixture definition $wanted not found in $file")
    return false
end

const SimulationConfiguration = NamedTuple
MissionConfiguration(; kwargs...) = (; kwargs...)
struct Earth end
struct Mars end
module SimulationModel
module SimConfig
_with_configuration(args; kwargs...) = merge(args, (; kwargs...))
end
end
module SimulationEngine
struct SolverIntegratorCache end
end
module ParallelProfiles
const events = NamedTuple[]
function record_outer_route_feedback!(state, features; kwargs...)
    push!(events, merge((scenario=features,), (; kwargs...)))
end
end
const _PERF_OUTER_ROUTE_STATE = Ref(nothing)
_outer_route_features(case; kwargs...) = case.name
_outer_route_tuning() = nothing
const HELPERS = joinpath(RUNTIME, "case_catalog", "profile_types_and_helpers.jl")
for name in (:ProfileSpec, :BenchmarkCase, :ParallelPriorityPlan, :_safe_div,
             :_safe_unique_join, :_parse_bool_env, :_perf_smoke_mode,
             :_perf_case_name_filter, :_exclude_entry_scenarios)
    load_definition(HELPERS, name)
end
const BACKEND = Ref(:none)
const MC_ENABLED = Ref(false)
const SCENARIOS = Ref([(name="mc_retry", variant=:retry)])
perf_parallel_backend() = BACKEND[]
parallel_priority_plan(case, route) = ParallelPriorityPlan(outer_route=route)
parallel_priority_env_pairs(plan) = Pair{String,String}[]
case_env_pairs(case, plan) = Pair{String,String}[]
_case_outer_threads_safe(case) = !startswith(case.name, "unsafe")
auto_backend_for_case(case; kwargs...) = startswith(case.name, "process") ? :process :
    startswith(case.name, "thread") ? :threads : :none
_split_rollout_benchmark_cases(cases) = cases
_multirate_rollout_benchmark_cases(cases) = cases
_include_montecarlo_scenarios() = MC_ENABLED[]
_active_montecarlo_scenarios() = SCENARIOS[]
_montecarlo_batch_mission_time_s(spec, variant) = spec.montecarlo_mission_s
_include_control_stress_per_orbit() = false
_perf_stream_orbit_logs() = false
make_spacecraft(planet; kwargs...) = nothing
orbital_period_seconds(sc, planet) = 10.0
perf_worker_planet() = Earth()
perf_worker_mars() = nothing
const WORKER_TRANSPORT = Ref(0)
ensure_perf_workers!() = (WORKER_TRANSPORT[] += 1; nothing)
workers() = [101, 102]
remotecall_wait(f, worker, args...) = f(args...)
pmap(f, inputs) = map(f, inputs)

function make_spec(; repeats=1, attempts=2, seeds=2, warmup=1)
    ProfileSpec(name="quick", repeats=repeats, warmup=warmup, max_attempts=attempts,
        mission_short_s=120.0, mission_long_s=600.0, montecarlo_samples=seeds,
        montecarlo_mission_s=120.0)
end
const SPEC = make_spec()
function fixture_case(name="retry_case"; behavior=:retry, category="baseline", mission_time=120.0, quick=true)
    mission = MissionConfiguration(mission_type=:fixture, keplerian=nothing,
        number_of_orbits=1, mission_time=mission_time, orientation_sim=false,
        num_steps_to_save=2)
    BenchmarkCase(name=name, category=category, description="Synthetic retry fixture",
        args_template=(mission_configuration=mission, behavior=behavior), run_in_quick=quick)
end
make_montecarlo_case(seed, mission_time, variant, planet; kwargs...) =
    fixture_case("mc_$variant"; behavior=variant, category="montecarlo", mission_time=mission_time)

const RECORD_LOCK = ReentrantLock()
const WARMUPS = NamedTuple[]
const MEASUREMENTS = NamedTuple[]
function run_warmup(case, count, profile; kwargs...)
    lock(RECORD_LOCK) do
        push!(WARMUPS, (scenario=case.name, count=count,
            mission_time=case.args_template.mission_configuration.mission_time,
            interface=case.entry_target_count_override))
    end
    nothing
end
# Reuse the established synthetic report schema, not an independent accounting
# implementation. These measurement values are deterministic and not timings.
const PERF_BASELINE_SCENARIO = "baseline_fixture"
load_definition(REPORT_FIXTURE, :raw_row)
function measure_case(case, profile, rep; seed=missing, attempt=1, plan=nothing, kwargs...)
    behavior = case.args_template.behavior
    success = behavior == :early || (behavior in (:retry, :metadata) && attempt == 2)
    elapsed = success ? 2.0 : 100.0
    lock(RECORD_LOCK) do
        push!(MEASUREMENTS, (scenario=case.name, repeat=rep, seed=seed, attempt=attempt,
            mission_time=case.args_template.mission_configuration.mission_time,
            interface=case.entry_target_count_override))
    end
    row = raw_row(category=case.category, scenario=case.name, description=case.description,
        mission_time_s=case.args_template.mission_configuration.mission_time,
        repeat=rep, seed=seed, attempt=attempt, solve_success=success,
        solve_retcode=success ? "Success" : "FixtureFailure", total_time_s=elapsed,
        copy_time_s=0.0, solve_time_s=elapsed, terminal_time_s=success ? 20.0 : 999.0,
        outer_route=string(plan.outer_route))
    return behavior == :metadata ? merge(row, (
        solver_mode=success ? "auto_stiff" : missing,
        solver_sequence=success ? "Tsit5->Rodas5P" : missing,
        solver_fallback_used=success, solver_fallback_count=success ? 1 : missing)) : row
end

const ROUTES = joinpath(RUNTIME, "case_catalog", "env_workers_and_routes.jl")
for name in (:_env_pairs_key, :_thread_plan_groups, :_split_threadsafe_indices,
             :_record_outer_route_feedback!)
    load_definition(ROUTES, name)
end
const BATCH = joinpath(RUNTIME, "measurement", "warmup_and_batch.jl")
for name in (:_perf_optional_nonnegative_int_env, :_case_sample_schedule,
             :_run_case_batch_core!, :run_case_batch!)
    load_definition(BATCH, name)
end
include(joinpath(RUNTIME, "measurement", "montecarlo.jl"))
include(joinpath(RUNTIME, "measurement", "per_orbit.jl"))
include(joinpath(RUNTIME, "measurement", "entry_duration.jl"))
include(joinpath(RUNTIME, "reporting", "summaries.jl"))
for name in (:_fmt, :_scenario_metric)
    load_definition(joinpath(RUNTIME, "reporting", "plots_and_reporting.jl"), name)
end
const FALLBACK = (hardware_class="fixture", machine_label="fixture", host_name="fixture",
    cpu_name="fixture", cpu_threads=2, julia_threads=2, os="fixture", arch="fixture")
for name in (:_runtime_hardware_snapshot, :_perf_default_solver_mode, :_perf_solver_mode_env,
    :_split_rollout_enabled, :_split_rollout_enforce, :_split_rollout_case_names,
    :_split_rollout_solver_variants, :_multirate_rollout_enabled, :_multirate_rollout_enforce,
    :_multirate_rollout_case_names, :_multirate_rollout_max_slowdown_ratio, :render_report)
    load_definition(REPORT_FIXTURE, name)
end
include(joinpath(RUNTIME, "reporting", "report_writing.jl"))

function reset!(backend=:none; mc=false)
    BACKEND[] = backend; MC_ENABLED[] = mc; WORKER_TRANSPORT[] = 0
    empty!(WARMUPS); empty!(MEASUREMENTS); empty!(ParallelProfiles.events)
end
function captured(f)
    mktemp() do _, io
        result = redirect_stdout(f, io)
        flush(io); seekstart(io)
        return result, read(io, String)
    end
end
function feedback_expected(successes, failures; route=nothing)
    @test length(ParallelProfiles.events) == length(successes)
    @test [e.successes for e in ParallelProfiles.events] == successes
    @test [e.failures for e in ParallelProfiles.events] == failures
    @test [e.elapsed_success_s for e in ParallelProfiles.events] == 2.0 .* successes
    @test [e.elapsed_success_sq_sum_s for e in ParallelProfiles.events] == 4.0 .* successes
    route === nothing || @test all(e.route == route for e in ParallelProfiles.events)
end
function accounting(df; requests, successes, failures, cost, mean_time)
    before = deepcopy(df)
    summary = summarize_results(df)
    @test isequal(df, before)
    @test sum(summary.samples_total) == successes + failures
    @test sum(summary.samples_success) == successes
    @test sum(summary.samples_failed) == failures
    @test sum(summary.requested_runs) == requests
    @test sum(summary.retries_total) == successes + failures - requests
    @test sum(summary.total_time_all_attempts_s) == cost
    @test all(v -> isequal(v, mean_time), summary.total_time_mean_s)
    @test all(summary.penalized_expected_wall_time_s .== cost / requests)
    mktempdir() do dir
        path = joinpath(dir, "attempts.csv")
        CSV.write(path, df)
        restored = CSV.read(path, DataFrame)
        @test isequal(restored, df)
        @test isequal(summarize_results(restored), summary)
        report = render_report(restored, summary)
        @test occursin("Failed attempts: `$failures/$(successes + failures)`", report)
        @test occursin("Counts and success/failure rates describe recorded attempts.", report)
        @test occursin("costs exclude warmups, worker startup, explicit GC and later processing", report)
        @test occursin("not complete campaign wall time", report)
        @test occursin("Robustness-adjusted expected wall time (all attempts): `$(cost / requests) s/requested run`", report)
        successes == 0 && @test occursin("No successful runs were recorded.", report)
    end
    summary
end

withenv("SPACEAGORA_PERF_WARMUP_OVERRIDE" => nothing,
        "SPACEAGORA_PERF_REPEATS_OVERRIDE" => nothing, "SPACEAGORA_PERF_CASES" => "",
        "SPACEAGORA_PERF_SMOKE" => "0", "SPACEAGORA_PERF_EXCLUDE_ENTRY_SCENARIOS" => "0") do
    @testset "deterministic returned attempts and unchanged terminal feedback" begin
        for (behavior, expected_attempts, success, failures, cost, mean_time) in
            ((:retry, [1, 2], 1, 1, 102.0, 2.0), (:fail, [1, 2], 0, 2, 200.0, missing),
             (:early, [1], 1, 0, 2.0, 2.0))
            reset!()
            df, log = captured(() -> run_benchmarks(SPEC, [fixture_case(; behavior=behavior)], Earth()))
            @test df.attempt == expected_attempts
            @test :is_terminal_attempt in propertynames(df) && df.is_terminal_attempt == (behavior == :early ? [true] : [false, true])
            @test length(MEASUREMENTS) == length(expected_attempts)
            @test length(WARMUPS) == 1 && only(WARMUPS).count == 1
            accounting(df; requests=1, successes=success, failures=failures, cost=cost, mean_time=mean_time)
            feedback_expected([success], [success == 0 ? 1 : 0]; route=:none)
            @test occursin(success == 0 ? "failed after 2 attempts" : "total=2.0 s", log)
        end
        reset!()
        limited, _ = captured(() -> run_benchmarks(make_spec(attempts=1), [fixture_case()], Earth()))
        @test limited.attempt == [1] && limited.solve_success == [false]
        @test length(MEASUREMENTS) == 1
        reset!()
        empty_df, _ = captured(() -> run_benchmarks(SPEC, BenchmarkCase[], Earth()))
        @test isempty(empty_df) && isempty(MEASUREMENTS) && isempty(ParallelProfiles.events)
        reset!()
        empty_df, _ = captured(() -> run_benchmarks(make_spec(repeats=0), [fixture_case()], Earth()))
        @test isempty(empty_df) && isempty(MEASUREMENTS)
        @test length(WARMUPS) == 1
    end

    @testset "solver metadata changes do not split one retried request" begin
        reset!()
        df, _ = captured(() -> run_benchmarks(SPEC, [fixture_case(; behavior=:metadata)], Earth()))
        @test isequal(df.solver_mode, [missing, "auto_stiff"])
        @test isequal(df.solver_sequence, [missing, "Tsit5->Rodas5P"])
        @test isequal(df.solver_fallback_count, [missing, 1])
        summary = accounting(df; requests=1, successes=1, failures=1, cost=102.0, mean_time=2.0)
        @test nrow(summary) == 1
        @test only(summary.solver_modes) == "auto_stiff"
        @test only(summary.solver_sequences) == "Tsit5->Rodas5P"
        @test only(summary.solver_fallback_any)
        @test only(summary.solver_fallback_count_mean) == 1.0
        @test only(summary.fallback_count_mean_all_attempts) == 0.5
        feedback_expected([1], [0]; route=:none)
    end

    @testset "collector routes retain order and terminal-only feedback" begin
        for backend in (:none, :threads, :process, :auto)
            reset!(backend)
            cases = [fixture_case("process_retry"), fixture_case("thread_retry"), fixture_case("serial_retry")]
            df, _ = captured(() -> run_benchmarks(make_spec(repeats=2), cases, Earth()))
            @test df.scenario == repeat([c.name for c in cases]; inner=4)
            @test df.repeat == repeat([1, 1, 2, 2], 3)
            @test df.attempt == repeat([1, 2], 6)
            @test :is_terminal_attempt in propertynames(df) && df.is_terminal_attempt == repeat([false, true], 6)
            feedback_expected([2, 2, 2], [0, 0, 0])
            @test sum(df.total_time_s) == 612.0
            @test length(MEASUREMENTS) == 12
            backend in (:process, :auto) && @test WORKER_TRANSPORT[] > 0
        end
    end

    @testset "ordinary orbit and entry sweep identities" begin
        for backend in (:none, :threads, :process, :auto)
            reset!(backend)
            cases = [fixture_case("process_retry"), fixture_case("thread_retry")]
            df, _ = captured(() -> run_per_orbit_for_scenarios(SPEC, cases, Earth()))
            @test df.scenario == repeat([c.name for c in cases]; inner=6)
            @test df.mission_time_multiplier == repeat([1, 1, 2, 2, 3, 3], 2)
            @test df.orbit_count == df.mission_time_multiplier
            @test df.mission_time_s == 10.0 .* df.mission_time_multiplier
            @test df.attempt == repeat([1, 2], 6)
            @test :is_terminal_attempt in propertynames(df) && df.is_terminal_attempt == repeat([false, true], 6)
            feedback_expected([3, 3], [0, 0])
            summary = summarize_per_orbit_results(df)
            @test all(summary.samples_total .== 2)
            @test all(summary.samples_failed .== 1)
            @test all(summary.total_time_mean_s .== 2.0)
            reset!(backend)
            entries = [fixture_case("process_entry"; category="entry"), fixture_case("thread_entry"; category="entry")]
            entry, _ = captured(() -> run_entry_duration_sweep(SPEC, entries, Earth()))
            @test entry.entry_run_role == vcat(fill("reference", 8), fill("measured", 8))
            @test entry.attempt == repeat([1, 2], 8)
            @test :is_terminal_attempt in propertynames(entry) && entry.is_terminal_attempt == repeat([false, true], 8)
            @test entry.entry_atmospheric_interface_count == repeat([1, 1, 2, 2], 4)
            @test all(entry.entry_reference_terminal_time_s .== 20.0)
            @test entry.entry_event_time_abs_error_s == repeat([979.0, 0.0], 8)
            feedback_expected([2, 2], [0, 0])
            entry_summary = summarize_entry_duration_results(entry)
            @test all(entry_summary.samples_total .== 2)
            @test all(entry_summary.samples_failed .== 1)
            @test all(entry_summary.total_time_mean_s .== 2.0)
            @test all(entry_summary.event_time_abs_error_mean_s .== 0.0)
            @test count(w -> w.count == 1, WARMUPS) == 8
        end
        reset!()
        failed = fixture_case("failed_entry"; category="entry", behavior=:fail)
        entry, _ = captured(() -> run_entry_duration_sweep(SPEC, [failed], Earth()))
        @test nrow(entry) == 8 && all(.!entry.solve_success)
        @test all(ismissing, entry.entry_reference_terminal_time_s)
        @test all(ismissing, entry.entry_event_time_abs_error_s)
        @test all(ismissing, summarize_entry_duration_results(entry).total_time_mean_s)
        feedback_expected([0], [2])
        @test isempty(measure_per_orbit_scenario(fixture_case(), SPEC, 10.0, Int[])[1])
        @test isempty(measure_entry_duration_scenario(fixture_case(category="entry"), SPEC, Int[])[1])
    end

    @testset "legacy feedback and Monte Carlo terminal-row compatibility" begin
        reset!()
        legacy = NamedTuple[(solve_success=false, total_time_s=100.0),
                            (solve_success=true, total_time_s=2.0)]
        _record_outer_route_feedback!(fixture_case(), legacy; route=:none)
        feedback_expected([1], [1]; route=:none)
        reset!()
        marked = NamedTuple[(solve_success=false, total_time_s=100.0, is_terminal_attempt=false),
                            (solve_success=true, total_time_s=2.0, is_terminal_attempt=true)]
        _record_outer_route_feedback!(fixture_case(), marked; route=:none)
        feedback_expected([1], [0]; route=:none)
        for behavior in (:retry, :fail, :early)
            reset!()
            terminal, err = measure_montecarlo_seed(SPEC, Earth(), 120.0, 1001; variant=behavior)
            @test terminal.attempt == (behavior == :early ? 1 : 2)
            @test terminal.is_terminal_attempt
            @test terminal.solve_success == (behavior != :fail)
            @test behavior == :fail ? occursin("failed after 2 attempts", err) : err === nothing
            worker_terminal, worker_err = perf_worker_measure_montecarlo_seed(SPEC, 120.0, 1001, behavior)
            @test worker_terminal.attempt == terminal.attempt
            @test worker_terminal.is_terminal_attempt
            @test isequal(worker_err, err)
            attempts, all_err = perf_worker_measure_montecarlo_seed(SPEC, 120.0, 1001, behavior; retain_attempts=true)
            @test [row.attempt for row in attempts] == (behavior == :early ? [1] : [1, 2])
            @test isequal(all_err, err)
        end
        reset!()
        custom, err = measure_montecarlo_seed(SPEC, Earth(), 120.0, 1001;
            variant=:retry, mars=Mars(), outer_route=:process,
            plan=ParallelPriorityPlan(outer_route=:threads), apply_env=false)
        @test err === nothing && custom.solve_success
        @test custom.outer_route == "threads"
        @test custom.attempt == 2
        reset!()
        terminal, err = measure_montecarlo_seed(make_spec(attempts=0), Earth(), 120.0, 1001;
            variant=:retry, apply_env=false)
        @test terminal === nothing
        @test err == "failed without attempt data"
        @test isempty(MEASUREMENTS)
    end

    @testset "Monte Carlo ordered attempts across serial thread and process collectors" begin
        for backend in (:none, :threads, :process), per_orbit in (false, true), behavior in (:retry, :fail, :early)
            reset!(backend; mc=true)
            SCENARIOS[] = [(name="mc_$behavior", variant=behavior)]
            rows = NamedTuple[]
            _, log = captured() do
                if per_orbit
                    run_montecarlo_per_orbit!(rows, SPEC, Earth(), 10.0, [1, 2])
                else
                    run_montecarlo_batch!(rows, SPEC, Earth())
                end
            end
            df = DataFrame(rows)
            attempts = behavior == :early ? 1 : 2
            sweeps = per_orbit ? 2 : 1
            seed_ids = per_orbit ? [1, 2] : [1001, 1002]
            @test df.seed == repeat(repeat(seed_ids; inner=attempts), sweeps)
            @test df.attempt == repeat(collect(1:attempts), 2sweeps)
            @test :is_terminal_attempt in propertynames(df) && df.is_terminal_attempt == repeat(behavior == :early ? [true] : [false, true], 2sweeps)
            @test length(MEASUREMENTS) == 2sweeps * attempts
            @test sum(df.solve_success) == (behavior == :fail ? 0 : 2sweeps)
            @test sum(df.total_time_s) == 2sweeps * (behavior == :retry ? 102.0 : behavior == :fail ? 200.0 : 2.0)
            feedback_expected(fill(behavior == :fail ? 0 : 2, sweeps), fill(behavior == :fail ? 2 : 0, sweeps); route=backend)
            if per_orbit
                @test df.orbit_count == repeat([1, 2]; inner=2attempts)
                @test df.mission_time_multiplier == df.orbit_count
                @test df.mission_time_s == 10.0 .* df.orbit_count
            else
                @test length(WARMUPS) == (backend == :process ? 2 : 1)
            end
            @test findfirst("seed 1/2=$(seed_ids[1])", log).start < findfirst("seed 2/2=$(seed_ids[2])", log).start
            backend == :process && @test WORKER_TRANSPORT[] > 0
        end
        reset!(:none; mc=false)
        rows = NamedTuple[]
        @test run_montecarlo_batch!(rows, SPEC, Earth()) === nothing
        @test run_montecarlo_per_orbit!(rows, SPEC, Earth(), 10.0, [1]) === nothing
        @test isempty(rows) && isempty(MEASUREMENTS) && isempty(WARMUPS)
    end
end
end # module
