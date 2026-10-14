module BenchmarkRolloutGateTests
using Test
using CSV
using DataFrames
using Dates
using LinearAlgebra

const RUNTIME = normpath(joinpath(@__DIR__, "..", "..", "..", "benchmarks", "studies", "performance_runtime_analysis"))

# Select actual definitions without main.jl's unconditional native imports or
# the catalog's process-policy state. Decisions, thresholds and writers below
# are production code; only warmup/simulation and their configuration are fake.
function load_definition(file::String, wanted::Symbol)
    for expr in Meta.parseall(read(file, String)).args
        expr isa Expr || continue
        node = expr.head == :macrocall ? expr.args[end] : expr
        node isa Expr || continue
        name = if node.head == :struct
            node.args[2]
        elseif node.head == :function
            signature = node.args[1]
            signature = signature.head == :(::) ? signature.args[1] : signature
            signature.args[1]
        else
            nothing
        end
        if name == wanted
            Core.eval(@__MODULE__, expr)
            return
        end
    end
    error("Required gate definition $wanted not found in $file")
end
struct SimulationConfiguration
    mission_configuration::NamedTuple
    dynamics_model::NamedTuple
end
const HELPERS = joinpath(RUNTIME, "case_catalog", "profile_types_and_helpers.jl")
for name in (:ProfileSpec, :BenchmarkCase, :_parse_bool_env, :_parse_positive_int_env,
             :_parse_nonnegative_float_env, :_perf_case_name_filter,
             :_exclude_entry_scenarios, :_perf_smoke_mode,
             :_multirate_env_float_or_missing, :_multirate_env_int_or_missing,
             :_multirate_rollout_setting_snapshot, :_split_rollout_solver_variants)
    load_definition(HELPERS, name)
end
for kind in ("split", "multirate"), suffix in ("enforce", "case_names", "max_slowdown_ratio",
        "pos_rel_tol", "vel_rel_tol", "q_angle_tol_rad", "omega_rel_tol", "sample_count")
    load_definition(HELPERS, Symbol("_", kind, "_rollout_", suffix))
end
load_definition(joinpath(RUNTIME, "measurement", "montecarlo.jl"), :selected_cases)
for name in (:_case_with_solver, :_solution_state_at, :_relative_vector_delta, :_quaternion_angle_delta_rad)
    load_definition(joinpath(RUNTIME, "reporting", "summaries.jl"), name)
end
load_definition(joinpath(RUNTIME, "reporting", "plots_and_reporting.jl"), :_fmt)
for name in (:_trajectory_delta_metrics, :_write_split_rollout_gate_report,
             :_write_multirate_rollout_gate_report, :evaluate_split_rollout_gate,
             :evaluate_multirate_rollout_gate)
    load_definition(joinpath(RUNTIME, "reporting", "rollout_gates.jl"), name)
end
const SPEC = ProfileSpec(name="quick", repeats=1, warmup=0, max_attempts=1,
    mission_short_s=120.0, mission_long_s=600.0, montecarlo_samples=0,
    montecarlo_mission_s=120.0)
function fixture_case(name; orientation=false, quick=true)
    args = SimulationConfiguration((orientation_sim=orientation,), (spacecraft=[nothing],))
    BenchmarkCase(name=name, category="fixture", description="No simulation",
        args_template=args, run_in_quick=quick)
end
const CASES = [fixture_case("good"), fixture_case("boundary"; orientation=true),
    fixture_case("slow"), fixture_case("trajectory_failed"), fixture_case("solve_failed"),
    fixture_case("zero_baseline"), fixture_case("full_only"; quick=false)]
struct FixtureSolution
    t::Vector{Float64}
    state::NamedTuple
end
(sol::FixtureSolution)(t) = sol.state
function fixture_solution(delta)
    state = (sc=[(pos=[1.0 + delta, 0.0, 0.0], vel=[1.0 + delta, 0.0, 0.0],
                  q=[1.0, 0.0, 0.0, 0.0], ω=[1.0, 0.0, 0.0])],)
    FixtureSolution([0.0, 1.0], state)
end
const WARMUPS = NamedTuple[]
const RUNS = NamedTuple[]
function run_warmup(case::BenchmarkCase, count::Int, profile::String)
    push!(WARMUPS, (scenario=case.name, count=count, profile=profile,
                   mode=case.solver_mode_override, solver=case.split_imex_solver_override))
    nothing
end
function _run_split_gate_solution(case::BenchmarkCase, profile::String;
                                 solver_mode::String, split_solver::Union{Nothing,String}=nothing)
    push!(RUNS, (scenario=case.name, mode=solver_mode, solver=split_solver))
    baseline = solver_mode == "auto_stiff"
    success = baseline || case.name != "solve_failed"
    elapsed = baseline ? (case.name == "zero_baseline" ? 0.0 : 1.0) :
        (case.name == "slow" ? nextfloat(1.25) : 1.25)
    delta = baseline ? 0.0 : (case.name == "boundary" ? 0.125 :
                            case.name == "trajectory_failed" ? 0.25 : 0.0)
    return (ok=success, elapsed_s=elapsed, success=success,
        retcode=success ? "Success" : "FixtureFailure", solver_mode=solver_mode,
        solver_sequence=something(split_solver, solver_mode),
        solution=success ? fixture_solution(delta) : nothing,
        error_text=success ? missing : "synthetic solve failure")
end

function run_gate(kind::Symbol, requested::Vector{String}; enforce=false,
                  cases=CASES, extra_env=Pair{String,String}[])
    prefix = "SPACEAGORA_PERF_$(uppercase(string(kind)))"
    settings = ["$(prefix)_ROLLOUT_CASES" => join(requested, ","),
        "$(prefix)_ROLLOUT_ENFORCE" => string(enforce),
        "$(prefix)_MAX_SLOWDOWN_RATIO" => "1.25", "$(prefix)_POS_REL_TOL" => "0.125",
        "$(prefix)_VEL_REL_TOL" => "0.125", "$(prefix)_Q_ANGLE_TOL_RAD" => "0",
        "$(prefix)_OMEGA_REL_TOL" => "0", "$(prefix)_TRAJ_SAMPLES" => "3",
        "SPACEAGORA_PERF_SPLIT_ROLLOUT_SOLVERS" => "kencarp58,kencarp4",
        "SPACEAGORA_PERF_CASES" => "", "SPACEAGORA_PERF_EXCLUDE_ENTRY_SCENARIOS" => "0",
        "SPACEAGORA_PERF_SMOKE" => "0", "SPACEAGORA_MULTIRATE_SLOW_SOLVER" => "tsit5",
        "SPACEAGORA_MULTIRATE_FAST_SOLVER" => "auto_stiff",
        "SPACEAGORA_MULTIRATE_SLOW_DT_S" => "0.5", "SPACEAGORA_MULTIRATE_FAST_SUBSTEPS" => "4"]
    settings = merge(Dict(settings), Dict(extra_env))
    empty!(RUNS); empty!(WARMUPS)
    return withenv(settings...) do
        mktempdir() do dir
            result = nothing
            caught = try
                result = (kind == :split ? evaluate_split_rollout_gate : evaluate_multirate_rollout_gate)(SPEC, cases, dir)
                nothing
            catch err
                err
            end
            # Read the artifacts even after an exception: enforcement must occur
            # after both writes, and those files must retain the failed rows.
            files = readdir(dir; join=true)
            csv_paths = filter(p -> endswith(p, ".csv"), files)
            report_paths = filter(p -> endswith(p, ".md"), files)
            @test length(csv_paths) == 1
            @test length(report_paths) == 1
            csv_path = only(csv_paths)
            report_path = only(report_paths)
            csv = isempty(strip(read(csv_path, String))) ? DataFrame() : CSV.read(csv_path, DataFrame)
            report = read(report_path, String)
            if result !== nothing
                @test result.csv_path == csv_path
                @test result.report_path == report_path
                @test isequal(csv, result.df)
            end
            return (error=caught, df=csv, report=report, result=result,
                    warmups=copy(WARMUPS), runs=copy(RUNS))
        end
    end
end

for kind in (:split, :multirate)
    per_case = kind == :split ? 2 : 1
    @testset "$kind rollout gate decisions and evidence" begin
        passing = run_gate(kind, ["boundary", "good"]; enforce=true)
        @test passing.error === nothing
        @test all(passing.df.pass_all)
        @test passing.df.scenario == repeat(["boundary", "good"]; inner=per_case)
        @test occursin("Gate pass count: `$(2per_case)/$(2per_case)`", passing.report)
        @test all(passing.df.runtime_ratio .== 1.25)
        @test all(passing.df.compared_samples .== 3)
        @test all(passing.df.pos_rel_max[1:per_case] .== 0.125)
        @test all(passing.df.vel_rel_max[1:per_case] .== 0.125)
        @test all(passing.df.pass_q) && all(passing.df.pass_omega)
        @test length(passing.runs) == 4per_case
        @test all(row.count == 1 && row.profile == "quick" for row in passing.warmups)
        if kind == :split
            @test passing.df.split_solver == ["kencarp58", "kencarp4", "kencarp58", "kencarp4"]
            @test length(passing.warmups) == 6
        else
            @test length(passing.warmups) == 4
            @test all(passing.df.multirate_slow_dt_s .== 0.5)
            @test all(passing.df.multirate_fast_substeps .== 4)
        end
        # The previous representable tolerance excludes the exact boundary;
        # the next representable runtime above the ceiling also fails.
        prefix = "SPACEAGORA_PERF_$(uppercase(string(kind)))"
        below = run_gate(kind, ["boundary"]; extra_env=["$(prefix)_POS_REL_TOL" => string(prevfloat(0.125))])
        @test below.error === nothing
        @test all(.!below.df.pass_pos)
        @test all(below.df.pass_vel)
        @test all(.!below.df.pass_all)
        for (names, pass_count) in ((["slow", "good", "solve_failed"], per_case),
                                   (["trajectory_failed", "slow"], 0))
            observed = run_gate(kind, names; enforce=false)
            @test observed.error === nothing
            @test count(observed.df.pass_all) == pass_count
            @test occursin("Gate pass count: `$pass_count/$(length(names) * per_case)`", observed.report)
            failed_names = names == ["slow", "good", "solve_failed"] ? ["slow", "solve_failed"] : names
            expected_failures = kind == :split ?
                ["$name:$solver" for name in failed_names for solver in ("kencarp58", "kencarp4")] : failed_names
            enforced = run_gate(kind, names; enforce=true)
            @test enforced.error isa ErrorException
            expected = "$(kind == :split ? "Split" : "Multirate") rollout gate failed for $(length(expected_failures)) configuration(s): $(join(expected_failures, ", "))"
            @test enforced.error isa ErrorException && sprint(showerror, enforced.error) == expected
            @test enforced.result === nothing
            @test isequal(enforced.df, observed.df)
            @test occursin("Gate pass count: `$pass_count/$(length(names) * per_case)`", enforced.report)
            if "solve_failed" in names
                rows = observed.df[observed.df.scenario .== "solve_failed", :]
                @test all(ismissing, rows.compared_samples)
                @test all(.!rows.pass_runtime)
                @test all(.!rows.pass_all)
                @test occursin("FixtureFailure", enforced.report)
            end
        end
        zero = run_gate(kind, ["zero_baseline"])
        @test zero.error === nothing
        @test all(isinf, zero.df.runtime_ratio)
        @test all(.!zero.df.pass_runtime)
        @test all(.!zero.df.pass_all)
    end
    @testset "$kind empty and filtered selections" begin
        empty_selection = run_gate(kind, String[]; enforce=true)
        @test empty_selection.error === nothing
        @test isempty(empty_selection.df)
        @test isempty(empty_selection.runs) && isempty(empty_selection.warmups)
        @test occursin("Gate pass count: `0/0`", empty_selection.report)
        for name in ("unknown", "full_only")
            result = @test_logs (:warn, Regex("requested scenario '$name'.*skipping")) run_gate(kind, [name]; enforce=true)
            @test result.error === nothing
            @test isempty(result.df)
            @test isempty(result.runs) && isempty(result.warmups)
            @test occursin("Gate pass count: `0/0`", result.report)
        end
        result = @test_logs (:warn, r"requested scenario 'good'.*skipping") run_gate(kind, ["good"]; enforce=true, cases=BenchmarkCase[])
        @test result.error === nothing
        @test isempty(result.df)
        @test isempty(result.runs) && isempty(result.warmups)
    end
end

@testset "rollout fixture avoids simulation and native imports" begin
    for name in (:SpaceAGORA, :SimulationModel, :GRAMSuite, :SPICE, :Plots)
        @test !isdefined(@__MODULE__, name)
    end
end
end # module
