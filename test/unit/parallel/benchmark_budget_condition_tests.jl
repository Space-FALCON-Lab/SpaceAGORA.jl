module BenchmarkBudgetConditionTests
using Test
using SpaceAGORA

const STUDIES = normpath(joinpath(@__DIR__, "..", "..", "..", "benchmarks", "studies"))
for file in ("cli.jl", "modes.jl", "cases.jl", "reporting.jl", "execution.jl")
    include(joinpath(STUDIES, "parallelization_performance", file))
end
for file in ("cli.jl", "main.jl", "reporting.jl")
    include(joinpath(STUDIES, "paper_parallelization_benchmarks", file))
end

cmd_env(cmd) = Dict(k => v for (k, v) in (split(pair, '='; limit=2) for pair in cmd.env))
worker_command(mode; threads=4) = ppc_worker_cmd(PPCConfig(process_workers=4);
    case="independent_1sat_1hr", mode=mode, threads=threads, repeat=1,
    seed=1, mc_samples=2, outfile="unused.csv", parity=false)

# Exercise the real result writer with synthetic samples. No integration,
# campaign, worker process, native model or plot is needed for provenance.
function ppc_run_sample_batch(case::PPCCaseSpec, cfg::PPCConfig, mode::PPCModeSpec, n::Int)
    sample = (success=true, wall_time_s=1.0, retcode=:Success,
              terminal=(terminal_time_s=1.0, pos_norm_m=2.0, vel_norm_mps=3.0, mass_kg=4.0))
    return (results=fill(sample, n), batch_wall_time_s=2.0, actual_backend="none",
            execution_scope="fixture", outer_tasks=1, mode=mode,
            env_extra=Pair{String,String}[], adaptive_allocation="", policy=(;))
end

ppc_single_config(case_name::String, cfg::PPCConfig; seed::Int=cfg.seed, mc_index::Int=1) = (; seed)
ppc_solve_once(args::NamedTuple, cfg::PPCConfig) = (solution=(retcode=:Success,),)
ppc_solve_success(sol::NamedTuple) = sol.retcode === :Success
ppc_compare_trajectories(ref, cmp, args; sample_count::Int=128) = (pass=true, samples=sample_count)

@testset "benchmark budget provenance and resume" begin
    try
        withenv("SPACEAGORA_CORE_BUDGET" => "32",
                "SPACEAGORA_PERF_HARDWARE_CLASS" => "large",
                "SPACEAGORA_PPC_BUDGET_CONDITION" => nothing,
                "SPACEAGORA_PPC_EQUAL_CORE_BUDGET" => nothing) do
            SpaceAGORA.ParallelProfiles.refresh_machine_topology!()
            uncapped = ppc_budget_condition()
            off_cmd = worker_command("predictive")
            @test cmd_env(off_cmd)["SPACEAGORA_CORE_BUDGET"] == "32"
            @test cmd_env(off_cmd)["SPACEAGORA_PERF_HARDWARE_CLASS"] == "large"
            @test cmd_env(off_cmd)["SPACEAGORA_PPC_BUDGET_CONDITION"] == uncapped
            withenv("SPACEAGORA_PPC_EQUAL_CORE_BUDGET" => "0") do
                @test ppc_budget_condition() == uncapped
                @test worker_command("predictive").exec == off_cmd.exec
            end
            withenv("SPACEAGORA_PPC_EQUAL_CORE_BUDGET" => "1") do
                capped = ppc_budget_condition()
                @test capped != uncapped
                on_cmd = worker_command("predictive")
                @test on_cmd.exec == off_cmd.exec
                @test cmd_env(on_cmd)["SPACEAGORA_CORE_BUDGET"] == "4"
                @test cmd_env(on_cmd)["SPACEAGORA_PERF_HARDWARE_CLASS"] == "large"
                @test cmd_env(on_cmd)["SPACEAGORA_PPC_BUDGET_CONDITION"] == capped
                for mode in ("serial", "outer_threads", "force_none@w1+l0+b4", "rhs_serial")
                    env = cmd_env(worker_command(mode))
                    @test env["SPACEAGORA_CORE_BUDGET"] == "32"
                    @test env["SPACEAGORA_PPC_BUDGET_CONDITION"] == capped
                end
                @test cmd_env(worker_command("predictive"; threads=2))["SPACEAGORA_CORE_BUDGET"] == "32"
                @test ppc_budget_condition(; cpu_pinning=[0, 1, 2, 3]) != capped

                mktempdir() do dir
                    file = joinpath(dir, "phase", "worker_rows", "perf_fixture.csv")
                    mkpath(dirname(file))
                    withenv("SPACEAGORA_CORE_BUDGET" => "4",
                            "SPACEAGORA_PPC_BUDGET_CONDITION" => capped) do
                        SpaceAGORA.ParallelProfiles.refresh_machine_topology!()
                        cfg = PPCConfig(worker_case="independent_1sat_1hr", worker_mode="predictive",
                                        worker_threads=4, process_workers=4, worker_mc_samples=2,
                                        worker_outfile=file)
                        ppc_run_worker_performance(cfg)
                        df = CSV.read(file, DataFrame)
                        @test df.budget_condition == [capped]
                        @test df.core_budget == [4]
                        @test df.hardware_class == ["large"]
                        @test occursin("SPACEAGORA_CORE_BUDGET=4", only(df.effective_env))
                        @test occursin("SPACEAGORA_PERF_HARDWARE_CLASS=large", only(df.effective_env))
                        @test ppc_hardware_snapshot().core_budget == 4
                        @test ppc_hardware_snapshot().budget_condition == capped
                        parity_file = joinpath(dirname(file), "parity_fixture.csv")
                        ppc_run_worker_parity(_ppc_with(cfg; worker_outfile=parity_file))
                        parity = CSV.read(parity_file, DataFrame)
                        @test parity.pass == [true]
                        @test parity.budget_condition == [capped]
                        @test parity.core_budget == [4]
                        @test _ppc_worker_already_done(parity_file; budget_condition=capped)
                        @test_throws ArgumentError _ppc_worker_already_done(parity_file; budget_condition=uncapped)
                        rm(parity_file)
                    end
                    SpaceAGORA.ParallelProfiles.refresh_machine_topology!()
                    @test _ppc_worker_already_done(file; budget_condition=capped)
                    @test ppc_validate_resume_budget(dir, capped) === nothing
                    before = read(file)
                    @test_throws ArgumentError _ppc_worker_already_done(file; budget_condition=uncapped)
                    @test_throws ArgumentError ppc_validate_resume_budget(dir, uncapped)
                    @test read(file) == before

                    # A different ambient cap or pinned machine class is also a
                    # different experiment, even with equal-core still enabled.
                    withenv("SPACEAGORA_CORE_BUDGET" => "16") do
                        SpaceAGORA.ParallelProfiles.refresh_machine_topology!()
                        @test_throws ArgumentError ppc_validate_resume_budget(dir, ppc_budget_condition())
                    end
                    SpaceAGORA.ParallelProfiles.refresh_machine_topology!()
                    withenv("SPACEAGORA_PERF_HARDWARE_CLASS" => "medium") do
                        @test_throws ArgumentError ppc_validate_resume_budget(dir, ppc_budget_condition())
                    end
                    good = CSV.read(file, DataFrame)
                    invalid = copy(good); invalid.core_budget .= 32
                    CSV.write(file, invalid)
                    @test_throws ArgumentError _ppc_worker_already_done(file; budget_condition=capped)
                    invalid = copy(good); invalid.hardware_class .= "small"
                    CSV.write(file, invalid)
                    @test_throws ArgumentError _ppc_worker_already_done(file; budget_condition=capped)
                    partial = vcat(good, good)
                    partial.success[2] = false
                    CSV.write(file, partial)
                    @test !_ppc_worker_already_done(file; budget_condition=capped)
                    @test_throws ArgumentError ppc_validate_resume_budget(dir, uncapped)
                    # A mismatched completed point in a later phase is found
                    # before any point of the resumed run starts.
                    later = joinpath(dir, "later", "worker_rows", "perf_later.csv")
                    mkpath(dirname(later))
                    other = copy(good); other.budget_condition .= uncapped
                    CSV.write(later, other)
                    @test_throws ArgumentError ppc_validate_resume_budget(dir, capped)
                    rm(later)
                    raw_file = joinpath(dir, "parallelization_performance_raw_old.csv")
                    CSV.write(raw_file, other)
                    @test_throws ArgumentError ppc_validate_resume_budget(dir, capped)
                    rm(raw_file)
                    CSV.write(file, DataFrame(success=[true]))
                    err = try
                        _ppc_worker_already_done(file; budget_condition=uncapped)
                        nothing
                    catch e
                        e
                    end
                    @test err isa ArgumentError
                    @test occursin("fresh output directory", sprint(showerror, err))
                    @test_throws ArgumentError _ppc_worker_already_done(file; budget_condition=capped)
                    CSV.write(file, DataFrame(success=[false]))
                    @test !_ppc_worker_already_done(file; budget_condition=capped)
                    CSV.write(file, DataFrame(success=Bool[]))
                    @test !_ppc_worker_already_done(file; budget_condition=capped)
                    @test !_ppc_worker_already_done(joinpath(dir, "missing.csv"); budget_condition=capped)
                end
            end
        end
    finally
        SpaceAGORA.ParallelProfiles.refresh_machine_topology!()
    end
end

function budget_rows()
    rows = NamedTuple[]
    for (condition, serial, static, adaptive) in (("uncapped", 100.0, 40.0, 50.0),
                                                 ("equal_core", 20.0, 10.0, 5.0)),
        (mode, elapsed) in (("serial", serial), ("outer_threads", static), ("predictive", adaptive)),
        repeat in 1:2
        push!(rows, (phase_id="B9", case="independent_1sat_1hr", family="fixture", mode=mode,
            thread_count=4, process_workers=4, mc_samples=8, success=true,
            wall_time_s=elapsed, throughput_samples_per_s=8 / elapsed,
            execution_scope="fixture", outer_backend_actual="threads",
            budget_condition=condition,
            core_budget=(condition == "equal_core" && mode == "predictive" ? 4 : 32),
            hardware_class="large"))
    end
    return DataFrame(rows)
end

@testset "aggregation keeps budget conditions and their own baselines" begin
    raw = budget_rows()
    for summarize in (df -> ppc_summarize(df, DataFrame()), _ppb_aggregate)
        result = summarize(raw)
        @test nrow(result) == 6
        @test Set(result.budget_condition) == Set(["uncapped", "equal_core"])
        @test :core_budget in propertynames(result)
        @test :hardware_class in propertynames(result)
        speedup = summarize === _ppb_aggregate ? :speedup : :speedup_vs_serial
        for (condition, expected) in (("uncapped", 2.0), ("equal_core", 4.0))
            row = only(eachrow(result[(result.mode .== "predictive") .&
                                      (result.budget_condition .== condition), :]))
            @test row[speedup] == expected
        end
        legacy = select(raw[raw.budget_condition .== "uncapped", :],
                        Not([:budget_condition, :core_budget, :hardware_class]))
        @test all(summarize(legacy).budget_condition .== "legacy_unrecorded")
        mixed = vcat(legacy, raw[raw.budget_condition .== "equal_core", :]; cols=:union)
        @test nrow(summarize(mixed)) == 6
        @test Set(summarize(mixed).budget_condition) == Set(["legacy_unrecorded", "equal_core"])
    end
    agg = _ppb_aggregate(raw)
    @test only(agg[(agg.mode .== "predictive") .& (agg.budget_condition .== "uncapped"), :regret_vs_best_static]) == 0.25
    @test only(agg[(agg.mode .== "predictive") .& (agg.budget_condition .== "equal_core"), :regret_vs_best_static]) == -0.5
    replace!(agg.mode, "predictive" => "full_smart")
    report = _ppb_router_regret_summary(agg)
    @test nrow(report) == 4
    @test Set(report.budget_condition) == Set(["uncapped", "equal_core"])
end
end # module
