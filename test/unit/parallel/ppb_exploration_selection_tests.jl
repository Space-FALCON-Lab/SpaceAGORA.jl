using Test

# Load the real catalog separately at portable and published budgets. Definitions
# only: no campaign, precompile workload, or worker is launched by this test.
function _ppb_selection_sandbox(budget::Int)
    root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    sandbox = Module(gensym(:PPBSelection))
    withenv("SPACEAGORA_PPB_PAPER_BUDGET" => string(budget)) do
        for path in (
            "parallelization_performance/cli.jl",
            "parallelization_performance/modes.jl",
            "parallelization_performance/cases.jl",
            "paper_parallelization_benchmarks/cli.jl",
            "paper_parallelization_benchmarks/reporting.jl",
            "paper_parallelization_benchmarks/main.jl",
        )
            Base.include(sandbox, joinpath(root, "benchmarks", "studies", path))
        end
    end
    return sandbox
end

@testset "route exploration requires explicit compatible selection" begin
    withenv("SPACEAGORA_PPB_MIN_REPEATS" => "", "SPACEAGORA_PPB_MIN_WARMUP" => "") do
        for budget in (8, 32)
            S = _ppb_selection_sandbox(budget)
            # Access newly loaded module bindings in their definition world.
            Base.invokelatest() do
                config(; kwargs...) = Base.invokelatest(S.PPBConfig; kwargs...)
                active(cfg) = Base.invokelatest(S._ppb_active_phases, cfg)
                # Every historical default phase remains; X phases are opt-in even
                # when the host is large enough to execute them.
                expected = [p.id for p in S.PAPER_BENCHMARK_PHASES if !startswith(p.id, "X")]
                @test [p.id for p in active(config(process_workers=budget))] == expected
                preview = active(config(preview=true, process_workers=budget))
                @test all(p -> !startswith(p.id, "X"), preview)
                @test Set(p.id for p in preview) == setdiff(Set(expected), S.PPB_PREVIEW_SKIP_PHASES)
                for id in ("X1", "X3", "X4", "X5a", "X5b")
                    cfg = config(phases=[id], threads=[32], process_workers=32)
                    if budget == 8
                        @test_throws ArgumentError active(cfg)
                    else
                        @test only(active(cfg)).id == id
                        @test_throws ArgumentError active(config(phases=[id], threads=[32],
                                                                 process_workers=32, preview=true))
                        if id == "X1"
                            @test_throws ArgumentError active(config(phases=[id], threads=[8]))
                        else
                            @test_throws ArgumentError active(config(phases=[id], process_workers=8))
                        end
                    end
                end
                @test [p.id for p in active(config(phases=["P4", "P3", "P4"]))] == ["P4", "P3"]
            end
        end
    end
end
