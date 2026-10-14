using Test, CSV, DataFrames

const COMPARE_SCRIPT = joinpath(@__DIR__, "compare.jl")

function comparison_fixture(dir; payloads, samples=[1, 1], omit_states=false, codes=("pre", "post"))
    for code in codes, mode in ("serial", "predictive")
        point = joinpath(dir, code, mode)
        mkpath(point)
        touch(joinpath(point, "done"))
        raw = DataFrame(
            success=fill(true, length(samples)), case=fill("synthetic", length(samples)),
            mode=fill(mode, length(samples)), thread_count=fill(mode == "serial" ? 1 : 2, length(samples)),
            mc_samples=samples, wall_time_s=fill(1.0, length(samples)),
            sample_alloc_mb_sum=fill(2.0, length(samples)), sample_gc_time_sum_s=fill(0.0, length(samples)),
            git_commit=fill("0123456789abcdef", length(samples)),
        )
        CSV.write(joinpath(point, "parallelization_performance_raw_fixture.csv"), raw)
        if !omit_states
            states = mkpath(joinpath(point, "states"))
            for (i, payload) in enumerate(payloads(code, mode))
                write(joinpath(states, "state_synthetic_i1_s1_c$(i).bin"), payload)
            end
        end
    end
    mod = Module(gensym(:SyntheticComparison))
    Core.eval(mod, :(const ARGS = [$dir]))
    output = joinpath(dir, "comparison.log")
    open(output, "w") do io
        redirect_stdout(io) do
            Base.include(mod, COMPARE_SCRIPT)
        end
    end
    return read(output, String), CSV.read(joinpath(dir, "summary.csv"), DataFrame)
end

@testset "state identity needs complete, consistent timed dumps" begin
    for (label, payloads, samples, omit_states, complete, identical) in (
            ("identical", (_, _) -> ["A", "A"], [1, 1], false, true, true),
            ("mixed", (code, _) -> code == "pre" ? ["A", "B"] : ["C", "D"], [1, 1], false, true, false),
            ("empty", (_, _) -> String[], [1, 1], false, false, false),
            ("missing directory", (_, _) -> String[], [1, 1], true, false, false),
            ("missing repeat", (_, _) -> ["A"], [1, 1], false, false, false),
            ("extra stale dump", (_, _) -> ["A", "A", "A"], [1, 1], false, false, false),
            ("multiple samples per repeat", (_, _) -> fill("A", 5), [2, 3], false, true, true),
        )
        @testset "$label" begin
            mktempdir() do dir
                output, summary = comparison_fixture(dir; payloads, samples, omit_states)
                @test all(==(sum(samples)), summary.expected_solves)
                @test all(==(complete), summary.state_dumps_complete)
                @test count("one state across routes/threads: $identical", output) == 2
                @test occursin("one identical state across pre/post routes/threads: $identical", output)
                @test !identical || all(==(1), summary.distinct_states)
            end
        end
    end
    @testset "different valid states do not match across codes" begin
        mktempdir() do dir
            output, summary = comparison_fixture(dir; payloads=(code, _) -> fill(code, 2))
            @test all(summary.state_dumps_complete)
            @test count("one state across routes/threads: true", output) == 2
            @test occursin("one identical state across pre/post routes/threads: false", output)
        end
    end
    @testset "matching multivalued sets cannot establish cross-code identity" begin
        for swapped in (false, true)
            mktempdir() do dir
                payloads = (code, mode) -> fill(swapped && code == "post" ?
                    (mode == "serial" ? "predictive" : "serial") : mode, 2)
                output, summary = comparison_fixture(dir; payloads)
                @test all(summary.state_dumps_complete)
                @test all(==(1), summary.distinct_states)  # each point is internally consistent
                @test count("one state across routes/threads: false", output) == 2
                @test occursin("one identical state across pre/post routes/threads: false", output)
            end
        end
    end
    @testset "unrelated tags are not a pre/post comparison" begin
        mktempdir() do dir
            output, _ = comparison_fixture(dir; payloads=(_, _) -> ["A", "A"], codes=("alpha", "beta"))
            @test !occursin("one identical state across pre/post routes/threads", output)
        end
    end
end
