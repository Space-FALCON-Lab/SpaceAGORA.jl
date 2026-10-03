# The coverage quality gate drops lines inside PrecompileTools
# `@setup_workload` / `@compile_workload` blocks (they run only at precompile
# time, never under `--code-coverage=user`) and nothing else.
using Test

module CoverageGateUnderTest
include(joinpath(@__DIR__, "..", "..", "gates", "ci_coverage_quality_gate.jl"))
end
const CGT = CoverageGateUnderTest

# A .cov line: a count or "-" in a 9-wide field, then the source text.
function _write_cov(path::String, source_lines::Vector{String}, counts::Dict{Int, Int})
    open(path, "w") do io
        for (i, text) in enumerate(source_lines)
            token = haskey(counts, i) ? string(counts[i]) : "-"
            println(io, lpad(token, 9), " ", text)
        end
    end
end

const FIXTURE = """
function helper_called_from_workload()
    return 1 + 1
end

@setup_workload begin
    x = 1
    @compile_workload begin
        helper_called_from_workload()
    end
    y = x
end

PrecompileTools.@compile_workload begin
    z = 3
end

@compile_workload helper_called_from_workload()

function after_workload()
    return 2
end
"""

@testset "coverage gate: precompile-only exclusion" begin
    @testset "span finder" begin
        spans = CGT.precompile_only_spans(FIXTURE; filename="fixture.jl")
        # Nested @compile_workload is inside the @setup_workload span, so only
        # the three top-most calls appear.
        @test spans == [5:11, 13:15, 17:17]
        @test isempty(CGT.precompile_only_spans("f(x) = x\n@foo begin\n  1\nend\n"))
    end

    @testset "summary excludes block lines only" begin
        mktempdir() do root
            src_dir = joinpath(root, "src")
            mkpath(src_dir)
            src_lines = split(chomp(FIXTURE), '\n') .|> String
            write(joinpath(src_dir, "pc.jl"), FIXTURE)
            # Executable lines: 2 (helper body, covered), 6, 8, 10, 14, 17
            # (workload bodies, never run under coverage), 20 (after, uncovered).
            counts = Dict(2 => 5, 6 => 0, 8 => 0, 10 => 0, 14 => 0, 17 => 0, 20 => 0)
            cov = joinpath(src_dir, "pc.jl.12345.cov")
            _write_cov(cov, src_lines, counts)

            plain = "g() = 1\nh() = 2\n"
            write(joinpath(src_dir, "plain.jl"), plain)
            _write_cov(joinpath(src_dir, "plain.jl.12345.cov"), String["g() = 1", "h() = 2"], Dict(1 => 1, 2 => 0))

            buf = IOBuffer()
            summaries = CGT.summarize_coverage([cov, joinpath(src_dir, "plain.jl.12345.cov")]; repo_root=root, io=buf)
            out = String(take!(buf))
            by_path = Dict(s.path => s for s in summaries)

            pc = by_path[joinpath("src", "pc.jl")]
            @test pc.executable == 2          # lines 2 and 20 remain
            @test pc.covered == 1             # the helper called from the workload still counts
            @test pc.percent == 50.0

            plain_summary = by_path[joinpath("src", "plain.jl")]
            @test (plain_summary.covered, plain_summary.executable) == (1, 2)

            @test occursin("precompile_only_excluded  $(joinpath("src", "pc.jl"))  lines=5 (5-11,13-15,17)", out)
            @test !occursin("precompile_only_excluded  $(joinpath("src", "plain.jl"))", out)
            @test occursin("precompile_only_excluded_total  files=1 lines=5", out)
        end
    end

    @testset "parse error fails loudly" begin
        mktempdir() do root
            src_dir = joinpath(root, "src")
            mkpath(src_dir)
            write(joinpath(src_dir, "broken.jl"), "function f(\n    1 +\n")
            cov = joinpath(src_dir, "broken.jl.1.cov")
            _write_cov(cov, String["function f(", "    1 +"], Dict(2 => 0))
            err = try
                CGT.summarize_coverage([cov]; repo_root=root, io=devnull)
                nothing
            catch e
                e
            end
            @test err isa ErrorException
            @test occursin("could not parse", err.msg)
            @test occursin("broken.jl", err.msg)
        end
    end

    @testset "current tree" begin
        # Precompile-only spans across the package's own sources: only the
        # workload block in src/precompile_workload.jl is expected.
        repo = normpath(joinpath(@__DIR__, "..", "..", ".."))
        found = Dict{String, Vector{UnitRange{Int}}}()
        for dir in ("src", "ext"), (root, _, files) in walkdir(joinpath(repo, dir)), f in files
            endswith(f, ".jl") || continue
            path = joinpath(root, f)
            spans = CGT.precompile_only_spans(read(path, String); filename=path)
            isempty(spans) || (found[relpath(path, repo)] = spans)
        end
        @test collect(keys(found)) == [joinpath("src", "precompile_workload.jl")]
        @test length(found[joinpath("src", "precompile_workload.jl")]) == 1
    end
end
