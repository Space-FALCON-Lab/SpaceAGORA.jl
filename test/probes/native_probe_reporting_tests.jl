# Native-free tests for the reporting path used by the actual coverage child.
using Test
if !isdefined(@__MODULE__, :NativeProbeReporting)
    include(joinpath(@__DIR__, "..", "helpers", "native_probe_reporting.jl"))
end

@testset "Native GRAM probe reporting" begin
    reporter = NativeProbeReporting
    @test !reporter.native_probe_required("0")
    @test reporter.native_probe_required("1")
    @test_throws ArgumentError reporter.native_probe_required("invalid")

    completed, skipped = reporter.COMPLETED * "\n", reporter.SKIPPED * "\n"
    for (output, child_ok, required, status, accepted) in (
        (completed, true, true, :completed, true),
        (completed, true, false, :completed, true),
        (skipped, true, false, :skipped, true),
        (skipped, true, true, :skipped, false),
        (completed, false, true, :completed, false),
        (completed, false, false, :completed, false),
        ("", true, true, :missing, false),
        ("", true, false, :missing, false),
        (completed * skipped, true, true, :invalid, false),
        (completed * completed, true, true, :invalid, false),
        (reporter.PREFIX * "unknown\n", true, true, :invalid, false),
    )
        io = IOBuffer()
        @test reporter.report_native_probes(io, output, child_ok; required) == accepted
        @test occursin(reporter.PREFIX * string(status), String(take!(io)))
    end

    io = IOBuffer()
    runs = Ref(0)
    callback() = (runs[] += 1)
    @test reporter.run_native_probes(callback; available=true, required=true, io) === :completed
    @test runs[] == 1
    @test String(take!(io)) == completed
    @test reporter.run_native_probes(callback; available=false, required=false, io) === :skipped
    @test runs[] == 1
    @test String(take!(io)) == skipped
    @test_throws ErrorException reporter.run_native_probes(callback; available=false, required=true, io)
    @test runs[] == 1
    @test String(take!(io)) == skipped
    @test_throws ErrorException reporter.run_native_probes(() -> error("child failed"); available=true, required=true, io)
    @test isempty(String(take!(io))) # A failed native test cannot emit completion.
end
