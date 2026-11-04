import Pkg
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(ROOT, "9_SpaceAGORA.jl-main_full integration"))
using Test, DataFrames, CSV, Arrow, LinearAlgebra
include(joinpath(ROOT, "Simulation_absolute_RN", "Code", "evaluate_encounters.jl"))
const EV = EncounterEvaluation

@testset "Distribution statistics and empty groups" begin
    stats = EV.describe([0.0, 10.0, 20.0, 30.0])
    @test stats.n == 4
    @test stats.mean == stats.median == 15
    @test stats.p05 ≈ 1.5
    @test stats.p25 == 7.5
    @test stats.p75 == 22.5
    @test stats.p95 ≈ 28.5
    @test stats.std_population ≈ sqrt(125)
    @test EV.describe([2.0]).std_population == 0
    @test EV.describe(Float64[]).n == 0
    @test ismissing(EV.describe(Float64[]).median)
    @test_throws ErrorException EV.describe([NaN])
end

@testset "Firing overlap accounting, not start encounter IDs" begin
    encounters = DataFrame(encounter_id=[10, 11, 12], sc_a=[1, 1, 2], sc_b=[3, 3, 3],
        start_time_s=[0., 8., 0.], end_time_s=[5., 12., 12.], duration_s=[5., 4., 12.])
    laser = DataFrame(laser_id=[20, 21], sc_a=[1, 1], sc_b=[3, 3],
        start_time_s=[3., 10.], end_time_s=[9., 11.], duration_s=[6., 1.])
    audit = DataFrame(laser_id=[20, 20, 20, 21], sc_a=fill(1, 4), sc_b=fill(3, 4),
        encounter_id=[10, missing, 11, 11], start_time_s=[3., 5., 8., 10.],
        end_time_s=[5., 8., 9., 11.], duration_s=[2., 3., 1., 1.],
        within_range=[true, false, true, true])
    EV.attach_firing!(encounters, laser, audit)
    @test encounters.laser_on_within_encounter_s == [2, 2, 0]
    @test laser.within_range_s == [3, 1]
    @test laser.outside_range_s == [3, 0]
    broken = copy(audit)
    broken.encounter_id[1] = 99
    @test_throws ErrorException EV.attach_firing!(encounters, laser, broken)
    @test_throws ErrorException EV.attach_firing!(encounters, laser, audit[2:end, :])
    @test_throws ErrorException EV.attach_firing!(encounters, laser, vcat(audit, audit[1:1, :]))
end

@testset "Interval validation" begin
    intervals = DataFrame(encounter_id=[1, 2], sc_a=[1, 1], sc_b=[2, 2],
        start_time_s=[0., 5.], end_time_s=[3., 10.], duration_s=[3., 5.],
        start_clipped=[true, false], end_clipped=[false, true])
    @test isnothing(EV.validate_intervals(intervals, :encounter_id, 10., 2))
    broken = copy(intervals)
    broken.duration_s[1] = 2
    @test_throws ErrorException EV.validate_intervals(broken, :encounter_id, 10., 2)
    broken = copy(intervals)
    broken.encounter_id[2] = 1
    @test_throws ErrorException EV.validate_intervals(broken, :encounter_id, 10., 2)
end

@testset "Real saved case and CSV-driven report" begin
    c = EV.StudyConfig.StudyCase(950, 0.5, 1, :gve_sma)
    e, l, row = EV.read_case(EV.StudyConfig.DEFAULT_OUTPUT, c)
    @test nrow(e) == 2
    @test nrow(l) == 2
    @test all(==(2), e.sc_b)
    @test sum(e.laser_on_within_encounter_s) ≈ sum(l.within_range_s)
    @test row.total_laser_on_s ≈ 3354.056842646578
    @test row.outside_range_laser_on_s ≈ 4.93032756325556
    raw = DataFrame(Arrow.Table(joinpath(row.source_directory, "encounters.feather")))
    for event in eachrow(e)
        boundary = only(eachrow(filter(r -> r.encounter_id == event.encounter_id, raw)))
        relative(field) = [getproperty(boundary, Symbol("start_sc_b_$(field)_$axis")) -
            getproperty(boundary, Symbol("start_sc_a_$(field)_$axis")) for axis in ("x", "y", "z")]
        r, v = relative("r"), relative("v")
        @test event.entry_relative_speed_mps ≈ norm(v)
        @test event.entry_range_rate_mps ≈ dot(r, v) / norm(r)
    end
    cases = DataFrame([row])
    # A zero-firing case remains in the denominator but supplies no firing interval.
    absent = merge(row, (run_id="zero-event-test", encounter_count=0, laser_interval_count=0))
    groups = EV.distributions(e, l, DataFrame([row, absent]))
    overall = filter(r -> r.dimension == "all" && r.population == "all_observed", groups)
    @test nrow(overall) == 7
    @test all(==(2), overall.cases)
    @test all(==(1), overall.contributing_cases)
    unclipped = filter(r -> r.dimension == "all" && r.population == "unclipped_only", groups)
    @test all(==(1), unclipped.n)
    @test all(iszero, unclipped.clipped_events)
    empty_stats = EV.distributions(e[1:0, :], l[1:0, :], cases)
    @test all(iszero, empty_stats.n)
    @test all(ismissing, empty_stats.mean)
    mktempdir() do output
        CSV.write(joinpath(output, "distributions.csv"), EV.distributions(e, l, cases))
        CSV.write(joinpath(output, "case_summary.csv"), cases)
        EV.write_report(output)
        text = read(joinpath(output, "summary.md"), String)
        @test all(occursin("## $k)", text) for k in 1:4)
        @test occursin("3354.06", text)
        stats = CSV.read(joinpath(output, "distributions.csv"), DataFrame)
        stats.mean .= 12345.0
        CSV.write(joinpath(output, "distributions.csv"), stats)
        EV.main(["--output", output, "--report-only"])
        @test occursin("12345", read(joinpath(output, "summary.md"), String))
    end
end

if "--campaign" in ARGS
    @testset "Full campaign output reconciliation" begin
        output = EV.OUTPUT
        e = CSV.read(joinpath(output, "target_helper_encounters.csv"), DataFrame)
        l = CSV.read(joinpath(output, "target_helper_laser_intervals.csv"), DataFrame)
        cases = CSV.read(joinpath(output, "case_summary.csv"), DataFrame)
        stats = CSV.read(joinpath(output, "distributions.csv"), DataFrame;
            types=Dict(:group => String))
        @test nrow(cases) == 270
        @test nrow(stats) == (1 + 3 + 5 + 9 + 2) * 2 * 7
        @test allunique(cases.run_id)
        @test all(e.sc_b .== e.helper_count .+ 1)
        @test all(l.sc_b .== l.helper_count .+ 1)
        @test all(e.min_relative_speed_mps .- 1e-8 .<= e.entry_relative_speed_mps .<=
            e.max_relative_speed_mps .+ 1e-8)
        @test all(e.min_range_rate_mps .- 1e-8 .<= e.entry_range_rate_mps .<= e.max_range_rate_mps .+ 1e-8)
        @test nrow(e) == sum(cases.encounter_count)
        @test nrow(l) == sum(cases.laser_interval_count)
        @test sum(e.laser_on_within_encounter_s) ≈ sum(l.within_range_s)
        @test sum(l.duration_s) ≈ sum(cases.total_laser_on_s)
        for (dimension, _) in EV.DIMENSIONS, (kind, field, _, _) in EV.METRICS
            selected = filter(r -> r.dimension == string(dimension) &&
                r.population == "all_observed" && r.event_kind == string(kind) &&
                r.metric == string(field), stats)
            table = kind == :encounter ? e : l
            @test sum(selected.n) == nrow(table)
            @test sum(selected.cases) == 270
            for row in eachrow(selected)
                values = table[[string(v) == row.group for v in table[!, dimension]], field]
                @test row.n == length(values)
                @test row.median ≈ EV.quantile(values, 0.5)
                @test row.p95 ≈ EV.quantile(values, 0.95)
            end
        end
        report = read(joinpath(output, "summary.md"), String)
        EV.main(["--report-only"])
        @test read(joinpath(output, "summary.md"), String) == report
    end
end
