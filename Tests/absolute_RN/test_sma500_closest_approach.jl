import Pkg
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(ROOT, "9_SpaceAGORA.jl-main_full integration"))
using Test, DataFrames, LinearAlgebra, CSV, Arrow, TOML, Printf, Statistics
include(joinpath(ROOT, "Simulation_absolute_RN", "Code", "write_sma500_closest_approach.jl"))
const CA = SMA500ClosestApproach

fresh() = DataFrame(closest_time_s=[NaN], closest_separation_m=[Inf],
    closest_relative_speed_mps=[NaN], closest_range_rate_mps=[NaN])

@testset "Speed at closest approach, not minimum speed" begin
    # A curved pass whose speed minimum and separation minimum occur at different times.
    function at(t)
        x = exp(t - 2) - 1
        r = [x, 1 + x^2, 0]
        v = exp(t - 2) .* [1, 2x, 0]
        (distance=norm(r), speed=norm(v), rate=dot(r, v) / norm(r))
    end
    rows = fresh()
    CA.scan_interval!(rows, 1, 0.0, 4.0, at)
    @test rows.closest_time_s[1] ≈ 2.0 atol=1e-7
    @test rows.closest_separation_m[1] ≈ 1.0
    @test rows.closest_relative_speed_mps[1] ≈ 1.0 atol=1e-7
    @test at(0.0).speed < rows.closest_relative_speed_mps[1]
    @test abs(rows.closest_range_rate_mps[1]) < 1e-7

    # A shrinking-radius spiral: speed increases while separation decreases.
    spiral(t) = (distance=10 - t, speed=1 + t, rate=-1.0)
    rows = fresh()
    CA.scan_interval!(rows, 1, 0.0, 4.0, spiral)
    @test rows.closest_time_s == [4.0]
    @test rows.closest_relative_speed_mps == [5.0]
end

@testset "Clipped, constant-range, and step-boundary minima" begin
    rows = fresh()
    CA.scan_interval!(rows, 1, 2.0, 3.0,
        t -> (distance=t, speed=5.0, rate=1.0))
    @test rows.closest_time_s == [2.0]
    @test rows.closest_relative_speed_mps == [5.0]
    rows = fresh()
    constant(t) = (distance=10.0, speed=2.0, rate=0.0)
    CA.scan_interval!(rows, 1, 2.0, 3.0, constant)
    CA.scan_interval!(rows, 1, 0.0, 2.0, constant)
    @test rows.closest_time_s == [0.0]
    rows = fresh()
    at(t) = (distance=1 + (t - 2)^2, speed=3 + t, rate=2(t - 2))
    CA.scan_interval!(rows, 1, 0.0, 2.0, at)
    CA.scan_interval!(rows, 1, 2.0, 4.0, at)
    @test rows.closest_time_s == [2.0]
    @test rows.closest_relative_speed_mps == [5.0]
    @test_throws ErrorException CA.scan_interval!(rows, 1, 3.0, 2.0, at)
end

@testset "Global minimum across repeated approaches" begin
    rows = fresh()
    for (left, right, center, offset) in ((0.0, 2.0, 1.0, 10.0),
            (2.0, 4.0, 3.0, 5.0), (4.0, 6.0, 5.0, 8.0))
        at(t) = (distance=offset + (t - center)^2, speed=center,
            rate=2(t - center))
        CA.scan_interval!(rows, 1, left, right, at)
    end
    @test rows.closest_time_s[1] ≈ 3.0
    @test rows.closest_separation_m[1] ≈ 5.0
    @test rows.closest_relative_speed_mps[1] ≈ 3.0
end

@testset "Relative state quantities" begin
    state = (sc=[(pos=[0., 0., 0.], vel=[1., 0., 0.]),
        (pos=[3., 4., 0.], vel=[1., 0., 2.])],)
    @test CA.quantities(state, 1, 2) == (speed=2.0, rate=0.0, distance=5.0)
    @test_throws ErrorException CA.quantities(state, 1, 1)
end

struct TestIntegrator{F,S}
    t::Float64
    uprev::S
    u::S
    at::F
end
(i::TestIntegrator)(t) = i.at(t)

@testset "Observer checks dense boundaries and step coverage" begin
    state(t) = (sc=[(pos=[0., 0., 0.], vel=[0., 0., 0.]),
        (pos=[t - 2, 1., 0.], vel=[1., 0., 0.])],)
    geometry = DataFrame(run_id=["test"], encounter_id=[1], sc_a=[1], sc_b=[2],
        start_time_s=[0.25], end_time_s=[3.75], start_clipped=[false], end_clipped=[false])
    for (boundary, time) in (("start", 0.25), ("end", 3.75)),
            (endpoint, index) in (("a", 1), ("b", 2)),
            (field, property) in (("r", :pos), ("v", :vel)), (axis, k) in zip("xyz", 1:3)
        geometry[!, Symbol("$(boundary)_sc_$(endpoint)_$(field)_$axis")] =
            [getproperty(state(time).sc[index], property)[k]]
    end

    if "--campaign" in ARGS
        @testset "Nine-case closest-approach report and preserved parameters" begin
            output = CA.SC.DEFAULT_OUTPUT
            original_events = CSV.read(joinpath(output, "encounter_evaluation",
                "target_helper_encounters.csv"), DataFrame)
            original_laser = CSV.read(joinpath(output, "encounter_evaluation",
                "target_helper_laser_intervals.csv"), DataFrame)
            records = CSV.read(joinpath(CA.OUTPUT, "encounters.csv"), DataFrame)
            text = read(CA.REPORT, String)
            table_rows = filter(line -> occursin(r"^\| (0\.0|0\.5|1\.0) \|", line),
                split(text, '\n'))
            @test length(table_rows) == 12
            @test occursin("Encounter speed at closest approach", text)
            @test occursin("Signed range rate at entry", text)
            specs = ((records, :closest_relative_speed_mps),
                (original_events, :entry_range_rate_mps), (original_events, :duration_s),
                (original_laser, :duration_s))
            for (page, (table, field)) in enumerate(specs),
                    (i, inclination) in enumerate(CA.INCLINATIONS),
                    (h, altitude) in enumerate(CA.ALTITUDES)
                c = CA.StudyCase(altitude, inclination, 500, :gve_sma)
                selected = filter(r -> r.run_id == CA.SC.case_signature(c), table)
                @test !isempty(selected)
                @test all(==(500), selected.helper_count)
                @test all(==("gve_sma"), selected.scheduler)
                values = selected[!, field]
                expected = [@sprintf("%.6g", v) for v in (minimum(values), mean(values), maximum(values))]
                cells = strip.(split(table_rows[3(page - 1) + i], '|')[2:end-1])
                @test cells[(3h - 1):(3h + 1)] == expected
            end
            for inclination in CA.INCLINATIONS, altitude in CA.ALTITUDES
                c = CA.StudyCase(altitude, inclination, 500, :gve_sma)
                id = CA.SC.case_signature(c)
                result = filter(r -> r.run_id == id, records)
                source = filter(r -> r.run_id == id, original_events)
                @test Set(result.encounter_id) == Set(source.encounter_id)
                @test allunique(result.encounter_id)
                @test all(result.start_time_s .<= result.closest_time_s .<= result.end_time_s)
                @test all(0 .< result.closest_separation_m .<= CA.SC.RANGE_M + 1e-3)
                verification = TOML.parsefile(CA.case_file(c) * ".verification.toml")
                @test verification["passed"]
                @test verification["boundary_checks"] == 2nrow(result)
                @test verification["result_sha256"] == CA.SC.filehash(CA.case_file(c))

                track = DataFrame(Arrow.Table(joinpath(CA.case_directory(output, c), "longest_encounter.feather")))
                track = sort(filter(r -> r.is_target_helper_longest, track), :time)
                encounter = only(eachrow(filter(r -> r.encounter_id == first(track.encounter_id), result)))
                time = encounter.closest_time_s
                before = searchsortedlast(track.time, time)
                after = searchsortedfirst(track.time, time)
                @test 1 <= before <= after <= nrow(track)
                track_state = index -> hcat([[track[index, Symbol("sc_$(endpoint)_$(field)_$axis")]
                    for field in ("pos", "vel") for axis in 1:3] for endpoint in ("a", "b")]...)
                interpolated = if before == after
                    track_state(before)
                else
                    duration = track.time[after] - track.time[before]
                    CA.OracleStudy.EncounterTimingRecorder.hermite_state(track_state(before), track_state(after),
                        duration, (time - track.time[before]) / duration)
                end
                speed, rate, distance = CA.encounter_quantities(interpolated[:, 1], interpolated[:, 2])
                @test speed ≈ encounter.closest_relative_speed_mps atol=1e-3
                @test distance ≈ encounter.closest_separation_m atol=1e-2
            end
        end
    end
    observer = CA.closest_observer(geometry)
    observer.observe(TestIntegrator(2.0, state(0.0), state(2.0), state), 0.0, [])
    observer.observe(TestIntegrator(4.0, state(2.0), state(4.0), state), 2.0, [])
    @test observer.rows.closest_time_s == [2.0]
    @test observer.rows.closest_relative_speed_mps == [1.0]
    @test observer.rows.closest_separation_m == [1.0]
    @test observer.checks().boundary_checks == 2
    @test observer.checks().max_position_error_m == 0
    @test observer.checks().max_velocity_error_mps == 0
    @test_throws ErrorException observer.observe(
        TestIntegrator(6.0, state(5.0), state(6.0), state), 5.0, [])
end
