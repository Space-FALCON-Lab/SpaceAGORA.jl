using Test
include("extract_feather_encounters.jl")
using OrdinaryDiffEq
include(joinpath(@__DIR__, "..", "2_SpaceAGORA.jl", "ORACLE", "functions", "FeatherRecording.jl"))
using .FeatherRecording

function write_encounter_fixture(path; active_helpers=[0, 2, 0, 2, 0], third_spacecraft=false)
    times = collect(0.0:10.0:40.0)
    columns = Dict{Symbol,AbstractVector}(:time => times)
    for spacecraft in 1:(third_spacecraft ? 3 : 2), field in ("pos", "vel"), component in 1:3
        columns[Symbol("sc$(spacecraft)_$(field)_$(component)")] = zeros(length(times))
    end
    columns[:sc1_pos_1] = 100 .+ times
    columns[:sc1_vel_1] = ones(length(times))
    columns[:sc2_pos_1] = columns[:sc1_pos_1] + [20, 0, 20, 0, 20]
    columns[:sc2_vel_1] = columns[:sc1_vel_1] + [-2, 2, 2, -2, 2]
    if third_spacecraft
        columns[:sc3_pos_1] = columns[:sc2_pos_1] .+ 5
        columns[:sc3_vel_1] = fill(10.0, length(times))
    end
    active_helpers !== nothing && (columns[:laser_active_helper] = active_helpers)
    Arrow.write(path, columns)
end

@testset "Feather encounter extraction" begin
    times = [0., 10., 20.]
    @test contact_windows(times, [20., 0., 20.], 10.) == [(5., 15.)]
    @test contact_windows(times, [0., 0., 20.], 10.) == [(0., 15.)]
    @test contact_windows(times, [20., 0., 0.], 10.) == [(5., 20.)]
    @test contact_windows(times, [0., 0., 0.], 10.) == [(0., 20.)]
    @test isempty(contact_windows(times, [20., 20., 20.], 10.))
    @test contact_windows(times, [20., 10., 20.], 10.) == [(10., 10.)]
    @test pair_motion(zeros(6), [2., 0., 0., -3., 0., 0.]) == (2., 3., -3.)
    @test pair_motion(zeros(6), [2., 0., 0., 3., 0., 0.]) == (2., 3., 3.)
    @test ismissing(last(pair_motion(zeros(6), zeros(6))))

    mktempdir() do directory
        path = joinpath(directory, "synthetic.feather")
        write_encounter_fixture(path)
        result = extract_feather_encounters(path; maximum_range_m=10., output_dir=joinpath(directory, "CSV"))
        @test result.encounters.start_time_s == [5., 25.]
        @test result.encounters.end_time_s == [15., 35.]
        @test result.encounters.geometry_duration_s == [10., 10.]
        @test result.encounters.laser_on_estimate_s == [5., 5.]
        @test !any(result.encounters.start_clipped)
        @test !any(result.encounters.end_clipped)
        @test result.samples.time_s == [5., 10., 15., 25., 30., 35.]
        @test result.samples.encounter_id == [1, 1, 1, 2, 2, 2]
        @test result.samples.interpolated == [true, false, true, true, false, true]
        @test result.samples.geometry_range_m == [10., 0., 10., 10., 0., 10.]
        @test result.samples.sc_a_r_x == 100 .+ result.samples.time_s
        @test result.samples.sc_b_r_x == result.samples.sc_a_r_x + [10., 0., 10., 10., 0., 10.]
        @test result.extrema.value == [2., 2., -2., 2., 10., 10., 5., 5.]
        @test result.extrema.sc_a == fill(1, 8)
        @test result.extrema.sc_b == fill(2, 8)
        for (csv_path, expected) in ((result.encounters_path, result.encounters),
                                     (result.samples_path, result.samples),
                                     (result.extrema_path, result.extrema))
            @test isequal(CSV.read(csv_path, DataFrame), expected)
        end

        write_encounter_fixture(path; active_helpers=nothing)
        unknown = extract_feather_encounters(path; maximum_range_m=10., output_dir=joinpath(directory, "CSV"))
        @test all(ismissing, unknown.encounters.laser_on_estimate_s)
        @test all(ismissing, unknown.samples.laser_on_held)

        write_encounter_fixture(path; active_helpers=[0, missing, 0, 2, 0])
        partial = extract_feather_encounters(path; maximum_range_m=10., output_dir=joinpath(directory, "CSV"))
        @test ismissing(partial.encounters.laser_on_estimate_s[1])
        @test partial.encounters.laser_on_estimate_s[2] == 5.

        write_encounter_fixture(path; third_spacecraft=true)
        all_pairs = extract_feather_encounters(path; maximum_range_m=10., output_dir=joinpath(directory, "CSV"))
        helper_pair = filter(row -> row.sc_a == 2 && row.sc_b == 3, all_pairs.encounters)
        @test nrow(helper_pair) == 1
        @test helper_pair.geometry_duration_s == [40.]
        @test helper_pair.laser_on_estimate_s == [0.]
        @test only(helper_pair.start_clipped) && only(helper_pair.end_clipped)
        @test all_pairs.extrema.value[2] == 11.

        empty_path = joinpath(directory, "empty.feather")
        columns = DataFrame(Arrow.Table(path); copycols=true)
        columns.sc2_pos_1 .= columns.sc1_pos_1 .+ 100
        columns.sc3_pos_1 .= columns.sc1_pos_1 .+ 200
        Arrow.write(empty_path, columns)
        empty_result = extract_feather_encounters(empty_path; maximum_range_m=10., output_dir=joinpath(directory, "CSV"))
        @test isempty(empty_result.encounters)
        @test isempty(empty_result.samples)
        @test all(ismissing, empty_result.extrema.value[5:8])
        @test names(CSV.read(empty_result.encounters_path, DataFrame)) == names(empty_result.encounters)
        @test_throws ArgumentError extract_feather_encounters(path; maximum_range_m=-1.)
        columns.time[2] = columns.time[1]
        Arrow.write(empty_path, columns)
        @test_throws ArgumentError extract_feather_encounters(empty_path; maximum_range_m=10.)
    end
end

@testset "Fixed-cadence solver endpoints" begin
    problem = ODEProblem((state, params, time)->1., 0., (0., 21.))
    for early_stop in (false, true)
        callback = early_stop ? CallbackSet(ContinuousCallback(
            (state, time, integrator)->time-15., integrator->terminate!(integrator);
            save_positions=(false, false))) : CallbackSet()
        solution = solve(problem, Vern9(); callback, saveat=10., save_end=true,
            reltol=1e-12, abstol=1e-12)
        @test solution.t ≈ (early_stop ? [0., 10., 15.] : [0., 10., 20., 21.])
    end
end

@testset "Runtime interval bundle" begin
    @test output_times(21., 10.) == [0., 10., 20., 21.]
    @test output_times(20., 10.) == [0., 10., 20.]
    @test_throws ArgumentError output_times(0.)
    initial = zeros(6, 2)
    initial[1, 2] = -20.
    initial[4, 2] = 2.
    states(state) = reshape(state, 6, 2)
    laser(state, time) = abs(state[1, 2]) < 10 && (time < 8 || time >= 12) ? [(1, 2)] : Tuple{Int,Int}[]
    gates(state, time) = [time-8., time-12.]
    recorder = IntervalRecorder(states, laser, vec(initial), 10.; gates)
    function straight_line!(derivative, state, params, time)
        fill!(derivative, 0.)
        derivative[1] = state[4]
        derivative[7] = state[10]
    end
    solution = solve(ODEProblem(straight_line!, vec(initial), (0., 21.)), Tsit5();
        callback=CallbackSet(recording_callbacks(recorder)...),
        saveat=output_times(21., 10.), save_everystep=false, save_end=true, dtmax=2.)
    finish_recording!(recorder, solution.u[end], solution.t[end])
    @test solution.t == output_times(21., 10.)
    @test recorder.geometry.start_time_s ≈ [5.]
    @test recorder.geometry.end_time_s ≈ [15.]
    @test recorder.geometry.duration_s ≈ [10.]
    @test recorder.laser_on.start_time_s ≈ [5., 12.]
    @test recorder.laser_on.end_time_s ≈ [8., 15.]
    @test sum(recorder.laser_on.duration_s) ≈ 6.
    @test recorder.laser_on.encounter_id == [1, 1]
    @test recorder.geometry.start_sc_b_r_x ≈ [-10.]
    @test recorder.geometry.end_sc_b_r_x ≈ [10.]
    @test !any(recorder.geometry.start_clipped)
    @test !any(recorder.geometry.end_clipped)
    mktempdir() do directory
        feather_dir = joinpath(directory, "smoke", "feather", "synthetic")
        mkpath(feather_dir)
        columns = Dict{Symbol,AbstractVector}(:time=>solution.t)
        for spacecraft in 1:2, (field, offset) in (("pos", 0), ("vel", 3)), component in 1:3
            columns[Symbol("sc$(spacecraft)_$(field)_$(component)")] =
                [state[6*(spacecraft-1)+offset+component] for state in solution.u]
        end
        path = joinpath(feather_dir, "trajectory.feather")
        Arrow.write(path, columns; metadata=Dict("scenario"=>"synthetic", "maximum_range_m"=>"10.0"))
        write_intervals(recorder, feather_dir; source="synthetic", scenario="synthetic", laser_timing="continuous test gates")
        result = extract_feather_encounters(feather_dir; output_dir=joinpath(directory, "CSV"))
        @test encounter_csv_directory(feather_dir) == joinpath(directory, "smoke", "CSV", "synthetic")
        @test encounter_csv_directory(path) == joinpath(directory, "smoke", "CSV", "synthetic")
        for project in ("prototype", "spaceagora")
            root = joinpath(directory, project, "output")
            @test encounter_csv_directory(joinpath(root, "feather", "synthetic", "trajectory.feather")) ==
                joinpath(root, "CSV", "synthetic")
        end
        @test result.encounters.laser_on_s ≈ [6.]
        @test result.encounters.geometry_duration_s ≈ [10.]
        @test result.samples.time_s ≈ [5., 10., 15.]
        @test result.samples.sc_b_r_x ≈ [-10., 0., 10.]
        @test !any(result.samples.interpolated)
        @test !any(result.samples.laser_on[2:3])
        @test result.extrema.metric[end] == "laser_on_s"
        for name in ("trajectory", "geometry_encounters", "laser_on", "encounters", "encounter_states", "extrema")
            folder = name in ("trajectory", "geometry_encounters", "laser_on") ? "sim_output" : "analysis"
            @test isfile(joinpath(directory, "CSV", folder, name*".csv"))
            @test !isfile(joinpath(directory, "CSV", name*".csv"))
        end
        @test dirname(result.extrema_path) == dirname(result.encounters_path) ==
            dirname(result.samples_path) == joinpath(directory, "CSV", "analysis")
        @test_throws ArgumentError extract_feather_encounters(feather_dir; maximum_range_m=11.)
    end
end

@testset "Simultaneous pair crossings" begin
    initial = zeros(6, 3)
    initial[1, 2:3] .= -20.
    initial[4, 2:3] .= 2.
    states(state) = reshape(state, 6, 3)
    laser(state, time) = [(1, helper) for helper in 2:3 if abs(state[1, helper]) < 10.]
    recorder = IntervalRecorder(states, laser, vec(initial), 10.)
    function straight_line!(derivative, state, params, time)
        fill!(derivative, 0.)
        for spacecraft in 1:3
            derivative[6*(spacecraft-1)+1] = state[6*(spacecraft-1)+4]
        end
    end
    solution = solve(ODEProblem(straight_line!, vec(initial), (0., 20.)), Tsit5();
        callback=CallbackSet(recording_callbacks(recorder)...), dtmax=2.)
    finish_recording!(recorder, solution.u[end], solution.t[end])
    @test sort(recorder.geometry.duration_s) ≈ [10., 10., 20.]
    contacts = filter(row -> row.sc_a == 1, recorder.geometry)
    @test contacts.start_time_s ≈ [5., 5.]
    @test contacts.end_time_s ≈ [15., 15.]
    @test !any(contacts.start_clipped) && !any(contacts.end_clipped)
    @test recorder.laser_on.start_time_s ≈ [5., 5.]
    @test recorder.laser_on.end_time_s ≈ [15., 15.]
    @test recorder.laser_on.duration_s ≈ [10., 10.]
    @test recorder.laser_on.encounter_id == contacts.encounter_id
end