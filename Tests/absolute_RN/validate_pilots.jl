include(joinpath(@__DIR__, "..", "..", "Simulation_absolute_RN", "Code", "VerifyOutput.jl"))
using .StudyVerification
using .StudyVerification.StudyConfig
using Arrow, DataFrames, LinearAlgebra, StaticArrays, Test, TOML

module BaselineVerification
include(joinpath(@__DIR__, "results", "source_snapshot_pre_absolute_RN", "project",
    "Simulation", "Code", "VerifyOutput.jl"))
end

const RESULTS = joinpath(@__DIR__, "results")
report_path(label) = joinpath(RESULTS, label * ".toml")
load_report(label) = TOML.parsefile(report_path(label))
trajectory(report) = Arrow.Table(joinpath(report["output_directory"], "trajectory.feather"))

function verify(label)
    report = load_report(label)
    report["fresh_file_verified"] = false
    atomic_toml(report_path(label), report)
    directory = report["output_directory"]
    if report["kind"] == "new"
        verify_output(directory; require_events=report["mode"] != "control")
    else
        BaselineVerification.StudyVerification.verify_output(directory)
    end
    actual = trajectory(report)
    @testset "Fresh pilot files: $label" begin
        @test length(actual.time) == 1001
        @test first(actual.time) == 0 && last(actual.time) == report["settings"]["duration_s"]
        schedule = Arrow.Table(joinpath(directory, "schedule.feather"))
        @test first(schedule.time) == 0 && last(schedule.time) == last(actual.time)
        @test length(schedule.time) == report["solver_accepted_steps"] + 1
        @test all(0 .< diff(schedule.time) .<= STEP_S + 2eps(last(actual.time)))
        @test all(0 .<= schedule.active_helper .<= report["settings"]["helper_count"])
        if report["mode"] == "diagnostic"
            @test report["passed"]
            reference = Arrow.Table(report["reference_file"])
            @test actual.time == reference.time
            maxima = Dict{String,Any}()
            for (axis, component) in enumerate(("R", "T", "N"))
                net = Symbol("target_dv_$(component)_mps")
                positive = Symbol("target_dv_$(component)_positive_mps")
                negative = Symbol("target_dv_$(component)_negative_mps")
                absolute = Symbol("target_dv_$(component)_absolute_mps")
                net_error = maximum(abs.(actual[net] - reference[net]))
                positive_error = maximum(abs.(actual[positive] - reference[positive]))
                negative_error = maximum(abs.(actual[negative] - reference[negative]))
                absolute_error = maximum(abs.(actual[positive] - actual[negative] - reference[absolute]))
                @test net_error <= report["net_tolerance_mps"]
                @test positive_error <= report["branch_tolerances_mps"][axis]
                @test negative_error <= report["branch_tolerances_mps"][axis]
                @test absolute_error <= report["branch_tolerances_mps"][axis]
                maxima[component] = Dict("net_mps" => net_error, "positive_mps" => positive_error,
                    "negative_mps" => negative_error, "absolute_mps" => absolute_error)
            end
            for name in propertynames(reference)
                endswith(string(name), "_laser_work_J") || continue
                @test maximum(abs.(actual[name] - reference[name])) <= report["work_tolerance_J"]
            end
            report["fresh_reference_errors"] = maxima
            report["fresh_reference_verified"] = true
        end
    end
    report["fresh_file_verified"] = true
    atomic_toml(report_path(label), report)
end

function schedule_disagreement(a, b)
    times = sort!(unique!(vcat(collect(a.time), collect(b.time))))
    ia, ib, disagreement = 1, 1, 0.0
    for index in 1:length(times)-1
        start, stop = times[index], times[index + 1]
        while ia < length(a.time) && a.time[ia + 1] <= start
            ia += 1
        end
        while ib < length(b.time) && b.time[ib + 1] <= start
            ib += 1
        end
        if a.active_helper[ia] != b.active_helper[ib]
            disagreement += stop - start
        end
    end
    disagreement
end

function compare(left_label, right_label, mode)
    mode in ("exact", "adaptive") || error("Comparison mode must be exact or adaptive")
    destination = joinpath(RESULTS, "comparison_$(left_label)_$(right_label).toml")
    atomic_toml(destination, Dict("passed" => false, "status" => "checking"))
    left, right = load_report(left_label), load_report(right_label)
    a, b = trajectory(left), trajectory(right)
    la = Arrow.Table(joinpath(left["output_directory"], "schedule.feather"))
    lb = Arrow.Table(joinpath(right["output_directory"], "schedule.feather"))
    duration = left["settings"]["duration_s"]
    disagreement = schedule_disagreement(la, lb)
    maximum_position, maximum_velocity = 0.0, 0.0
    for sc in 1:left["settings"]["helper_count"]+1
        position2, velocity2 = zeros(length(a.time)), zeros(length(a.time))
        for axis in 1:3
            position = Symbol("sc$(sc)_pos_$axis")
            velocity = Symbol("sc$(sc)_vel_$axis")
            position2 .+= (a[position] - b[position]) .^ 2
            velocity2 .+= (a[velocity] - b[velocity]) .^ 2
        end
        maximum_position = max(maximum_position, sqrt(maximum(position2)))
        maximum_velocity = max(maximum_velocity, sqrt(maximum(velocity2)))
    end
    acceleration = left["settings"]["force_N"] / MASS_KG
    # Held-link differences are at most 2a. Allow near-circular position response,
    # rotating-frame velocity response, and n*delta-r when comparing ECI velocities.
    position_response = 10duration * acceleration * disagreement
    mean_motion = sqrt(MU / (EARTH_RADIUS_M + 1000 * left["settings"]["target_altitude_km"])^3)
    position_tolerance = 0.02 + position_response
    velocity_tolerance = 5e-5 + 14acceleration * disagreement + 1.01mean_motion * position_response
    @testset "$mode physical/scheduler comparison: $left_label / $right_label" begin
        @test left["fresh_file_verified"] && right["fresh_file_verified"]
        @test a.time == b.time
        for key in intersect(keys(left["settings"]), keys(right["settings"]))
            key in ("schema_version", "diagnostics") && continue
            @test left["settings"][key] == right["settings"][key]
        end
        left_info = TOML.parsefile(joinpath(left["output_directory"], "run_info.toml"))
        right_info = TOML.parsefile(joinpath(right["output_directory"], "run_info.toml"))
        @test left_info["initial_states"] == right_info["initial_states"]
        if mode == "exact"
            @test collect(propertynames(a)) == collect(propertynames(b))
            @test all(a[name] == b[name] for name in propertynames(a))
            @test left["solver_accepted_steps"] == right["solver_accepted_steps"]
            @test left["solver_rejected_steps"] == right["solver_rejected_steps"]
            @test left["schedule_sha256"] == right["schedule_sha256"]
            @test disagreement == 0
        else
            @test left["kind"] != right["kind"]
            @test disagreement / duration <= 0.02
            @test maximum_position <= position_tolerance
            @test maximum_velocity <= velocity_tolerance
        end
    end
    report = Dict("passed" => true, "left" => left_label, "right" => right_label,
        "mode" => mode, "maximum_position_difference_m" => maximum_position,
        "maximum_velocity_difference_mps" => maximum_velocity,
        "schedule_disagreement_s" => disagreement,
        "schedule_disagreement_fraction" => disagreement / duration,
        "position_tolerance_m" => mode == "exact" ? 0.0 : position_tolerance,
        "velocity_tolerance_mps" => mode == "exact" ? 0.0 : velocity_tolerance,
        "schedule_disagreement_fraction_limit" => 0.02)
    atomic_toml(destination, report)
    TOML.print(stdout, report; sorted=true)
end

function coverage(labels)
    atomic_toml(joinpath(RESULTS, "sign_coverage.toml"), Dict("passed" => false, "status" => "checking"))
    reports = load_report.(labels)
    values = Dict{String,Any}()
    @testset "Meaningful positive and negative R/N coverage for both schedulers" begin
        @test all(report["mode"] == "diagnostic" && report["passed"] &&
            report["fresh_file_verified"] && report["fresh_reference_verified"] for report in reports)
        for policy in ("gve_sma", "gve_inc"), (component, axis) in (("R", 1), ("N", 3))
            selected = filter(report -> report["settings"]["scheduler"] == policy &&
                report["settings"]["target_inclination_deg"] > 0, reports)
            @test !isempty(selected)
            positive = maximum(report["final_independent_positive_mps"][axis] for report in selected)
            negative = minimum(report["final_independent_negative_mps"][axis] for report in selected)
            @test positive > 1e-7
            @test negative < -1e-7
            values["$(policy)_$component"] = Dict("positive_mps" => positive, "negative_mps" => negative)
        end
    end
    atomic_toml(joinpath(RESULTS, "sign_coverage.toml"),
        Dict("passed" => true, "pilots" => labels, "coverage" => values))
end

isempty(ARGS) && error("Use verify LABEL..., compare LABEL LABEL exact|adaptive, or coverage LABEL...")
if ARGS[1] == "verify"
    foreach(verify, ARGS[2:end])
elseif ARGS[1] == "compare"
    length(ARGS) == 4 || error("compare LABEL LABEL exact|adaptive")
    compare(ARGS[2], ARGS[3], ARGS[4])
elseif ARGS[1] == "coverage"
    coverage(ARGS[2:end])
else
    error("Unknown pilot validation command")
end
