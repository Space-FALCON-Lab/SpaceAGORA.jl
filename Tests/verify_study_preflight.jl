include(joinpath(@__DIR__, "preflight_study_checks.jl"))

function load_checked(label)
    directory = joinpath(PREFLIGHT_ROOT, label)
    summary = TOML.parsefile(joinpath(directory, "summary.toml"))
    table = Arrow.Table(joinpath(directory, "simulation_results.feather"))
    diagnostics = DataFrame(Arrow.Table(joinpath(directory, "diagnostics.feather")))
    schedule = DataFrame(Arrow.Table(joinpath(directory, "schedule.feather")))
    geometry = DataFrame(Arrow.Table(joinpath(directory, "geometry.feather")))
    laser = DataFrame(Arrow.Table(joinpath(directory, "laser.feather")))
    segments = DataFrame(Arrow.Table(joinpath(directory, "range_segments.feather")))
    expected = collect(0.0:summary["save_interval_s"]:summary["duration_s"])
    last(expected) < summary["duration_s"] && push!(expected, summary["duration_s"])
    @test table.time == expected
    @test diagnostics.time_s == expected
    @test summary["retcode"] == "Success"
    @test all(all(isfinite, table[name]) for name in propertynames(table))
    @test all(all(isfinite, diagnostics[!, name]) for name in names(diagnostics))
    @test Vector(diagnostics[end, 2:4]) == summary["dv_RTN_mps"]
    @test isapprox(sum(Vector(diagnostics[end, 7:ncol(diagnostics)])),
        summary["optical_work_total_J"]; atol=1e-8)
    @test first(schedule.time_s) == 0.0 && last(schedule.time_s) == last(expected)
    @test all(diff(schedule.time_s) .> 0)
    @test length(schedule.time_s) == summary["accepted_steps"] + 1
    @test maximum(diff(schedule.time_s)) <= summary["dtmax_s"] + 2eps(last(expected))
    @test all(0 <= helper <= 20 for helper in schedule.active_helper)
    @test all(laser.sc_b .== 21)
    @test all(1 .<= geometry.sc_a .< geometry.sc_b .<= 21)
    @test all(laser.start_time_s .∈ Ref(Set(schedule.time_s)))
    @test all(laser.end_time_s .∈ Ref(Set(schedule.time_s)))
    @test length(unique(laser.laser_id)) == nrow(laser)
    @test length(unique(geometry.encounter_id)) == nrow(geometry)
    @test nrow(geometry) == summary["encounter_count"]
    @test nrow(laser) == summary["firing_count"]
    @test isapprox(sum(laser.duration_s), summary["firing_s"]; atol=1e-8)
    @test isapprox(sum(diff(schedule.time_s) .* (schedule.active_helper[1:end-1] .> 0)),
        summary["firing_s"]; atol=1e-8)
    @test summary["max_energy_work_error_J"] < summary["energy_work_tolerance_J"]
    @test summary["max_relative_angular_momentum_error"] < 1e-9
    @test summary["min_radius_m"] - RANGE_M > summary["earth_radius_m"]
    check_segments(geometry, laser, segments)
    # Segment flags must agree with the recorded trajectory at every regular sample.
    for segment in eachrow(segments)
        left = searchsortedfirst(table.time, segment.start_time_s + 1e-6)
        right = searchsortedlast(table.time, segment.end_time_s - 1e-6)
        for row in left:right
            range = norm(readvec(table, segment.sc_b, "pos", row) -
                readvec(table, segment.sc_a, "pos", row))
            @test segment.within_range ? range <= RANGE_M + 1e-3 : range >= RANGE_M - 1e-3
        end
    end
    return (; summary, table, diagnostics, schedule, geometry, laser, segments)
end

function encounter_extrema(run)
    isempty(run.geometry) && return (speed_min=missing, speed_max=missing,
        range_rate_min=missing, range_rate_max=missing)
    speeds, rates = Float64[], Float64[]
    function add(r, v)
        push!(speeds, norm(v))
        norm(r) > 0 || error("Zero separation in encounter kinematics")
        push!(rates, dot(r, v) / norm(r))
    end
    for event in eachrow(run.geometry)
        firstrow = searchsortedfirst(run.table.time, event.start_time_s)
        lastrow = searchsortedlast(run.table.time, event.end_time_s)
        for row in firstrow:lastrow
            add(readvec(run.table, event.sc_b, "pos", row) - readvec(run.table, event.sc_a, "pos", row),
                readvec(run.table, event.sc_b, "vel", row) - readvec(run.table, event.sc_a, "vel", row))
        end
        for boundary in ("start", "end")
            vector(sc, field) = SVector{3,Float64}(ntuple(k ->
                event[Symbol("$(boundary)_sc_$(sc)_$(field)_$(("x","y","z")[k])")], 3))
            add(vector("b", "r") - vector("a", "r"), vector("b", "v") - vector("a", "v"))
        end
    end
    return (speed_min=minimum(speeds), speed_max=maximum(speeds),
        range_rate_min=minimum(rates), range_rate_max=maximum(rates))
end

function paired_changes(a, b)
    left = sort(a.geometry, [:sc_a, :sc_b, :start_time_s])
    right = sort(b.geometry, [:sc_a, :sc_b, :start_time_s])
    same_pairs = nrow(left) == nrow(right) &&
        left.sc_a == right.sc_a && left.sc_b == right.sc_b
    entry_error = same_pairs && !isempty(left) ?
        maximum(abs, right.start_time_s - left.start_time_s) : missing
    exit_error = same_pairs && !isempty(left) ?
        maximum(abs, right.end_time_s - left.end_time_s) : missing
    duration_error = same_pairs && !isempty(left) ?
        maximum(abs, right.duration_s - left.duration_s) : missing
    last_a, last_b = length(a.table.time), length(b.table.time)
    return (same_pair_occurrences=same_pairs,
        encounter_count_change=nrow(right)-nrow(left),
        max_entry_difference_s=entry_error, max_exit_difference_s=exit_error,
        max_geometry_duration_difference_s=duration_error,
        delta_sma_difference_m=b.summary["delta_sma_m"]-a.summary["delta_sma_m"],
        delta_inc_difference_deg=b.summary["delta_inc_deg"]-a.summary["delta_inc_deg"],
        dv_R_difference_mps=b.summary["dv_RTN_mps"][1]-a.summary["dv_RTN_mps"][1],
        dv_T_difference_mps=b.summary["dv_RTN_mps"][2]-a.summary["dv_RTN_mps"][2],
        dv_N_difference_mps=b.summary["dv_RTN_mps"][3]-a.summary["dv_RTN_mps"][3],
        firing_difference_s=b.summary["firing_s"]-a.summary["firing_s"],
        final_position_difference_m=norm(readvec(b.table,21,"pos",last_b)-readvec(a.table,21,"pos",last_a)),
        final_velocity_difference_mps=norm(readvec(b.table,21,"vel",last_b)-readvec(a.table,21,"vel",last_a)))
end

function verify_all()
    rows = NamedTuple[]
    differences = NamedTuple[]
    checks = Dict{String,Any}()
    for case in ("sma850", "inc1050")
        base_label = "$(case)_dt10.0_save10.0"
        base = load_checked(base_label)
        reconstructed = Recorder(table_state(base.table, 1), RANGE_M)
        previous_state = table_state(base.table, 1)
        for row in 2:length(base.table.time)
            current = table_state(base.table, row)
            start, stop = base.table.time[row - 1], base.table.time[row]
            dense(t) = hermite_state(previous_state, current, stop-start, (t-start)/(stop-start))
            observe_geometry!(reconstructed, previous_state, current, start, stop, dense)
            previous_state = current
        end
        finish!(reconstructed, previous_state, last(base.table.time))
        @test select(reconstructed.geometry, :sc_a, :sc_b) ==
            select(base.geometry, :sc_a, :sc_b)
        @test nrow(reconstructed.geometry) == nrow(base.geometry)
        boundary_difference = maximum(vcat(abs.(reconstructed.geometry.start_time_s - base.geometry.start_time_s),
            abs.(reconstructed.geometry.end_time_s - base.geometry.end_time_s)))
        @test boundary_difference < 1e-3
        checks["$(case)_independent_geometry_max_boundary_difference_s"] = boundary_difference
        function record_row(run)
            s = run.summary
            push!(rows, (; label=s["label"], dtmax_s=s["dtmax_s"], save_interval_s=s["save_interval_s"],
                delta_sma_m=s["delta_sma_m"], delta_inc_deg=s["delta_inc_deg"],
                dv_R_mps=s["dv_RTN_mps"][1], dv_T_mps=s["dv_RTN_mps"][2], dv_N_mps=s["dv_RTN_mps"][3],
                encounter_count=s["encounter_count"], firing_s=s["firing_s"],
                outside_firing_s=s["outside_firing_s"],
                max_energy_work_error_J=s["max_energy_work_error_J"],
                max_relative_momentum_error=s["max_relative_angular_momentum_error"],
                encounter_extrema(run)...))
        end
        record_row(base)
        previous = base
        for step in (5.0, 2.5)
            run = load_checked("$(case)_dt$(step)_save10.0")
            record_row(run)
            push!(differences, (; left=previous.summary["label"], right=run.summary["label"],
                paired_changes(previous, run)...))
            previous = run
        end
        for interval in (5.0, 1.0)
            run = load_checked("$(case)_dt10.0_save$(interval)")
            record_row(run)
            @test run.schedule == base.schedule
            @test run.geometry == base.geometry
            @test run.laser == base.laser
            common = indexin(base.table.time, run.table.time)
            @test !any(isnothing, common)
            indices = Int[i for i in common]
            @test all(base.table[name] == run.table[name][indices] for name in propertynames(base.table))
            @test all(base.diagnostics[!, name] == run.diagnostics[indices, name]
                for name in names(base.diagnostics))
            checks["$(case)_save$(interval)_schedule_and_common_states_identical"] = true
            push!(differences, (; left=base_label, right=run.summary["label"],
                paired_changes(base, run)...))
        end
        control = joinpath(PREFLIGHT_ROOT, base_label * "_control")
        control_table = Arrow.Table(joinpath(control, "simulation_results.feather"))
        control_schedule = DataFrame(Arrow.Table(joinpath(control, "schedule.feather")))
        @test base.schedule == control_schedule
        @test all(base.table[name] == control_table[name] for name in propertynames(base.table))
        @test base.diagnostics == DataFrame(Arrow.Table(joinpath(control, "diagnostics.feather")))
        checks["$(case)_event_recorder_invariance"] = true
    end
    no_j2 = load_checked("inc1050_dt10.0_save10.0_noJ2")
    checks["noJ2_full_vector_relative_momentum_error"] = no_j2.summary["max_relative_angular_momentum_error"]
    checks["noJ2_max_energy_work_error_J"] = no_j2.summary["max_energy_work_error_J"]
    empty = DataFrame(Arrow.Table(joinpath(PREFLIGHT_ROOT, "no_firing_segments.feather")))
    @test isempty(empty)
    @test eltype(empty.duration_s) == Float64
    multiple = DataFrame(Arrow.Table(joinpath(PREFLIGHT_ROOT, "multiple_crossing_segments.feather")))
    @test nrow(multiple) == 7
    @test isapprox(sum(multiple.duration_s), 1.0; atol=1e-9)
    @test multiple.end_time_s[1:end-1] == multiple.start_time_s[2:end]
    CSV.write(joinpath(PREFLIGHT_ROOT, "sensitivity_runs.csv"), DataFrame(rows))
    CSV.write(joinpath(PREFLIGHT_ROOT, "sensitivity_differences.csv"), DataFrame(differences))
    open(joinpath(PREFLIGHT_ROOT, "verification.toml"), "w") do io
        TOML.print(io, checks; sorted=true)
    end
    show(stdout, DataFrame(rows); allrows=true, allcols=true)
    println()
    show(stdout, DataFrame(differences); allrows=true, allcols=true)
    println()
    TOML.print(stdout, checks; sorted=true)
end

@testset "Persisted preflight results" begin
    verify_all()
end
