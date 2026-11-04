using LinearAlgebra, Statistics, TOML, SHA, Test

const INTEGRATION_ROOT = normpath(joinpath(@__DIR__, "..", "9_SpaceAGORA.jl-main_full integration"))
include(joinpath(INTEGRATION_ROOT, "examples", "common.jl"))
using Arrow, CSV, DataFrames, StaticArrays, DiffEqBase
include(joinpath(@__DIR__, "EncounterTimingRecorder.jl"))
using .EncounterTimingRecorder

const COMPARISON_ROOT = normpath(joinpath(@__DIR__, "..", "..", "Summary_Table_for_Other_Labs"))
const OUTPUT_ROOT = joinpath(@__DIR__, "results", "N20_encounter_recording_validation")
const TARGET = 21
const RANGE_M = 200e3
const MASS_KG = 227.0
const FORCE_N = 100 * 10_000 / 299_792_458.0
scenario(inc) = "N20_h1000km_t1050km_ih0_it$(inc)_e0_nu0_T318567.042s_R200km_P10000W_B100_m227kg_gve_sma_J2true_dt10s"
filehash(path) = open(io -> bytes2hex(sha256(io)), path)
readvec(table, sc, field, row) =
    SVector{3, Float64}(ntuple(k -> table[Symbol("sc$(sc)_$(field)_$(k)")][row], 3))
table_state(table, row) = hcat([vcat(readvec(table, i, "pos", row),
    readvec(table, i, "vel", row)) for i in 1:TARGET]...)

function reference(inc)
    inc in (0, 5) || error("Only the requested 0 and 5 degree cases are supported.")
    bundle = joinpath(COMPARISON_ROOT, "2_SpaceAGORA.jl", "output", "feather", "spaceagora_" * scenario(inc))
    deliverable = joinpath(COMPARISON_ROOT, "6_Deliverables", scenario(inc))
    table = Arrow.Table(joinpath(bundle, "trajectory.feather"))
    geometry = DataFrame(Arrow.Table(joinpath(bundle, "geometry_encounters.feather")))
    laser = DataFrame(Arrow.Table(joinpath(bundle, "laser_on.feather")))
    published_geometry = CSV.read(joinpath(deliverable, "encounter_record.csv"), DataFrame)
    published_laser = CSV.read(joinpath(deliverable, "laser_link_record.csv"), DataFrame)
    @test select(geometry, names(published_geometry)) == published_geometry
    @test laser == published_laser
    @test Arrow.getmetadata(table)["target_id"] == string(TARGET)
    @test first(table.time) == 0.0
    return (; bundle, deliverable, table, geometry, laser)
end

function configuration(old, output; duration=last(old.time))
    planet = make_no_gram_planet(:earth)
    spacecraft = map(1:TARGET) do i
        bus = Link(root=true, m=MASS_KG)
        initial = CartesianInitialCondition(readvec(old, i, "pos", 1), readvec(old, i, "vel", 1))
        SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia,
            initial_condition=initial, id=i, n_terminal=1)
    end
    links = ConstellationInterLinkModel()
    for helper in 1:TARGET-1
        register_allowed_link!(links, spacecraft, (helper, 1), (TARGET, 1);
            parameters=InterLinkParameters(P=10_000.0, B=100.0, range=RANGE_M))
    end
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=true, generate_plots=false,
            save_csv=false, results_directory=output),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            mission_time=duration, orientation_sim=false, data_rate=10.0),
        environment_model=make_no_gram_environment(planet=planet),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredJ2GravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(year=2026, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-12,
            abstol_orbit=1e-12, dt_max_orbit=10.0),
        solver_config=SM.SolverConfig(solver_mode=:tsit5),
        interlink_model=links, scheduling_policy_model=SchedulingPolicyModel(:gve_sma; target_idx=TARGET))
end

sample_pairs(sample) = [minmax(edge[1][1], edge[2][1]) for edge in sample.active]

function passive_callback(initial, duration)
    recorder = Recorder(initial, RANGE_M)
    previous = copy(initial)
    previous_time = 0.0
    previous_dv = fill(zero(SVector{3, Float64}), TARGET)
    max_force_residual = Ref(0.0)
    max_inactive_impulse = Ref(0.0)
    samples = Ref(0)
    function affect!(integrator)
        time = Float64(integrator.t)
        current = snapshot(integrator.u)
        dt = time - previous_time
        history = integrator.p.args.interlink_model.history
        index = history[end].time > previous_time ? length(history) - 1 : length(history)
        sample = history[index]
        sample.time == previous_time || error("Missing scheduler state for the completed step")
        previous_active = sample_pairs(sample)
        observe_firing!(recorder, previous_active, previous, previous_time; initial=samples[] == 0)
        dense(t) = snapshot(integrator(t))
        observe_geometry!(recorder, previous, current, previous_time, time, dense)
        expected = zero(SVector{3, Float64})
        helper = 0
        if !isempty(previous_active)
            length(previous_active) == 1 || error("Multiple simultaneous target links")
            pair = only(previous_active)
            pair[2] == TARGET || error("Unexpected non-target link")
            helper = pair[1]
            middle = dense((previous_time + time) / 2)
            direction(state) = normalize(SVector{3, Float64}(state[1:3, TARGET] - state[1:3, helper]))
            expected = FORCE_N / MASS_KG * dt / 6 *
                (direction(previous) + 4direction(middle) + direction(current))
        end
        for i in 1:TARGET
            dv = SVector{3, Float64}(integrator.u.sc[i].laser_dv)
            increment = dv - previous_dv[i]
            prediction = i == TARGET ? expected : i == helper ? -expected : zero(expected)
            max_force_residual[] = max(max_force_residual[], norm(increment - prediction))
            if i != helper && (i != TARGET || helper == 0)
                max_inactive_impulse[] = max(max_inactive_impulse[], norm(increment))
            end
            previous_dv[i] = dv
        end
        previous, previous_time = current, time
        samples[] += 1
    end
    callback = DiffEqBase.DiscreteCallback((u, t, integrator) -> true, affect!;
        save_positions=(false, false))
    return (; callback, recorder, max_force_residual, max_inactive_impulse, samples)
end

function output_times(duration)
    times = collect(0.0:10.0:duration)
    last(times) < duration && push!(times, duration)
    return times
end

function run_case(inc; baseline=false, smoke=false)
    old = reference(inc)
    duration = smoke ? 3500.0 : last(old.table.time)
    label = smoke ? "smoke" : baseline ? "baseline" : "recorded"
    output = joinpath(OUTPUT_ROOT, "i$(inc)", label)
    isdir(output) && !isempty(readdir(output)) && error("Refusing to overwrite $output")
    mkpath(output)
    initial = table_state(old.table, 1)
    args = configuration(old.table, output; duration)
    observer = baseline ? nothing : passive_callback(initial, duration)
    fields = filter(f -> f.name in (:position, :velocity, :laser_dv, :laser_delta_sma),
        SM.default_save_fields(args))
    println("Running i=$inc, $label, duration=$duration s, all-pair range recording=$(observer !== nothing)")
    flush(stdout)
    solution = withenv("SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "false") do
        run_simulation(args; return_solution=true, save_fields=fields,
            extra_callbacks=observer === nothing ? () : (observer.callback,))
    end
    @test string(solution.retcode) == "Success"
    @test unique(solution.t) == [0.0, duration]
    @test snapshot(solution.u[1]) == initial
    @test all(isfinite, solution.u[end])
    history = solution.prob.p.args.interlink_model.history
    @test length(history) == solution.stats.naccept + 1
    @test all(length(s.active) <= 1 for s in history)
    schedule = DataFrame(time_s=[s.time for s in history],
        active_helper=[isempty(s.active) ? 0 : only(s.active)[1][1] for s in history])
    Arrow.write(joinpath(output, "accepted_step_schedule.feather"), schedule)
    table = Arrow.Table(joinpath(output, "simulation_results.feather"))
    @test table.time == output_times(duration)
    @test table_state(table, 1) == initial
    @test table_state(table, length(table.time)) == snapshot(solution.u[end])
    summary = Dict{String,Any}(
        "retcode" => string(solution.retcode), "duration_s" => duration, "inclination_deg" => inc,
        "accepted_steps" => solution.stats.naccept, "rejected_steps" => solution.stats.nreject,
        "max_step_s" => maximum(diff(schedule.time_s)), "target_id" => TARGET, "all_pair_count" => 210,
        "J2" => args.environment_model.planet.J2, "earth_radius_m" => args.environment_model.planet.Rp_e,
        "mu_m3_s2" => args.environment_model.planet.μ, "range_m" => RANGE_M,
        "force_N" => FORCE_N, "mass_kg" => MASS_KG, "initial_target_nu_deg" => 0.0,
        "reference_trajectory_sha256" => filehash(joinpath(old.bundle, "trajectory.feather")),
        "runner_sha256" => filehash(@__FILE__),
        "recorder_sha256" => filehash(joinpath(@__DIR__, "EncounterTimingRecorder.jl")),
        "schedule_applied_firing_s" => sum(diff(schedule.time_s) .* (schedule.active_helper[1:end-1] .> 0)),
        "recording_method" => baseline ? "No passive recorder" :
            "Passive accepted-step observer; range roots on current-step solver dense interpolation; no solver stops added")
    if observer !== nothing
        @test observer.samples[] == solution.stats.naccept
        @test observer.max_force_residual[] < 1e-9
        @test observer.max_inactive_impulse[] < 1e-12
        finish!(observer.recorder, snapshot(solution.u[end]), duration)
        @test observer.recorder.max_boundary_residual_m < 1e-3
        geometry, laser = observer.recorder.geometry, observer.recorder.laser
        encounters, audit = analyze_intervals(geometry, laser)
        for (name, data) in (("geometry_encounters", geometry), ("laser_on", laser),
            ("encounter_analysis", encounters), ("firing_range_audit", audit))
            Arrow.write(joinpath(output, "$name.feather"), data)
            CSV.write(joinpath(output, "$name.csv"), data)
        end
        @test isapprox(sum(laser.duration_s), summary["schedule_applied_firing_s"]; atol=1e-8)
        @test isapprox(sum(encounters.laser_on_s) + sum(audit.outside_range_s),
            sum(laser.duration_s); atol=1e-8)
        @test all(audit.outside_range_s .>= -1e-8)
        summary["encounter_count"] = nrow(geometry)
        summary["firing_interval_count"] = nrow(laser)
        summary["total_geometry_s"] = sum(geometry.duration_s)
        summary["total_actual_firing_s"] = sum(laser.duration_s)
        summary["total_in_range_firing_s"] = sum(encounters.laser_on_s)
        summary["total_out_of_range_firing_s"] = sum(audit.outside_range_s)
        summary["max_range_root_residual_m"] = observer.recorder.max_boundary_residual_m
        summary["max_step_impulse_simpson_residual_mps"] = observer.max_force_residual[]
        summary["max_inactive_step_impulse_mps"] = observer.max_inactive_impulse[]
    end
    open(joinpath(output, "summary.toml"), "w") do io
        TOML.print(io, summary; sorted=true)
    end
    println("Saved ", output)
    TOML.print(stdout, summary; sorted=true)
end

function synthetic_tests()
    fixture(y=0.0) = t -> [0.0 -2.0+4t; 0.0 y; 0.0 0.0; 0.0 4.0; 0.0 0.0; 0.0 0.0]
    @testset "Analytical range and firing fixtures" begin
        state = fixture()
        recorder = Recorder(state(0.0), 1.0)
        observe_geometry!(recorder, state(0.0), state(1.0), 0.0, 1.0, state)
        finish!(recorder, state(1.0), 1.0)
        @test nrow(recorder.geometry) == 1
        @test isapprox(recorder.geometry.start_time_s[1], 0.25; atol=1e-8)
        @test isapprox(recorder.geometry.end_time_s[1], 0.75; atol=1e-8)
        # Both endpoints outside, including a much shorter interior encounter.
        short = fixture(sqrt(1 - 1e-6))
        r = Recorder(short(0.0), 1.0)
        observe_geometry!(r, short(0.0), short(1.0), 0.0, 1.0, short)
        finish!(r, short(1.0), 1.0)
        @test nrow(r.geometry) == 1
        @test isapprox(r.geometry.duration_s[1], 0.0005; atol=1e-8)
        for y in (1.0, 2.0)
            outside = fixture(y)
            r = Recorder(outside(0.0), 1.0)
            observe_geometry!(r, outside(0.0), outside(1.0), 0.0, 1.0, outside)
            finish!(r, outside(1.0), 1.0)
            @test isempty(r.geometry)
            @test isempty(r.laser)
        end
        state = t -> [0.0 t; 0.0 0.0; 0.0 0.0; 0.0 1.0; 0.0 0.0; 0.0 0.0]
        r = Recorder(state(0.0), 0.5)
        observe_firing!(r, [(1, 2)], state(0.0), 0.0; initial=true)
        observe_geometry!(r, state(0.0), state(1.0), 0.0, 1.0, state)
        observe_firing!(r, Tuple{Int,Int}[], state(1.0), 1.0)
        finish!(r, state(1.0), 1.0)
        encounters, audit = analyze_intervals(r.geometry, r.laser)
        @test isapprox(only(encounters.duration_s), 0.5; atol=1e-8)
        @test only(r.laser.duration_s) == 1.0
        @test isapprox(only(audit.within_range_s), 0.5; atol=1e-8)
        @test isapprox(only(audit.outside_range_s), 0.5; atol=1e-8)
        @test only(r.geometry.start_clipped)
        initial = [0.0 0.2 0.4; zeros(5, 3)]
        r = Recorder(initial, 1.0)
        observe_firing!(r, [(1, 2)], initial, 0.0; initial=true)
        observe_firing!(r, [(1, 3)], initial, 0.4)
        finish!(r, initial, 1.0)
        @test r.laser.duration_s == [0.4, 0.6]
        @test all(r.geometry.end_clipped)
        @test r.laser.end_clipped == [false, true]
    end
end

function compare_case(inc)
    old = reference(inc)
    root = joinpath(OUTPUT_ROOT, "i$(inc)")
    recorded = joinpath(root, "recorded")
    baseline = joinpath(root, "baseline")
    table = Arrow.Table(joinpath(recorded, "simulation_results.feather"))
    baseline_table = Arrow.Table(joinpath(baseline, "simulation_results.feather"))
    schedule = DataFrame(Arrow.Table(joinpath(recorded, "accepted_step_schedule.feather")))
    baseline_schedule = DataFrame(Arrow.Table(joinpath(baseline, "accepted_step_schedule.feather")))
    @test schedule == baseline_schedule
    @test propertynames(table) == propertynames(baseline_table)
    @test all(table[name] == baseline_table[name] for name in propertynames(table))
    geometry = DataFrame(Arrow.Table(joinpath(recorded, "geometry_encounters.feather")))
    laser = DataFrame(Arrow.Table(joinpath(recorded, "laser_on.feather")))
    encounters, audit = analyze_intervals(geometry, laser)
    @test encounters == DataFrame(Arrow.Table(joinpath(recorded, "encounter_analysis.feather")))
    @test audit == DataFrame(Arrow.Table(joinpath(recorded, "firing_range_audit.feather")))
    @test all(geometry.duration_s .== geometry.end_time_s .- geometry.start_time_s)
    @test all(laser.duration_s .== laser.end_time_s .- laser.start_time_s)
    @test all(geometry.duration_s .> 0)
    @test all(laser.duration_s .> 0)
    @test all(laser.sc_b .== TARGET)
    @test all(laser.start_time_s .∈ Ref(Set(schedule.time_s)))
    @test all(laser.end_time_s .∈ Ref(Set(schedule.time_s)))
    @test isapprox(sum(laser.duration_s),
        sum(diff(schedule.time_s) .* (schedule.active_helper[1:end-1] .> 0)); atol=1e-8)
    # An independent reconstruction uses the regular 10-second output, not accepted steps.
    reconstructed = Recorder(table_state(table, 1), RANGE_M)
    previous = table_state(table, 1)
    for row in 2:length(table.time)
        current = table_state(table, row)
        start, stop = table.time[row - 1], table.time[row]
        dense(t) = hermite_state(previous, current, stop - start, (t - start) / (stop - start))
        observe_geometry!(reconstructed, previous, current, start, stop, dense)
        previous = current
    end
    finish!(reconstructed, previous, last(table.time))
    rebuilt = reconstructed.geometry
    @test select(geometry, :sc_a, :sc_b) == select(rebuilt, :sc_a, :sc_b)
    @test nrow(geometry) == nrow(rebuilt)
    reconstruction_error = maximum(vcat(abs.(geometry.start_time_s .- rebuilt.start_time_s),
        abs.(geometry.end_time_s .- rebuilt.end_time_s)))
    @test reconstruction_error < 1e-3
    CSV.write(joinpath(root, "independent_10s_reconstruction.csv"), rebuilt)
    old_encounters, _ = analyze_intervals(old.geometry, old.laser)
    sort!(old_encounters, [:sc_a, :sc_b, :start_time_s])
    sort!(encounters, [:sc_a, :sc_b, :start_time_s])
    comparison = DataFrame(sc_a=Int[], sc_b=Int[], occurrence=Int[],
        old_start_s=Union{Missing,Float64}[], new_start_s=Union{Missing,Float64}[],
        entry_difference_s=Union{Missing,Float64}[], exit_difference_s=Union{Missing,Float64}[],
        old_geometry_s=Union{Missing,Float64}[], new_geometry_s=Union{Missing,Float64}[],
        geometry_difference_s=Union{Missing,Float64}[], old_laser_on_s=Union{Missing,Float64}[],
        new_in_range_laser_s=Union{Missing,Float64}[], in_range_laser_difference_s=Union{Missing,Float64}[])
    pairs = sort!(union([(r.sc_a, r.sc_b) for r in eachrow(old_encounters)],
        [(r.sc_a, r.sc_b) for r in eachrow(encounters)]))
    unmatched = 0
    for pair in pairs
        left = filter(r -> (r.sc_a, r.sc_b) == pair, old_encounters)
        right = filter(r -> (r.sc_a, r.sc_b) == pair, encounters)
        for i in 1:max(nrow(left), nrow(right))
            a = i <= nrow(left) ? left[i, :] : nothing
            b = i <= nrow(right) ? right[i, :] : nothing
            value(row, key) = row === nothing ? missing : row[key]
            a === nothing || b === nothing ? (unmatched += 1) : nothing
            push!(comparison, (pair..., i, value(a, :start_time_s), value(b, :start_time_s),
                value(b, :start_time_s) - value(a, :start_time_s),
                value(b, :end_time_s) - value(a, :end_time_s),
                value(a, :duration_s), value(b, :duration_s),
                value(b, :duration_s) - value(a, :duration_s),
                value(a, :laser_on_s), value(b, :laser_on_s),
                value(b, :laser_on_s) - value(a, :laser_on_s)))
        end
    end
    CSV.write(joinpath(root, "encounter_comparison.csv"), comparison)
    firing_tables = DataFrame[]
    for (prefix, source) in (("old", old.laser), ("new", laser))
        frame = select(source, :sc_a, :sc_b, :start_time_s, :end_time_s, :duration_s)
        sort!(frame, [:sc_a, :sc_b, :start_time_s])
        transform!(groupby(frame, [:sc_a, :sc_b]),
            :start_time_s => (times -> collect(eachindex(times))) => :occurrence)
        rename!(frame, :start_time_s => Symbol(prefix * "_start_s"),
            :end_time_s => Symbol(prefix * "_end_s"),
            :duration_s => Symbol(prefix * "_duration_s"))
        push!(firing_tables, frame)
    end
    firing_comparison = outerjoin(firing_tables...;
        on=[:sc_a, :sc_b, :occurrence], validate=(true, true))
    sort!(firing_comparison, [:sc_a, :sc_b, :occurrence])
    firing_comparison.start_difference_s = firing_comparison.new_start_s .- firing_comparison.old_start_s
    firing_comparison.end_difference_s = firing_comparison.new_end_s .- firing_comparison.old_end_s
    firing_comparison.duration_difference_s = firing_comparison.new_duration_s .- firing_comparison.old_duration_s
    CSV.write(joinpath(root, "laser_interval_comparison.csv"), firing_comparison)
    @test all(audit.outside_range_s .>= -1e-8)
    @test maximum(audit.outside_range_s; init=0.0) <= 10.0 + 1e-8
    @test isapprox(sum(audit.within_range_s) + sum(audit.outside_range_s), sum(laser.duration_s); atol=1e-8)
    summary = Dict{String,Any}(
        "old_encounter_count" => nrow(old_encounters), "new_encounter_count" => nrow(encounters),
        "unmatched_pair_occurrences" => unmatched, "old_firing_interval_count" => nrow(old.laser),
        "new_firing_interval_count" => nrow(laser),
        "old_total_geometry_s" => sum(old_encounters.duration_s),
        "new_total_geometry_s" => sum(encounters.duration_s),
        "old_recorded_firing_s" => sum(old.laser.duration_s),
        "new_actual_firing_s" => sum(laser.duration_s),
        "new_in_range_firing_s" => sum(encounters.laser_on_s),
        "new_out_of_range_firing_s" => sum(audit.outside_range_s),
        "max_out_of_range_firing_per_interval_s" => maximum(audit.outside_range_s; init=0.0),
        "old_startup_firing_s" => sum(old.laser.duration_s[old.laser.start_time_s .== 0.0]),
        "new_startup_firing_s" => sum(laser.duration_s[laser.start_time_s .== 0.0]),
        "max_abs_entry_difference_s" => maximum(abs, skipmissing(comparison.entry_difference_s)),
        "max_abs_exit_difference_s" => maximum(abs, skipmissing(comparison.exit_difference_s)),
        "max_abs_geometry_duration_difference_s" => maximum(abs, skipmissing(comparison.geometry_difference_s)),
        "max_abs_in_range_firing_difference_s" => maximum(abs, skipmissing(comparison.in_range_laser_difference_s)),
        "independent_reconstruction_max_boundary_difference_s" => reconstruction_error,
        "baseline_schedule_identical" => schedule == baseline_schedule,
        "baseline_trajectory_identical" => all(table[name] == baseline_table[name] for name in propertynames(table)))
    open(joinpath(root, "comparison.toml"), "w") do io
        TOML.print(io, summary; sorted=true)
    end
    println("Comparison i=$inc")
    TOML.print(stdout, summary; sorted=true)
end

if abspath(PROGRAM_FILE) == @__FILE__
@testset "Encounter recording validation" begin
    if ARGS == ["fixtures"]
        synthetic_tests()
    elseif length(ARGS) == 2 && ARGS[1] in ("run", "baseline", "smoke")
        run_case(parse(Int, ARGS[2]); baseline=ARGS[1] == "baseline", smoke=ARGS[1] == "smoke")
    elseif length(ARGS) == 2 && ARGS[1] == "compare"
        compare_case(parse(Int, ARGS[2]))
    else
        error("Usage: validate_encounter_recording.jl fixtures | run INC | baseline INC | smoke INC | compare INC")
    end
    end
end
