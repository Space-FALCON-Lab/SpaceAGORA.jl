if !isdefined(@__MODULE__, :REPO_ROOT)
    include(joinpath(@__DIR__, "common.jl"))
end
using Random
using LinearAlgebra
using Statistics
using Printf
using CSV
using DataFrames
using DiffEqBase

edge_key(graph::InterLinkGraph, edge::Int) =
    (graph.vertices[graph.edges[edge][1]], graph.vertices[graph.edges[edge][2]])

function walker_constellation(; satellites::Int=4000, planes::Int=80,
    seed::Int=20261003, altitude::Float64=550e3, inclination::Float64=53.0,
    target_idx::Int=1, range_m::Float64=2e6, all_pairs::Bool=false,
    policy_name::Symbol=:gve_sma)
    satellites > 1 || throw(ArgumentError("satellites must exceed one."))
    1 <= planes <= satellites && satellites % planes == 0 ||
        throw(ArgumentError("planes must be a positive divisor of satellites."))
    1 <= target_idx <= satellites || throw(ArgumentError("target_idx is outside the constellation."))
    isfinite(altitude) && altitude > 0 || throw(ArgumentError("altitude must be finite and positive."))
    isfinite(inclination) && 0 <= inclination <= 180 ||
        throw(ArgumentError("inclination must be in [0, 180] degrees."))
    parameters = InterLinkParameters(range=range_m)
    policy = SchedulingPolicyModel(policy_name; target_idx)
    generator = MersenneTwister(seed)
    phasing = rand(generator, 0:(planes - 1))
    raan_offset, anomaly_offset = 360rand(generator), 360rand(generator)
    planet = make_no_gram_planet(:earth)
    per_plane = satellites ÷ planes
    spacecraft = map(1:satellites) do index
        plane, slot = divrem(index - 1, per_plane)
        bus = Link(root=true, m=227.0)
        initial_condition = InitialCondition(a=planet.Rp_e + altitude, i=inclination,
            Ω=mod(raan_offset + 360plane / planes, 360),
            ν=mod(anomaly_offset + 360slot / per_plane + 360phasing * plane / satellites, 360))
        SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia,
            initial_condition=initial_condition, id=index, n_terminal=1)
    end
    model = ConstellationInterLinkModel()
    if all_pairs
        for first_satellite in 1:(satellites - 1), second_satellite in (first_satellite + 1):satellites
            register_allowed_link!(model, spacecraft, (first_satellite, 1), (second_satellite, 1); parameters)
        end
    else
        for partner in eachindex(spacecraft)
            partner == target_idx && continue
            register_allowed_link!(model, spacecraft, (target_idx, 1), (partner, 1); parameters)
        end
    end
    metadata = (; satellites, planes, per_plane, phasing, seed, altitude, inclination,
        raan_offset, anomaly_offset, range_m, target_idx, all_pairs, policy_name)
    return (; spacecraft, planet, model, policy, metadata)
end

function walker_snapshot(constellation)
    return (sc=map(constellation.spacecraft) do vehicle
        pos, vel = SimulationEngine.orbitalelemtorv(vehicle.initial_condition, constellation.planet)
        (; pos, vel, mass=227.0)
    end,)
end

"""Evaluate circular two-body motion from the initial state; laser forces do not alter these snapshots."""
function advance_walker_snapshot!(state, initial, mean_motion, time)
    sine, cosine = sincos(mean_motion * time)
    for index in eachindex(state.sc)
        current, origin = state.sc[index], initial.sc[index]
        @. current.pos = cosine * origin.pos + (sine / mean_motion) * origin.vel
        @. current.vel = -mean_motion * sine * origin.pos + cosine * origin.vel
    end
    return state
end

function validate_walker_selection(model, policy, selected)
    occupied = Set{vertex}()
    selected_set = Set(selected)
    for key in selected
        connection = model.graph.edge_index_lookup[key]
        @assert model.graph.available[connection] && model.graph.active[connection] && policy.matching.score[connection] > 0.0
        for endpoint in key
            @assert !(endpoint in occupied)
            push!(occupied, endpoint)
        end
    end
    @assert all(model.graph.active[connection] == (key in selected_set) for (connection, key) in enumerate([edge_key(model.graph, e) for e in eachindex(model.graph.edges)]))
    return nothing
end

function schedule_table(history)
    table = DataFrame(time_s=Float64[], first_satellite=Int[], first_terminal=Int[],
        second_satellite=Int[], second_terminal=Int[])
    for sample in history, key in sample.active
        push!(table, (sample.time, key[1]..., key[2]...))
    end
    return table
end

function benchmark_walker(constellation; duration::Float64=3600.0, interval::Float64=60.0)
    isfinite(duration) && duration > 0 || throw(ArgumentError("duration must be finite and positive."))
    isfinite(interval) && 0 < interval <= duration ||
        throw(ArgumentError("interval must be positive and no larger than duration."))
    (; spacecraft, planet, model, policy) = constellation
    initial = walker_snapshot(constellation)
    state = deepcopy(initial)
    mean_motion = sqrt(planet.μ / spacecraft[1].initial_condition.a^3)
    schedule_interlinks!(model, policy, spacecraft, state, planet.μ; candidate_force=force_on_endpoint)
    metrics = DataFrame(time_s=Float64[], available=Int[], selected=Int[], entered=Int[], left=Int[],
        availability_s=Float64[], scoring_s=Float64[], matching_s=Float64[], total_s=Float64[],
        availability_bytes=Int[], scoring_bytes=Int[], matching_bytes=Int[], allocated_bytes=Int[])
    history = InterLinkHistory[]
    for time in 0.0:interval:duration
        advance_walker_snapshot!(state, initial, mean_motion, time)
        previous_available = [model.graph.available[connection] for connection in eachindex(model.graph.edges)]
        availability = @timed update_availability!(model, spacecraft, state)
        entered = count(index -> !previous_available[index] && availability.value[index], eachindex(previous_available))
        left = count(index -> previous_available[index] && !availability.value[index], eachindex(previous_available))
        scoring = @timed score_candidates!(model, policy, state, planet.μ; candidate_force=force_on_endpoint)
        matching = @timed select_interlinks!(model, policy)
        validate_walker_selection(model, policy, matching.value)
        @assert all(all(isfinite, current.pos) && all(isfinite, current.vel) for current in state.sc)
        @assert length(matching.value) <= (constellation.metadata.all_pairs ? length(spacecraft) ÷ 2 : 1)
        if !constellation.metadata.all_pairs
            @assert sum((policy.matching.score[model.graph.edge_index_lookup[key]] for key in matching.value); init=0.0) ==
                maximum((model.graph.available[connection] ? max(0.0, policy.matching.score[connection]) : 0.0 for connection in eachindex(model.graph.edges)); init=0.0)
        end
        push!(history, InterLinkHistory(time, matching.value))
        push!(metrics, (time, count(availability.value), length(matching.value), entered, left,
            availability.time, scoring.time, matching.time, availability.time + scoring.time + matching.time,
            availability.bytes, scoring.bytes, matching.bytes, availability.bytes + scoring.bytes + matching.bytes))
    end
    @assert any(metrics.selected .> 0)
    return (; metrics, history)
end

function run_walker_engine(constellation; duration::Float64=20.0, dt_max::Float64=5.0,
    wall_limit_s::Float64=Inf)
    isfinite(duration) && duration > 0 || throw(ArgumentError("duration must be finite and positive."))
    isfinite(dt_max) && dt_max > 0 || throw(ArgumentError("dt_max must be finite and positive."))
    wall_limit_s > 0 || throw(ArgumentError("wall_limit_s must be positive."))
    args = SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, generate_plots=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            mission_time=duration, orientation_sim=false, data_rate=duration),
        environment_model=make_no_gram_environment(planet=constellation.planet),
        dynamics_model=DynamicsModel(constellation.spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-11,
            dt_max_orbit=dt_max),
        solver_config=SM.SolverConfig(solver_mode=:tsit5),
        interlink_model=constellation.model, scheduling_policy_model=constellation.policy)
    initial = SimulationEngine.build_initial_conditions(args)
    mu = constellation.planet.μ
    initial_sma = [inv(2 / norm(current.pos) - dot(current.vel, current.vel) / mu) for current in initial.sc]
    wall_limit_reached = Ref(false)
    wall_start = time()
    extra_callbacks = if isfinite(wall_limit_s)
        guard = DiscreteCallback((state, current_time, integrator) -> time() - wall_start >= wall_limit_s,
            integrator -> begin
                wall_limit_reached[] = true
                terminate!(integrator)
            end; save_positions=(false, false))
        (guard,)
    else
        ()
    end
    measured = @timed withenv("SPACEAGORA_RHS_CALIBRATE" => "off",
        "SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "false", "SPACEAGORA_SOLVER_DENSE" => "false") do
        run_simulation(args; isolate_state=false, return_solution=true, save_fields=SaveField[], extra_callbacks)
    end
    solution = measured.value
    history = constellation.model.history
    completed = solution.t[end] == duration
    @assert completed || wall_limit_reached[]
    @assert length(solution.u[end].sc) == constellation.metadata.satellites
    @assert all(isfinite, solution.u[end])
    @assert length(history) == solution.destats.naccept + 1
    @assert history[1].time == 0.0 && history[end].time == solution.t[end]
    @assert all(diff([sample.time for sample in history]) .> 0)
    @assert any(!isempty(sample.active) for sample in history)
    for sample in history
        endpoints = [endpoint for key in sample.active for endpoint in key]
        @assert length(unique(endpoints)) == length(endpoints)
    end
    validate_walker_selection(constellation.model, constellation.policy, history[end].active)
    final_sma = [inv(2 / norm(current.pos) - dot(current.vel, current.vel) / mu) for current in solution.u[end].sc]
    target_index = constellation.policy.target_idx
    target = solution.u[end].sc[target_index]
    target_delta_sma = final_sma[target_index] - initial_sma[target_index]
    sma_residual = maximum(abs(final_sma[index] - initial_sma[index] - current.laser_delta_sma)
        for (index, current) in enumerate(solution.u[end].sc))
    @assert target.laser_delta_sma > 0 && target_delta_sma > 0
    @assert sma_residual < 1e-3
    initial_period = 2pi * sqrt(constellation.spacecraft[1].initial_condition.a^3 / mu)
    summary = (satellites=length(solution.u[end].sc), candidates=length(constellation.model.graph.edges),
        duration_s=duration, dtmax_s=dt_max, gravity="two_body", solver="Tsit5",
        altitude_m=constellation.metadata.altitude, all_pairs=constellation.metadata.all_pairs,
        scoring_policy=String(constellation.policy.policy), initial_period_s=initial_period,
        requested_initial_orbits=duration / initial_period, simulated_duration_s=solution.t[end],
        completed=completed, termination_reason=completed ? "completed" : "wall_time_limit",
        wall_limit_s=isfinite(wall_limit_s) ? wall_limit_s : missing,
        reltol=1e-9, abstol=1e-11, elapsed_s=measured.time, allocated_bytes=measured.bytes,
        gc_s=measured.gctime, accepted_steps=solution.destats.naccept, rejected_steps=solution.destats.nreject,
        history_samples=length(history), nonempty_schedules=count(sample -> !isempty(sample.active), history),
        target_delta_sma_m=target_delta_sma, target_integrated_laser_delta_sma_m=target.laser_delta_sma,
        target_integrated_laser_dv_norm_m_s=norm(target.laser_dv), max_sma_diagnostic_residual_m=sma_residual,
        retcode=string(solution.retcode), julia_threads=Threads.nthreads(), source_file=pathof(SpaceAGORA))
    return (; summary, history, solution, initial, initial_sma, final_sma)
end

function run_walker_engine_case(; seed::Int=20261003, duration::Float64=20.0,
    altitude::Float64=550e3, all_pairs::Bool=false, policy_name::Symbol=:gve_sma,
    wall_limit_s::Float64=Inf)
    constellation = walker_constellation(; seed, altitude, all_pairs, policy_name)
    @assert length(constellation.spacecraft) == 4000
    @assert length(constellation.model.graph.edges) == (all_pairs ? 4000 * 3999 ÷ 2 : 3999)
    println("Running 4,000-satellite force-coupled integration for ", duration,
        " s; policy=", policy_name, "; all_pairs=", all_pairs,
        "; scheduling at every accepted step; package=", pathof(SpaceAGORA))
    flush(stdout)
    engine = run_walker_engine(constellation; duration, wall_limit_s)
    suffix = altitude == 550e3 && !all_pairs && policy_name === :gve_sma ? "" :
        "_$(altitude / 1e3)km_$(all_pairs ? "all_pairs" : "target_star")_$(policy_name)"
    output = joinpath(REPO_ROOT, "III_output", "walker_interlinks", "seed_$seed", "force_coupled_$(duration)s$suffix")
    mkpath(output)
    CSV.write(joinpath(output, "constellation.csv"), DataFrame([constellation.metadata]))
    CSV.write(joinpath(output, "engine_summary.csv"), DataFrame([engine.summary]))
    CSV.write(joinpath(output, "accepted_step_schedule.csv"), schedule_table(engine.history))
    CSV.write(joinpath(output, "accepted_step_counts.csv"),
        DataFrame(time_s=[sample.time for sample in engine.history],
            selected=[length(sample.active) for sample in engine.history]))
    trajectory = DataFrame(time_s=Float64[], satellite=Int[], x_m=Float64[], y_m=Float64[], z_m=Float64[],
        vx_m_s=Float64[], vy_m_s=Float64[], vz_m_s=Float64[], mass_kg=Float64[], sma_m=Float64[],
        delta_sma_m=Float64[], laser_dvx_m_s=Float64[], laser_dvy_m_s=Float64[], laser_dvz_m_s=Float64[],
        laser_delta_sma_m=Float64[])
    mu = constellation.planet.μ
    for sample in eachindex(engine.solution.t)
        sample < lastindex(engine.solution.t) && engine.solution.t[sample] == engine.solution.t[sample + 1] && continue
        for (satellite, current) in enumerate(engine.solution.u[sample].sc)
            sma = inv(2 / norm(current.pos) - dot(current.vel, current.vel) / mu)
            push!(trajectory, (engine.solution.t[sample], satellite, current.pos..., current.vel..., current.mass,
                sma, sma - engine.initial_sma[satellite], current.laser_dv..., current.laser_delta_sma))
        end
    end
    CSV.write(joinpath(output, "trajectory.csv"), trajectory)
    @printf("Force-coupled Walker: %d accepted / %d rejected steps, target delta SMA %.9f m, wall %.3f s\n",
        engine.summary.accepted_steps, engine.summary.rejected_steps, engine.summary.target_delta_sma_m,
        engine.summary.elapsed_s)
    println("Completed requested duration: ", engine.summary.completed,
        "; simulated seconds: ", engine.summary.simulated_duration_s,
        "; termination: ", engine.summary.termination_reason)
    println("Saved force-coupled integration results to ", output)
    return (; constellation, engine, output)
end

function run_walker_case(; seed::Int=20261003)
    constellation = walker_constellation(; seed)
    @assert length(constellation.spacecraft) == 4000
    @assert length(constellation.model.graph.edges) == 3999
    output = joinpath(REPO_ROOT, "III_output", "walker_interlinks", "seed_$seed")
    mkpath(output)
    CSV.write(joinpath(output, "constellation.csv"), DataFrame([constellation.metadata]))
    elements = DataFrame([(satellite=index, plane=(index - 1) ÷ constellation.metadata.per_plane + 1,
        slot=mod(index - 1, constellation.metadata.per_plane) + 1, mass_kg=227.0, terminals=vehicle.n_terminal,
        semimajor_axis_m=vehicle.initial_condition.a, eccentricity=vehicle.initial_condition.e,
        inclination_deg=vehicle.initial_condition.i, raan_deg=vehicle.initial_condition.Ω,
        true_anomaly_deg=vehicle.initial_condition.ν) for (index, vehicle) in enumerate(constellation.spacecraft)])
    CSV.write(joinpath(output, "spacecraft.csv"), elements)
    benchmark = benchmark_walker(constellation)
    @assert benchmark.metrics.time_s == collect(0.0:60.0:3600.0)
    @assert all(benchmark.metrics.selected .== 1)
    CSV.write(joinpath(output, "snapshot_metrics.csv"), benchmark.metrics)
    CSV.write(joinpath(output, "snapshot_schedule.csv"), schedule_table(benchmark.history))
    @printf("Walker %d/%d/%d: %d snapshots, median scheduler %.6f s, median allocation %.0f bytes\n",
        constellation.metadata.satellites, constellation.metadata.planes, constellation.metadata.phasing,
        nrow(benchmark.metrics), median(benchmark.metrics.total_s), median(benchmark.metrics.allocated_bytes))
    println("Saved analytic snapshot results to ", output)
    return (; constellation, benchmark)
end

if abspath(PROGRAM_FILE) == @__FILE__
    mode = isempty(ARGS) ? "snapshot" : ARGS[1]
    if mode == "engine"
        duration = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 20.0
        run_walker_engine_case(; duration)
    elseif mode == "all_pairs"
        orbits = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 20.0
        wall_limit_s = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 120.0
        planet = make_no_gram_planet(:earth)
        duration = orbits * 2pi * sqrt((planet.Rp_e + 1000e3)^3 / planet.μ)
        run_walker_engine_case(; duration, altitude=1000e3, all_pairs=true,
            policy_name=:gve_sma_first_index, wall_limit_s)
    elseif mode == "snapshot"
        run_walker_case()
    else
        throw(ArgumentError("Mode must be snapshot, engine, or all_pairs."))
    end
end