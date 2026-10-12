if !isdefined(@__MODULE__, :REPO_ROOT)
    include(joinpath(@__DIR__, "common.jl"))
end
using Random
using Statistics
using Printf
using CSV
using DataFrames

function walker_constellation(; satellites::Int=4000, planes::Int=80,
    seed::Int=20261003, altitude::Float64=550e3, inclination::Float64=53.0,
    target_idx::Int=1, range_m::Float64=2e6)
    satellites > 1 || throw(ArgumentError("satellites must exceed one."))
    1 <= planes <= satellites && satellites % planes == 0 ||
        throw(ArgumentError("planes must be a positive divisor of satellites."))
    1 <= target_idx <= satellites || throw(ArgumentError("target_idx is outside the constellation."))
    isfinite(altitude) && altitude > 0 || throw(ArgumentError("altitude must be finite and positive."))
    isfinite(inclination) && 0 <= inclination <= 180 ||
        throw(ArgumentError("inclination must be in [0, 180] degrees."))
    parameters = InterLinkParameters(range=range_m)
    rng = MersenneTwister(seed)
    phasing = rand(rng, 0:(planes - 1))
    raan_offset, anomaly_offset = 360rand(rng), 360rand(rng)
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
    model = InterLinkModel()
    for partner in eachindex(spacecraft)
        partner == target_idx && continue
        register_candidate!(model, spacecraft, (target_idx, 1), (partner, 1); parameters)
    end
    policy = SchedulingPolicyModel(:gve_sma; target_idx)
    metadata = (; satellites, planes, per_plane, phasing, seed, altitude, inclination,
        raan_offset, anomaly_offset, range_m, target_idx)
    return (; spacecraft, planet, model, policy, metadata)
end

function walker_snapshot(constellation)
    return (sc=map(constellation.spacecraft) do vehicle
        pos, vel = SimulationEngine.orbitalelemtorv(vehicle.initial_condition, constellation.planet)
        (; pos, vel, mass=227.0)
    end,)
end

function advance_walker_snapshot!(state, initial, mean_motion, time)
    sine, cosine = sincos(mean_motion * time)
    for index in eachindex(state.sc)
        current, origin = state.sc[index], initial.sc[index]
        @. current.pos = cosine * origin.pos + (sine / mean_motion) * origin.vel
        @. current.vel = -mean_motion * sine * origin.pos + cosine * origin.vel
    end
    return state
end

function validate_walker_selection(model, selected)
    occupied = Set{TerminalEndpoint}()
    for key in selected
        connection = model.linkgraph[key]
        connection.state.available && connection.state.active ||
            error("Selected Walker link is not available and active: $key")
        connection.state.score > model.active_link_penalty ||
            error("Selected Walker link has nonpositive net benefit: $key")
        for endpoint in key
            endpoint in occupied && error("Walker schedule reuses terminal $endpoint")
            push!(occupied, endpoint)
        end
    end
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
    schedule_interlinks!(model, policy, spacecraft, state, planet.μ)
    metrics = DataFrame(time_s=Float64[], available=Int[], selected=Int[],
        availability_s=Float64[], scoring_s=Float64[], matching_s=Float64[],
        total_s=Float64[], allocated_bytes=Int[])
    history = InterLinkScheduleSample[]
    for time in 0.0:interval:duration
        advance_walker_snapshot!(state, initial, mean_motion, time)
        availability = @timed update_availability!(model, spacecraft, state)
        scoring = @timed score_candidates!(model, policy, state, planet.μ)
        matching = @timed select_interlinks!(model)
        validate_walker_selection(model, matching.value)
        push!(history, InterLinkScheduleSample(time, matching.value))
        push!(metrics, (time, count(model.available), length(matching.value),
            availability.time, scoring.time, matching.time,
            availability.time + scoring.time + matching.time,
            availability.bytes + scoring.bytes + matching.bytes))
    end
    any(metrics.selected .> 0) || error("Walker snapshot experiment produced no selected links.")
    return (; metrics, history)
end

function walker_configuration(constellation; duration::Float64=20.0, dt_max::Float64=5.0)
    isfinite(duration) && duration > 0 || throw(ArgumentError("duration must be finite and positive."))
    isfinite(dt_max) && dt_max > 0 || throw(ArgumentError("dt_max must be finite and positive."))
    return SimulationConfiguration(
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
end

function run_walker_engine(constellation; duration::Float64=20.0)
    args = walker_configuration(constellation; duration)
    measured = @timed run_simulation(args; return_solution=true)
    solution = measured.value
    model = solution.prob.p.args.interlink_model
    history = model.history
    solution.t[end] == duration || error("Walker integration did not reach its final time.")
    length(solution.u[end].sc) == constellation.metadata.satellites ||
        error("Walker integration has the wrong spacecraft count.")
    all(isfinite, solution.u[end]) || error("Walker integration produced nonfinite state.")
    length(history) == solution.destats.naccept + 1 ||
        error("Walker scheduler history does not match accepted steps.")
    all(diff([sample.time for sample in history]) .> 0) ||
        error("Walker scheduler history times are not strictly increasing.")
    any(!isempty(sample.active) for sample in history) ||
        error("Walker integration produced no selected links.")
    for sample in history
        endpoints = [endpoint for key in sample.active for endpoint in key]
        length(unique(endpoints)) == length(endpoints) ||
            error("Walker accepted-step schedule reuses a terminal.")
    end
    target = solution.u[end].sc[constellation.policy.target_idx]
    target.laser_delta_sma > 0 || error("Walker target has no positive integrated laser benefit.")
    validate_walker_selection(model, history[end].active)
    summary = (satellites=length(solution.u[end].sc), duration_s=duration,
        elapsed_s=measured.time, allocated_bytes=measured.bytes,
        accepted_steps=solution.destats.naccept, rejected_steps=solution.destats.nreject,
        history_samples=length(history), target_delta_sma_m=target.laser_delta_sma)
    return (; summary, history)
end

function run_walker_case(; seed::Int=20261003)
    constellation = walker_constellation(; seed)
    output = joinpath(REPO_ROOT, "III_output", "walker_interlinks", "seed_$seed")
    mkpath(output)
    CSV.write(joinpath(output, "constellation.csv"), DataFrame([constellation.metadata]))
    benchmark = benchmark_walker(constellation)
    CSV.write(joinpath(output, "snapshot_metrics.csv"), benchmark.metrics)
    CSV.write(joinpath(output, "snapshot_schedule.csv"), schedule_table(benchmark.history))
    @printf("Walker %d/%d/%d: median scheduler %.6f s, median allocation %.0f bytes\n",
        constellation.metadata.satellites, constellation.metadata.planes, constellation.metadata.phasing,
        median(benchmark.metrics.total_s), median(benchmark.metrics.allocated_bytes))
    engine = run_walker_engine(constellation)
    CSV.write(joinpath(output, "engine_summary.csv"), DataFrame([engine.summary]))
    CSV.write(joinpath(output, "accepted_step_schedule.csv"), schedule_table(engine.history))
    CSV.write(joinpath(output, "accepted_step_counts.csv"),
        DataFrame(time_s=[sample.time for sample in engine.history],
            selected=[length(sample.active) for sample in engine.history]))
    println("Walker full-engine validation: ", engine.summary)
    println("Saved Walker results to ", output)
    return (; constellation, benchmark, engine)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_walker_case()
end
