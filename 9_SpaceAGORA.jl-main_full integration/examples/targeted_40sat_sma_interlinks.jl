if !isdefined(@__MODULE__, :REPO_ROOT)
    include(joinpath(@__DIR__, "common.jl"))
end
using LinearAlgebra
using Printf
import Plots

# The graph is schematic, not orbital geometry; animation sampling can skip brief links.
function targeted_40sat_configuration(; save_results::Bool=true, orbits::Float64=20.0)
    isfinite(orbits) && orbits > 0 || throw(ArgumentError("orbits must be finite and positive."))
    planet = make_no_gram_planet(:earth)
    helper_indices = 1:30
    target_indices = 31:40
    helper_altitude = 1000e3
    target_altitude = 1100e3
    target_phase_offset = 1.0
    spacecraft = map(1:40) do index
        ring = index in helper_indices ? helper_indices : target_indices
        ring_index = index - first(ring)
        altitude = index in helper_indices ? helper_altitude : target_altitude
        phase_offset = index in target_indices ? target_phase_offset : 0.0
        bus = Link(root=true, m=227.0)
        initial_condition = InitialCondition(a=planet.Rp_e + altitude, e=0.0,
            i=0.0, Ω=0.0, ν=360.0 * ring_index / length(ring) + phase_offset)
        SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia,
            initial_condition=initial_condition, id=index, n_terminal=1)
    end

    interlinks = ConstellationInterLinkModel()
    parameters = InterLinkParameters(P=10_000.0, B=100.0, range=200e3)
    for helper in helper_indices, target in target_indices
        register_allowed_link!(interlinks, spacecraft, (helper, 1), (target, 1); parameters)
    end
    policy = SchedulingPolicyModel(:gve_sma; target_idxs=collect(target_indices))

    target_period = 2pi * sqrt((planet.Rp_e + target_altitude)^3 / planet.μ)
    duration = orbits * target_period
    output = joinpath(REPO_ROOT, "output", "targeted_40sat_sma_interlinks")
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=save_results, generate_plots=false,
            results_directory=output, save_csv=true),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            mission_time=duration, orientation_sim=false, data_rate=10.0),
        environment_model=make_no_gram_environment(planet=planet),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-11,
            dt_max_orbit=10.0),
        solver_config=SM.SolverConfig(solver_mode=:tsit5),
        interlink_model=interlinks,
        scheduling_policy_model=policy)
end

function interlink_animation_frames(solution; frame_count::Int=600)
    frame_count >= 2 || throw(ArgumentError("frame_count must be at least two."))
    args = solution.prob.p.args
    model = deepcopy(args.interlink_model)
    history = model.history
    isempty(history) && error("No interlink schedule was recorded.")
    history_times = [sample.time for sample in history]
    indices = unique(searchsortedlast.(Ref(history_times),
        range(first(history_times), last(history_times); length=frame_count)))
    times = history_times[indices]
    endpoints = [map(v -> model.graph.vertices[v][1], edge) for edge in model.graph.edges]
    status = map(indices) do index
        sample = history[index]
        saved = searchsortedlast(solution.t, sample.time)
        saved > 0 && solution.t[saved] == sample.time ||
            error("Animation needs saved accepted steps; rerun with animate=true.")
        available = update_availability!(model, args.dynamics_model.spacecraft, solution.u[saved])
        active = Set(sample.active)
        [minmax(model.graph.vertices[edge[1]], model.graph.vertices[edge[2]]) in active ?
            2 : Int(available[k]) for (k, edge) in enumerate(model.graph.edges)]
    end
    target = first(args.scheduling_policy_model.target_idxs)
    period = 2pi * sqrt(args.dynamics_model.spacecraft[target].initial_condition.a^3 /
        args.environment_model.planet.μ)
    return (; times, status, endpoints, period)
end

function save_interlink_animation(solution; output::String, frame_count::Int=600, fps::Int=20)
    fps > 0 || throw(ArgumentError("fps must be positive."))
    frames = interlink_animation_frames(solution; frame_count)
    targets = solution.prob.p.args.scheduling_policy_model.target_idxs
    helpers = setdiff(eachindex(solution.prob.p.args.dynamics_model.spacecraft), targets)
    x, y = zeros(40), zeros(40)
    for (ring, row) in ((helpers, 350.0), (targets, 190.0))
        for (k, index) in enumerate(ring)
            x[index] = 70.0 + 940.0 * (k - 1) / (length(ring) - 1)
            y[index] = row + (isodd(k) ? -30 : 30)
        end
    end
    colors = ["#cbd5e1", "#16a34a", "#f97316"]
    labels = ["Unavailable", "Available (idle)", "Active"]
    mkpath(output)
    gif_path = joinpath(output, "interlinks.gif")
    mktempdir() do directory
        animation = Plots.Animation(directory, String[])
        for (time, status) in zip(frames.times, frames.status)
            figure = Plots.plot(; size=(1080, 560), xlims=(0, 1080), ylims=(0, 560),
                yflip=true, axis=false, ticks=false, legend=false, grid=false,
                background_color=:white)
            for kind in 0:2
                xs, ys = Float64[], Float64[]
                for (edge, (first, second)) in enumerate(frames.endpoints)
                    status[edge] == kind || continue
                    append!(xs, (x[first], x[second], NaN))
                    append!(ys, (y[first], y[second], NaN))
                end
                Plots.plot!(figure, xs, ys; color=colors[kind + 1],
                    linewidth=kind == 0 ? 0.4 : 2.0, alpha=kind == 0 ? 0.3 : 1.0)
                Plots.annotate!(figure, 200 + 340kind, 505,
                    Plots.text("$(labels[kind + 1]): $(count(==(kind), status))", 11, colors[kind + 1]))
            end
            Plots.scatter!(figure, x, y; markersize=11, markerstrokewidth=0,
                color=[index in helpers ? "#0369a1" : "#4338ca" for index in 1:40])
            for index in 1:40
                Plots.annotate!(figure, x[index], y[index], Plots.text(string(index), 8, :white))
            end
            Plots.annotate!(figure, 540, 40, Plots.text(
                @sprintf("t = %.1f s | %.3f target orbits", time, time / frames.period), 14))
            Plots.annotate!(figure, 540, 100, Plots.text("1100 km targets (31-40)", 12))
            Plots.annotate!(figure, 540, 440, Plots.text("1000 km helpers (1-30)", 12))
            Plots.frame(animation, figure)
        end
        Plots.gif(animation, gif_path; fps)
    end
    println("Saved animation to ", gif_path)
    return gif_path
end

function run_targeted_40sat_example(; save_results::Bool=true, orbits::Float64=20.0,
    animate::Bool=true, frame_count::Int=600, fps::Int=20)
    frame_count >= 2 || throw(ArgumentError("frame_count must be at least two."))
    fps > 0 || throw(ArgumentError("fps must be positive."))
    args = targeted_40sat_configuration(; save_results, orbits)
    solution = if animate
        withenv("SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => "true") do
            run_simulation(args; return_solution=true)
        end
    else
        run_simulation(args; return_solution=true)
    end
    solution.t[end] == args.mission_configuration.mission_time ||
        error("Simulation stopped before the requested duration: $(solution.retcode)")
    model = solution.prob.p.args.interlink_model
    mu = args.environment_model.planet.μ
    initial_sma = [vehicle.initial_condition.a for vehicle in args.dynamics_model.spacecraft]
    targets = args.scheduling_policy_model.target_idxs
    target_delta_sma = [
        inv(2 / norm(solution.u[end].sc[index].pos) -
            dot(solution.u[end].sc[index].vel, solution.u[end].sc[index].vel) / mu) -
        initial_sma[index] for index in targets
    ]
    target_orbit_count = solution.t[end] / (2pi * sqrt(initial_sma[first(targets)]^3 / mu))
    active_samples = count(sample -> !isempty(sample.active), model.history)
    @printf("Simulated %.2f target orbits; %d accepted steps; %d/%d samples with active links\n",
        target_orbit_count, solution.stats.naccept, active_samples, length(model.history))
    @printf("Target SMA changes [m] (satellites 31-40): %s\n", target_delta_sma)
    if animate
        save_interlink_animation(solution; output=args.simulation_settings.results_directory,
            frame_count, fps)
    end
    return solution
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_targeted_40sat_example()
end
