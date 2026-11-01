include(joinpath(@__DIR__, "common.jl"))
using LinearAlgebra
using Printf

# Edge scheduling for controllability on a four-spacecraft formation: propagate
# the formation, build the feasible-edge graph of its laser links at each step,
# linearize the Cartesian two-body motion along the propagated states, and
# compare a constant link set with greedy schedules by their Gramians on the
# formation's ranges and range rates, and by the least input energy of one
# reconfiguration. The model is linearized about zero link force, so its inputs
# are the link forces themselves: the energies are those of a two-way actuator,
# and the sign check says whether a laser that only pushes could follow them.

const SMOKE = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"

function formation_configuration(; duration::Float64)
    # Step 1: Four 100 kg spacecraft near one 550 km circular orbit, spread along
    # track, the fourth also shifted across track by its RAAN.
    planet = make_no_gram_planet(:earth)
    offsets = [(ν=0.0, Ω=0.0), (ν=0.4, Ω=0.0), (ν=-0.4, Ω=0.0), (ν=0.15, Ω=0.25)]
    spacecraft = map(enumerate(offsets)) do (index, offset)
        bus = Link(root=true, m=100.0)
        initial_condition = InitialCondition(a=planet.Rp_e + 550e3, e=1e-4, i=53.0, Ω=offset.Ω, ν=offset.ν)
        SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia, initial_condition=initial_condition,
            id=index, n_terminal=1)
    end
    # Step 2: Two-body motion, saved densely enough to sample the geometry.
    args = SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, generate_plots=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime, mission_time=duration,
            orientation_sim=false, data_rate=10.0),
        environment_model=make_no_gram_environment(planet=planet),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]), navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]), initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10, abstol_orbit=1e-8, dt_max_orbit=10.0),
        solver_config=SM.SolverConfig(solver_mode=:tsit5))
    return args, planet, [100.0 for _ in spacecraft]
end

# The laser links' feasibility: a Gaussian-beam laser that the Earth and the
# atmosphere below 100 km block, closing only while it pushes with at least
# 2 mN; fore and aft terminals that steer within 70° of the along-track axis;
# and each spacecraft lasing at most four link-steps in any six steps.
function formation_feasibility(count)
    beam = GaussianBeam(wavelength_m=1.064e-6, waist_radius_m=0.2, aperture_radius_m=0.2)
    link = LaserBeamLink(power_w=10_000.0, magnification=100.0, range_m=200e3, beam=beam,
        occultation=(; min_ray_altitude_m=100e3))
    fore = OpticalTerminal(boresight=[0.0, 1.0, 0.0], field_of_regard_rad=deg2rad(70.0))
    aft = OpticalTerminal(boresight=[0.0, -1.0, 0.0], field_of_regard_rad=deg2rad(70.0))
    return EdgeFeasibilityModel(link, [[fore, aft] for _ in 1:count];
        max_range=MinimumForceRange(min_force_n=2e-3), duty_cycle=DutyCycleLimit(on_steps=4, window_steps=6))
end

# The point-mass gravity gradient ∂a/∂r at `r`.
gravity_gradient(r, μ) = (μ / norm(r)^3) .* (3.0 .* (r * r') ./ dot(r, r) - I)

# One step of the stacked Cartesian model at the step's midpoint state: each
# spacecraft's [r; v] with its own gravity gradient, and the link forces of the
# feasible-edge graph there, held over the step.
function cartesian_step(state, graph, masses, μ, dt)
    count = length(masses)
    Ac = zeros(6 * count, 6 * count)
    for i in 1:count
        block = 6 * (i - 1)
        Ac[block .+ (1:3), block .+ (4:6)] .= Matrix(1.0I, 3, 3)
        Ac[block .+ (4:6), block .+ (1:3)] .= gravity_gradient(collect(state.sc[i].pos), μ)
    end
    A, B = zoh_discretize(Ac, edge_input_matrix(graph, masses), dt)
    return EdgeLinearStep(A, B, graph.pairs)
end

function run_constellation_edge_scheduling(; steps::Int=SMOKE ? 8 : 24, dt::Float64=SMOKE ? 15.0 : 60.0,
                                           max_active::Int=2)
    # Step 1: Propagate the formation over the window.
    args, planet, masses = formation_configuration(duration=steps * dt)
    solution = run_simulation(args; return_solution=true)

    # Step 2: The feasible-edge graph at each step's midpoint.
    feasibility = formation_feasibility(length(masses))
    midpoints = [(k - 0.5) * dt for k in 1:steps]
    graphs = feasible_edge_graphs(feasibility, solution, midpoints; planet=planet)
    println("Feasible links per step: ", join([count(graph.feasible) for graph in graphs], " "))

    # Step 3: The linear model of the window, projected on the ranges and range
    # rates of every pair at its end, the rates in metres through the mean
    # motion: outputs that ignore the formation's rigid motion, which the links
    # cannot drive. Internal forces leave the barycenter alone; in orbit,
    # unequal gravity couples it to the shape over a step, so the barycenter
    # rows of B are small, not zero.
    model_steps = [cartesian_step(solution(t), graph, masses, planet.μ, dt) for (t, graph) in zip(midpoints, graphs)]
    final = solution(steps * dt)
    mean_motion = sqrt(planet.μ / norm(final.sc[1].pos)^3)
    pairs = feasibility.pairs
    model = EdgeLinearModel(model_steps, range_output_matrix(final, pairs; rate_scale=1 / mean_motion))
    transform = barycentric_transform(masses)
    barycenter_share = maximum(norm(transform.T[transform.barycenter_rows, :] * step.B) / norm(step.B)
        for step in model_steps)
    @printf("Largest share of B in the barycenter rows: %.1e\n", barycenter_share)

    # Step 4: The baselines: a constant set, greedy λ₁, greedy log det and greedy λ_min.
    epsilon = 1e-6 * tr(controllability_gramian(model)) / size(model.projection, 1)
    volume = LogDetMetric(epsilon=epsilon)
    strategies = [
        "constant set" => ConstantEdgeSet(metric=volume),
        "greedy λ₁" => GreedyEdgeSchedule(LargestEigenvalueMetric()),
        "greedy log det" => GreedyEdgeSchedule(volume),
        "greedy λ_min" => GreedyEdgeSchedule(SmallestEigenvalueMetric(); tie_break=volume),
    ]

    # Step 5: One reconfiguration: open the range between spacecraft 1 and 2 by
    # 10 m at the window's end, with every other range and range rate unchanged.
    target = zeros(size(model.projection, 1))
    target[findfirst(==((1, 2)), pairs)] = 10.0

    println(rpad("strategy", 16), rpad("link-steps", 12), rpad("tied", 6), rpad("λ₁", 12), rpad("λ_min", 12),
        rpad("log det", 10), rpad("two-way energy (N²)", 24), "push-only")
    rows = map(strategies) do (name, strategy)
        schedule = schedule_edges(model, strategy; max_active, graphs)
        for (k, active) in enumerate(schedule.active)
            admissible_edge_set(graphs[k], active; max_active, previous=schedule.active[max(1, k - 5):k-1]) ||
                error("The $name schedule broke a link limit at step $k.")
        end
        metrics = gramian_metrics(schedule.gramian; epsilon)
        energy = minimum_input_energy(model, schedule, target)
        reach = energy.reachable ? @sprintf("%.3e", energy.energy) : @sprintf("%.3e (unreached)", energy.energy)
        sign = energy.nonnegative ? "yes" : @sprintf("no (min %.1e N)", energy.most_negative_input)
        @printf("%-16s%-12d%-6d%-12.3e%-12.3e%-10.2f%-24s%s\n", name, sum(length, schedule.active),
            schedule.tied_rounds, metrics.largest_eigenvalue, metrics.smallest_eigenvalue, metrics.log_det, reach, sign)
        (; name, schedule, metrics, energy)
    end
    return (; model, graphs, rows)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_constellation_edge_scheduling()
end
