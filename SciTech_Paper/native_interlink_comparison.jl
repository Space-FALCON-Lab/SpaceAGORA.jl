using LinearAlgebra
using TOML

const PAPER_DIR = @__DIR__
length(ARGS) in (0, 2, 4) ||
    error("Usage: native runner [helper_count orbit_count [schedule target_inclination_deg]]")
const N_HELPERS = isempty(ARGS) ? 20 : parse(Int, ARGS[1])
const ORBITS = isempty(ARGS) ? 5 : parse(Int, ARGS[2])
const SCHEDULE = length(ARGS) == 4 ? Symbol(ARGS[3]) : :gve_sma
const TARGET_INC_DEG = length(ARGS) == 4 ? parse(Float64, ARGS[4]) : 0.0
N_HELPERS > 0 && ORBITS > 0 || error("Helper and orbit counts must be positive.")
SCHEDULE in (:gve_sma, :gve_inc) || error("Supported comparison schedules: gve_sma, gve_inc.")
isfinite(TARGET_INC_DEG) && 0 <= TARGET_INC_DEG <= 180 ||
    error("Target inclination must be finite and between 0 and 180 degrees.")
const CASE_DIR = isempty(ARGS) ? "" : "$(N_HELPERS)_helpers_$(ORBITS)_orbits" *
    (length(ARGS) == 4 ? "_$(SCHEDULE)_target_inc_$(TARGET_INC_DEG)_deg" : "")
const MASS_KG = 227.0
const LASER_POWER_W = 10_000.0
const LASER_MAGNIFICATION = 100.0
const LASER_RANGE_M = 200e3
const DT_MAX_S = 10.0
const EARTH_RADIUS_M = 6_378_136.6

include(joinpath(INTEGRATION_ROOT, EXAMPLES_DIR, "common.jl"))
normpath(pathof(SpaceAGORA)) == normpath(joinpath(INTEGRATION_ROOT, PACKAGE_ENTRY)) ||
    error("Loaded SpaceAGORA from the wrong project: $(pathof(SpaceAGORA))")

function sma_from_state(state, mu)
    return inv(2 / norm(state.pos) - dot(state.vel, state.vel) / mu)
end

function inclination_deg_from_state(state)
    h = cross(state.pos, state.vel)
    return rad2deg(atan(hypot(h[1], h[2]), h[3]))
end

function run_case(; output_root=PAPER_DIR, helper_count::Int=N_HELPERS)
    helper_count > 0 || error("Helper count must be positive.")
    case_dir = helper_count == N_HELPERS ? CASE_DIR :
        "$(helper_count)_helpers_$(ORBITS)_orbits_$(SCHEDULE)_target_inc_$(TARGET_INC_DEG)_deg"
    planet = make_no_gram_planet(:earth)
    helper_radius = EARTH_RADIUS_M + 1000e3
    target_radius = EARTH_RADIUS_M + 1050e3
    target_period = 2pi * sqrt(target_radius^3 / planet.μ)
    duration = ORBITS * target_period
    spacecraft = map(1:(helper_count + 1)) do index
        is_target = index == 1
        radius = is_target ? target_radius : helper_radius
        phase_deg = is_target ? 0.0 : 360.0 * (index - 2) / helper_count
        bus = Link(root=true, m=MASS_KG)
        initial_condition = InitialCondition(a=radius, e=0.0, i=is_target ? TARGET_INC_DEG : 0.0,
            Ω=0.0, ω=0.0, ν=phase_deg)
        SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia,
            initial_condition=initial_condition, id=index, n_terminal=1)
    end

    links = ConstellationInterLinkModel()
    parameters = InterLinkParameters(P=LASER_POWER_W, B=LASER_MAGNIFICATION,
        range=LASER_RANGE_M)
    for helper in 2:(helper_count + 1)
        register_allowed_link!(links, spacecraft, (1, 1), (helper, 1); parameters)
    end
    policy = SchedulingPolicyModel(SCHEDULE; target_idx=1)
    args = SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, generate_plots=false,
            results_directory=joinpath(output_root, "runs", case_dir, RESULT_CODE)),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            mission_time=duration, orientation_sim=false, data_rate=10.0),
        environment_model=make_no_gram_environment(planet=planet),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9,
            abstol_orbit=1e-11, dt_max_orbit=DT_MAX_S),
        solver_config=SM.SolverConfig(solver_mode=:tsit5),
        interlink_model=links,
        scheduling_policy_model=policy)

    elapsed = @elapsed solution = run_simulation(args; return_solution=true)
    string(solution.retcode) == "Success" || error("Simulation failed: $(solution.retcode)")
    solution.t[end] == duration || error("Simulation ended before the requested duration.")
    model = solution.prob.p.args.interlink_model
    history = model.history
    active_helper(sample) = isempty(sample.active) ? 0 : sample.active[1][2][1]
    activations = count(index -> isempty(history[index - 1].active) &&
        !isempty(history[index].active),
        2:length(history))
    final_state = solution.u[end].sc[1]
    initial_inc = inclination_deg_from_state(first(solution.u).sc[1])
    final_inc = inclination_deg_from_state(final_state)
    isapprox(initial_inc, TARGET_INC_DEG; atol=1e-10) ||
        error("Initial target inclination does not match the requested value.")
    result = (
        code=basename(INTEGRATION_ROOT),
        helper_count=helper_count,
        schedule=string(SCHEDULE),
        initial_target_altitude_m=1050e3,
        initial_helper_altitude_m=1000e3,
        initial_helper_inclination_deg=0.0,
        initial_inclination_deg=initial_inc,
        final_inclination_deg=final_inc,
        delta_inclination_deg=final_inc - initial_inc,
        requested_target_orbits=ORBITS,
        duration_s=duration,
        target_period_s=target_period,
        initial_sma_m=target_radius,
        final_sma_m=sma_from_state(final_state, planet.μ),
        delta_sma_m=sma_from_state(final_state, planet.μ) - target_radius,
        final_eccentricity=norm(cross(final_state.vel,
            cross(final_state.pos, final_state.vel)) / planet.μ -
            final_state.pos / norm(final_state.pos)),
        scheduler_updates=length(history),
        link_activations=activations,
        active_updates=count(sample -> !isempty(sample.active), history),
        final_helper=isempty(history) ? 0 : active_helper(last(history)),
        runtime_s=elapsed,
        retcode=string(solution.retcode),
    )

    result_dir = joinpath(output_root, "results", case_dir, RESULT_CODE)
    mkpath(result_dir)
    open(joinpath(result_dir, "summary.toml"), "w") do io
        TOML.print(io, Dict{String, Any}(string(k) => v for (k, v) in pairs(result)))
    end
    open(joinpath(result_dir, "schedule.csv"), "w") do io
        println(io, "time_s,target_orbits,active_helper")
        for sample in history
            println(io, sample.time, ",", sample.time / target_period, ",", active_helper(sample))
        end
    end
    println(result)
    return result
end
