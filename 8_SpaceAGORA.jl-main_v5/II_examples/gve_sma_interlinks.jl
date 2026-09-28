include(joinpath(@__DIR__, "common.jl"))
using LinearAlgebra
using Statistics
using Printf
using CSV
using DataFrames

function gve_sma_configuration(; dt_max_orbit::Float64=10.0, duration::Float64=3600.0,
    power::Float64=10_000.0, save_results::Bool=false, with_interlinks::Bool=true)
    planet = make_no_gram_planet(:earth)
    radii = planet.Rp_e .+ [1000e3, 1050e3, 1050e3]
    phases = [0.0, 0.018, -0.018]
    masses = [227.0, 227.0, 454.0]
    spacecraft = map(eachindex(radii)) do index
        bus = Link(root=true, m=masses[index])
        initial_condition = InitialCondition(a=radii[index], ν=rad2deg(phases[index]))
        SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia,
            initial_condition=initial_condition, id=index, n_terminal=1)
    end
    model = InterLinkModel()
    for partner in 2:length(spacecraft)
        register_candidate!(model, spacecraft, (1, 1), (partner, 1);
            parameters=InterLinkParameters(P=power, B=100.0, range=200e3))
    end
    output = joinpath(REPO_ROOT, "III_output", "gve_sma_interlinks", "dt_$(dt_max_orbit)s")
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
            dt_max_orbit=dt_max_orbit),
        solver_config=SM.SolverConfig(solver_mode=:tsit5),
        interlink_model=with_interlinks ? model : nothing,
        scheduling_policy_model=SchedulingPolicyModel(:gve_sma))
end

function interlink_switches(history)
    return [history[index] for index in 2:length(history)
        if history[index].active != history[index - 1].active]
end

function run_gve_sma_case(; dt_max_orbit::Float64=10.0, duration::Float64=3600.0,
    save_results::Bool=true)
    args = gve_sma_configuration(; dt_max_orbit, duration, save_results)
    solution = run_simulation(args; return_solution=true)
    model = solution.prob.p.args.interlink_model
    history = model.history
    intervals = diff([sample.time for sample in history])
    delta_sma = [current.laser_delta_sma for current in solution.u[end].sc]
    @printf("dtmax=%.1f s: %d accepted steps, max/median dt=%.6f/%.6f s, %d switches\n",
        dt_max_orbit, length(intervals), maximum(intervals), median(intervals), length(interlink_switches(history)))
    println("  Integrated laser semimajor-axis changes [m]: ", delta_sma)
    if save_results
        output = args.simulation_settings.results_directory
        mkpath(output)
        schedule = DataFrame(time_s=[sample.time for sample in history],
            active_connections=[join(string.(sample.active), ";") for sample in history])
        CSV.write(joinpath(output, "accepted_step_schedule.csv"), schedule)
    end
    return solution
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_gve_sma_case(dt_max_orbit=10.0)
    run_gve_sma_case(dt_max_orbit=5.0)
end