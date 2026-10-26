# Worked example D of "Writing a GNC Algorithm" (docs/src/user/gnc_planning.md):
# a planner of the burn-plan family, run through the generic binding
# `AlgorithmGuidance`. Each region between a `# snippet:` line and `# end snippet`
# appears on that page as it is here; the HPC and extensibility docs gate compares them.
using SpaceAGORA
using SpaceAGORA.Conformance: test_algorithm, burn_plan_conformance_requests, burn_plan_conformance_simulation
using Test

# snippet: planner
# A prograde burn of a set Δv that starts a set time after the request.
Base.@kwdef struct BurnSoonPlanner <: AbstractBurnPlanner
    delta_v_mps::Float64 = 5.0      # the burn's Δv, m/s
    lead_s::Float64 = 10.0          # time from the request to the burn's start, s
end
SpaceAGORA.consumes_rng(::BurnSoonPlanner) = false

function SpaceAGORA.propose!(state, planner::BurnSoonPlanner, request::BurnPlanRequest, rng)
    # The scheduler's own size check, which the family's validation repeats.
    sizing = burn_plan_sizing(request.mass_kg, planner.delta_v_mps, 0.0, request.thrust_n, request.isp_s;
                              available_propellant_kg=request.available_propellant_kg)
    sizing.reason === :accepted || return AlgorithmResult(request_id=request.request_id,
        status=:infeasible, termination=sizing.reason)
    start_s = request.time_s + planner.lead_s
    plan = PropulsiveBurnPlan(valid=true, delta_v_mps=planner.delta_v_mps, direction_rad=0.0,
        source_orbit=-1, start_burn_s=start_s, stop_burn_s=start_s + sizing.burn_duration_s,
        thrust_n=request.thrust_n, isp_s=request.isp_s,
        commanded_impulse_n_s=sizing.commanded_impulse_n_s,
        propellant_required_kg=sizing.propellant_required_kg)
    return AlgorithmResult(request_id=request.request_id, status=:candidate,
        termination=:completed, product=plan)
end
# end snippet

# snippet: run
# The family installs into the burn scheduler, which plans no burn of its own here
# (zero Δv, no window) and fires a 50 N engine of 300 s specific impulse.
planet = make_no_gram_planet(:earth)
root = Link(root=true, m=100.0, ref_area=2.0)
orbit = InitialCondition(ra=planet.Rp_e + 2000e3, rp=planet.Rp_e + 400e3, i=30.0, ω=0.0, Ω=0.0, ν=10.0)
spacecraft = SpacecraftModel(joints=Joint[], links=[root], root=root, prop_mass=50.0,
                             inertia_tensor=root.inertia, initial_condition=orbit, id=1)
scheduler = BaseThrusterModel(thrust=[50.0], direction=[0.0], Δv=[0.0], start_burn_time=[-1.0],
                              stop_burn_time=[-1.0], Isp=[300.0])
binding = AlgorithmGuidance(algorithm=BurnSoonPlanner(), spacecraft_index=1,
                            events=[AlgorithmEvent(0.0)])      # one request, at the start
config = SimulationConfiguration(
    simulation_settings=SimulationSettings(
        results=false, verbose=false, generate_plots=false, normalize=false),
    mission_configuration=MissionConfiguration(
        mission_type=MissionTime, mission_time=60.0, number_of_orbits=1,
        keplerian=true, orientation_sim=false, num_steps_to_save=100),
    environment_model=EnvironmentModel(
        planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false, wind=false),
    dynamics_model=DynamicsModel([spacecraft], (InverseSquaredGravityModel(),)),
    guidance_model=GuidanceModel(guidance_effectors=(binding,), guidance_rates=[60.0]),
    navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
    control_model=ControlModel(control_effectors=(scheduler,), control_rates=[1.0]),
    initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
    integration_tolerances=IntegrationTolerances())

solution = run_simulation(config; return_solution=true)
records = only(algorithm_run_report(solution).bindings).records
plan = only(filter(record -> record.event === :installed, records)).product

# By hand: the engine takes the propellant the family validated, and that is the
# rocket equation's mass for the burn's Δv.
used_kg = (100.0 + 50.0) - solution.u[end].sc[1].mass
rocket_kg = (100.0 + 50.0) * (1 - exp(-plan.delta_v_mps / (plan.isp_s * 9.80665)))
println("burn from $(plan.start_burn_s) s to $(round(plan.stop_burn_s, digits=3)) s; propellant used ",
        "$(round(used_kg, digits=6)) kg, planned $(round(plan.propellant_required_kg, digits=6)) kg, ",
        "rocket equation $(round(rocket_kg, digits=6)) kg")
@test used_kg ≈ plan.propellant_required_kg rtol=1e-6
@test plan.propellant_required_kg ≈ rocket_kg rtol=1e-9
# end snippet

# snippet: conformance
@testset "BurnSoonPlanner conformance" begin
    test_algorithm(BurnSoonPlanner(), burn_plan_conformance_requests();
                   simulation=burn_plan_conformance_simulation())
end
# end snippet
