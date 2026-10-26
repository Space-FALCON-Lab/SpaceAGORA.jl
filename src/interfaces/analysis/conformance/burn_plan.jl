# The burn-plan family: its representative requests, its conformance link and a short
# simulation through the generic binding (`AlgorithmGuidance`). Its validation needs no
# trusted queries.

const _BURN = BurnPlanning

"""
    burn_plan_conformance_requests() -> Vector{BurnPlanRequest}

Representative requests of the burn-plan family for [`check_algorithm`](@ref), for
one 150 kg spacecraft at perigee of a 400 by 2000 km Earth orbit with 50 kg of
propellant, a 50 N thruster and a 300 s specific impulse: an initial request at
0 s; a forced request at 100 s while an installed burn is in progress, which the
family refuses (`:burn_in_progress`) and, by its `:keep_previous` policy, keeps the
installed plan for; and a replanning request at 300 s, after that burn has ended.
The family's validation needs no trusted queries.
"""
function burn_plan_conformance_requests()
    planet = SimulationModel.make_no_gram_planet(:earth)
    state = (spacecraft_index=1, position_m=(planet.Rp_e + 400e3, 0.0, 0.0), velocity_mps=(0.0, 8200.0, 0.0),
        mass_kg=150.0, available_propellant_kg=50.0, thrust_n=50.0, isp_s=300.0, mu_m3ps2=planet.μ,
        planet_radius_m=planet.Rp_e, entry_interface_altitude_m=120e3)
    return [_BURN.BurnPlanRequest(; state..., request_id=1, reason=:initial, time_s=0.0),
        _BURN.BurnPlanRequest(; state..., request_id=2, reason=:forced, time_s=100.0, installed_window_s=(90.0, 200.0)),
        _BURN.BurnPlanRequest(; state..., request_id=3, reason=:replan, time_s=300.0, installed_window_s=(90.0, 200.0))]
end

conformance_link(::_BURN.BurnPlanFamily) =
    _BURN.BurnPlanSlot(SimulationModel.PropulsiveBurnPlan[SimulationModel.PropulsiveBurnPlan()], 1)

"""
    burn_plan_conformance_simulation(; mission_time_s=10.0) -> Function

A `simulation` for [`check_algorithm`](@ref): a function of a burn planner that
runs it through `AlgorithmGuidance`, with one request at the run's start, on one
spacecraft in a 400 by 2000 km Earth orbit whose burn scheduler (a
`BaseThrusterModel`) plans no burn of its own. The run is made twice in sequence
and, with two or more threads, twice concurrently, and the binding's run records of
the first run are returned. It throws when a run fails or the runs differ.
"""
burn_plan_conformance_simulation(; mission_time_s::Real=10.0) =
    planner -> _burn_plan_simulation_records(planner, Float64(mission_time_s))

function _burn_plan_configuration(planner, mission_time_s::Float64)
    S = SimulationModel
    planet = S.make_no_gram_planet(:earth)
    root = S.Link(root=true, m=100.0, ref_area=2.0)
    initial = S.InitialCondition(ra=planet.Rp_e + 2000e3, rp=planet.Rp_e + 400e3, i=30.0, ω=0.0, Ω=0.0, ν=10.0)
    spacecraft = S.SpacecraftModel(joints=S.Joint[], links=[root], root=root, prop_mass=50.0,
        inertia_tensor=root.inertia, initial_condition=initial, id=1)
    scheduler = S.ThrusterModels.BaseThrusterModel(thrust=[50.0], direction=[0.0], Δv=[0.0], start_burn_time=[-1.0],
        stop_burn_time=[-1.0], Isp=[300.0])
    binding = AlgorithmLifecycle.AlgorithmGuidance(; algorithm=planner, spacecraft_index=1,
        events=[AlgorithmLifecycle.AlgorithmEvent(0.0)])
    environment = S.EnvironmentModel(planet=planet, EI=120.0, density_model=S.NoAtmosphereModel(),
        ephemerides_model=S.SimpleEphemeridesModel(),
        thermal_model=S.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet), topography=false, wind=false)
    return S.SimulationConfiguration(
        simulation_settings=S.SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=S.MissionConfiguration(mission_type=S.MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=mission_time_s, orientation_sim=false, num_steps_to_save=100),
        environment_model=environment,
        dynamics_model=S.DynamicsModel(S.SpacecraftModel[spacecraft], (S.InverseSquaredGravityModel(),)),
        guidance_model=S.GuidanceModel((binding,), [1.0e6]),
        navigation_model=S.NavigationModel((), Float64[]),
        control_model=S.ControlModel(control_effectors=(scheduler,), control_rates=[1.0]),
        initial_time=S.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
end

function _burn_plan_simulation_records(planner, mission_time_s::Float64)
    configuration = _burn_plan_configuration(planner, mission_time_s)
    base = _solve(configuration)
    _failed(base) && throw(base.error)
    again = _solve(configuration)
    gap = _run_difference(base, again)
    gap === nothing || error("two runs of the burn-plan configuration differ: $(gap)")
    concurrency = _concurrency_check(configuration, base, again, FIRST_RELEASE)
    concurrency.status === :fail && error(concurrency.detail)
    return only(AlgorithmLifecycle.algorithm_run_report(base.solution).bindings).records
end
