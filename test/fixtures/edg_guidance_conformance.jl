# The EDG guidance family's conformance fixture: requests for one three-body
# spacecraft above the 160 km entry interface of an exponential Earth atmosphere,
# and a short run through `AlgorithmGuidance`. Proposals need the SpaceAGORAEDG
# provider. Included by conformance_instances.jl.

# The analytic EDG run of the optional-package tests: a circular 120 km orbit inside
# the atmosphere, guidance and control every 1/3 s. `guidance` is `:model` for
# `AerobrakingEnergyDepletionGuidanceModel` or `:binding` for an `AlgorithmGuidance`
# of the EDG guidance family with `algorithm`.
function edg_guidance_configuration(; guidance::Symbol=:binding, algorithm=SpaceAGORA.EDGGuidanceAlgorithm(),
        duration_s::Real=1.0, guidance_modes=(:max_energy_depletion,), results::Bool=false)
    planet = SpaceAGORA.make_no_gram_planet(:earth)
    spacecraft = SpaceAGORA.make_three_body_spacecraft(
        bus_dims=(2.0, 2.0, 2.0), panel_dims=(0.02, 2.0, 1.0),
        bus_mass=500.0, panel_mass_each=10.0, panel_offset_y=2.0,
        ic=SpaceAGORA.InitialCondition(ra=planet.Rp_e + 120e3, rp=planet.Rp_e + 120e3,
            i=0.0, ω=0.0, Ω=0.0, ν=0.0), prop_mass=0.0, id=1)
    config = SpaceAGORA.AerobrakingEnergyDepletionConfig(guidance_modes=guidance_modes,
        max_energy_submodes=(:heat_rate, :structural_load), controlled_panel_links=(2, 3),
        heat_rate_limit_w_cm2=0.15, structural_load_limit_pa=0.5, min_alpha_rad=1e-4, max_alpha_rad=pi / 2)
    state = SpaceAGORA.AerobrakingEnergyDepletionState(num_sats=1)
    model = guidance === :model ? SpaceAGORA.AerobrakingEnergyDepletionGuidanceModel(config, state) :
        SpaceAGORA.AlgorithmGuidance(; algorithm=algorithm, spacecraft_index=1,
            settings=SpaceAGORA.EDGGuidanceSettings(config=config, state=state), replan_interval_s=1 / 3)
    control = SpaceAGORA.AerobrakingEnergyDepletionControlModel(config, state)
    return SpaceAGORA.SimulationConfiguration(
        simulation_settings=SpaceAGORA.SimulationSettings(results=results, verbose=false, generate_plots=false),
        mission_configuration=SpaceAGORA.MissionConfiguration(mission_type=SpaceAGORA.MissionTime,
            mission_time=Float64(duration_s), keplerian=false, orientation_sim=false, num_steps_to_save=25),
        environment_model=SpaceAGORA.EnvironmentModel(planet=planet, EI=160.0,
            density_model=SpaceAGORA.ExponentialAtmosphereModel(1e-8, 120e3, 20e3; temperature_k=200.0),
            ephemerides_model=SpaceAGORA.SimpleEphemeridesModel(), wind=false, topography=false,
            thermal_model=SpaceAGORA.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet)),
        dynamics_model=SpaceAGORA.DynamicsModel([spacecraft],
            (SpaceAGORA.InverseSquaredGravityModel(), SpaceAGORA.AerodynamicCoefficientfM())),
        guidance_model=SpaceAGORA.GuidanceModel((model,), [1 / 3]),
        control_model=SpaceAGORA.ControlModel((control,), [1 / 3]),
        initial_time=SpaceAGORA.InitialTime(year=2024, month=1, day=1),
        solver_config=SpaceAGORA.SolverConfig(solver_mode=:tsit5, parallel=false))
end

"""
Requests of the EDG guidance family for the conformance checks, each from a fresh
EDG state of the fixture above with its spacecraft at 300 km, outside the
atmosphere: a maximum-depletion request at 10 s, whose decision selects that mode,
and a targeting request at 20 s, which leaves the spacecraft inactive because the
targeting plan is made only inside a drag passage. Each installed decision
therefore changes the state it is installed in.
"""
function edg_guidance_conformance_requests()
    args = edg_guidance_configuration()
    p = SM.ODEParams(n_sats=1, args=args)
    planet = args.environment_model.planet
    u = (sc=[(pos=[planet.Rp_e + 300e3, 0.0, 0.0], vel=[0.0, 7700.0, 0.0], mass=520.0, heat_loads=zeros(3))],)
    settings(modes) = SpaceAGORA.EDGGuidanceSettings(
        config=SpaceAGORA.AerobrakingEnergyDepletionConfig(guidance_modes=modes, controlled_panel_links=(2, 3)),
        state=SpaceAGORA.AerobrakingEnergyDepletionState(num_sats=1))
    family = SpaceAGORA.EDGGuidanceFamily()
    return [SpaceAGORA.family_request(family, settings((:max_energy_depletion,)), u, p, 10.0, 1, :interval, 1),
        SpaceAGORA.family_request(family, settings((:targeting, :max_energy_depletion)), u, p, 20.0, 1, :interval, 2)]
end

"""
A `simulation` for `check_algorithm`: the fixture's run with `algorithm` through
`AlgorithmGuidance`, returning the binding's run records.
"""
function edg_guidance_conformance_simulation()
    return algorithm -> begin
        result = SpaceAGORA.run_simulation(edg_guidance_configuration(; algorithm=algorithm);
            return_results=true, return_solution=true)
        only(SpaceAGORA.algorithm_run_report(result.solution).bindings).records
    end
end
