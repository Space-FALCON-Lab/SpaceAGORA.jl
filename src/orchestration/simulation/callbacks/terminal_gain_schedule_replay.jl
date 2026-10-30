module TerminalGainScheduleReplay
using DiffEqCallbacks: PresetTimeCallback
using ..SimulationModel
using ..ConstellationScheduling: TerminalGainPolicy

export terminal_gain_schedule_callback, constellation_replay_configuration

"""Factory for run_simulation(...; extra_callback_factories=(terminal_gain_schedule_callback,)).

Land on every activation boundary; the existing native interlink callback then
installs that interval's matching. The existing force effectors do all forcing.
"""
function terminal_gain_schedule_callback(args)
    policy=args.scheduling_policy_model.objective
    policy isa TerminalGainPolicy || throw(ArgumentError("Expected a TerminalGainPolicy."))
    return PresetTimeCallback(policy.times[2:end],integrator->nothing;save_positions=(false,false))
end

"""Configure the native Cartesian SpaceAGORA engine to replay a fixed terminal-gain plan.

Pass the returned configuration to run_simulation with terminal_gain_schedule_callback
as an extra_callback_factory. This uses the same initial states, constants and
native links as the KS reference/replay. No orbital-rate policy is substituted.
"""
function constellation_replay_configuration(case,policy::TerminalGainPolicy)
    first(policy.times)==0 || throw(ArgumentError("The constellation example starts at t=0."))
    policy.keys==case.keys || throw(ArgumentError("Policy keys must match the configured constellation."))
    c=case.config
    planet=SimulationModel.Planets.Earth(Rp_e=c.planet.Re,μ=c.planet.μ,J2=c.planet.J2)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false,generate_plots=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            mission_time=last(policy.times),orientation_sim=false,data_rate=minimum(diff(policy.times))),
        environment_model=make_no_gram_environment(planet=planet),
        dynamics_model=DynamicsModel(deepcopy(case.spacecraft),(InverseSquaredJ2GravityModel(),)),
        guidance_model=GuidanceModel((),Float64[]),navigation_model=NavigationModel((),Float64[]),
        control_model=ControlModel((),Float64[]),initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=c.tolerance,
            abstol_orbit=c.tolerance,dt_max_orbit=minimum(diff(policy.times))/4),
        solver_config=SimulationModel.SolverConfig(solver_mode=:tsit5),
        interlink_model=deepcopy(case.interlinks),scheduling_policy_model=policy)
end
end
