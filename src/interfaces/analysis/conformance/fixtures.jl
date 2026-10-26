# The fixture: a short no-GRAM run built from public constructors, and the runs
# and comparisons the behavioural checks make on it.

"""
    ConformanceFixture

A configuration the checks run, with the concrete argument types its run passes
to the hooks. Build it with [`conformance_fixture`](@ref).

- `configuration`: the `SimulationConfiguration` every check runs. Checks run it
  with `isolate_state=true`, except `deepcopy_equivalent`, which runs it once with
  `isolate_state=false` after the other runs.
- `arguments`: the configuration the engine's pre-run hook check reads, after the
  same per-spacecraft GNC flattening and density epoch alignment `run_simulation`
  applies first.
- `call`: the run's state and parameter types, as the pre-run check derives them
  from `build_initial_conditions` and `ODEParams`.
- `solver_mode`: the solver mode of the run.
"""
struct ConformanceFixture{C, A}
    configuration::C
    arguments::A
    call::NamedTuple{(:state, :params)}
    solver_mode::Symbol
end

function ConformanceFixture(configuration::SimulationModel.SimulationConfiguration)
    arguments = SimulationEngine._with_density_model_epoch(SimulationEngine._flatten_spacecraft_gnc(configuration))
    solver = arguments.solver_config === nothing ? SimulationEngine._solver_config_from_env() : arguments.solver_config
    u0 = SimulationEngine.build_initial_conditions(arguments)
    params = SimulationModel.ODEParams(n_sats=length(arguments.dynamics_model.spacecraft), args=arguments)
    return ConformanceFixture(configuration, arguments, SimulationEngine._hook_call_types(u0, params), solver.solver_mode)
end

function Base.show(io::IO, fixture::ConformanceFixture)
    dynamics = fixture.configuration.dynamics_model
    print(io, "ConformanceFixture(", length(dynamics.spacecraft), " spacecraft, ",
        length(dynamics.dynamic_effectors), " force/torque models, ",
        fixture.configuration.mission_configuration.mission_time, " s, ", fixture.solver_mode, ")")
end

"""
    conformance_fixture(; configuration=nothing, spacecraft_count=1, mission_time_s=60.0,
        guidance=(), navigation=(), control=(), effectors=(), sensors=(), guidance_rate_s=10.0,
        navigation_rate_s=10.0, control_rate_s=1.0, orientation=true) -> ConformanceFixture

The run the checks use. Without `configuration` it is a short Earth orbit that needs
no native GRAM, SPICE kernels or private data: `make_no_gram_planet(:earth)`,
`NoAtmosphereModel()`, `SimpleEphemeridesModel()`, `spacecraft_count` single-link
spacecraft (100 kg, 2 m², on a 650 km by 700 km orbit at 45°, spaced 5° apart in true
anomaly), a `MissionTime` of `mission_time_s` seconds, attitude propagation when
`orientation` is `true`, `InverseSquaredJ2GravityModel()` followed by `effectors` as
the force/torque models, `sensors` mounted on the first spacecraft's root link with
`add_sensor!`, and the given guidance, navigation and control models at their rates
in seconds.

A model that needs more (an atmosphere for drag, a second spacecraft with a
particular orbit, articulated links) is checked on a `configuration` the author
builds, which must already contain the models under test; the component keywords
are then refused.
"""
function conformance_fixture(; configuration::Union{Nothing, SimulationModel.SimulationConfiguration}=nothing,
        spacecraft_count::Integer=1, mission_time_s::Real=60.0, guidance::Tuple=(), navigation::Tuple=(),
        control::Tuple=(), effectors::Tuple=(), sensors::Tuple=(), guidance_rate_s::Real=10.0,
        navigation_rate_s::Real=10.0, control_rate_s::Real=1.0, orientation::Bool=true)
    if configuration !== nothing
        isempty(guidance) && isempty(navigation) && isempty(control) && isempty(effectors) && isempty(sensors) || throw(ArgumentError(
            "conformance_fixture: pass the models under test either inside `configuration` or as keywords, not both."))
        return ConformanceFixture(configuration)
    end
    return ConformanceFixture(_fixture_configuration(; spacecraft_count, mission_time_s, guidance, navigation,
        control, effectors, sensors, guidance_rate_s, navigation_rate_s, control_rate_s, orientation))
end

# The built-in fixture configuration. `density_model`, `planet` and
# `ephemerides_model` are for core's own instances of the built-in models that need
# an atmosphere or SPICE ephemerides.
function _fixture_configuration(; spacecraft_count::Integer=1, mission_time_s::Real=60.0, guidance::Tuple=(),
        navigation::Tuple=(), control::Tuple=(), effectors::Tuple=(), sensors::Tuple=(), guidance_rate_s::Real=10.0,
        navigation_rate_s::Real=10.0, control_rate_s::Real=1.0, orientation::Bool=true,
        density_model=SimulationModel.NoAtmosphereModel(), planet=SimulationModel.make_no_gram_planet(:earth),
        ephemerides_model=SimulationModel.SimpleEphemeridesModel())
    spacecraft_count >= 1 || throw(ArgumentError("conformance_fixture: spacecraft_count must be at least 1."))
    mission_time_s > 0 || throw(ArgumentError("conformance_fixture: mission_time_s must be positive."))
    spacecraft = map(1:spacecraft_count) do id
        root = SimulationModel.Link(root=true, m=100.0, ref_area=2.0)
        orbit = SimulationModel.InitialCondition(ra=planet.Rp_e + 700e3, rp=planet.Rp_e + 650e3, i=45.0,
            ω=0.0, Ω=0.0, ν=10.0 + 5.0 * (id - 1))
        SimulationModel.SpacecraftModel(joints=SimulationModel.Joint[], links=[root], root=root, prop_mass=0.0,
            inertia_tensor=root.inertia, initial_condition=orbit, id=id)
    end
    first_spacecraft = first(spacecraft)
    for sensor in sensors
        SimulationModel.add_sensor!(first_spacecraft, first_spacecraft.root, sensor)
    end
    environment = SimulationModel.EnvironmentModel(planet=planet, EI=120.0, density_model=density_model,
        ephemerides_model=ephemerides_model,
        thermal_model=SimulationModel.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false, wind=false)
    return SimulationModel.SimulationConfiguration(
        simulation_settings=SimulationModel.SimulationSettings(results=false, verbose=false, generate_plots=false,
            normalize=false),
        mission_configuration=SimulationModel.MissionConfiguration(mission_type=SimulationModel.MissionTime,
            keplerian=true, number_of_orbits=1, mission_time=Float64(mission_time_s), orientation_sim=orientation,
            num_steps_to_save=100),
        environment_model=environment,
        dynamics_model=SimulationModel.DynamicsModel(SimulationModel.SpacecraftModel[spacecraft...],
            (SimulationModel.InverseSquaredJ2GravityModel(), effectors...)),
        solver_config=SimulationEngine.SolverConfig(),
        guidance_model=SimulationModel.GuidanceModel(guidance, fill(Float64(guidance_rate_s), length(guidance))),
        navigation_model=SimulationModel.NavigationModel(navigation, fill(Float64(navigation_rate_s), length(navigation))),
        control_model=SimulationModel.ControlModel(control, fill(Float64(control_rate_s), length(control))),
        initial_time=SimulationModel.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
end

# A finished run: its solution, or the error it stopped with.
struct _RunOutcome
    solution::Any
    error::Any
end

_failed(run::_RunOutcome) = run.error !== nothing

function _solve(configuration; isolate_state::Bool=true, kwargs...)
    try
        return _RunOutcome(SimulationEngine.run_simulation(configuration; isolate_state, return_solution=true, kwargs...), nothing)
    catch err
        return _RunOutcome(nothing, err)
    end
end

# How two runs differ, or `nothing` when their saved times and states (or, for a
# sensor, its measurement times and measurements) are equal bit for bit (`isequal`
# on every element).
function _run_difference(a::_RunOutcome, b::_RunOutcome; noun::AbstractString="state")
    _failed(a) && return "the first run failed: " * _summary(a.error)
    _failed(b) && return "the second run failed: " * _summary(b.error)
    ta, tb = a.solution.t, b.solution.t
    length(ta) == length(tb) || return "the runs recorded $(length(ta)) and $(length(tb)) $(noun)s"
    for k in eachindex(ta, tb)
        isequal(ta[k], tb[k]) || return "the time of $(noun) $k differs: $(ta[k]) s and $(tb[k]) s"
        ua, ub = a.solution.u[k], b.solution.u[k]
        if !isequal(ua, ub)
            gap = length(ua) == length(ub) ? maximum(abs.(collect(ua) .- collect(ub))) : NaN
            return "the $(noun) at t = $(ta[k]) s differs (largest element difference $(gap))"
        end
    end
    return nothing
end
