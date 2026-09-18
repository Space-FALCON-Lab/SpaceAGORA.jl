# Accepted-step diagnostics for the optional plume effector. The RHS never
# updates this integral, so rejected stages cannot advance diagnostic mass.
const _PSI = PlumeSurfaceInteraction

function _plume_effector(args::SimulationConfiguration)
    models = filter(model -> model isa _PSI.PlumeSurfaceInteractionModel,
        args.dynamics_model.dynamic_effectors)
    length(models) <= 1 || throw(ArgumentError("At most one plume effector may publish plume diagnostics in a simulation."))
    return isempty(models) ? nothing : only(models)
end

function _plume_sample_state(model, u, p, i, t)
    engine = _simulation_engine_module()
    x = _simulation_model_module.StateSample(engine._state_position_ii(u, i),
        engine._state_velocity_ii(u, i), engine._state_mass_kg(u, p.args, i);
        q_ib=engine._state_quaternion(u, i),
        spacecraft=p.args.dynamics_model.spacecraft[i])
    env = engine.sample_environment(_PSI.environment_requirements(model), model,
        x, p, i, Float64(t))
    return x, env
end

function _accept_plume_step!(model, integrator)
    n = length(integrator.p.args.dynamics_model.spacecraft)
    _PSI._control_spacecraft_count(model.control) == n ||
        throw(ArgumentError("Plume thrust vector must match the simulation spacecraft count."))
    _PSI._validate_plume_state(model.state, n)
    for i in 1:n
        x, env = _plume_sample_state(model, integrator.u, integrator.p, i, integrator.t)
        _PSI._accept_plume_sample!(model, x, env, Float64(integrator.t), i;
            active=integrator.p.is_active[i])
    end
    return nothing
end

function get_plume_callback(args::SimulationConfiguration)
    model = _plume_effector(args)
    model === nothing && return nothing
    n = length(args.dynamics_model.spacecraft)
    _PSI._control_spacecraft_count(model.control) == n ||
        throw(ArgumentError("Plume thrust vector must match the simulation spacecraft count."))
    _PSI._validate_plume_state(model.state, n)
    # A fresh callback is built for every run. Its initialization also runs at
    # checkpoint-segment boundaries; those preserve this run's accumulated mass.
    initialized = Ref(false)
    function initialize(c, u, t, integrator)
        if !initialized[]
            _PSI._reset_plume_state!(model.state)
            initialized[] = true
        end
        _accept_plume_step!(model, integrator)
        return nothing
    end
    return DiscreteCallback((u, t, integrator) -> true,
        integrator -> _accept_plume_step!(model, integrator);
        initialize=initialize, save_positions=(false, false))
end

function _save_plume_field(num_sats, field, u, t, integrator)
    # Use the model actually owned by this solve, including isolate_state=true.
    model = _plume_effector(integrator.p.args)
    model === nothing && throw(ArgumentError("Plume save field needs a plume effector."))
    return [begin
        if field === :eroded_kg
            _PSI._plume_mass_at(model.state, i, Float64(t))
        else
            x, env = _plume_sample_state(model, u, integrator.p, i, t)
            held = _PSI._plume_saved_thrust(model.state, i, Float64(t))
            sample = _PSI._plume_sample(model, x, env, Float64(t), i; thrust_n=held)
            getproperty(sample, field)
        end
    end for i in 1:num_sats]
end

function plume_save_fields(args::SimulationConfiguration)
    _plume_effector(args) === nothing && return SaveField[]
    n = length(args.dynamics_model.spacecraft)
    return [SaveField(Symbol("plume_", name),
        (u, t, integrator) -> _save_plume_field(n, name, u, t, integrator);
        per_satellite=true, column_prefix="plume_" * String(name))
        for name in (:height_m, :pressure_pa, :shear_pa, :erosion_kg_s,
                     :eroded_kg, :ejecta_mps, :ground_effect_n)]
end
