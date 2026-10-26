function get_navigation_callbacks(num_sats::Int, args::SimulationConfiguration)::Vector{DiscreteCallback}
    navigation_models = args.navigation_model.navigation_effectors
    navigation_rates = args.navigation_model.navigation_rates
    callbacks = Vector{DiscreteCallback}(undef, length(navigation_models))
    for i in eachindex(navigation_models)
        navigation_model = navigation_models[i]
        navigation_rate = navigation_rates[i]
        navigation_func = integrator -> _gnc_tick_scope!(update_navigation!, navigation_model, integrator.u,
            integrator.p, integrator.t, num_sats, _run_scope_slot(integrator.p, :navigation, i))
        callbacks[i] = _run_clock_periodic_callback(navigation_func, navigation_rate)
    end
    return callbacks
end

function get_guidance_callbacks(num_sats::Int, args::SimulationConfiguration)::Vector{DiscreteCallback}
    guidance_models = args.guidance_model.guidance_effectors
    guidance_rates = args.guidance_model.guidance_rates
    callbacks = Vector{DiscreteCallback}(undef, length(guidance_models))
    for i in eachindex(guidance_models)
        guidance_model = guidance_models[i]
        guidance_rate = guidance_rates[i]
        guidance_func = integrator -> _gnc_tick_scope!(update_guidance!, guidance_model, integrator.u,
            integrator.p, integrator.t, num_sats, _run_scope_slot(integrator.p, :guidance, i))
        callbacks[i] = _run_clock_periodic_callback(guidance_func, guidance_rate)
    end
    return callbacks
end
