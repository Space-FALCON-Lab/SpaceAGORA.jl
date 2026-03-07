function get_navigation_callbacks(num_sats::Int, args::SimulationConfiguration)::Vector{DiscreteCallback}
    navigation_models = args.navigation_model.navigation_effectors
    navigation_rates = args.navigation_model.navigation_rates
    use_invokelatest = callback_use_invokelatest()
    callbacks = Vector{DiscreteCallback}(undef, length(navigation_models))
    for i in eachindex(navigation_models)
        navigation_model = navigation_models[i]
        navigation_rate = navigation_rates[i]
        navigation_func = (integrator) -> begin
            if use_invokelatest
                @inbounds for sat_idx in 1:num_sats
                    Base.invokelatest(calcNavigationEffect!, navigation_model, integrator.u, integrator.p, integrator.t, sat_idx)
                end
            else
                @inbounds for sat_idx in 1:num_sats
                    calcNavigationEffect!(navigation_model, integrator.u, integrator.p, integrator.t, sat_idx)
                end
            end
        end
        callbacks[i] = PeriodicCallback(navigation_func, navigation_rate)
    end
    return callbacks
end

function get_guidance_callbacks(num_sats::Int, args::SimulationConfiguration)::Vector{DiscreteCallback}
    # Implement a callback to calculate guidance commands at each time step based on the current state and the guidance model defined in the simulation configuration
    guidance_models = args.guidance_model.guidance_effectors
    guidance_rates = args.guidance_model.guidance_rates
    use_invokelatest = callback_use_invokelatest()
    callbacks = Vector{DiscreteCallback}(undef, length(guidance_models))
    for i in eachindex(guidance_models)
        guidance_model = guidance_models[i]
        guidance_rate = guidance_rates[i]
        guidance_func = (integrator) -> begin
            if use_invokelatest
                @inbounds for sat_idx in 1:num_sats
                    # Dev mode: keep Revise/hot-reload workflows free of world-age errors.
                    Base.invokelatest(calcGuidanceEffect!, guidance_model, integrator.u, integrator.p, integrator.t, sat_idx)
                end
            else
                @inbounds for sat_idx in 1:num_sats
                    # Production mode: direct dispatch avoids invokelatest overhead.
                    calcGuidanceEffect!(guidance_model, integrator.u, integrator.p, integrator.t, sat_idx)
                end
            end
        end
        callbacks[i] = PeriodicCallback(guidance_func, guidance_rate)
    end
    return callbacks
end
