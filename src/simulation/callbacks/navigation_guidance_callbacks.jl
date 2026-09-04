function get_navigation_callbacks(num_sats::Int, args::SimulationConfiguration)::Vector{DiscreteCallback}
    navigation_models = args.navigation_model.navigation_effectors
    navigation_rates = args.navigation_model.navigation_rates
    callbacks = Vector{DiscreteCallback}(undef, length(navigation_models))
    for i in eachindex(navigation_models)
        navigation_model = navigation_models[i]
        navigation_rate = navigation_rates[i]
        navigation_func = (integrator) -> begin
            @inbounds for sat_idx in 1:num_sats
                calcNavigationEffect!(navigation_model, integrator.u, integrator.p, integrator.t, sat_idx)
            end
        end
        callbacks[i] = PeriodicCallback(navigation_func, navigation_rate)
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
        guidance_func = (integrator) -> begin
            @inbounds for sat_idx in 1:num_sats
                calcGuidanceEffect!(guidance_model, integrator.u, integrator.p, integrator.t, sat_idx)
            end
        end
        callbacks[i] = PeriodicCallback(guidance_func, guidance_rate)
    end
    return callbacks
end
