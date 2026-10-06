# Construct the nominal KS trajectory from atmospheric entry through exit.
Base.@kwdef struct AerobrakingMPCReferenceConfig
    h_cut_m::Float64
    delta_s::Float64 = 10.0e-7
    max_coast_steps::Int
    max_pass_steps::Int
end

function _altitude_from_ks_state(state, params::AerobrakingMPCParams)
    cart = ks_state_to_cartesian(state)
    return norm(cart.position_ii_m) - params.Re
end

function _radial_velocity_from_ks_state(state)
    cart = ks_state_to_cartesian(state)
    return dot(cart.position_ii_m, cart.velocity_ii_m) / norm(cart.position_ii_m)
end

function build_reference_drag_pass(
    params::AerobrakingMPCParams,
    position_ii_m,
    velocity_ii_m;
    config::AerobrakingMPCConfig,
    reference::AerobrakingMPCReferenceConfig,
    nominal_area_m2::Real=config.bus_reference_area_m2 + config.controllable_area_m2,
    density::Function=(altitude_m, elapsed_time_s) -> 0.0,
)
    area = Float64(nominal_area_m2)
    state = cartesian_to_ks_state(position_ii_m, velocity_ii_m, params)
    h_prev = _altitude_from_ks_state(state, params)
    reached = h_prev <= reference.h_cut_m && _radial_velocity_from_ks_state(state) < 0.0

    for _ in 1:(reached ? 0 : reference.max_coast_steps)
        next_state = ks_rk4_step(
            state,
            params,
            area,
            reference.delta_s,
        )
        h_next = _altitude_from_ks_state(next_state, params)
        vr_next = _radial_velocity_from_ks_state(next_state)
        if h_prev > reference.h_cut_m && h_next <= reference.h_cut_m && vr_next < 0.0
            α = (h_prev - reference.h_cut_m) / (h_prev - h_next)
            state = state .+ α .* (next_state .- state)
            reached = true
            break
        end
        state = next_state
        h_prev = h_next
    end
    reached || throw(ErrorException("Reference propagation did not reach inbound cutoff altitude."))

    rows = Vector{Vector{Float64}}()
    push!(rows, copy(state))
    exited = false
    for _ in 1:reference.max_pass_steps
        h_now = _altitude_from_ks_state(state, params)
        next_state = ks_rk4_step(
            state,
            params,
            area,
            reference.delta_s;
            density_kg_m3=density,
            drag_coefficient=config.drag_coefficient,
            mass_kg=config.mass_kg,
            use_drag=true,
        )
        h_next = _altitude_from_ks_state(next_state, params)
        vr_next = _radial_velocity_from_ks_state(next_state)
        if h_now <= reference.h_cut_m && h_next > reference.h_cut_m && vr_next > 0.0
            α = (reference.h_cut_m - h_now) / (h_next - h_now)
            push!(rows, copy(state .+ α .* (next_state .- state)))
            exited = true
            break
        end
        all(isfinite, next_state) || throw(ErrorException("Reference propagation became non-finite during the atmospheric pass."))
        state = next_state
        push!(rows, copy(state))
    end
    exited || throw(ErrorException("Reference propagation did not reach the outbound cutoff altitude."))

    states = reduce(vcat, transpose.(rows))
    times_s = states[:, 10]
    altitude_m = [_altitude_from_ks_state(view(states, i, :), params) for i in axes(states, 1)]
    return (
        states=states,
        time_s=times_s,
        altitude_m=altitude_m,
        nominal_area_m2=area,
        cutoff_altitude_m=reference.h_cut_m,
        delta_s=reference.delta_s,
    )
end
