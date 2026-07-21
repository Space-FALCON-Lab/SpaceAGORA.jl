using ComponentArrays
using RecursiveArrayTools: ArrayPartition
using StaticArrays

@inline function _is_flat_xyz_state(u)::Bool
    return !hasproperty(u, :sc) &&
        u isa AbstractVector &&
        eltype(u) <: Real &&
        length(u) >= 3
end

@inline function _is_flat_rv_state(u)::Bool
    return _is_flat_xyz_state(u) && length(u) >= 6
end

@inline function _check_flat_single_sat(u, sat_idx::Int)
    sat_idx == 1 || throw(BoundsError(u, sat_idx))
    return nothing
end

@inline function _is_gravity_backbone_state(u)::Bool
    return u isa ArrayPartition &&
        length(u.x) == 2 &&
        ((hasproperty(u.x[1], :sc) && hasproperty(u.x[2], :sc)) ||
         (_is_flat_xyz_state(u.x[1]) && _is_flat_xyz_state(u.x[2])))
end

@inline function _gravity_backbone_velocity_state(u)
    _is_gravity_backbone_state(u) || throw(ArgumentError("Expected gravity-backbone second-order state."))
    return u.x[1]
end

@inline function _gravity_backbone_position_state(u)
    _is_gravity_backbone_state(u) || throw(ArgumentError("Expected gravity-backbone second-order state."))
    return u.x[2]
end

@inline function _gravity_backbone_spacecraft_state(u)
    if _is_gravity_backbone_state(u)
        return _gravity_backbone_position_state(u)
    end
    if hasproperty(u, :sc)
        return getproperty(u, :sc)
    end
    return u
end

@inline function _state_position_ii(u, sat_idx::Int)::SVector{3, Float64}
    if _is_gravity_backbone_state(u)
        spacecraft_state = _gravity_backbone_spacecraft_state(u)
        if hasproperty(spacecraft_state, :sc)
            sc_state = spacecraft_state.sc[sat_idx]
            return SVector{3, Float64}(sc_state[1], sc_state[2], sc_state[3])
        end
        if _is_flat_xyz_state(spacecraft_state)
            _check_flat_single_sat(spacecraft_state, sat_idx)
            return SVector{3, Float64}(spacecraft_state[1], spacecraft_state[2], spacecraft_state[3])
        end
        sc_state = spacecraft_state[sat_idx]
        return SVector{3, Float64}(sc_state[1], sc_state[2], sc_state[3])
    end
    if _is_flat_rv_state(u)
        _check_flat_single_sat(u, sat_idx)
        return SVector{3, Float64}(u[1], u[2], u[3])
    end
    sc_state = u.sc[sat_idx]
    return SVector{3, Float64}(sc_state[1], sc_state[2], sc_state[3])
end

@inline function _state_velocity_ii(u, sat_idx::Int)::SVector{3, Float64}
    if _is_gravity_backbone_state(u)
        spacecraft_state = _gravity_backbone_velocity_state(u)
        if hasproperty(spacecraft_state, :sc)
            sc_state = spacecraft_state.sc[sat_idx]
            return SVector{3, Float64}(sc_state[1], sc_state[2], sc_state[3])
        end
        if _is_flat_xyz_state(spacecraft_state)
            _check_flat_single_sat(spacecraft_state, sat_idx)
            return SVector{3, Float64}(spacecraft_state[1], spacecraft_state[2], spacecraft_state[3])
        end
        sc_state = spacecraft_state[sat_idx]
        return SVector{3, Float64}(sc_state[1], sc_state[2], sc_state[3])
    end
    if _is_flat_rv_state(u)
        _check_flat_single_sat(u, sat_idx)
        return SVector{3, Float64}(u[4], u[5], u[6])
    end
    sc_state = u.sc[sat_idx]
    return SVector{3, Float64}(sc_state[4], sc_state[5], sc_state[6])
end

@inline function _state_mass_kg(u, args, sat_idx::Int)::Float64
    if _is_gravity_backbone_state(u) || _is_flat_rv_state(u)
        spacecraft = args.dynamics_model.spacecraft[sat_idx]
        return spacecraft.dry_mass + spacecraft.prop_mass
    end
    return Float64(u.sc[sat_idx].mass)
end

@inline function _state_has_mass(u, args, sat_idx::Int)::Bool
    return true
end

@inline function _state_heat_loads(u, args, sat_idx::Int)
    if _is_gravity_backbone_state(u)
        return zeros(Float64, length(args.dynamics_model.spacecraft[sat_idx].links))
    end
    if _is_flat_rv_state(u)
        return Float64[]
    end
    return u.sc[sat_idx].heat_loads
end

@inline function _state_has_heat_loads(u, args, sat_idx::Int)::Bool
    return !(_is_gravity_backbone_state(u) || _is_flat_rv_state(u))
end

@inline function _state_quaternion(u, sat_idx::Int)
    if _is_gravity_backbone_state(u) || _is_flat_rv_state(u) || !hasproperty(u.sc[sat_idx], :q)
        return nothing
    end
    return SVector{4, Float64}(u.sc[sat_idx].q)
end

@inline function _state_has_quaternion(u, sat_idx::Int)::Bool
    return !isnothing(_state_quaternion(u, sat_idx))
end

@inline function _gravity_backbone_initial_states(u0, args)
    if _is_gravity_backbone_state(u0)
        return deepcopy(_gravity_backbone_position_state(u0)), deepcopy(_gravity_backbone_velocity_state(u0))
    end

    q_shapes = map(args.dynamics_model.spacecraft) do _
        (pos=zeros(3),)
    end
    dq_shapes = map(args.dynamics_model.spacecraft) do _
        (vel=zeros(3),)
    end
    q0 = ComponentVector(sc=q_shapes)
    dq0 = ComponentVector(sc=dq_shapes)
    @inbounds for i in eachindex(args.dynamics_model.spacecraft)
        q0.sc[i].pos .= u0.sc[i].pos
        dq0.sc[i].vel .= u0.sc[i].vel
    end
    return q0, dq0
end
