module LaserLinkEffectors

using LinearAlgebra
using StaticArrays
using DiffEqBase

using ...AbstractTypes: AbstractForceTorqueModel
import ..DynamicEffectors: calcForceTorque, solver_partition

export OpenCavityLaserLinkModel, laser_link_scheduler_callback
export laser_link_force_magnitude, laser_link_pair_force, laser_link_active_pair
export update_laser_link_schedule!, accumulate_laser_link_forces!

const SPEED_OF_LIGHT_MPS = 299_792_458.0
const SUPPORTED_LASER_SCHEDULES = (:naive_next_entering, :positive_along_track)

"""
    OpenCavityLaserLinkModel(; kwargs...)

Open-cavity laser interlink force model for one target satellite and a helper
population. The model applies at most one active target-helper link at a time,
with equal-and-opposite inertial-frame forces along the instantaneous
inter-satellite line of sight.
"""
mutable struct OpenCavityLaserLinkModel <: AbstractForceTorqueModel
    target_idx::Int
    helper_indices::Vector{Int}
    range_m::Float64
    power_w::Float64
    magnification::Float64
    beta::Float64
    eta::Float64
    schedule::Symbol
    active_helper_idx::Int
    previous_in_range::Vector{Bool}
    link_activation_count::Int
    active_link_step_count::Int
end

function OpenCavityLaserLinkModel(;
    target_idx::Integer=1,
    helper_indices::AbstractVector{<:Integer}=Int[],
    range_m::Real=200e3,
    power_w::Real=10_000.0,
    magnification::Real=100.0,
    beta::Real=1.0,
    eta::Real=2.0,
    schedule::Symbol=:naive_next_entering,
    active_helper_idx::Integer=0,
    previous_in_range::AbstractVector{Bool}=Bool[],
    link_activation_count::Integer=0,
    active_link_step_count::Integer=0,
)
    model = OpenCavityLaserLinkModel(
        Int(target_idx),
        collect(Int, helper_indices),
        Float64(range_m),
        Float64(power_w),
        Float64(magnification),
        Float64(beta),
        Float64(eta),
        schedule,
        Int(active_helper_idx),
        collect(Bool, previous_in_range),
        Int(link_activation_count),
        Int(active_link_step_count),
    )
    _validate_laser_link_model!(model)
    _ensure_laser_link_state!(model)
    return model
end

function OpenCavityLaserLinkModel(
    target_idx::Integer,
    helper_indices::AbstractVector{<:Integer};
    range_m::Real=200e3,
    power_w::Real=10_000.0,
    magnification::Real=100.0,
    beta::Real=1.0,
    eta::Real=2.0,
    schedule::Symbol=:naive_next_entering,
)
    model = OpenCavityLaserLinkModel(
        target_idx=Int(target_idx),
        helper_indices=collect(Int, helper_indices),
        range_m=Float64(range_m),
        power_w=Float64(power_w),
        magnification=Float64(magnification),
        beta=Float64(beta),
        eta=Float64(eta),
        schedule=schedule,
    )
    _validate_laser_link_model!(model)
    _ensure_laser_link_state!(model)
    return model
end

function _validate_laser_link_model!(model::OpenCavityLaserLinkModel)::Nothing
    model.target_idx > 0 || throw(ArgumentError("OpenCavityLaserLinkModel.target_idx must be positive."))
    isempty(model.helper_indices) && throw(ArgumentError("OpenCavityLaserLinkModel requires at least one helper index."))
    any(==(model.target_idx), model.helper_indices) &&
        throw(ArgumentError("OpenCavityLaserLinkModel helper_indices cannot include target_idx."))
    length(unique(model.helper_indices)) == length(model.helper_indices) ||
        throw(ArgumentError("OpenCavityLaserLinkModel helper_indices must be unique."))
    model.range_m >= 0.0 || throw(ArgumentError("OpenCavityLaserLinkModel.range_m must be >= 0.0."))
    model.power_w >= 0.0 || throw(ArgumentError("OpenCavityLaserLinkModel.power_w must be >= 0.0."))
    model.magnification >= 0.0 || throw(ArgumentError("OpenCavityLaserLinkModel.magnification must be >= 0.0."))
    model.beta >= 0.0 || throw(ArgumentError("OpenCavityLaserLinkModel.beta must be >= 0.0."))
    model.eta >= 0.0 || throw(ArgumentError("OpenCavityLaserLinkModel.eta must be >= 0.0."))
    model.schedule in SUPPORTED_LASER_SCHEDULES ||
        throw(ArgumentError("Unsupported OpenCavityLaserLinkModel.schedule=$(repr(model.schedule)). Supported schedules: $(SUPPORTED_LASER_SCHEDULES)."))
    return nothing
end

@inline function _ensure_laser_link_state!(model::OpenCavityLaserLinkModel)::Nothing
    if length(model.previous_in_range) != length(model.helper_indices)
        resize!(model.previous_in_range, length(model.helper_indices))
        fill!(model.previous_in_range, false)
    end
    return nothing
end

@inline solver_partition(::OpenCavityLaserLinkModel) = :explicit

function calcForceTorque(
    ::OpenCavityLaserLinkModel,
    x,
    p,
    i::Int64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end

@inline laser_link_force_magnitude(model::OpenCavityLaserLinkModel)::Float64 =
    model.eta * model.beta * model.magnification * model.power_w / SPEED_OF_LIGHT_MPS

@inline function _extract_pos_vel_mass(x)
    pos = SVector{3, Float64}(x[1], x[2], x[3])
    vel = SVector{3, Float64}(x[4], x[5], x[6])
    mass = length(x) >= 7 ? Float64(x[7]) : NaN
    return pos, vel, mass
end

function _state_vectors(u)::Tuple{Vector{SVector{3, Float64}}, Vector{SVector{3, Float64}}}
    sc_state = u.sc
    pos = Vector{SVector{3, Float64}}(undef, length(sc_state))
    vel = Vector{SVector{3, Float64}}(undef, length(sc_state))
    @inbounds for idx in eachindex(sc_state)
        pos[idx], vel[idx], _ = _extract_pos_vel_mass(sc_state[idx])
    end
    return pos, vel
end

@inline _helper_slot(model::OpenCavityLaserLinkModel, helper_idx::Int)::Union{Nothing, Int} =
    findfirst(==(helper_idx), model.helper_indices)

function _in_range_flags!(
    flags::Vector{Bool},
    model::OpenCavityLaserLinkModel,
    pos::AbstractVector{SVector{3, Float64}},
)::Vector{Bool}
    _ensure_laser_link_state!(model)
    if length(flags) != length(model.helper_indices)
        resize!(flags, length(model.helper_indices))
    end
    target_pos = pos[model.target_idx]
    @inbounds for slot in eachindex(model.helper_indices)
        helper_idx = model.helper_indices[slot]
        flags[slot] = norm(pos[helper_idx] - target_pos) <= model.range_m
    end
    return flags
end

@inline function _rtn_basis(pos_t::SVector{3, Float64}, vel_t::SVector{3, Float64})
    rnorm = norm(pos_t)
    rhat = rnorm > 0.0 ? pos_t / rnorm : SVector{3, Float64}(1.0, 0.0, 0.0)
    h = cross(pos_t, vel_t)
    hnorm = norm(h)
    nhat = hnorm > 0.0 ? h / hnorm : SVector{3, Float64}(0.0, 0.0, 1.0)
    that = cross(nhat, rhat)
    return rhat, that, nhat
end

@inline function _along_track_projection(
    model::OpenCavityLaserLinkModel,
    helper_idx::Int,
    pos::AbstractVector{SVector{3, Float64}},
    vel::AbstractVector{SVector{3, Float64}},
)::Float64
    target_pos = pos[model.target_idx]
    rel = pos[helper_idx] - target_pos
    rho = norm(rel)
    rho > 0.0 || return 0.0
    _, that, _ = _rtn_basis(target_pos, vel[model.target_idx])
    return dot(rel / rho, that)
end

function _closest_in_range_helper(
    model::OpenCavityLaserLinkModel,
    pos::AbstractVector{SVector{3, Float64}},
    in_range::AbstractVector{Bool},
)::Int
    target_pos = pos[model.target_idx]
    best_idx = 0
    best_range = Inf
    @inbounds for slot in eachindex(model.helper_indices)
        in_range[slot] || continue
        helper_idx = model.helper_indices[slot]
        rho = norm(pos[helper_idx] - target_pos)
        if rho < best_range
            best_idx = helper_idx
            best_range = rho
        end
    end
    return best_idx
end

function _closest_entering_helper(
    model::OpenCavityLaserLinkModel,
    pos::AbstractVector{SVector{3, Float64}},
    in_range::AbstractVector{Bool},
)::Int
    target_pos = pos[model.target_idx]
    best_idx = 0
    best_range = Inf
    @inbounds for slot in eachindex(model.helper_indices)
        in_range[slot] || continue
        model.previous_in_range[slot] && continue
        helper_idx = model.helper_indices[slot]
        rho = norm(pos[helper_idx] - target_pos)
        if rho < best_range
            best_idx = helper_idx
            best_range = rho
        end
    end
    return best_idx
end

function _activate_helper!(model::OpenCavityLaserLinkModel, helper_idx::Int)::Nothing
    if helper_idx > 0 && model.active_helper_idx != helper_idx
        model.link_activation_count += 1
    end
    model.active_helper_idx = helper_idx
    return nothing
end

function update_laser_link_schedule!(
    model::OpenCavityLaserLinkModel,
    pos::AbstractVector{SVector{3, Float64}},
    vel::AbstractVector{SVector{3, Float64}},
)::Int
    _validate_laser_link_model!(model)
    in_range = _in_range_flags!(Bool[], model, pos)

    if model.schedule === :naive_next_entering
        if model.active_helper_idx > 0
            slot = _helper_slot(model, model.active_helper_idx)
            if slot === nothing || !in_range[slot]
                model.active_helper_idx = 0
                helper_idx = _closest_entering_helper(model, pos, in_range)
                helper_idx > 0 && _activate_helper!(model, helper_idx)
            end
        else
            helper_idx = any(model.previous_in_range) ?
                _closest_entering_helper(model, pos, in_range) :
                _closest_in_range_helper(model, pos, in_range)
            helper_idx > 0 && _activate_helper!(model, helper_idx)
        end
    elseif model.schedule === :positive_along_track
        if model.active_helper_idx > 0
            slot = _helper_slot(model, model.active_helper_idx)
            valid = slot !== nothing &&
                in_range[slot] &&
                _along_track_projection(model, model.active_helper_idx, pos, vel) > 0.0
            valid || (model.active_helper_idx = 0)
        end
        if model.active_helper_idx == 0
            best_idx = 0
            best_projection = 0.0
            @inbounds for slot in eachindex(model.helper_indices)
                in_range[slot] || continue
                helper_idx = model.helper_indices[slot]
                projection = _along_track_projection(model, helper_idx, pos, vel)
                if projection > best_projection
                    best_projection = projection
                    best_idx = helper_idx
                end
            end
            best_idx > 0 && _activate_helper!(model, best_idx)
        end
    end

    model.previous_in_range .= in_range
    model.active_helper_idx > 0 && (model.active_link_step_count += 1)
    return model.active_helper_idx
end

function update_laser_link_schedule!(model::OpenCavityLaserLinkModel, u)::Int
    pos, vel = _state_vectors(u)
    return update_laser_link_schedule!(model, pos, vel)
end

@inline laser_link_active_pair(model::OpenCavityLaserLinkModel)::Tuple{Int, Int} =
    (model.target_idx, model.active_helper_idx)

function laser_link_pair_force(
    model::OpenCavityLaserLinkModel,
    target_pos::SVector{3, Float64},
    helper_pos::SVector{3, Float64},
)::SVector{3, Float64}
    rel = helper_pos - target_pos
    rho = norm(rel)
    (rho > 0.0 && rho <= model.range_m) || return SVector{3, Float64}(0.0, 0.0, 0.0)
    return laser_link_force_magnitude(model) * rel / rho
end

function accumulate_laser_link_forces!(
    totals::AbstractMatrix{Float64},
    model::OpenCavityLaserLinkModel,
    pos::AbstractVector{SVector{3, Float64}},
    active_flags,
)::Nothing
    helper_idx = model.active_helper_idx
    helper_idx > 0 || return nothing
    target_idx = model.target_idx
    target_idx <= length(pos) && helper_idx <= length(pos) || return nothing
    active_flags[target_idx] && active_flags[helper_idx] || return nothing

    force = laser_link_pair_force(model, pos[target_idx], pos[helper_idx])
    @inbounds begin
        totals[1, target_idx] += force[1]
        totals[2, target_idx] += force[2]
        totals[3, target_idx] += force[3]
        totals[1, helper_idx] -= force[1]
        totals[2, helper_idx] -= force[2]
        totals[3, helper_idx] -= force[3]
    end
    return nothing
end

function _update_matching_laser_models!(template::OpenCavityLaserLinkModel, integrator)::Nothing
    for effector in integrator.p.args.dynamics_model.dynamic_effectors
        effector isa OpenCavityLaserLinkModel || continue
        if effector === template ||
           (effector.target_idx == template.target_idx && effector.helper_indices == template.helper_indices)
            update_laser_link_schedule!(effector, integrator.u)
        end
    end
    return nothing
end

function laser_link_scheduler_callback(model::OpenCavityLaserLinkModel)
    condition(u, t, integrator) = true
    affect!(integrator) = _update_matching_laser_models!(model, integrator)
    initialize = (cb, u, t, integrator) -> _update_matching_laser_models!(model, integrator)
    return DiffEqBase.DiscreteCallback(condition, affect!; initialize=initialize)
end

end
