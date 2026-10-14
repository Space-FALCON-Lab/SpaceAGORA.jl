
# Core-owned value projections and borrowed prediction queries for the EDG facade.
# No physical state, mutable vehicle geometry or ODEParams is passed to the provider.
using ..AbstractTypes: AbstractPlanet
using ..SpacecraftModels: SpacecraftModel
using ..CommandTypes: AerobrakingControlCommand
using ..GravityEffectors: _inverse_squared_gravity_accel, _inverse_squared_j2_gravity_accel,
    NBodyGravityModel, nbody_acceleration_ii_at_epoch

const CONTRACT_VERSION = v"1.0.0"
struct CompatibilityError <: Exception
    message::String
end
Base.showerror(io::IO, e::CompatibilityError) = print(io, e.message)
struct EDGUnavailableError <: Exception
    message::String
end
Base.showerror(io::IO, e::EDGUnavailableError) = print(io, e.message)
function require_version(expected::VersionNumber)
    expected == CONTRACT_VERSION || throw(CompatibilityError(
        "EDG requires SpaceAGORA service contract $expected; this build provides $CONTRACT_VERSION."))
    return nothing
end

const _CORE = parentmodule(parentmodule(@__MODULE__))
function require_provider()
    extension = Base.get_extension(_CORE, :SpaceAGORAEDGExt)
    extension === nothing && throw(EDGUnavailableError(
        "Energy-depletion guidance requires the optional SpaceAGORAEDG package. Install the supported package pair and load SpaceAGORAEDG before running EDG; no alternate guidance has been selected."))
    extension.require_version(CONTRACT_VERSION)
    return extension
end
function available()
    Base.get_extension(_CORE, :SpaceAGORAEDGExt) === nothing && return false
    require_provider()
    return true
end

"""A core-owned command result. Diagnostics contain copied values, not live state."""
struct EDGControlDecision{D}
    command::AerobrakingControlCommand
    apply::Bool
    switch_action::Symbol
    diagnostics::D
end

"""Borrowed prediction capability, valid only during one synchronous EDG call."""
abstract type AbstractPredictionContext end
struct _PredictionContext{F} <: AbstractPredictionContext
    query::F
    active::Base.RefValue{Bool}
end
Base.propertynames(::_PredictionContext, private::Bool=false) = ()
Base.getproperty(::_PredictionContext, ::Symbol) = throw(ArgumentError(
    "EDG prediction contexts are opaque; use the versioned query functions."))
function _context_query(context::_PredictionContext, kind, args...)
    return getfield(context, :query)(kind, args...)
end
function _borrow_context(p)
    active = Ref(true)
    query = function (kind, args...)
        active[] || throw(ArgumentError("An EDG prediction context cannot be used after its callback returns."))
        return _prediction_query(p, kind, args...)
    end
    return _PredictionContext(query, active)
end

# Named tuples and StaticArrays copy only the values consumed by the kernels.
function _planet_values(planet)
    fields = (:name, :μ, :Rp_e, :Rp_p, :T, :R, :γ, :g_ref)
    values = NamedTuple{fields}(Tuple(getproperty(planet, name) for name in fields))
    return hasproperty(planet, :H) ? merge(values, (H=Float64(planet.H),)) : values
end
_link_values(link) = (root=Bool(link.root), m=Float64(link.m), ref_area=Float64(link.ref_area),
    reflection_coefficient=Float64(link.reflection_coefficient),
    dims=SVector{3,Float64}(link.dims), α=Float64(link.α))
_geometry_values(spacecraft) = (links=Tuple(_link_values(link) for link in spacecraft.links),
    root=_link_values(spacecraft.root), prop_mass=Float64(spacecraft.prop_mass))
function _state_values(state)
    pos, vel, mass = _edg_control_pos_vel_mass(state)
    base = (;pos, vel, mass)
    return hasproperty(state, :heat_loads) ? merge(base, (heat_loads=Tuple(Float64.(state.heat_loads)),)) : base
end
_observation_values(u) = hasproperty(u, :sc) ? (sc=Tuple(_state_values(sc) for sc in u.sc),) : _state_values(u)

_prediction_query(p, ::Val{:planet}) = _planet_values(p.args.environment_model.planet)
_prediction_query(p, ::Val{:thermal}) = let thermal = p.args.environment_model.thermal_model
    hasproperty(thermal, :thermal_accomodation_factor) ?
        (thermal_accomodation_factor=Float64(thermal.thermal_accomodation_factor),) : NamedTuple()
end
_prediction_query(p, ::Val{:entry_altitude}) = 1e3 * Float64(p.args.environment_model.EI)
_prediction_query(p, ::Val{:geometry}, i) = _geometry_values(p.args.dynamics_model.spacecraft[i])
_prediction_query(p, ::Val{:environment}, u, t, i) = _edg_environment_state(u, p, t, i)
_prediction_query(p, ::Val{:trajectory_environment}, r, v, t) = _edg_targeting_prediction_environment(p, r, v, t)
_prediction_query(p, ::Val{:atmosphere}, altitude, t) = _edg_sample_prediction_atmosphere(p, altitude, t)
_prediction_query(p, ::Val{:ephemeris_time}, t) = _edg_ephemeris_time(p, t)
_prediction_query(p, ::Val{:frame}, t) = SMatrix{3,3,Float64}(_edg_planet_frame_lpi(p, t))
_prediction_query(p, ::Val{:clock}, config) = _edg_targeting_control_clock(config, p)
_prediction_query(p, ::Val{:gravity}, r, v, mass, t) = _edg_prediction_gravity_query(p, r, v, mass, t)
function _prediction_query(p, ::Val{:frame_state}, r, v, t)
    env = p.args.environment_model
    return r_intor_p!(r, v, env.planet, _edg_ephemeris_time(p, t), env.ephemerides_model)
end
function _prediction_query(p, ::Val{:scale_height})
    density_model = p.args.environment_model.density_model
    planet = p.args.environment_model.planet
    if hasproperty(density_model, :H)
        return max(Float64(density_model.H), 1.0)
    end
    return hasproperty(planet, :H) ? max(Float64(planet.H), 1.0) : 7_000.0
end

planet_constants(p::AbstractPredictionContext) = _context_query(p, Val(:planet))
thermal_parameters(p::AbstractPredictionContext) = _context_query(p, Val(:thermal))
entry_altitude_m(p::AbstractPredictionContext) = _context_query(p, Val(:entry_altitude))
spacecraft_geometry(p::AbstractPredictionContext, i::Int) = _context_query(p, Val(:geometry), i)
density_scale_height(p::AbstractPredictionContext) = _context_query(p, Val(:scale_height))
frame_state(p::AbstractPredictionContext, r, v, t) = _context_query(p, Val(:frame_state), r, v, t)
gravity_acceleration(p::AbstractPredictionContext, r, v, mass, t) = _context_query(p, Val(:gravity), r, v, mass, t)
_edg_environment_state(u, p::AbstractPredictionContext, t::Float64, i::Int) = _context_query(p, Val(:environment), u, t, i)
_edg_targeting_prediction_environment(p::AbstractPredictionContext, r, v, t::Float64) = _context_query(p, Val(:trajectory_environment), r, v, t)
_edg_sample_prediction_atmosphere(p::AbstractPredictionContext, altitude::Float64, t::Float64) = _context_query(p, Val(:atmosphere), altitude, t)
_edg_ephemeris_time(p::AbstractPredictionContext, t::Float64) = _context_query(p, Val(:ephemeris_time), t)
_edg_planet_frame_lpi(p::AbstractPredictionContext, t::Float64) = _context_query(p, Val(:frame), t)
_edg_targeting_control_clock(config, p::AbstractPredictionContext) = _context_query(p, Val(:clock), config)
_edg_in_drag_passage(p::AbstractPredictionContext, env)::Bool =
    isfinite(env.altitude_m) && isfinite(entry_altitude_m(p)) && env.altitude_m <= entry_altitude_m(p)

function _project_argument(x, contexts)
    if x isa ODEParams
        context = _borrow_context(x)
        push!(contexts, context)
        return context
    elseif x isa SpacecraftModel
        return _geometry_values(x)
    elseif x isa AbstractPlanet
        return _planet_values(x)
    elseif hasproperty(x, :sc) || (hasproperty(x, :pos) && hasproperty(x, :vel))
        return _observation_values(x)
    elseif x isa AbstractArray && !(x isa StaticArray)
        # Legacy helper arguments are copied before crossing the package boundary.
        return copy(x)
    end
    return x
end
function _invoke_kernel(name::Symbol, args...; kwargs...)
    extension = require_provider()
    contexts = AbstractPredictionContext[]
    try
        projected = map(x -> _project_argument(x, contexts), args)
        return extension.call_kernel(name, projected...; kwargs...)
    finally
        for context in contexts
            getfield(context, :active)[] = false
        end
    end
end

function _edg_prediction_gravity_query(
    p::ODEParams,
    r::SVector{3, Float64},
    v::SVector{3, Float64},
    mass::Float64,
    t_abs::Float64,
)
    planet = p.args.environment_model.planet
    ephemerides_model = p.args.environment_model.ephemerides_model
    et = _edg_ephemeris_time(p, t_abs)
    pos_pp, _ = r_intor_p!(r, v, planet, et, ephemerides_model)
    dynamic_effectors = p.args.dynamics_model.dynamic_effectors
    use_j2 = any(nameof(typeof(effector)) == :InverseSquaredJ2GravityModel for effector in dynamic_effectors)
    acceleration_ii = if use_j2
        l_pi = _edg_planet_frame_lpi(p, t_abs)
        SVector{3, Float64}(l_pi' * _inverse_squared_j2_gravity_accel(pos_pp, planet))
    else
        _inverse_squared_gravity_accel(r, planet)
    end
    et = _edg_ephemeris_time(p, t_abs)
    for effector in dynamic_effectors
        if effector isa NBodyGravityModel
            acceleration_ii += nbody_acceleration_ii_at_epoch(effector, r, p, et)
        end
    end
    return acceleration_ii
end

