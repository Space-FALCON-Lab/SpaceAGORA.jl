# Extracted from Evan's PR121 (cc891c76), with construction validation,
# pure wrench evaluation, and accepted-step diagnostic integration.
# This is a phenomenological plume/regolith closure, not a validated Apollo
# calibration or a resolved gas/soil simulation.
module PlumeSurfaceInteraction
using StaticArrays, LinearAlgebra
using ...AbstractTypes: AbstractForceTorqueModel, AbstractTerrainModel
using ...TerrainModels: NoTerrainModel, DEMTerrainModel, terrain_radius
using ...EffectorSampling: StateSample, EnvironmentSample, EffectorEnvironmentRequirements
using ...SimulationModel: rot
import ..DynamicEffectors: wrench, wrench_caching!, environment_requirements
export PlumeSurfaceConfig, PlumeSurfaceInteractionModel, PlumeSurfaceState
export plume_surface_footprint, plume_erosion_onset_height, plume_ground_effect_force, plume_quantities

"""
    PlumeSurfaceConfig(; ...)

Phenomenological Gaussian pressure footprint, excess-shear erosion and near-ground
thrust augmentation. SI units are used except `plume_half_angle_deg`. All values
must be finite; geometric/material scales are positive and coefficients are
nonnegative. The defaults are illustrative model choices, not a validated flight
calibration. `expansion_ratio`, `chamber_pressure_pa`, `exit_mach`,
`bulk_density_kg_m3` and `cohesion_pa` are recorded metadata, not closure inputs.

`nozzle_offset_m` locates the exit plane along body +z (the exhaust direction)
from the integrated vehicle reference point. The augmentation acts along the
opposite direction on this same line, so its torque about that reference point
is zero. No off-axis nozzle or terrain-normal torque is modeled.

The footprint radius is the larger of the nozzle radius and exit-plane height
multiplied by `tan(plume_half_angle_deg)`. Wall shear is proportional to the
Gaussian pressure, and excess shear integrated to three footprint radii sets
erosion through `erosion_efficiency` and the bounded grain speed. The spatial
quadrature uses 128 midpoint cells. `max_height_m` is the exit-plane height cutoff.
"""
Base.@kwdef struct PlumeSurfaceConfig
    nozzle_exit_radius_m::Float64 = 0.75
    nozzle_offset_m::Float64 = 1.5
    expansion_ratio::Float64 = 47.5
    chamber_pressure_pa::Float64 = 7.2e5
    exit_mach::Float64 = 4.8
    plume_half_angle_deg::Float64 = 25.0
    friction_coefficient::Float64 = 0.01
    bulk_density_kg_m3::Float64 = 1_500.0
    particle_density_kg_m3::Float64 = 3_100.0
    particle_diameter_m::Float64 = 70.0e-6
    cohesion_pa::Float64 = 1.0e-3
    threshold_shear_pa::Float64 = 0.15
    erosion_efficiency::Float64 = 10.0
    particle_drag_coefficient::Float64 = 2.0
    ejecta_speed_min_mps::Float64 = 5.0
    ejecta_speed_max_mps::Float64 = 150.0
    ground_effect_max_fraction::Float64 = 0.03
    ground_effect_scale::Float64 = 0.8
    ground_effect_cutoff::Float64 = 2.0
    max_height_m::Float64 = 250.0
    function PlumeSurfaceConfig(nozzle_exit_radius_m,
        nozzle_offset_m,
        expansion_ratio,
        chamber_pressure_pa,
        exit_mach,
        plume_half_angle_deg,
        friction_coefficient,
        bulk_density_kg_m3,
        particle_density_kg_m3,
        particle_diameter_m,
        cohesion_pa,
        threshold_shear_pa,
        erosion_efficiency,
        particle_drag_coefficient,
        ejecta_speed_min_mps,
        ejecta_speed_max_mps,
        ground_effect_max_fraction,
        ground_effect_scale,
        ground_effect_cutoff,
        max_height_m)
        values = Float64.((nozzle_exit_radius_m, nozzle_offset_m, expansion_ratio, chamber_pressure_pa, exit_mach, plume_half_angle_deg, friction_coefficient, bulk_density_kg_m3, particle_density_kg_m3, particle_diameter_m, cohesion_pa, threshold_shear_pa, erosion_efficiency, particle_drag_coefficient, ejecta_speed_min_mps, ejecta_speed_max_mps, ground_effect_max_fraction, ground_effect_scale, ground_effect_cutoff, max_height_m,))
        all(isfinite, values) || throw(ArgumentError("Plume config values must be finite."))
        config = NamedTuple{(:nozzle_exit_radius_m, :nozzle_offset_m, :expansion_ratio, :chamber_pressure_pa, :exit_mach, :plume_half_angle_deg, :friction_coefficient, :bulk_density_kg_m3, :particle_density_kg_m3, :particle_diameter_m, :cohesion_pa, :threshold_shear_pa, :erosion_efficiency, :particle_drag_coefficient, :ejecta_speed_min_mps, :ejecta_speed_max_mps, :ground_effect_max_fraction, :ground_effect_scale, :ground_effect_cutoff, :max_height_m,)}(values)
        _validate_plume_config(config)
        return new(values...)
    end
end

function _validate_plume_config(c)
    for name in (:nozzle_exit_radius_m, :chamber_pressure_pa, :exit_mach,
                 :bulk_density_kg_m3, :particle_density_kg_m3, :particle_diameter_m,
                 :ejecta_speed_min_mps, :ground_effect_scale, :ground_effect_cutoff, :max_height_m)
        getproperty(c, name) > 0 || throw(ArgumentError("Plume $name must be positive."))
    end
    for name in (:nozzle_offset_m, :friction_coefficient, :cohesion_pa,
                 :threshold_shear_pa, :erosion_efficiency, :particle_drag_coefficient)
        getproperty(c, name) >= 0 || throw(ArgumentError("Plume $name must be nonnegative."))
    end
    c.expansion_ratio >= 1 || throw(ArgumentError("Plume expansion_ratio must be >= 1."))
    0 < c.plume_half_angle_deg < 90 || throw(ArgumentError("Plume half-angle must lie strictly between 0 and 90 degrees."))
    c.ejecta_speed_max_mps >= c.ejecta_speed_min_mps || throw(ArgumentError("Plume ejecta speed bounds are reversed."))
    0 <= c.ground_effect_max_fraction <= 1 || throw(ArgumentError("Plume ground-effect fraction must lie in [0, 1]."))
    return nothing
end

"""
    PlumeSurfaceState(num_sats)

Per-spacecraft accepted-state diagnostics and diagnostic erosion integral.
`eroded_kg` uses trapezoidal quadrature over accepted solver steps; it is not an
ODE state controlled by solver tolerances. Check convergence by reducing the
maximum step. Initialization resets each new run and samples its initial state;
checkpoint continuation within that run preserves the accumulated diagnostic.
A separately resumed run begins a new diagnostic integral at its resume time.
"""
mutable struct PlumeSurfaceState
    height_m::Vector{Float64}
    pressure_pa::Vector{Float64}
    shear_pa::Vector{Float64}
    erosion_kg_s::Vector{Float64}
    eroded_kg::Vector{Float64}
    ejecta_mps::Vector{Float64}
    ground_effect_n::Vector{Float64}
    last_time_s::Vector{Float64}
    last_rate_kg_s::Vector{Float64}
    last_thrust_n::Vector{Float64}
    previous_time_s::Vector{Float64}
    previous_rate_kg_s::Vector{Float64}
    previous_eroded_kg::Vector{Float64}
    previous_thrust_n::Vector{Float64}
    interval_end_rate_kg_s::Vector{Float64}
end
function PlumeSurfaceState(num_sats::Integer)
    num_sats >= 1 || throw(ArgumentError("PlumeSurfaceState needs at least one spacecraft."))
    n = Int(num_sats)
    return PlumeSurfaceState(zeros(n), zeros(n), zeros(n), zeros(n), zeros(n), zeros(n), zeros(n), fill(NaN, n), zeros(n), zeros(n), fill(NaN, n), zeros(n), zeros(n), zeros(n), zeros(n))
end

"""
    PlumeSurfaceInteractionModel(control, terrain=NoTerrainModel(); config=PlumeSurfaceConfig(), num_sats=..., reference_radius_m=nothing)

Add only modeled ground-effect augmentation to an existing engine's thrust.
`control.actuators.thrust_n` must be a finite, nonnegative vector in simulation
spacecraft order and remain held between accepted control updates. This effector
does not apply primary engine thrust. Its length must match the simulation.

Geometry uses a local radial tangent plane at the subspacecraft terrain sample,
with planetocentric latitude/longitude in degrees and the terrain reference
sphere. Supply `reference_radius_m` when touchdown/guidance uses an explicit
sphere; it must agree with a DEM's declared radius. When omitted, DEM terrain
uses its declared radius and other terrain uses the planet's equatorial radius.
Slant height is radial clearance divided by the downward axis cosine;
curvature, the tilted impact point, DEM slopes and ray tracing are not modeled.
Horizontal or upward exhaust gives no interaction. Without attitude the exhaust
is radial inward; with attitude it is body +z. The scalar nozzle offset lies on
that same axis, so augmentation returns no reference-point torque.

`wrench` and `wrench_caching!` do not mutate plume diagnostics. The engine's
accepted-step callback owns the diagnostic state and mass quadrature. Default
save fields sample the supplied saved state and the held interval thrust, with
mass interpolated from accepted-step quadrature. The four-argument pure wrench
is valid for one spacecraft; pass an explicit fifth spacecraft index otherwise.
"""
struct PlumeSurfaceInteractionModel{C, T <: AbstractTerrainModel} <: AbstractForceTorqueModel
    config::PlumeSurfaceConfig
    control::C
    terrain::T
    state::PlumeSurfaceState
    reference_radius_m::Union{Nothing,Float64}
    function PlumeSurfaceInteractionModel(config::PlumeSurfaceConfig, control,
            terrain::T, state::PlumeSurfaceState, reference_radius_m=nothing) where {T <: AbstractTerrainModel}
        n = _control_spacecraft_count(control)
        _validate_plume_state(state, n)
        radius = if reference_radius_m === nothing
            nothing
        else
            reference_radius_m isa Real && isfinite(reference_radius_m) && reference_radius_m > 0 ||
                throw(ArgumentError("Plume reference_radius_m must be finite and positive."))
            value = Float64(reference_radius_m)
            isfinite(value) && value > 0 ||
                throw(ArgumentError("Plume reference_radius_m must be representable as a finite positive Float64."))
            terrain isa DEMTerrainModel && value != terrain.reference_radius_m &&
                throw(ArgumentError("Plume reference_radius_m differs from the DEM reference radius."))
            value
        end
        return new{typeof(control), T}(config, control, terrain, state, radius)
    end
end
function PlumeSurfaceInteractionModel(control, terrain::AbstractTerrainModel=NoTerrainModel();
        config::PlumeSurfaceConfig=PlumeSurfaceConfig(), num_sats::Integer=_control_spacecraft_count(control),
        reference_radius_m=nothing)
    return PlumeSurfaceInteractionModel(config, control, terrain, PlumeSurfaceState(num_sats), reference_radius_m)
end
function _control_spacecraft_count(control)
    hasproperty(control, :actuators) && hasproperty(control.actuators, :thrust_n) ||
        throw(ArgumentError("Plume control must provide actuators.thrust_n."))
    values = control.actuators.thrust_n
    values isa AbstractVector || throw(ArgumentError("Plume thrust_n must be a vector."))
    !isempty(values) && all(v -> v isa Real && isfinite(v) && v >= 0, values) ||
        throw(ArgumentError("Plume thrust values must be finite, nonnegative and nonempty."))
    return length(values)
end
function _validate_plume_state(state, n)
    n >= 1 || throw(ArgumentError("Plume needs at least one spacecraft."))
    all(name -> length(getfield(state, name)) == n, fieldnames(PlumeSurfaceState)) ||
        throw(ArgumentError("All plume state vectors and thrust_n must have matching spacecraft counts."))
    return nothing
end
function _engine_thrust_n(model, i::Int)
    1 <= i <= length(model.state.height_m) || throw(BoundsError(model.state.height_m, i))
    thrust = Float64(model.control.actuators.thrust_n[i])
    isfinite(thrust) && thrust >= 0 || throw(ArgumentError("Plume thrust must be finite and nonnegative."))
    return thrust
end
function _plume_inputs(thrust_n, height_m)
    F, h = Float64(thrust_n), Float64(height_m)
    isfinite(F) && F >= 0 || throw(ArgumentError("Plume thrust must be finite and nonnegative."))
    isfinite(h) && h >= 0 || throw(ArgumentError("Plume height must be finite and nonnegative."))
    return F, h
end

const _PLUME_QUADRATURE_POINTS = 128
const _PLUME_QUADRATURE_LIMIT = 3.0        # footprint radii; exp(-9) is already negligible

"Peak wall shear stress of the Gaussian footprint, at r = R_p / sqrt(2)."
@inline _plume_peak_shear(cfg::PlumeSurfaceConfig, p0::Float64)::Float64 = cfg.friction_coefficient * p0 * sqrt(2.0) * exp(-0.5)

"""
Peak Gaussian pressure (Pa) and footprint radius (m), normalized over the infinite plane to thrust (N). Inputs must be finite and nonnegative.
"""
function plume_surface_footprint(cfg::PlumeSurfaceConfig, thrust_n::Real, height_m::Real)
    F, h = _plume_inputs(thrust_n, height_m)
    R = max(h * tand(cfg.plume_half_angle_deg), cfg.nozzle_exit_radius_m)
    p0 = F / (pi * R * R)
    return p0, R
end

"""
Largest modeled exit-plane height (m) at which excess shear can occur, bounded by max_height_m and the nozzle-radius floor; zero when erosion is absent.
"""
function plume_erosion_onset_height(cfg::PlumeSurfaceConfig, thrust_n::Real)::Float64
    F, _ = _plume_inputs(thrust_n, 0.0)
    F > 0 && cfg.friction_coefficient > 0 || return 0.0
    peak_at_nozzle = _plume_peak_shear(cfg, F / (pi * cfg.nozzle_exit_radius_m^2))
    peak_at_nozzle > cfg.threshold_shear_pa || return 0.0
    cfg.threshold_shear_pa == 0 && return cfg.max_height_m
    k = sqrt(2.0) * exp(-0.5) * cfg.friction_coefficient * F /
        (pi * tand(cfg.plume_half_angle_deg)^2 * cfg.threshold_shear_pa)
    return min(sqrt(k), cfg.max_height_m)
end

"""
Instantaneous phenomenological pressure, peak shear, erosion rate, bounded ejecta speed, sampled annulus and ground-effect force. Erosion uses 128-cell midpoint quadrature to three footprint radii. This function is pure.
"""
function plume_quantities(cfg::PlumeSurfaceConfig, thrust_n::Real, height_m::Real)
    F, h = _plume_inputs(thrust_n, height_m)
    zero_out = (pressure_pa=0.0, shear_pa=0.0, erosion_kg_s=0.0, ejecta_mps=0.0, inner_m=0.0, outer_m=0.0,
                ground_effect_n=0.0)
    (isfinite(F) && F > 0.0 && isfinite(h) && h >= 0.0 && h <= cfg.max_height_m) || return zero_out
    p0, R = plume_surface_footprint(cfg, F, h)
    tau_peak = _plume_peak_shear(cfg, p0)
    ge = plume_ground_effect_force(cfg, F, h)
    tau_t = cfg.threshold_shear_pa
    if !(tau_peak > tau_t)
        return (pressure_pa=p0, shear_pa=tau_peak, erosion_kg_s=0.0, ejecta_mps=0.0, inner_m=0.0, outer_m=0.0,
                ground_effect_n=ge)
    end
    # Ejecta speed: a grain dragged across one footprint radius by the gas
    # dynamic pressure q = τ_peak / c_f, starting from rest.
    q_gas = tau_peak / cfg.friction_coefficient
    accel = 3.0 * cfg.particle_drag_coefficient * q_gas / (4.0 * cfg.particle_density_kg_m3 * cfg.particle_diameter_m)
    v_ej = clamp(sqrt(2.0 * accel * R), cfg.ejecta_speed_min_mps, cfg.ejecta_speed_max_mps)
    # Excess shear force over the annulus, by midpoint quadrature in x = r / R_p.
    A = cfg.friction_coefficient * p0
    dx = _PLUME_QUADRATURE_LIMIT / _PLUME_QUADRATURE_POINTS
    excess_n = 0.0
    inner = Inf
    outer = 0.0
    @inbounds for k in 1:_PLUME_QUADRATURE_POINTS
        x = (k - 0.5) * dx
        tau = A * 2.0 * x * exp(-x * x)
        tau > tau_t || continue
        excess_n += (tau - tau_t) * 2.0 * pi * x * dx * R * R
        inner = min(inner, x * R)
        outer = max(outer, x * R)
    end
    excess_n > 0.0 || return (pressure_pa=p0, shear_pa=tau_peak, erosion_kg_s=0.0, ejecta_mps=0.0,
                              inner_m=0.0, outer_m=0.0, ground_effect_n=ge)
    rate = cfg.erosion_efficiency * excess_n / v_ej
    return (pressure_pa=p0, shear_pa=tau_peak, erosion_kg_s=rate, ejecta_mps=v_ej,
            inner_m=(isfinite(inner) ? inner : 0.0), outer_m=outer, ground_effect_n=ge)
end

"""
Phenomenological axial augmentation (N), decaying exponentially from the configured fraction at contact to zero at ground_effect_cutoff nozzle diameters. This is an illustrative closure, not a calibrated flight prediction.
"""
function plume_ground_effect_force(cfg::PlumeSurfaceConfig, thrust_n::Real, height_m::Real)::Float64
    F, h = _plume_inputs(thrust_n, height_m)
    (isfinite(F) && F > 0.0 && isfinite(h) && h >= 0.0) || return 0.0
    D = 2.0 * cfg.nozzle_exit_radius_m
    D > 0.0 || return 0.0
    x = h / D
    x0 = cfg.ground_effect_cutoff
    x < x0 || return 0.0
    s = cfg.ground_effect_scale
    (s > 0.0 && x0 > 0.0) || return 0.0
    ratio = exp(-x / s) * (-expm1(-(x0 - x) / s)) / (-expm1(-x0 / s))
    return F * cfg.ground_effect_max_fraction * ratio
end

@inline environment_requirements(::PlumeSurfaceInteractionModel) = EffectorEnvironmentRequirements(planet_frame=true)
const _ZERO3 = SVector(0.0, 0.0, 0.0)
const _ZERO_QUANTITIES = (pressure_pa=0.0, shear_pa=0.0, erosion_kg_s=0.0,
    ejecta_mps=0.0, inner_m=0.0, outer_m=0.0, ground_effect_n=0.0)

function plume_engine_axis(x::StateSample)
    radius = norm(x.pos_ii)
    isfinite(radius) && radius > 0 || throw(ArgumentError("Plume position must have finite positive radius."))
    x.q_ib === nothing && return -x.pos_ii / radius
    qnorm = norm(x.q_ib)
    isfinite(qnorm) && qnorm > 0 || throw(ArgumentError("Plume attitude quaternion must be finite and nonzero."))
    A = rot(x.q_ib / qnorm)
    return SVector(A[3, 1], A[3, 2], A[3, 3])
end

function _plume_height_along_axis(model, x, env, axis)
    pf = env.planet_frame
    pf === nothing && throw(ArgumentError("Plume interaction requires a sampled planet frame."))
    radius = norm(pf.pos_pp)
    isfinite(radius) && radius > 0 || throw(ArgumentError("Plume planet-frame radius must be finite and positive."))
    # Terrain uses planetocentric coordinates, not pf.lat_rad (geodetic).
    lat = atand(pf.pos_pp[3], hypot(pf.pos_pp[1], pf.pos_pp[2]))
    lon = atand(pf.pos_pp[2], pf.pos_pp[1])
    reference = model.reference_radius_m === nothing ?
        (model.terrain isa DEMTerrainModel ? model.terrain.reference_radius_m : Float64(env.planet.Rp_e)) :
        model.reference_radius_m
    ground = terrain_radius(model.terrain, lat, lon, reference)
    down_cos = dot(axis, -x.pos_ii / norm(x.pos_ii))
    down_cos > 0 || return Inf
    return max(0.0, radius - ground) / min(down_cos, 1.0)
end

function _plume_sample(model, x, env, t::Float64, i::Int;
        thrust_n::Real=_engine_thrust_n(model, i))
    1 <= i <= length(model.state.height_m) || throw(BoundsError(model.state.height_m, i))
    isfinite(t) || throw(ArgumentError("Plume sample time must be finite."))
    thrust, _ = _plume_inputs(thrust_n, 0.0)
    axis = plume_engine_axis(x)
    height = _plume_height_along_axis(model, x, env, axis)
    quantities = isfinite(height) ? plume_quantities(model.config, thrust,
        max(0.0, height - model.config.nozzle_offset_m)) : _ZERO_QUANTITIES
    return (; height_m=height, quantities..., force_ii=-quantities.ground_effect_n * axis)
end

function wrench(model::PlumeSurfaceInteractionModel, x::StateSample, env::EnvironmentSample, t::Float64)
    length(model.state.height_m) == 1 || throw(ArgumentError(
        "A multi-spacecraft plume needs wrench(model, state, environment, time, sat_idx)."))
    return wrench(model, x, env, t, 1)
end
function wrench(model::PlumeSurfaceInteractionModel, x::StateSample, env::EnvironmentSample, t::Float64, i::Int)
    sample = _plume_sample(model, x, env, t, i)
    return sample.force_ii, _ZERO3
end
function wrench_caching!(model::PlumeSurfaceInteractionModel, x::StateSample,
        env::EnvironmentSample, t::Float64, p, i::Int)
    return wrench(model, x, env, t, i)
end

function _reset_plume_state!(state::PlumeSurfaceState)
    for name in fieldnames(PlumeSurfaceState)
        fill!(getfield(state, name), name in (:last_time_s, :previous_time_s) ? NaN : 0.0)
    end
    return nothing
end
function _cache_plume_sample!(state, i, sample)
    for name in (:height_m, :pressure_pa, :shear_pa, :erosion_kg_s, :ejecta_mps, :ground_effect_n)
        getfield(state, name)[i] = getproperty(sample, name)
    end
    return nothing
end

# Called only by the accepted-step callback. No RHS path calls this helper.
function _accept_plume_sample!(model, x, env, t::Float64, i::Int; active::Bool=true)
    state = model.state
    thrust = active ? _engine_thrust_n(model, i) : 0.0
    sample = _plume_sample(model, x, env, t, i; thrust_n=thrust)
    last = state.last_time_s[i]
    if isfinite(last) && t < last
        throw(ArgumentError("Plume accepted sample time moved backwards."))
    end
    if isfinite(last) && t > last
        # Control callbacks have already updated the actuator at t. Integrate
        # the interval with the previously held thrust, then retain the new
        # held value for the next interval and for the endpoint diagnostic.
        endpoint = _plume_sample(model, x, env, t, i; thrust_n=state.last_thrust_n[i])
        state.previous_time_s[i] = last
        state.previous_rate_kg_s[i] = state.last_rate_kg_s[i]
        state.previous_eroded_kg[i] = state.eroded_kg[i]
        state.previous_thrust_n[i] = state.last_thrust_n[i]
        state.interval_end_rate_kg_s[i] = endpoint.erosion_kg_s
        state.eroded_kg[i] += (t - last) * (state.last_rate_kg_s[i] + endpoint.erosion_kg_s) / 2
    elseif !isfinite(last)
        state.previous_time_s[i] = t
        state.previous_rate_kg_s[i] = sample.erosion_kg_s
        state.previous_eroded_kg[i] = state.eroded_kg[i]
        state.previous_thrust_n[i] = thrust
        state.interval_end_rate_kg_s[i] = sample.erosion_kg_s
    end
    state.last_time_s[i] = t
    state.last_rate_kg_s[i] = sample.erosion_kg_s
    state.last_thrust_n[i] = thrust
    _cache_plume_sample!(state, i, sample)
    return nothing
end

function _plume_saved_thrust(state, i, t)
    t == state.last_time_s[i] && return state.last_thrust_n[i]
    state.previous_time_s[i] <= t < state.last_time_s[i] && return state.previous_thrust_n[i]
    throw(ArgumentError("Plume saves must lie in the current accepted-step interval."))
end
function _plume_mass_at(state, i, t)
    t == state.last_time_s[i] && return state.eroded_kg[i]
    before, after = state.previous_time_s[i], state.last_time_s[i]
    before <= t < after || throw(ArgumentError("Plume mass saves must lie in the current accepted-step interval."))
    dt = t - before
    fraction = dt / (after - before)
    rate0 = state.previous_rate_kg_s[i]
    rate = rate0 + fraction * (state.interval_end_rate_kg_s[i] - rate0)
    return state.previous_eroded_kg[i] + dt * (rate0 + rate) / 2
end
end # module PlumeSurfaceInteraction
