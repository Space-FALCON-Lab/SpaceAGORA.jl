# Plume-surface interaction of a descent engine with the regolith: the surface
# pressure and shear stress the exhaust lays down under the vehicle, the soil
# it erodes, the speed the grains leave at, and the small thrust augmentation
# the reflected plume gives the vehicle in ground effect.
#
# The gas-dynamic part follows Roberts' treatment of a hypersonic jet acting on
# a dust layer (L. Roberts, "The action of a hypersonic jet on a dust layer",
# IAS Paper 63-50, 1963) in the form used for the Moon by P. T. Metzger and
# co-workers (Metzger, Immer, Donahue, Vu, Latta, Deyo-Svendsen, "Jet-induced
# cratering of a granular surface with application to lunar spaceports",
# J. Aerosp. Eng. 22, 2009; Metzger, Smith, Lane, "Phenomenology of soil
# erosion due to rocket exhaust on the Moon and the Mauna Kea lunar test site",
# J. Geophys. Res. 116, E06005, 2011): a momentum-conserving surface-pressure
# footprint, a wall shear stress proportional to it, and viscous erosion driven
# by the shear stress in excess of the soil's threshold.
#
# Two parameters are calibrated, not derived, and both are anchored on Apollo 11
# observables (see `docs/src/user/lunar_landing.md`):
#   * `threshold_shear_pa` fixes the height at which erosion starts;
#   * `erosion_efficiency` fixes how much soil the descent moves in total.
module PlumeSurfaceInteraction

using StaticArrays
using LinearAlgebra
using ...AbstractTypes: AbstractForceTorqueModel, AbstractTerrainModel
using ...TerrainModels: NoTerrainModel, DEMTerrainModel, terrain_height
using ...EffectorSampling: StateSample, EnvironmentSample, EffectorEnvironmentRequirements
import ...SimulationModel
using ...SimulationModel: rot
import ..DynamicEffectors: wrench, wrench_caching!, environment_requirements

export PlumeSurfaceConfig, PlumeSurfaceInteractionModel, PlumeSurfaceState
export plume_surface_footprint, plume_erosion_onset_height, plume_ground_effect_force, plume_quantities

"""
    PlumeSurfaceConfig(; ...)

Nozzle, plume and regolith properties of a [`PlumeSurfaceInteractionModel`](@ref).
The defaults are the Apollo lunar module's descent propulsion system (DPS) over
mare regolith: a 1.5 m exit diameter, an area ratio of 47.5, a chamber pressure
of 7.2 bar and 45.04 kN at full throttle (throttleable to 10 percent, so about
4.5 kN to 45 kN), over soil of 1500 kg/m³ bulk density made of 70 µm grains.

Gas dynamics:

- `nozzle_exit_radius_m`, `expansion_ratio`, `chamber_pressure_pa` and
  `exit_mach` describe the engine. Only the exit radius enters the force model
  (it sets the ground-effect length scale and floors the plume footprint); the
  rest are carried so a scenario records the engine it flew.
- `nozzle_offset_m` is how far the nozzle exit plane sits below the vehicle's
  reference point along the engine axis (1.5 m for the LM, where the state
  point is level with the descent stage deck). The plume geometry is measured
  from the exit plane; the `height_m` the model records is measured from the
  reference point, which is the point the trajectory is integrated at and the
  point the viewer draws the plume from.
- `plume_half_angle_deg` is the half-angle of the momentum-carrying core of the
  vacuum plume. The surface pressure is spread over a Gaussian footprint of
  radius `R_p = h tan(θ_p)` normalized so its integral is the engine thrust.
- `friction_coefficient` is the wall skin-friction coefficient that turns the
  local surface pressure into wall shear stress.

Regolith:

- `bulk_density_kg_m3`, `particle_density_kg_m3` and `particle_diameter_m` are
  the soil properties; `cohesion_pa` is recorded for reference.
- `threshold_shear_pa` is the wall shear stress below which nothing moves. It
  is the parameter that sets the erosion onset height and defaults to 0.15 Pa,
  which puts the onset at about 31 m for the Apollo 11 approach thrust — the
  height at which the crew first reported blowing dust. It is three to four
  orders of magnitude below the bulk cohesion of lunar regolith (0.1–1 kPa),
  as the mobile surface layer must be.
- `erosion_efficiency` multiplies the momentum-balance erosion rate to account
  for the saltation cascade (each impacting grain splashes several more), which
  a pure momentum balance cannot produce. It defaults to 10.
- `particle_drag_coefficient`, `ejecta_speed_min_mps` and
  `ejecta_speed_max_mps` bound the ejecta speed.

Ground effect (see [`plume_ground_effect_force`](@ref)):
`ground_effect_max_fraction`, `ground_effect_scale` and
`ground_effect_cutoff` give the thrust augmentation inside
`ground_effect_cutoff` exit diameters of the ground.

`max_height_m` short-circuits the whole model: above it every plume quantity
is zero.
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
end

"""
    PlumeSurfaceState(num_sats)

Per-spacecraft record the effector keeps as it runs, so a run with
`isolate_state=false` can read it afterwards and the save fields can publish
it: the height of the engine above the ground along the engine axis, the peak
surface pressure and wall shear stress, the mass erosion rate and its time
integral, the characteristic ejecta speed and the ground-effect force. The
last two vectors carry the state of the cumulative integral.
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
end

function PlumeSurfaceState(num_sats::Integer)
    n = Int(num_sats)
    n >= 1 || throw(ArgumentError("PlumeSurfaceState needs at least one spacecraft"))
    return PlumeSurfaceState(zeros(n), zeros(n), zeros(n), zeros(n), zeros(n), zeros(n), zeros(n), fill(NaN, n), zeros(n))
end

"""
    PlumeSurfaceInteractionModel(control, terrain; config=PlumeSurfaceConfig(), num_sats=...)

Dynamic effector for the plume-surface interaction of a descent engine.

`control` is the descent control effector that owns the engine: the model reads
the engine's actual thrust from its `actuators.thrust_n` vector every
evaluation, so the plume follows the throttle the controller commanded rather
than the guidance demand. `terrain` is the terrain model the height above
ground is measured against (the same one the guidance uses as its altimeter).

The force it returns is the ground-effect thrust augmentation along the engine
axis; it applies no torque. The erosion, pressure, shear and ejecta quantities
are diagnostics kept in [`PlumeSurfaceState`](@ref) and published as result
columns `sc{i}_plume_*` when the effector is present.

```julia
control = ApolloDescentControlModel(ApolloDescentControlConfig(), gcfg, state, terrain)
plume = PlumeSurfaceInteractionModel(control, terrain)
# ... dynamic_effectors = (gravity..., plume)
```
"""
struct PlumeSurfaceInteractionModel{C, T <: AbstractTerrainModel} <: AbstractForceTorqueModel
    config::PlumeSurfaceConfig
    control::C
    terrain::T
    state::PlumeSurfaceState
end

function PlumeSurfaceInteractionModel(control, terrain::AbstractTerrainModel=NoTerrainModel();
                                      config::PlumeSurfaceConfig=PlumeSurfaceConfig(),
                                      num_sats::Integer=_control_spacecraft_count(control))
    return PlumeSurfaceInteractionModel(config, control, terrain, PlumeSurfaceState(num_sats))
end

@inline function _control_spacecraft_count(control)::Int
    hasproperty(control, :actuators) && hasproperty(control.actuators, :thrust_n) ||
        throw(ArgumentError("PlumeSurfaceInteractionModel needs a descent control effector with an `actuators.thrust_n` vector."))
    return length(control.actuators.thrust_n)
end

@inline function _engine_thrust_n(model::PlumeSurfaceInteractionModel, i::Int)::Float64
    thrusts = model.control.actuators.thrust_n
    return (1 <= i <= length(thrusts)) ? max(0.0, Float64(thrusts[i])) : 0.0
end

# ---- the plume on the ground ---------------------------------------------------------

const _PLUME_QUADRATURE_POINTS = 128
const _PLUME_QUADRATURE_LIMIT = 3.0        # footprint radii; exp(-9) is already negligible

"Peak wall shear stress of the Gaussian footprint, at r = R_p / sqrt(2)."
@inline _plume_peak_shear(cfg::PlumeSurfaceConfig, p0::Float64)::Float64 = cfg.friction_coefficient * p0 * sqrt(2.0) * exp(-0.5)

"""
    plume_surface_footprint(config, thrust_n, height_m) -> (p0_pa, radius_m)

Peak (stagnation) surface pressure and the footprint radius `R_p = h tan(θ_p)`
of the plume of `thrust_n` newtons standing `height_m` above the ground. The
pressure is `p(r) = p0 exp(-(r/R_p)^2)`, normalized so `∫ p dA` is the thrust:
all of the engine's axial momentum is turned by the surface. The radius never
falls below the nozzle exit radius.
"""
@inline function plume_surface_footprint(cfg::PlumeSurfaceConfig, thrust_n::Real, height_m::Real)
    R = max(Float64(height_m) * tand(cfg.plume_half_angle_deg), cfg.nozzle_exit_radius_m)
    p0 = Float64(thrust_n) / (pi * R * R)
    return p0, R
end

"""
    plume_erosion_onset_height(config, thrust_n) -> Float64

Height (m) at which the peak wall shear stress of a `thrust_n` plume falls to
the soil's threshold: erosion starts below it and is identically zero above it.
With the defaults and the Apollo 11 approach thrust (about 11.5 kN of lunar
weight near the end of the descent) this is about 31 m, matching the roughly
30 m at which the Apollo 11 crew first reported blowing dust.
"""
function plume_erosion_onset_height(cfg::PlumeSurfaceConfig, thrust_n::Real)::Float64
    F = Float64(thrust_n)
    (F > 0.0 && cfg.threshold_shear_pa > 0.0) || return 0.0
    # τ_peak = 0.8578 c_f F / (π h² tan²θ) = τ_t
    k = sqrt(2.0) * exp(-0.5) * cfg.friction_coefficient * F / (pi * tand(cfg.plume_half_angle_deg)^2 * cfg.threshold_shear_pa)
    return sqrt(k)
end

"""
    plume_quantities(config, thrust_n, height_m) -> NamedTuple

The whole plume-surface state at one instant: `pressure_pa` and `shear_pa` are
the peaks of the surface distributions, `erosion_kg_s` is the mass erosion rate
integrated over the eroded annulus, `ejecta_mps` the characteristic speed the
grains leave at, `inner_m` and `outer_m` the edges of the annulus where the
shear stress exceeds the threshold, and `ground_effect_n` the thrust
augmentation. Everything but the pressure and the shear stress is zero above
the erosion onset height.

The erosion closure is Roberts' viscous erosion written as a momentum balance:
the shear stress in excess of the threshold, integrated over the annulus, is
the force available to accelerate grains, and a mass flux `ṁ` leaving at
`v_ej` carries `ṁ v_ej` of momentum, so

```math
\\dot m = \\frac{\\eta}{v_{ej}} \\int \\max(\\tau(r) - \\tau_t,\\, 0)\\, \\mathrm{d}A
```

with `η = erosion_efficiency` covering the saltation cascade. The ejecta speed
comes from a drag balance on one grain accelerated across the footprint by the
gas dynamic pressure at the shear peak.
"""
function plume_quantities(cfg::PlumeSurfaceConfig, thrust_n::Real, height_m::Real)
    F = Float64(thrust_n)
    h = Float64(height_m)
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
    plume_ground_effect_force(config, thrust_n, height_m) -> Float64

Thrust augmentation (N) the reflected plume gives the vehicle in ground effect,
along the engine axis. The correlation is an exponential in the height over the
nozzle exit diameter,

```math
\\Delta F = f_{max} F \\frac{e^{-x/s} - e^{-x_0/s}}{1 - e^{-x_0/s}},\\qquad x = h/D_e,
```

zero at and above `x_0 = ground_effect_cutoff` exit diameters and rising
monotonically to `f_max = ground_effect_max_fraction` of the thrust at contact.
The form follows the exponential decay of the base-pressure rise measured for
nozzles near a plate; the magnitude (3 percent at contact, vanishing by two
exit diameters) is a modeling choice sized to the few-percent effect reported
for lunar-lander-class plumes (Roberts 1963; Metzger et al. 2011), not an
Apollo flight measurement.
"""
function plume_ground_effect_force(cfg::PlumeSurfaceConfig, thrust_n::Real, height_m::Real)::Float64
    F = Float64(thrust_n)
    h = Float64(height_m)
    (isfinite(F) && F > 0.0 && isfinite(h) && h >= 0.0) || return 0.0
    D = 2.0 * cfg.nozzle_exit_radius_m
    D > 0.0 || return 0.0
    x = h / D
    x0 = cfg.ground_effect_cutoff
    x < x0 || return 0.0
    s = cfg.ground_effect_scale
    (s > 0.0 && x0 > 0.0) || return 0.0
    tail = exp(-x0 / s)
    return F * cfg.ground_effect_max_fraction * (exp(-x / s) - tail) / (1.0 - tail)
end

# ---- effector interface --------------------------------------------------------------

@inline environment_requirements(::PlumeSurfaceInteractionModel) = EffectorEnvironmentRequirements(planet_frame=true)

@inline function _plume_reference_radius(model::PlumeSurfaceInteractionModel, planet)::Float64
    return model.terrain isa DEMTerrainModel ? model.terrain.reference_radius_m : Float64(planet.Rp_e)
end

"""
    plume_engine_axis(x) -> SVector{3}

Unit vector of the engine axis in the inertial frame: the descent engine
thrusts along body `-z`, so the plume travels along body `+z`. Without an
attitude state the axis is taken straight down, along the local vertical.
"""
@inline function plume_engine_axis(x::StateSample)::SVector{3, Float64}
    q = x.q_ib
    if q === nothing
        r = x.pos_ii
        n = norm(r)
        return n > 0.0 ? -r / n : SVector{3, Float64}(0.0, 0.0, -1.0)
    end
    A = rot(SVector{4, Float64}(q))          # rows: body axes in inertial coordinates
    return SVector{3, Float64}(A[3, 1], A[3, 2], A[3, 3])
end

function wrench(model::PlumeSurfaceInteractionModel, x::StateSample, env::EnvironmentSample, t::Float64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    return _plume_wrench(model, x, env, t, 1)
end

function wrench_caching!(model::PlumeSurfaceInteractionModel, x::StateSample, env::EnvironmentSample, t::Float64, p, sat_idx::Int)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    return _plume_wrench(model, x, env, t, sat_idx)
end

function _plume_wrench(model::PlumeSurfaceInteractionModel, x::StateSample, env::EnvironmentSample, t::Float64, i::Int)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    zero3 = SVector{3, Float64}(0.0, 0.0, 0.0)
    st = model.state
    (1 <= i <= length(st.height_m)) || return zero3, zero3
    axis = plume_engine_axis(x)
    thrust = _engine_thrust_n(model, i)
    height = _plume_height_along_axis(model, x, env, axis)
    q = plume_quantities(model.config, thrust, max(0.0, height - model.config.nozzle_offset_m))
    st.height_m[i] = height
    st.pressure_pa[i] = q.pressure_pa
    st.shear_pa[i] = q.shear_pa
    st.erosion_kg_s[i] = q.erosion_kg_s
    st.ejecta_mps[i] = q.ejecta_mps
    st.ground_effect_n[i] = q.ground_effect_n
    _accumulate_eroded_mass!(st, i, t, q.erosion_kg_s)
    # The augmentation pushes the vehicle away from the ground, along the
    # thrust direction (the engine axis is the direction the plume travels).
    return -q.ground_effect_n * axis, zero3
end

"""
    _plume_height_along_axis(model, x, env, axis) -> Float64

Slant distance (m) from the vehicle to the ground along the engine axis: the
height of the vehicle above the terrain divided by the cosine of the angle
between the engine axis and the local vertical. The footprint is treated as if
the plume met the surface square on, so a tilted vehicle only moves the
impingement point, it does not skew the pressure distribution.
"""
function _plume_height_along_axis(model::PlumeSurfaceInteractionModel, x::StateSample, env::EnvironmentSample, axis::SVector{3, Float64})::Float64
    pf = env.planet_frame
    pf === nothing && return NaN
    radius = norm(pf.pos_pp)
    ground = _plume_reference_radius(model, env.planet) + terrain_height(model.terrain, rad2deg(pf.lat_rad), rad2deg(pf.lon_rad))
    vertical = radius - ground
    vertical > 0.0 || return 0.0
    r = x.pos_ii
    n = norm(r)
    n > 0.0 || return vertical
    cosang = clamp(dot(axis, -r / n), 0.1, 1.0)     # engine axis versus straight down
    return vertical / cosang
end

"""
    _accumulate_eroded_mass!(state, i, t, rate)

Trapezoidal time integral of the erosion rate, guarded so it only ever advances
in time. The solver evaluates the right-hand side at stage times inside a trial
step and repeats them when a step is rejected; integrating only when the
evaluation time exceeds every time already integrated keeps the cumulative mass
monotone and free of the double counting a naive integral would pick up.
"""
@inline function _accumulate_eroded_mass!(st::PlumeSurfaceState, i::Int, t::Float64, rate::Float64)
    isfinite(t) || return nothing
    last_t = st.last_time_s[i]
    if !isfinite(last_t)
        st.last_time_s[i] = t
        st.last_rate_kg_s[i] = rate
        return nothing
    end
    if t > last_t
        st.eroded_kg[i] += 0.5 * (st.last_rate_kg_s[i] + rate) * (t - last_t)
        st.last_time_s[i] = t
        st.last_rate_kg_s[i] = rate
    end
    return nothing
end

end # module PlumeSurfaceInteraction
