# Plume-surface interaction of a descent engine with the regolith: the surface
# pressure and shear stress the exhaust lays down under the vehicle, the soil
# it erodes, the crater it digs, where the grains land, and the small thrust
# augmentation the reflected plume gives the vehicle in ground effect.
#
# The model is assembled from four pieces, each of which is its own file and
# each of which can be swapped or evaluated on its own:
#
#   * the GAS ON THE GROUND comes from a plume field
#     (`plume_gas_field.jl`), which the configuration carries and which answers
#     `plume_gas_state(field, config, thrust, height, radius)`.
#     `PlumeAnalyticField`, the default, is the momentum-conserving Gaussian
#     pressure footprint with a skin-friction shear law; `PlumeFieldTable` is
#     the Simons source-flow plume of a named engine, tabulated by
#     `scripts/dev/psi/build_plume_field.jl`. Supply one with
#     `PlumeSurfaceConfig(field=load_plume_field("data/psi/apollo_lmde.json"))`.
#   * the SOIL MOVED at each radius comes from the erosion regimes
#     (`regolith_erosion.jl`): viscous erosion in Metzger's energy-flux form,
#     diffusion-driven flow and bearing-capacity failure, each with its own
#     onset criterion derived from the soil rather than fitted. The effector
#     integrates their local rate over the footprint and records which regime
#     moved the most mass. This is the default (`erosion_model = :regimes`).
#   * the OLD LAW is still reachable as `erosion_model = :roberts_fitted`: the
#     shear-excess momentum balance of L. Roberts, "The action of a hypersonic
#     jet on a dust layer", IAS Paper 63-50, 1963, with the two calibrated
#     constants `threshold_shear_pa = 0.15` (onset at 31 m with the analytic
#     field) and `erosion_efficiency = 10` (the saltation cascade). It is kept
#     only so the two laws can be compared on one gas state; Metzger (2024a,
#     conclusions) holds that its form is wrong.
#   * WHERE THE GRAINS GO comes from the ejecta transport module
#     (`ejecta_transport.jl`), evaluated once per saved sample rather than per
#     right-hand-side evaluation, weighted by the regimes' own local erosion
#     rate and by a lognormal grain-mass distribution fitted to the soil's D50
#     and D84/D50.
#
# The CRATER is this file's own: a radial grid of eroded depth under the
# vehicle, advanced with the same monotone time guard as the cumulative eroded
# mass, whose depth under the stagnation point feeds back into the height the
# field is queried at. Its deepest point and its edge are saved quantities, and
# the whole profile is readable through `plume_crater_profile`.
#
# WHAT IS CALIBRATED. With the default `:regimes` law, nothing in the erosion
# path is fitted inside this file: every coefficient belongs to
# `RegolithProperties` and is sourced there. What remains a modeling choice here
# is the ground-effect correlation (see `plume_ground_effect_force`), the crater
# grid's geometry, and the gas residence time that sets the pore-pressure
# diffusion depth; all three are named configuration fields.
# `docs/src/user/lunar_landing.md` states the same list.
module PlumeSurfaceInteraction

using StaticArrays
using LinearAlgebra
using ...AbstractTypes: AbstractForceTorqueModel, AbstractTerrainModel
using ...TerrainModels: NoTerrainModel, DEMTerrainModel, terrain_height
using ...EffectorSampling: StateSample, EnvironmentSample, EffectorEnvironmentRequirements
import ...SimulationModel
using ...SimulationModel: rot
using ..PlumeGasField: PlumeGasState, PlumeAnalyticField, PlumeFieldTable, plume_gas_state,
                       plume_field_footprint, plume_field_shear_coefficient, _plume_shear_at,
                       _plume_quadrature_limit
using ..RegolithErosion: RegolithProperties, lunar_mare_regolith, ErosionEnvironment, ViscousErosionRoberts,
                         ErosionRegimeKind, NoErosion, ViscousErosion,
                         DiffusionDrivenFlowRegime, BearingCapacityFailureRegime,
                         DiffusionDrivenFlow, BearingCapacityFailure,
                         default_erosion_regimes, regolith_erosion_rate
using ..EjectaTransport: EjectaTransportConfig, ejecta_distribution, ejecta_lognormal_mass_weights
import ..DynamicEffectors: wrench, wrench_caching!, environment_requirements

export PlumeSurfaceConfig, PlumeSurfaceInteractionModel, PlumeSurfaceState
export plume_surface_footprint, plume_erosion_onset_height, plume_ground_effect_force, plume_quantities
export plume_crater_profile, plume_ejecta_summary, plume_refresh_ejecta!, plume_regime_code
export plume_default_regimes

"""
    PlumeSurfaceConfig(; ...)

Nozzle, plume and regolith properties of a [`PlumeSurfaceInteractionModel`](@ref).
The defaults are the Apollo lunar module's descent propulsion system (DPS) over
mare regolith: a 1.5 m exit diameter, an area ratio of 47.5, a chamber pressure
of 7.2 bar and 45.04 kN at full throttle (throttleable to 10 percent, so about
4.5 kN to 45 kN), over soil of 1500 kg/m³ bulk density made of 70 µm grains.

Plume field:

- `field` is the object every surface quantity is read from: a
  `PlumeAnalyticField` (the default, the Gaussian footprint described below) or
  a `PlumeFieldTable` loaded with `load_plume_field`. Swapping it changes the
  pressure, shear stress, footprint radius and everything derived from them,
  and nothing else; the default keeps earlier runs bit for bit.

Gas dynamics:

- `nozzle_exit_radius_m`, `expansion_ratio`, `chamber_pressure_pa` and
  `exit_mach` describe the engine. Only the exit radius enters the force model
  (it sets the ground-effect length scale and floors the plume footprint); the
  rest are carried so a scenario records the engine it flew. `exit_mach`
  defaults to 5.03, the isentropic exit Mach number of an area ratio of 47.5 at
  the exhaust's ratio of specific heats, which is also the value Morris 2012
  Table 3.1.1 reports for this engine (see `PlumeNozzle`); it replaces an
  earlier unsourced 4.8, and no computed quantity depends on it.
- `nozzle_offset_m` is how far the nozzle exit plane sits below the vehicle's
  reference point along the engine axis (1.5 m for the LM, where the state
  point is level with the descent stage deck). The plume geometry is measured
  from the exit plane; the `height_m` the model records is measured from the
  reference point, which is the point the trajectory is integrated at and the
  point the viewer draws the plume from.
- `plume_half_angle_deg` is the half-angle of the momentum-carrying core of the
  vacuum plume, used by `PlumeAnalyticField` only: the surface pressure is
  spread over a Gaussian footprint of radius `R_p = h tan(θ_p)` normalized so
  its integral is the engine thrust.
- `friction_coefficient` is the wall skin-friction coefficient that turns the
  local surface pressure into wall shear stress, again for `PlumeAnalyticField`
  only. A table carries its own surface drag coefficient.

Erosion:

- `erosion_model` selects the law the local mass flux comes from.
  `:regimes` (the default) integrates
  [`regolith_erosion_rate`](@ref)'s regimes over the footprint on the gas state
  the field reports at every radius, so the active regime is whichever one's
  own onset criterion is met; nothing in it is fitted here.
  `:roberts_fitted` is the shear-excess momentum balance this model used
  before, with the two calibrated constants below, kept so the two laws can be
  compared on one gas state.
- `soil` is the [`RegolithProperties`](@ref) the regimes read (grain sizes,
  strength, permeability, and Metzger's sourced erosion coefficients); it
  defaults to `lunar_mare_regolith()`.
- `regimes` is the tuple of regime singletons evaluated at each radius,
  `default_erosion_regimes()` by default. It is a type parameter, so the loop
  over it unrolls and allocates nothing.
- `gravity_m_s2` is the surface gravity the regimes weigh the soil against when
  a caller evaluates the pure functions directly; inside a run the effector
  overrides it with the local `mu/r^2` of the planet it is flying over. The
  default 1.625 m/s² is the Moon's mean surface gravity.
- `gas_residence_time_s` is how long the plume is taken to have been loading the
  patch of ground being evaluated, which sets the pore-pressure diffusion depth
  of the diffusion-driven-flow regime. ASSUMPTION, not a measurement: 1 s is the
  timescale over which the surface pressure at a fixed radius changes at the
  meter-per-second descent rates of the last few meters.

Regolith (the legacy `:roberts_fitted` path and the ground geometry):

- `bulk_density_kg_m3`, `particle_density_kg_m3` and `particle_diameter_m` are
  the soil properties; `cohesion_pa` is recorded for reference. The `:regimes`
  law reads `soil` instead, and `bulk_density_kg_m3` is what the crater's depth
  is computed from.
- `threshold_shear_pa` is the wall shear stress below which nothing moves under
  the `:roberts_fitted` law. It defaults to 0.15 Pa, which puts that law's onset
  at about 31 m for the Apollo 11 approach thrust with the analytic field — the
  height at which the crew first reported blowing dust. **It is a calibrated
  constant of that law only**; the `:regimes` law derives its thresholds from
  the soil and ignores this field entirely.
- `erosion_efficiency` multiplies the momentum-balance erosion rate of the
  `:roberts_fitted` law to account for the saltation cascade. It defaults to 10
  and, like `threshold_shear_pa`, is unused under `:regimes`.
- `particle_drag_coefficient`, `ejecta_speed_min_mps` and
  `ejecta_speed_max_mps` bound the characteristic ejecta speed both laws report.

Crater (see [`plume_crater_profile`](@ref)):

- `crater_bins` (64), `crater_min_radius_m` (0.05 m) and `crater_max_radius_m`
  (40 m) size the radial grid of eroded depth kept under the vehicle. The nodes
  are logarithmic in radius, because the crater is toroidal and its inner wall
  sits inside the first meter while its outer edge reaches tens of meters; a
  uniform grid that covers the second cannot resolve the first. The grid is
  axisymmetric about the current impingement point and does not follow it as the
  vehicle translates:
  the model assumes the vehicle is nearly stationary over the ground while it
  is eroding, which is true of the last tens of meters of a landing and is
  reported in the demo's output as the horizontal travel below the onset height.
- `crater_edge_fraction` (0.1) and `crater_min_depth_m` (1 mm) define the
  crater's edge, and so `crater_radius_m`: the outermost radius whose depth is
  still a tenth of the deepest point, and at least the floor. Both are
  definitions, not measurements. The floor exists because a pure ratio test
  returns a wide radius the instant the plume touches the ground, when the
  surface has lost a few grain layers everywhere and nothing that could be
  called a crater exists; one millimeter is about fourteen median grain
  diameters. Below the floor no radius is reported at all. Neither field changes
  any other quantity.
- `crater_height_feedback` (true) adds the depth eroded under the stagnation
  point to the height the plume field is queried at, so a deepening crater moves
  the ground away from the nozzle.

Ejecta (see [`plume_ejecta_summary`](@ref)):

- `ejecta_diagnostic` (true) computes the grain-transport distribution once per
  saved sample while erosion is active. It is off the right-hand side: no force
  depends on it.
- `ejecta` is the [`EjectaTransportConfig`](@ref) the transport model uses.
- `ejecta_sizes` are the grain diameters swept, and `ejecta_radii` the number of
  launch radii spread over the eroding annulus. Both are cost knobs: the default
  5 sizes by 8 radii is 40 trajectories per saved sample.

Ground effect (see [`plume_ground_effect_force`](@ref)):
`ground_effect_max_fraction`, `ground_effect_scale` and
`ground_effect_cutoff` give the thrust augmentation inside
`ground_effect_cutoff` exit diameters of the ground.

`max_height_m` short-circuits the whole model: above it every plume quantity
is zero.
"""
Base.@kwdef struct PlumeSurfaceConfig{F, R}
    nozzle_exit_radius_m::Float64 = 0.75
    nozzle_offset_m::Float64 = 1.5
    expansion_ratio::Float64 = 47.5
    chamber_pressure_pa::Float64 = 7.2e5
    exit_mach::Float64 = 5.03
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
    field::F = PlumeAnalyticField()
    erosion_model::Symbol = :regimes
    soil::RegolithProperties = lunar_mare_regolith()
    regimes::R = plume_default_regimes()
    gravity_m_s2::Float64 = 1.625
    gas_residence_time_s::Float64 = 1.0
    crater_bins::Int = 64
    crater_min_radius_m::Float64 = 0.05
    crater_max_radius_m::Float64 = 40.0
    crater_edge_fraction::Float64 = 0.1
    crater_min_depth_m::Float64 = 1.0e-3
    crater_height_feedback::Bool = true
    ejecta_diagnostic::Bool = true
    ejecta::EjectaTransportConfig = EjectaTransportConfig()
    ejecta_sizes::Vector{Float64} = [5.0e-6, 2.0e-5, 7.0e-5, 2.0e-4, 5.0e-4]
    ejecta_radii::Int = 8
end

"""
    PlumeSurfaceState(num_sats; crater_bins=64, crater_min_radius_m=0.05, crater_max_radius_m=40.0)

Per-spacecraft record the effector keeps as it runs, so a run with
`isolate_state=false` can read it afterwards and the save fields can publish
it.

Per-evaluation quantities: the height of the engine above the ground along the
engine axis, the peak surface pressure and wall shear stress, the mass erosion
rate and its time integral, the characteristic ejecta speed, the ground-effect
force, `regime`, the [`ErosionRegimeKind`](@ref) that moved the most mass over
the footprint as its integer code (see [`plume_regime_code`](@ref)), and
`erosion_radius_m`, the outer edge of the region that moved anything.

Crater: `crater_profile_m` is the eroded depth on the radial grid
`crater_radii_m` (logarithmic in radius, bins down the rows, spacecraft across
the columns),
`crater_depth_m` its deepest point and `crater_radius_m` the outermost radius
still at `crater_edge_fraction` of that depth.

Ejecta, refreshed once per saved sample by [`plume_refresh_ejecta!`](@ref) and
not by the right-hand side: the mass-weighted mean ejection angle above the
local horizontal, the mass-weighted mean deposition radius, and the mass
fraction leaving faster than escape speed.

The trailing vectors carry the state the time integrals need: the last time and
rate integrated, the last local rate on the crater grid, the thrust, height,
gravity, body radius and eroding annulus the ejecta refresh re-uses, and the
time the ejecta summary was last computed at.
"""
mutable struct PlumeSurfaceState
    height_m::Vector{Float64}
    pressure_pa::Vector{Float64}
    shear_pa::Vector{Float64}
    erosion_kg_s::Vector{Float64}
    eroded_kg::Vector{Float64}
    ejecta_mps::Vector{Float64}
    ground_effect_n::Vector{Float64}
    regime::Vector{Float64}
    erosion_radius_m::Vector{Float64}
    crater_depth_m::Vector{Float64}
    crater_radius_m::Vector{Float64}
    ejecta_angle_deg::Vector{Float64}
    ejecta_range_m::Vector{Float64}
    ejecta_escape_frac::Vector{Float64}
    crater_radii_m::Vector{Float64}
    crater_profile_m::Matrix{Float64}
    crater_rate_kg_m2_s::Matrix{Float64}
    crater_scratch_kg_m2_s::Matrix{Float64}
    last_time_s::Vector{Float64}
    last_rate_kg_s::Vector{Float64}
    last_thrust_n::Vector{Float64}
    last_query_height_m::Vector{Float64}
    last_gravity_m_s2::Vector{Float64}
    last_body_radius_m::Vector{Float64}
    last_inner_m::Vector{Float64}
    last_outer_m::Vector{Float64}
    ejecta_time_s::Vector{Float64}
end

function PlumeSurfaceState(num_sats::Integer; crater_bins::Integer=64,
                           crater_min_radius_m::Real=0.05, crater_max_radius_m::Real=40.0)
    n = Int(num_sats)
    n >= 1 || throw(ArgumentError("PlumeSurfaceState needs at least one spacecraft"))
    nb = max(Int(crater_bins), 2)
    rmin = Float64(crater_min_radius_m)
    rmax = Float64(crater_max_radius_m)
    (rmin > 0.0 && rmax > rmin) ||
        throw(ArgumentError("PlumeSurfaceState needs 0 < crater_min_radius_m < crater_max_radius_m"))
    radii = collect(exp.(range(log(rmin), log(rmax); length=nb)))
    return PlumeSurfaceState(zeros(n), zeros(n), zeros(n), zeros(n), zeros(n), zeros(n), zeros(n),
                             zeros(n), zeros(n), zeros(n), zeros(n), zeros(n), zeros(n), zeros(n),
                             radii, zeros(nb, n), zeros(nb, n), zeros(nb, n),
                             fill(NaN, n), zeros(n), zeros(n), zeros(n), fill(1.625, n),
                             fill(1.7374e6, n), zeros(n), zeros(n), fill(NaN, n))
end

"""
    plume_default_regimes() -> Tuple

The erosion regimes a `PlumeSurfaceConfig` evaluates unless it is given others:
`(ViscousErosionRoberts(), DiffusionDrivenFlow(), BearingCapacityFailure())`.

This differs in one place from `default_erosion_regimes()` in
`regolith_erosion.jl`, and the reason is a measurement, not a preference. That
module's default viscous closure is `ViscousErosionEnergyFlux`, Metzger's
energy-flux law (2024a equation 16), whose threshold `E_th = 0.123 J/(m^2 s)`
was fitted to the Apollo 16 landing video **through Metzger's own plume model**.
Evaluated on the surface gas state either of this repository's plume fields
reports, that threshold is crossed by three to four orders of magnitude: at the
Apollo 12 approach thrust the energy flux at the shear peak is 4.3 W/m^2 at
31.9 m and 41 W/m^2 at 10 m, against an `E_th` of 0.123, and integrating the
resulting local rate over the footprint gives 4e4 kg/s on the analytic field and
6e6 kg/s on the table, at every height, against the 11 to 99 kg/s Lane and
Metzger measured. The rate also comes out nearly independent of height, which is
the degeneracy `ViscousErosionEnergyFlux`'s own docstring predicts for a density
that tracks the pressure.

The law is not being rejected; the pairing is. `E = 3 tau^2/(rho_s vbar)` needs
the density in the laminar sublayer at the wall, and what both fields report is
the post-shock free-stream density, which is the wrong quantity by a factor
these fields cannot supply. Until a field computes the sublayer state, the
energy-flux closure cannot be driven from one, and using it anyway would be
reporting a number that is wrong by 10^4 because two sourced pieces were bolted
together without checking that they meet.

`ViscousErosionRoberts` with the derived Shields threshold is used instead: it
needs only the wall shear stress, which both fields do compute, and it lands
within a factor of 3.3 of every one of the eleven Apollo 12 altitudes Lane and
Metzger published (see `benchmarks/studies/psi_validation/README.md`). It still
carries the soil's unsourced `saltation_efficiency`; that is the model's
remaining fitted constant and is reported as such.

Pass `regimes=default_erosion_regimes()` to run the energy-flux closure anyway.
"""
@inline plume_default_regimes() = (ViscousErosionRoberts(), DiffusionDrivenFlow(), BearingCapacityFailure())

"""
    plume_regime_code(kind) -> Float64

Integer code of an [`ErosionRegimeKind`](@ref) as the `sc{i}_plume_regime`
column carries it: 0 no erosion, 1 viscous erosion, 2 diffusion-driven flow,
3 bearing-capacity failure.
"""
@inline plume_regime_code(kind::ErosionRegimeKind)::Float64 = Float64(Int(kind))

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

# the same descent with the tabulated plume of the LM descent engine
cfg = PlumeSurfaceConfig(field=load_plume_field("data/psi/apollo_lmde.json"))
plume = PlumeSurfaceInteractionModel(control, terrain; config=cfg)
```
"""
struct PlumeSurfaceInteractionModel{F, R, C, T <: AbstractTerrainModel} <: AbstractForceTorqueModel
    config::PlumeSurfaceConfig{F, R}
    control::C
    terrain::T
    state::PlumeSurfaceState
end

function PlumeSurfaceInteractionModel(control, terrain::AbstractTerrainModel=NoTerrainModel();
                                      config::PlumeSurfaceConfig=PlumeSurfaceConfig(),
                                      num_sats::Integer=_control_spacecraft_count(control))
    config.erosion_model in (:regimes, :roberts_fitted) ||
        throw(ArgumentError("PlumeSurfaceConfig.erosion_model must be :regimes or :roberts_fitted, got $(repr(config.erosion_model))"))
    state = PlumeSurfaceState(num_sats; crater_bins=config.crater_bins,
                              crater_min_radius_m=config.crater_min_radius_m,
                              crater_max_radius_m=config.crater_max_radius_m)
    return PlumeSurfaceInteractionModel(config, control, terrain, state)
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

const _PLUME_QUADRATURE_POINTS = 128       # radial nodes of the erosion quadrature

"""
    _plume_peak_shear(cfg, thrust_n, height_m, p0, radius_m) -> Float64

Peak wall shear stress over the footprint. The analytic field peaks at
`r = R_p/sqrt(2)` with the closed-form value `sqrt(2) e^(-1/2) c_f p_0`; a
tabulated field has no closed form, so its profile is scanned on the same grid
the erosion quadrature uses.
"""
@inline function _plume_peak_shear(cfg::PlumeSurfaceConfig{<:PlumeAnalyticField}, thrust_n::Float64,
                                   height_m::Float64, p0::Float64, R::Float64)::Float64
    return cfg.friction_coefficient * p0 * sqrt(2.0) * exp(-0.5)
end

function _plume_peak_shear(cfg::PlumeSurfaceConfig, thrust_n::Float64, height_m::Float64,
                           p0::Float64, R::Float64)::Float64
    peak = 0.0
    dx = _plume_quadrature_limit(cfg.field, height_m, R) / _PLUME_QUADRATURE_POINTS
    @inbounds for k in 1:_PLUME_QUADRATURE_POINTS
        peak = max(peak, _plume_shear_at(cfg.field, cfg, thrust_n, height_m, p0, R, (k - 0.5) * dx))
    end
    return peak
end

"""
    plume_surface_footprint(config, thrust_n, height_m) -> (p0_pa, radius_m)

Peak (stagnation) surface pressure and the footprint radius of the plume of
`thrust_n` newtons standing `height_m` above the ground, read from the
configuration's field. The radius is the one containing `1 - 1/e` of the
integral of the surface pressure over the ground, which for the default
analytic field is exactly `R_p = h tan(θ_p)`: there `p(r) = p0 exp(-(r/R_p)^2)`
normalized so `∫ p dA` is the thrust, so all of the engine's axial momentum is
turned by the surface. The radius never falls below the nozzle exit radius.
"""
@inline function plume_surface_footprint(cfg::PlumeSurfaceConfig, thrust_n::Real, height_m::Real)
    return plume_field_footprint(cfg.field, cfg, thrust_n, height_m)
end

# ---- the erosion regimes on the footprint ---------------------------------------------

"""
    _erosion_environment(config, radius_m) -> ErosionEnvironment

The scenario geometry the erosion regimes need, built from the plume's own
footprint: the footprint radius sets the lateral scale the pore pressure has to
beat for diffusion-driven flow, and the loaded width used by the
bearing-capacity factors is the footprint diameter. The residence time is the
configuration's `gas_residence_time_s`; the viscosity is left `NaN` so the
regimes compute it from the gas temperature.
"""
@inline function _erosion_environment(cfg::PlumeSurfaceConfig, radius_m::Float64)
    R = max(radius_m, cfg.nozzle_exit_radius_m)
    return ErosionEnvironment(footprint_radius_m=R, residence_time_s=cfg.gas_residence_time_s,
                              bearing_width_m=2.0 * R, gas_viscosity_pa_s=NaN)
end

"""
    _local_erosion(config, thrust_n, height_m, radius_m, g, env)

The regimes' total local mass flux and dominant regime at one radius on the
ground, on the gas state the configuration's field reports there. Pure and
allocation-free: `config.regimes` is a tuple of singletons the compiler unrolls.
"""
@inline function _local_erosion(cfg::PlumeSurfaceConfig, F::Float64, h::Float64, r::Float64,
                                g::Float64, env::ErosionEnvironment)
    gas = plume_gas_state(cfg.field, cfg, F, h, r)
    return regolith_erosion_rate(cfg.regimes, gas, cfg.soil, g, env)
end

"""
    _regime_erosion(config, thrust_n, height_m, p0, radius_m, g)

Integrate the regimes' local mass flux over the ground,
`ṁ = ∫ ρ̇(r) 2π r dr`, by the same midpoint quadrature in `x = r/R` the shear
integral uses, so the nodes follow the footprint as the vehicle descends.

Returns `(rate_kg_s, dominant, inner_m, outer_m)`. `dominant` is the regime that
moved the most mass: each quadrature ring is attributed whole to the regime that
dominates *at that radius*, which is exact wherever one regime is alone (the
usual case on the Moon) and is a share, not a decomposition, where two overlap.
`inner_m` and `outer_m` are the edges of the region that moved anything.
"""
function _regime_erosion(cfg::PlumeSurfaceConfig, F::Float64, h::Float64, p0::Float64,
                         R::Float64, g::Float64)
    env = _erosion_environment(cfg, R)
    dx = _plume_quadrature_limit(cfg.field, h, R) / _PLUME_QUADRATURE_POINTS
    total = 0.0
    m_viscous = 0.0
    m_diffusion = 0.0
    m_bearing = 0.0
    inner = Inf
    outer = 0.0
    @inbounds for k in 1:_PLUME_QUADRATURE_POINTS
        x = (k - 0.5) * dx
        r = x * R
        out = _local_erosion(cfg, F, h, r, g, env)
        out.rate_kg_m2_s > 0.0 || continue
        mass = out.rate_kg_m2_s * 2.0 * pi * r * dx * R
        total += mass
        inner = min(inner, r)
        outer = max(outer, r)
        if out.dominant === ViscousErosion
            m_viscous += mass
        elseif out.dominant === DiffusionDrivenFlowRegime
            m_diffusion += mass
        elseif out.dominant === BearingCapacityFailureRegime
            m_bearing += mass
        end
    end
    kind = NoErosion
    best = 0.0
    m_viscous > best && (best = m_viscous; kind = ViscousErosion)
    m_diffusion > best && (best = m_diffusion; kind = DiffusionDrivenFlowRegime)
    m_bearing > best && (best = m_bearing; kind = BearingCapacityFailureRegime)
    return (rate_kg_s=total, dominant=kind, inner_m=(isfinite(inner) ? inner : 0.0), outer_m=outer)
end

"""
    _roberts_erosion(config, thrust_n, height_m, p0, radius_m, v_ej)

The legacy `:roberts_fitted` closure, unchanged: the wall shear stress in excess
of the fitted `threshold_shear_pa`, integrated over the annulus where it is
positive, divided by the ejecta speed and multiplied by the fitted
`erosion_efficiency`. Kept only so the two laws can be compared on one gas
state; see the module header for why it is not the default.
"""
function _roberts_erosion(cfg::PlumeSurfaceConfig, F::Float64, h::Float64, p0::Float64,
                          R::Float64, v_ej::Float64)
    tau_t = cfg.threshold_shear_pa
    dx = _plume_quadrature_limit(cfg.field, h, R) / _PLUME_QUADRATURE_POINTS
    excess_n = 0.0
    inner = Inf
    outer = 0.0
    @inbounds for k in 1:_PLUME_QUADRATURE_POINTS
        x = (k - 0.5) * dx
        tau = _plume_shear_at(cfg.field, cfg, F, h, p0, R, x)
        tau > tau_t || continue
        excess_n += (tau - tau_t) * 2.0 * pi * x * dx * R * R
        inner = min(inner, x * R)
        outer = max(outer, x * R)
    end
    (excess_n > 0.0 && v_ej > 0.0) || return (rate_kg_s=0.0, dominant=NoErosion, inner_m=0.0, outer_m=0.0)
    return (rate_kg_s=cfg.erosion_efficiency * excess_n / v_ej, dominant=ViscousErosion,
            inner_m=(isfinite(inner) ? inner : 0.0), outer_m=outer)
end

"""
    _characteristic_ejecta_speed(config, tau_peak, radius_m) -> Float64

The single characteristic ejecta speed both laws report as `ejecta_mps`: one
grain dragged across one footprint radius from rest by the gas dynamic pressure
at the shear peak, recovered from the peak shear stress through the field's own
shear coefficient, clamped to the configured band. It is unchanged from the
model this file has always carried, and it is the quantity the validation
study's Surveyor III case scores. The full population is
[`plume_ejecta_summary`](@ref).
"""
@inline function _characteristic_ejecta_speed(cfg::PlumeSurfaceConfig, tau_peak::Float64, R::Float64)::Float64
    q_gas = tau_peak / plume_field_shear_coefficient(cfg.field, cfg)
    accel = 3.0 * cfg.particle_drag_coefficient * q_gas / (4.0 * cfg.particle_density_kg_m3 * cfg.particle_diameter_m)
    return clamp(sqrt(max(2.0 * accel * R, 0.0)), cfg.ejecta_speed_min_mps, cfg.ejecta_speed_max_mps)
end

"""
    plume_erosion_onset_height(config, thrust_n[; gravity_m_s2]) -> Float64

Height (m) above which the plume of `thrust_n` newtons moves no soil at all:
the erosion rate is identically zero above it and positive below.

Under the fitted `:roberts_fitted` law this is the height at which the peak wall
shear stress falls to `threshold_shear_pa`, which has a closed form for the
analytic field (about 31 m at the Apollo 11 approach thrust of 11.5 kN, the
height at which the crew first reported blowing dust) and is bisected for a
table.

Under the default `:regimes` law there is no single threshold to invert — three
regimes each have their own criterion, applied to the gas state at every radius
— so the height is bisected on whether the integrated rate is positive. The
bisection runs between the nozzle exit radius and `max_height_m`, and returns
`max_height_m` when the plume is still eroding there.
"""
function plume_erosion_onset_height(cfg::PlumeSurfaceConfig{<:PlumeAnalyticField}, thrust_n::Real;
                                    gravity_m_s2::Real=cfg.gravity_m_s2)::Float64
    F = Float64(thrust_n)
    F > 0.0 || return 0.0
    if cfg.erosion_model === :roberts_fitted
        cfg.threshold_shear_pa > 0.0 || return 0.0
        # τ_peak = 0.8578 c_f F / (π h² tan²θ) = τ_t
        k = sqrt(2.0) * exp(-0.5) * cfg.friction_coefficient * F /
            (pi * tand(cfg.plume_half_angle_deg)^2 * cfg.threshold_shear_pa)
        return sqrt(k)
    end
    return _bisect_onset_height(cfg, F, Float64(gravity_m_s2))
end

function plume_erosion_onset_height(cfg::PlumeSurfaceConfig, thrust_n::Real;
                                    gravity_m_s2::Real=cfg.gravity_m_s2)::Float64
    F = Float64(thrust_n)
    F > 0.0 || return 0.0
    cfg.erosion_model === :roberts_fitted && cfg.threshold_shear_pa <= 0.0 && return 0.0
    return _bisect_onset_height(cfg, F, Float64(gravity_m_s2))
end

"Bisection of the height at which the erosion rate of whichever law is selected first vanishes."
function _bisect_onset_height(cfg::PlumeSurfaceConfig, F::Float64, g::Float64)::Float64
    erodes(h) = plume_quantities(cfg, F, h; gravity_m_s2=g).erosion_kg_s > 0.0
    lo = cfg.nozzle_exit_radius_m
    hi = cfg.max_height_m
    erodes(lo) || return 0.0
    erodes(hi) && return hi
    for _ in 1:60
        mid = 0.5 * (lo + hi)
        if erodes(mid)
            lo = mid
        else
            hi = mid
        end
    end
    return 0.5 * (lo + hi)
end

"""
    plume_quantities(config, thrust_n, height_m; gravity_m_s2=config.gravity_m_s2) -> NamedTuple

The whole plume-surface state at one instant: `pressure_pa` and `shear_pa` are
the peaks of the surface distributions, `erosion_kg_s` is the mass erosion rate
integrated over the ground, `regime` the [`ErosionRegimeKind`](@ref) that moved
the most of it, `ejecta_mps` the characteristic speed the grains leave at,
`inner_m` and `outer_m` the edges of the region that moved anything, and
`ground_effect_n` the thrust augmentation. Everything but the pressure, the
shear stress and the ground effect is zero above the erosion onset height.

With the default `erosion_model = :regimes` the rate is
[`regolith_erosion_rate`](@ref) integrated over the ground on the gas state the
field reports at each radius, so the active regime is whichever one's onset
criterion the gas meets — viscous erosion in Metzger's energy-flux form,
diffusion-driven flow, or bearing-capacity failure. Every coefficient of those
laws is sourced in `regolith_erosion.jl`; this function fits nothing.

With `erosion_model = :roberts_fitted` it is instead the shear-excess momentum
balance,

```math
\\dot m = \\frac{\\eta}{v_{ej}} \\int \\max(\\tau(r) - \\tau_t,\\, 0)\\, \\mathrm{d}A
```

with the two calibrated constants `erosion_efficiency` and
`threshold_shear_pa`, kept so the two can be compared on one gas state.

`gravity_m_s2` is the surface gravity the soil is weighed against; the effector
passes the local `mu/r^2` of the planet it is flying over.
"""
function plume_quantities(cfg::PlumeSurfaceConfig, thrust_n::Real, height_m::Real;
                          gravity_m_s2::Real=cfg.gravity_m_s2)
    F = Float64(thrust_n)
    h = Float64(height_m)
    g = Float64(gravity_m_s2)
    zero_out = (pressure_pa=0.0, shear_pa=0.0, erosion_kg_s=0.0, ejecta_mps=0.0, inner_m=0.0,
                outer_m=0.0, ground_effect_n=0.0, regime=NoErosion)
    (isfinite(F) && F > 0.0 && isfinite(h) && h >= 0.0 && h <= cfg.max_height_m) || return zero_out
    p0, R = plume_surface_footprint(cfg, F, h)
    tau_peak = _plume_peak_shear(cfg, F, h, p0, R)
    ge = plume_ground_effect_force(cfg, F, h)
    v_ej = _characteristic_ejecta_speed(cfg, tau_peak, R)
    er = cfg.erosion_model === :roberts_fitted ?
        _roberts_erosion(cfg, F, h, p0, R, v_ej) :
        _regime_erosion(cfg, F, h, p0, R, g)
    er.rate_kg_s > 0.0 || return (pressure_pa=p0, shear_pa=tau_peak, erosion_kg_s=0.0, ejecta_mps=0.0,
                                  inner_m=0.0, outer_m=0.0, ground_effect_n=ge, regime=NoErosion)
    return (pressure_pa=p0, shear_pa=tau_peak, erosion_kg_s=er.rate_kg_s, ejecta_mps=v_ej,
            inner_m=er.inner_m, outer_m=er.outer_m, ground_effect_n=ge, regime=er.dominant)
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
    cfg = model.config
    (1 <= i <= length(st.height_m)) || return zero3, zero3
    axis = plume_engine_axis(x)
    thrust = _engine_thrust_n(model, i)
    height = _plume_height_along_axis(model, x, env, axis)
    ground_radius = _plume_ground_radius(model, env)
    g = _plume_surface_gravity(env, ground_radius, cfg)
    # The nozzle exit plane above the *current* ground: the plume has already
    # taken some soil out from under it, so the crater's own depth at the
    # stagnation point adds to the height the field is queried at.
    nozzle_height = max(0.0, height - cfg.nozzle_offset_m)
    query_height = cfg.crater_height_feedback ?
        nozzle_height + (isempty(st.crater_profile_m) ? 0.0 : st.crater_profile_m[1, i]) : nozzle_height
    q = plume_quantities(cfg, thrust, query_height; gravity_m_s2=g)
    st.height_m[i] = height
    st.pressure_pa[i] = q.pressure_pa
    st.shear_pa[i] = q.shear_pa
    st.erosion_kg_s[i] = q.erosion_kg_s
    st.ejecta_mps[i] = q.ejecta_mps
    st.ground_effect_n[i] = q.ground_effect_n
    st.regime[i] = plume_regime_code(q.regime)
    st.erosion_radius_m[i] = q.outer_m
    st.last_thrust_n[i] = thrust
    st.last_query_height_m[i] = query_height
    st.last_gravity_m_s2[i] = g
    st.last_body_radius_m[i] = ground_radius
    st.last_inner_m[i] = q.inner_m
    st.last_outer_m[i] = q.outer_m
    _advance_integrals!(model, i, t, q.erosion_kg_s, thrust, query_height, g)
    # The augmentation pushes the vehicle away from the ground, along the
    # thrust direction (the engine axis is the direction the plume travels).
    return -q.ground_effect_n * axis, zero3
end

"Distance from the planet's center to the ground under the vehicle."
@inline function _plume_ground_radius(model::PlumeSurfaceInteractionModel, env::EnvironmentSample)::Float64
    pf = env.planet_frame
    pf === nothing && return Float64(env.planet.Rp_e)
    return _plume_reference_radius(model, env.planet) +
           terrain_height(model.terrain, rad2deg(pf.lat_rad), rad2deg(pf.lon_rad))
end

"""
    _plume_surface_gravity(env, ground_radius_m, config) -> Float64

Local surface gravity the erosion regimes weigh the soil against, `mu/r^2` at
the ground under the vehicle. Falls back to the configuration's
`gravity_m_s2` when the planet carries no gravitational parameter.
"""
@inline function _plume_surface_gravity(env::EnvironmentSample, ground_radius_m::Float64,
                                        cfg::PlumeSurfaceConfig)::Float64
    planet = env.planet
    mu = hasproperty(planet, :μ) ? Float64(getproperty(planet, :μ)) : 0.0
    (mu > 0.0 && ground_radius_m > 0.0) || return cfg.gravity_m_s2
    return mu / (ground_radius_m * ground_radius_m)
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
    _advance_integrals!(model, i, t, rate, thrust_n, height_m, g)

Trapezoidal time integrals of everything the model accumulates: the eroded mass
and the crater's depth profile.

Both are guarded so they only ever advance in time. The solver evaluates the
right-hand side at stage times inside a trial step and repeats them when a step
is rejected; integrating only when the evaluation time exceeds every time
already integrated keeps the integrals monotone and free of the double counting
a naive integral would pick up. That guard is what makes this an accepted-step
advance without an accepted-step callback, and it survives `isolate_state=false`
because the state lives on the effector.
"""
function _advance_integrals!(model::PlumeSurfaceInteractionModel, i::Int, t::Float64,
                             rate::Float64, thrust_n::Float64, height_m::Float64, g::Float64)
    st = model.state
    isfinite(t) || return nothing
    last_t = st.last_time_s[i]
    if !isfinite(last_t)
        st.last_time_s[i] = t
        st.last_rate_kg_s[i] = rate
        _crater_rates!(view(st.crater_rate_kg_m2_s, :, i), model, thrust_n, height_m, g)
        return nothing
    end
    t > last_t || return nothing
    dt = t - last_t
    st.eroded_kg[i] += 0.5 * (st.last_rate_kg_s[i] + rate) * dt
    _advance_crater!(model, i, dt, thrust_n, height_m, g)
    st.last_time_s[i] = t
    st.last_rate_kg_s[i] = rate
    return nothing
end

"""
    _crater_rates!(out, model, thrust_n, height_m, g)

Write the local mass flux the selected erosion law gives at every node of the
crater's radial grid into `out`. The grid is in absolute meters from the
impingement point, unlike the erosion quadrature, which is in footprint radii:
the crater is a fixed patch of ground, the footprint is not.
"""
function _crater_rates!(out, model::PlumeSurfaceInteractionModel, thrust_n::Float64,
                        height_m::Float64, g::Float64)
    cfg = model.config
    radii = model.state.crater_radii_m
    if !(thrust_n > 0.0 && isfinite(height_m) && height_m >= 0.0 && height_m <= cfg.max_height_m)
        @inbounds for k in eachindex(radii)
            out[k] = 0.0
        end
        return nothing
    end
    p0, R = plume_surface_footprint(cfg, thrust_n, height_m)
    if cfg.erosion_model === :roberts_fitted
        # The legacy law is a momentum balance over the whole annulus and has no
        # local rate of its own; its surface flux is the integrand that balance
        # is built from, the excess shear at this radius over the ejecta speed.
        v_ej = _characteristic_ejecta_speed(cfg, _plume_peak_shear(cfg, thrust_n, height_m, p0, R), R)
        tau_t = cfg.threshold_shear_pa
        @inbounds for k in eachindex(radii)
            tau = plume_gas_state(cfg.field, cfg, thrust_n, height_m, radii[k]).shear_pa
            out[k] = (tau > tau_t && v_ej > 0.0) ? cfg.erosion_efficiency * (tau - tau_t) / v_ej : 0.0
        end
        return nothing
    end
    env = _erosion_environment(cfg, R)
    @inbounds for k in eachindex(radii)
        out[k] = _local_erosion(cfg, thrust_n, height_m, radii[k], g, env).rate_kg_m2_s
    end
    return nothing
end

"""
    _advance_crater!(model, i, dt, thrust_n, height_m, g)

Advance the crater's depth profile over `dt` by the trapezoidal rule on the
local mass flux, `Δz = (ρ̇_prev + ρ̇_now) dt / (2 ρ_b)`, then update the derived
depth and radius. The bulk density that converts mass flux to depth is
`bulk_density_kg_m3`: the soil leaves at its in-situ density, so the depth is
the depth of the hole, not the loose volume thrown out of it.
"""
function _advance_crater!(model::PlumeSurfaceInteractionModel, i::Int, dt::Float64,
                          thrust_n::Float64, height_m::Float64, g::Float64)
    st = model.state
    rho_b = model.config.bulk_density_kg_m3
    radii = st.crater_radii_m
    (rho_b > 0.0 && !isempty(radii) && dt > 0.0) || return nothing
    scratch = view(st.crater_scratch_kg_m2_s, :, i)
    _crater_rates!(scratch, model, thrust_n, height_m, g)
    @inbounds for k in eachindex(radii)
        previous = st.crater_rate_kg_m2_s[k, i]
        now = scratch[k]
        (previous > 0.0 || now > 0.0) &&
            (st.crater_profile_m[k, i] += 0.5 * (previous + now) * dt / rho_b)
        st.crater_rate_kg_m2_s[k, i] = now
    end
    _update_crater_summary!(model, i)
    return nothing
end

"Deepest point and edge of the crater, from its depth profile."
function _update_crater_summary!(model::PlumeSurfaceInteractionModel, i::Int)
    st = model.state
    cfg = model.config
    radii = st.crater_radii_m
    isempty(radii) && return nothing
    depth = 0.0
    @inbounds for k in eachindex(radii)
        depth = max(depth, st.crater_profile_m[k, i])
    end
    st.crater_depth_m[i] = depth
    if depth < cfg.crater_min_depth_m
        st.crater_radius_m[i] = 0.0
        return nothing
    end
    # The edge is a tenth of the deepest point, but never shallower than the
    # reporting floor: on a profile that is a fraction of a millimeter deep and
    # nearly flat, a pure ratio test puts the "edge" tens of meters out, wider
    # than the region the plume is eroding at all.
    edge = max(cfg.crater_edge_fraction * depth, cfg.crater_min_depth_m)
    outer = 0.0
    @inbounds for k in eachindex(radii)
        st.crater_profile_m[k, i] >= edge && (outer = radii[k])
    end
    st.crater_radius_m[i] = outer
    return nothing
end

"""
    plume_crater_profile(model, i=1) -> (radius_m, depth_m)

The crater under spacecraft `i` as two vectors: the radial grid in meters from
the impingement point, and the depth eroded at each node. The radii are the
model's own grid (`crater_bins` nodes, logarithmic from `crater_min_radius_m`
to `crater_max_radius_m`) and the depths are a view into the live state, so
reading them during a run is free.

The profile is the quantity a viewer draws and a reference case compares
against; `crater_depth_m` and `crater_radius_m` in [`PlumeSurfaceState`](@ref)
are its deepest point and its edge.
"""
function plume_crater_profile(model::PlumeSurfaceInteractionModel, i::Integer=1)
    st = model.state
    k = Int(i)
    (1 <= k <= size(st.crater_profile_m, 2)) ||
        throw(ArgumentError("spacecraft index $k is outside the plume state"))
    return st.crater_radii_m, view(st.crater_profile_m, :, k)
end

# ---- ejecta transport, once per saved sample -----------------------------------------

"""
    _PlumeEjectaField(config)

Adapter that lets the ejecta transport module read this model's plume field.
`ejecta_gas_state` calls `field(ejecta_config, thrust, height, radius)`, and the
plume field needs the *plume* configuration instead, so this callable carries
it. A struct rather than a closure so the call stays concrete.
"""
struct _PlumeEjectaField{C}
    config::C
end

@inline (f::_PlumeEjectaField)(_ejecta_config, thrust_n, height_m, radius_m) =
    plume_gas_state(f.config.field, f.config, thrust_n, height_m, radius_m)

"""
    _PlumeEjectaWeight(config, gravity, environment, dr)

The mass each launch radius contributes to the ejecta distribution: the
regimes' own local erosion rate there, times the annulus area `2π r Δr`. This
replaces the transport module's default weighting, which is proportional to the
wall shear stress with no threshold and therefore describes where the shear is
rather than where soil actually moves.
"""
struct _PlumeEjectaWeight{C}
    config::C
    gravity_m_s2::Float64
    environment::ErosionEnvironment
    dr_m::Float64
end

@inline function (w::_PlumeEjectaWeight)(gas, radius_m, _diameter_m)
    cfg = w.config
    rate = regolith_erosion_rate(cfg.regimes, gas, cfg.soil, w.gravity_m_s2, w.environment).rate_kg_m2_s
    return rate * 2.0 * pi * Float64(radius_m) * w.dr_m
end

"""
    plume_ejecta_summary(model, i=1) -> Union{Nothing, NamedTuple}

The full grain-transport distribution under spacecraft `i` at the state the
effector last evaluated: what [`ejecta_distribution`](@ref) returns, including
the speed, angle and deposition histograms, or `nothing` when nothing is
eroding.

The launch radii span the region the erosion regimes moved soil in, the mass at
each radius is that region's own local erosion rate, and the mass in each grain
size is [`ejecta_lognormal_mass_weights`](@ref) of the soil's `D50` and
`D84/D50`. Those two weightings are what separate this from a sweep of the shear
profile: without them the micron fines, which carry almost none of the soil's
mass, dominate the deposition radius.

This is expensive by the standards of a right-hand side (one trajectory
integration per size and radius), so it is called once per saved sample by
[`plume_refresh_ejecta!`](@ref) and never from the dynamics.
"""
function plume_ejecta_summary(model::PlumeSurfaceInteractionModel, i::Integer=1)
    cfg = model.config
    st = model.state
    k = Int(i)
    (1 <= k <= length(st.height_m)) || return nothing
    F = st.last_thrust_n[k]
    h = st.last_query_height_m[k]
    (st.erosion_kg_s[k] > 0.0 && F > 0.0) || return nothing
    inner = st.last_inner_m[k]
    outer = st.last_outer_m[k]
    outer > 0.0 || return nothing
    n = max(cfg.ejecta_radii, 2)
    lo = max(inner, 0.25 * (outer - inner) / n, cfg.nozzle_exit_radius_m * 0.1)
    lo < outer || (lo = 0.5 * outer)
    radii = collect(range(lo, outer; length=n))
    dr = (outer - lo) / (n - 1)
    _, R = plume_surface_footprint(cfg, F, h)
    weight = _PlumeEjectaWeight(cfg, st.last_gravity_m_s2[k], _erosion_environment(cfg, R), dr)
    return ejecta_distribution(_PlumeEjectaField(cfg), cfg.ejecta, cfg.soil, F, h;
                               sizes=cfg.ejecta_sizes, radii=radii,
                               size_weights=ejecta_lognormal_mass_weights(cfg.ejecta_sizes, cfg.soil),
                               gravity_m_s2=st.last_gravity_m_s2[k],
                               body_radius_m=st.last_body_radius_m[k], weight=weight)
end

"""
    plume_refresh_ejecta!(model, t, i=1)

Recompute the ejecta summary for spacecraft `i` if it has not already been
computed at time `t`, and publish its three scalars into the state
(`ejecta_angle_deg`, `ejecta_range_m`, `ejecta_escape_frac`). The save fields
call this, so the distribution is evaluated once per saved sample however many
columns read it, and not at all when `ejecta_diagnostic` is off or nothing is
eroding.
"""
function plume_refresh_ejecta!(model::PlumeSurfaceInteractionModel, t::Real, i::Integer=1)
    cfg = model.config
    st = model.state
    k = Int(i)
    (1 <= k <= length(st.height_m)) || return nothing
    cfg.ejecta_diagnostic || return nothing
    time = Float64(t)
    (isfinite(time) && st.ejecta_time_s[k] == time) && return nothing
    st.ejecta_time_s[k] = time
    dist = plume_ejecta_summary(model, k)
    if dist === nothing
        st.ejecta_angle_deg[k] = 0.0
        st.ejecta_range_m[k] = 0.0
        st.ejecta_escape_frac[k] = 0.0
        return nothing
    end
    st.ejecta_angle_deg[k] = dist.mean_angle_deg
    st.ejecta_range_m[k] = dist.mean_deposition_radius_m
    st.ejecta_escape_frac[k] = dist.escape_fraction
    return nothing
end

end # module PlumeSurfaceInteraction
