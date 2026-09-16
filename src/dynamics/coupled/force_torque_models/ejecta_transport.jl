# Transport of the grains a descent engine's plume lifts off the surface: how
# fast the local wall jet can accelerate a grain of a given size, where that
# grain lands, and what the population of grains looks like as a distribution of
# speed, ejection angle and deposition radius.
#
# This module replaces the single characteristic ejecta speed the plume-surface
# effector carries today (`plume_surface_interaction.jl`, `plume_quantities`,
# one grain dragged across one footprint radius at a fixed drag coefficient of
# 2) with a drag balance in a named drag law and a trajectory integration. It is
# pure: it owns no state, mutates nothing, and takes the gas on the ground as a
# NamedTuple with the six fields the plume-field module publishes
# (`pressure_pa`, `shear_pa`, `density_kg_m3`, `speed_mps`, `temperature_k`,
# `mach`), so it can be called from a per-step effector loop.
#
# ---- flow regime ----------------------------------------------------------
#
# The grain is NOT in continuum flow. Worked from the reference field below for
# the Apollo lunar module descent engine at full thrust (45 kN) and 10 m above
# the ground, at the radius where the wall shear stress peaks:
#
#   surface static pressure   p  = 4.0e2 Pa
#   static temperature        T  = 8.5e2 K   (configuration, see below)
#   gas density               rho = p/(R_s T) = 1.2e-3 kg/m^3
#   wall-jet speed            u  = 2.2e3 m/s
#   Chapman-Enskog viscosity  mu = 2.2e-5 Pa s
#   grain diameter            d  = 70e-6 m
#   particle Reynolds number  Re = rho u d / mu    ~ 8
#   particle Mach number      Ma = u / a           ~ 3.4
#   particle Knudsen number   Kn = sqrt(pi gamma/2) Ma/Re ~ 0.6
#
# Kn ~ 0.6 puts the grain in the TRANSITIONAL regime on the Schaaf and Chambre
# classification (continuum Kn <= 0.01, slip 0.01-0.1, transitional 0.1-10, free
# molecular > 10) as tabulated in Capecelatro (2022) Table 2; at the larger
# radii and lower thrusts of the approach the same grain crosses into free
# molecular flow. A continuum drag law is therefore wrong everywhere in the
# problem, and a purely free-molecular one is wrong near the impingement point.
#
# The drag law used is Henderson's, which is a single correlation spanning
# continuum through free-molecular flow and is valid over Re < 2e4 and Ma < 6,
# which brackets every condition above:
#
#   C. B. Henderson, "Drag coefficients of spheres in continuum and rarefied
#   flows", AIAA Journal 14(6), 1976, pp. 707-708.
#
# It is transcribed here from the Matlab listing in Appendix A.2 of
#
#   J. Capecelatro, "Modeling high-speed gas-particle flows relevant to
#   spacecraft landings: a review and perspectives", Journal of Fluid Mechanics
#   (review), arXiv:2109.02523, 2021/2022 (Appendix A.2 for the Henderson
#   listing, Eq. 8 for Kn = sqrt(pi gamma/2) Ma/Re, Table 2 for the regimes).
#
# Henderson's rarefied limit is checked in the tests against the closed-form
# free-molecular drag on a sphere with diffuse reflection (Schaaf and Chambre,
# "Flow of Rarefied Gases", Princeton Aeronautical Paperbacks, 1961): the two
# agree to under one percent at Ma = 6, Re -> 0, which is what "this correlation
# really does reach free-molecular flow" means.
#
# ---- ejection angle -------------------------------------------------------
#
# The ejection angle is NOT computed here; it is an input, and it is the one
# ejecta quantity for which there is a direct Apollo measurement:
#
#   C. D. Immer, J. E. Lane, P. T. Metzger, S. Clements, "Apollo video
#   photogrammetry estimation of plume impingement effects", Icarus 214, 2011,
#   pp. 46-52.
#
# Their Table 1 gives the dust ejection angle above the local horizontal,
# measured at launch from the shadow the lunar module casts on the ejecta sheet:
# Apollo 11 2.6 deg, Apollo 14 2.4 deg, Apollo 15 8.1 deg, Apollo 16 1.4 deg,
# Apollo 17 2.0 deg, average 3.3 deg (Apollo 12 could not be measured). They
# state the angle is 1-3 deg for the majority of the landing and attribute the
# Apollo 15 outlier to an 11 deg surface slope. They also note (their section 3) that on Roberts'
# theory the angle is set by the wall angle of the scour crater, so it is only
# weakly coupled to thrust -- which is why treating it as a configured range
# rather than a computed quantity is defensible. Metzger (2023), below, makes
# the same assumption and takes the angle uniformly distributed over 1-3 deg.
#
# ---- what the speeds should come out near ---------------------------------
#
#   P. T. Metzger, "The damage to lunar orbiting spacecraft caused by the ejecta
#   of lunar landers", Proc. ASCE Earth and Space 2022/2023,
#   arXiv:2305.12234, section "Ejecta trajectories":
#   the maximum ejecta velocity for the finest dust equals the exhaust velocity
#   of the propellant, about 3,100 m/s for the Apollo lunar module
#   (Aerozine-50/N2O4); particles of a given size are distributed between 30 and
#   100 percent of that maximum; smaller particles go faster because the ratio
#   of drag force to inertia scales as 1/d.
#
#   Immer et al. (2011) section 1, summarizing the Surveyor 3 pitting analyses:
#   40 m/s (Nickle and Carroll 1972), > 70 m/s (Jaffe 1972), 100 m/s
#   (Cour-Palais et al. 1972), and the estimate they call the most reliable,
#   300 to 2000 m/s from the surface structure of the pits (Brownlee, Bucher et
#   al. 1972). They also quote the lunar escape velocity as about 2373 m/s.
#
# These are the numbers the distribution's mean and maximum speed are tested
# against; they are bounds on the fast tail that hit Surveyor 3, not a measured
# mean of the whole population, and the test tolerances are correspondingly wide.
#
# ---- soil -----------------------------------------------------------------
#
# Lunar mare values are from the Lunar Sourcebook (G. Heiken, D. Vaniman,
# B. M. French, eds., "Lunar Sourcebook: A User's Guide to the Moon", Cambridge
# University Press, 1991), chapter 9 "Physical Properties of the Lunar Surface"
# by W. D. Carrier III, G. R. Olhoeft and W. Mendell:
#   * median particle size 40 to 130 um, average 70 um (section 9.1.1);
#   * specific gravity 2.3 to > 3.2, "we recommend a value of 3.1 for general
#     scientific and engineering analyses" (section 9.1.3, Table 9.3), so a
#     particle density of 3100 kg/m^3;
#   * in situ bulk density about 1.5 g/cm^3 in the top 15 cm (section 9.1.4,
#     Tables 9.4 and 9.5), so 1500 kg/m^3.
#
# ---- what is assumed, not sourced -----------------------------------------
#
# Every one of these is a named field of `EjectaTransportConfig` or
# `EjectaReferenceGasField` with its reasoning on the field, and each is listed
# again in `docs/src/user/lunar_landing.md`:
#   * `entrainment_length_factor` / `entrainment_length_min_m` -- the distance
#     over which the wall jet accelerates a grain before it leaves the surface;
#   * `wall_jet_growth_rate` / `wall_jet_thickness_min_m` -- how fast the drag
#     decays with height above the surface;
#   * `molecular_collision_diameter_m` -- the hard-sphere diameter that sets the
#     Chapman-Enskog viscosity;
#   * `grain_temperature_k` -- enters the drag law only through T_p/T_gas;
#   * `gas_molar_mass_kg_mol`, `gas_gamma` -- exhaust composition, overridden by
#     whatever the plume field supplies;
#   * `EjectaReferenceGasField.exit_static_temperature_k` and
#     `wall_jet_decay_radii` -- the reference field is a stand-in for testing and
#     for use before a real plume field is available, not a plume model.
module EjectaTransport

using SpecialFunctions: erf

export EjectaSoil, EjectaTransportConfig, EjectaReferenceGasField
export ejecta_gas_state, ejecta_vacuum_gas_state, ejecta_gas_viscosity
export ejecta_particle_flow_numbers, ejecta_flow_regime, ejecta_regime_name
export ejecta_drag_coefficient, ejecta_free_molecular_drag_coefficient
export ejecta_launch_speed, ejecta_trajectory, ejecta_distribution, ejecta_escape_speed

# SI defining constants (CODATA/SI 2019 exact values).
const EJECTA_BOLTZMANN_J_PER_K = 1.380649e-23
const EJECTA_AVOGADRO_PER_MOL = 6.02214076e23
const EJECTA_GAS_CONSTANT_J_PER_MOL_K = EJECTA_BOLTZMANN_J_PER_K * EJECTA_AVOGADRO_PER_MOL

"""Flow-regime codes returned by [`ejecta_flow_regime`](@ref), after the
Knudsen-number classification of Schaaf and Chambre (1958) as tabulated in
Capecelatro (2022), Table 2."""
const EJECTA_REGIME_CONTINUUM = 0x01
const EJECTA_REGIME_SLIP = 0x02
const EJECTA_REGIME_TRANSITIONAL = 0x03
const EJECTA_REGIME_FREE_MOLECULAR = 0x04

"""
    ejecta_regime_name(code) -> Symbol

Name of a flow-regime code: `:continuum`, `:slip`, `:transitional`,
`:free_molecular`, or `:unknown`.
"""
function ejecta_regime_name(code::Integer)::Symbol
    c = UInt8(code)
    c == EJECTA_REGIME_CONTINUUM && return :continuum
    c == EJECTA_REGIME_SLIP && return :slip
    c == EJECTA_REGIME_TRANSITIONAL && return :transitional
    c == EJECTA_REGIME_FREE_MOLECULAR && return :free_molecular
    return :unknown
end

# ---- soil and configuration ----------------------------------------------------------

"""
    EjectaSoil(; particle_density_kg_m3=3100.0, bulk_density_kg_m3=1500.0,
               median_particle_diameter_m=70.0e-6)

The three soil properties grain transport needs. Defaults are lunar mare values
from the Lunar Sourcebook chapter 9 (see the module header): a particle density
of 3100 kg/m³ (the recommended specific gravity of 3.1, section 9.1.3), an in
situ bulk density of 1500 kg/m³ in the top 15 cm (section 9.1.4) and a median
particle diameter of 70 µm (section 9.1.1, whose range across the returned soils
is 40 to 130 µm).

Only `particle_density_kg_m3` enters the drag balance; the other two are carried
so a call site records the soil it flew over and so the distribution can be
weighted by a size distribution. Every function here reaches the soil through
property access, so the regolith module's `RegolithProperties` -- or any other
object exposing `particle_density_kg_m3` -- can be passed instead of this type.
"""
Base.@kwdef struct EjectaSoil
    particle_density_kg_m3::Float64 = 3_100.0
    bulk_density_kg_m3::Float64 = 1_500.0
    median_particle_diameter_m::Float64 = 70.0e-6
end

"""
    EjectaTransportConfig(; ...)

Everything the transport model needs that is not the gas state, the soil or the
grain size. Fields marked ASSUMPTION below have no measurement behind them; they
are the model's free parameters and are reported as such.

Exhaust gas composition (used for the viscosity, the Knudsen number and, when
the gas state carries no Mach number, the speed of sound):

- `gas_molar_mass_kg_mol` (0.0215) and `gas_gamma` (1.24) are representative of
  N₂O₄/Aerozine-50 combustion products. ASSUMPTION: they are order-of-magnitude
  right for a storable hypergolic exhaust but are not taken from an equilibrium
  chemistry calculation. A plume field that carries its own composition should
  override them.
- `molecular_collision_diameter_m` (4.0e-10) is the hard-sphere diameter in the
  Chapman-Enskog viscosity. ASSUMPTION: 4 Å is between the kinetic diameters of
  N₂ (3.6 Å) and CO₂ (4.5 Å), the bulk of a hypergolic exhaust. The viscosity
  enters only through the Reynolds number, and the drag coefficient is a weak
  function of Re in the rarefied branch, so a 30 percent error here moves the
  launch speed by a few percent.
- `grain_temperature_k` (250.0) is the grain's own temperature, which enters
  Henderson's correlation only as the ratio T_p/T_gas. ASSUMPTION: lunar surface
  soil at the low sun elevations (about 11°) of the Apollo landings.

Drag law:

- `drag_model` is `:henderson` (the default; Henderson 1976, valid Re < 2e4,
  Ma < 6, continuum through free molecular) or `:constant`, which uses
  `constant_drag_coefficient` unchanged. `:constant` exists so the launch
  integration can be checked against its closed-form solution; it is not a
  physical model.
- `constant_drag_coefficient` (2.0) is the hypersonic free-molecular limit for a
  sphere, and is what `plume_surface_interaction.jl` uses today.

Launch (see [`ejecta_launch_speed`](@ref)):

- `launch_steps` (256) is the fixed number of RK4 steps of the launch
  integration.
- `entrainment_length_factor` (1.0) and `entrainment_length_min_m` (0.5) set the
  default entrainment length used by [`ejecta_distribution`](@ref) for a grain
  starting at radius `r`: `max(factor * r, min_m)`. ASSUMPTION, and the single
  most influential free parameter of the launch model. The reasoning for scaling
  it with `r` is that the wall jet's own properties vary on the scale of the
  radius itself (the surface pressure and shear fall like exp(-(r/R)²) with
  R ∝ h), so a grain is inside gas of roughly its launch condition for a
  distance of order `r`; the floor keeps the innermost radii from being frozen
  at zero speed. `ejecta_launch_speed` takes the length explicitly, so a caller
  with a better estimate is not bound by this.

Ejection angle (see the module header for the source):

- `ejection_angle_min_deg` (1.0) and `ejection_angle_max_deg` (3.0) bound the
  angle above the local horizontal, measured at launch, at which grains leave.
  SOURCED: Immer, Lane, Metzger and Clements (Icarus 214, 2011; Earth and Space
  2008) measure 1-3 degrees from the Apollo landing films (their Table 1: 2.6,
  2.4, 8.1, 1.4 and 2.0 degrees for Apollo 11, 14, 15, 16 and 17, mean 3.3, the
  8.1 attributed to an 11 degree surface slope). The distribution over the range
  is taken uniform, following Metzger (2023). `ejecta_distribution` reports the
  mass-weighted mean as `mean_angle_deg` and the spread as `angle_edges_deg` /
  `angle_fraction`, in the same units and convention.

Flight (see [`ejecta_trajectory`](@ref)):

- `trajectory_steps` (512) is the nominal step count over the drag-free flight
  time; `trajectory_max_steps` (8192) caps a flight that drag lengthens.
- `wall_jet_growth_rate` (0.1) and `wall_jet_thickness_min_m` (0.05) set the
  height `delta = max(growth * r, min_m)` over which the gas density seen by a
  grain in flight decays, as `exp(-z/delta)`. SOURCED in form only: the outer
  layer of a turbulent radial wall jet grows linearly with radius (M. Poreh,
  Y. G. Tsuei, J. E. Cermak, "Investigation of a turbulent radial wall jet",
  Journal of Applied Mechanics 34(2), 1967, pp. 457-463). ASSUMPTION: the
  growth rate of 0.1 and the exponential (rather than measured) profile shape.
  This is what "the plume drag decaying away from the impingement point" means
  here: the drag falls both because the field's own gas state decays with radius
  and because the grain climbs out of the jet.
- `max_flight_radius_m` (1.0e5) stops a grain that is neither landing nor
  escaping.
"""
Base.@kwdef struct EjectaTransportConfig
    gas_molar_mass_kg_mol::Float64 = 0.0215
    gas_gamma::Float64 = 1.24
    molecular_collision_diameter_m::Float64 = 4.0e-10
    grain_temperature_k::Float64 = 250.0
    drag_model::Symbol = :henderson
    constant_drag_coefficient::Float64 = 2.0
    launch_steps::Int = 256
    entrainment_length_factor::Float64 = 1.0
    entrainment_length_min_m::Float64 = 0.5
    ejection_angle_min_deg::Float64 = 1.0
    ejection_angle_max_deg::Float64 = 3.0
    trajectory_steps::Int = 512
    trajectory_max_steps::Int = 8_192
    wall_jet_growth_rate::Float64 = 0.1
    wall_jet_thickness_min_m::Float64 = 0.05
    max_flight_radius_m::Float64 = 1.0e5
end

# ---- gas state -----------------------------------------------------------------------

"""
    ejecta_gas_state(field, config, thrust_n, height_m, radius_m) -> NamedTuple

The gas on the ground at `radius_m` from the stagnation point, under an engine
of `thrust_n` newtons whose exit plane is `height_m` above the surface, as the
six-field NamedTuple the plume-field module publishes: `pressure_pa`,
`shear_pa`, `density_kg_m3`, `speed_mps`, `temperature_k`, `mach`.

The generic method forwards to `field(config, thrust_n, height_m, radius_m)`, so
`field` may be the plume-field module's query closure
`(cfg, F, h, r) -> plume_gas_state(plume_field, cfg, F, h, r)` -- this module
never imports the field module, it only calls what it is handed. A method for
[`EjectaReferenceGasField`](@ref) is provided so the transport model is usable
and testable on its own.
"""
@inline function ejecta_gas_state(field::F, config, thrust_n::Real, height_m::Real,
                                  radius_m::Real) where {F}
    return field(config, Float64(thrust_n), Float64(height_m), Float64(radius_m))
end

"""
    ejecta_vacuum_gas_state() -> NamedTuple

A gas state with nothing in it. A grain flown through this feels only gravity,
which is how the ballistic limit of [`ejecta_trajectory`](@ref) is tested.
"""
@inline ejecta_vacuum_gas_state() = (pressure_pa=0.0, shear_pa=0.0, density_kg_m3=0.0,
                                     speed_mps=0.0, temperature_k=1.0, mach=0.0)

"""
    EjectaReferenceGasField(; ...)

A minimal surface gas field, for testing this module and for running it before a
real plume field is available. It is a STAND-IN, not a plume model: the plume
field is owned by `plume_gas_field.jl`, and a caller with one should pass that
instead.

Its pressure and shear stress are the closure `plume_surface_interaction.jl`
already uses -- a Gaussian footprint of radius `R = max(h tan(θ_p), r_exit)`
normalized so `∫p dA` is the thrust, and `τ = c_f p₀ (2r/R) exp(-(r/R)²)` -- which
follows L. Roberts, "The action of a hypersonic jet on a dust layer", IAS Paper
63-50, 1963, in the form used by Metzger and co-workers. On top of that:

- `exhaust_velocity_mps` (3100.0) is the speed the wall jet carries at the
  stagnation point. SOURCED: Metzger (2023, arXiv:2305.12234) gives about
  3,100 m/s as the Apollo lunar module exhaust velocity; the descent engine's
  tabulated vacuum specific impulse of 311 s is an independent 3,050 m/s.
- `wall_jet_decay_radii` (2.0) is the e-folding length of that speed in
  footprint radii. ASSUMPTION.
- `exit_static_temperature_k` (850.0) is the static temperature of the gas.
  ASSUMPTION: for a γ ≈ 1.24 exhaust leaving at 3,100 m/s, energy conservation
  from a hypergolic flame temperature of order 3,000 K leaves several hundred
  kelvin of static temperature, and the isentropic relation at the descent
  engine's exit Mach number of about 4.8 gives the same order.
- `friction_coefficient` (0.01), `plume_half_angle_deg` (25.0) and
  `nozzle_exit_radius_m` (0.75) are the values `PlumeSurfaceConfig` carries for
  the Apollo descent engine.

Density follows from the ideal gas law at the local static pressure and this
temperature, and the Mach number from the speed and `sqrt(γ R_s T)`.
"""
Base.@kwdef struct EjectaReferenceGasField
    plume_half_angle_deg::Float64 = 25.0
    nozzle_exit_radius_m::Float64 = 0.75
    friction_coefficient::Float64 = 0.01
    exhaust_velocity_mps::Float64 = 3_100.0
    wall_jet_decay_radii::Float64 = 2.0
    exit_static_temperature_k::Float64 = 850.0
end

@inline function ejecta_gas_state(field::EjectaReferenceGasField, config::EjectaTransportConfig,
                                  thrust_n::Real, height_m::Real, radius_m::Real)
    F = Float64(thrust_n)
    h = Float64(height_m)
    r = Float64(radius_m)
    (isfinite(F) && F > 0.0 && isfinite(h) && h >= 0.0 && isfinite(r) && r >= 0.0) ||
        return ejecta_vacuum_gas_state()
    R = max(h * tand(field.plume_half_angle_deg), field.nozzle_exit_radius_m)
    p0 = F / (pi * R * R)
    x = r / R
    e = exp(-x * x)
    p = p0 * e
    tau = field.friction_coefficient * p0 * 2.0 * x * e
    u = field.exhaust_velocity_mps * exp(-x / max(field.wall_jet_decay_radii, 1.0e-12))
    T = field.exit_static_temperature_k
    Rs = EJECTA_GAS_CONSTANT_J_PER_MOL_K / config.gas_molar_mass_kg_mol
    rho = p / (Rs * T)
    a = sqrt(config.gas_gamma * Rs * T)
    return (pressure_pa=p, shear_pa=tau, density_kg_m3=rho, speed_mps=u, temperature_k=T,
            mach=u / a)
end

# ---- gas transport properties --------------------------------------------------------

"""
    ejecta_gas_viscosity(config, temperature_k) -> Float64

Dynamic viscosity (Pa s) of the exhaust gas from Chapman-Enskog hard-sphere
kinetic theory,

```math
\\mu = \\frac{5}{16}\\frac{\\sqrt{\\pi m k_B T}}{\\pi \\sigma^2},
```

with `m` the molecular mass from `gas_molar_mass_kg_mol` and `σ` the
`molecular_collision_diameter_m` (S. Chapman and T. G. Cowling, "The
Mathematical Theory of Non-Uniform Gases", 3rd ed., Cambridge, 1970, ch. 10).
A hard-sphere model is used rather than a Sutherland fit because no Sutherland
constants exist for this exhaust mixture; it is the honest first-principles
option, and the collision diameter it needs is a documented assumption.
"""
@inline function ejecta_gas_viscosity(config::EjectaTransportConfig, temperature_k::Real)::Float64
    T = Float64(temperature_k)
    T > 0.0 || return 0.0
    m = config.gas_molar_mass_kg_mol / EJECTA_AVOGADRO_PER_MOL
    sigma = config.molecular_collision_diameter_m
    return 0.3125 * sqrt(pi * m * EJECTA_BOLTZMANN_J_PER_K * T) / (pi * sigma * sigma)
end

@inline function _ejecta_sound_speed(config::EjectaTransportConfig, gas)::Float64
    u = Float64(gas.speed_mps)
    M = Float64(gas.mach)
    (M > 0.0 && u > 0.0 && isfinite(M)) && return u / M
    T = Float64(gas.temperature_k)
    T > 0.0 || return 0.0
    Rs = EJECTA_GAS_CONSTANT_J_PER_MOL_K / config.gas_molar_mass_kg_mol
    return sqrt(config.gas_gamma * Rs * T)
end

"""
    ejecta_particle_flow_numbers(config, gas, relative_speed_mps, diameter_m)
        -> (reynolds, mach, knudsen)

The three dimensionless numbers the drag law needs for a grain of `diameter_m`
moving at `relative_speed_mps` through `gas`: the particle Reynolds number
`ρ u d / μ`, the particle Mach number `u / a`, and the particle Knudsen number.

The Knudsen number uses the ideal-gas identity `Kn = sqrt(πγ/2) Ma / Re`
(Clift, Grace and Weber, "Bubbles, Drops, and Particles", as quoted in
Capecelatro 2022 Eq. 8), which keeps it consistent with the Reynolds number the
drag correlation is evaluated at instead of introducing a second, independent
mean-free-path model.
"""
@inline function ejecta_particle_flow_numbers(config::EjectaTransportConfig, gas,
                                              relative_speed_mps::Real, diameter_m::Real)
    w = Float64(relative_speed_mps)
    d = Float64(diameter_m)
    rho = Float64(gas.density_kg_m3)
    (w > 0.0 && d > 0.0 && rho > 0.0) || return (0.0, 0.0, Inf)
    mu = ejecta_gas_viscosity(config, gas.temperature_k)
    mu > 0.0 || return (0.0, 0.0, Inf)
    re = rho * w * d / mu
    a = _ejecta_sound_speed(config, gas)
    ma = a > 0.0 ? w / a : 0.0
    kn = re > 0.0 ? sqrt(pi * config.gas_gamma / 2.0) * ma / re : Inf
    return (re, ma, kn)
end

"""
    ejecta_flow_regime(knudsen) -> UInt8

Flow regime of a grain at a given Knudsen number, on the classification of
Schaaf and Chambre (1958) as tabulated in Capecelatro (2022) Table 2:
continuum for `Kn <= 0.01`, slip to `0.1`, transitional to `10`, free molecular
above. Returns one of `EJECTA_REGIME_CONTINUUM`, `EJECTA_REGIME_SLIP`,
`EJECTA_REGIME_TRANSITIONAL`, `EJECTA_REGIME_FREE_MOLECULAR`; name it with
`ejecta_regime_name`.
"""
@inline function ejecta_flow_regime(knudsen::Real)::UInt8
    kn = Float64(knudsen)
    isfinite(kn) || return EJECTA_REGIME_FREE_MOLECULAR
    kn <= 0.01 && return EJECTA_REGIME_CONTINUUM
    kn <= 0.1 && return EJECTA_REGIME_SLIP
    kn <= 10.0 && return EJECTA_REGIME_TRANSITIONAL
    return EJECTA_REGIME_FREE_MOLECULAR
end

# ---- drag law ------------------------------------------------------------------------

@inline function _henderson_subsonic(re::Float64, ma::Float64, gamma::Float64, tr::Float64)::Float64
    s = ma * sqrt(0.5 * gamma)
    sqrt_re = sqrt(re)
    denom_a = re + ma * sqrt(0.5 * gamma) *
              (4.33 + (3.65 - 1.53 * tr) / (1.0 + 0.353 * tr)) *
              exp(-0.247 * re / max(s, 1.0e-300))
    term_a = 24.0 / max(denom_a, 1.0e-300)
    poly = 0.03 * re + 0.48 * sqrt_re
    term_b = exp(-0.5 * ma / max(sqrt_re, 1.0e-300)) *
             ((4.5 + 0.38 * poly) / (1.0 + poly) + 0.1 * ma^2 + 0.2 * ma^8)
    term_c = (1.0 - exp(-ma / max(re, 1.0e-300))) * 0.6 * s
    return term_a + term_b + term_c
end

@inline function _henderson_supersonic(re::Float64, ma::Float64, gamma::Float64, tr::Float64)::Float64
    s = ma * sqrt(0.5 * gamma)
    s2 = s * s
    root = 1.86 * sqrt(ma / max(re, 1.0e-300))
    bracket = 2.0 + 2.0 / s2 + (1.058 / s) * sqrt(tr) - 1.0 / (s2 * s2)
    return (0.9 + 0.34 / (ma * ma) + root * bracket) / (1.0 + root)
end

"""
    ejecta_drag_coefficient(config, gas, relative_speed_mps, diameter_m) -> Float64

Drag coefficient of a spherical grain of `diameter_m` moving at
`relative_speed_mps` through `gas`.

With `config.drag_model === :constant` this is `constant_drag_coefficient`, for
checking the integrators against closed forms. Otherwise it is Henderson's
correlation (C. B. Henderson, "Drag coefficients of spheres in continuum and
rarefied flows", AIAA Journal 14(6), 1976, 707-708), transcribed from the
listing in Capecelatro (2022) Appendix A.2: the subsonic branch below Ma = 1,
the supersonic branch above Ma = 1.75, and Henderson's linear interpolation
`C_D(1) + (4/3)(Ma - 1)[C_D(1.75) - C_D(1)]` between them.

Henderson's correlation distinguishes the flow relative to the particle from the
undisturbed freestream, through `Re_inf` and `Ma_inf`. A grain sitting in a wall
jet has no separate freestream, so `Re_inf = Re_p` and `Ma_inf = Ma_p` here.
That is a modeling choice, not part of the correlation.

The correlation is validated over `Re < 2e4` and `Ma < 6`, which covers the wall
jet; outside that range it is still evaluated, and the caller is responsible for
knowing it is extrapolating.
"""
@inline function ejecta_drag_coefficient(config::EjectaTransportConfig, gas,
                                         relative_speed_mps::Real, diameter_m::Real)::Float64
    config.drag_model === :constant && return config.constant_drag_coefficient
    re, ma, _ = ejecta_particle_flow_numbers(config, gas, relative_speed_mps, diameter_m)
    (re > 0.0 && ma > 0.0) || return 0.0
    g = config.gas_gamma
    tr = config.grain_temperature_k / max(Float64(gas.temperature_k), 1.0e-12)
    ma <= 1.0 && return _henderson_subsonic(re, ma, g, tr)
    ma >= 1.75 && return _henderson_supersonic(re, ma, g, tr)
    cd1 = _henderson_subsonic(re, 1.0, g, tr)
    cd175 = _henderson_supersonic(re, 1.75, g, tr)
    return cd1 + (4.0 / 3.0) * (ma - 1.0) * (cd175 - cd1)
end

"""
    ejecta_free_molecular_drag_coefficient(speed_ratio, temperature_ratio) -> Float64

Closed-form drag coefficient of a sphere in free molecular flow with diffuse
reflection and complete accommodation,

```math
C_D = \\frac{2S^2+1}{\\sqrt{\\pi}S^3}e^{-S^2}
    + \\frac{4S^4+4S^2-1}{2S^4}\\mathrm{erf}(S)
    + \\frac{2\\sqrt{\\pi}}{3S}\\sqrt{T_w/T_\\infty},
```

with `S = Ma sqrt(γ/2)` the molecular speed ratio (S. A. Schaaf and
P. L. Chambre, "Flow of Rarefied Gases", Princeton Aeronautical Paperbacks,
Princeton University Press, 1961). This is not used by the transport model; it
is the asymptote [`ejecta_drag_coefficient`](@ref) is checked against as
`Re -> 0`, which is how the claim that Henderson's correlation really reaches
free molecular flow is tested rather than asserted.
"""
@inline function ejecta_free_molecular_drag_coefficient(speed_ratio::Real, temperature_ratio::Real)::Float64
    s = Float64(speed_ratio)
    s > 0.0 || return Inf
    s2 = s * s
    s4 = s2 * s2
    return (2.0 * s2 + 1.0) / (sqrt(pi) * s2 * s) * exp(-s2) +
           (4.0 * s4 + 4.0 * s2 - 1.0) / (2.0 * s4) * erf(s) +
           2.0 * sqrt(pi) / (3.0 * s) * sqrt(Float64(temperature_ratio))
end

# ---- launch --------------------------------------------------------------------------

"""
    ejecta_launch_speed(gas, soil, diameter_m, entrainment_length_m; config) -> NamedTuple

Speed a grain of `diameter_m` reaches while the local wall jet drags it over
`entrainment_length_m` of surface, starting from rest.

The balance is the drag on a sphere against its own inertia,

```math
\\frac{\\mathrm{d}v}{\\mathrm{d}t}
 = \\frac{3\\rho_g C_D(|u-v|)}{4\\rho_p d}\\,(u-v)\\,|u-v|,
\\qquad \\frac{\\mathrm{d}x}{\\mathrm{d}t} = v,
```

integrated with fixed-step RK4 until `x` reaches `entrainment_length_m`; the
partial last step is resolved on the local quadratic in `x`. `C_D` is
[`ejecta_drag_coefficient`](@ref), so at the densities under a descent engine it
is Henderson's transitional/free-molecular correlation and not a continuum law
(see the module header for the worked Knudsen number). The grain never overtakes
the gas: `u` is the asymptote of this equation, which is why the fastest and
finest ejecta approach the exhaust velocity itself.

Gravity and cohesion do not appear. This function answers "how fast does a grain
that is already mobile get going"; whether it is mobile at all is the erosion
model's question (`regolith_erosion.jl`), and the caller is expected to have
asked it first.

Returns `(speed_mps, drag_coefficient, reynolds, mach, knudsen, regime,
distance_m, reached_length)`, where `regime` is a
[`ejecta_flow_regime`](@ref) code for the conditions at launch and
`reached_length` says whether the integration actually covered the entrainment
length before running out of steps. Allocation-free.
"""
function ejecta_launch_speed(gas, soil, diameter_m::Real, entrainment_length_m::Real;
                             config::EjectaTransportConfig=EjectaTransportConfig())
    d = Float64(diameter_m)
    len = Float64(entrainment_length_m)
    u = Float64(gas.speed_mps)
    rho = Float64(gas.density_kg_m3)
    rho_p = Float64(soil.particle_density_kg_m3)
    re0, ma0, kn0 = ejecta_particle_flow_numbers(config, gas, u, d)
    regime = ejecta_flow_regime(kn0)
    dead = (speed_mps=0.0, drag_coefficient=0.0, reynolds=re0, mach=ma0, knudsen=kn0,
            regime=regime, distance_m=0.0, reached_length=false)
    (d > 0.0 && len > 0.0 && u > 0.0 && rho > 0.0 && rho_p > 0.0) || return dead
    cd0 = ejecta_drag_coefficient(config, gas, u, d)
    cd0 > 0.0 || return dead
    # k has units 1/m: the reciprocal of the distance over which a grain would
    # reach the gas speed if the drag stayed at its initial value.
    k0 = 3.0 * rho * cd0 / (4.0 * rho_p * d)
    t_span = 4.0 * (1.0 / (k0 * u) + len / u)
    steps = max(config.launch_steps, 8)
    dt = t_span / steps
    x = 0.0
    v = 0.0
    reached = false
    @inbounds for _ in 1:steps
        x_prev = x
        v_prev = v
        a_prev = _launch_accel(config, gas, soil, d, v)
        k1x, k1v = v, a_prev
        k2x, k2v = v + 0.5 * dt * k1v, _launch_accel(config, gas, soil, d, v + 0.5 * dt * k1v)
        k3x, k3v = v + 0.5 * dt * k2v, _launch_accel(config, gas, soil, d, v + 0.5 * dt * k2v)
        k4x, k4v = v + dt * k3v, _launch_accel(config, gas, soil, d, v + dt * k3v)
        x += (dt / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x)
        v += (dt / 6.0) * (k1v + 2.0 * k2v + 2.0 * k3v + k4v)
        if x >= len
            # Cubic Hermite interpolation across the step just taken -- the
            # positions and speeds at both ends are known, so the crossing is
            # resolved to the order of the integrator itself rather than to the
            # order of a linear back-off.
            a_end = _launch_accel(config, gas, soil, d, v)
            theta = _hermite_crossing(x_prev, v_prev * dt, x, v * dt, len)
            v = _hermite(v_prev, a_prev * dt, v, a_end * dt, theta)
            x = len
            reached = true
            break
        end
    end
    v = min(v, u)
    cd = ejecta_drag_coefficient(config, gas, max(u - v, 0.0), d)
    re, ma, kn = ejecta_particle_flow_numbers(config, gas, max(u - v, 0.0), d)
    return (speed_mps=v, drag_coefficient=cd, reynolds=re, mach=ma, knudsen=kn,
            regime=ejecta_flow_regime(kn), distance_m=x, reached_length=reached)
end

@inline function _launch_accel(config::EjectaTransportConfig, gas, soil, d::Float64, v::Float64)::Float64
    u = Float64(gas.speed_mps)
    w = u - v
    w > 0.0 || return 0.0
    cd = ejecta_drag_coefficient(config, gas, w, d)
    cd > 0.0 || return 0.0
    return 3.0 * Float64(gas.density_kg_m3) * cd * w * w / (4.0 * Float64(soil.particle_density_kg_m3) * d)
end

"""
    _hermite(p0, m0, p1, m1, theta) -> Float64

Cubic Hermite interpolation on a unit interval, with `m0`/`m1` the derivatives
already scaled by the step length.
"""
@inline function _hermite(p0::Float64, m0::Float64, p1::Float64, m1::Float64, theta::Float64)::Float64
    t2 = theta * theta
    t3 = t2 * theta
    return (2.0t3 - 3.0t2 + 1.0) * p0 + (t3 - 2.0t2 + theta) * m0 +
           (-2.0t3 + 3.0t2) * p1 + (t3 - t2) * m1
end

"""
    _hermite_crossing(p0, m0, p1, m1, target) -> Float64

Position in `[0, 1]` where the cubic Hermite through `(p0, m0)` and `(p1, m1)`
reaches `target`, by bisection. The interpolant is monotone over a step of an
integration whose speed never changes sign, which is every use here.
"""
@inline function _hermite_crossing(p0::Float64, m0::Float64, p1::Float64, m1::Float64,
                                   target::Float64)::Float64
    p1 > p0 || return 1.0
    lo = 0.0
    hi = 1.0
    @inbounds for _ in 1:60
        mid = 0.5 * (lo + hi)
        if _hermite(p0, m0, p1, m1, mid) < target
            lo = mid
        else
            hi = mid
        end
    end
    return 0.5 * (lo + hi)
end

"""
    _quadratic_crossing(f0, f1, f2, dt) -> Float64

Smallest root in `(0, dt]` of `f0 + f1 s + f2 s²/2`, falling back to linear
interpolation when the quadratic term is degenerate. Used to land the last step
of an integration exactly on a crossing (the entrainment length, or the ground).
"""
@inline function _quadratic_crossing(f0::Float64, f1::Float64, f2::Float64, dt::Float64)::Float64
    if abs(f2) > 1.0e-30
        disc = f1 * f1 - 2.0 * f2 * f0
        if disc >= 0.0
            root = sqrt(disc)
            s1 = (-f1 + root) / f2
            s2 = (-f1 - root) / f2
            lo, hi = minmax(s1, s2)
            lo > 0.0 && lo <= dt && return lo
            hi > 0.0 && hi <= dt && return hi
        end
    end
    abs(f1) > 1.0e-30 && return clamp(-f0 / f1, 0.0, dt)
    return dt
end

# ---- flight --------------------------------------------------------------------------

"""
    ejecta_trajectory(gas_source, soil, diameter_m, launch_speed_mps, angle_rad,
                      gravity_m_s2; config, start_radius_m, escape_speed_mps) -> NamedTuple

Ballistic flight of one grain, launched at `launch_speed_mps` and `angle_rad`
above the local horizontal from `start_radius_m`, under constant gravity
`gravity_m_s2` and the drag of the plume as it decays away from the impingement
point.

`gas_source(radius_m, height_m)` returns the gas the grain is flying through as
the usual six-field NamedTuple, with the gas velocity taken radially outward at
`speed_mps`. [`ejecta_distribution`](@ref) builds one from a plume field, fading
the surface density as `exp(-z/δ)` with `δ = max(wall_jet_growth_rate * r,
wall_jet_thickness_min_m)`; pass a source returning
`ejecta_vacuum_gas_state` to get the pure ballistic problem, whose range
is the textbook `v² sin(2θ)/g`.

The surface is flat: at these ranges (hundreds of meters against a 1737 km lunar
radius) the curvature correction is below a part in 10⁴, and the terrain, not the
figure of the body, is what would actually matter.

Returns `(range_m, apex_m, flight_time_s, impact_speed_mps, final_speed_mps,
escaped, steps)`. `escaped` is set when the grain's speed reaches
`escape_speed_mps` (see [`ejecta_escape_speed`](@ref)), in which case `range_m`
is `Inf`. Allocation-free, given an allocation-free `gas_source`.
"""
function ejecta_trajectory(gas_source::G, soil, diameter_m::Real, launch_speed_mps::Real,
                           angle_rad::Real, gravity_m_s2::Real;
                           config::EjectaTransportConfig=EjectaTransportConfig(),
                           start_radius_m::Real=0.0, escape_speed_mps::Real=Inf) where {G}
    d = Float64(diameter_m)
    v0 = Float64(launch_speed_mps)
    ang = Float64(angle_rad)
    g = Float64(gravity_m_s2)
    r = Float64(start_radius_m)
    v_esc = Float64(escape_speed_mps)
    (d > 0.0 && v0 > 0.0 && g > 0.0 && isfinite(ang)) ||
        return (range_m=r, apex_m=0.0, flight_time_s=0.0, impact_speed_mps=0.0,
                final_speed_mps=0.0, escaped=false, steps=0)
    v0 >= v_esc && return (range_m=Inf, apex_m=Inf, flight_time_s=Inf, impact_speed_mps=0.0,
                           final_speed_mps=v0, escaped=true, steps=0)
    sn, cs = sincos(ang)
    vr = v0 * cs
    vz = v0 * sn
    z = 0.0
    t = 0.0
    apex = 0.0
    # The drag-free flight time is the natural step scale; drag can only stretch
    # the flight, which the step cap absorbs.
    t_ballistic = 2.0 * v0 * max(sn, 1.0e-6) / g
    steps_nominal = max(config.trajectory_steps, 8)
    dt = t_ballistic / steps_nominal
    max_steps = max(config.trajectory_max_steps, steps_nominal)
    used = 0
    @inbounds for step in 1:max_steps
        used = step
        r0, z0, vr0, vz0 = r, z, vr, vz
        a1r, a1z = _flight_accel(gas_source, config, soil, d, r0, z0, vr0, vz0, g)
        r2, z2, vr2, vz2 = r0 + 0.5dt * vr0, z0 + 0.5dt * vz0, vr0 + 0.5dt * a1r, vz0 + 0.5dt * a1z
        a2r, a2z = _flight_accel(gas_source, config, soil, d, r2, z2, vr2, vz2, g)
        r3, z3, vr3, vz3 = r0 + 0.5dt * vr2, z0 + 0.5dt * vz2, vr0 + 0.5dt * a2r, vz0 + 0.5dt * a2z
        a3r, a3z = _flight_accel(gas_source, config, soil, d, r3, z3, vr3, vz3, g)
        r4, z4, vr4, vz4 = r0 + dt * vr3, z0 + dt * vz3, vr0 + dt * a3r, vz0 + dt * a3z
        a4r, a4z = _flight_accel(gas_source, config, soil, d, r4, z4, vr4, vz4, g)
        r = r0 + (dt / 6.0) * (vr0 + 2.0vr2 + 2.0vr3 + vr4)
        z = z0 + (dt / 6.0) * (vz0 + 2.0vz2 + 2.0vz3 + vz4)
        vr = vr0 + (dt / 6.0) * (a1r + 2.0a2r + 2.0a3r + a4r)
        vz = vz0 + (dt / 6.0) * (a1z + 2.0a2z + 2.0a3z + a4z)
        t += dt
        apex = max(apex, z)
        speed = hypot(vr, vz)
        if speed >= v_esc
            return (range_m=Inf, apex_m=Inf, flight_time_s=Inf, impact_speed_mps=0.0,
                    final_speed_mps=speed, escaped=true, steps=used)
        end
        if z <= 0.0
            if !(z0 > 0.0) && step == 1
                # Launched into the ground: no flight at all.
                return (range_m=r0, apex_m=0.0, flight_time_s=0.0, impact_speed_mps=v0,
                        final_speed_mps=v0, escaped=false, steps=used)
            end
            # z(s) = z0 + vz0 s + a1z s^2 / 2 on the step just taken; exact when
            # the only vertical force is gravity, which is the ballistic limit.
            s = _quadratic_crossing(z0, vz0, a1z, dt)
            frac = dt > 0.0 ? clamp(s / dt, 0.0, 1.0) : 1.0
            r_land = r0 + frac * (r - r0)
            vr_land = vr0 + frac * (vr - vr0)
            vz_land = vz0 + frac * (vz - vz0)
            t += (frac - 1.0) * dt
            return (range_m=max(r_land, 0.0), apex_m=apex, flight_time_s=max(t, 0.0),
                    impact_speed_mps=hypot(vr_land, vz_land),
                    final_speed_mps=hypot(vr_land, vz_land), escaped=false, steps=used)
        end
        r >= config.max_flight_radius_m && break
    end
    return (range_m=r, apex_m=apex, flight_time_s=t, impact_speed_mps=hypot(vr, vz),
            final_speed_mps=hypot(vr, vz), escaped=false, steps=used)
end

@inline function _flight_accel(gas_source::G, config::EjectaTransportConfig, soil, d::Float64,
                               r::Float64, z::Float64, vr::Float64, vz::Float64, g::Float64) where {G}
    gas = gas_source(r, z)
    rho = Float64(gas.density_kg_m3)
    rho > 0.0 || return (0.0, -g)
    wr = Float64(gas.speed_mps) - vr
    wz = -vz
    w = hypot(wr, wz)
    w > 0.0 || return (0.0, -g)
    cd = ejecta_drag_coefficient(config, gas, w, d)
    cd > 0.0 || return (0.0, -g)
    k = 3.0 * rho * cd * w / (4.0 * Float64(soil.particle_density_kg_m3) * d)
    return (k * wr, k * wz - g)
end

"""
    ejecta_escape_speed(gravity_m_s2, body_radius_m) -> Float64

Escape speed `sqrt(2 g R)` at the surface of a body of surface gravity `g` and
radius `R`. For the Moon (1.62 m/s², 1737.4 km) this is 2373 m/s, the value
Immer et al. (2011) quote when bounding the Surveyor 3 impact speeds.
"""
@inline ejecta_escape_speed(gravity_m_s2::Real, body_radius_m::Real)::Float64 =
    sqrt(2.0 * Float64(gravity_m_s2) * Float64(body_radius_m))

# ---- distribution --------------------------------------------------------------------

"""
Default grain sizes for [`ejecta_distribution`](@ref): 1, 5, 20, 70, 200 and
500 µm. The 70 µm entry is the Lunar Sourcebook's average median size (chapter
9, section 9.1.1); the range spans the lunar soil distribution, whose fines
reach far below a micrometre (Park, Liu, Kihm and Taylor, "Characterization of
lunar dust for toxicological studies I: particle size distribution",
J. Aerospace Engineering 21(4), 2008, give 0.019 µm as the smallest size).
"""
const EJECTA_DEFAULT_SIZES_M = (1.0e-6, 5.0e-6, 2.0e-5, 7.0e-5, 2.0e-4, 5.0e-4)

"""
    ejecta_distribution(field, config, soil, thrust_n, height_m; kwargs...) -> NamedTuple

The population of grains a plume of `thrust_n` newtons standing `height_m` above
the ground throws off the surface: histograms of launch speed, ejection angle
and deposition radius, plus the summary scalars a results table would carry.

`field` is anything [`ejecta_gas_state`](@ref) accepts -- a
[`EjectaReferenceGasField`](@ref), or the closure
`(cfg, F, h, r) -> plume_gas_state(plume_field, cfg, F, h, r)` over the plume
field module's table.

For every launch radius in `radii` and grain size in `sizes` the model
(1) reads the gas on the ground there, (2) accelerates the grain over
`max(config.entrainment_length_factor * r, config.entrainment_length_min_m)` of
surface with [`ejecta_launch_speed`](@ref), (3) assigns it an ejection angle
spread uniformly over `config.ejection_angle_min_deg` to
`ejection_angle_max_deg` (the Apollo film range; see the module header), and
(4) flies it with [`ejecta_trajectory`](@ref) through the decaying plume.

Keywords: `sizes` (`EJECTA_DEFAULT_SIZES_M`), `radii` (16 points spanning
three footprint radii of the analytic footprint `h tan 25°`, floored at 0.25 m),
`size_weights` (`nothing`, meaning equal mass per size bin),
`gravity_m_s2` (1.62, the Moon), `body_radius_m` (1.7374e6, the Moon),
`speed_bins`/`angle_bins`/`deposition_bins` (24/12/24), and `weight`, a callable
`(gas, radius_m, diameter_m) -> Float64` giving the mass each cell contributes.

The default `weight` is `gas.shear_pa * 2π r Δr` times the size weight. THIS IS
A PROXY, not an erosion model: it follows Roberts' viscous erosion in making the
entrained mass proportional to the local wall shear stress, but it does not
apply a threshold, so the histograms describe "where the shear stress is" rather
than "where soil actually moves". Pass the regolith erosion module's rate as
`weight` to fix that; the summary scalars change, the trajectories do not.

Returns a NamedTuple whose scalar members are the candidates for result columns:
`mean_speed_mps`, `median_speed_mps`, `max_speed_mps`, `mean_angle_deg`,
`mean_deposition_radius_m`, `p90_deposition_radius_m`, `max_deposition_radius_m`,
`escape_fraction`, `escape_speed_mps`, `mean_knudsen`, `dominant_regime`,
`sample_count` and `total_weight`, alongside the histogram edges and the
mass fractions in each bin (`speed_edges_mps`/`speed_fraction`,
`angle_edges_deg`/`angle_fraction`,
`deposition_edges_m`/`deposition_fraction`). This builder allocates; the
per-grain functions it calls do not.
"""
function ejecta_distribution(field, config::EjectaTransportConfig, soil, thrust_n::Real,
                             height_m::Real;
                             sizes=EJECTA_DEFAULT_SIZES_M,
                             radii=nothing,
                             size_weights=nothing,
                             gravity_m_s2::Real=1.62,
                             body_radius_m::Real=1.7374e6,
                             speed_bins::Integer=24,
                             angle_bins::Integer=12,
                             deposition_bins::Integer=24,
                             weight=nothing)
    F = Float64(thrust_n)
    h = Float64(height_m)
    g = Float64(gravity_m_s2)
    v_esc = ejecta_escape_speed(g, body_radius_m)
    size_vec = collect(Float64, sizes)
    radius_vec = radii === nothing ? _default_radii(h) : collect(Float64, radii)
    n_sizes = length(size_vec)
    n_radii = length(radius_vec)
    (n_sizes >= 1 && n_radii >= 1) ||
        throw(ArgumentError("ejecta_distribution needs at least one size and one radius"))
    w_size = size_weights === nothing ? fill(1.0 / n_sizes, n_sizes) : collect(Float64, size_weights)
    length(w_size) == n_sizes ||
        throw(ArgumentError("size_weights must have one entry per size"))

    n = n_sizes * n_radii
    speeds = zeros(n)
    angles = zeros(n)
    ranges = zeros(n)
    weights = zeros(n)
    knudsens = zeros(n)
    regimes = zeros(Int, 5)
    escaped_weight = 0.0
    total_weight = 0.0
    idx = 0
    gas_source = _EjectaFlightGas(field, config, F, h)
    for (ir, r) in enumerate(radius_vec)
        dr = _radial_width(radius_vec, ir)
        gas = ejecta_gas_state(field, config, F, h, r)
        entrain = max(config.entrainment_length_factor * r, config.entrainment_length_min_m)
        for (is, d) in enumerate(size_vec)
            idx += 1
            launch = ejecta_launch_speed(gas, soil, d, entrain; config=config)
            # The ejection angle is an input, spread uniformly over the measured
            # Apollo range. The van der Corput sequence samples that range
            # deterministically (no RNG, reproducible run to run) without letting
            # the angle correlate with the radius or the size the sweep happens
            # to be walking.
            frac = _radical_inverse_base2(idx)
            ang_deg = config.ejection_angle_min_deg +
                      frac * (config.ejection_angle_max_deg - config.ejection_angle_min_deg)
            traj = ejecta_trajectory(gas_source, soil, d, launch.speed_mps, deg2rad(ang_deg), g;
                                     config=config, start_radius_m=r, escape_speed_mps=v_esc)
            w = weight === nothing ?
                Float64(gas.shear_pa) * 2.0 * pi * r * dr * w_size[is] :
                Float64(weight(gas, r, d)) * w_size[is]
            w = max(w, 0.0)
            speeds[idx] = launch.speed_mps
            angles[idx] = ang_deg
            ranges[idx] = traj.escaped ? Inf : traj.range_m
            weights[idx] = w
            knudsens[idx] = isfinite(launch.knudsen) ? launch.knudsen : 0.0
            code = Int(launch.regime)
            1 <= code <= 4 && (regimes[code] += 1)
            total_weight += w
            traj.escaped && (escaped_weight += w)
        end
    end

    finite_ranges = [isfinite(ranges[i]) ? ranges[i] : 0.0 for i in 1:n]
    ground_weights = [isfinite(ranges[i]) ? weights[i] : 0.0 for i in 1:n]
    speed_edges, speed_fraction = _mass_histogram(speeds, weights, speed_bins)
    angle_edges, angle_fraction = _mass_histogram(angles, weights, angle_bins)
    dep_edges, dep_fraction = _mass_histogram(finite_ranges, ground_weights, deposition_bins)
    dominant = sum(regimes) > 0 ? argmax(view(regimes, 1:4)) : 0

    return (
        speed_edges_mps=speed_edges, speed_fraction=speed_fraction,
        angle_edges_deg=angle_edges, angle_fraction=angle_fraction,
        deposition_edges_m=dep_edges, deposition_fraction=dep_fraction,
        mean_speed_mps=_weighted_mean(speeds, weights),
        median_speed_mps=_weighted_quantile(speeds, weights, 0.5),
        max_speed_mps=maximum(speeds),
        mean_angle_deg=_weighted_mean(angles, weights),
        mean_deposition_radius_m=_weighted_mean(finite_ranges, ground_weights),
        p90_deposition_radius_m=_weighted_quantile(finite_ranges, ground_weights, 0.9),
        max_deposition_radius_m=maximum(finite_ranges),
        escape_fraction=total_weight > 0.0 ? escaped_weight / total_weight : 0.0,
        escape_speed_mps=v_esc,
        mean_knudsen=_weighted_mean(knudsens, weights),
        dominant_regime=ejecta_regime_name(UInt8(dominant)),
        sample_count=n,
        total_weight=total_weight,
    )
end

"""
    _EjectaFlightGas(field, config, thrust_n, height_m)

Callable that gives [`ejecta_trajectory`](@ref) the gas a grain is flying
through: the surface state at the grain's radius, with the density faded as
`exp(-z/δ)`, `δ = max(wall_jet_growth_rate * r, wall_jet_thickness_min_m)`. A
callable struct rather than a closure so the call is concrete and allocation-free.
"""
struct _EjectaFlightGas{F, C}
    field::F
    config::C
    thrust_n::Float64
    height_m::Float64
end

_EjectaFlightGas(field, config, thrust_n::Real, height_m::Real) =
    _EjectaFlightGas(field, config, Float64(thrust_n), Float64(height_m))

@inline function (gs::_EjectaFlightGas)(radius_m::Float64, height_above_surface_m::Float64)
    gas = ejecta_gas_state(gs.field, gs.config, gs.thrust_n, gs.height_m, max(radius_m, 0.0))
    z = max(height_above_surface_m, 0.0)
    delta = max(gs.config.wall_jet_growth_rate * max(radius_m, 0.0),
                gs.config.wall_jet_thickness_min_m)
    fade = exp(-z / delta)
    return (pressure_pa=gas.pressure_pa, shear_pa=gas.shear_pa,
            density_kg_m3=gas.density_kg_m3 * fade, speed_mps=gas.speed_mps,
            temperature_k=gas.temperature_k, mach=gas.mach)
end

"""
    _radical_inverse_base2(i) -> Float64

The `i`-th term of the van der Corput sequence in base 2, in `[0, 1)`: the bits
of `i` reversed and read as a binary fraction. Used to spread the ejection angle
over its measured range without an RNG and without correlating it with the order
the (radius, size) grid is walked in.
"""
@inline function _radical_inverse_base2(i::Integer)::Float64
    bits = UInt32(i) & 0xFFFFFFFF
    bits = (bits << 16) | (bits >> 16)
    bits = ((bits & 0x55555555) << 1) | ((bits & 0xAAAAAAAA) >> 1)
    bits = ((bits & 0x33333333) << 2) | ((bits & 0xCCCCCCCC) >> 2)
    bits = ((bits & 0x0F0F0F0F) << 4) | ((bits & 0xF0F0F0F0) >> 4)
    bits = ((bits & 0x00FF00FF) << 8) | ((bits & 0xFF00FF00) >> 8)
    return Float64(bits) * 2.3283064365386963e-10
end

function _default_radii(height_m::Float64)
    R = max(height_m * tand(25.0), 0.75)
    return collect(range(max(0.25, 0.05 * R), 3.0 * R; length=16))
end

@inline function _radial_width(radii::Vector{Float64}, i::Int)::Float64
    n = length(radii)
    n == 1 && return max(radii[1], 1.0e-6)
    i == 1 && return radii[2] - radii[1]
    i == n && return radii[n] - radii[n - 1]
    return 0.5 * (radii[i + 1] - radii[i - 1])
end

function _weighted_mean(values::Vector{Float64}, weights::Vector{Float64})::Float64
    total = 0.0
    acc = 0.0
    @inbounds for i in eachindex(values)
        w = weights[i]
        w > 0.0 || continue
        total += w
        acc += w * values[i]
    end
    return total > 0.0 ? acc / total : 0.0
end

function _weighted_quantile(values::Vector{Float64}, weights::Vector{Float64}, q::Float64)::Float64
    total = sum(weights)
    total > 0.0 || return 0.0
    order = sortperm(values)
    target = q * total
    acc = 0.0
    @inbounds for i in order
        acc += weights[i]
        acc >= target && return values[i]
    end
    return values[order[end]]
end

function _mass_histogram(values::Vector{Float64}, weights::Vector{Float64}, nbins::Integer)
    nb = max(Int(nbins), 1)
    lo = Inf
    hi = -Inf
    @inbounds for i in eachindex(values)
        weights[i] > 0.0 || continue
        lo = min(lo, values[i])
        hi = max(hi, values[i])
    end
    (isfinite(lo) && isfinite(hi)) || return (collect(range(0.0, 1.0; length=nb + 1)), zeros(nb))
    hi > lo || (hi = lo + 1.0)
    edges = collect(range(lo, hi; length=nb + 1))
    counts = zeros(nb)
    total = 0.0
    width = (hi - lo) / nb
    @inbounds for i in eachindex(values)
        w = weights[i]
        w > 0.0 || continue
        b = clamp(Int(floor((values[i] - lo) / width)) + 1, 1, nb)
        counts[b] += w
        total += w
    end
    total > 0.0 && (counts ./= total)
    return (edges, counts)
end

end # module EjectaTransport
