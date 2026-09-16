# The gas a descent engine's plume lays on the ground: the field the
# plume-surface interaction model reads its surface state from.
#
# Two fields answer the same query,
# `plume_gas_state(field, config, thrust_n, height_m, radius_m) -> PlumeGasState`,
# with `height_m` the nozzle exit plane above the ground and `radius_m` measured
# from the stagnation point on the ground:
#
#   * `PlumeAnalyticField` is the Gaussian pressure footprint the effector has
#     always used. It owns no thermodynamics, so its density, temperature and
#     Mach number are reconstructed from its own pressure through the same
#     impingement relations the table uses. It stays the default, so a run that
#     does not ask for a table reproduces earlier results exactly.
#   * `PlumeFieldTable` is a gridded evaluation of a far-field plume computed
#     from the nozzle, built once by `scripts/dev/psi/build_plume_field.jl` and
#     loaded with `load_plume_field`.
#
# WHERE THE PHYSICS COMES FROM
#
# Plume: the Simons source-flow model of a rocket exhausting into vacuum with
# Boynton's nozzle-boundary-layer correction,
#   G. A. Simons, "Effect of nozzle boundary layers on rocket exhaust plumes",
#     AIAA Journal 10(11), 1972, pp. 1534-1535;
#   F. P. Boynton, "Exhaust plumes from nozzles with wall boundary layers",
#     J. Spacecraft and Rockets 5(10), 1968, pp. 1143-1147,
# in the explicit form collected by
#   P. J. Herraiz, J. M. Fernandez and J. R. Villa, "Development of a MATLAB
#     plume impingement tool for fast system analysis", 8th European Conference
#     for Aeronautics and Aerospace Sciences (EUCASS), 2019,
#     DOI 10.13009/EUCASS2019-661,
# whose section 2.1 equations (3) to (14) are the ones implemented here and are
# cited below equation by equation.
#
# Impingement: classical Newtonian impact theory for the wall pressure
# (`p = rho_inf U^2 cos^2(theta)`, the `C_p,max = 2` limit) and the oblique
# normal-shock relations for the state just above the wall. Newtonian is chosen
# over modified Newtonian because it satisfies the global axial-momentum balance
# on the plane exactly: the integral of the wall pressure over the ground
# recovers the axial momentum flux of the plume identically, which is the
# property `plume_surface_footprint` and the momentum test rely on. Modified
# Newtonian would scale it by `C_p,max/2`, which is 0.95 for this exhaust.
#
# Wall shear: the wall-jet dynamic pressure times a rough-surface drag
# coefficient, the treatment L. Roberts used for a meteorite-gardened lunar
# surface (Roberts 1966, quoted with its value of 0.2 by A. B. Morris,
# "Simulation of rocket plume impingement and dust dispersal on the lunar
# surface", PhD dissertation, University of Texas at Austin, 2012, section
# 4.4.2). The `PlumeAnalyticField` keeps its own skin-friction closure instead,
# unchanged.
#
# ASSUMPTIONS. Values that no source pinned down are named fields of
# `PlumeNozzle` with the reasoning for their magnitude on the field: the nozzle
# exit lip half-angle `exit_lip_angle_deg`, the mean limit-speed ratio
# `limit_speed_ratio` of equation (14), and the divergent length
# `divergent_length_m`. They are documented as assumptions here, in
# `docs/src/user/lunar_landing.md`, and in `data/psi/README.md`.
module PlumeGasField

using JSON

export PlumeGasState, PlumeNozzle, PlumeAnalyticField, PlumeFieldTable
export plume_gas_state, plume_field_footprint, plume_field_shear_coefficient
export plume_wall_shear, plume_mean_shear, plume_scour_radius
export plume_field_name, build_plume_field_table, save_plume_field, load_plume_field
export plume_limit_speed, plume_limit_angle, plume_angular_mass_flux

"""
    PlumeGasState

The gas on the ground under an impinging plume, at one radius from the
stagnation point: `pressure_pa` (static pressure on the surface), `shear_pa`
(wall shear stress), `density_kg_m3` and `speed_mps` (the wall jet just above
the surface, the speed being parallel to it), `temperature_k` and `mach` (the
wall jet's local Mach number).
"""
const PlumeGasState = @NamedTuple{
    pressure_pa::Float64,
    shear_pa::Float64,
    density_kg_m3::Float64,
    speed_mps::Float64,
    temperature_k::Float64,
    mach::Float64,
}

const _PLUME_ZERO_STATE = PlumeGasState((0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
const _G0 = 9.80665                    # standard gravity, the definition of specific impulse
const _R_UNIVERSAL = 8314.462618       # J/(kmol K), CODATA 2018

# ---- engine and propellant -----------------------------------------------------------

"""
    PlumeNozzle(; ...)

The engine the plume is computed from. The defaults are the Apollo lunar module
descent engine (LMDE) in the Apollo 11 configuration, with every number either
sourced or derived from sourced numbers:

SOURCED

- `specific_impulse_s` = 305 s and `area_ratio` = 47.5: A. B. Morris,
  "Simulation of rocket plume impingement and dust dispersal on the lunar
  surface", PhD dissertation, University of Texas at Austin, 2012, Table 2.2.1;
  the same two numbers appear in "Rocket Propulsion Evolution" 9.42
  (enginehistory.org), which quotes the TRW LMDE data sheet.
- `exit_mach` = 5.03 and `chamber_temperature_k` = 2730 K: Morris 2012
  Table 3.1.1 ("Estimated LMDE exit plane and throat properties"), computed there
  from the engine's rated performance with the JANNAF exhaust composition.
- `exit_radius_m` = 0.75 m: the 59 inch exit diameter of the area-ratio 47.5
  nozzle Apollo 11 flew ("Rocket Propulsion Evolution" 9.42; Wikipedia's
  "Descent propulsion system" gives the same 59 to 63 inch range for the two
  nozzles). Morris 2012 Table 2.2.1 lists 1.62 m instead; the exit diameter only
  sets the table's height axis, so the difference is cosmetic.
- `surface_drag_coefficient` = 0.2: the rough-surface drag coefficient L.
  Roberts applied to the wall dynamic pressure for the meteorite-gardened lunar
  surface (Roberts 1966), quoted as 0.2 by Morris 2012 section 4.4.2.

DERIVED (the derivation is in the comment, so it can be checked)

- `gamma` = 1.30 and `gas_constant_j_kg_k` = 461 J/(kg K): from Morris 2012
  Table 3.1.1's exit state (M = 5.03, rho = 1.09e-3 kg/m^3, p = 296 Pa,
  T = 589 K, V = 2992 m/s). `R = p/(rho T)` = 461 J/(kg K), which is a mean
  exhaust molar mass of 18.0 kg/kmol, matching that dissertation's statement
  that the mean molecular weight of the LMDE exhaust is close to water's;
  `gamma = (V/M)^2/(R T)` = 1.30. The pair is self-consistent with the area
  ratio: the isentropic area-Mach relation at `gamma = 1.30` gives `M = 5.03`
  for `A_e/A_t = 47.5` to three digits.
- `divergent_length_m` = 1.91 m: an 80 percent bell for `A_e/A_t` = 47.5 and a
  1.50 m exit, i.e. 0.8 times the 15 degree conical length
  `(D_e - D_t)/(2 tan 15 deg)` with `D_t = D_e/sqrt(47.5)`.
- `exit_boundary_layer_ratio` = 0.0595, the exit boundary-layer thickness over
  the exit radius: from EUCASS 2019 equation (9), `delta = 10.82 l/sqrt(Re_l)`,
  with `Re_l` scaled from the exit Reynolds number 182,000 that Morris 2012
  Table 3.1.1 reports for the LMDE at 13.34 kN (the thrust Apollo flew in the
  last minute of the descent). The 10.82 in that correlation is itself a fit
  reported by the EUCASS authors, not a universal boundary-layer constant, and
  the thickness scales as `thrust^-1/2`, so at full throttle it would be 0.032.

ASSUMPTIONS (no source found; each is a named field so a study can sweep it)

- `exit_lip_angle_deg` = 10 deg, the flow angle at the nozzle lip. Rao-contoured
  bells of this area ratio exit between about 5 and 15 degrees; 10 degrees is
  the middle of that band. It enters the limiting turn angle of equation (7) and
  the near-field terms of equations (5) and (6).
- `limit_speed_ratio` = 0.75, the `vbar_lim/v_lim` of EUCASS 2019 equation (14),
  which that paper bounds only as `v_lim/2 < vbar_lim < v_lim`. 0.75 is the
  midpoint of the stated bound. It sets how fast the boundary-layer wing of the
  plume decays outboard of `theta_0`.

NEAR FIELD. Two optional treatments of the region where the plume is not yet a
point source, both off by default, so the shipped table is the plain
source-flow model with its source at the center of the exit plane:

- `near_field_correction` switches on the `r^2 - a r + b` denominator of EUCASS
  2019 equation (3), with `a = 3 sqrt(theta_e) M_e r_e` (equation 5) and
  `b = 5 theta_e M_e^2 r_e^2` (equation 6). It is off because it is a density
  multiplier, not a geometric change, so it does not conserve the plume's axial
  momentum against a surface: with these LMDE numbers the integral of the wall
  pressure over the ground comes to 1.81 times the thrust at a 5 m height, 1.25
  at 20 m and 1.10 at 50 m, tending to 1 only in the far field. It stays
  available because it is the published form.
- `virtual_source_offset_m` moves the point source that far downstream of the
  exit plane along the axis, which is geometric and therefore conserves
  momentum exactly. Morris 2012 section 4.6 measured this offset at about 1.8 m
  for the LMDE, where the nozzle's internal compression wave meets the axis,
  and identified the point source sitting at the exit plane as the reason
  Roberts' theory under-predicts the surface stress below about ten nozzle
  diameters. It defaults to zero because it is only meaningful while the source
  stays above the ground: the source-to-ground distance is floored at one exit
  radius, and below that the table reports the frozen value.

`reference_thrust_n` is only the normalization the table is written at: every
tabulated quantity that carries thrust is stored divided by it, so the choice
cancels (see `PlumeFieldTable`).
"""
Base.@kwdef struct PlumeNozzle
    name::String = "apollo_lmde"
    gamma::Float64 = 1.30
    gas_constant_j_kg_k::Float64 = 461.0
    chamber_temperature_k::Float64 = 2730.0
    exit_mach::Float64 = 5.03
    area_ratio::Float64 = 47.5
    exit_radius_m::Float64 = 0.75
    exit_lip_angle_deg::Float64 = 10.0
    divergent_length_m::Float64 = 1.91
    exit_boundary_layer_ratio::Float64 = 0.0595
    limit_speed_ratio::Float64 = 0.75
    specific_impulse_s::Float64 = 305.0
    reference_thrust_n::Float64 = 45_040.0
    surface_drag_coefficient::Float64 = 0.2
    near_field_correction::Bool = false
    virtual_source_offset_m::Float64 = 0.0
end

"""
    plume_limit_speed(nozzle) -> Float64

Limiting (fully expanded) exhaust speed, `sqrt(2 gamma R T_0/(gamma - 1))`: all
of the chamber's enthalpy turned into directed kinetic energy. The Simons source
flow carries this speed along every streamline. For the LMDE defaults it is
3303 m/s, 10 percent above the 2992 m/s exit speed of Morris 2012 Table 3.1.1,
as the remaining expansion into vacuum should give.
"""
@inline function plume_limit_speed(n::PlumeNozzle)::Float64
    return sqrt(2.0 * n.gamma * n.gas_constant_j_kg_k * n.chamber_temperature_k / (n.gamma - 1.0))
end

"Stagnation density of the chamber, from the choked throat and the mass flow of `thrust_n`."
@inline function _chamber_density(n::PlumeNozzle, thrust_n::Float64)::Float64
    mdot = thrust_n / (n.specific_impulse_s * _G0)
    a_throat = sqrt(n.gamma * n.gas_constant_j_kg_k * n.chamber_temperature_k * 2.0 / (n.gamma + 1.0))
    area_throat = pi * n.exit_radius_m^2 / n.area_ratio
    rho_throat = mdot / (area_throat * a_throat)
    return rho_throat * ((n.gamma + 1.0) / 2.0)^(1.0 / (n.gamma - 1.0))
end

"Prandtl-Meyer function (radians) of a Mach number, and its `M -> Inf` limit."
@inline function _prandtl_meyer(mach::Float64, gamma::Float64)::Float64
    mach <= 1.0 && return 0.0
    b = sqrt((gamma + 1.0) / (gamma - 1.0))
    m = sqrt(mach * mach - 1.0)
    return b * atan(m / b) - atan(m)
end

@inline _prandtl_meyer_limit(gamma::Float64)::Float64 = 0.5 * pi * (sqrt((gamma + 1.0) / (gamma - 1.0)) - 1.0)

"""
    plume_limit_angle(nozzle) -> (theta_lim_rad, theta_0_rad, decay_per_rad)

The three angles of the Simons/Boynton angular distribution:

- `theta_lim = theta_e + (nu(Inf) - nu(M_e)) - theta_BL` is the limiting
  streamline angle, EUCASS 2019 equation (7), with the boundary-layer deflection
  `theta_BL = atan(delta cos(theta_e)/l)` of equation (11);
- `theta_0`, the edge of the isentropic core, from equation (12),
  `theta_0/theta_lim = (2/pi) acos[(2 delta/r_e - (delta/r_e)^2)^((g-1)/(g+1))]`;
- `decay_per_rad` is the `c_rho` of equation (14), the exponential decay rate of
  the boundary-layer wing outboard of `theta_0`.

For the LMDE defaults these are 78.4 and 35.7 degrees and about 9.9 per radian.
"""
function plume_limit_angle(n::PlumeNozzle)
    g = n.gamma
    theta_e = deg2rad(n.exit_lip_angle_deg)
    zeta = n.exit_boundary_layer_ratio
    delta = zeta * n.exit_radius_m
    theta_bl = atan(delta * cos(theta_e) / n.divergent_length_m)
    theta_lim = theta_e + (_prandtl_meyer_limit(g) - _prandtl_meyer(n.exit_mach, g)) - theta_bl
    theta_lim = clamp(theta_lim, 1e-3, pi)
    arg = clamp((2.0 * zeta - zeta * zeta)^((g - 1.0) / (g + 1.0)), -1.0, 1.0)
    theta_0 = clamp((2.0 / pi) * acos(arg), 0.0, 1.0) * theta_lim
    # A_P of equation (4): the closed-form normalization of the isentropic core.
    core(t) = cos(0.5 * pi * t / theta_lim)^(2.0 / (g - 1.0))
    integral = _simpson(t -> core(t) * sin(t), 0.0, theta_lim, 512)
    a_p = integral > 0.0 ? 0.5 * sqrt((g - 1.0) / (g + 1.0)) / integral : 0.0
    c_rho = a_p * sqrt((g + 1.0) / (g - 1.0)) * (2.0 * n.limit_speed_ratio) *
            (n.exit_radius_m / (2.0 * delta))^((g - 1.0) / (g + 1.0))
    return theta_lim, theta_0, c_rho
end

"Composite Simpson rule on `n` (even) intervals."
function _simpson(f, a::Float64, b::Float64, n::Int)::Float64
    n = iseven(n) ? n : n + 1
    h = (b - a) / n
    s = f(a) + f(b)
    for k in 1:(n - 1)
        s += (isodd(k) ? 4.0 : 2.0) * f(a + k * h)
    end
    return s * h / 3.0
end

"""
    plume_angular_mass_flux(nozzle, thrust_n) -> (theta -> dmdot/dOmega)

The plume's mass flux per steradian as a function of the polar angle from the
engine axis: the `g(theta)` of EUCASS 2019 equation (13) -- a
`cos^(2/(gamma-1))` core out to `theta_0`, an exponential wing decaying at
`c_rho` out to `theta_lim`, and that wing's end value held beyond it -- scaled
so that the integral over the sphere is exactly the engine's mass flow
`thrust/(Isp g0)`. Normalizing numerically rather than with the closed-form
`A_P` of equation (4) makes mass conservation exact for the piecewise shape,
including the boundary-layer wing that equation (4) leaves out.
"""
function plume_angular_mass_flux(n::PlumeNozzle, thrust_n::Real)
    theta_lim, theta_0, c_rho = plume_limit_angle(n)
    g = n.gamma
    core(t) = cos(0.5 * pi * t / theta_lim)^(2.0 / (g - 1.0))
    core_0 = core(theta_0)
    shape(t) = t <= theta_0 ? core(t) :
               (t <= theta_lim ? core_0 * exp(-c_rho * (t - theta_0)) :
                core_0 * exp(-c_rho * (theta_lim - theta_0)))
    total = 2.0 * pi * _simpson(t -> shape(t) * sin(t), 0.0, Float64(pi), 2048)
    mdot = Float64(thrust_n) / (n.specific_impulse_s * _G0)
    scale = total > 0.0 ? mdot / total : 0.0
    return t -> scale * shape(t)
end

"""
    _near_field_denominator(nozzle, slant_m) -> Float64

The `r^2 - a r + b` of EUCASS 2019 equation (3), with `a = 3 sqrt(theta_e) M_e r_e`
(equation 5) and `b = 5 theta_e M_e^2 r_e^2` (equation 6), `theta_e` in radians.
It replaces the `1/r^2` of a bare point source close to the nozzle: the
discriminant `a^2 - 4b` is negative for these engines, so the denominator never
vanishes and the density stays finite at the source. It is used only when
`near_field_correction` is set; `PlumeNozzle` records why it is off by default.
"""
@inline function _near_field_denominator(n::PlumeNozzle, slant_m::Float64)::Float64
    s2 = slant_m * slant_m
    n.near_field_correction || return s2
    theta_e = deg2rad(n.exit_lip_angle_deg)
    a = 3.0 * sqrt(theta_e) * n.exit_mach * n.exit_radius_m
    b = 5.0 * theta_e * n.exit_mach^2 * n.exit_radius_m^2
    return max(s2 - a * slant_m + b, 1e-6)
end

"""
    _source_height(nozzle, height_m) -> Float64

Height of the plume's point source above the ground when the nozzle exit plane
is `height_m` above it: the exit height less `virtual_source_offset_m`, floored
at one exit radius so the source never reaches the ground.
"""
@inline function _source_height(n::PlumeNozzle, height_m::Float64)::Float64
    return max(height_m - n.virtual_source_offset_m, n.exit_radius_m)
end

# ---- impingement ---------------------------------------------------------------------

"""
    _impinge(gamma, r_gas, rho_inf, speed, cos_inc, t_total, drag_coefficient) -> PlumeGasState

Turn a free-stream plume state meeting a flat surface at incidence `cos_inc`
(the cosine of the angle between the streamline and the surface normal) into the
gas state on the wall.

- Pressure: classical Newtonian, `p = rho_inf U^2 cos^2(inc)`, the normal
  momentum flux the surface has to absorb.
- Density: the normal-shock density ratio at the normal Mach number
  `M_n = M_inf cos(inc)`, `rho_w/rho_inf = (g+1)M_n^2/((g-1)M_n^2 + 2)`, which
  tends to `(g+1)/(g-1)` = 7.67 for this exhaust.
- Speed: the tangential component `U sin(inc)` survives the shock, so the wall
  jet leaves the stagnation point at rest and accelerates outward.
- Temperature: `p/(rho_w R)`, capped at the plume's total temperature, which the
  gas behind a shock cannot exceed. For this exhaust the cap is active inside
  about 20 degrees of the stagnation point, where the classical Newtonian
  pressure combined with the strong-shock density ratio would otherwise give a
  recovery temperature above the chamber's.
- Shear: `drag_coefficient` times the wall jet's dynamic pressure.
"""
@inline function _impinge(gamma::Float64, r_gas::Float64, rho_inf::Float64, speed::Float64,
                          cos_inc::Float64, mach_inf::Float64, t_total::Float64,
                          drag_coefficient::Float64)::PlumeGasState
    (rho_inf > 0.0 && speed > 0.0) || return _PLUME_ZERO_STATE
    c = clamp(cos_inc, 0.0, 1.0)
    s = sqrt(max(0.0, 1.0 - c * c))
    pressure = rho_inf * speed * speed * c * c
    mn2 = (mach_inf * c)^2
    ratio = mn2 > 0.0 ? (gamma + 1.0) * mn2 / ((gamma - 1.0) * mn2 + 2.0) : 1.0
    rho_w = rho_inf * max(ratio, 1.0)
    u_t = speed * s
    temperature = min(pressure / max(rho_w * r_gas, eps()), t_total)
    a_w = sqrt(max(gamma * r_gas * temperature, eps()))
    return PlumeGasState((pressure, drag_coefficient * 0.5 * rho_w * u_t * u_t, rho_w, u_t, temperature, u_t / a_w))
end

# ---- radial profile accessors --------------------------------------------------------
#
# The surface state is a radial profile, not a single peak, and a measurement of
# a plume on the ground is usually an average or an edge rather than a peak: a
# dust-cloud inversion gives the shear averaged over the region it sees, and
# post-landing photography gives the radius the scour reaches. These three
# accessors make those comparisons like-for-like without reaching into the field.

"""
    plume_wall_shear(field, config, thrust_n, height_m, radius_m) -> Float64

Wall shear stress (Pa) at one radius. The scalar form of
`plume_gas_state(...).shear_pa`.
"""
@inline function plume_wall_shear(field, cfg, thrust_n::Real, height_m::Real, radius_m::Real)::Float64
    return plume_gas_state(field, cfg, thrust_n, height_m, radius_m).shear_pa
end

"""
    plume_mean_shear(field, config, thrust_n, height_m, radius_m; points=256) -> Float64

Wall shear stress averaged over the ground inside `radius_m`,

```math
\\bar\\tau(a) = \\frac{2}{a^2}\\int_0^a \\tau(r)\\, r\\, \\mathrm{d}r
```

by the midpoint rule. This is the quantity a dust-cloud inversion reports, and
it is well below the peak for any impinging plume, whose shear vanishes at the
stagnation point.
"""
function plume_mean_shear(field, cfg, thrust_n::Real, height_m::Real, radius_m::Real; points::Int=256)::Float64
    a = Float64(radius_m)
    (isfinite(a) && a > 0.0 && points >= 1) || return 0.0
    dr = a / points
    total = 0.0
    for k in 1:points
        r = (k - 0.5) * dr
        total += plume_gas_state(field, cfg, thrust_n, height_m, r).shear_pa * r * dr
    end
    return 2.0 * total / (a * a)
end

"""
    plume_scour_radius(field, config, thrust_n, height_m, threshold_pa; max_radius_over_height=8.0, points=512) -> Float64

Outer radius (m) at which the wall shear stress last exceeds `threshold_pa`:
the edge of the region the plume can move, the scour radius that post-landing
surface photography measures. Zero when the plume never reaches the threshold.
The search runs out to `max_radius_over_height` times the height.
"""
function plume_scour_radius(field, cfg, thrust_n::Real, height_m::Real, threshold_pa::Real;
                            max_radius_over_height::Float64=8.0, points::Int=512)::Float64
    h = Float64(height_m)
    (isfinite(h) && h > 0.0 && points >= 1) || return 0.0
    dr = max_radius_over_height * h / points
    outer = 0.0
    for k in 1:points
        r = k * dr
        plume_gas_state(field, cfg, thrust_n, h, r).shear_pa > threshold_pa && (outer = r)
    end
    return outer
end

# ---- the analytic (Gaussian) field ---------------------------------------------------

"""
    PlumeAnalyticField(; gamma, gas_constant_j_kg_k, chamber_temperature_k, limit_speed_mps)

The Gaussian surface-pressure footprint the plume-surface effector has always
used, and the default field of `PlumeSurfaceConfig`. Its pressure and shear
stress are unchanged:

```math
p(r) = \\frac{F}{\\pi R_p^2} e^{-(r/R_p)^2}, \\qquad
\\tau(r) = c_f\\, p(r)\\, \\frac{2r}{R_p}, \\qquad R_p = \\max(h \\tan\\theta_p,\\, r_e)
```

with `R_p`, `c_f` and `theta_p` read from the configuration passed to
`plume_gas_state`, so a run that keeps the default field reproduces earlier
results bit for bit.

The closure has no thermodynamic content, so the other three quantities of
`PlumeGasState` are reconstructed from its own pressure by inverting the same
Newtonian relation the table uses: the free-stream density that would deliver
`p(r)` at the local incidence, then the strong-shock density ratio, the
tangential speed `v_lim sin(theta)` and `T = p/(rho_w R)` capped at the total
temperature. The four thermodynamic constants default to the LMDE values of
`PlumeNozzle`; they affect only the reconstructed density, temperature and Mach
number, never the pressure, the shear stress or any force.
"""
Base.@kwdef struct PlumeAnalyticField
    gamma::Float64 = 1.30
    gas_constant_j_kg_k::Float64 = 461.0
    chamber_temperature_k::Float64 = 2730.0
    limit_speed_mps::Float64 = 3302.6
end

plume_field_name(::PlumeAnalyticField) = "analytic"

"""
    plume_field_shear_coefficient(field, config) -> Float64

The coefficient that turns the wall jet's dynamic pressure into the wall shear
stress for this field, so a caller can recover the dynamic pressure from
`shear_pa` without a second query. It is the configuration's
`friction_coefficient` for the analytic field and the table's own
`surface_drag_coefficient` (Roberts' 0.2) for a table.
"""
@inline plume_field_shear_coefficient(::PlumeAnalyticField, cfg) = Float64(cfg.friction_coefficient)

"Footprint radius of the analytic field: `R_p = max(h tan(theta_p), r_e)`."
@inline function _analytic_radius(cfg, height_m::Float64)::Float64
    return max(height_m * tand(cfg.plume_half_angle_deg), Float64(cfg.nozzle_exit_radius_m))
end

"""
    plume_field_footprint(field, config, thrust_n, height_m) -> (p0_pa, radius_m)

Peak (stagnation) surface pressure and the footprint radius, the radius that
contains `1 - 1/e` of the integral of the surface pressure over the ground. For
the analytic field that definition returns `R_p` exactly, which is what
`plume_surface_footprint` has always returned.
"""
@inline function plume_field_footprint(::PlumeAnalyticField, cfg, thrust_n::Real, height_m::Real)
    R = _analytic_radius(cfg, Float64(height_m))
    return Float64(thrust_n) / (pi * R * R), R
end

"""
    plume_gas_state(field, config, thrust_n, height_m, radius_m) -> PlumeGasState

The gas on the ground at `radius_m` from the stagnation point, with the nozzle
exit plane `height_m` above it and the engine producing `thrust_n`.
"""
function plume_gas_state(f::PlumeAnalyticField, cfg, thrust_n::Real, height_m::Real, radius_m::Real)::PlumeGasState
    F = Float64(thrust_n)
    h = Float64(height_m)
    r = Float64(radius_m)
    (isfinite(F) && F > 0.0 && isfinite(h) && h >= 0.0 && isfinite(r) && r >= 0.0) || return _PLUME_ZERO_STATE
    R = _analytic_radius(cfg, h)
    p0 = F / (pi * R * R)
    x = r / R
    e = exp(-x * x)
    pressure = p0 * e
    shear = (Float64(cfg.friction_coefficient) * p0) * 2.0 * x * e
    pressure > 0.0 || return _PLUME_ZERO_STATE
    slant = sqrt(h * h + r * r)
    cos_inc = slant > 0.0 ? h / slant : 1.0
    u = f.limit_speed_mps
    rho_inf = pressure / max(u * u * cos_inc * cos_inc, eps())
    # Free-stream Mach number from the isentropic relation to the chamber; a
    # strong-shock ratio is reached by M_n of a few, so the exact value only
    # matters near the footprint edge.
    st = _impinge(f.gamma, f.gas_constant_j_kg_k, rho_inf, u, cos_inc, 1.0e3, f.chamber_temperature_k, 1.0)
    return PlumeGasState((pressure, shear, st.density_kg_m3, st.speed_mps, st.temperature_k, st.mach))
end

"Wall shear stress at the normalized radius `x = r/R`, without the `x*R/R` round trip."
@inline function _plume_shear_at(f::PlumeAnalyticField, cfg, thrust_n::Float64, height_m::Float64,
                                 p0::Float64, R::Float64, x::Float64)::Float64
    return (Float64(cfg.friction_coefficient) * p0) * 2.0 * x * exp(-x * x)
end

"""
    _plume_quadrature_limit(field, height_m, radius_m) -> Float64

How far out, in footprint radii, the effector's radial quadrature has to run to
capture the whole shear profile. The Gaussian is down by `exp(-9)` at three
footprint radii and needs no more. A source-flow table decays far more slowly
outboard, so its quadrature runs to the last radius node of the table, which is
where it stops carrying information anyway.
"""
@inline _plume_quadrature_limit(::PlumeAnalyticField, height_m::Float64, R::Float64)::Float64 = 3.0

# ---- the tabulated field -------------------------------------------------------------

"""
    PlumeFieldTable

A plume field evaluated once on a grid and stored, so a run pays an
interpolation instead of a source-flow integration. The axes are
`height_over_diameter` (`h/D_e`, the nozzle exit plane above the ground over the
exit diameter) and `radius_over_height` (`r/h`, the ground radius over that
height), both ascending; a query outside either axis clamps to its end.

Every stored quantity is dimensionless and carries the thrust scaling
explicitly, which is what makes one table serve the whole throttle band:

| stored | recovered as |
|---|---|
| `pressure_hat` | `p = pressure_hat * F / h^2` |
| `shear_hat` | `tau = shear_hat * F / h^2` |
| `density_hat` | `rho = density_hat * F / (v_lim^2 h^2)` |
| `speed_hat` | `u = speed_hat * v_lim` |
| `temperature_hat` | `T = temperature_hat * T_0` |
| `mach` | as stored |

THE THRUST SCALING. The LMDE throttles by moving propellant flow at a fixed
mixture ratio and a fixed nozzle, so the chamber pressure tracks the thrust
(Wikipedia's "Descent propulsion system" quotes 110 psi at 100 percent thrust
and 11 psi at 10 percent, a ratio of 10 for a thrust ratio of 10) while the
exhaust velocity, the exit Mach number and the plume's angular shape stay put.
Mass flux, and with it density, pressure and shear stress, is then proportional
to thrust and speed, temperature and Mach number are not. That is the scaling
above, and `reference_thrust_n` cancels out of it. What the scaling misses is
that the nozzle boundary layer thickens as the chamber pressure falls -- the
exit Reynolds number drops with thrust, so `delta ~ F^-1/2` -- which widens the
real plume at low throttle. The table is frozen at the boundary layer of its
build thrust; treat the wings of the footprint at deep throttle as approximate.

`source` records the generator and the engine the table was built from.
"""
struct PlumeFieldTable
    name::String
    reference_thrust_n::Float64
    exit_diameter_m::Float64
    limit_speed_mps::Float64
    chamber_temperature_k::Float64
    gamma::Float64
    gas_constant_j_kg_k::Float64
    surface_drag_coefficient::Float64
    height_over_diameter::Vector{Float64}
    radius_over_height::Vector{Float64}
    pressure_hat::Matrix{Float64}
    shear_hat::Matrix{Float64}
    density_hat::Matrix{Float64}
    speed_hat::Matrix{Float64}
    temperature_hat::Matrix{Float64}
    mach::Matrix{Float64}
    footprint_over_height::Vector{Float64}
    source::String
end

plume_field_name(t::PlumeFieldTable) = t.name
@inline plume_field_shear_coefficient(t::PlumeFieldTable, cfg) = t.surface_drag_coefficient

"Smallest height the table was built for; queries below it clamp to it."
@inline _table_min_height(t::PlumeFieldTable)::Float64 = first(t.height_over_diameter) * t.exit_diameter_m

"Index and weight of `v` in the ascending axis `ax`, clamped to its ends."
@inline function _axis_weight(ax::Vector{Float64}, v::Float64)
    n = length(ax)
    n > 1 || return 1, 0.0
    v <= ax[1] && return 1, 0.0
    v >= ax[n] && return n - 1, 1.0
    i = searchsortedlast(ax, v)
    i = clamp(i, 1, n - 1)
    return i, (v - ax[i]) / (ax[i + 1] - ax[i])
end

@inline function _bilinear(m::Matrix{Float64}, i::Int, wi::Float64, j::Int, wj::Float64)::Float64
    a = m[i, j] * (1.0 - wj) + m[i, j + 1] * wj
    b = m[i + 1, j] * (1.0 - wj) + m[i + 1, j + 1] * wj
    return a * (1.0 - wi) + b * wi
end

function plume_gas_state(t::PlumeFieldTable, cfg, thrust_n::Real, height_m::Real, radius_m::Real)::PlumeGasState
    F = Float64(thrust_n)
    h = Float64(height_m)
    r = Float64(radius_m)
    (isfinite(F) && F > 0.0 && isfinite(h) && h >= 0.0 && isfinite(r) && r >= 0.0) || return _PLUME_ZERO_STATE
    h = max(h, _table_min_height(t))
    i, wi = _axis_weight(t.height_over_diameter, h / t.exit_diameter_m)
    j, wj = _axis_weight(t.radius_over_height, r / h)
    scale = F / (h * h)
    u = t.limit_speed_mps
    return PlumeGasState((
        _bilinear(t.pressure_hat, i, wi, j, wj) * scale,
        _bilinear(t.shear_hat, i, wi, j, wj) * scale,
        _bilinear(t.density_hat, i, wi, j, wj) * scale / (u * u),
        _bilinear(t.speed_hat, i, wi, j, wj) * u,
        _bilinear(t.temperature_hat, i, wi, j, wj) * t.chamber_temperature_k,
        _bilinear(t.mach, i, wi, j, wj),
    ))
end

function plume_field_footprint(t::PlumeFieldTable, cfg, thrust_n::Real, height_m::Real)
    F = Float64(thrust_n)
    h = max(Float64(height_m), _table_min_height(t))
    (isfinite(F) && F > 0.0 && isfinite(h)) || return 0.0, Float64(cfg.nozzle_exit_radius_m)
    i, wi = _axis_weight(t.height_over_diameter, h / t.exit_diameter_m)
    p0 = _bilinear(t.pressure_hat, i, wi, 1, 0.0) * F / (h * h)
    y = t.footprint_over_height[i] * (1.0 - wi) + t.footprint_over_height[i + 1] * wi
    return p0, max(y * h, Float64(cfg.nozzle_exit_radius_m))
end

@inline function _plume_shear_at(t::PlumeFieldTable, cfg, thrust_n::Float64, height_m::Float64,
                                 p0::Float64, R::Float64, x::Float64)::Float64
    return plume_gas_state(t, cfg, thrust_n, height_m, x * R).shear_pa
end

@inline function _plume_quadrature_limit(t::PlumeFieldTable, height_m::Float64, R::Float64)::Float64
    R > 0.0 || return 3.0
    return clamp(last(t.radius_over_height) * max(height_m, _table_min_height(t)) / R, 3.0, 12.0)
end

# ---- building the table --------------------------------------------------------------

"""
    build_plume_field_table(nozzle; height_nodes, radius_nodes, max_incidence_deg, source) -> PlumeFieldTable

Evaluate the Simons/Boynton plume of `nozzle` on the ground and tabulate it.

The point source sits at the center of the nozzle exit plane; a ground point at
radius `r` under a nozzle `h` above the ground is a slant distance
`s = sqrt(h^2 + r^2)` away along a streamline making the angle
`theta = atan(r/h)` with the axis, and meets the ground at that same angle from
its normal. The free-stream density there is the plume's mass flux per steradian
divided by the limiting speed and by the near-field denominator of equation (3),
the free-stream temperature follows the isentropic relation to the chamber
density, and `_impinge` turns the result into the wall state.

The radius axis is uniform in the incidence angle out to `max_incidence_deg`,
which puts the nodes where the physics varies rather than where the radius does;
the height axis is logarithmic.
"""
function build_plume_field_table(n::PlumeNozzle; height_nodes::Int=24, radius_nodes::Int=64,
                                 height_over_diameter_min::Float64=0.5,
                                 height_over_diameter_max::Float64=200.0,
                                 max_incidence_deg::Float64=76.0,
                                 source::String="")
    height_nodes >= 2 && radius_nodes >= 2 || throw(ArgumentError("build_plume_field_table needs at least 2 nodes per axis"))
    D = 2.0 * n.exit_radius_m
    F = n.reference_thrust_n
    u_lim = plume_limit_speed(n)
    flux = plume_angular_mass_flux(n, F)
    rho_0 = _chamber_density(n, F)
    g = n.gamma
    r_gas = n.gas_constant_j_kg_k
    t0 = n.chamber_temperature_k
    xs = exp.(range(log(height_over_diameter_min), log(height_over_diameter_max); length=height_nodes))
    ys = tand.(range(0.0, max_incidence_deg; length=radius_nodes))
    P = zeros(height_nodes, radius_nodes); S = similar(P); Rho = similar(P)
    U = similar(P); T = similar(P); M = similar(P)
    for i in 1:height_nodes
        h = xs[i] * D
        hs = _source_height(n, h)
        for j in 1:radius_nodes
            r = ys[j] * h
            slant = sqrt(hs * hs + r * r)
            theta = atan(r, hs)
            cos_inc = hs / slant
            rho_inf = flux(theta) / (u_lim * _near_field_denominator(n, slant))
            t_inf = t0 * (rho_inf / rho_0)^(g - 1.0)
            mach_inf = u_lim / sqrt(max(g * r_gas * t_inf, eps()))
            st = _impinge(g, r_gas, rho_inf, u_lim, cos_inc, mach_inf, t0, n.surface_drag_coefficient)
            P[i, j] = st.pressure_pa * h * h / F
            S[i, j] = st.shear_pa * h * h / F
            Rho[i, j] = st.density_kg_m3 * u_lim * u_lim * h * h / F
            U[i, j] = st.speed_mps / u_lim
            T[i, j] = st.temperature_k / t0
            M[i, j] = st.mach
        end
    end
    foot = [_footprint_fraction(view(P, i, :), ys) for i in 1:height_nodes]
    src = isempty(source) ? "scripts/dev/psi/build_plume_field.jl, engine $(n.name)" : source
    return PlumeFieldTable(n.name, F, D, u_lim, t0, g, r_gas, n.surface_drag_coefficient,
                           collect(xs), collect(ys), P, S, Rho, U, T, M, foot, src)
end

"""
    _footprint_fraction(pressure_row, ys) -> Float64

The `r/h` at which the running integral of the surface pressure over the ground
reaches `1 - 1/e` of its value over the whole tabulated footprint, by the
trapezoidal rule on `p y dy`. It is the definition of footprint radius that
returns `R_p` exactly for a Gaussian, so the table and the analytic field report
the same thing.
"""
function _footprint_fraction(p, ys::Vector{Float64})::Float64
    n = length(ys)
    cum = zeros(n)
    for k in 2:n
        cum[k] = cum[k - 1] + 0.5 * (p[k] * ys[k] + p[k - 1] * ys[k - 1]) * (ys[k] - ys[k - 1])
    end
    total = cum[n]
    total > 0.0 || return ys[n]
    target = (1.0 - exp(-1.0)) * total
    for k in 2:n
        if cum[k] >= target
            w = (target - cum[k - 1]) / max(cum[k] - cum[k - 1], eps())
            return ys[k - 1] + w * (ys[k] - ys[k - 1])
        end
    end
    return ys[n]
end

# ---- storage -------------------------------------------------------------------------

const _PLUME_FIELD_SCHEMA = "spaceagora_plume_field_v1"

"Round to `digits` significant digits so the stored table stays small."
_sig(m, digits::Int) = [round(v; sigdigits=digits) for v in m]

"""
    save_plume_field(path, table; sigdigits=6) -> path

Write a `PlumeFieldTable` as JSON. Six significant digits keeps an Apollo-sized
table under 150 kB while staying far finer than the model's own accuracy.
"""
function save_plume_field(path::AbstractString, t::PlumeFieldTable; sigdigits::Int=6)::String
    doc = Dict{String, Any}(
        "schema" => _PLUME_FIELD_SCHEMA,
        "name" => t.name,
        "reference_thrust_n" => t.reference_thrust_n,
        "exit_diameter_m" => t.exit_diameter_m,
        "limit_speed_mps" => t.limit_speed_mps,
        "chamber_temperature_k" => t.chamber_temperature_k,
        "gamma" => t.gamma,
        "gas_constant_j_kg_k" => t.gas_constant_j_kg_k,
        "surface_drag_coefficient" => t.surface_drag_coefficient,
        "height_over_diameter" => _sig(t.height_over_diameter, sigdigits),
        "radius_over_height" => _sig(t.radius_over_height, sigdigits),
        "footprint_over_height" => _sig(t.footprint_over_height, sigdigits),
        "source" => t.source,
    )
    for (key, m) in ("pressure_hat" => t.pressure_hat, "shear_hat" => t.shear_hat,
                     "density_hat" => t.density_hat, "speed_hat" => t.speed_hat,
                     "temperature_hat" => t.temperature_hat, "mach" => t.mach)
        doc[key] = [_sig(view(m, i, :), sigdigits) for i in 1:size(m, 1)]
    end
    mkpath(dirname(abspath(String(path))))
    open(String(path), "w") do io
        JSON.print(io, doc)
    end
    return String(path)
end

"""
    load_plume_field(path) -> PlumeFieldTable

Read a table written by [`save_plume_field`](@ref), for example
`load_plume_field("data/psi/apollo_lmde.json")`.
"""
function load_plume_field(path::AbstractString)::PlumeFieldTable
    doc = try
        JSON.parsefile(String(path))
    catch err
        throw(ArgumentError("$(path) is not a plume field table: $(sprint(showerror, err))"))
    end
    doc isa AbstractDict && get(doc, "schema", "") == _PLUME_FIELD_SCHEMA ||
        throw(ArgumentError("$(path) is not a plume field table (schema $(repr(get(doc, "schema", nothing)))))."))
    mat(key) = Matrix{Float64}(reduce(vcat, [Float64.(row)' for row in doc[key]]))
    return PlumeFieldTable(String(doc["name"]), Float64(doc["reference_thrust_n"]), Float64(doc["exit_diameter_m"]),
        Float64(doc["limit_speed_mps"]), Float64(doc["chamber_temperature_k"]), Float64(doc["gamma"]),
        Float64(doc["gas_constant_j_kg_k"]), Float64(doc["surface_drag_coefficient"]),
        Float64.(doc["height_over_diameter"]), Float64.(doc["radius_over_height"]),
        mat("pressure_hat"), mat("shear_hat"), mat("density_hat"), mat("speed_hat"),
        mat("temperature_hat"), mat("mach"), Float64.(doc["footprint_over_height"]),
        String(get(doc, "source", "")))
end

end # module PlumeGasField
