# Regolith erosion and cratering regimes under a rocket plume.
#
# The literature on plume-surface interaction does not describe one erosion
# process; it describes several, each with its own onset criterion, and only
# one of them (viscous erosion) is thought to have operated through most of the
# Apollo landings. Metzger (2024a, section 1) lists them:
#
#   "This effort has identified several distinct regimes in which gas jets move
#    the soil. One of these has been called viscous erosion, surface erosion, or
#    simply erosion. It is the process where gas flows across the soil's surface,
#    lifting grains and accelerating them away. The other regimes are
#    collectively called deep cratering because they move soil in bulk with the
#    potential to quickly form deep craters. These regimes include bearing
#    capacity failure, diffused gas eruption, diffusion driven shearing, and
#    diffused gas explosive erosion."
#
# This module implements three of them as separate, individually testable
# closures over a soil description and a gas state on the ground, plus a
# dispatcher that reports which one dominates:
#
#   * viscous erosion, in two forms -- the shear-excess momentum balance the
#     repository already uses (Roberts 1963) and the laminar-sublayer energy-flux
#     law that replaces it in the modern literature (Metzger 2024a/2024b);
#   * diffusion-driven flow, gas percolating into the pores and lifting the
#     surface layer in bulk (Scott and Ko 1968, as summarized in the Lunar
#     Sourcebook section 9.1.8);
#   * bearing capacity failure, the stagnation pressure exceeding the soil's
#     ultimate bearing capacity (Durgunoglu and Mitchell 1975, as used in the
#     Lunar Sourcebook section 9.1.9).
#
# Everything here is a pure function of a gas state, a `RegolithProperties`, the
# local gravity and an `ErosionEnvironment`. Nothing allocates, so the whole set
# can be evaluated per satellite inside a per-step effector loop.
#
# ---------------------------------------------------------------------------
# Sources (author, title, venue, year -- and the equation or table used)
# ---------------------------------------------------------------------------
#
# [LSB]  G. H. Heiken, D. T. Vaniman, B. M. French (eds.), "Lunar Sourcebook: A
#        User's Guide to the Moon", Cambridge University Press, 1991; chapter 9,
#        "Physical Properties of the Lunar Surface", W. D. Carrier III,
#        G. R. Olhoeft, W. Mendell.
#          - section 9.1.2: median particle size 40-130 um, average 70 um.
#          - section 9.1.3: recommended specific gravity 3.1 for lunar soils.
#          - section 9.1.4 / Table 9.4: in situ bulk density 1.50 +/- 0.05 g/cm^3
#            for the top 15 cm, 1.66 +/- 0.05 g/cm^3 for the top 60 cm.
#          - Table 9.5: average in situ porosity 49 percent in the top 30 cm.
#          - Table 9.12: recommended cohesion and friction angle, intercrater
#            areas: 0-15 cm, c = 0.52 kPa (0.44-0.62), phi = 42 deg (41-43);
#            0-30 cm, 0.90 kPa, 46 deg; 30-60 cm, 3.0 kPa, 54 deg; 0-60 cm,
#            1.6 kPa (1.3-1.9), 49 deg (48-51).
#          - section 9.1.8: permeability of the top ~25 cm deduced as
#            1-7 x 10^-12 m^2 from the Surveyor 5 vernier engine firing
#            (Choate et al., 1968); and the statement that rocket-exhaust
#            inflow was studied "to estimate the amount of erosion that could
#            occur when the accumulated gas pressure in the pores exceeded the
#            weight of the overlying soil" (cf. Scott and Ko, 1968) -- the
#            diffusion-driven-flow criterion implemented below.
#          - section 9.1.9: ultimate bearing capacity after Durgunoglu and
#            Mitchell (1975); approximately 6000 kPa for a 1 m footing, and
#            3000-11000 kPa for the Apollo 11 LM footpad (nearly 1 m wide).
#
# [SL00] Y. Shao, H. Lu, "A simple expression for wind erosion threshold
#        friction velocity", Journal of Geophysical Research 105(D17),
#        22437-22443, 2000; equation (22),
#            u_*t^2 = f(Re_*t) (sigma_p g d + gamma / (rho d)),
#        fitted with f(Re_*t) = A_N = 0.0123 and gamma between 1.65e-4 and
#        5e-4 kg/s^2 (see the paragraph following equation (22) and the summary
#        on the paper's last page).
#
# [M24a] P. T. Metzger, "Erosion rate of lunar soil under a landing rocket,
#        part 1: identifying the rate-limiting physics", Icarus, 2024
#        (arXiv:2403.18583); equation (16), the laminar-sublayer energy-flux
#        erosion law, with the mean lift height <D> = 1.5 D84 and the molecular
#        relations tau = rho_s v_s vbar / 6 and E = v_s tau / 2 in section 4;
#        D50 = 77.5 um and D84 = 2.3 D50 for the baseline lunar soil (section 2).
#
# [M24b] P. T. Metzger, "Erosion rate of lunar soil under a landing rocket,
#        part 2: benchmarking and predictions", Icarus, 2024
#        (arXiv:2403.18584); section 2.4, dust blowing first observed in the
#        Apollo 16 landing video at 31.5 m altitude, giving the erosion
#        threshold E_th = 0.123 J/m^2/s; section 3.2, best fit erosion
#        efficiency eps = 0.0029 for the baseline soil (0.0029-0.0042 across the
#        three soil models of Table 1); section 4.2, cohesive energy density
#        alpha_0 = 0.289 J/m^3 at a surface bulk density of 1000 kg/m^3;
#        section 6, total eroded mass 11-26 t for an Apollo landing.
#
# [R63]  L. Roberts, "The action of a hypersonic jet on a dust layer",
#        IAS Paper 63-50, 1963 -- the shear-excess erosion closure the
#        repository's `plume_surface_interaction.jl` currently uses. Metzger
#        (2024a, conclusions) states this closure "is incorrect and should not
#        be used"; it is kept here only so the old and new laws can be compared
#        on the same gas state.
#
# [V73]  A. S. Vesic, "Analysis of ultimate loads of shallow foundations",
#        Journal of the Soil Mechanics and Foundations Division (ASCE) 99(SM1),
#        45-73, 1973; N_gamma = 2 (N_q + 1) tan(phi), and the circular-footing
#        shape factors s_c = 1 + N_q/N_c, s_gamma = 0.6.
#        N_q = exp(pi tan phi) tan^2(45 deg + phi/2) is Reissner (1924) and
#        N_c = (N_q - 1) cot(phi) is Prandtl (1920), both standard.
#
# [W06]  F. M. White, "Viscous Fluid Flow", 3rd ed., McGraw-Hill, 2006;
#        equation (1-36), Sutherland's law for air,
#        mu = 1.458e-6 T^1.5 / (T + 110.4) Pa s. Used here only as a stand-in
#        for the viscosity of the exhaust mixture (see
#        `gas_dynamic_viscosity_pa_s`); this is an ASSUMPTION, not a
#        measurement of Aerozine-50/N2O4 combustion products.
#
# ---------------------------------------------------------------------------
# What is sourced, what is derived, what is assumed
# ---------------------------------------------------------------------------
#
# Sourced, with the citation above each field or function: every default in
# `lunar_mare_regolith`, the Shao-Lu threshold coefficients, Metzger's
# E_th / eps / alpha_0 / <D>, the bearing-capacity factor formulas, the Lunar
# Sourcebook's own bearing-capacity anchor.
#
# Derived in this file, from the sourced framework: the wall-shear form of
# Metzger's energy threshold (`energy_flux_threshold_shear_pa`); the pressure
# diffusion depth from unsteady compressible Darcy flow
# (`pressure_diffusion_depth_m`); the net uplift available to
# diffusion-driven flow after the plume's own downward pressure on the lifting
# plug is subtracted; the Mohr-Coulomb tensile cutoff that resists it; the
# bulk-flow speeds that turn both deep-cratering criteria into mass fluxes.
#
# Assumed, and therefore exposed as documented configuration rather than buried
# literals: `saltation_efficiency` (the repository's existing fitted factor of
# 10 in the Roberts closure, kept only for comparison), the exhaust viscosity
# stand-in, and the `ErosionEnvironment` geometry and timing (footprint radius,
# gas residence time, bearing-failure width), which are properties of the
# scenario rather than of the soil.
module RegolithErosion

export RegolithProperties, lunar_mare_regolith
export ErosionEnvironment, erosion_environment
export ErosionRegimeKind, NoErosion, ViscousErosion, DiffusionDrivenFlowRegime, BearingCapacityFailureRegime
export AbstractErosionRegime, ViscousErosionRoberts, ViscousErosionEnergyFlux
export DiffusionDrivenFlow, BearingCapacityFailure
export erosion_rate, erosion_onset, regime_kind, default_erosion_regimes, regolith_erosion_rate
export shields_threshold_shear_pa, energy_flux_threshold_shear_pa, soil_bearing_capacity_pa
export mean_thermal_speed_mps, mean_lift_height_m, gas_dynamic_viscosity_pa_s, pressure_diffusion_depth_m
export soil_tensile_strength_pa

# ---- soil ----------------------------------------------------------------------------

"""
    RegolithProperties(; ...)

Mechanical, granulometric and transport properties of a soil, as the erosion
regimes of this module need them. The defaults are lunar mare regolith in the
top few centimeters (see [`lunar_mare_regolith`](@ref) for the sources of every
number); construct one with different values for a different site or body.

Granulometry and density (Lunar Sourcebook chapter 9):

- `bulk_density_kg_m3` -- in situ bulk density. 1500 kg/m³ is the Sourcebook's
  best estimate for the top 15 cm (Table 9.4).
- `particle_density_kg_m3` -- mineral density of the grains, 3100 kg/m³ from
  the recommended specific gravity of 3.1 (section 9.1.3).
- `median_diameter_m` -- median grain size `D50`, 70 µm (section 9.1.2; the
  measured range is 40-130 µm).
- `d84_over_d50` -- ratio of the 84th-percentile to the median grain size.
  2.3 for Metzger's baseline lunar soil (2024a, section 2). It only enters the
  mean lift height (see [`mean_lift_height_m`](@ref)).
- `porosity` -- pore volume fraction. 0.52 follows from 1500/3100 by the
  Sourcebook's own relation `n = 1 - rho / (G rho_w)`; Table 9.5 quotes 49
  percent averaged over the top 30 cm, where the soil is denser.

Strength (Lunar Sourcebook Table 9.12, intercrater areas, 0-15 cm):

- `cohesion_pa` -- Mohr-Coulomb cohesion, 520 Pa (range 440-620 Pa).
- `friction_angle_deg` -- Mohr-Coulomb friction angle, 42° (range 41-43°).

Transport (Lunar Sourcebook section 9.1.8):

- `permeability_m2` -- absolute permeability. 3e-12 m² is the middle of the
  1-7 × 10⁻¹² m² deduced from the Surveyor 5 vernier-engine firing over the top
  25 cm (Choate et al., 1968). No direct measurement on returned samples exists.

Erosion-law coefficients:

- `shields_coefficient` and `cohesion_parameter_kg_s2` are Shao and Lu's `A_N`
  and `gamma` (2000, equation 22): 0.0123, and 3e-4 kg/s² out of a fitted range
  of 1.65e-4 to 5e-4. They set the threshold shear stress of
  [`shields_threshold_shear_pa`](@ref).
- `erosion_energy_threshold_w_m2` = 0.123 J/(m² s), the energy-flux erosion
  threshold Metzger (2024b, section 2.4) obtains from the 31.5 m altitude at
  which dust first blows in the Apollo 16 landing video.
- `erosion_efficiency` = 0.0029, the fraction of the energy crossing the lift
  height that does mechanical work, fitted to the Apollo 16 dust opacity
  (Metzger 2024b, section 3.2; 0.0029-0.0042 across the three soil models).
- `cohesive_energy_density_j_m3` = 0.289 J/m³, the van der Waals cohesive
  energy density of the baseline lunar soil at a surface bulk density of
  1000 kg/m³ (Metzger 2024b, section 4.2). It rises with compaction by a law
  whose coefficient is unconstrained, so this is a lower bound for denser soil.
- `lift_height_over_d84` = 1.5, from `<D> = 1.5 D84` (Metzger 2024a, section 4).
- `saltation_efficiency` is NOT sourced. It is the factor of 10 the
  repository's existing Roberts closure uses to stand in for the saltation
  cascade, kept here so the old law can be evaluated on the same gas state as
  the new one. Metzger (2024a, conclusions) holds that the Roberts closure is
  wrong in form, not just in this coefficient.
"""
Base.@kwdef struct RegolithProperties
    bulk_density_kg_m3::Float64 = 1_500.0
    particle_density_kg_m3::Float64 = 3_100.0
    median_diameter_m::Float64 = 70.0e-6
    d84_over_d50::Float64 = 2.3
    porosity::Float64 = 0.52
    cohesion_pa::Float64 = 520.0
    friction_angle_deg::Float64 = 42.0
    permeability_m2::Float64 = 3.0e-12
    shields_coefficient::Float64 = 0.0123
    cohesion_parameter_kg_s2::Float64 = 3.0e-4
    erosion_energy_threshold_w_m2::Float64 = 0.123
    erosion_efficiency::Float64 = 0.0029
    cohesive_energy_density_j_m3::Float64 = 0.289
    lift_height_over_d84::Float64 = 1.5
    saltation_efficiency::Float64 = 10.0
end

"""
    lunar_mare_regolith(; kwargs...) -> RegolithProperties

Lunar mare regolith as the Lunar Sourcebook recommends it for the top few
centimeters, which is the layer a descent plume touches. Identical to
`RegolithProperties()`; the named constructor exists so a scenario records which
soil it flew, and so the sources sit next to the numbers:

| Property | Value | Source |
|---|---|---|
| bulk density | 1500 kg/m³ | Lunar Sourcebook Table 9.4, top 15 cm (1.50 ± 0.05 g/cm³) |
| particle density | 3100 kg/m³ | Lunar Sourcebook section 9.1.3, recommended specific gravity 3.1 |
| median diameter `D50` | 70 µm | Lunar Sourcebook section 9.1.2 (range 40-130 µm) |
| `D84/D50` | 2.3 | Metzger 2024a section 2, baseline lunar soil |
| porosity | 0.52 | `1 - 1500/3100`, the Sourcebook's relation; Table 9.5 gives 49 percent over the top 30 cm |
| cohesion | 520 Pa | Lunar Sourcebook Table 9.12, 0-15 cm, average of 440-620 Pa |
| friction angle | 42° | Lunar Sourcebook Table 9.12, 0-15 cm, average of 41-43° |
| permeability | 3 × 10⁻¹² m² | Lunar Sourcebook section 9.1.8, 1-7 × 10⁻¹² m² from Surveyor 5 (Choate et al. 1968) |

Pass keyword arguments to override any of them, for instance
`lunar_mare_regolith(bulk_density_kg_m3=1660.0, cohesion_pa=1600.0, friction_angle_deg=49.0)`
for the Sourcebook's 0-60 cm column.
"""
lunar_mare_regolith(; kwargs...) = RegolithProperties(; kwargs...)

"""
    mean_lift_height_m(soil) -> Float64

Height a grain must rise before the horizontal gas flow can sweep it away,
`<D> = 1.5 D84` (Metzger 2024a, section 4). With the lunar defaults
`D84 = 2.3 D50 = 161 µm` and `<D> = 241 µm`. This is the length over which the
energy-flux erosion law does its work, and the denominator of that law is the
energy needed to lift a unit volume of soil through it. Metzger's own baseline
soil has `D50 = 77.5 µm` and so `<D> = 267 µm`; the default here uses the Lunar
Sourcebook's 70 µm median instead, an 11 percent difference in `<D>` and in the
rate's gravitational resisting term, which is in any case the smaller of the two
terms in the denominator.
"""
@inline function mean_lift_height_m(soil::RegolithProperties)::Float64
    return soil.lift_height_over_d84 * soil.d84_over_d50 * soil.median_diameter_m
end

"""
    soil_tensile_strength_pa(soil) -> Float64

Mohr-Coulomb tensile cutoff of the soil,

```math
\\sigma_t = \\frac{2 c \\cos\\varphi}{1 + \\sin\\varphi},
```

the isotropic tension at which a Mohr circle of zero confining stress first
touches the failure envelope. Derived here from `cohesion_pa` and
`friction_angle_deg`; with the lunar 0-15 cm values (520 Pa, 42°) it is 463 Pa.
It is the strength a lifting plug of surface soil must break in
[`DiffusionDrivenFlow`](@ref).
"""
@inline function soil_tensile_strength_pa(soil::RegolithProperties)::Float64
    phi = deg2rad(soil.friction_angle_deg)
    s = sin(phi)
    return 2.0 * soil.cohesion_pa * cos(phi) / (1.0 + s)
end

# ---- gas-state helpers ---------------------------------------------------------------

const _SUTHERLAND_C1 = 1.458e-6        # kg / (m s K^0.5), White 2006 eq. (1-36)
const _SUTHERLAND_S_K = 110.4          # K, White 2006 eq. (1-36)

"""
    gas_dynamic_viscosity_pa_s(temperature_k) -> Float64

Dynamic viscosity of the gas at the surface, from Sutherland's law for air
(White, "Viscous Fluid Flow", 3rd ed., 2006, equation 1-36):

```math
\\mu = \\frac{1.458\\times10^{-6}\\, T^{3/2}}{T + 110.4}\\ \\mathrm{Pa\\,s}.
```

**This is an assumption, not a measurement of rocket exhaust.** Air is used as a
stand-in because the viscosity of Aerozine-50/N₂O₄ combustion products (mostly
N₂, H₂O, CO₂, CO and H₂) is not tabulated here; for a mixture of that
composition the true value is within roughly 30 percent of air's at the same
temperature, and the one place it is used -- the pressure diffusion depth --
depends on it only as `1/sqrt(mu)`. It returns 2.67 × 10⁻⁵ Pa s at 500 K and
5.26 × 10⁻⁵ Pa s at 1500 K. Pass an override through
[`ErosionEnvironment`](@ref) when a better value is available.
"""
@inline function gas_dynamic_viscosity_pa_s(temperature_k::Real)::Float64
    T = Float64(temperature_k)
    (isfinite(T) && T > 0.0) || return 0.0
    return _SUTHERLAND_C1 * T * sqrt(T) / (T + _SUTHERLAND_S_K)
end

"""
    mean_thermal_speed_mps(gas) -> Float64

Mean molecular thermal speed of the gas at the surface,

```math
\\bar v = \\sqrt{\\frac{8 R T}{\\pi M}} = \\sqrt{\\frac{8 p}{\\pi \\rho}},
```

the second form following from the ideal-gas law `p = ρ R T / M`, so no
molecular weight has to be supplied: the gas state's pressure and density carry
it. `v̄` is the transport speed of Metzger's molecular-diffusion picture of the
laminar sublayer (2024a, section 4).
"""
@inline function mean_thermal_speed_mps(gas)::Float64
    p = Float64(gas.pressure_pa)
    rho = Float64(gas.density_kg_m3)
    (isfinite(p) && p > 0.0 && isfinite(rho) && rho > 0.0) || return 0.0
    return sqrt(8.0 * p / (pi * rho))
end

"""
    _lift_height_speed_mps(gas, soil) -> Float64

Gas speed at the mean lift height `<D>`, from Metzger's molecular relation
`tau = rho_s v_s vbar / 6` (2024a, section 4) inverted for `v_s`:

```math
v_s = \\frac{6\\,\\tau}{\\rho_s \\bar v}.
```

This is the velocity that appears in the downward energy flux
`E = v_s \\tau / 2`, and it is far below the free-stream wall-jet speed: with
`τ ≈ 1 Pa` under an Apollo LM at 10 m it is of order 10 m/s, not the kilometer
per second of the jet itself. Using the gas state's own `speed_mps` here would
overstate the energy flux by two orders of magnitude.
"""
@inline function _lift_height_speed_mps(gas, soil::RegolithProperties)::Float64
    tau = Float64(gas.shear_pa)
    rho = Float64(gas.density_kg_m3)
    vbar = mean_thermal_speed_mps(gas)
    (isfinite(tau) && tau > 0.0 && rho > 0.0 && vbar > 0.0) || return 0.0
    return 6.0 * tau / (rho * vbar)
end

"""
    _lift_height_energy_flux_w_m2(gas, soil) -> Float64

Downward kinetic-energy flux across the mean lift height,
`E = v_s τ / 2 = 3 τ² / (ρ_s v̄)` (Metzger 2024a, section 4). This is the
quantity the modern erosion law is linear in, and the quantity whose threshold
`E_th` Metzger fixes from the Apollo 16 video.
"""
@inline function _lift_height_energy_flux_w_m2(gas, soil::RegolithProperties)::Float64
    tau = Float64(gas.shear_pa)
    return 0.5 * _lift_height_speed_mps(gas, soil) * max(tau, 0.0)
end

# ---- thresholds ----------------------------------------------------------------------

"""
    shields_threshold_shear_pa(soil, g_m_s2) -> Float64

Threshold wall shear stress for entraining a grain, from the soil's own
properties rather than from a fit. Shao and Lu (2000, equation 22) write the
threshold friction velocity of a loose bed of spherical grains as

```math
u_{*t}^2 = A_N\\left(\\frac{\\rho_p}{\\rho} g d + \\frac{\\gamma}{\\rho d}\\right),
```

which, multiplied through by the gas density, gives the threshold shear stress
free of any gas property at all:

```math
\\tau_t = \\rho u_{*t}^2 = A_N\\left(\\rho_p g d + \\frac{\\gamma}{d}\\right).
```

The first term is the grain's weight, the second its interparticle cohesion.
Shao and Lu fit `A_N = 0.0123` and `γ` in 1.65e-4 to 5e-4 kg/s² to the
threshold data of Iversen and White (1982) for 50-1800 µm grains.

With the lunar defaults (70 µm grains, ρ_p = 3100 kg/m³, g = 1.625 m/s²) the
weight term contributes 0.0043 Pa and the cohesion term 0.0527 Pa: **the lunar
threshold is cohesion-dominated by a factor of 12**, because one sixth of
Earth's gravity suppresses the weight term while the cohesion term is
unchanged. The result is

    tau_t = 0.057 Pa   (0.033 Pa at gamma = 1.65e-4, 0.092 Pa at gamma = 5e-4)

against the **0.15 Pa fitted constant** in `PlumeSurfaceConfig`, which was
chosen to place the erosion onset at 31 m. The derived value is 2.6 times
lower, and 1.6 times below even the upper end of Shao and Lu's `γ` range. Two
honest readings of that gap: Shao and Lu's `γ` is fitted to terrestrial dust
under air, where adsorbed water dominates the cohesion, so extrapolating it to
airless, electrostatically charged lunar fines is an extrapolation; and the
0.15 Pa constant absorbs whatever error the repository's assumed Gaussian shear
law makes in the shear stress itself. See
[`energy_flux_threshold_shear_pa`](@ref) for a second, independent estimate,
which lands much closer to 0.15 Pa.
"""
@inline function shields_threshold_shear_pa(soil::RegolithProperties, g_m_s2::Real)::Float64
    d = soil.median_diameter_m
    g = Float64(g_m_s2)
    (isfinite(d) && d > 0.0 && isfinite(g) && g >= 0.0) || return Inf
    return soil.shields_coefficient * (soil.particle_density_kg_m3 * g * d + soil.cohesion_parameter_kg_s2 / d)
end

"""
    energy_flux_threshold_shear_pa(gas, soil) -> Float64

The same threshold read off Metzger's energy-flux law instead of a Shields
criterion, and therefore an independent check on it. Metzger (2024b, section
2.4) fixes the erosion threshold at `E_th = 0.123 J/(m² s)` from the 31.5 m
altitude at which dust first blows in the Apollo 16 landing video. Because
`E = 3 τ² / (ρ_s v̄)` (2024a, section 4), that energy threshold is a shear
threshold:

```math
\\tau_t = \\sqrt{\\frac{E_{th}\\, \\rho_s \\bar v}{3}}.
```

Unlike [`shields_threshold_shear_pa`](@ref) this depends on the local gas state,
not on the soil alone, which is the more physical statement: the same shear
stress carried by a denser, slower gas moves more soil.

Evaluated at the surface state an 11.5 kN Apollo plume lays down at 10 m
altitude under the repository's Gaussian footprint (p = 168 Pa) with an assumed
500 K, 21.5 g/mol exhaust (ρ = 8.7 × 10⁻⁴ kg/m³, v̄ = 705 m/s) this returns
**0.159 Pa** -- six percent from the 0.15 Pa fitted constant in
`PlumeSurfaceConfig`, and reached from an entirely different observable (Apollo
16 dust opacity rather than Apollo 11 crew reports).

Read that agreement with the two caveats it deserves. Substituting
`v̄ = sqrt(8p/(πρ_s))` makes the threshold `∝ (ρ_s p)^{1/4}`, or `∝ sqrt(p)` for
an ideal gas at fixed temperature, so it is not a soil constant at all: at
31.5 m, where the same footprint gives a tenth the pressure, it falls to
0.05 Pa. And the density and temperature it needs are assumed here, because the
repository's analytic plume field supplies pressure and shear only; a factor of
two in surface density at fixed pressure moves the threshold by 19 percent. The
comparison is therefore suggestive, not a validation, until a plume field that
computes density and temperature is available.
"""
@inline function energy_flux_threshold_shear_pa(gas, soil::RegolithProperties)::Float64
    rho = Float64(gas.density_kg_m3)
    vbar = mean_thermal_speed_mps(gas)
    (isfinite(rho) && rho > 0.0 && vbar > 0.0) || return Inf
    return sqrt(soil.erosion_energy_threshold_w_m2 * rho * vbar / 3.0)
end

"""
    soil_bearing_capacity_pa(soil, g_m_s2, width_m) -> Float64

Static ultimate bearing capacity of the soil under a circular loaded area of
diameter `width_m`. The Lunar Sourcebook (section 9.1.9, after Durgunoglu and
Mitchell 1975) writes it as

```math
q_{ult} = \\tfrac12 \\rho g_m B N_{\\gamma} \\xi_{\\gamma} + c N_c \\xi_c,
```

and this implementation uses the classical bearing-capacity factors:
`N_q = e^{\\pi \\tan\\varphi}\\tan^2(45° + \\varphi/2)` (Reissner 1924),
`N_c = (N_q - 1)\\cot\\varphi` (Prandtl 1920) and
`N_{\\gamma} = 2 (N_q + 1)\\tan\\varphi` (Vesic 1973), with Vesic's
circular-footing shape factors `ξ_c = 1 + N_q/N_c` and `ξ_γ = 0.6`.

Known gap, stated rather than tuned away: with the Sourcebook's own 0-60 cm
soil (c = 1.6 kPa, φ = 49°, ρ = 1660 kg/m³) and `width_m = 1.0` this returns
about 1.3 MPa, whereas the Sourcebook quotes roughly 6 MPa for a 1 m footing
and 3-11 MPa for the Apollo 11 LM footpad. The difference is the factors:
Durgunoglu and Mitchell's wedge-penetration factors, calibrated to
penetrometer data, are several times the classical shallow-footing ones used
here. Treat this function as a lower bound on the soil's strength, which makes
the [`BearingCapacityFailure`](@ref) onset criterion conservative (it triggers
sooner than the Sourcebook numbers would).
"""
function soil_bearing_capacity_pa(soil::RegolithProperties, g_m_s2::Real, width_m::Real)::Float64
    g = Float64(g_m_s2)
    B = Float64(width_m)
    phi = deg2rad(soil.friction_angle_deg)
    (isfinite(g) && g >= 0.0 && isfinite(B) && B >= 0.0) || return Inf
    (phi > 0.0 && phi < 0.5 * pi) || return Inf
    t = tan(phi)
    n_q = exp(pi * t) * tan(0.25 * pi + 0.5 * phi)^2
    n_c = (n_q - 1.0) / t
    n_gamma = 2.0 * (n_q + 1.0) * t
    xi_c = 1.0 + n_q / n_c
    xi_gamma = 0.6
    return soil.cohesion_pa * n_c * xi_c + 0.5 * soil.bulk_density_kg_m3 * g * B * n_gamma * xi_gamma
end

"""
    pressure_diffusion_depth_m(gas, soil, residence_time_s, viscosity_pa_s) -> Float64

Depth to which the plume's gas has percolated into the pores after
`residence_time_s` of loading. Isothermal Darcy flow of a compressible gas
through a porous medium obeys

```math
n \\frac{\\partial p}{\\partial t} = \\frac{k}{\\mu}\\nabla\\cdot(p\\,\\nabla p),
```

a nonlinear diffusion with pressure diffusivity `D = k p / (μ n)`; the front
therefore advances as `δ = sqrt(D t)`. Derived here from Darcy's law (the Lunar
Sourcebook states it in section 9.1.8) and mass conservation; the numerical
coefficient of order one that a similarity solution would supply is dropped.

With the lunar permeability (3 × 10⁻¹² m²), porosity 0.52, the 500 K air
viscosity of 2.67 × 10⁻⁵ Pa s and one second of loading, `δ` is 6 mm under the
168 Pa an Apollo plume lays down at 10 m, 3.7 cm under the 6.5 kPa it reaches at
contact on the approach throttle, and 7.4 cm under the 25 kPa of a full-throttle
descent engine at contact. Because `δ ∝ sqrt(t)`, a ten times longer hover
deepens the front only threefold.
"""
@inline function pressure_diffusion_depth_m(gas, soil::RegolithProperties, residence_time_s::Real,
                                            viscosity_pa_s::Real)::Float64
    p = Float64(gas.pressure_pa)
    mu = Float64(viscosity_pa_s)
    t = Float64(residence_time_s)
    n = soil.porosity
    k = soil.permeability_m2
    (isfinite(p) && p > 0.0 && isfinite(t) && t > 0.0) || return 0.0
    (isfinite(mu) && mu > 0.0 && n > 0.0 && k > 0.0) || return 0.0
    return sqrt(k * p * t / (mu * n))
end

# ---- scenario geometry and timing ----------------------------------------------------

"""
    ErosionEnvironment(; footprint_radius_m, residence_time_s, bearing_width_m, gas_viscosity_pa_s)

The parts of the problem that belong to the scenario rather than to the soil or
to the gas state at one point. None of these are sourced physical constants;
each is a documented configuration field.

- `footprint_radius_m` -- radius of the plume's surface pressure footprint. It
  sets the lateral scale over which the surface pressure varies, which is what
  makes [`DiffusionDrivenFlow`](@ref) possible at all (see that regime's
  docstring). Default 1 m, which is the order of an Apollo LM's footprint at
  touchdown; a caller inside the effector loop passes the real radius.
- `residence_time_s` -- how long the plume has been loading the patch of ground
  being evaluated, which sets the pressure diffusion depth. Default 1 s: the
  surface pressure at a fixed radius changes on the timescale over which the
  vehicle's height changes appreciably, one second at the meter-per-second
  descent rates of the last few meters.
- `bearing_width_m` -- diameter of the loaded area used for the bearing-capacity
  factors. Default 1 m; `NaN` means "use the footprint diameter".
- `gas_viscosity_pa_s` -- dynamic viscosity of the exhaust at the surface.
  `NaN` (the default) means "compute it from the gas temperature with
  [`gas_dynamic_viscosity_pa_s`](@ref)", which is itself a documented
  assumption.
"""
Base.@kwdef struct ErosionEnvironment
    footprint_radius_m::Float64 = 1.0
    residence_time_s::Float64 = 1.0
    bearing_width_m::Float64 = 1.0
    gas_viscosity_pa_s::Float64 = NaN
end

"""
    erosion_environment(; kwargs...) -> ErosionEnvironment

Keyword constructor for [`ErosionEnvironment`](@ref), for symmetry with
[`lunar_mare_regolith`](@ref).
"""
erosion_environment(; kwargs...) = ErosionEnvironment(; kwargs...)

const DEFAULT_EROSION_ENVIRONMENT = ErosionEnvironment()

@inline function _viscosity(gas, env::ErosionEnvironment)::Float64
    mu = env.gas_viscosity_pa_s
    return isfinite(mu) && mu > 0.0 ? mu : gas_dynamic_viscosity_pa_s(gas.temperature_k)
end

@inline function _bearing_width(env::ErosionEnvironment)::Float64
    B = env.bearing_width_m
    return isfinite(B) && B > 0.0 ? B : 2.0 * env.footprint_radius_m
end

# ---- regimes -------------------------------------------------------------------------

"""
    AbstractErosionRegime

Supertype of the erosion and cratering regimes. Each concrete subtype is a
singleton carrying no data, so a tuple of them dispatches without allocating,
and each implements [`erosion_onset`](@ref) and [`erosion_rate`](@ref).
"""
abstract type AbstractErosionRegime end

"""
    ViscousErosionRoberts()

Viscous erosion in the form Roberts (1963) gave it and
`plume_surface_interaction.jl` currently uses: the wall shear stress in excess
of a threshold is the force available to accelerate grains, and a mass flux
`ṁ` leaving at speed `v` carries `ṁ v` of momentum, so

```math
\\dot m = \\frac{\\eta\\,\\max(\\tau - \\tau_t,\\ 0)}{v}.
```

The threshold `τ_t` is [`shields_threshold_shear_pa`](@ref), derived from the
soil, replacing the fitted 0.15 Pa. `η = saltation_efficiency` is unchanged and
still unsourced. The grain speed is taken as the gas speed at the lift height,
[`_lift_height_speed_mps`](@ref), rather than the free-stream speed.

**This closure is retained for comparison only.** Metzger (2024a, conclusions)
states plainly: "Roberts' hypothesis for lunar soil erosion, that it is a
shearing process linear with gas shear stress and governed both by shear
strength of the soil and by acceleration of the eroded particles reducing shear
stress of the gas, is incorrect and should not be used." Prefer
[`ViscousErosionEnergyFlux`](@ref); this type exists so a study can quantify
what changes when the old law is swapped out.
"""
struct ViscousErosionRoberts <: AbstractErosionRegime end

"""
    ViscousErosionEnergyFlux()

Viscous erosion as the laminar-sublayer energy-flux process Metzger (2024a,
equation 16) argues it actually is:

```math
\\dot m = \\rho_b\\,\\frac{\\varepsilon\\,(E - E_{th})}{\\rho_b\\, g \\langle D\\rangle + \\alpha},
\\qquad E > E_{th},
```

and zero otherwise. The numerator is the rate at which kinetic energy diffuses
downward across the mean lift height `<D> = 1.5 D84`, of which a fraction
`ε` does mechanical work; the denominator is the energy a unit volume of soil
costs to disassemble -- the potential energy to raise it through `<D>` plus the
cohesive energy density `α` that binds it. `E = 3 τ² / (ρ_s v̄)` follows from
the molecular relation `τ = ρ_s v_s v̄ / 6` and `E = v_s τ / 2` (2024a, section
4), with `v̄` from [`mean_thermal_speed_mps`](@ref).

Note what this law demands of the gas state: it is linear in `τ²/ρ_s`, so it
needs the gas density and temperature at the surface, not only the pressure.
Fed a density synthesized from an analytic pressure footprint at a fixed
temperature (`ρ ∝ p`), it degenerates -- `E ∝ τ²/ρ ∝ p`, whose integral over a
momentum-conserving footprint is the thrust, so the total erosion rate comes out
independent of altitude, which is wrong. Use it with a plume field that computes
density and temperature.

Every coefficient is sourced: `ε = 0.0029` and `E_th = 0.123 J/(m² s)` from the
Apollo 16 landing video (Metzger 2024b, sections 2.4 and 3.2), `α = 0.289 J/m³`
from the van der Waals cohesion of the baseline soil (2024b, section 4.2), and
`<D>` from `D50 = 77.5 µm`, `D84 = 2.3 D50` (2024a, section 2). Note that `α`
is quoted at a surface bulk density of 1000 kg/m³ and rises with compaction by
a law whose coefficient Metzger treats as the dominant remaining uncertainty,
so this rate is an upper bound for denser soil.
"""
struct ViscousErosionEnergyFlux <: AbstractErosionRegime end

"""
    DiffusionDrivenFlow()

Gas percolating into the pores and lifting the surface layer in bulk. The Lunar
Sourcebook (section 9.1.8) states the criterion directly: rocket-exhaust inflow
was studied "to estimate the amount of erosion that could occur when the
accumulated gas pressure in the pores exceeded the weight of the overlying
soil" (cf. Scott and Ko, 1968).

Taken literally that criterion never fires, and saying why is the whole model.
Under a laterally uniform footprint the same pressure that pressurizes the pores
presses down on the plug of soil above them, so the net uplift is zero. Uplift
exists only where the surface pressure varies laterally faster than the pore
pressure can follow: the pore pressure under a plug is the surface pressure
smoothed over the diffusion length `δ`, so over a footprint whose pressure
falls off on a scale `R_p` the net uplift available is of order

```math
\\Delta p = p \\frac{\\delta}{R_p},
\\qquad \\delta = \\sqrt{\\frac{k\\,p\\,t}{\\mu\\,n}},
```

derived here from unsteady Darcy flow ([`pressure_diffusion_depth_m`](@ref))
and a first-order expansion of the surface pressure. The plug lifts when that
uplift exceeds its own weight plus the tensile strength holding it to the bed:

```math
\\Delta p > \\rho_b g \\delta + \\sigma_t,
\\qquad \\sigma_t = \\frac{2 c \\cos\\varphi}{1 + \\sin\\varphi}.
```

Once it lifts, a layer `δ` deep is removed in the time it took to pressurize,
so the mass flux is `ṁ = ρ_b δ / t`. That rate is a derived order of magnitude,
not a measured one; no published measurement of a diffusion-driven-flow mass
flux on lunar soil is used here.

With the lunar defaults, the repository's Gaussian footprint and the 11.5 kN
the Apollo LM carried through its final approach, this regime never fires: at
contact the uplift reaches 325 Pa against 554 Pa of weight and tensile strength.
It fires at contact above about 20 kN, and a full-throttle 45 kN descent engine
crosses the criterion from roughly 2 m down. So the model says a nominal Apollo
landing, which throttles deep for touchdown, saw no diffusion-driven flow, while
a lander that arrives on a harder throttle would. That is the behavior Metzger
(2024a, section 1) reports: deep cratering "generally did not occur during lunar
landings due to the high mechanical strength and low permeability of the soil
and because rocket exhaust in vacuum is not collimated into a jet", yet "during
the final moments of the Apollo landings when the engine nozzle was close to the
surface, an abrupt blast of soil sometimes occurred". The model was not tuned to
produce that; it follows from the Sourcebook's permeability and Table 9.12's
strength.
"""
struct DiffusionDrivenFlow <: AbstractErosionRegime end

"""
    BearingCapacityFailure()

The stagnation pressure exceeds the soil's ultimate bearing capacity and the
surface fails structurally, digging a deep crater rather than being scoured
grain by grain (Metzger et al., "Jet-induced cratering of a granular surface
with application to lunar spaceports", Journal of Aerospace Engineering 22(1),
24-32, 2009, which names bearing capacity failure as one of the cratering
mechanisms; the strength model is the Lunar Sourcebook's, section 9.1.9).

Onset is `p > q_ult(B)` with `q_ult` from [`soil_bearing_capacity_pa`](@ref)
and `B` the loaded width from the [`ErosionEnvironment`](@ref). The mass flux
once it fails is derived from a plug momentum balance -- the excess pressure
accelerates the failing soil, which reaches `v = sqrt(2 (p - q_{ult}) / ρ_b)`,
carrying

```math
\\dot m = \\rho_b v = \\sqrt{2 \\rho_b (p - q_{ult})}.
```

This is a derived upper bound, not a measured rate.

With the lunar defaults and an Apollo LM this never fires: the stagnation
pressure reaches about 25 kPa at contact against a bearing capacity of several
hundred kPa even on the conservative factors used here. Metzger (2024a, section
1) reaches the same conclusion for the Apollo landings while noting "it is
possible that larger lunar landers could dig deep craters" -- in this model that
takes about 450 kN through a 1.5 m footprint, ten times the Apollo descent
engine at full throttle.
"""
struct BearingCapacityFailure <: AbstractErosionRegime end

"""
    ErosionRegimeKind

Compact tag for the regime a dispatcher found dominant. It is an `Int8`-backed
enumeration so the dispatcher's return value is `isbits` and costs no
allocation in a per-step loop. `ViscousErosion` covers both viscous closures;
they are alternative laws for the same physical regime, not different regimes.
"""
@enum ErosionRegimeKind::Int8 begin
    NoErosion = 0
    ViscousErosion = 1
    DiffusionDrivenFlowRegime = 2
    BearingCapacityFailureRegime = 3
end

@doc """
    NoErosion

[`ErosionRegimeKind`](@ref) tag reported when no regime is above its onset.
""" NoErosion

@doc """
    ViscousErosion

[`ErosionRegimeKind`](@ref) tag of the viscous erosion regime, whichever of the
two closures ([`ViscousErosionRoberts`](@ref), [`ViscousErosionEnergyFlux`](@ref))
produced the rate.
""" ViscousErosion

@doc """
    DiffusionDrivenFlowRegime

[`ErosionRegimeKind`](@ref) tag of [`DiffusionDrivenFlow`](@ref). It is spelled
differently from the regime type so the tag and the singleton can both be
exported.
""" DiffusionDrivenFlowRegime

@doc """
    BearingCapacityFailureRegime

[`ErosionRegimeKind`](@ref) tag of [`BearingCapacityFailure`](@ref). It is
spelled differently from the regime type so the tag and the singleton can both
be exported.
""" BearingCapacityFailureRegime

"""
    regime_kind(regime) -> ErosionRegimeKind

The [`ErosionRegimeKind`](@ref) tag of a regime singleton.
"""
@inline regime_kind(::ViscousErosionRoberts) = ViscousErosion
@inline regime_kind(::ViscousErosionEnergyFlux) = ViscousErosion
@inline regime_kind(::DiffusionDrivenFlow) = DiffusionDrivenFlowRegime
@inline regime_kind(::BearingCapacityFailure) = BearingCapacityFailureRegime

# ---- rates ---------------------------------------------------------------------------

@inline function _gas_is_usable(gas)::Bool
    p = Float64(gas.pressure_pa)
    tau = Float64(gas.shear_pa)
    rho = Float64(gas.density_kg_m3)
    T = Float64(gas.temperature_k)
    return isfinite(p) && p >= 0.0 && isfinite(tau) && tau >= 0.0 &&
           isfinite(rho) && rho >= 0.0 && isfinite(T) && T >= 0.0
end

"""
    erosion_rate(regime, gas, soil, g_m_s2, env=ErosionEnvironment()) -> Float64

Mass erosion rate in kg/(m² s) that `regime` produces at one point on the
surface, given the gas state there, the soil, the local gravity and the
scenario's [`ErosionEnvironment`](@ref). Identically zero below the regime's
onset (see [`erosion_onset`](@ref)), and zero for a non-finite or empty gas
state.

`gas` is any object with the fields of the shared plume gas state:
`pressure_pa`, `shear_pa`, `density_kg_m3`, `speed_mps`, `temperature_k` and
`mach`. Pure, non-allocating, and safe to call from a right-hand-side
evaluation.
"""
function erosion_rate(regime::AbstractErosionRegime, gas, soil::RegolithProperties, g_m_s2::Real,
                      env::ErosionEnvironment=DEFAULT_EROSION_ENVIRONMENT)::Float64
    return _erosion_rate(regime, gas, soil, Float64(g_m_s2), env)
end

function _erosion_rate(::ViscousErosionRoberts, gas, soil::RegolithProperties, g::Float64,
                       env::ErosionEnvironment)::Float64
    _gas_is_usable(gas) || return 0.0
    tau = Float64(gas.shear_pa)
    tau_t = shields_threshold_shear_pa(soil, g)
    excess = tau - tau_t
    excess > 0.0 || return 0.0
    v = _lift_height_speed_mps(gas, soil)
    v > 0.0 || return 0.0
    return soil.saltation_efficiency * excess / v
end

function _erosion_rate(::ViscousErosionEnergyFlux, gas, soil::RegolithProperties, g::Float64,
                       env::ErosionEnvironment)::Float64
    _gas_is_usable(gas) || return 0.0
    E = _lift_height_energy_flux_w_m2(gas, soil)
    excess = E - soil.erosion_energy_threshold_w_m2
    excess > 0.0 || return 0.0
    resist = soil.bulk_density_kg_m3 * g * mean_lift_height_m(soil) + soil.cohesive_energy_density_j_m3
    resist > 0.0 || return 0.0
    return soil.bulk_density_kg_m3 * soil.erosion_efficiency * excess / resist
end

"""
    _diffusion_uplift(gas, soil, g, env) -> (uplift_pa, resistance_pa, depth_m)

The three quantities [`DiffusionDrivenFlow`](@ref)'s criterion compares: the net
uplift the lateral pressure gradient makes available, the weight plus tensile
strength of the plug that resists it, and the depth of the plug.
"""
@inline function _diffusion_uplift(gas, soil::RegolithProperties, g::Float64, env::ErosionEnvironment)
    mu = _viscosity(gas, env)
    delta = pressure_diffusion_depth_m(gas, soil, env.residence_time_s, mu)
    R = env.footprint_radius_m
    (delta > 0.0 && isfinite(R) && R > 0.0) || return (0.0, Inf, 0.0)
    uplift = Float64(gas.pressure_pa) * min(delta / R, 1.0)
    resistance = soil.bulk_density_kg_m3 * g * delta + soil_tensile_strength_pa(soil)
    return (uplift, resistance, delta)
end

function _erosion_rate(::DiffusionDrivenFlow, gas, soil::RegolithProperties, g::Float64,
                       env::ErosionEnvironment)::Float64
    _gas_is_usable(gas) || return 0.0
    uplift, resistance, delta = _diffusion_uplift(gas, soil, g, env)
    (uplift > resistance && delta > 0.0) || return 0.0
    t = env.residence_time_s
    (isfinite(t) && t > 0.0) || return 0.0
    return soil.bulk_density_kg_m3 * delta / t
end

function _erosion_rate(::BearingCapacityFailure, gas, soil::RegolithProperties, g::Float64,
                       env::ErosionEnvironment)::Float64
    _gas_is_usable(gas) || return 0.0
    q_ult = soil_bearing_capacity_pa(soil, g, _bearing_width(env))
    excess = Float64(gas.pressure_pa) - q_ult
    excess > 0.0 || return 0.0
    return sqrt(2.0 * soil.bulk_density_kg_m3 * excess)
end

"""
    erosion_onset(regime, gas, soil, g_m_s2, env=ErosionEnvironment()) -> Bool

Whether `regime`'s onset criterion is met at this gas state. Every regime
satisfies `erosion_onset(...) == (erosion_rate(...) > 0)`; the predicate exists
so a caller can ask which regimes are live without evaluating their rates, and
so the criteria can be read and tested on their own.

The criteria are: wall shear above the soil's threshold for the viscous
regimes ([`shields_threshold_shear_pa`](@ref) for
[`ViscousErosionRoberts`](@ref), the energy threshold `E > E_th` for
[`ViscousErosionEnergyFlux`](@ref)); pore-pressure uplift above the overburden
weight plus tensile strength for [`DiffusionDrivenFlow`](@ref); stagnation
pressure above the ultimate bearing capacity for
[`BearingCapacityFailure`](@ref).
"""
function erosion_onset(regime::AbstractErosionRegime, gas, soil::RegolithProperties, g_m_s2::Real,
                       env::ErosionEnvironment=DEFAULT_EROSION_ENVIRONMENT)::Bool
    return _erosion_onset(regime, gas, soil, Float64(g_m_s2), env)
end

function _erosion_onset(::ViscousErosionRoberts, gas, soil::RegolithProperties, g::Float64,
                        env::ErosionEnvironment)::Bool
    _gas_is_usable(gas) || return false
    return Float64(gas.shear_pa) > shields_threshold_shear_pa(soil, g) && _lift_height_speed_mps(gas, soil) > 0.0
end

function _erosion_onset(::ViscousErosionEnergyFlux, gas, soil::RegolithProperties, g::Float64,
                        env::ErosionEnvironment)::Bool
    _gas_is_usable(gas) || return false
    return _lift_height_energy_flux_w_m2(gas, soil) > soil.erosion_energy_threshold_w_m2
end

function _erosion_onset(::DiffusionDrivenFlow, gas, soil::RegolithProperties, g::Float64,
                        env::ErosionEnvironment)::Bool
    _gas_is_usable(gas) || return false
    uplift, resistance, delta = _diffusion_uplift(gas, soil, g, env)
    return uplift > resistance && delta > 0.0 && env.residence_time_s > 0.0
end

function _erosion_onset(::BearingCapacityFailure, gas, soil::RegolithProperties, g::Float64,
                        env::ErosionEnvironment)::Bool
    _gas_is_usable(gas) || return false
    return Float64(gas.pressure_pa) > soil_bearing_capacity_pa(soil, g, _bearing_width(env))
end

# ---- dispatcher ----------------------------------------------------------------------

"""
    default_erosion_regimes() -> Tuple

The regimes a plume-surface model should evaluate by default:
[`ViscousErosionEnergyFlux`](@ref), [`DiffusionDrivenFlow`](@ref) and
[`BearingCapacityFailure`](@ref).

[`ViscousErosionRoberts`](@ref) is deliberately left out. It is an alternative
law for the same physical regime as the energy-flux closure, so including both
would count the same mass twice, and Metzger (2024a) holds that it is the wrong
law; pass it explicitly when comparing the two.
"""
@inline default_erosion_regimes() = (ViscousErosionEnergyFlux(), DiffusionDrivenFlow(), BearingCapacityFailure())

"""
    regolith_erosion_rate(regimes, gas, soil, g_m_s2, env=ErosionEnvironment())

Evaluate every regime in `regimes` at one point on the surface and report the
total mass flux and which regime dominates it. Returns the `isbits` named tuple

```julia
(rate_kg_m2_s, dominant::ErosionRegimeKind, dominant_rate_kg_m2_s, active_count::Int)
```

`dominant` is [`NoErosion`](@ref ErosionRegimeKind) when nothing is above its
onset. `regimes` is a tuple of regime singletons, which Julia unrolls at compile
time, so the whole call is allocation-free and callable from a per-step
effector loop; `default_erosion_regimes()` is the usual argument.

The total is the sum of the regimes' rates: they are distinct mechanisms rather
than alternative descriptions of one flux, except for the two viscous closures,
which must not both appear in `regimes` (see
[`default_erosion_regimes`](@ref)). Which regime comes out dominant is a model
output, not a foregone conclusion -- under a lander large enough to drive all
three, whether the deep-cratering flux beats the viscous one depends on the gas
state, and on the synthesized states available today viscous erosion still wins.
Report `dominant` rather than assuming it.

```julia
soil = lunar_mare_regolith()
gas  = (pressure_pa=170.0, shear_pa=1.0, density_kg_m3=8.7e-4,
        speed_mps=1000.0, temperature_k=500.0, mach=3.0)
env  = erosion_environment(footprint_radius_m=4.66)
out  = regolith_erosion_rate(default_erosion_regimes(), gas, soil, 1.625, env)
out.dominant === ViscousErosion
```
"""
function regolith_erosion_rate(regimes::Tuple, gas, soil::RegolithProperties, g_m_s2::Real,
                               env::ErosionEnvironment=DEFAULT_EROSION_ENVIRONMENT)
    g = Float64(g_m_s2)
    total = 0.0
    best = 0.0
    kind = NoErosion
    active = 0
    @inbounds for regime in regimes
        r = _erosion_rate(regime, gas, soil, g, env)
        r > 0.0 || continue
        total += r
        active += 1
        if r > best
            best = r
            kind = regime_kind(regime)
        end
    end
    return (rate_kg_m2_s=total, dominant=kind, dominant_rate_kg_m2_s=best, active_count=active)
end

"""
    regolith_erosion_rate(gas, soil, g_m_s2, env=ErosionEnvironment())

Convenience method over [`default_erosion_regimes`](@ref).
"""
@inline function regolith_erosion_rate(gas, soil::RegolithProperties, g_m_s2::Real,
                                       env::ErosionEnvironment=DEFAULT_EROSION_ENVIRONMENT)
    return regolith_erosion_rate(default_erosion_regimes(), gas, soil, g_m_s2, env)
end

end # module RegolithErosion
