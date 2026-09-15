---
id: module.environment
label: EnvironmentModels
kind: module
source:
  file: src/environment/physical_models.jl
  symbol: EnvironmentModels
inputs: []
outputs:
- id: api
  type: Module
  units: n/a
  description: The exported environment surface — the density-model structs `NoAtmosphereModel`,
    `ExponentialAtmosphereModel`, `PiecewiseExponentialAtmosphereModel`, `TabulatedFlightAtmosphereModel`,
    `TimeTabulatedAtmosphereModel`, `NRLMSISE00AtmosphereModel`, `GRAMAtmosphereModel`,
    `GRAMAtmosphereModelSurrogate`, `ConstantDensityModel`, together with `getDensity`,
    `getDensityBatch!`, `init_nrlmsise_space_indices!`, `precompute_gram_static_grids!`
    and `clear_gram_static_grid_cache!`.
tags:
- module
charts:
- master
origin: agent
---

# EnvironmentModels

## Purpose
`EnvironmentModels` is the planetary-environment layer: it supplies the atmospheric state
(density, temperature, wind), the gravity field models, the planet constant tables and the
ephemeris/frame machinery that every dynamic effector queries. The module body in
`src/environment/physical_models.jl` includes `atmosphere/density_models.jl` and exports the
`getDensity` interface; the sibling files `gravity/gravity_models.jl`,
`ephemerides/simple_ephemerides.jl`, `ephemerides/planets.jl` and `ephemerides/planet_shapes.jl`
supply gravity, epoch handling and body-fixed frame rotations. `gravity/gravity_effectors.jl`
and `aerodynamics/aerodynamic_effectors.jl` are thin re-export shims that hand the same
symbols back through the environment namespace.

## Theory & Math
Atmospheric density uses a barometric (isothermal-layer) law:

$$
\rho(h) = \rho_{\text{ref}}\,\exp\!\left(\frac{h_{\text{ref}} - h}{H}\right)
$$

where $\rho$ is mass density [kg/m³], $h$ geometric altitude above the reference ellipsoid [m],
$\rho_{\text{ref}}$ the reference density [kg/m³] at reference altitude $h_{\text{ref}}$ [m],
and $H$ the scale height [m].

Gravity is a point mass plus the $J_2$ zonal term, evaluated in planet-fixed coordinates:

$$
\mathbf{g} = -\frac{\mu}{r^2}\hat{\mathbf{r}}
+ \frac{3}{2}\,\frac{J_2\,\mu\,R_e^2}{r^4}
\begin{bmatrix}
\frac{x}{r}\left(\frac{5z^2}{r^2}-1\right)\\[2pt]
\frac{y}{r}\left(\frac{5z^2}{r^2}-1\right)\\[2pt]
\frac{z}{r}\left(\frac{5z^2}{r^2}-3\right)
\end{bmatrix}
$$

with $\mu$ the gravitational parameter [m³/s²], $r = \lVert(x,y,z)\rVert$ the planet-centred
radius [m], $R_e$ the equatorial radius [m] (`planet.Rp_e`), and $J_2$ the second zonal
harmonic [dimensionless]. The induced gravity-gradient torque on an extended body is

$$
\boldsymbol{\tau}_{gg} = \frac{3\mu}{r^3}\,\hat{\mathbf{r}}_b \times \left(\mathbf{J}\,\hat{\mathbf{r}}_b\right)
$$

with $\hat{\mathbf{r}}_b$ the unit nadir vector in body axes and $\mathbf{J}$ the inertia
tensor [kg·m²], giving $\boldsymbol{\tau}_{gg}$ in N·m.

Earth prime-meridian orientation under the simple ephemeris uses IAU-82 GMST:

$$
\theta_{\text{GMST}} = \mathrm{mod}_{2\pi}\!\left(
\frac{2\pi}{86400}\;\mathrm{rem}\!\big(P(T_u),\,86400\big)\right),\quad
T_u = \frac{t_{\text{UT1}}}{36525\cdot 86400}
$$

where $P$ is the polynomial with coefficients `67310.54841, 3.164400184812866e9, 0.093104, -6.2e-6`
in seconds and $t_{\text{UT1}}$ is seconds past the J2000 epoch.

## Model & Assumptions
- Exponential and piecewise-exponential atmospheres carry a *constant* placeholder
  `temperature_k` (default 200.0 K) and return an identically zero wind vector; only the
  tabulated, GRAM and NRLMSISE-00 paths return real thermospheric temperature and winds.
- `valid_min_altitude_m` and `valid_max_altitude_m` on `ExponentialAtmosphereModel` are
  advisory only — `_exponential_density` extrapolates the same exponential outside them.
- The $J_2$ coefficient is normalised to the *equatorial* radius `planet.Rp_e`, matching
  `j2_secular_rates` and the spherical-harmonics path, not to a mean radius.
- The simple ephemeris treats leap-second-free UTC seconds past J2000 as UT1. The source
  records the resulting error as $|UT1 - UTC| < 0.9$ s, under $4\times10^{-3}$ deg of Earth
  rotation.
- With J2000 as the internal inertial frame, a single spin-axis rotation is taken as the
  complete J2000→planet-centred-planet-fixed transform: precession, nutation and polar motion
  are omitted from `SimpleEphemeridesModel`.

## Design & Implementation
`getDensity` is a multiple-dispatch interface keyed on the density-model type, with a uniform
signature `(model, h, lat, lon, el_time, wind, p)` returning `Tuple{Float64, Float64, SVector{3,Float64}}`
= (density, temperature, wind). Models with no epoch dependence define the six-argument method
and a one-line `@inline` seven-argument forwarder, so the hot solver path never pays for the
unused `ODEParams`. `PiecewiseExponentialAtmosphereModel` selects its layer with
`searchsortedlast(model.h_breaks_m, h)` clamped into `1:length(model.ρ_refs)`, so an altitude
below the first break reuses the lowest layer instead of erroring. `getDensityBatch!` writes
into caller-supplied `rhos`, `Ts`, `winds`, `lats`, `lons` vectors and validates that all
lengths match `hs` before any evaluation. On the gravity side,
`_inverse_squared_j2_gravity_accel` is `@inline` over `SVector{3,Float64}` and
`_gravity_gradient_torque_body` short-circuits to the zero vector whenever
`model.gravity_gradient` is false, `orientation_sim` is off, the index `i` is out of range, or
the state carries no `:q` field. `planet_frame_lpi` dispatches on the ephemeris model:
`SpiceEphemeridesModel` calls SPICE `pxform` under `SPICE_LOCK` with `ITRF93` for Earth and an
`IAU_EARTH` fallback, while `SimpleEphemeridesModel` builds `_rotation_about_spin_axis(θ)`
analytically. `ephemerides_cache_key` encodes `NaN` (the planet-true prime-meridian sentinel)
as `typemin(Int64)` because `round` would throw on it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `api` | Module | n/a | — | The exported environment surface — the density-model structs `NoAtmosphereModel`, `ExponentialAtmosphereModel`, `PiecewiseExponentialAtmosphereModel`, `TabulatedFlightAtmosphereModel`, `TimeTabulatedAtmosphereModel`, `NRLMSISE00AtmosphereModel`, `GRAMAtmosphereModel`, `GRAMAtmosphereModelSurrogate`, `ConstantDensityModel`, together with `getDensity`, `getDensityBatch!`, `init_nrlmsise_space_indices!`, `precompute_gram_static_grids!` and `clear_gram_static_grid_cache!`. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `api` → [[envana.env_aerodynamic_effectors_aerodynamiceffectors|AerodynamicEffectors]] · `module_api` · call · `src/environment/aerodynamics/aerodynamic_effectors.jl`
- `api` → [[envana.env_ephemerides_models_ephemeridesmodels|EphemeridesModels]] · `module_api` · call · `src/environment/ephemerides/ephemerides_models.jl`
- `api` → [[envana.env_gravity_effectors_gravityeffectors|GravityEffectors]] · `module_api` · call · `src/environment/gravity/gravity_effectors.jl`
- `api` → [[envana.env_gravity_models_inversesquaredj2gravitymodel|InverseSquaredJ2GravityModel]] · `module_api` · call · `src/environment/gravity/gravity_models.jl`
- `api` → [[envana.env_physical_models_environmentmodels|EnvironmentModels]] · `module_api` · call · `src/environment/physical_models.jl`
- `api` → [[envana.env_planets_topographyharmonicsworkspace_bang|TopographyHarmonicsWorkspace!]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.density_models__batch_elapsed_time|_batch_elapsed_time]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__density_scalar_for_batch|_density_scalar_for_batch]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__gram_core_density_state|_gram_core_density_state]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__gram_default_surrogate_file|_gram_default_surrogate_file]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__gram_lock_scope|_gram_lock_scope]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__gram_not_loaded_error|_gram_not_loaded_error]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__gram_point_density|_gram_point_density]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__nrlmsise_ap_bins|_nrlmsise_ap_bins]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__nrlmsise_ap_slot_index|_nrlmsise_ap_slot_index]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__nrlmsise_ap_value|_nrlmsise_ap_value]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__nrlmsise_density_state|_nrlmsise_density_state]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__nrlmsise_eval_datetime|_nrlmsise_eval_datetime]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__nrlmsise_flux_value|_nrlmsise_flux_value]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__nrlmsise_provider_indices|_nrlmsise_provider_indices]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__nrlmsise_resolved_indices|_nrlmsise_resolved_indices]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__nrlmsise_space_indices_ap_bin|_nrlmsise_space_indices_ap_bin]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__nrlmsise_space_indices_ap_vector|_nrlmsise_space_indices_ap_vector]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__nrlmsise_space_indices_f107|_nrlmsise_space_indices_f107]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__nrlmsise_space_indices_f107a|_nrlmsise_space_indices_f107a]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__nrlmsise_space_indices_indices|_nrlmsise_space_indices_indices]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__nrlmsise_space_indices_lookup|_nrlmsise_space_indices_lookup]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__piecewise_layer_index|_piecewise_layer_index]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__planet_polyfit_valid_max_altitude_m|_planet_polyfit_valid_max_altitude_m]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__planet_polyfit_valid_min_altitude_m|_planet_polyfit_valid_min_altitude_m]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__polyfit_density|_polyfit_density]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__polyfit_eval_altitude_km|_polyfit_eval_altitude_km]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__polyfit_log_density|_polyfit_log_density]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models__validate_density_batch_lengths|_validate_density_batch_lengths]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models_clear_gram_offline_surrogate_cache_bang|clear_gram_offline_surrogate_cache!]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models_clear_gram_static_grid_cache_bang|clear_gram_static_grid_cache!]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models_constantdensitymodel|ConstantDensityModel]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models_exponentialatmospheremodel|ExponentialAtmosphereModel]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models_gramatmospheremodel|GRAMAtmosphereModel]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models_gramatmospheremodelsurrogate|GRAMAtmosphereModelSurrogate]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models_init_nrlmsise_space_indices_bang|init_nrlmsise_space_indices!]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models_noatmospheremodel|NoAtmosphereModel]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models_nrlmsise00atmospheremodel|NRLMSISE00AtmosphereModel]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models_nrlmsise00spaceindicesprovider|NRLMSISE00SpaceIndicesProvider]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models_piecewiseexponentialatmospheremodel|PiecewiseExponentialAtmosphereModel]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models_polynomialfitatmospheremodel|PolynomialFitAtmosphereModel]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models_precompute_gram_static_grids_bang|precompute_gram_static_grids!]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.density_models_tabulatedflightatmospheremodel|TabulatedFlightAtmosphereModel]] · `module_api` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[environment.gravity_models__gravity_gradient_torque_body|_gravity_gradient_torque_body]] · `module_api` · call · `src/environment/gravity/gravity_models.jl`
- `api` → [[environment.gravity_models__gravity_runtime_field|_gravity_runtime_field]] · `module_api` · call · `src/environment/gravity/gravity_models.jl`
- `api` → [[environment.gravity_models__inverse_squared_gravity_accel|_inverse_squared_gravity_accel]] · `module_api` · call · `src/environment/gravity/gravity_models.jl`
- `api` → [[environment.gravity_models__inverse_squared_j2_gravity_accel|_inverse_squared_j2_gravity_accel]] · `module_api` · call · `src/environment/gravity/gravity_models.jl`
- `api` → [[environment.gravity_models_calcforcetorque|calcForceTorque]] · `module_api` · call · `src/environment/gravity/gravity_models.jl`
- `api` → [[environment.gravity_models_constantgravitymodel|ConstantGravityModel]] · `module_api` · call · `src/environment/gravity/gravity_models.jl`
- `api` → [[environment.gravity_models_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `module_api` · call · `src/environment/gravity/gravity_models.jl`
- `api` → [[environment.gravity_models_gravity_backbone_structure|gravity_backbone_structure]] · `module_api` · call · `src/environment/gravity/gravity_models.jl`
- `api` → [[environment.gravity_models_gravity_gradient|gravity_gradient]] · `module_api` · call · `src/environment/gravity/gravity_models.jl`
- `api` → [[environment.gravity_models_inversesquaredgravitymodel|InverseSquaredGravityModel]] · `module_api` · call · `src/environment/gravity/gravity_models.jl`
- `api` → [[environment.gravity_models_j2_secular_rates|j2_secular_rates]] · `module_api` · call · `src/environment/gravity/gravity_models.jl`
- `api` → [[environment.planet_shapes_earth_elevation_bang|Earth_elevation!]] · `module_api` · call · `src/environment/ephemerides/planet_shapes.jl`
- `api` → [[environment.planet_shapes_mars_elevation_bang|Mars_elevation!]] · `module_api` · call · `src/environment/ephemerides/planet_shapes.jl`
- `api` → [[environment.planet_shapes_venus_elevation_bang|Venus_elevation!]] · `module_api` · call · `src/environment/ephemerides/planet_shapes.jl`
- `api` → [[environment.planets__furnsh_first_existing|_furnsh_first_existing]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets__furnsh_first_existing_if_available|_furnsh_first_existing_if_available]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets__furnsh_mars_pck|_furnsh_mars_pck]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets__furnsh_mars_system_kernel|_furnsh_mars_system_kernel]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets__furnsh_once|_furnsh_once]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets__furnsh_planetary_kernel|_furnsh_planetary_kernel]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets__furnsh_required|_furnsh_required]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets__gravity_constants_kernel_if_available|_gravity_constants_kernel_if_available]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets__planetary_kernel_override_relpath|_planetary_kernel_override_relpath]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets__reset_furnished_kernels_bang|_reset_furnished_kernels!]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets__spice_backed_planet_kwargs|_spice_backed_planet_kwargs]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets__spice_body_fixed_frame|_spice_body_fixed_frame]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets__spice_body_gm_m3s2|_spice_body_gm_m3s2]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets__spice_body_pool_name|_spice_body_pool_name]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets__spice_body_radii_m|_spice_body_radii_m]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets_earth|Earth]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets_mars|Mars]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets_moon|Moon]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets_planets|Planets]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets_titan|Titan]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets_topographyharmonicsworkspace|TopographyHarmonicsWorkspace]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.planets_venus|Venus]] · `module_api` · call · `src/environment/ephemerides/planets.jl`
- `api` → [[environment.simple_ephemerides__earth_gmst_iau82_rad|_earth_gmst_iau82_rad]] · `module_api` · call · `src/environment/ephemerides/simple_ephemerides.jl`
- `api` → [[environment.simple_ephemerides__initial_time_datetime|_initial_time_datetime]] · `module_api` · call · `src/environment/ephemerides/simple_ephemerides.jl`
- `api` → [[environment.simple_ephemerides__rotation_about_spin_axis|_rotation_about_spin_axis]] · `module_api` · call · `src/environment/ephemerides/simple_ephemerides.jl`
- `api` → [[environment.simple_ephemerides__spice_body_fixed_frame|_spice_body_fixed_frame]] · `module_api` · call · `src/environment/ephemerides/simple_ephemerides.jl`
- `api` → [[environment.simple_ephemerides__spice_planet_frame_lpi|_spice_planet_frame_lpi]] · `module_api` · call · `src/environment/ephemerides/simple_ephemerides.jl`
- `api` → [[environment.simple_ephemerides__spice_position_j2000_m_unlocked|_spice_position_j2000_m_unlocked]] · `module_api` · call · `src/environment/ephemerides/simple_ephemerides.jl`
- `api` → [[environment.simple_ephemerides_ephemerides_cache_key|ephemerides_cache_key]] · `module_api` · call · `src/environment/ephemerides/simple_ephemerides.jl`
- `api` → [[environment.simple_ephemerides_planet_frame_lpi|planet_frame_lpi]] · `module_api` · call · `src/environment/ephemerides/simple_ephemerides.jl`
- `api` → [[environment.simple_ephemerides_spiceephemeridesmodel|SpiceEphemeridesModel]] · `module_api` · call · `src/environment/ephemerides/simple_ephemerides.jl`
- `api` → [[grp.src_environment_atmosphere|environment/atmosphere/]] · `members_in` · call · `src/environment/atmosphere/density_models.jl`
- `api` → [[grp.src_environment_ephemerides|environment/ephemerides/]] · `members_in` · call · `src/environment/ephemerides/ephemerides_models.jl`
- `api` → [[grp.src_environment_gravity|environment/gravity/]] · `members_in` · call · `src/environment/gravity/gravity_effectors.jl`
- `api` → [[grp.src_environment_small|environment/ (small files)]] · `members_in` · call · `src/environment/aerodynamics/aerodynamic_effectors.jl`
- `api` → [[module.core|SimulationModel]] · `environment` · call · `src/core/simulation_model.jl:17-50`
<!-- vulcan:connections:end -->

## Limitations
- `getDensity(::NRLMSISE00AtmosphereModel, h, lat, lon, el_time, wind)` without `p` throws:
  the six-argument method has no scenario epoch and directs the caller to the seven-argument
  form that maps `el_time` through `p.args.initial_time`.
- `_polyfit_log_density` clamps its exponent into `[_POLYFIT_LOG_DENSITY_MIN, _POLYFIT_LOG_DENSITY_MAX]`,
  so densities saturate rather than under/overflowing at extreme altitudes.
- `ExponentialAtmosphereModel` rejects `H <= 0.0` and `valid_min_altitude_m > valid_max_altitude_m`
  at construction, but performs no altitude-range check at evaluation time.
- `gravity_gradient` returns the zero vector for non-finite or non-positive radius, so a
  corrupted position silently produces no torque instead of raising.
- SPICE-backed frame and position queries serialise on the global `SPICE_LOCK`, capping
  multi-threaded throughput for `SpiceEphemeridesModel` scenarios.

## Provenance
Mapped from `src/environment/physical_models.jl`, with density physics from
`src/environment/atmosphere/density_models.jl`, gravity from
`src/environment/gravity/gravity_models.jl`, and epoch/frame handling from
`src/environment/ephemerides/simple_ephemerides.jl` and
`src/environment/ephemerides/ephemerides_models.jl`.
