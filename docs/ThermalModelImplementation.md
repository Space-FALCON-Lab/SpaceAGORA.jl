# Thermal Model Implementation

This document describes the aerodynamic heating implementation in SpaceAGORA.jl, covering the physical model, software architecture, simulation integration, and GNC guidance usage. It is intended as a reference for adding thermal modeling detail to the JAIS paper.

---

## Physical Model: Free-Molecular Maxwellian Heat Transfer

The primary thermal model is `MaxwellianHeat`, which computes surface heat flux in the free-molecular (FM) regime using the Maxwellian velocity distribution of incident gas molecules. This is the appropriate model for high-altitude atmospheric flight where the Knudsen number is large (mean free path much larger than the vehicle scale), a regime characteristic of aerobraking corridors at Mars and other bodies.

### Governing Equation

The heat rate per unit area (W/cm²) is:

$$
\dot{q} = \alpha_{TAF} \, \rho \, R \, T_p \sqrt{\frac{R T_p}{2\pi}} \left[ \left( S^2 + \frac{\gamma}{\gamma-1} - \frac{\gamma+1}{2(\gamma-1)} \frac{T_w}{T_p} \right) \mathcal{F}(S, \alpha) - \frac{1}{2} e^{-(S\sin\alpha)^2} \right] \times 10^{-4}
$$

where:

- $\alpha_{TAF}$: thermal accommodation factor (0–1); represents the fraction of kinetic energy transferred to the surface per molecular collision. A value of 1.0 indicates complete accommodation (full energy transfer).
- $\rho$: atmospheric density (kg/m³)
- $R$: specific gas constant of the atmosphere (J/kg·K)
- $T_p$: free-stream temperature (K)
- $T_w$: wall temperature (K); in the current implementation, $T_w = T_p$ (thermal equilibrium assumption)
- $S$: speed ratio (ratio of bulk flow speed to most-probable molecular thermal speed), computed from Mach number as $S = \sqrt{\gamma/2} \cdot M$
- $\alpha$: angle of attack of the surface panel (rad)
- $\gamma$: ratio of specific heats of the atmosphere

The flux factor $\mathcal{F}(S, \alpha)$ depends on the thermal contact boundary condition:

**No thermal contact** (open surface, `thermal_contact=false`):

$$
\mathcal{F}(S,\alpha) = e^{-(S\sin\alpha)^2} + \sqrt{\pi}\, S\sin\alpha \cdot \mathrm{erf}(S\sin\alpha) \cdot e^{(S\sin\alpha)^2}
$$

**Full thermal contact** (`thermal_contact=true`, default):

$$
\mathcal{F}(S,\alpha) = e^{-(S\sin\alpha)^2} + \sqrt{\pi}\, S\sin\alpha \cdot \left(1 + \mathrm{erf}(S\sin\alpha)\right)
$$

The thermal contact flag determines whether reflected molecules (those departing the surface) contribute to the total heat flux. With `thermal_contact=true` the re-emitted flux is included, corresponding to a wall that has reached thermal equilibrium with impinging molecules.

### Derived Intermediate Quantities

The model also computes recovery and stagnation temperatures internally, though they are not currently used to modify $T_w$:

- Stagnation temperature: $T_0 = T_p \left(1 + \frac{\gamma-1}{\gamma} S^2\right)$
- Recovery temperature: $T_r = T_p + \frac{\gamma}{\gamma+1} r' (T_0 - T_p)$

where $r'$ is an energy accommodation factor that differs between the thermal contact cases, encoding the rotational energy exchange fraction via an erf-based function of $S\sin\alpha$.

### Stagnation Heating Models (Legacy/Entry Corridor)

Two additional blunt-body models are present in `thermal_models.jl` but are not currently wired into the main simulation dispatch:

- `heatrate_convective(S, T, m, ρ, v, α)`: convective stagnation-point heating using a continuum-regime correlation: $\dot{q}_{conv} = k \sqrt{\rho / r_n} \, v^3 \times 10^{-4}$, where $k$ is a planetary constant and $r_n$ is the nose radius.
- `heatrate_radiative(S, T, m, ρ, v, α)`: radiative stagnation-point heating using the Tauber–Sutton formulation with a velocity-dependent $f$-factor interpolated from a tabulated lookup (velocities 0–16,000 m/s), exponent $a = 1.072 \times 10^6 v^{-1.88} \rho^{-0.325}$, and constant $C = 4.736 \times 10^4$.
- `heatrate_convective_radiative`: sum of both; currently returns only the convective term (radiative term commented out).

These are retained for entry/descent scenarios but are inactive in the aerobraking simulation path.

---

## Software Architecture

### Type Hierarchy

```
AbstractThermalModel          (abstract_types.jl)
└── MaxwellianHeat{P}         (thermal_models.jl)
        thermal_accomodation_factor :: Float64
        planet                      :: P <: AbstractPlanet
        thermal_contact             :: Bool = false
```

`AbstractThermalModel` is a stable extension interface exported from `SpaceAGORA.jl`. Users can implement custom thermal models by subtyping it and defining `getHeatRate(model, S, T, ρ, v, α)::Float64`. Startup validation (`_validate_thermal_model_support!` in `setup.jl`) checks for this method dispatch at simulation initialization and throws a descriptive error if it is absent.

### Module Layout

```
src/vehicle/thermal/
    thermal_models_module.jl   — VehicleThermalModels module; exports MaxwellianHeat, getHeatRate
    thermal_models.jl          — struct definition and method implementations
```

`VehicleThermalModels` is re-exported from `SimulationModel` and exposed at the top-level `SpaceAGORA` namespace.

### Configuration

The thermal model is a typed field of `EnvironmentModel{P, D, E, T}`:

```julia
@kwdef struct EnvironmentModel{..., T <: AbstractThermalModel}
    thermal_model::T
    ...
end
```

The type parameter `T` is propagated through `SimulationConfiguration`, enabling full compile-time specialization of the thermal dispatch path — no runtime type checks or branching during integration.

Default construction (used in the no-GRAM preset and verification scenarios):

```julia
MaxwellianHeat(thermal_accomodation_factor = 1.0, planet = planet_model)
```

---

## Simulation Integration

### Thermal Callback

Heat rates are computed outside the ODE right-hand side (RHS) in a dedicated `DiscreteCallback` that fires at every integration step. This decouples the potentially expensive thermal evaluation from the ODE stiffness structure and allows independent threading.

The callback (`get_thermal_callback` in `thermal_callbacks.jl`) executes the following per-satellite:

1. Read atmospheric state (`ρ`, `T`, wind vector) from `shared_buffers` populated by the density callback.
2. Transform the inertial velocity to the planet-fixed frame and add the wind vector to form the relative velocity $v_{rel}$.
3. Compute the local speed of sound: $a = \sqrt{\gamma R T}$, Mach number $M = |v_{rel}|/a$, and speed ratio $S = \sqrt{\gamma/2} \cdot M$.
4. Loop over all aerodynamic links (surface panels) of the spacecraft. Each link carries a precomputed angle of attack $\alpha_j$.
5. Call `getHeatRate(thermal_model, S, T, ρ, v, α_j)` and store the result in `shared_buffers.heat_rates[i][j]`.
6. Non-finite or non-positive results are clamped to zero.

The callback is initialized to fire immediately at $t=0$ and runs as a `DiscreteCallback(condition=always_true, affect!)`.

### Heat Load State Integration

Heat flux values from `shared_buffers.heat_rates` are promoted into the ODE state vector via `_assign_heat_rate_derivative!` (in `dynamics_rhs.jl`). This means the ODE integrates cumulative heat load directly as an additional state:

$$
\frac{d}{dt}Q_j(t) = \dot{q}_j(t)
$$

so $Q_j(t)$ (J/cm²) at each panel is a native ODE state subject to the integration tolerances. The tolerances applied to heat load states are:

```
reltol_heat_load = 1e-7
abstol_heat_load = 1e-9
```

### Saved Fields

Two quantities are saved to the output feather/CSV files at each data step:

- **`heat_rate`**: instantaneous total heat rate per satellite (W/cm²), summed across panels via `_save_heat_rate` (sum of `shared_buffers.heat_rates[i]`).
- **`heat_load`**: cumulative heat load per satellite (J/cm²), taken directly from the ODE state `u.sc[i].heat_loads`.

### Multi-Satellite Parallelism

When simulating constellation scenarios the thermal callback supports adaptive threading via `ParallelPolicy.threaded_foreach_persistent(:thermal_callback, num_sats, decision.allotment)`. The threading decision is governed by `ParallelConfig.thermal_callback_parallel_mode` (default `"auto"`), which applies the same adaptive policy framework used by the density and control callbacks. Timing observations are recorded to feed back into the adaptive policy.

---

## GNC Guidance Usage

The aerobraking guidance system uses a vectorized version of the same Maxwellian heat rate formula for constraint-based trajectory optimization. The function `heat_rate_calc` (in `tracking_executor.jl`) evaluates the heat rate over arrays of trajectory states:

```julia
function heat_rate_calc(taf, ρ, T_w, T_p, R, γ, S, angle)
    first_term = ρ .* (1e-4 * taf * R * T_p * sqrt.(R * T_p / (2π)))
    term_a = exp.(-(S .* sin.(angle)).^2) .+ sqrt(π) * S .* sin.(angle) .* (1 .+ erf.(S .* sin.(angle)))
    term_b = (S.^2 .+ γ/(γ-1) .- (γ+1)/(2(γ-1)) * (T_w ./ T_p)) .* term_a
    return (term_b .- 0.5 * exp.(-(S .* sin.(angle)).^2)) .* first_term
end
```

This is algebraically identical to `getHeatRate` in the Maxwellian model (with `thermal_contact=true`), operating on trajectory state vectors rather than scalar values.

### Guidance Integration Points

The heat rate function is called in the following guidance contexts:

- **T-EDG (time-based targeting solver)** (`targeting_solver.jl`): heat rate evaluated over the predicted trajectory; angle of attack $\alpha$ is adjusted via co-state (Lagrangian multiplier) integration to minimize heat load accumulation while meeting orbit-lowering constraints.
- **E-EDG (energy-based solvers)** (`energy_profile_solver.jl`, `second_switch_solver.jl`): heat rate limits inform the switch-time computation between high-drag (banked) and low-drag (ballistic) flight segments.
- **Heat load convergence loop** (`heat_rate_models.jl:func`): iterates angle-of-attack profiles to convergence on the total integrated heat load $Q = \int \dot{q}\,dt$, stopping when $|Q_k - Q_{k-1}| < 10^{-3}$.
- **Heat load constraint** (`heat_rate_models.jl`): $\delta Q = Q - Q_{limit}$ drives a root-finding step (via `heat_load_sol` parameter) that determines the maximum-depth drag passage.

The guidance path uses a configurable multiplicative factor `args[:multiplicative_factor_heatload]` applied to the thermal accommodation factor, providing a safety margin between the predicted and maximum allowable heat load.

---

## Key Parameters Summary

| Parameter | Symbol | Location | Notes |
|-----------|--------|----------|-------|
| Thermal accommodation factor | $\alpha_{TAF}$ | `MaxwellianHeat.thermal_accomodation_factor` | 0–1; default 1.0 |
| Thermal contact mode | — | `MaxwellianHeat.thermal_contact` | `false` = open surface, `true` = wall equilibrium |
| Speed ratio | $S$ | computed in thermal callback | $\sqrt{\gamma/2} \cdot M$ |
| Planet gas constant | $R$ | `planet.R` | Specific gas constant (J/kg·K) |
| Ratio of specific heats | $\gamma$ | `planet.γ` | Used in $S$, stagnation temperature, heat rate |
| Heat load tolerances | — | `IntegrationTolerances` | reltol=1e-7, abstol=1e-9 |
| Guidance heat load margin | — | `args[:multiplicative_factor_heatload]` | Scales TAF in predictors |

---

## Extension Interface

Custom thermal models can be registered by:

1. Subtyping `AbstractThermalModel`.
2. Defining `getHeatRate(model::MyModel, S::Float64, T::Float64, ρ::Float64, v::Float64, α::Float64)::Float64`.
3. Passing an instance as `thermal_model` in `EnvironmentModel`.

The simulation engine validates this interface at startup and will error with a descriptive message if the method is missing.
