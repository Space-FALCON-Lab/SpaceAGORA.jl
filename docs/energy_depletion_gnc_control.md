# Energy-Depletion Aerobraking Guidance and Control

## Purpose and Scope

This document describes the energy-depletion guidance and control implementation exactly as
it currently exists. It covers every live mode decision, constraint decision, switch-planning
branch, numerical fallback, state reset, and solar-panel angle-of-attack (AoA) command.

The implementation ports the Odyssey control laws from `deprecated_control_code/src/control`
into the current callback architecture without changing the simulation right-hand side (RHS).
Guidance chooses between maximum energy depletion, targeting, and safe low drag. Control
computes switch times, applies instantaneous heat-rate and structural constraints, and rotates
the configured panel links. The control effector itself returns zero force and torque; the
existing aerodynamic dynamics observe the changed link orientations.

The authoritative implementation files are:

- `src/gnc/guidance/target_energy_bracketing.jl`
- `src/gnc/control/heat_rate_control.jl`
- `src/gnc/control/struct_load_control.jl`
- `src/gnc/control/heat_load_control.jl`
- `src/gnc/control/targeting_control.jl`
- `src/gnc/control/control_hooks.jl`

The Odyssey executable example is
`examples/AGORA_Odyssey_Energy_Depletion_Control_Test.jl`.

## Units and Symbols

| Symbol | Meaning | Units |
| --- | --- | --- |
| `alpha` or `α` | Controlled panel AoA | rad |
| `rho` or `ρ` | Atmospheric density | kg/m3 |
| `V` | Wind-relative speed | m/s |
| `S` | Molecular speed ratio | dimensionless |
| `q` | Dynamic pressure, `rho*V^2/2` | Pa |
| `qdot` | Surface heat rate | W/cm2 |
| `Q` | Integrated surface heat load | J/cm2 |
| `CL`, `CD` | Lift and drag coefficients | dimensionless |
| `A` | Reference area | m2 |
| `epsilon` | Specific orbital energy | J/kg |
| `mu` | Planet gravitational parameter | m3/s2 |
| `Rp` | Planet equatorial radius | m |
| `gamma` | Flight-path angle, except where used for gas specific-heat ratio | rad |

`alpha_max` normally means high drag and `alpha_min` means low drag. The Odyssey run uses
`alpha_max = pi/2` and `alpha_min = 0`.

## Callback Rates and Execution Ownership

The Odyssey example configures:

| Operation | Period | Rate |
| --- | ---: | ---: |
| Guidance callback | `1/3 s` | 3 Hz |
| Control callback | `0.1 s` | 10 Hz |
| CSV save callback | `1.0 s` | 1 Hz |

Guidance owns only mode selection, energy bracketing, and target-energy construction. Heat-load
and targeting switch times are solved in control. Heat-rate and structural AoA constraints are
recomputed by the live control callback at 10 Hz. Targeting is planned once. Initial heat-load
switch planning is performed once at drag-passage entry, followed only by the explicitly
described second-switch correction and security logic.

## Configuration Validation

`AerobrakingEnergyDepletionConfig` accepts the following mode symbols:

```julia
guidance_modes = (:max_energy_depletion,)
guidance_modes = (:targeting,)
guidance_modes = (:targeting, :max_energy_depletion)

max_energy_submodes = (:heat_rate, :structural_load, :heat_load)
heat_load_switch_solver = :closed_form       # or :tpbvp_integration
```

Construction applies these decisions:

1. A single symbol is converted to a one-element tuple. Other iterables are converted element
   by element with `Symbol(value)`.
2. `guidance_modes` must be nonempty and contain only `:max_energy_depletion` or `:targeting`.
3. `max_energy_submodes` must be nonempty and contain only `:heat_rate`, `:structural_load`, or
   `:heat_load`.
4. `heat_load_switch_solver` must be `:closed_form` or `:tpbvp_integration`.
5. `controlled_panel_links` must be nonempty and every index must be positive. Bounds and root
   link checks occur later during actuation.
6. AoA bounds must be finite and satisfy `0 <= min_alpha_rad <= max_alpha_rad`.
7. `planning_horizon_s` and `switch_recompute_interval_s` must be finite and positive.

The current defaults are a 5000 s planning horizon, 20 s switch-recompute interval, enabled
second-switch re-evaluation, and enabled heat-load security mode. The
`switch_recompute_interval_s` field is validated and stored but is not read by the current
switch scheduler; the scheduler uses the explicit 10 s, 3 s, and 0.8 s thresholds documented
below.

## Shared State

`AerobrakingEnergyDepletionState` stores one element per spacecraft. New state is initialized
as follows:

| Field group | Initial value |
| --- | --- |
| Selected mode | `:inactive` |
| Targeting, safe-low-drag, bracket flags | `false` |
| Bracket evaluation count | `0` |
| Target and bracket energies | `NaN` |
| Heat-load switches | `(Inf, Inf)` |
| Heat-load solved/passage/security flags | `false` |
| Heat-load entry time | `NaN` |
| Last heat-load re-evaluation | `-Inf` |
| Previous heat load | `0` |
| Targeting switch | `Inf` |
| Last switch solve time | `-Inf` |
| Last AoAs and telemetry | `NaN` |

An invalid spacecraft index causes either callback to return immediately without changing
state.

## Environment Calculation

The live control and numerical predictors transform inertial position and velocity into the
planet-fixed frame. Density, temperature, and wind are requested at the resulting altitude,
latitude, longitude, and simulation time. Wind in North-East-Up order is transformed to a
planet-fixed vector. Relative velocity is

```math
\mathbf V_{rel}=\mathbf V_{planet-fixed}-\mathbf V_{wind}.
```

The implementation then computes

```math
V=\lVert\mathbf V_{rel}\rVert,
\qquad a_s=\sqrt{\gamma_g R T},
\qquad S=\sqrt{\frac{\gamma_g}{2}}\frac{V}{a_s},
\qquad q=\frac{1}{2}\rho V^2.
```

Density is clamped to at least zero and temperature to at least machine epsilon. The vehicle
is considered inside the drag passage when altitude is finite and

```math
h \le 1000\,EI_{km}.
```

## End-to-End Guidance Flow

The guidance callback evaluates top-level modes in this strict priority order:

1. If `:targeting` is present anywhere in `guidance_modes`, run targeting energy bracketing.
   The presence of `:max_energy_depletion` does not bypass this branch.
2. Otherwise, if `:max_energy_depletion` is present, set mode directly to
   `:max_energy_depletion`, clear targeting and safe-low-drag flags, and leave
   `energy_bracketing_evaluated = false`.
3. Otherwise, select `:safe_low_drag`, disable targeting, and set the safe-low-drag flag.

Therefore a maximum-energy-only run never performs an energy bracket.

### Targeting bracket scheduling

When targeting is configured, bracketing follows these decisions:

1. If the bracket was already evaluated, return immediately. It is never recomputed.
2. Build the current environment. If outside the drag passage, return without setting any
   mode or bracket state. The initial `:inactive` mode remains until entry.
3. Read the maximum nonnegative finite heat load among the configured controlled links. A
   missing heat-load state is treated as zero.
4. Predict a complete low-drag passage and a maximum-energy-depletion passage.
5. Set the reachable interval to the extrema of those two predicted exit energies.
6. Construct the desired exit energy from target apoapsis and the vacuum J2/N-body
   correction.
7. Mark the bracket evaluated and increment its count exactly once.
8. If the target lies in the interval within tolerance, select `:targeting`.
9. If unreachable and `:max_energy_depletion` is available, select it as fallback.
10. If unreachable and maximum energy depletion is unavailable, select `:safe_low_drag`.

The reachability tolerance is

```math
epsilon_{tol}=10^{-6}\max(|epsilon_{min}|,|epsilon_{max}|,1).
```

The target must be finite and satisfy

```math
epsilon_{min}-epsilon_{tol}\le epsilon_{target}
\le epsilon_{max}+epsilon_{tol}.
```

`_edg_set_targeting_fallback!` implements the same fallback state assignment, but the current
live bracketing function performs the assignment inline and does not call that helper.

## End-to-End Control Callback

Every 10 Hz control callback performs this sequence:

1. Return immediately for an invalid spacecraft index.
2. If mode is `:inactive` and targeting is not configured, select maximum energy depletion if
   available, otherwise safe low drag. If targeting is configured, leave the mode inactive so
   guidance can perform entry bracketing.
3. Read spacecraft state and calculate the current environment.
4. Read the maximum controlled-panel heat load.
5. Run targeting or heat-load switch scheduling as applicable.
6. Determine whether the heat-load low-drag interval is active.
7. Select a base AoA from mode and switch state.
8. Apply enabled instantaneous constraints unless heat-load low drag is active.
9. Clamp and store the final command and telemetry.
10. Rotate only the configured panel links.

While a targeting-configured vehicle remains `:inactive` outside EI, the base-AoA fallthrough
is `alpha_max`. Thus the panels are held high until entry guidance selects targeting or a
fallback. Atmospheric forces are negligible outside the drag passage, but the link command is
still applied.

### Base-AoA decision

The base command is selected in this order:

1. If mode is `:safe_low_drag` or the safe-low-drag flag is true, use `alpha_min`.
2. If mode is active targeting, use `alpha_max` before the targeting switch and `alpha_min`
   at and after the switch. With an infinite switch, targeting stays at `alpha_max`.
3. If mode is maximum energy depletion and the current time is inclusively inside both finite
   heat-load switches, use `alpha_min`.
4. Otherwise use `alpha_max`.

Heat-load low drag is active only when mode is `:max_energy_depletion`, `:heat_load` is an
enabled submode, both switches are finite, and

```math
t_{s1}\le t\le t_{s2}.
```

## Heat-Rate Control

### Heat-rate model

The thermal accommodation factor is read from the configured thermal model and defaults to
one if the field is absent. Wall temperature is set equal to atmospheric temperature,
`T_w = T_p`. Define

```math
s_n=S\sin\alpha,
\qquad E=e^{-s_n^2},
\qquad A=E+\sqrt{\pi}s_n(1+\operatorname{erf}(s_n)),
```

```math
B=S^2+\frac{\gamma_g}{\gamma_g-1}
-\frac{\gamma_g+1}{2(\gamma_g-1)}\frac{T_w}{T_p},
```

```math
L=10^{-4}\,TAF\,\rho R T_p\sqrt{\frac{R T_p}{2\pi}}.
```

The controlled heat rate is

```math
\dot Q(\alpha)=L\left(BA-\frac{1}{2}E\right).
```

The `1e-4` factor converts the legacy area units to W/cm2. A nonfinite result is treated as
infinite by the root evaluator. The telemetry evaluator instead returns zero for invalid
environment inputs and clamps a nonfinite or negative final rate to zero.

### Heat-rate AoA decision

The constraint solver uses a minimum constraint angle

```math
\alpha_{constraint,min}=\max(\alpha_{configured,min},10^{-4}\text{ rad}).
```

This is intentionally distinct from a bang-bang `alpha_min = 0` command.

Given a base upper bound `alpha_base`, the solver follows these decisions:

1. If the limit is nonfinite or nonpositive, return `alpha_base` unchanged.
2. If density, temperature, or speed ratio is invalid or nonpositive, return `alpha_base`.
3. Set the root target to `qdot_limit - 1e-5 W/cm2`. If this is nonpositive, return the
   constraint minimum.
4. Evaluate heat rate at the base maximum and constraint minimum.
5. If heat rate at the maximum is strictly below the target, return the maximum.
6. If heat rate at the minimum is strictly above the target, no feasible angle exists; return
   the minimum protective command.
7. Otherwise solve `qdot(alpha) - (qdot_limit - 1e-5) = 0` with Brent's method on
   `[alpha_min, alpha_base]`.
8. If root finding throws, return the minimum.
9. If the root is finite and lies in `[0, alpha_base]`, return it; otherwise return zero.

The `alpha_past` argument is passed through live and predicted paths for legacy interface
parity but is not used by the current root equation.

The live wrapper also returns the base command without solving when the base is below the
`1e-4 rad` constraint floor. This preserves a commanded zero-AoA heat-load or targeting phase.

## Structural-Load Control

### Legacy flat-plate coefficients

For speed ratio `s`, reflection coefficient `sigma`, and AoA `alpha`, define

```math
s_n=s\sin\alpha,
\qquad E=e^{-s_n^2},
\qquad F=1+\operatorname{erf}(s_n),
```

```math
C_N=\frac{
\left(\frac{2-\sigma}{\sqrt\pi}s_n+\frac\sigma2\right)E+
\left((2-\sigma)(s_n^2+\tfrac12)+\frac{\sigma\sqrt\pi s_n}{2}\right)F
}{s^2},
```

```math
C_A=\frac{\sigma\cos\alpha}{\sqrt\pi s}
\left(E+\sqrt\pi s_nF\right),
```

```math
C_L=C_N\cos\alpha-C_A\sin\alpha,
\qquad
C_D=C_A\cos\alpha+C_N\sin\alpha.
```

The controlled links are evaluated at panel AoA. The remaining bus area is evaluated at
`pi/2`. Coefficients are area weighted:

```math
C_{D,sc}=\frac{C_{D,panel}A_{panel}+C_{D,bus}A_{bus}}{A_{total}},
```

with the same expression for lift. The root link reflection coefficient is used for both.

### Structural AoA decision

The configured structural limit has units of equivalent maximum-AoA dynamic pressure. The
reference drag area is

```math
(C_DA)_{ref}=C_D(\alpha_{configured,max})A_{total},
```

and the force limit is

```math
D_{limit}=q_{limit}(C_DA)_{ref}.
```

At the current dynamic pressure, the controller follows these decisions:

1. Outside the drag passage, return the base AoA.
2. If the base is below the `1e-4 rad` constraint floor, return the base AoA.
3. If limit or dynamic pressure is nonfinite or nonpositive, return the base AoA.
4. Clamp the root maximum to `[alpha_constraint_min, alpha_configured_max]`.
5. If drag at that maximum is at or below the force limit, return the maximum.
6. If drag at the constraint minimum is above the force limit, no feasible angle exists;
   return the minimum protective command.
7. Otherwise solve
   `q*CD(alpha)*A_total - D_limit = 0` with bisection on `[0, pi/2]`.
8. If the solve throws, return the constraint minimum.
9. If the result is finite and in `[0, alpha_base]`, return it; otherwise return zero.

The stored structural telemetry is the equivalent normalized pressure

```math
q_{struct}=q\frac{(C_DA)_{command}}{(C_DA)_{ref}}.
```

The numerical trajectory predictor separately uses the current per-link
`aerodynamic_coefficient_fM` model, area weighted over all links. Controlled links receive the
candidate AoA; every other link receives `pi/2`.

## Constraint Composition

The live command starts from the mode/switch base AoA. Heat-rate control is active when its
submode is selected and heat-load low drag is not active. Structural control uses the same
gate. The composition is exactly:

```math
\alpha_c=
\begin{cases}
\min(\alpha_{heat},\alpha_{struct}), & \text{both active},\\
\alpha_{heat}, & \text{heat rate only},\\
\alpha_{struct}, & \text{structural only},\\
\alpha_{base}, & \text{neither active}.
\end{cases}
```

The result is clamped to configured AoA bounds. During the heat-load minimum-AoA interval,
both instantaneous controls are disabled, so neither may raise the command. During a high-AoA
heat-load phase or the pre-switch targeting phase, both remain active and the smaller AoA wins.

## Heat-Load Control

### Requested command shape

Heat-load control is maximum-energy-depletion-only and enforces one high-low-high interval:

```math
\alpha(t)=
\begin{cases}
\alpha_{high}, & t<t_{s1},\\
\alpha_{min}, & t_{s1}\le t\le t_{s2},\\
\alpha_{high}, & t>t_{s2}.
\end{cases}
```

Any raw costate profile with more than one low interval is reduced to its first contiguous low
interval. Therefore only the first two transitions are applied.

### Planning inputs and fallbacks

Total area sums finite positive link areas and never returns less than machine epsilon.
Controlled area sums valid configured links. The current heat-load state is the maximum
nonnegative finite value on the controlled links. If the state has no `heat_loads` property,
it is zero. A valid positive state mass is used; otherwise link masses plus propellant are
summed and clamped to at least 1 kg.

The closed-form aerodynamic slope uses legacy aggregate coefficients at zero and maximum AoA:

```math
C_{D,slope}=\frac{C_D(\alpha_{max})-C_D(0)}{\pi/2}.
```

If the slope is nonfinite or effectively zero, the fallback is

```math
C_{D,slope}=\frac{2.2-0.8}{\pi/2},\qquad C_{D,0}=0.8.
```

The closed-form aggregate treats every non-root link as a panel for coefficient weighting,
whereas controlled reference area and live actuation use only `controlled_panel_links`.

### Drag-passage duration

Orbital elements are computed from current position and velocity. Invalid, circular,
unbound, or nonpositive-semimajor-axis states use the full planning horizon.

For an inbound state, outbound true anomaly is the mirror of the current true anomaly. This
preserves the legacy same-radius passage duration. For an outbound state, the code solves true
anomaly at the configured entry-interface radius:

```math
\cos\nu_{exit}=\frac{a(1-e^2)/r_{EI}-1}{e}.
```

Eccentric and mean anomalies are

```math
E=\operatorname{atan2}(\sqrt{1-e^2}\sin\nu,e+\cos\nu),
\qquad M=E-e\sin E,
```

wrapped to `[0,2pi)`. Duration is

```math
Delta t=\frac{(M_{exit}-M_0)\bmod 2\pi}{\sqrt{\mu/a^3}}.
```

A nonfinite duration or one not greater than 1 s becomes `min(planning_horizon,1000 s)`.
Otherwise duration is capped by the planning horizon. Time grids contain at least 64 and at
most 20,000 points. Closed form normally requests 0.1 s spacing; numerical TPBVP requests 1 s.

### Active Chebyshev corrections

For Chebyshev polynomials `T0` through `T4`, the implementation evaluates

```math
\hat V=\frac{V_0-V_{center}}{V_{half}},
\qquad
\hat\gamma=\frac{\gamma_{0,deg}-\gamma_{center}}{\gamma_{half}},
```

```math
f_1=\sum_{i=0}^{4}c_iT_i(\hat V),
```

and a total-degree-four bivariate polynomial

```math
\log f_{2,base}=\sum_{p+q\le4}d_{pq}T_p(\hat\gamma)T_q(\hat V),
\qquad
f_2(t)=e^{\log f_{2,base}}\frac{t}{2t_p}.
```

The coefficient order is `(0,0),(1,0),(0,1),(2,0),(1,1),(0,2),(3,0),(2,1),`
`(1,2),(0,3),(4,0),(3,1),(2,2),(1,3),(0,4)`.

| Planet | `V_center`, `V_half` | `gamma_center`, `gamma_half` (deg) | `f1` coefficients |
| --- | --- | --- | --- |
| Mars | 4350, 450 | -5.25, 2.25 | 4.560850, -1.715925, -0.375150, 0.376025, -0.254500 |
| Venus | 9000, 600 | -4.375, 1.375 | 16.3939125, -0.6316500, -0.3870500, 0.3974500, -0.2209625 |
| Earth | 9075, 725 | -5.125, 2.125 | 19.807700, -3.164075, 0.592100, -0.113625, 0.043900 |
| Titan | 2050, 150 | -12.5, 6.5 | 148.480400, 10.480925, 0.030500, -0.003725, 0.000600 |

The complete `f2` vectors in the bivariate order above are:

```text
Mars  = [-4.74628450, -3.69824500, -6.06348525,  1.07688200,  3.54740275,
          1.27053200, -0.26357900, -1.12508600, -0.91570600, -0.17288575,
          0.045299125, 0.280380250, 0.224278500, -0.004521000, -0.025270125]
Venus = [-11.578333625, -6.131747500, -3.504510250, 1.440936250, 2.801654000,
          0.536784250, 0.049912500, -0.281389000, -0.520214000, -0.099717750,
         -0.094991375, -0.214584750, -0.037281250, 0.078927750, 0.016912750]
Earth = [-1.989471500, 0.523520500, -1.992902250, 0.462889750, 1.965380250,
          0.591760750, -0.182282500, -0.886610000, -0.871433000, -0.200818750,
          0.024188625, 0.286623750, 0.468237250, 0.332824000, 0.043238125]
Titan = [-8.397098625, -4.739801250, -1.473748250, 0.671454000, 1.024378750,
          0.184619500, -0.031104750, -0.118313000, -0.145654000, -0.027513750,
         -0.034697125, -0.048001000, 0.005169000, 0.019228250, 0.003789250]
```

An unsupported planet throws an `ArgumentError`; there is no generic correction fallback.

Density is evaluated as

```math
\rho(h)=\exp\left(\sum_{i=0}^{N}b_iT_i(\hat h)\right),
\qquad
\hat h=\frac{h_{km}-h_c}{h_w}.
```

The active normalized density domains are Mars `(175,125)`, Venus `(150,100)`, Earth
`(275.5,225)`, and Titan `(1025,975)` km for center and half-width. The complete coefficient
vectors are:

```text
Mars = [-22.904887634115337, -11.841971600850178, 2.449273531108809,
         0.07122528331622979, -0.3143478095124491, 0.1734745257301972,
        -0.14396272246410283, 0.06104818577407287, 0.0738268254164405,
        -0.11983534975043726, 0.052506209635111954, 0.0022260853200571835,
        -0.0245022530854842]
Venus = [-17.17603495187312, -15.971487996237625, 3.1916728583530163,
          0.8296122739608919, -0.891116626156514, -0.018054118127076735,
          0.2016004439406335, 0.0537350005662433, -0.06184791826147179,
         -0.002221184567720138, -0.03852447322489234, 0.11389766896205009,
         -0.023720859333601994]
Earth = [-20.5897273251208688, -8.56737992269037107, 3.45019353462444567,
         -1.76485096058359936, 0.601801952591328626, 0.0724429484383713046,
         -0.330361927785729759, 0.307509524503914389, -0.158408358148278777,
          0.00777052279665210299, 0.0781913274877799047, -0.0903365513831131534,
          0.0560317850649417070, -0.0122200602431319916, -0.0147378888840894624,
          0.0175437183414442650, -0.00428819900201210177, -0.0110920129900025574,
          0.0182659646193858177, -0.0151388155198109841, 0.00619700872461139817]
Titan = [-19.4278696179288133, -14.4426058071846857, 2.69604850426823095,
         -0.336162464509102088, 0.0900778126111024535, -0.136262683746282226,
          0.232875124056659416, -0.226385345471503546, 0.118549197857537245,
         -0.0289659248379974324, 0.0165054791585107323, -0.0478606382329083910,
          0.0650600435627339962, -0.0503265704035375031, 0.0226198266668309125,
         -0.00645321729948887882, 0.00509967575456080676, -0.00943672721749253159,
          0.00948561230043162319, -0.00620859370479269826, 0.00120768474131884517]
```

Unsupported planets throw an `ArgumentError`.

### Closed-form trajectory equations

Let predicted passage duration be `Delta t`, `t_p = Delta t/2`, and
`c_3 = V_0*gamma_0`. Then

```math
h(t)=h_0+c_3\left(t-\frac{t^2}{2t_p}\right),
\qquad
\gamma(t)=\frac{c_3(1-t/t_p)}{V(t)}.
```

AoA profiles are normalized to the grid: an empty profile becomes all maximum AoA, a long
profile is truncated, and a short profile is extended with its final value.

The code linearly interpolates aggregate coefficients with AoA:

```math
C_D(t)=C_{D,0}+\alpha(t)\frac{C_{D,max}-C_{D,0}}{\pi/2},
\qquad
C_L(t)=C_{L,0}+\alpha(t)\frac{C_{L,max}-C_{L,0}}{\pi/2}.
```

It then forms

```math
c_1=\frac{\rho C_DA}{2m}\alpha,
\qquad
c_2=\frac{\rho C_LA}{2m},
```

```math
epsilon=f_1+f_2\alpha_{max}\frac{A_{panels}}{A}
+f_2\frac{\pi}{2}\frac{A_{bus}}{A},
```

```math
k_1=c_2+\frac{1}{R_p+h},
\qquad
k_2=c_1c_3(1-t/t_p),
\qquad
k_3=-g_{ref}-epsilon.
```

Velocity uses the negative quadratic root plus an offset that forces the first sample to equal
the actual entry speed:

```math
V_{root}=\frac{k_2/k_1-\sqrt{(k_2/k_1)^2-4k_3/k_1}}{2},
\qquad
V(t)=\max(V_{root}(t)+V_0-V_{root}(0),\epsilon_{machine}).
```

### Constrained high-AoA branch

Heat-load control has a high branch and a low branch. The low branch is always
`alpha_min`; neither instantaneous controller may raise it. At every high-branch sample, the
same heat-rate and structural roots used by the live controller are evaluated. Disabled
constraints contribute `alpha_max`, so the applied high command is exactly

```math
\alpha_{high}(t)=\operatorname{clamp}\left(
\min\left(\alpha_{heat\ rate}(t),\alpha_{structural}(t)\right),
\alpha_{min},\alpha_{max}\right).
```

Consequently, the phase-2 command is

```math
\alpha(t)=\begin{cases}
\alpha_{high}(t),&t<t_{s1},\\
\alpha_{min},&t_{s1}\le t\le t_{s2},\\
\alpha_{high}(t),&t>t_{s2}.
\end{cases}
```

This composition is used in the closed-form predictor, numerical shooting integration,
full-trajectory replay, and live control. A unit test independently computes both constraint
profiles and checks that the combined profile equals their pointwise minimum.

### Numerical TPBVP shooting method

The numerical option is a forward shooting method. It does not infer the initial costates by
repeated backward profile updates. The three unknown shooting variables are
`lambda_V(0)`, `lambda_gamma(0)`, and `lambda_h(0)`. For each guess, fixed-grid RK4 propagates
the nine-state system

```math
y=(\mathbf r,\mathbf v,\lambda_V,\lambda_\gamma,\lambda_h)
```

forward over the predicted drag passage. Its trajectory equations are

```math
\dot{\mathbf r}=\mathbf v,
\qquad
\dot{\mathbf v}=\mathbf a_g+\mathbf a_{aero}.
```

The configured inverse-square or J2 gravity predictor supplies `a_g`. The current per-link
free-molecular aerodynamic predictor supplies `a_aero`. Define

```math
r=\lVert\mathbf r\rVert,
\qquad
V=\lVert\mathbf v\rVert,
\qquad
\gamma=\sin^{-1}\left(\frac{\mathbf r\cdot\mathbf v}{rV}\right),
\qquad
g=\lVert\mathbf a_g\rVert.
```

With `CD=CD0+alpha*CDslope`, the costate equations are

```math
\dot\lambda_V=-\frac{3k\rho V^2\alpha}{\pi}
+\lambda_V\frac{\rho A C_DV}{m}
-\lambda_\gamma\left(\frac{\rho A C_{L,0}}{2m}+\frac{g}{V^2}+\frac{1}{r}\right)
-\lambda_h\gamma,
```

```math
\dot\lambda_\gamma=\lambda_Vg-\lambda_hV,
```

```math
\dot\lambda_h=\frac{k\rho V^3\alpha}{\pi H}
-\lambda_V\left(\frac{\rho AC_DV^2}{2mH}+\frac{2g\gamma}{r}\right)
+\lambda_\gamma\left(\frac{\rho AC_{L,0}V}{2mH}
-\frac{2g}{rV}+\frac{V}{r^2}\right).
```

Scale height comes from density-model field `H`, otherwise planet field `H`, otherwise 7000 m,
and is clamped to at least 1 m.

Every RK4 stage converts its currently propagated inertial `r` and `v` into the planet-fixed
frame, computes altitude, latitude, and longitude, and calls the configured `getDensity` method
at absolute time `t_entry + tau`. Therefore, when the Odyssey run configures
`GRAMAtmosphereModel`, every shooting residual and every finite-difference Jacobian
perturbation queries GRAM directly. Density and temperature are not frozen or interpolated from
the seed trajectory. GRAM wind is also included in the relative velocity used by the
instantaneous heat-rate and structural constraints and by the aerodynamic acceleration.

The terminal residual solved by `NLsolve.nlsolve` is

```math
R=\begin{bmatrix}
\lambda_V(t_f)-V(t_f)\\
\lambda_\gamma(t_f)\\
\lambda_h(t_f)-\mu/r_f^2
\end{bmatrix}=0.
```

Residuals and shooting variables are scaled by the predicted initial costate magnitudes, with
the legacy magnitudes `(500,80000,9)` as lower bounds. The Jacobian uses central relative finite
differences whose perturbation shrinks with the switching homotopy width. For the first outer
`k` sample, two initial guesses are attempted: costates obtained from a backward legacy pass
and the legacy numerical seed `(-500,80000,9)`. Later `k` samples first use the converged
costates from the nearest previous `k`, while still querying GRAM afresh. The solve is rejected
if the infinity norm of the scaled terminal residual is above `1e-3`, which is the supported
accuracy of the 1-second bang-bang shooting grid.

The discontinuous switching condition is

```math
\lambda_{switch}=\frac{2kmV}{|AC_{D,slope}\pi|},
\qquad
\alpha_{raw}=\begin{cases}
\alpha_{max},&\lambda_V\ge\lambda_{switch},\\
\alpha_{min},&\lambda_V<\lambda_{switch}.
\end{cases}
```

Newton continuation smooths only this discrete decision with

```math
w_{high}=\frac{1}{2}\left[1+\tanh\left(
\frac{\lambda_V-\lambda_{switch}}{\epsilon_s}\right)\right]
```

and successively reduces `epsilon_s` through
`(1,0.1,0.01)` times the `lambda_V` scale. This is numerical homotopy for the
shooting solve, not a heat-load prediction multiplier. Once converged, switch extraction uses
the unsmoothed bang-bang sign. If more than two crossings exist, only the first contiguous
minimum-AoA interval is retained, yielding exactly two switches. The resulting profile is then
replayed with the full inertial trajectory predictor and current atmosphere/aerodynamics to
evaluate heat load for the outer `k` root.

### Closed-form fixed-point method

For each candidate `k`, the closed-form path uses the Chebyshev `f1` and `f2` trajectory model
described above. It evaluates the costates backward from

```math
\lambda_V(t_f)=V(t_f),\qquad
\lambda_\gamma(t_f)=0,\qquad
\lambda_h(t_f)=\frac{\mu}{(R_p+h_f)^2},
```

applies the bang-bang switching condition, retains the first contiguous low interval, applies
`alpha_high` everywhere else, and regenerates the closed-form trajectory. Iteration stops when
the absolute predicted heat-load change is at most `1e-3 J/cm^2` or after 100 iterations.
There are no empirical switch-time shifts or prediction multipliers.

Both solvers use the legacy rectangular heat-load quadrature

```math
Q_{pred}=\left(\sum_{j=1}^{N}\dot Q_j\right)\frac{t_f}{N}.
```

### Outer `k` root and switch extraction

If the heat-load limit is nonfinite or nonpositive, switches are `(Inf,Inf)`. Candidate root
brackets are:

| Solver | Mars | Venus | Earth | Titan/other |
| --- | --- | --- | --- | --- |
| Numerical root variable | 3.3 to 10 | 0.1 to 100 | 1 to 10 | 1e-5 to 10 |
| Closed-form root variable | 0 to 10 | 0 to 100 | 0 to 10 | 0 to 10 |

For numerical prediction, the root variable is divided by 100 before entering the costates.
The root residual always targets the configured physical limit directly:

```math
F(k)=Q_{current}+Q_{predicted}(k)-Q_{limit}.
```

No `0.92` prediction factor, switch-time scale factor, or balanced-window fallback is present.
The decision flow is:

1. Evaluate both configured bracket endpoints.
2. If their residuals have opposite signs, use Brent and re-evaluate the converged root to
   store its trajectory and profile. Closed form uses relative tolerance `1e-5`; the numerical
   TPBVP uses `1e-3`, matching its 1-second switch grid.
3. If the high-end residual is negative, return `(t,t)` because the requested heat load cannot
   be reached even at that endpoint.
4. If the low-end residual is positive, use the legacy unreachable-limit fallback: half the
   predicted passage for closed form or 1000 s for numerical integration.
5. Otherwise return `(Inf,Inf)`.
6. Convert only the first contiguous low-AoA interval with at least two samples into
   `(t_s1,t_s2)`; an invalid interval returns `(Inf,Inf)`.

### Live heat-load scheduling and reset

Scheduling runs only in selected maximum energy depletion with the heat-load submode. At drag
entry the state is initialized, but planning waits until the vehicle is 500 m below the entry
interface so the predictor has a meaningful in-atmosphere state. Failed plans may retry no
faster than once per second; a finite ordered pair marks the passage solved.

The numerical shooting plan is calculated once for that drag passage. It is not re-planned and
has no online correction factor. Closed form may retain the configured legacy second-switch
re-evaluation and security behavior. Closed-form re-evaluation keeps switch 1 fixed, roots
switch 2 over `[t,t+200 s]`, and targets the physical heat-load limit with no time scaling.
Nonfinite endpoints, an unbracketed residual, or a Brent failure preserve the existing pair.

Outside the drag passage, active-passage state is reset: the solved and passage flags become
false, switches become `(Inf,Inf)`, entry time becomes `NaN`, last re-evaluation becomes
`-Inf`, security clears, and previous heat load becomes zero.

### Closed-form security mode

Security is considered only as an `elseif` after re-evaluation, and only when all are true:

- security mode is enabled;
- solver is closed form;
- security is not already active;
- current heat load exceeds 98% of the limit;
- heat-load increase since the previous 10 Hz callback is below 2 J/cm2.

It predicts an all-minimum-AoA closed-form remainder. The current implementation explicitly
sets predicted heat-rate samples with `track.time <= planet.T` to zero; `planet.T` is the
planet temperature value numerically interpreted as seconds in this comparison. Remaining
load then uses the rectangular expression. If current plus remaining load exceeds the limit,
the switch pair becomes `(current_time,predicted_outbound_end)` and security remains active for
the rest of the passage. This is the physical two-transition interval reported by telemetry.

## Energy Bracketing and Targeting

### Reachable endpoints

The low-drag endpoint calls the targeting predictor with switch time equal to current time, so
the entire predicted passage uses minimum AoA. The maximum-energy-depletion endpoint solves
heat-load switches if that submode is enabled, then predicts the same heat-load, heat-rate, and
structural control combination used by maximum-energy-depletion mode. Without heat-load
control, switches are `(Inf,Inf)` and high phases are constrained only by selected instantaneous
submodes.

Both trajectories extend predicted passage duration by 120 s, capped by the planning horizon,
and use a grid with nominal maximum 1 s spacing, at least 64 points, and at most 20,000 points.
Outbound exit is the first sample at or above EI after the minimum-altitude sample. If none is
found, the final sample is used.

### Target-energy construction

Specific orbital energy and bound-orbit radii are

```math
epsilon=\frac{V^2}{2}-\frac{\mu}{r},
\qquad r_p=a(1-e),
\qquad r_a=a(1+e).
```

Nonnegative or nonfinite energy returns `NaN` periapsis and infinite apoapsis. Invalid orbital
elements return `NaN` radii.

The desired two-body energy for target apoapsis is

```math
epsilon_{2body,target}=-\frac{\mu}{r_{a,target}+r_p}.
```

Invalid target or periapsis radius returns `NaN`.

To reproduce the deprecated target-planning sequence, guidance performs two vacuum
central/J2/N-body propagations with 1 s RK4:

1. Propagate from current entry state to outbound EI, with a 2000 s cap. Periapsis is considered
   passed once radial velocity is nonnegative. Stop when outbound altitude reaches EI.
2. Propagate that exit state toward apoapsis. The time cap is the larger of 20,000 s and 1.25
   times the osculating two-body half-period, so high-apoapsis cases such as VEX can reach the
   event. Stop when radial velocity crosses from positive to nonpositive.

The gravity predictor starts with the configured inverse-square or inverse-square/J2 model and
adds every configured `NBodyGravityModel`. Third-body positions come from the simulation's
ephemeris cache when available and otherwise use the same memoized SPICE fallback as the RHS.
The second propagation returns the resulting change in two-body specific energy. If its
periapsis is finite, that periapsis is used; otherwise the low-drag predicted periapsis is used.
The atmospheric-exit target is

```math
epsilon_{target}=epsilon_{2body,target}-Delta epsilon_{pert}.
```

A positive perturbation energy change therefore makes the required exit target more negative.
The same configured N-body acceleration is also included in low-drag, maximum-depletion,
heat-load shooting, and targeting switch predictions.

### Targeting trajectory and force model

For a candidate absolute switch time, the base profile is maximum before the switch and
minimum at and after it. At minimum, constraints return immediately and cannot raise AoA. At
maximum, enabled heat-rate and structural controls are evaluated and their minimum is used.
The predictor's heat-rate root uses the configured minimum directly; the early minimum-AoA
return still protects an exact zero command.

Aerodynamic directions in planet-fixed coordinates are

```math
\hat d=-\frac{\mathbf V_{rel}}{V},
\qquad
\hat l=\operatorname{normalize}(\hat h\times\hat V_{rel}),
\qquad
\mathbf F=qA(C_D\hat d+C_L\hat l).
```

Force is rotated to inertial and divided by mass. A degenerate lift vector becomes zero lift.
RK4 holds the sample's AoA constant across its four substeps.

### One-switch root

Targeting switch planning runs in control only when selected mode is targeting, targeting is
active, the vehicle is inside the drag passage, and the stored switch is not already finite.
Thus it runs once. A nonfinite target returns an infinite switch.

For delay `tau`, define

```math
F(\tau)=\frac{epsilon_{exit}(\tau)-epsilon_{target}}{10^6}.
```

The implemented root interval is `[0,2000 s]`. Both endpoint predictions overwrite the stored
energy bracket with their extrema. Delay selection is:

1. If the zero-delay residual is exactly zero, choose zero.
2. Else, if the 2000 s residual is exactly zero, choose 2000 s.
3. Else, if both are finite and have opposite signs, solve with Brent using `rtol=atol=1e-8`.
4. Otherwise choose whichever endpoint has smaller absolute residual.

Unlike the instantaneous constraint roots and second-switch correction, this Brent call is not
wrapped in `try/catch`; a root-solver exception propagates to the caller.

For an interior solution, evaluate residuals one second to either side. If both are finite and
unequal, linearly interpolate a zero and clamp it between those neighboring delays. The stored
absolute switch is `current_time + delay`.

The targeting-switch function accepts `heat_load_j_cm2` for interface consistency, but the
current one-switch residual does not read it. Targeting itself has no two-switch heat-load
phase; heat-load control enters only the maximum-energy endpoint used by guidance bracketing.

After the switch becomes finite, all later control callbacks skip targeting planning. The
current targeting state and switch are not reset at outbound EI.

`_edg_targeting_switch_outcomes` is a support helper that predicts immediate low drag and a
high-drag switch one second after the low-drag predicted duration. It is not called by the live
guidance or control path.

## Present but Inactive Helpers

The source contains three additional helpers that do not affect the current callback flow:

- `_edg_set_targeting_fallback!` assigns the same unreachable-target fallback that guidance
  currently assigns inline.
- `_edg_sample_prediction_atmosphere` samples density and temperature at zero latitude and
  longitude after clamping altitude to zero, but no current predictor calls it.
- `_edg_padded_heat_load_window` clamps a finite switch window to track start and end, but no
  current switch solver calls it. A nonfinite window would be returned unchanged.

`_edg_targeting_switch_outcomes`, described above, is likewise not called by the live path.

## Solar-Panel Actuation

The AoA effector defaults to links 2 and 3. For every configured link it enforces these checks:

1. Index must be within the spacecraft link array; otherwise throw `ArgumentError`.
2. Link must not be the root; otherwise throw `ArgumentError`.
3. Rotation axis is `abs(link.r)`.
4. Apply the absolute link rotation `rotate_link(link,axis,pi/2-alpha)`.
5. Store the command in `link.α`.

Only configured links rotate; the bus remains unchanged. Both the combined energy-depletion
effector and the panel-AoA effector return exactly zero force and torque. Aerodynamic forces are
generated later by the normal dynamics effector from the updated orientations.

## Stored Telemetry

After every control command, state records final AoA, heat-rate candidate AoA, structural
candidate AoA, controlled heat rate, maximum controlled-link heat load, dynamic pressure, and
equivalent structural load. The Odyssey CSV additionally saves mode code, targeting flag and
switch, target and bracket energies, both heat-load switches, security state, and both panel
AoAs.

For targeting, guidance first propagates a vacuum trajectory from atmospheric exit to the next
apoapsis with the configured central/J2 and N-body gravity models. If the resulting change in
two-body specific orbital energy is `Delta epsilon_pert`, the drag-passage target is

```math
\epsilon_{target,exit}=-\frac{\mu}{r_{a,target}+r_p}-\Delta\epsilon_{pert}.
```

The correction and bracket are calculated once at drag-passage entry. N-body positions use the
same ephemeris cache and SPICE fallback as the simulation force model. The correction is stored
in `target_perturbation_energy_change_jkg` for validation.

Mode codes are `2` for targeting, `1` for maximum energy depletion, `-1` for safe low drag,
and `0` for inactive. Controlled-panel heat telemetry is intentionally separate from the
standard all-link maximum because the bus is not actuated.

## Current-Code Integration Boundary

Relative wind is spacecraft planet-fixed velocity minus atmospheric wind. Some deprecated
predictor paths used the opposite sign. Numerical full-trajectory replay uses the current
per-link aerodynamic model and configured planet-fixed gravity model so the predicted switch
profile follows the unchanged current RHS. The closed-form option intentionally retains the
deprecated aggregate flat-plate model and updated Chebyshev fits.

## Odyssey Configuration and Validation

The validation example uses the deprecated Odyssey values: 411 kg dry spacecraft mass, 50 kg
propellant, a `2.2 x 2.6 x 1.7 m` bus, two half-panels forming a `3.76 x 1.93 m` planform,
reflection coefficient 0.9, 9800 km initial apoapsis radius, 100 km periapsis altitude, 93.6 deg
inclination, and the 2001-11-06 19:00 epoch. Limits are 0.15 W/cm2 heat rate, 30 J/cm2 heat load,
and 0.5 Pa equivalent structural load.

```bash
SPACEAGORA_EDG_PHASE=constraints julia --project=. examples/AGORA_Odyssey_Energy_Depletion_Control_Test.jl
SPACEAGORA_EDG_PHASE=heat_load SPACEAGORA_EDG_HEAT_LOAD_SOLVER=closed_form julia --project=. examples/AGORA_Odyssey_Energy_Depletion_Control_Test.jl
SPACEAGORA_EDG_PHASE=heat_load SPACEAGORA_EDG_HEAT_LOAD_SOLVER=tpbvp_integration julia --project=. examples/AGORA_Odyssey_Energy_Depletion_Control_Test.jl
SPACEAGORA_EDG_PHASE=targeting julia --project=. examples/AGORA_Odyssey_Energy_Depletion_Control_Test.jl
```

| Case | Maximum heat rate | Maximum heat load | Structural load | Switch/target result |
| --- | ---: | ---: | ---: | --- |
| Constraints only | 0.149990 W/cm2 | not controlled | 0.45767 Pa | no low-AoA interval |
| Heat load, closed form | 0.149990 W/cm2 | 30.04747 J/cm2 | 0.45724 Pa | 8314.000 s, 8498.735 s; security active |
| Heat load, numerical shooting with live GRAM | validation incomplete | validation incomplete | validation incomplete | one-time planner exceeded 40 minutes |
| Targeting | 0.101729 W/cm2 | 4.18560 J/cm2 | 0.26411 Pa | one switch at 8131.706 s; 9754.927 km for 9755 km target |

All CSVs use 1 s sample spacing. Each heat-load case contains exactly one contiguous low-AoA
interval. The closed-form run reaches 30.04747 J/cm2, 0.16% above the 30 J/cm2 limit, without
any switch-time or heat-load factor. The previous numerical result used the obsolete frozen
atmosphere track and is intentionally no longer reported. A full live-GRAM Odyssey shooting
run was stopped after 40 minutes while still evaluating the one-time planner, so this table
does not claim an unverified numerical result. The targeting apoapsis miss in the validated
run is approximately 72.66 m.

### Venus N-body validation

`examples/AGORA_Vex_Energy_Depletion_Control_Test.jl` uses the deprecated VEX control-test
epoch, orbit, mass, geometry, limits, Venus GRAM atmosphere, inverse-square/J2 gravity, and Sun
`NBodyGravityModel`. The legacy spacecraft reference areas are 5.74 m2 for the bus and 5.70 m2
for the complete solar array. The target can be changed from the run file with
`SPACEAGORA_VEX_TARGET_RA_M`.

The original 72,480 km target remains just above the current reachable interval. Its corrected
target energy is `-4.129880367e6 J/kg`, while the bracket is
`[-4.133222565e6,-4.130329559e6] J/kg`; guidance therefore selects maximum energy depletion.
This case exposed and corrected an earlier VEX example bus-area mismatch.

The reachable validation command is

```bash
SPACEAGORA_VEX_TARGET_RA_M=72450000 julia --project=. examples/AGORA_Vex_Energy_Depletion_Control_Test.jl
```

For that run, the exit-to-apoapsis J2/N-body energy change is `95.778186 J/kg`. Targeting is
active, the single switch is at `43104.264115 s`, and the final osculating apoapsis radius is
`72450.407647 km`, 407.647 m above the 72,450 km target.
