# Lunar Landing

A powered descent from orbit to the surface, flown 6-DOF: an Apollo-style
quadratic guidance commands the descent engine's thrust and the attitude that
points it, an RCS attitude controller tracks that attitude, and the terrain
comes from digital elevation maps (DEMs) that serve as the radar altimeter,
stop the run at touchdown, and shape the surface the viewer draws.

`scripts/dev/viewer_demos/apollo11_landing.jl` flies Apollo 11 from powered
descent initiation (PDI, 1969-07-20 20:05:05 UTC, 15.24 km above Tranquility
Base and 480 km uprange on the descent orbit) to touchdown with NASA's lunar
module model and the LROC terrain of the site.

## Terrain models

`DEMTerrainModel` holds one or more `DEMGrid`s (regular latitude/longitude
grids of heights above the body's reference sphere, finest first) and answers
`terrain_height(model, lat_deg, lon_deg)` by bilinear interpolation in the
first grid that covers the point; `NoTerrainModel()` is the flat sphere.
`load_site_terrain("data/terrain/moon/apollo11/site.json")` reads a site
directory written by `scripts/dev/terrain/fetch_moon_site.py`, which fetches:

- the LOLA 128 pixel-per-degree global grid (237 m/px) in a window around the
  site, by HTTP range requests straight out of the PDS raster;
- the LROC NAC digital terrain model of the site (2 m/px, resampled to 4 m)
  where one exists (Apollo 11 does);
- nested imagery patches from NASA Moon Trek, from the LROC WAC global mosaic
  (83 m/px over a 4° window) down to the 26 cm NAC mosaic (0.65 m/px over the
  last 500 m).

```
python3 scripts/dev/terrain/fetch_moon_site.py --site 0.67416 23.47314 --out data/terrain/moon/apollo11
```

It also writes an albedo-normalized copy of every imagery level
(`level_k_albedo.jpg`), which is what the page drapes when it lights the ground
itself; `--reuse DIR` derives a second copy of a site from an existing one
without touching the network, and the demo draws the copy that `--site-json
PATH` or `SPACEAGORA_TERRAIN_SITE` names. See
[Visualization](visualization.md) for what the page does with them.

The data stay under `data/` (not tracked). Grids are little-endian Float32
files with a JSON header, so other DEMs can be dropped in the same layout.

## Powered descent guidance

`ApolloDescentGuidanceModel` implements the LM guidance computer's quadratic
law (Klumpp, *Apollo Lunar Descent Guidance*, 1974) in a site-fixed frame
(x uprange, y crossrange, z up; velocities relative to the rotating body):

```math
a_c = a_T - \frac{6\,(v_T + v)}{T} + \frac{12\,(r_T - r)}{T^2}
```

with the time-to-go `T` re-solved every cycle from the uprange cubic that adds
a jerk target. Two quadratic phases follow the Apollo programs, P63 braking
(targets beyond the high gate, the phase ends at the high-gate altitude) and
P64 approach (a hover point above the site), then a P66 rate-of-descent phase
nulls horizontal velocity and descends at 1 m/s, 0.5 m/s below 10 m, until
the engine cutoff altitude. `apollo11_descent_targets()` returns targets in
the spirit of Apollo 11's; `DescentPhaseTargets` lets you set your own.

The thrust command is `m |a_c - g|` with the rotating-frame terms included;
the descent engine's throttle follows the DPS envelope (`throttle_min` to
`throttle_max`, 10 to 60 percent, and full thrust above it). The attitude
command points body `-z` (the engine axis) along the thrust with the windows
(body `+x`) along the projection of up plus the flight direction: windows up
during braking, facing the site once the vehicle pitches back.

`ApolloDescentControlModel` throttles the engine toward the command with a
slew limit, runs a rate-limited attitude loop whose torques are bounded by
the RCS authority, applies the thrust along the vehicle's actual engine axis,
and drains propellant through both. The two effectors share an
`ApolloDescentState` (phase, time-to-go, thrust, commanded attitude, radar
altitude, site-frame state, touchdown record). Run the simulation with
`isolate_state=false` when you want to read that state afterwards; the
scenario also saves thrust, throttle, time-to-go, radar altitude, phase and
attitude error as result columns.

When a control effector reports a landing (`touchdown_spec`), the engine
replaces its impact event with a touchdown event on the terrain: the run
ends when the vehicle's reference point reaches `touchdown_height_m` above
the ground and the touchdown time, ground-relative velocity and miss distance
are recorded.

```julia
terrain, site = load_site_terrain("data/terrain/moon/apollo11/site.json")
braking, approach = apollo11_descent_targets()
gcfg = ApolloDescentConfig(site_lat_deg=site.lat_deg, site_lon_deg=site.lon_deg, approach_azimuth_deg=270.0,
    braking=braking, approach=approach)
state = ApolloDescentState(1)
guidance = ApolloDescentGuidanceModel(gcfg, state, terrain)
control = ApolloDescentControlModel(ApolloDescentControlConfig(touchdown_height_m=2.6), gcfg, state, terrain)
# ... GuidanceModel(guidance_effectors=(guidance,), guidance_rates=[2.0]),
#     ControlModel(control_effectors=(control,), control_rates=[0.05]), orientation_sim=true
```

The control effector also reports what its thrusters are doing, through
`control_thruster_levels`: the descent engine's firing level is its actual
thrust over the thruster's rating, and the sixteen RCS jets' levels come from
a least-norm allocation of the commanded body torque over their torque arms
about the spacecraft reference point, clipped into 0 to 1 and refreshed every
control cycle. The run writes them as `sc1_thruster_level_1..17` and the
viewer draws a plume on each firing thruster, so the descent engine burns
throughout the braking phase and the jets puff at the phase handovers. A
steady phase such as P66 commands a small fraction of the RCS authority, so
those puffs are small; the levels are in the selection panel and plot like
any other quantity.

## Plume-surface interaction

`PlumeSurfaceInteractionModel` is a dynamic effector for what the descent
engine's exhaust does to the ground: the pressure and shear stress it lays
down, the regolith it erodes and by which mechanism, the crater it digs, where
the grains land, and the small thrust augmentation the reflected plume gives the
vehicle in the last couple of nozzle diameters. It is constructed with the descent control effector
(whose actuator state carries the engine's actual thrust) and the terrain,
so it always sees the throttle the controller flew and the ground the
altimeter measured:

```julia
control = ApolloDescentControlModel(ApolloDescentControlConfig(), gcfg, state, terrain)
plume = PlumeSurfaceInteractionModel(control, terrain)
# ... dynamic_effectors = (gravity..., plume)
```

`PlumeSurfaceConfig` carries the engine (1.5 m exit diameter, area ratio 47.5,
7.2 bar chamber pressure, 4.5 to 45 kN throttle band — the Apollo LM's descent
propulsion system by default), the soil (a `RegolithProperties`, lunar mare by
default), the erosion regimes evaluated on it, the plume field the surface state
is read from, and the crater and ejecta settings.

**What is computed and what is calibrated.** Everything below comes from a
sourced model or is derived in the file that owns it, except for the rows marked
ASSUMPTION, every one of which is a named configuration field with its reasoning
on the field rather than a buried literal:

| Quantity | Where it comes from |
|---|---|
| surface pressure, wall shear, gas density, temperature, Mach | the plume field (a Gaussian closure, or the Simons source flow of the nozzle) |
| erosion rate and the active regime | `regolith_erosion.jl`, every coefficient sourced — except the soil's `saltation_efficiency` |
| erosion threshold | Shao and Lu (2000) applied to the Lunar Sourcebook's soil; **derived, not fitted** |
| crater depth and radius | the regimes' rate integrated in time on a radial grid, over the soil's in-situ bulk density |
| ejecta speed, angle, deposition radius, escape fraction | `ejecta_transport.jl`; the angle is an **input** from the Apollo films |
| **`soil.saltation_efficiency` = 10** | ASSUMPTION: the saltation cascade, unsourced. The best fit to Lane and Metzger's eleven Apollo 12 altitudes is 17.5 (analytic field) and 13.7 (table); it is left at 10, see the validation study for why |
| **`gas_residence_time_s` = 1 s** | ASSUMPTION: sets the pore-pressure diffusion depth of the diffusion-driven-flow regime |
| **the ground-effect correlation** | ASSUMPTION: 3 percent of thrust at contact over two exit diameters; the literature has not characterized the near-ground base-pressure rise, and this is the only force the model applies to the vehicle |
| **`crater_edge_fraction`, `crater_min_depth_m`** | DEFINITIONS, not measurements: where the crater's edge is taken to be. They affect the reported radius and nothing else |
| **the ejecta transport free parameters** | ASSUMPTIONS: the entrainment length, the wall-jet thickness, the collision diameter and the lognormal shape of the grain-mass distribution; tabulated in "What is assumed" under Ejecta transport below |

The two constants the previous model calibrated — `threshold_shear_pa = 0.15 Pa`
and `erosion_efficiency = 10` — are still on the configuration but are used only
by `erosion_model = :roberts_fitted`.

### The model

The model is assembled from four separable pieces, each of which is its own
file and each of which can be swapped or evaluated on its own.

**The gas on the ground** comes from a *plume field* that `PlumeSurfaceConfig`
carries in its `field`. `PlumeAnalyticField` is the Gaussian footprint described
below and is the default, so a run that asks for nothing else keeps the gas
dynamics this model has always had. `PlumeFieldTable` is a plume computed from
the nozzle and tabulated; it is described in "The tabulated plume field" below.

The plume's momentum-carrying core spreads over a footprint of radius
`R_p = h tan(θ_p)` at the nozzle exit's height `h` above the ground, with a
Gaussian surface pressure normalized so all of the engine's axial momentum is
turned by the surface:

```math
p(r) = p_0 \exp\!\left(-\left(\tfrac{r}{R_p}\right)^2\right), \qquad p_0 = \frac{F}{\pi R_p^2}
```

The wall shear stress is the skin-friction coefficient times that pressure,
vanishing at the stagnation point and peaking at `r = R_p/\sqrt{2}`, as an
impinging jet's does:

```math
\tau(r) = c_f\, p_0\, \frac{2r}{R_p} \exp\!\left(-\left(\tfrac{r}{R_p}\right)^2\right)
```

**The soil moved at each radius** comes from the erosion regimes of
`regolith_erosion.jl`, documented in their own section below. The effector
evaluates `regolith_erosion_rate` at 128 radii across the footprint on the gas
state the field reports there, integrates the result over the ground, and
records which regime moved the most mass as `sc{i}_plume_regime`. Nothing in
that path is fitted inside this file; every coefficient belongs to
`RegolithProperties` and is sourced there.

Which regimes are evaluated is `PlumeSurfaceConfig.regimes`, which defaults to
`plume_default_regimes()` — viscous erosion in Roberts' shear-excess form with
the *derived* Shields threshold, plus diffusion-driven flow and bearing-capacity
failure. That differs in one place from `default_erosion_regimes()` in the
regolith module, whose viscous closure is Metzger's energy-flux law, and the
reason is a measurement rather than a preference:

!!! warning "Metzger's energy-flux law cannot be driven from these plume fields"
    `ViscousErosionEnergyFlux` computes `E = 3τ²/(ρ_s v̄)` and compares it with
    `E_th = 0.123 J/(m² s)`, a threshold fitted to the Apollo 16 landing video
    **through Metzger's own plume model**. Evaluated on the surface gas state
    either of this repository's fields reports, that threshold is crossed by
    three to four orders of magnitude: at the Apollo 12 approach thrust the
    energy flux at the shear peak is 4.3 W/m² at 31.9 m and 41 W/m² at 10 m, and
    integrating the resulting local rate over the footprint gives about
    4 × 10⁴ kg/s on the analytic field and 6 × 10⁶ kg/s on the table, at every
    height, against the 11 to 99 kg/s Lane and Metzger measured. The rate also
    comes out nearly independent of height, exactly the degeneracy that law's
    own docstring predicts for a density that tracks the pressure. `E` needs the
    density in the laminar sublayer at the wall; what both fields report is the
    post-shock free-stream density. Until a field computes the sublayer state,
    the law and the fields do not meet, and bolting them together would report a
    number wrong by 10⁴. Pass `regimes=default_erosion_regimes()` to run it
    anyway.

The old single law is still reachable as `erosion_model = :roberts_fitted`, so
the two can be compared on one gas state:

```julia
PlumeSurfaceInteractionModel(control, terrain;
    config=PlumeSurfaceConfig(erosion_model=:roberts_fitted))
```

It is the shear in excess of a fitted threshold, closed as a momentum balance
over the annulus where the excess is positive:

```math
\dot m = \frac{\eta}{v_{ej}} \int \max(\tau(r) - \tau_t,\, 0)\, \mathrm{d}A
```

with `τ_t = threshold_shear_pa = 0.15 Pa` chosen to put the onset at 31 m and
`η = erosion_efficiency = 10` standing in for the saltation cascade. **Both
constants are used by that law only.** The default law ignores them; its
threshold is `shields_threshold_shear_pa(soil, g)`, 0.057 Pa on lunar mare at
1.625 m/s², which is a property of the soil and does not move when the plume
field changes.

Both laws report the same `ejecta_mps`: one grain dragged across one footprint
radius from rest by the gas dynamic pressure at the shear peak. The cumulative
eroded mass is the time integral of the rate, advanced only when the solver
moves past every time already integrated, so rejected steps neither double count
nor run it backwards.

### The crater

The effector keeps a radial grid of eroded depth under the vehicle and advances
it with the same monotone time guard as the cumulative mass: at each accepted
time it evaluates the regimes' local mass flux at every node, trapezoidally
integrates it, and divides by the soil's in-situ bulk density, so the depth is
the depth of the hole rather than the loose volume thrown out of it.

The grid is logarithmic in radius, 64 nodes from 5 cm to 40 m by default,
because the crater is toroidal: its inner wall sits inside the first meter while
its outer edge reaches tens of meters, and a uniform grid that covers the second
cannot resolve the first. `plume_crater_profile(model, i)` returns the radii and
the live depths; `sc{i}_plume_crater_depth_m` is the deepest point and
`sc{i}_plume_crater_radius_m` the outermost radius still at
`crater_edge_fraction` (a tenth) of it.

The depth under the stagnation point feeds back into the height the plume field
is queried at, so a deepening crater moves the ground away from the nozzle. On
an Apollo descent that feedback is 1.6 cm against a 50 m onset height and moves
nothing measurable; it exists because the model would otherwise be inconsistent,
not because it matters here. Turn it off with `crater_height_feedback=false`.

Where the crater stops is a definition rather than a measurement, and it is two
configuration fields: the edge is the outermost radius still at
`crater_edge_fraction` (a tenth) of the deepest point *and* at least
`crater_min_depth_m` (1 mm, about fourteen median grain diameters) deep. The
floor matters: a pure ratio test returns a wide radius the instant the plume
touches the ground, when the surface has lost a few grain layers everywhere and
nothing that could be called a crater exists, and the viewer would draw a scour
mark wider than the region the plume is eroding. Below the floor no radius is
reported at all, and neither field changes any other quantity.

!!! warning "The crater assumes a stationary vehicle, and Apollo 11 was not one"
    The grid is axisymmetric about the *current* impingement point and does not
    follow it as the vehicle translates, so every second of erosion is credited
    to the same patch of ground. On the Apollo 11 demo the vehicle covers 310 m
    of horizontal distance below the erosion onset height, and 63 percent of the
    eroded mass is laid down over the last 92 m of it — against a crater whose
    own edge radius is 4.3 m. The profile the model reports is therefore the
    crater a **hovering** vehicle would dig over the same erosion history, not
    the trench a translating one leaves, and the depth is an upper bound for any
    vehicle that moves. That cuts the same way in the validation study: the
    3.2 cm it scores against Metzger's 6 cm is an upper bound under-predicting a
    measurement, so the real gap is wider than the ratio says. Making the crater
    follow the ground needs a two-dimensional grid in the planet frame, which
    this model does not have.

### Ejecta as a per-sample diagnostic

The grain-transport model of `ejecta_transport.jl` is evaluated once per **saved
sample**, never on the right-hand side: no force depends on it, and one
trajectory integration per grain size and launch radius is far too expensive for
a solver stage. `plume_refresh_ejecta!` recomputes it at most once per sample
however many columns read it, and not at all when nothing is eroding or when
`ejecta_diagnostic=false`.

Two weightings separate this from a sweep of the shear profile, and both were
the fix the ejecta work identified for its kilometer-scale deposition radii:

- each launch radius contributes the **regimes' own local erosion rate** there
  times its annulus area, instead of the raw wall shear stress with no
  threshold;
- each grain size carries the mass a **lognormal fitted to the soil's own `D50`
  and `D84/D50`** gives it (`ejecta_lognormal_mass_weights`). A lognormal has
  `ln D84 − ln D50 = σ`, so the one sourced ratio fixes the whole width; with
  the lunar values the 70 µm median carries 51 percent of the mass and the 5 µm
  bin under 1 percent, where equal weights per bin would have given it a fifth.
  The lognormal *shape* is an assumption; the two percentiles it is fitted to
  are sourced.

Together they halve the deposition radius and no more: at the Apollo approach
thrust 10 m up the mass-weighted mean goes from 1553 m to 755 m and the mean
speed from 143 to 107 m/s. The rest is not an artifact — a 70 µm grain leaving
at 107 m/s and 2° above the horizontal has a ballistic range of 424 m in lunar
gravity. See the validation study for why that is not the same quantity as the
visible blast zone.

`plume_ejecta_summary(model, i)` returns the whole distribution, histograms
included. Three scalars are saved as columns: the mass-weighted mean ejection
angle, the mass-weighted mean deposition radius and the mass fraction leaving
faster than escape speed.

### The ground effect

The force the effector returns is the ground-effect thrust augmentation along
the engine axis, and no torque. The correlation is an exponential in the height
over the nozzle exit diameter `D_e`,

```math
\Delta F = f_{max} F\, \frac{e^{-x/s} - e^{-x_0/s}}{1 - e^{-x_0/s}}, \qquad x = h/D_e
```

zero at and above `x_0 = 2` exit diameters and rising monotonically to
`f_{max} = 3` percent of the thrust at contact. The exponential form follows the
decay of the base-pressure rise measured for nozzles near a plate; the
magnitude is a modeling choice sized to the few-percent effect reported for
lunar-lander-class plumes, not an Apollo flight measurement. Change it through
`ground_effect_max_fraction`, `ground_effect_scale` and `ground_effect_cutoff`.

This is the only wrench the model applies to the vehicle, and it is the one
quantity in the whole model for which the validation study could find no
published reference at all.

### The tabulated plume field

The Gaussian footprint is a closure, not a plume: its width is an assumed
half-angle and its shear stress an assumed fraction of the static pressure.
`PlumeFieldTable` replaces both with a plume computed from the engine. Ask for
one by loading it into the configuration:

```julia
cfg = PlumeSurfaceConfig(field=load_plume_field("data/psi/apollo_lmde.json"))
plume = PlumeSurfaceInteractionModel(control, terrain; config=cfg)
```

The plume itself is the Simons source-flow model of a rocket exhausting into
vacuum with Boynton's nozzle-boundary-layer correction (G. A. Simons, *Effect of
nozzle boundary layers on rocket exhaust plumes*, AIAA Journal 10(11), 1972;
F. P. Boynton, *Exhaust plumes from nozzles with wall boundary layers*,
J. Spacecraft and Rockets 5(10), 1968), in the explicit form of section 2.1 of
P. J. Herráiz, J. M. Fernández and J. R. Villa, *Development of a MATLAB plume
impingement tool for fast system analysis*, EUCASS 2019, DOI
10.13009/EUCASS2019-661. It meets the ground through classical Newtonian impact
theory for the wall pressure and the oblique normal-shock relations for the gas
behind the surface shock, and the wall shear stress is the wall jet's dynamic
pressure times the rough-surface drag coefficient of 0.2 that Roberts used for
the meteorite-gardened lunar surface. Classical Newtonian is chosen over
modified Newtonian because it conserves the plume's axial momentum on the plane
exactly: the built table's surface pressure integrates to 1.008 times the engine
thrust at every height.

Every engine parameter is sourced or derived from a sourced number, and the
three that are neither — the nozzle exit lip half-angle, the divergent length
and the mean limit-speed ratio of the boundary-layer decay — are named fields of
`PlumeNozzle` marked as assumptions there. `data/psi/README.md` tabulates all of
them with their sources, and `scripts/dev/psi/build_plume_field.jl` rebuilds the
tables.

What changes when a table is used, at the Apollo 11 approach thrust of 11.5 kN:

| | analytic | table | ratio |
|---|---|---|---|
| stagnation pressure at 10 m | 168 Pa | 205 Pa | 1.22 |
| footprint radius at 10 m | 4.66 m | 4.40 m | 0.94 |
| peak wall shear at 10 m | 1.44 Pa | 9.18 Pa | 6.4 |
| erosion onset height | 50.3 m | 130 m | 2.6 |

The pressure and the footprint agree to about 20 percent over the whole descent,
which is the real verdict on the Gaussian footprint: as a *pressure* closure it
was close. The shear stress does not agree, because Roberts' drag coefficient
acting on the wall jet is a far stronger coupling than a 0.01 skin-friction
coefficient acting on the static pressure, and the erosion onset follows it.

!!! note "The calibration that did not transfer, and no longer has to"
    With the old fitted law, swapping the analytic field for the table moved the
    Apollo 11 descent from 429 kg of eroded soil to 29 tonnes, a factor of 68,
    because `threshold_shear_pa` and `erosion_efficiency` had been calibrated
    against the analytic field's much weaker wall shear. **The default law does
    not use either constant**, and its threshold is a property of the soil rather
    than of the field, so the same soil erodes under both plumes: over Lane and
    Metzger's Apollo 12 descent profile the two fields now give 1200 kg and
    1396 kg, a ratio of 1.16 instead of 27. Nothing was recalibrated to achieve
    that. The one constant that is still unsourced is the soil's
    `saltation_efficiency` (10); the validation study's README records the value
    that best fits Lane and Metzger's eleven Apollo 12 altitudes (17.5 on the
    analytic field, 13.7 on the table) and why it was left alone.

One published number for this engine at a stated condition is a
peak laminar smooth-wall shear stress of 92 Pa for the LMDE hovering 5 m above
the surface at 13.34 kN, from the DSMC simulations of A. B. Morris, *Simulation
of rocket plume impingement and dust dispersal on the lunar surface*, PhD
dissertation, University of Texas at Austin, 2012, section 4.4.2. The shipped
table gives 41.4 Pa there, 2.2 times low. That is the deficiency the same work
diagnoses in section 4.6: a point source at the nozzle exit plane under-predicts
the surface stress below about ten nozzle diameters, because the plume really
originates about 1.8 m further downstream, where the nozzle's internal
compression wave meets the axis. `data/psi/apollo_lmde_virtual_source.json` is
the same table built with that measured offset; it gives 100.0 Pa at the same
condition, 9 percent above the DSMC value, and has to freeze within about 2.6 m
of the ground where the source would reach the surface. The conservative one
ships as the default, and neither was tuned to match the 92 Pa.

### The plume field against the Apollo measurements

Lane and Metzger, *Estimation of Apollo lunar dust transport using optical
extinction measurements*, Acta Geophysica 63(2), 2015, inverted the Apollo 12
descent films and put the radius of the eroding region at 1.3 to 2.2 times the
nozzle height, an effective half-angle of 52 to 66 degrees, against the 25
degrees `plume_half_angle_deg` assumes (both results are taken from the
validation study's reading of that paper in `benchmarks/studies/psi_validation/`
rather than checked against the primary source here). The computed field reproduces that
without widening the plume: at 11.5 kN its eroding radius is 1.9, 1.5 and 1.3
times the height at 5, 10 and 15 m, while its momentum-carrying pressure core
stays at 24 degrees. The two decouple because the wall shear stress is the skin
friction of the radial wall jet rather than a fraction of the local static
pressure, and the wall jet keeps accelerating outward while its density thins,
so the shear reaches far beyond the pressure footprint.

With the erosion regimes reading that field, `apollo12_scour_radius` passes on
both fields for the first time: 0.34 to 0.91 of the measured radius on the
analytic field and 0.57 to 1.76 on the table.

Use `plume_wall_shear`, `plume_mean_shear` and `plume_scour_radius` to compare
against measurements of a profile rather than of a peak; comparing a model peak
against an area-averaged measurement flatters the analytic field considerably.
`data/psi/README.md` gives the full comparison, the two published references
that disagree with each other by a factor of 25, and the erosion-onset residual
that is left after both are taken into account.

### What it records

The effector keeps a `PlumeSurfaceState` per spacecraft and
`default_save_fields` publishes it as fourteen result columns whenever the
effector is in the run. Ten are updated at every right-hand-side evaluation:

| Column | Meaning |
|---|---|
| `sc{i}_plume_height_m` | height above the terrain along the engine axis, from the vehicle's reference point |
| `sc{i}_plume_shear_pa` | peak wall shear stress under the plume |
| `sc{i}_plume_pressure_pa` | peak surface pressure under the plume |
| `sc{i}_plume_erosion_kg_s` | mass erosion rate, integrated over the ground |
| `sc{i}_plume_eroded_kg` | its time integral |
| `sc{i}_plume_regime` | erosion regime moving the most mass: 0 none, 1 viscous, 2 diffusion-driven flow, 3 bearing-capacity failure |
| `sc{i}_plume_erosion_radius_m` | outer edge of the region moving soil right now |
| `sc{i}_plume_crater_depth_m` | deepest point of the crater eroded so far |
| `sc{i}_plume_crater_radius_m` | its edge |
| `sc{i}_plume_ejecta_mps` | characteristic ejecta speed |
| `sc{i}_plume_ground_effect_n` | ground-effect thrust augmentation |

and three come from the per-sample ejecta diagnostic:

| Column | Meaning |
|---|---|
| `sc{i}_plume_ejecta_angle_deg` | mass-weighted mean ejection angle above the local horizontal |
| `sc{i}_plume_ejecta_range_m` | mass-weighted mean deposition radius |
| `sc{i}_plume_ejecta_escape_frac` | mass fraction leaving faster than escape speed |

The viewer's dust module reads all of them: the sheet's radius follows
`erosion_radius_m`, its elevation follows `ejecta_angle_deg`, the scour mark
follows `crater_radius_m` and the haze follows `ejecta_range_m`, so the three
numbers that used to be hard-coded in `viewer/src/dust.js` (a 1 to 3 degree
sheet, a 45 m sheet radius and a 14 m scour) are now model outputs.

Run with `isolate_state=false` to read `plume.state`, `plume_crater_profile` and
`plume_ejecta_summary` directly afterwards. `plume_quantities`,
`plume_surface_footprint`, `plume_erosion_onset_height` and
`plume_ground_effect_force` are the same model as plain functions, for sizing a
scenario without running one.

On the Apollo 11 descent (`scripts/dev/viewer_demos/apollo11_landing.jl`, the
analytic field and the shipped defaults) erosion begins 67 s before touchdown
with the engine 54.5 m above the ground, peaks at 42 kg/s, and moves 905 kg of
regolith in all — against Lane and Metzger's 2.6 t for Apollo 12 and Metzger's
2024 estimate of 11 to 26 t. The crater reaches 1.6 cm deep with a 4.3 m edge (read with the caveat above:
the vehicle is still translating through most of that);
the peak surface pressure reaches 3.3 kPa and the peak wall shear 28 Pa; the
ground-effect augmentation peaks at 24 N, well under a percent of the engine's
thrust. Viscous erosion is the dominant regime throughout, which is what the
literature reports for the Apollo landings.

The onset moved up from the fitted law's 31 m, and that is the honest cost of
deriving the threshold instead of fitting it. The 31.5 m the Apollo 16 film
gives is the height at which dust becomes *visible*; the physical onset must be
above it, and the same body of work proves it — Lane and Metzger still measure
10.7 kg/s at 36.6 m, the top of their usable range.

### Validation against published measurements

`benchmarks/studies/psi_validation/` runs the model against every published
plume-surface measurement that could be sourced precisely — Apollo descent-film
erosion rates and eroding radii, the Apollo post-landing scour depth, the
Surveyor III ejecta speeds, NASA's subscale vacuum-chamber crater tests, the
Apollo film ejection angles — with one frozen manifest per case naming its
source and its tolerance. Run it on either field:

```bash
julia --project=. benchmarks/studies/psi_validation/run_validation.jl --field=analytic
julia --project=. benchmarks/studies/psi_validation/run_validation.jl --field=table
```

In short: six of seven gate-eligible cases pass on the analytic field and five
of seven on the table, against four of six for the single fitted law. The
eroding radius and the scour depth are the two that changed, and the erosion
rate's variation with height is what still fails. Read that README before
trusting any of these quantities quantitatively; it also records the three
published laws for how the surface stress varies with height, which disagree
with each other by a factor of 30 over an Apollo descent and which this model
agrees with none of.

## Erosion and cratering regimes

`src/dynamics/coupled/force_torque_models/regolith_erosion.jl` is the module the
effector's default erosion law comes from: the separate regimes the
plume-surface literature names, each with its own onset criterion and each
sourced, replacing the single law with two fitted constants the model used to
carry. It is also usable on its own — every entry point is a pure function of a
soil description, a surface gas state, the local gravity and an
`ErosionEnvironment`, and nothing in it allocates, which is what lets the
effector call it at 128 radii inside a right-hand-side evaluation.

### The soil

`RegolithProperties`, built by `lunar_mare_regolith()`, carries the Lunar
Sourcebook's recommended lunar values (Carrier, Olhoeft and Mendell, "Physical
Properties of the Lunar Surface", chapter 9 of Heiken, Vaniman and French, *Lunar
Sourcebook*, Cambridge University Press, 1991):

| Property | Value | Source |
|---|---|---|
| bulk density | 1500 kg/m³ | Table 9.4, top 15 cm (1.50 ± 0.05 g/cm³) |
| particle density | 3100 kg/m³ | section 9.1.3, recommended specific gravity 3.1 |
| median grain size `D50` | 70 µm | section 9.1.2, average of a 40-130 µm range |
| porosity | 0.52 | `1 - 1500/3100`; Table 9.5 gives 49 percent over the top 30 cm |
| cohesion | 520 Pa | Table 9.12, 0-15 cm (range 440-620 Pa) |
| friction angle | 42° | Table 9.12, 0-15 cm (range 41-43°) |
| permeability | 3 × 10⁻¹² m² | section 9.1.8, 1-7 × 10⁻¹² m² from the Surveyor 5 vernier firing (Choate et al., 1968) |

The erosion-law coefficients it also carries come from Metzger's two 2024 *Icarus*
papers, "Erosion rate of lunar soil under a landing rocket", parts 1 and 2: the
energy threshold `E_th = 0.123 J/(m² s)`, the erosion efficiency `ε = 0.0029`,
the cohesive energy density `α = 0.289 J/m³` and the mean lift height
`⟨D⟩ = 1.5 D84`.

### The threshold that replaced the fitted 0.15 Pa

`threshold_shear_pa = 0.15` in `PlumeSurfaceConfig` was chosen to put the erosion
onset at 31 m, and is now used only by `erosion_model = :roberts_fitted`. Two
derivations sit beside it, and the first is what the default law uses:

- `shields_threshold_shear_pa(soil, g)` applies Shao and Lu's threshold
  expression ("A simple expression for wind erosion threshold friction
  velocity", *JGR* 105(D17), 2000, equation 22), which for a wall shear stress
  reads `τ_t = A_N (ρ_p g d + γ/d)` with `A_N = 0.0123` and `γ` fitted to
  1.65e-4 to 5e-4 kg/s². With the lunar soil above it gives **0.057 Pa**
  (0.033-0.092 Pa over the `γ` range), a factor of 2.6 below the fitted value.
  On the Moon the cohesive term beats the weight term by twelve to one, because
  a sixth of Earth's gravity suppresses the latter and not the former.
- `energy_flux_threshold_shear_pa(gas, soil)` recasts Metzger's energy
  threshold, itself read off the 31.5 m altitude at which dust first blows in
  the Apollo 16 landing video, as a wall shear stress:
  `τ_t = sqrt(E_th ρ_s v̄ / 3)`. At the surface state an 11.5 kN Apollo plume
  lays down at 10 m it gives **0.159 Pa**, six percent from the fitted constant.

That six percent is not independent confirmation, and the module does not claim
it is. Both numbers are anchored to the same observable -- the onset of visible
dust at roughly 31 m under an LM-class engine. What the agreement shows is that
two different plume models put comparable wall shear on the ground there.

The independent target is Stubbs and Mehta's data-derived threshold ("A
Data-Derived Scaling Approach for Plume-Surface Interaction Crater Formation",
AIAA SciTech 2026, NTRS 20250011216, section IV, Fig. 14): **0.25 Pa**. The
Shields estimate is 4.4 times below it; matching it would need a cohesion
parameter 2.8 times the top of Shao and Lu's terrestrial range, which is
plausible for airless lunar fines but unmeasured, so the default stays at their
value and the gap is reported.

The threshold and the shear law cannot be separated on today's plume field.
Applied to the Gaussian shear law, the three candidate thresholds put the
erosion onset at 24.0 m (0.25 Pa), 31.0 m (0.15 Pa) and 50.3 m (0.057 Pa) for
the Apollo approach thrust -- and Lane and Metzger (*Acta Geophysica* 63(2),
2015, Table 2) still measure 10.7 kg/s of Apollo 12 erosion at 31.9 m. The
current effector returns exactly zero there; integrating this module's local
rate over the footprint with the derived threshold gives 28 kg/s, non-zero
where the measurement is non-zero and 2.7 times high rather than infinitely
low. Whether the remaining gap belongs to the threshold or to the wall shear
the Gaussian law puts on the ground at height is a question for the plume
field, not for the soil.

### The three regimes

| Regime | Onset criterion | Source |
|---|---|---|
| `ViscousErosionEnergyFlux` | downward energy flux across the lift height above `E_th` | Metzger 2024a equation 16 |
| `ViscousErosionRoberts` | wall shear above `shields_threshold_shear_pa` | Roberts, IAS Paper 63-50, 1963; kept for comparison only |
| `DiffusionDrivenFlow` | pore-pressure uplift above the overburden weight plus tensile strength | Scott and Ko 1968, as stated in Lunar Sourcebook section 9.1.8 |
| `BearingCapacityFailure` | stagnation pressure above the ultimate bearing capacity | Lunar Sourcebook section 9.1.9, after Durgunoglu and Mitchell 1975 |

`erosion_rate(regime, gas, soil, g, env)` returns kg/(m² s) and is identically
zero below the regime's onset; `erosion_onset(...)` is the predicate alone.
`regolith_erosion_rate(default_erosion_regimes(), gas, soil, g, env)` evaluates
all three and reports the total and which one dominates:

```julia
soil = lunar_mare_regolith()
gas  = (pressure_pa=168.0, shear_pa=1.44, density_kg_m3=8.7e-4,
        speed_mps=1000.0, temperature_k=500.0, mach=3.0)
out  = regolith_erosion_rate(gas, soil, 1.625, erosion_environment(footprint_radius_m=4.66))
out.dominant === ViscousErosion
```

With the lunar soil and the Apollo LM the model puts viscous erosion alone in
play throughout the descent. Diffusion-driven flow needs the plume on the
ground at better than about 20 kN, which the LM is below by touchdown; bearing
capacity failure needs about 450 kN through the same nozzle, ten times the
descent engine at full throttle. That ordering is what Metzger reports for the
Apollo landings, and it was not tuned in: it follows from the Sourcebook's
permeability and Table 9.12's strength.

### What is assumed rather than sourced

Three things, all documented configuration fields rather than buried literals:
`saltation_efficiency` (the fitted factor of 10 in the Roberts closure — which
the effector's default law *does* use, and which is therefore the model's one
remaining fitted constant; the validation study records the value that best fits
Lane and Metzger's Apollo 12 profile, 17.5 on the analytic field and 13.7 on the
table, and why it was left at 10); the exhaust viscosity, for which
`gas_dynamic_viscosity_pa_s` uses Sutherland's law for air as a stand-in; and
the `ErosionEnvironment` geometry and timing -- the footprint radius, the gas
residence time that sets the pressure diffusion depth, and the width used for
the bearing-capacity factors. The effector fills the first and third of those
from the plume's own footprint at every evaluation and carries the residence
time as `gas_residence_time_s`.

One known gap is stated rather than tuned away. `soil_bearing_capacity_pa` uses
the classical Prandtl/Reissner/Vesic factors and returns about 1.3 MPa for a 1 m
footing on the Sourcebook's 0-60 cm soil, where the Sourcebook itself quotes
roughly 6 MPa (3-11 MPa for the Apollo 11 footpad) from Durgunoglu and
Mitchell's larger wedge-penetration factors. The function is therefore a lower
bound on the soil's strength, which makes the bearing-capacity onset
conservative: the model fires that regime earlier than the Sourcebook's own
numbers would. It changes nothing for an Apollo-class lander, which reaches
about 25 kPa of stagnation pressure at contact against 1.3 MPa even on the
conservative factors, but a lander an order of magnitude larger would cross the
model's threshold before it crossed the Sourcebook's, and this is the reason.

## Ejecta transport

`ejecta_transport.jl` replaces the single characteristic ejecta speed above with
grain trajectories. The effector calls it once per saved sample (see "Ejecta as
a per-sample diagnostic" above), and it is also a standalone, pure module, so it
can be used to size a landing without running one.

```julia
field = EjectaReferenceGasField()            # or the plume field module's table
cfg = EjectaTransportConfig()
soil = EjectaSoil()                           # lunar mare, Lunar Sourcebook ch. 9

# one grain: how fast the wall jet gets it going over 3.3 m of surface
gas = ejecta_gas_state(field, cfg, 45_040.0, 10.0, 3.3)
launch = ejecta_launch_speed(gas, soil, 70e-6, 3.3; config=cfg)

# the whole population under a 45 kN engine 10 m above the ground
dist = ejecta_distribution(field, cfg, soil, 45_040.0, 10.0)
dist.mean_speed_mps, dist.mean_angle_deg, dist.escape_fraction
```

### The flow the grain is in

A 70 µm grain under a descent engine is **not in continuum flow**. At the radius
where the wall shear stress peaks, with the engine at full thrust 10 m above the
ground, the reference field gives a gas density of 1.2e-3 kg/m³ and a wall-jet
speed of 2.2 km/s, so the particle Reynolds number is about 8, the particle Mach
number about 3.4 and the particle Knudsen number about 0.6 -- the transitional
regime on the Schaaf and Chambré classification (continuum below Kn = 0.01, slip
to 0.1, transitional to 10, free molecular above), as tabulated in Capecelatro's
review of spacecraft-landing gas-particle flows. Higher on the approach, and at
throttled thrust, the same grain is in free molecular flow.

The drag law is therefore Henderson's correlation (C. B. Henderson, "Drag
coefficients of spheres in continuum and rarefied flows", *AIAA Journal* 14(6),
1976, 707-708), which spans continuum through free-molecular flow and is valid
over Re < 2e4 and Ma < 6. `ejecta_drag_coefficient` implements it; the tests
check that as Re → 0 it converges on the closed-form free-molecular drag of a
sphere with diffuse reflection (Schaaf and Chambré, *Flow of Rarefied Gases*,
1961) -- to 5 percent at Ma = 2 and better than 1 percent at Ma = 6.

### Launch, flight and the distribution

`ejecta_launch_speed` integrates the drag balance
`dv/dt = 3 ρ C_D(|u-v|) (u-v)|u-v| / (4 ρ_p d)` from rest over an entrainment
length. The gas speed is the asymptote, which is why the finest grains approach
the exhaust velocity itself -- about 3,100 m/s for the Apollo lunar module,
the figure Metzger uses when estimating ejecta damage to lunar orbiters. Larger
grains fall behind as 1/d.

`ejecta_trajectory` flies the grain ballistically under the body's gravity
through gas that decays both with radius (the field's own profile) and with
height above the surface (an exponential over a wall-jet thickness that grows
linearly with radius). With no gas it reduces exactly to `v² sin 2θ / g`.

`ejecta_distribution` sweeps launch radii and grain sizes and returns speed,
angle and deposition histograms plus the summary scalars a results table would
carry: `mean_speed_mps`, `median_speed_mps`, `max_speed_mps`, `mean_angle_deg`,
`mean_deposition_radius_m`, `p90_deposition_radius_m`, `escape_fraction` and
`escape_speed_mps`.

The ejection angle is an **input**, not a computed quantity: Immer, Lane,
Metzger and Clements ("Apollo video photogrammetry estimation of plume
impingement effects", *Icarus* 214, 2011) measured 1-3 degrees above the local
horizontal, measured at launch, from the Apollo landing films (their Table 1:
2.6°, 2.4°, 8.1°, 1.4° and 2.0° for Apollo 11, 14, 15, 16 and 17, mean 3.3°, the
outlier attributed to an 11° surface slope; Apollo 12 could not be measured),
and they note that on Roberts' theory the angle follows the scour crater's wall
slope, so it is only weakly coupled to thrust.

The model now has a crater, so that slope is computable, and it does not yet
supply the angle. Measured off the crater profile the model digs over Lane and
Metzger's Apollo 12 descent, the inner wall stands at 0.80° and the outer at
0.44°, with a steepest local outer slope of 1.38° (1.01°, 0.47° and 1.60° with
the tabulated field). Against a measured 1.4 to 8.1° that is the right order of
magnitude and at or below the bottom of the range — but the crater's deepest
point sits at 0.18 m where the reference puts it at 1 to 2 m, so its radial
scale has to be right before its slope can be trusted to set anything. The angle
stays an input, and the validation study records the residual.

`ejecta_distribution` publishes the angle in exactly that convention -- degrees
above the local horizontal at launch -- as the mass-weighted `mean_angle_deg`
and as the `angle_edges_deg`/`angle_fraction` histogram. Because the range is an
input, the mean comes out at 1.99°, inside the 1.4-8.1° spread the five measured
Apollo landings cover and just below their 3.3° mean; the model cannot currently
be said to *predict* the angle, only to carry it consistently.

### What the model gives on the Apollo descent envelope

Evaluated against the reference gas field, with equal mass weight per size bin
and the default 1 to 500 µm sizes:

| Engine | Height | Mean speed | Max speed | Mean angle | Escape fraction |
|---|---|---|---|---|---|
| 4.5 kN (minimum throttle) | 40 m | 149 m/s | 653 m/s | 1.99° | 0 |
| 11.5 kN (approach) | 30 m | 253 m/s | 1055 m/s | 1.99° | 0.004 |
| 45 kN (full) | 10 m | 610 m/s | 2145 m/s | 1.99° | 0.21 |
| 45 kN (full) | 2 m | 895 m/s | 2710 m/s | 1.98° | 0.53 |

These sit inside the only measured bounds there are: Immer et al. summarize the
Surveyor 3 pitting analyses as 40 m/s (Nickle and Carroll 1972), > 70 m/s
(Jaffe 1972), 100 m/s (Cour-Palais et al. 1972) and, from the pit structure and
called the most reliable, 300 to 2000 m/s (Brownlee, Bucher et al. 1972); lunar
escape is 2373 m/s. Those bounds are on the fast tail that reached Surveyor 3,
not on the mean of the population, so agreement here is a sanity check, not a
validation.

Those rows are the module on its own: the reference gas field, equal mass per
size bin, no erosion threshold. The effector runs it differently — the plume
field's gas state, the regimes' local erosion rate as the radial weight and the
lognormal grain-mass distribution as the size weight — and that is what the
saved columns carry. The difference, at the Apollo approach thrust with the
analytic field 10 m above the ground:

| Weighting | Mean speed | Mean deposition radius |
|---|---|---|
| wall shear, equal mass per size bin | 143 m/s | 1553 m |
| erosion rate, lognormal grain mass | 107 m/s | 755 m |

The deposition radii are still hundreds of meters to kilometers, far outside the
tens of meters of visible blast zone around an Apollo site, and the weighting is
no longer the reason. A 70 µm grain leaving at 107 m/s and 2° above the
horizontal has a ballistic range of 424 m in lunar gravity; that follows from
the measured ejection angle and the measured speeds and nothing else. The
visible scour is where the *coarse mass* lands, not where the population lands,
and the model has no separate output for it — so the deposition radius should be
read as what it is, the mass-weighted mean range of the whole population, and
not compared with a blast-zone radius.

### What is assumed

Every free parameter is a named configuration field, with its reasoning on the
field in `ejecta_transport.jl`:

| Field | What it absorbs |
|---|---|
| `entrainment_length_factor`, `entrainment_length_min_m` | the distance the wall jet accelerates a grain over before it leaves the surface; the most influential free parameter of the launch model |
| `wall_jet_growth_rate`, `wall_jet_thickness_min_m` | how quickly the drag falls off as the grain climbs out of the jet (linear growth is the measured behavior of a radial wall jet, the rate of 0.1 and the exponential profile are not) |
| `molecular_collision_diameter_m` | the hard-sphere diameter behind the Chapman-Enskog viscosity, hence the Reynolds number |
| `grain_temperature_k` | enters only as T_p/T_gas in the drag law |
| `gas_molar_mass_kg_mol`, `gas_gamma` | exhaust composition; a plume field carrying its own should override them |
| `EjectaReferenceGasField.exit_static_temperature_k`, `wall_jet_decay_radii` | the reference gas field is a stand-in for testing, not a plume model |

| the lognormal *shape* of `ejecta_lognormal_mass_weights` | the two percentiles it is fitted to (`D50` and `D84/D50`) are sourced; that the distribution between them is lognormal is not |

The mass weighting in `ejecta_distribution` still *defaults* to a proxy
proportional to the local wall shear stress with no threshold, for callers using
the module on its own. It is not an erosion model; pass the regolith erosion
rate as the `weight` keyword to replace it, which is what the effector does.

## Surface in the viewer

`export_visualization(prefix; terrain="data/terrain/moon/apollo11/site.json")`
embeds the site's grids and imagery. The page drapes each imagery level over a
patch displaced by the DEM, nested from coarse to fine and drawn all the time,
so the surface sharpens by itself as the camera closes in on the site; the
globe is cut open under the outermost patch. The selection panel gains a
"height above terrain" quantity, and the info panel names the terrain and its
finest resolution.
