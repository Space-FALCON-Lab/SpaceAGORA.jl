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
down, the regolith it erodes, the speed the grains leave at, and the small
thrust augmentation the reflected plume gives the vehicle in the last couple
of nozzle diameters. It is constructed with the descent control effector
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
propulsion system by default), the soil (1500 kg/m³ bulk density, 70 µm grains,
cohesion, threshold shear stress) and the plume field the surface state is read
from.

### The model

Every surface quantity is read from a *plume field* that `PlumeSurfaceConfig`
carries in its `field`, so the gas dynamics and the erosion closure are
separable. Two fields exist. `PlumeAnalyticField` is the Gaussian footprint
described in this section and is the default, so a run that asks for nothing
else is unchanged. `PlumeFieldTable` is a plume computed from the nozzle and
tabulated; it is described in "The tabulated plume field" below.

The erosion closure on top of either field follows Roberts' treatment of a
hypersonic jet acting on a dust layer (L. Roberts, *The action of a hypersonic jet on a dust layer*, IAS
Paper 63-50, 1963) in the form used for the Moon by Metzger and co-workers
(Metzger, Immer, Donahue, Vu, Latta and Deyo-Svendsen, *Jet-induced cratering
of a granular surface with application to lunar spaceports*, J. Aerosp. Eng.
22, 2009; Metzger, Smith and Lane, *Phenomenology of soil erosion due to rocket
exhaust on the Moon and the Mauna Kea lunar test site*, J. Geophys. Res. 116,
E06005, 2011).

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

Soil moves where `τ(r)` exceeds the threshold shear stress `τ_t`, which happens
in an annulus, not a disk. Roberts' viscous erosion is closed as a momentum
balance over that annulus: the excess shear is the force available to
accelerate grains, and a mass flux `ṁ` leaving at the ejecta speed `v_ej`
carries `ṁ v_ej` of momentum, so

```math
\dot m = \frac{\eta}{v_{ej}} \int \max(\tau(r) - \tau_t,\, 0)\, \mathrm{d}A
```

The ejecta speed comes from a drag balance on a single grain accelerated across
one footprint radius by the gas dynamic pressure at the shear peak; it runs
from tens of m/s at the erosion onset to the configured cap near the ground,
the range Metzger and Immer measured from the Apollo landing films. The
cumulative eroded mass is the time integral of `ṁ`, advanced only when the
solver moves past every time already integrated, so rejected steps neither
double-count nor run it backwards.

Two parameters are calibrated rather than derived, and both are anchored on
Apollo 11 observables:

- `threshold_shear_pa` (0.15 Pa) sets the erosion onset height. With the
  approach thrust near the end of the descent — about 11.5 kN of lunar weight —
  `plume_erosion_onset_height` returns 31 m, and the Apollo 11 crew first
  reported blowing dust at about 30 m. The value is three to four orders of
  magnitude below the bulk cohesion of lunar regolith (0.1 to 1 kPa), as the
  mobile surface layer must be. Above the onset every plume quantity is zero.
- `erosion_efficiency` (10) multiplies the momentum-balance rate to cover the
  saltation cascade, in which each impacting grain splashes several more — a
  mechanism a pure momentum balance cannot produce.

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
| erosion onset height | 31 m | 80 m | 2.6 |

The pressure and the footprint agree to about 20 percent over the whole descent,
which is the real verdict on the Gaussian footprint: as a *pressure* closure it
was close. The shear stress does not agree, because Roberts' drag coefficient
acting on the wall jet is a far stronger coupling than a 0.01 skin-friction
coefficient acting on the static pressure, and the erosion onset follows it.
`threshold_shear_pa` and `erosion_efficiency` were calibrated against Apollo 11
observables with the analytic field and are *not* retuned for the table, so a
run with a table erodes more soil, earlier: re-evaluating the Apollo 11 descent
on the tabulated field moves 29 tonnes of regolith against the analytic field's
429 kg, which is far above the tonne-scale estimates from the landings. The
calibration does not transfer, and it is not adjusted here to hide that; the
numbers quoted below are the analytic ones. Recalibrating the erosion closure
against the literature is separate work.

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
so the shear reaches far beyond the pressure footprint. A shear proportional to
the static pressure, which is what the analytic field uses, forces the eroding
region and the pressure footprint to be the same width, and that single
geometric error is what the wide measured region exposes.

Use `plume_wall_shear`, `plume_mean_shear` and `plume_scour_radius` to compare
against measurements of a profile rather than of a peak; comparing a model peak
against an area-averaged measurement flatters the analytic field considerably.
`data/psi/README.md` gives the full comparison, the two published references
that disagree with each other by a factor of 25, and the erosion-onset residual
that is left after both are taken into account.

### What it records

The effector keeps a `PlumeSurfaceState` per spacecraft, updated at every
right-hand-side evaluation, and `default_save_fields` publishes it as seven
result columns whenever the effector is in the run:

| Column | Meaning |
|---|---|
| `sc{i}_plume_height_m` | height above the terrain along the engine axis, from the vehicle's reference point |
| `sc{i}_plume_shear_pa` | peak wall shear stress under the plume |
| `sc{i}_plume_pressure_pa` | peak surface pressure under the plume |
| `sc{i}_plume_erosion_kg_s` | mass erosion rate |
| `sc{i}_plume_eroded_kg` | its time integral |
| `sc{i}_plume_ejecta_mps` | characteristic ejecta speed |
| `sc{i}_plume_ground_effect_n` | ground-effect thrust augmentation |

Run with `isolate_state=false` to read `plume.state` directly afterwards.
`plume_quantities`, `plume_surface_footprint`, `plume_erosion_onset_height` and
`plume_ground_effect_force` are the same model as plain functions, for sizing a
scenario without running one.

On the Apollo 11 descent
(`scripts/dev/viewer_demos/apollo11_landing.jl`) erosion begins 33 s before
touchdown with the engine 35 m above the ground, peaks at 18 kg/s, and moves
about 430 kg of regolith in all — the same order as the tonne-scale estimates
Metzger and co-workers derive from the Apollo landings, and low by roughly a
factor of two, which is where a model with one lumped cascade multiplier should
be expected to sit. The peak surface pressure reaches 3.3 kPa and the peak wall
shear 29 Pa; the ground-effect augmentation peaks at 25 N, well under a
percent of the engine's thrust.

## Surface in the viewer

`export_visualization(prefix; terrain="data/terrain/moon/apollo11/site.json")`
embeds the site's grids and imagery. The page drapes each imagery level over a
patch displaced by the DEM, nested from coarse to fine and drawn all the time,
so the surface sharpens by itself as the camera closes in on the site; the
globe is cut open under the outermost patch. The selection panel gains a
"height above terrain" quantity, and the info panel names the terrain and its
finest resolution.
