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
propulsion system by default) and the soil (1500 kg/m³ bulk density, 70 µm
grains, cohesion, threshold shear stress).

### The model

The gas-dynamic part follows Roberts' treatment of a hypersonic jet acting on a
dust layer (L. Roberts, *The action of a hypersonic jet on a dust layer*, IAS
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

## Ejecta transport

`ejecta_transport.jl` replaces the single characteristic ejecta speed above with
grain trajectories. It is a standalone, pure module -- nothing in the descent
run calls it yet -- so it can be used to size a landing or to generate ejecta
result columns without touching the effector.

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

The deposition radii the model gives are kilometers, far outside the tens of
meters of visible blast zone around an Apollo site. That is not obviously wrong
-- the visible scour is where the *coarse mass* lands, while the kilometer-scale
figures are driven by the micron fines that carry almost no mass but travel at
hundreds of meters per second -- but until the weighting comes from a real
erosion model and a real size distribution the deposition numbers should be read
as an upper envelope, not a prediction.

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

The mass weighting in `ejecta_distribution` defaults to a proxy proportional to
the local wall shear stress with no threshold. It is not an erosion model; pass
the regolith erosion rate as the `weight` keyword to replace it.

## Surface in the viewer

`export_visualization(prefix; terrain="data/terrain/moon/apollo11/site.json")`
embeds the site's grids and imagery. The page drapes each imagery level over a
patch displaced by the DEM, nested from coarse to fine and drawn all the time,
so the surface sharpens by itself as the camera closes in on the site; the
globe is cut open under the outermost patch. The selection panel gains a
"height above terrain" quantity, and the info panel names the terrain and its
finest resolution.
