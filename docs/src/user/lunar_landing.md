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
  site, by HTTP range requests straight out of the PDS raster, and a second,
  box-averaged LOLA window over the whole imagery region so the ground has
  relief out to the horizon;
- the LROC NAC digital terrain model of the site (2 m/px, resampled to 4 m)
  where one exists (Apollo 11 does);
- an imagery quadtree that follows the descent corridor, from the LROC WAC
  global mosaic through Kaguya's TC ortho mosaic to the Apollo 11 NAC mosaics
  on NASA Moon Trek, and below Moon Trek's own cap from the Planetary Data
  System (0.33 m/px at the site).

Those mosaics are not served on the same control, and since a level of the
quadtree takes its detail from whichever source reaches its zoom, an
uncorrected offset would put the same crater in two places either side of a
level change. `LAYER_REGISTRATION` in the script holds each source's measured
offset in meters and the sampling window is shifted to match: the two Apollo
11 NAC mosaics agree with each other to within a pixel but sit 24 m south and
4 m east of the site's PDS NAC digital terrain model, `A11_60x60km` sits 36 m
north of them, and the PDS browse raster that feeds the deepest level is 0.6 m
south and 1.9 m west of the corrected NAC mosaics. The comment in the script
records how each number was measured, and
`scripts/dev/terrain/check_site_registration.py` re-measures them:

```
python3 scripts/dev/terrain/check_site_registration.py --site 0.67416 23.47314 \
    --dem data/terrain/moon/apollo11/dem_nac.json \
    --cache data/terrain/moon/apollo11/pds_cache \
    --trek-cache data/terrain/moon/apollo11/trek_cache --sun-azimuth 90
```

It reports three independent readings of where a source draws the ground: a
phase correlation against a hillshade of the site's NAC digital terrain model
(the approach illumination is a sun azimuth of 90°), a correlation against the
Moon Trek layer feeding the level above, and the position of the lunar module
itself against the Wagner et al. (2017) coordinate the simulation targets.

```
python3 scripts/dev/terrain/fetch_moon_site.py --site 0.67416 23.47314 --out data/terrain/moon/apollo11
```

The data stay under `data/` (not tracked); `SPACEAGORA_TERRAIN_SITE` points
the demo at a `site.json` somewhere else, such as a regenerated site under
`output/terrain`. Grids are little-endian Float32 files with a JSON header, so
other DEMs can be dropped in the same layout.

### The imagery quadtree

The root of the quadtree is a square region one Moon Trek tile wide at
`--root-zoom` (22.5°, 682 km, by default), its corner snapped to that zoom's
pixel grid and centered on the ground track, so a node at level `L` is exactly
one Trek tile at zoom `root_zoom + L` and every node is an integer-pixel crop
of Trek's own tiles rather than a resample. `--approach-azimuth` and
`--uprange-km` (270° and 480 km, the Apollo 11 descent) place the track, and
`--profile` gives the height above the site against the distance still to go.

Nodes are built only where they are worth having, so the coverage is a funnel,
not a full pyramid: a level's half-extent is `--lod-factor` times its own node
side — the radius within which a view-dependent quadtree still asks for that
level — intersected with the part of the track from which the camera can see
that far. The coarse levels therefore span the whole 480 km corridor and each
finer level narrows toward the site. A node that was not built is not an error:
the viewer textures it from the nearest present ancestor and that ancestor's
matching UV sub-rectangle, so every point of the root is textured at the best
resolution that exists for it and the sharpness falls off gradually with
distance from the site instead of ending at a seam. `--quality-coarse` and
`--quality-fine` ramp the JPEG quality from the wide coarse tiles to the ones
at the site, which is where the page's byte budget is worth spending.

Re-tiling is cheap: Moon Trek tiles are cached under `<out>/trek_cache` and
archive tiles under `<out>/pds_cache`, and `--reuse DIR` borrows another site
directory's DEM files, raw downloads and both caches, so `--retile` with
different extents costs no new downloads.

### Below Moon Trek: the archive level

Moon Trek declares a deepest tile matrix per layer, and for the Apollo 11 NAC
mosaics it is zoom 15, 0.651 m/px — about 2.5 times coarser than the mosaic
Trek is serving. Level 13 of the quadtree (0.325 m/px) therefore comes from
the Planetary Data System instead, through `ARCHIVE_SOURCES` in the fetch
script: the map-projected browse raster
`NAC_PHO_E010N0230_M175124932R.PYR.TIF` of LROC NAC frame **M175124932R**,
imaged 2011-11-05 from a 24 km low-altitude pass, which is the finest frame
over this site and the frame Trek's own mosaic was made from. The 915 MB
raster is read in place: it is uncompressed and internally tiled 256 × 256, so
the reader takes the tile table out of its TIFF directory and fetches a whole
tile row per HTTP range request, caching every 64 kB tile so a re-run or a
resumed run costs nothing. Its georeferencing is read from the file
(`ModelTiepointTag`, `ModelPixelScaleTag` and the GeoTIFF keys: equirectangular
about 180° with a standard parallel of 1° on a 1737.4 km sphere) and checked
against the PDS4 label's bounding coordinates, which agree to 0.17 m. The
source grid is 0.2302 m/px, not a power-of-two relative of Trek's, so the
archive level is the one place the pipeline area-averages rather than cropping.

The archive rasters in the PDS are public domain; images served through the
LROC project's own non-PDS interfaces are not, so only the PDS product ships.
`imagery/tiles.json` carries the product, its URL and the credit, and the
attribution is **NASA/GSFC/Arizona State University**, with the instrument
reference Robinson, M. S., et al. (2010), "Lunar Reconnaissance Orbiter Camera
(LROC) Instrument Overview", *Space Science Reviews* **150**, 81–124.

**What that resolution is worth.** 0.325 m/px is the grid the tiles are drawn
on, not the size of the smallest thing in them. The frame samples the ground
0.24 m across the detector line but 0.56 m along track (`SCALED_PIXEL_WIDTH`
and `SCALED_PIXEL_HEIGHT` of M175124932R in the PDS volume index), and power
spectra of the archive raster lose their signal at about 0.5 m of feature
size, so the imagery resolves about half a meter and the popular "26 cm"
figure is the cross-track sampling alone. `tiles.json` records this in its
`resolution` block and the viewer's info panel reports it beside the sampling.
Level 14 is not built: it would cost another 0.25 MB of the page's budget to
interpolate detail the data does not contain.
`docs/reference/lunar_imagery_sources.md` has the survey and the measurements
behind all of this.

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

## Surface in the viewer

`export_visualization(prefix; terrain="data/terrain/moon/apollo11/site.json")`
embeds the site's grids and its imagery quadtree. The page walks the quadtree
against the camera and drapes each node it needs over geometry displaced by the
DEM, so the surface sharpens by itself as the camera closes in on the site and
stays covered out to the horizon while the lander is still uprange. The
selection panel gains a "height above terrain" quantity, and the info panel
names the terrain and its finest resolution.
