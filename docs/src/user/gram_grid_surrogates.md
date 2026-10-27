# GRAM Grid Surrogates

A grid surrogate is a GRAM-derived atmosphere frozen at one instant and stored
as data, so a simulation can use it without installing native GRAM. The public
GRAMSuite Julia package evaluates it. That package is the lab's open wrapper,
not NASA's GRAM Suite; the Odyssey example environment
(`examples/odyssey_surrogate_env`) installs it.

Related pages:

- [Atmosphere Models](atmosphere_models.md)
- [Odyssey guidance with a frozen atmosphere](../tutorials/odyssey_surrogate.md)
- [Assets & Modes](../assets.md)
- [GRAMSuite Setup](gramsuite_setup.md)

## Fixed GRAM grid snapshot

For a named atmosphere that SpaceAGORA retrieves for you, start with the
[Odyssey surrogate workflow](../tutorials/odyssey_surrogate.md). It loads
`GRAMSuite` and then calls
`surrogate_preset_model("odyssey_p20_frozen_v1"; version="1.0.0")`. Each named
preset is a lazy Julia artifact declared in the repository's `Artifacts.toml`
and published as a release asset of GRAMSuite.jl. The first call downloads it
into the Julia artifact store. Every call then checks the payload's size and
SHA256 against the catalog, `data/surrogate_presets.toml`, before it builds the
model, so retrieval and verification are finished before the solver starts.
With `offline=true`, a preset that is not installed yet is refused instead of
downloaded. `available_surrogate_presets()` lists the catalog without
downloading anything.

Four named presets are published, six versions in all. Every call names an
exact version; there is no "latest". All of them are frozen at the Odyssey P20
instant, 2001-11-07T11:51:04.794789Z:

| Preset | Domain | Validated use |
| --- | --- | --- |
| `odyssey_p20_frozen_v1` 1.0.0 | 100 to 260 km, 40 to 90 degrees north | Odyssey P20 passages within the tutorial's envelope |
| `mars_global_upper_p20_frozen_v1` 1.0.0 | 80 to 365 km, all latitudes and longitudes | Pointwise within 225 s of the frozen instant, down to 80 km; propagated passes designed for periapsis from 85 to 130 km (see below) |
| `mars_global_near_surface_p20_frozen_v1` 1.0.0 | 5 m above the local surface to 75 km areoid height, planetocentric latitudes within 85 degrees, surface below 9 km; density, temperature and pressure, no winds | Pointwise at the frozen instant |
| `mars_global_near_surface_p20_frozen_v1` 1.1.0 | As 1.0.0, up to 81 km areoid height, above the global upper preset's floor everywhere | Pointwise at the frozen instant |
| `mars_global_near_surface_p20_frozen_v1` 1.2.0 | As 1.1.0, with stored east, north and vertical winds; refused where a stored wind it needs is unavailable | Pointwise, on fixed descent paths and in propagated descents at the frozen instant; ballistic entries in the engine, with the global upper preset above 80 km |
| `mars_global_descent_p20_frozen_v1` 1.0.0 | A composite with no data of its own: near-surface 1.2.0 below 80 km above the ellipsoid, the global upper preset 1.0.0 at and above it, up to 365 km | Ballistic descents from 125 km to 100 m above the terrain, within the tested envelope (see [Mars descent preset](#Mars-descent-preset)) |

The global upper preset was validated against native Mars-GRAM under the lab's
release limits for pointwise density and wind error and for per-pass drag and
heating; its archive README lists the limits. The three kinds of check reach
different heights. Pointwise checks cover the whole domain, down to its 80 km
floor. Propagated passes were designed for periapsis from 85 to 130 km, and
the lowest altitude they reached was 84.56 km. Fixed-path arcs, which follow
prescribed paths without trajectory feedback, reach 80.25 km and check drag
and heating proxies along the path.

No complete propagated pass was validated below 84.56 km. Between 80 and
84.56 km the evidence is the pointwise density and wind checks and the
fixed-path drag and heating proxies, and these do not establish
propagated-pass accuracy in that interval. The archive README gives the pass
range as periapsis from 80 to 130 km; the figures above state the coverage
precisely.

The global upper grid's height and latitude spacing is not uniform: there are
extra nodes where Mars-GRAM has structure. The catalog lists every node, and
`surrogate_preset_model` checks the loaded grid against that list. The grid
stops at 80 km and queries below it fail. The near-surface preset below covers
heights up to 75 km (1.0.0) or 81 km (1.1.0 and 1.2.0).

```julia
using SpaceAGORA
import GRAMSuite
density_model = surrogate_preset_model("mars_global_upper_p20_frozen_v1"; version="1.0.0")
```

### Near-surface Mars preset

`mars_global_near_surface_p20_frozen_v1` covers the lower atmosphere down to
the ground. It is not a grid. It follows native Mars-GRAM's own near-surface
rule, which is driven by the local terrain. The query's own MOLA surface
height gives the first exposed table level, and that level selects the
regime: interpolation between table levels above it, and a surface-layer law
from 30 m above the surface up to it and from 5 m to 30 m. The stored
components sit on a 1.5-degree lattice and are interpolated bilinearly.

For this preset `surrogate_preset_model` returns a
`GRAMNearSurfaceAtmosphereModel`, whose `getDensity` returns density,
temperature and the wind vector. Versions 1.0.0 and 1.1.0 have format 1
payloads, which store no winds, so the wind vector is zero and a simulation
using them has no atmospheric wind. Version 1.2.0 has a format 2 payload and
returns its stored east, north and vertical winds, with the horizontal
components clipped at 0.7 times the speed of sound, as native clips them; the
vertical wind is not clipped. A simulation with its wind setting off
(`EnvironmentModel.wind = false`) uses zero wind, as it does with the grid
presets' stored winds. The catalog declares which format each version has, and
the loaded payload must agree.

Two GRAMSuite functions return more than `getDensity` does. Both take
latitude and longitude in degrees and height in meters.
`GRAMSuite.near_surface_state(model.core, lat_deg, lon_deg, h_m)` returns
pressure, the regime and the status of the surface-layer model used. For a
format 2 payload,
`GRAMSuite.near_surface_wind_state(model.core, lat_deg, lon_deg, h_m)` adds
the winds, the unclipped horizontal winds, the speed of sound and which
components were clipped.

```julia
using SpaceAGORA
import GRAMSuite
density_model = surrogate_preset_model("mars_global_near_surface_p20_frozen_v1"; version="1.2.0")
rho, T, wind = getDensity(density_model, 250.0, deg2rad(-4.5), deg2rad(137.4), 0.0, true)
```

All three published versions stay available. Version 1.0.0 is unchanged.
Version 1.1.0 adds two table levels, at 80.0323 and 85.0323 km areoid height,
where Mars-GRAM places the first two levels of its upper table at this
instant, and raises the top from 75 to 81 km. Between 75 km and the first
added level, native joins its lower and upper tables, and the preset follows
the same rule. Below 75 km, 1.1.0 gives the same results as 1.0.0. Version
1.2.0 adds the stored winds and the speed of sound that limits them; its
density, temperature and pressure equal those of 1.1.0 wherever it serves a
query.

A query fails with a `DomainError` that names the reason when its
planetocentric latitude is beyond 85 degrees, the surface height is 9 km or
more (volcano flanks), it is less than 5 m above the surface, it is above the
version's top (75 km areoid height in 1.0.0, 81 km in 1.1.0 and 1.2.0), a
component it needs is unavailable at that position, or no surface-layer model
exists for its cell. In 1.2.0 a query also fails when a stored wind value it
needs is unavailable or its speed of sound is not determined.

Within 85 degrees of the equator, about 98% of the area is served down to 5 m
with a qualified surface-layer model and about 1.7% with a provisional model,
which each result's status flags. The rest is refused at the lowest heights or
on volcano flanks. From 75 to 81 km (1.1.0 and 1.2.0), every location with a
surface below 9 km is served. Over about 99.65% of the area, 1.2.0 serves
winds wherever it serves the scalars. Over about 0.07%, in deep, narrow
terrain, a stored wind is unavailable at one or more heights between −5 and
7 km areoid height and those queries are refused;
`GRAMSuite.near_surface_state` still serves their scalars, as 1.1.0 does. The
archive's support map lists the cells.

Each version was validated pointwise against native Mars-GRAM at the frozen
instant, under the lab's release limits. Version 1.0.0 was validated
separately for qualified-model areas, provisional-model areas and all served
queries. For 1.1.0, the added heights were validated with every point's
expected outcome specified in advance; the 1.0.0 validation carries over below
75 km, where results are unchanged. For 1.2.0, the stored winds and the speed
of sound were validated against native Mars-GRAM pointwise, on fixed descent
paths and in propagated descents, and in SpaceAGORA's engine: ballistic
entries flown with 1.2.0 below 80 km and the global upper preset above, against
native Mars-GRAM in the same engine loop. The scalar validations carry over.

The 1.2.0 archive README lists the limits, what was not established and three
exemptions. At the exact wind jumps at 78.75 and 82.5 degrees north and on
exact first-level contours, native can resolve a query to the other side, and
the values 1.2.0 returns there were validated only where its own side matched
native. Native's extrapolation band at the contours and the 80 km composition
switch are the other two.

Some limits apply across versions. In 1.0.0, the 75 to 80 km interval lies
between this preset and the upper one, and neither covers it. Versions 1.1.0
and 1.2.0 reach 81 km areoid height, which is above the upper preset's 80 km
ellipsoidal floor everywhere, so both presets are defined between that floor
and 81 km. `CombinedAtmosphereModel` hands over between them; see
[Height-based composition of two frozen atmospheres](#Height-based-composition-of-two-frozen-atmospheres)
and the [Mars descent preset](#Mars-descent-preset). The 1.2.0 validation does
not establish real-world landing accuracy, other vehicles, epochs or entries,
or the engine's use of the preset's speed of sound: the engine computes Mach
number from the planet's own heat-capacity ratio and gas constant.

The near-surface payloads are published under CC BY 4.0. Each payload's
terrain component holds Mars-GRAM's MOLA terrain values at its lattice nodes;
credit NASA MOLA as the archive README states.

### Load your own grid

`GRAMGridAtmosphereModel` connects GRAMSuite's offline interpolation kernel to
SpaceAGORA. It needs the GRAMSuite package with its native-free grid API and a
trusted serialized grid payload. Neither construction nor density evaluation
uses a native GRAM installation or a native fallback. The example is a
template; put the path of your own file in place of the placeholder:

```julia
using SpaceAGORA
import GRAMSuite

density_model = GRAMGridAtmosphereModel(
    planet="earth",
    surrogate_file="/path/to/authorized/earth_surrogate.jls",
    above_grid=:error,
)
```

The constructor passes grid-loader options through to GRAMSuite, including
`search_roots` and `expected_sha256`. Without `surrogate_file`, the loader
looks for `<planet>_surrogate.jls` in `search_roots` and then in a list of
known locations. When you give no `search_roots`, that list includes this
repository's `data/planets/gram_surrogate`, where the per-planet grids are Git
LFS files. The first file found is used: a missing file, a Git LFS pointer
that was never downloaded or an invalid payload is an error, not a reason to
try the next location.

A checksum pins a file's contents. It does not establish scientific
provenance or permission to distribute the file. For a custom grid you provide
the payload yourself and should record its source, access conditions and
checksum. The named presets have their own retrieval and verification route,
and that route does not make an arbitrary local grid publicly reproducible.

Queries take altitude in meters and latitude and longitude in radians, on the
grid's documented reference surfaces. Results are density in kg/m³,
temperature in kelvin, and local east/north/up winds in m/s. The stored
atmosphere is frozen: neither elapsed time nor the wind argument changes its
density or stored winds. Choose and validate the atmospheric epoch,
coordinates, forcing and coverage for your mission before you rely on a grid.
Older payload metadata may leave some of these unknown, and loading a file
cannot recover them or establish physical accuracy.

By default, altitudes and latitudes outside the grid are rejected; longitude
is periodic. With `above_grid=:vacuum`, a query above the ceiling at a valid
latitude returns zero density and wind and the temperature
`vacuum_temperature` (default 200 K). Below the floor there is still no
extrapolation. `above_grid` is an option of the generic
`GRAMGridAtmosphereModel` only. A named preset fixes the default policy, and
`surrogate_preset_model` rejects grid options with an `ArgumentError`. To apply
another policy to a preset's grid, build the generic model from the resolved
file; the result carries no named-preset contract. In a session with
SpaceAGORA and GRAMSuite loaded:

```julia
preset = resolve_surrogate_preset("odyssey_p20_frozen_v1"; version="1.0.0")
density_model = GRAMGridAtmosphereModel(
    planet="Mars",
    surrogate_file=preset.file,
    expected_sha256=preset.expected_sha256,
    above_grid=:vacuum,
)
```

Queries below the grid floor fail under every policy. The adapter does not
apply the native model's entry-interface polynomial, its fixed 2000 km cutoff
or its lower-altitude clamp.

The grid model ignores `SPACEAGORA_VACUUM_GRAM_CACHE` and
`SPACEAGORA_DENSITY_FREEZE_PER_STEP`, and it evaluates the current coordinates
even when a buffered sample has the same timestamp. Native GRAM track caches,
isolated pools and per-satellite native copies do not apply to it. Threads can
query one model at the same time because the queries only read it; keep its
arrays and metadata unchanged during a run. An ordinary `deepcopy` produces
independent arrays, including when configuration isolation copies a whole
run. Managed process workers have their own native GRAM warm-up at startup,
which this model does not use.

### Height-based composition of two frozen atmospheres

`CombinedAtmosphereModel` joins two frozen, native-free atmosphere snapshots by
height above the reference ellipsoid: below the handover the lower model
answers, and at or above it the upper model does. With the published Mars
presets it serves a descent from the upper atmosphere into the lower one:

```julia
using SpaceAGORA
import GRAMSuite
lower = surrogate_preset_model("mars_global_near_surface_p20_frozen_v1"; version="1.2.0")
upper = surrogate_preset_model("mars_global_upper_p20_frozen_v1"; version="1.0.0")
density_model = CombinedAtmosphereModel(lower, upper; handover_height_m=80e3)
```

This pairing is also available as a named preset, which carries its validated
use; see the [Mars descent preset](#Mars-descent-preset).

Hand over at 80 km above the ellipsoid, the upper preset's floor. Near-surface
versions 1.1.0 and 1.2.0 reach 81 km areoid height, above that floor
everywhere. Version 1.0.0 stops at 75 km and leaves a gap below the floor.

The composition covers exactly what its components cover. There is no
blending, extrapolation or fallback between the two, so a query that the
selected component refuses fails with that component's `DomainError`. Below
the handover, coverage stops 5 m above supported terrain and excludes
latitudes beyond 85 degrees, volcano flanks and positions where a component is
unavailable, even though the upper preset covers those columns above 80 km.

Each component returns its own wind. Near-surface versions 1.0.0 and 1.1.0
store none, so with them the wind is zero below the handover and the upper
preset's stored wind at and above it, and it changes abruptly there. Version
1.2.0, a format 2 payload, returns its stored winds below the handover, but
the two components' winds still differ at the handover. Density and
temperature can also step by the difference between the presets.

Both presets are frozen at the same instant. A composition built with
`CombinedAtmosphereModel` adds no time evolution and has no validated accuracy
of its own; each preset's documented validation applies on its side of the
handover.

The constructor rejects components whose recorded planet, frozen instant or
reference ellipsoid differ, and a grid component whose height range does not
contain the handover. It cannot check that a near-surface lower component
reaches the handover, because that component's top is an areoid height. That
check, and choosing components whose validation suits your use, are up to
you. Only `GRAMNearSurfaceAtmosphereModel` and `GRAMGridAtmosphereModel` can
be composed; any other model raises an `ArgumentError`.

The engine treats the composition like its components. To propagate below
50 km, give the spacecraft a touchdown specification (see
[Terrain contact for a landing controller](stop_conditions.md#Terrain-contact-for-a-landing-controller));
otherwise the engine's default 50 km impact stop ends the descent there. Leave
a margin between the touchdown height and the lower component's minimum
clearance (5 m for the near-surface preset): the solver can evaluate a state
beyond the touchdown surface before it locates the crossing, and a state below
that clearance fails with the component's `DomainError`.

`atmosphere_provenance(density_model)` records the handover rule and both
presets' provenance.

### Mars descent preset

`mars_global_descent_p20_frozen_v1` 1.0.0 names the composition validated for
descents: near-surface 1.2.0 below 80 km above the ellipsoid and the global
upper preset 1.0.0 at and above it, frozen at the Odyssey P20 instant. It has
no data of its own. `surrogate_preset_model` installs and verifies each
component through its own catalog entry, exactly as when that preset is
selected alone, and composes them as described above:

```julia
using SpaceAGORA
import GRAMSuite
density_model = surrogate_preset_model("mars_global_descent_p20_frozen_v1"; version="1.0.0")
```

From 365 km down to 80 km above the ellipsoid, the domain is the upper
preset's: all latitudes and longitudes. Below 80 km it is the near-surface
preset's: planetocentric latitudes within 85 degrees, surface below 9 km, at
least 5 m above the surface, and where its stored winds are available. A
descent that drops below 80 km over a polar cap or a volcano flank fails there
with the near-surface preset's `DomainError`. At the 80 km handover, density,
temperature and winds change abruptly from one preset's values to the
other's; nothing is blended.

The composition was validated against native Mars-GRAM at the frozen instant,
under the lab's release limits: on 48 fixed descent paths and in 24
propagated descents in a reference integrator, and in 23 ballistic entries in
SpaceAGORA's engine with native Mars-GRAM in the same engine loop. The
descents started 125 km above the ellipsoid at 5.5 km/s relative to the
planet, with a ballistic coefficient of 100 kg/m² and flight-path angles from
−10 to −20 degrees, and ended 100 m above the terrain. Outside that envelope,
each preset's own validation applies on its side of the handover.

The validation does not establish passes that skim or recross the handover,
as in aerobraking; other entry speeds, ballistic coefficients and vehicles;
lifting or guided entries; real-world landing accuracy; the engine's use of
the near-surface speed of sound; or other epochs.

`atmosphere_provenance(density_model)` names the composite, its version, its
components and the catalog, with both presets' own provenance, and the
results-bundle manifest records it under `atmosphere`.

The composite takes no `file`, and `resolve_surrogate_preset` refuses it. To
compose local payload files, build each preset from its file and use
`CombinedAtmosphereModel`; that composition carries no named contract.

## Native GRAM and surrogates

The two workflows in this section need native GRAM, and so the NASA GRAM Suite
distribution; see [GRAMSuite Setup](gramsuite_setup.md).

### Compare a surrogate with native GRAM

To see how a surrogate differs from native GRAM in your own scenario, run the
same case twice and change only the density model: once with the named
preset, and once with a `GRAMAtmosphereModel`. Two native references answer
different questions. Native GRAM evaluated at the preset's frozen instant
isolates the grid's interpolation error. Native GRAM run with the actual time
along the trajectory adds the error of freezing the atmosphere.

The Odyssey preset has accepted diagnostic and propagation comparisons against
the second reference, along bounded P20 passages. These comparisons are
evidence, not numerical release limits: the catalog leaves application
accuracy requirements unset. At arbitrary points away from those passages,
pointwise differences can be larger, particularly near the 260 km ceiling at
high northern latitudes.

`atmosphere_provenance(model)` lists how the preset was generated: planet,
frozen UTC instant, Mars-GRAM configuration and input identities. Check each
setting against the native model you build. A plain `GRAMAtmosphereModel` is
not configured the same way automatically, and the Odyssey preset's diagnostic
comparisons used a dedicated native adapter matched to the preset's settings.
Compare pointwise density, temperature and wind along the native trajectory,
and the per-pass quantities your algorithm depends on, such as drag delta-v,
heat load and peak heat rate.

### Prepare another frozen snapshot

A named preset covers one planet, one frozen instant and one bounded domain.
For another scenario, generate a new frozen grid with GRAMSuite's
recorded-recipe generator, following its
[grid generation guide](https://github.com/Space-FALCON-Lab/GRAMSuite.jl/blob/main/docs/grid_generation.md).
The guide starts with a dry run of the recipe that needs no native GRAM; only
the final generation does. Validate the new grid against native GRAM for the
intended scenario before relying on it, then load it with the generic model.
This is a template; fill in your file's path and SHA256:

```julia
density_model = GRAMGridAtmosphereModel(
    planet="Mars",
    surrogate_file="/path/to/new_grid.jls",
    expected_sha256="<sha256 of the file>",
)
```

A grid loaded this way carries no named-preset contract: its domain, epoch and
accuracy are whatever your own generation and validation established.
