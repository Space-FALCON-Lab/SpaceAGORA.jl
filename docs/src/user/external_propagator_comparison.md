# Comparing Against an External Propagator

`examples/Earth_Harmonics_LatLon_Comparison.jl` is the reference case for
checking SpaceAGORA against a propagator written somewhere else — the MATLAB
one this page was written for, or any other. It starts from a geodetic
latitude and longitude, flies a spherical-harmonics gravity field while the
planet rotates underneath it, and writes both the trajectory and a set of plots
to a local folder.

Run it with:

```text
julia --project=. examples/Earth_Harmonics_LatLon_Comparison.jl
```

Everything lands under `output/harmonics_latlon_comparison/`: the trajectory as
`simulation_results.csv` and `simulation_results.feather`, and one PNG per saved
value under `plots/`.

## The case block

The top of the script is a single block of constants, and nothing outside it
needs touching for an ordinary comparison:

```julia
const INITIAL_LATITUDE_DEG = 0.0        # geodetic latitude of the start point
const INITIAL_LONGITUDE_DEG = 0.0       # planet-fixed longitude of the start point
const INITIAL_ALTITUDE_M = 500e3        # above the reference ellipsoid
const INCLINATION_DEG = 45.0            # orbit inclination through that point
const DESCENDING = false                # true selects the southbound crossing
const FLIGHT_PATH_ANGLE_DEG = 0.0       # 0 starts at an apsis
const ECCENTRICITY = 0.0                # 0 is circular; the start is periapsis

const HARMONICS_DEGREE = 40
const HARMONICS_ORDER = 40

const INITIAL_EPOCH = InitialTime(year=2026, month=9, day=8, hour=5, minute=0, second=0.0)
const MISSION_TIME_S = 86400.0          # one day
const DATA_RATE_S = 10.0                # seconds between saved rows

const SAVED_VALUES = (:orbital_elements, :gravity_accel)
```

Change `INITIAL_LATITUDE_DEG` and `INITIAL_LONGITUDE_DEG` and rerun; the
initial condition, the plots and the CSV all follow. The script derives the
speed itself, from the start radius and `ECCENTRICITY`, so a change of altitude
or latitude does not need a hand-computed velocity to go with it.

## Matching the two sides

Longitude is planet-fixed, so the start point only means something together
with an epoch — the same latitude and longitude an hour later is a different
inertial state. Three inputs have to agree for the comparison to be about the
propagators rather than about their setup:

1. **The epoch.** `INITIAL_EPOCH` here, and whatever seeds the planet rotation
   on the other side.
2. **The geodetic point.** Latitude, longitude and altitude above the reference
   ellipsoid — not above a sphere. SpaceAGORA places the point on the oblate
   ellipsoid, so a propagator using a spherical Earth starts tens of kilometers
   away at mid-latitudes.
3. **The field.** `HARMONICS_DEGREE` and `HARMONICS_ORDER`, against the same
   coefficient set. This example reads GGM05C from
   `data/Gravity_harmonics_data/EarthGGM05C.csv`, whose reference radius and GM
   are in the file's header and are not the same as the generic Earth radius
   used elsewhere in the run.

The script prints the inertial state the geodetic point resolved to, so the
other propagator can be seeded from the identical vector rather than from the
geodetic inputs:

```text
Comparison case
  epoch                 : InitialTime(2026, 9, 8, 5, 0, 0.0)
  geodetic start        : lat 0.0 deg, lon 0.0 deg, alt 500000.0 m
  inclination           : 45.0 deg (ascending)
  gravity field         : GGM05C 40x40
  inertial position (m) : [3.1995804387021046e6, 6.08863269581393e6, 0.0]
  inertial velocity (m/s): [-4765.050059355565, 2504.0369030361185, 5382.92698073559]
```

Seeding both sides from those two vectors isolates the propagators. Seeding
both from latitude and longitude instead also tests that the two agree about
the shape of the Earth and where the prime meridian was pointing — useful, but
a different question, and worth deciding which one you are asking.

## What gets plotted

Every name in `SAVED_VALUES` becomes one figure, with a panel per column. The
orbital elements come out as six panels: semimajor axis, eccentricity,
inclination, RAAN, argument of periapsis and true anomaly.

Adding a value is one edit:

```julia
const SAVED_VALUES = (:orbital_elements, :gravity_accel, :total_accel)
```

It is saved and plotted with no further changes — a field the script has no
labels for still gets a panel per component, titled by column name. To give it
proper axis labels, add a row to `PANEL_LABELS`, which maps a field name to one
`(title, ylabel, scale)` per column:

```julia
const PANEL_LABELS = Dict{Symbol, Vector{Tuple{String, String, Float64}}}(
    :orbital_elements => [
        ("Semimajor Axis", "a (km)", 1e-3),
        ...
    ],
)
```

The scale multiplies the column on its way to the plot, which is how the
semimajor axis is stored in meters and drawn in kilometers.

See [Simulation Outputs](outputs.md) for the full list of names
`SAVED_VALUES` accepts and for writing a `SaveField` of your own.

## Reading the results elsewhere

The CSV is the portable artifact; the feather file beside it holds the same
columns and loads faster. In MATLAB:

```matlab
T = readtable("output/harmonics_latlon_comparison/simulation_results.csv");
a_km = T.sc1_orbital_elements_1 / 1e3;
inc_deg = T.sc1_orbital_elements_3;
```

Column names and units are in [Simulation Outputs](outputs.md).

## What to expect

At degree and order 40 over a day, from a 500 km circular start at 45°
inclination, the secular signature is the one to check first: the node
regresses a few degrees per day, scaling with the cosine of the inclination,
while the semimajor axis and inclination oscillate about fixed means at twice
per orbit without drifting. Two propagators that agree on the field and the
epoch agree on that regression rate closely; a mismatch there usually means the
planet rotation or the epoch differs, not the gravity model.
