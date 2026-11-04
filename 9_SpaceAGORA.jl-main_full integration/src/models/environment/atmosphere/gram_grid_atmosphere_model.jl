"""
    GRAMGridAtmosphereModel(; planet, surrogate_file="", search_roots=String[],
                           expected_sha256=nothing, above_grid=:error,
                           vacuum_temperature=200.0)

SpaceAGORA adapter for a fixed GRAM surrogate grid. Load `GRAMSuite` with its
native-free grid API before keyword construction. Its grid loader validates the
payload and owns the contained arrays; no native GRAM installation is needed.
`core` holds the `GRAMSuite.GRAMGridAtmosphereModel` snapshot.

`getDensity` takes altitude in metres, latitude/longitude in radians, and elapsed
time in seconds. It returns density in kg/m^3, temperature in kelvin, and the
stored local east/north/up wind in m/s. Elapsed time and the wind selector do not
change a fixed snapshot. Match the grid's documented coordinates and atmospheric
epoch to the mission before use; loading a grid does not establish its physical
accuracy or recover missing provenance.

Altitude/latitude queries outside grid coverage throw `DomainError`. The explicit
`above_grid=:vacuum` option permits zero density above the grid ceiling only.
At an exact latitude pole, stored components follow the requested longitude
label; this API does not establish a unique physical wind vector there.
The adapter applies no entry-interface polynomial, fixed 2000 km cutoff, or
lower-altitude clamp. Trajectory-density caches and per-step density freezing
are bypassed for this model so current coordinates reach the grid evaluator.

Treat `core`, its arrays, and metadata as read-only during a run. Concurrent
queries share the snapshot; ordinary `deepcopy` makes an independent copy.
This model does not use native GRAM pools or per-satellite native instances.
"""
struct GRAMGridAtmosphereModel{C} <: AbstractDensityModel
    core::C
end
