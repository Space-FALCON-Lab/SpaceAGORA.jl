"""
    GRAMNearSurfaceAtmosphereModel(; planet="Mars", surrogate_file, expected_sha256=nothing)

SpaceAGORA adapter for a native-free near-surface Mars atmosphere: density and temperature from 5 m above the local
terrain to the payload's areoid-height top (75 km in version 1.0.0 of the published preset, 81 km in 1.1.0), following
native Mars-GRAM's own near-surface rule. Load `GRAMSuite` with its near-surface API before keyword construction. Its loader validates the
payload and owns the contained arrays; no native GRAM installation is needed. `core` holds the
`GRAMSuite.GRAMNearSurfaceAtmosphereModel`.

`getDensity` takes height above the reference ellipsoid in metres, geodetic latitude and east longitude in radians,
and elapsed time in seconds. It returns density in kg/m^3, temperature in kelvin and the east, north and up wind in
m/s.
- A format 1 payload (the published versions 1.0.0 and 1.1.0) stores no winds: the wind vector is zero, so a simulation
  with this model has no atmospheric wind.
- A format 2 payload returns its stored winds, with the horizontal components clipped at 0.7 times the speed of sound
  as native clips them. As for the stored winds of grid presets, the `wind` argument does not suppress them; the engine
  masks winds when its environment disables them.

Elapsed time is not applied to the frozen snapshot. Pressure, the regime and the surface-layer model status are
available from `GRAMSuite.near_surface_state(model.core, lat_deg, lon_deg, h_m)`. With a wind layer,
`GRAMSuite.near_surface_wind_state` adds the winds before clipping, the speed of sound and the wind regime.

Queries outside the supported domain throw `DomainError` naming the reason: planetocentric latitude beyond the
payload's limit, volcano flanks, clearance below the minimum, areoid height above the top, or an unavailable component.
There is no extrapolation and no native fallback. Trajectory-density caches and per-step density freezing are
bypassed for this model, as for fixed grids, so current coordinates always reach the evaluator. Treat `core`, its
arrays and metadata as read-only during a run; concurrent queries share the snapshot.
"""
struct GRAMNearSurfaceAtmosphereModel{C} <: AbstractDensityModel
    core::C
end

# Whether a near-surface model's payload has a wind layer. The GRAMSuite extension answers for real cores.
_near_surface_model_has_winds(model) = false

# Native-free snapshots that must see every current coordinate (their domain checks and spatial variation make a
# trajectory spline or a per-step frozen sample unsafe): fixed grids and near-surface payloads. `CombinedAtmosphereModel`
# joins two of them; `_NativeFreeSnapshotModel` (defined with it) covers the components and their combinations.
const _NativeFreeSnapshotComponent = Union{GRAMGridAtmosphereModel, GRAMNearSurfaceAtmosphereModel}
