"""
    GRAMNearSurfaceAtmosphereModel(; planet="Mars", surrogate_file, expected_sha256=nothing)

SpaceAGORA adapter for a native-free near-surface Mars atmosphere: density and temperature from 5 m above the local
terrain to the payload's areoid-height top (75 km for the published preset), following native Mars-GRAM's own
near-surface rule. Load `GRAMSuite` with its near-surface API before keyword construction. Its loader validates the
payload and owns the contained arrays; no native GRAM installation is needed. `core` holds the
`GRAMSuite.GRAMNearSurfaceAtmosphereModel`.

`getDensity` takes height above the reference ellipsoid in metres, geodetic latitude and east longitude in radians,
and elapsed time in seconds. It returns density in kg/m^3, temperature in kelvin and a zero wind vector: the payload
stores no winds, so a simulation with this model has no atmospheric wind. Elapsed time is not applied to the frozen
snapshot. Pressure, the regime and the surface-layer model status are available from
`GRAMSuite.near_surface_state(model.core, lat_deg, lon_deg, h_m)`.

Queries outside the supported domain throw `DomainError` naming the reason: planetocentric latitude beyond the
payload's limit, volcano flanks, clearance below the minimum, areoid height above the top, or an unavailable component.
There is no extrapolation and no native fallback. Trajectory-density caches and per-step density freezing are
bypassed for this model, as for fixed grids, so current coordinates always reach the evaluator. Treat `core`, its
arrays and metadata as read-only during a run; concurrent queries share the snapshot.
"""
struct GRAMNearSurfaceAtmosphereModel{C} <: AbstractDensityModel
    core::C
end

# Native-free snapshots that must see every current coordinate (their domain checks and spatial variation make a
# trajectory spline or a per-step frozen sample unsafe): fixed grids and near-surface payloads. `CombinedAtmosphereModel`
# joins two of them; `_NativeFreeSnapshotModel` (defined with it) covers the components and their combinations.
const _NativeFreeSnapshotComponent = Union{GRAMGridAtmosphereModel, GRAMNearSurfaceAtmosphereModel}
