# Three equatorial ground stations

Generated: 2026-11-09T01:33:52.455 UTC.

This follow-up supplements, rather than replaces, [the original results](REPORT.md). Data now reaches actual ground-station endpoints, not a destination satellite.

## Configuration

- Latitude 0 degrees; longitudes 0, 120 E and 120 W; altitude 0 m. Equal 120-degree longitude spacing was checked numerically.
- One receiving antenna per station, 10-degree minimum elevation, 5000 km range. Stations rotate with Earth using the simulation ephemerides frame. The mask enforces horizon clearance; station links do not use the spacecraft-to-spacecraft occultation test.
- Same 40 spacecraft orbital geometry as the original example: 30 at 1000 km and 10 at 1100 km, circular equatorial orbits; 227 kg each. Each has one virtual data terminal, a 1 Gbit store initially holding 20 Mbit, and 0.1 Mbit/s continuous production.
- All 120 satellite-to-station downlinks are allowed; no intersatellite links or momentum forces in this separate data experiment. Fixed rate 1 Mbit/s; backlog horizon 60 s; no radio link-budget or battery limitation is asserted.
- Duration 1608.9445727393636 s (quarter target orbit), Tsit5, dtmax 5 s, reltol 1e-10, abstol 1e-8. No movie is generated.

## Measured delivery

| Station | Longitude | Received (Mbit) | Hand lower (Mbit) | Hand upper (Mbit) |
|---|---:|---:|---:|---:|
| Equator_0 | 0.0 | 1471.1580125333287 | 1471.1580094348135 | 1471.1580307337301 |
| Equator_120E | 120.0 | 1454.207105938022 | 1454.2070982390965 | 1454.2071473397873 |
| Equator_120W | -120.0 | 1463.0785622114495 | 1463.0785576969145 | 1463.0785993580544 |

Total Earth delivery: **4388.443680682801 Mbit**. 3924 scheduling samples; 3771 selected-set changes. All three stations receive data. No production loss occurred to the 1-bit tolerance.

## Checks and interpretation

- Independently computed range/elevation for all 120 candidates at every accepted step. Separate tests at -5, 9, 11 and 90 degrees verify blocking and the 10-degree mask at every station.
- Actual availability and selected scores agree with those calculations. A spacecraft can see at most one station at once here, so feasible components are stars. The best matching is one maximum positive-score edge per station. This is checked without a Python cyclic matcher. Visible does not mean selected.
- Independent interval upper transfer: min(rate * dt, starting stored bits + production * dt), summed only over actually selected accepted-step intervals. When a store starts within its documented 1e-9-capacity empty threshold (1 bit here), the lower transfer is production * dt: a tiny positive residual may already be flagged empty. These interval uncertainties are accumulated explicitly, not hidden in solver tolerances. Every station counter lies between the resulting lower and upper integrals, with 1 bit numerical allowance. Final uncertainty widths: [21.29891668067939, 49.10069089001681, 41.66113997084538] bits.
- At every saved step: per-spacecraft balance, bounds, zero loss, and sum(spacecraft sent) = sum(ground received). Network balance is initial + produced = remaining spacecraft storage + Earth delivery + loss. Maximum node residual: 4.3259933590888977e-7 bit. Ground delivery is now outside the spacecraft storage total, unlike the earlier destination-satellite fixture.
- Contact boundaries and switches have accepted-step resolution, not exact continuously located visibility crossings. A 5-second step cap is not a claim that this scheduling behavior is timestep-converged. Busy-relay comparison and its refinement remain in the original separate fixture.

## Viewer

Open [viewer.html](results/ground_stations/viewer.html). This is a self-contained, Earth-fixed north-polar **2D diagnostic viewer**, not the optional SpaceAGORAVisualizer 3D application. Green markers are stations; orange lines are active downlinks; dashed green lines are available but idle. The slider uses every accepted scheduling sample, so it cannot silently skip short recorded windows. Playback advances samples, not uniform physical time. At the final sample it shows the final recorded selection, without extrapolating beyond the run.

The repository-native [viewer link payload](results/ground_stations/interlink_viewer_payload.json) and [scene](results/ground_stations/scene.json) are also saved. Their station coordinates and contact/active intervals were checked against history. The optional 3D renderer, terrain, weather, antenna pointing, atmospheric losses, power shedding and resumed histories are not validated by this fixture.

The separate browser check writes [VIEWER_CHECK.md](VIEWER_CHECK.md) and a SHA-256-bound validation record after checking the actual rendered links, coordinates, delivery counters and controls at every recorded sample. Rerun that check after regenerating this viewer; see [README.md](README.md).

## Evidence and rerun

- [Stations](results/ground_stations/stations.csv)
- [Candidates before propagation](results/ground_stations/candidates_before_run.csv)
- [Geometry probes](results/ground_stations/geometry_probes.csv)
- [All candidates and selections](results/ground_stations/candidate_scores.csv)
- [Recorded windows](results/ground_stations/link_schedule.csv)
- [Station delivery](results/ground_stations/station_delivery.csv)
- [Node accounting](results/ground_stations/node_accounting.csv)
- [Network accounting](results/ground_stations/network_accounting.csv)
- [Settings and hashes](results/ground_stations/settings.toml)

From the repository root:

```sh
julia --startup-file=no --project=. kuang_later_tasks/ground_stations.jl
```

