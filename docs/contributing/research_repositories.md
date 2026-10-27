# Research repository ownership

Paper and mission studies are independent consumers of SpaceAGORA. Small
teaching examples, regression fixtures and the reusable benchmark and
verification tools stay in core; the [user guide](../src/user/studies_benchmarks.md)
explains which entry point to use.

The four research repositories are private under Space-FALCON-Lab, and access
is by invitation. Public core documentation names them without linking to
pages most readers cannot open. Reusable replication capabilities stay in the
public core.

| Repository | Owns |
|---|---|
| SpaceAGORAAtmosphereStudies | Atmosphere and entry-density research |
| SpaceAGORAParallelizationStudies | Paper campaigns and scaling investigations |
| SpaceAGORAMissionReconstructions | Mission reconstruction and telemetry fitting |
| SpaceAGORARPOStudies | RPO paper reconstruction |

## Source and pins

Each repository was created by importing study source from core. Its import
records keep the original Git blobs, the license and every intentional change
of path or environment, and its first commits identify the transferred
source. Each repository pins the core revision it runs against in its own
`SpaceAGORASource.toml`, and the repositories may pin different revisions.
Core never depends on a research repository.

The study environments still use some private core names. Their source and
environment manifests document that boundary; using a name there does not make
it public API, and it does not make a newer core commit compatible. Moving a
study's core pin needs the affected scientific checks run on the old and the
new pin, with any native or data case that could not be run stated
explicitly.

## Tests and callers

- Core keeps the generic simulation, planner, actuator, telemetry-parser and
  parallel-execution tests. Small attributed fixtures cover the generic
  contracts that used to be exercised by importing a mission driver.
- Paper-only gates and tests live with their study. The parallelization
  repository keeps the paper workload coverage gate and its controls; the RPO
  repository keeps the saved paper inputs, the swarm-policy checks and the
  refused-plan evidence.
- The mission repository tests the actual study producers. Core's smaller
  fixtures do not show that a later edit to an external study still produces
  the same data.
- One example shard of core CI (`.github/workflows/julia-ci.yml`, step
  "Extracted pinned Mars RAAN native smoke", run by
  `scripts/ci/run_mission_native_smoke.py`) checks the Mars RAAN study from
  SpaceAGORAMissionReconstructions at a fixed revision, when the
  `SPACEFALCON_REPO_READ_TOKEN` secret can read that repository. It uses the
  study's pinned core and the public GRAM wrapper, and links the native assets
  the job already has. Without a credential the job reports that native
  coverage was not exercised. With one, access, asset, setup and native
  failures fail the job; they never count as a skip or as passing coverage.
- The core CLI's `benchmark runtime-analysis`, `benchmark smart-parallel-ladder`
  and `telemetry` commands use core-owned tools. Paper workload preparation and
  the paper's remote recipes belong to the parallelization repository.

## Maintenance rules

- Make the replacement in the research repository reachable before deleting a
  core caller.
- Keep the recorded source commits reachable from each research repository's
  main history. Merge imports with merge commits; before rewriting history,
  update every recorded pin and qualify it again.
- The hosted CI checks are still required. Native qualification depends on
  private-repository credentials, so report its result separately from the
  local deterministic checks.
