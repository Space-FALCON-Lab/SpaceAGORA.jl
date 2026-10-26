# Research repository ownership

Paper and mission studies are independent consumers of SpaceAGORA. Small teaching
examples, regression fixtures and the reusable benchmark and verification tools
stay in core. The [user guide](../src/user/studies_benchmarks.md) selects the right
entry point.

All four research repositories are private under Space-FALCON-Lab. Lab members
retain access; outside collaborators and students need an invitation. Public core
documentation names the repositories without presenting inaccessible links.
Reusable replication capabilities remain in the public core.

## Source transfer

The study snapshot comes from core
`bfa32f6a44c5f242714b8cf182a451e863119a05`. Each destination records original Git
blobs, the license and intentional path/environment adaptations. The initial
import commits identify the transferred source. Later commits on each import adapt
it (study folder names and wording, and the parallelization studies' core pin, now
the v0.6.0 release); the table records the reviewed import heads below. Neither is a
floating dependency pin for core runtime.

| Repository | Owns | Imported source | Current head |
|---|---|---|---|
| SpaceAGORAAtmosphereStudies | Atmosphere and entry-density research | `32ca5be9d6b6fe4e865365a470d03a2a56b3a1fe` | `32ca5be9d6b6fe4e865365a470d03a2a56b3a1fe` |
| SpaceAGORAParallelizationStudies | Paper campaigns and scaling investigations | `0071dd9b6e0f94828fb4da8cfcddc8425c403abd` | `f5c57d52a7d0662771e38e2e29b3cd0148d6ce52` |
| SpaceAGORAMissionReconstructions | Mission reconstruction and telemetry fitting | `ff67a512face495afa6dec10171806209a62033d` | `64bb98ff16ce0610589fd24bfd3e4e44fb805f51` |
| SpaceAGORARPOStudies | RPO paper reconstruction | `9356eacc6de660947133c8574fc652bae3b8e58b` | `10a1518153678bb7347103bca6360f7a610d3cc7` |

## Tests and callers

- Core retains generic simulation, planner, actuator, telemetry-parser and
  parallel-execution assertions. Minimal attributed fixtures preserve the generic
  contracts formerly exercised by importing a mission driver.
- Paper-only gates and tests move with the complete study. The parallelization
  project keeps the paper workload coverage gate and its controls. The RPO
  project keeps saved-paper-input, swarm-policy and refused-plan evidence.
- The mission project tests the actual study producers. Core's smaller fixtures
  do not claim that a later edit to an external study still produces the same data.
- The automatic native example job checks the transferred Mars RAAN study at an
  immutable revision when `SPACEFALCON_REPO_READ_TOKEN` can read its private
  repository. It uses the study's pinned core and public GRAM wrapper and reuses
  provisioned assets. Without a credential, the job explicitly reports that native
  coverage was not exercised. With a credential, access, asset, setup and native
  failures fail the job; they never count as a skip or passing coverage.
- The core CLI's runtime-analysis, smart-parallel-ladder and telemetry commands
  continue to use core-owned tools. Paper workload preparation and paper remote
  recipes belong to the parallelization project.

The study environments intentionally retain some private core dependencies.
Their source and environment manifests document that boundary. This extraction
does not promote those names to public APIs or make arbitrary newer core commits
compatible. A pin update needs the affected original-versus-new scientific
checks, with unavailable native/data cases stated explicitly.

## Merge and maintenance order

Make the replacement repository source reachable before deleting a core caller.
Keep the documented source commits reachable from each research repository's
main history. Use merge commits for these imports, or update all recorded pins
and requalify them before choosing a history rewrite.

The core relocation is composed with the examples environment, optional plotting,
HYPR-loader retirement and public documentation work. Those changes must be
integrated before this relocation is merged into main. Hosted checks remain
required. Native qualification is conditional on private-repository credentials;
its result must be reported separately from local deterministic checks.
