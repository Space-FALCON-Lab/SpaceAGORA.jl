# Kuang's later-task fixtures

Self-contained follow-up to the small known-answer scheduler check. All new code,
generated Markdown and numerical results stay in this directory. The parent
SpaceAGORA package is used without changing its code or dependency manifests.

## Scheduling policy explained

Read [SCHEDULING_POLICY.md](SCHEDULING_POLICY.md) for the backlog scoring formula,
terminal matching, the two-sender switching timeline, relay outflow credit,
score-versus-transfer-rate distinctions, and policy limitations.

## Saved animations of the six backlog cases

Open [animations.html](animations.html) in a browser to play all six saved MP4s.
[ANIMATIONS.md](ANIMATIONS.md) links the individual videos, initial/final PNGs,
measured totals and rendering details. Each case's `animation/` subdirectory
contains its video, snapshots, frame-to-simulation-time CSV and validation record.
These are the original **satellite-to-satellite** cases, not the ground-station
experiment below.

The renderer adapts the GLMakie `Observable`, `Axis3`, `record` and headless-Xvfb
approach from `Summary_Table_for_Other_Labs`. It replays the saved data without
rerunning or changing any simulation. The enlarged directed network is a
schematic; the separate 3D panel uses the saved orbital positions. Brief review
holds make every recorded selection change visible, including the busy relay's
very short initial direct connection. Holds are labelled and do not lengthen
the simulated contacts.

From the repository root, render all six using the existing reference
environment and Earth texture, which are read without modification:

```sh
set -o pipefail
xvfb-run -a julia --startup-file=no \
  --project="../Summary_Table_for_Other_Labs/2_SpaceAGORA.jl" \
  kuang_later_tasks/animate_backlog.jl 2>&1 |
  tee kuang_later_tasks/animations.log
```

Append a case name for a focused rerender; run the full command to refresh the
gallery and index. Render serially. Unlike the numerical runner, this command
requires the reference project's plotting environment and `xvfb-run`.
`results/backlog/animations_status.toml` records the render status. The renderer
checks the displayed states, unchanged inputs, H.264 format, dimensions, frame
count and duration, and fully decodes each saved video with FFmpeg.

Using the local browser installation described below, independently check that
all six videos load, play and seek in Chromium:

```sh
node kuang_later_tasks/browser_check/check_animations.cjs
```

This saves [browser evidence](results/backlog/animation_browser_validation.json)
and a [gallery screenshot](results/backlog/animation_gallery.png). Rerun this
check after rendering; the evidence includes video, frame-map and gallery hashes.
An edge label such as `2.2->3.2` means spacecraft 2's terminal 2 to spacecraft 3's
terminal 2. Link scores and maximum rates are not measurements of carried data.

## Follow-up: three equatorial ground stations

Run the separate Earth-downlink experiment from the repository root:

```sh
julia --startup-file=no --project=. kuang_later_tasks/ground_stations.jl
```

It uses the original forty-spacecraft orbital geometry with data stores and
120 allowed downlinks to three ground stations at latitude 0 degrees and
longitudes 0, 120 E and 120 W. Each station has one antenna and a 10-degree
elevation mask. The data now reaches Earth, not a destination spacecraft.

Read [GROUND_REPORT.md](GROUND_REPORT.md) and open the self-contained
[ground viewer](results/ground_stations/viewer.html). The viewer is an
Earth-fixed 2D diagnostic display, not the optional SpaceAGORAVisualizer 3D
renderer. Its source template is [ground_viewer.html](ground_viewer.html);
open the generated page under results, which includes the actual data.
Geometry, matching, per-station delivery, accounting, scene and native viewer
link-payload checks are included. Results and a separate status file are under
[results/ground_stations/](results/ground_stations/).

The original report below describes the earlier experiments only; the follow-up
report records the subsequently added ground-station work. Original momentum and
busy-relay fixtures are retained unchanged for comparison.

### Browser validation

The simulation needs only the root Julia environment. Browser tests use an
isolated Playwright installation under this folder; the root dependency
manifests are not changed. From the repository root:

```sh
cd kuang_later_tasks/browser_check
npm ci
PLAYWRIGHT_BROWSERS_PATH="$PWD/.browsers" npx playwright install chromium --only-shell
npm test
```

The test opens the generated HTML directly from disk, with no server. It checks
all recorded samples, line endpoints, counters and controls, and saves
[VIEWER_CHECK.md](VIEWER_CHECK.md), a screenshot and a JSON validation record.
The JSON contains the viewer's hash so a later simulation cannot silently reuse
an unrelated browser result. The integrated browser could not access this host,
so verification was performed in local headless Chromium instead.

## Run

From the repository root:

```sh
julia --startup-file=no --project=. kuang_later_tasks/run.jl
```

Or, from this directory:

```sh
julia --startup-file=no --project=.. run.jl
```

To retain the console output from the repository root:

```sh
set -o pipefail
julia --startup-file=no --project=. kuang_later_tasks/run.jl 2>&1 |
  tee kuang_later_tasks/run.log
```

The existing root Julia environment is sufficient. There is no plotting,
PythonCall, NetworkX or native GRAM dependency for these fixtures. No installation
or environment activation is performed by the runner. Run serially: another
invocation would overwrite the same output files.

Successful execution writes:

- [REPORT.md](REPORT.md): measured results, assumptions, numerical tolerances,
  interpretation and the handoff to Nihal.
- [results/status.toml](results/status.toml): `running`, `failed` or `passed`.
  Treat an older report as stale if this status is not `passed`.
- [results/provenance.toml](results/provenance.toml): Julia/package versions,
  checkout revision, existing worktree status and code/project/manifest hashes.
- [results/constellation/](results/constellation/): a quarter-orbit forty-spacecraft
  run, an otherwise identical no-link control, candidate scores, independently
  checked matchings, recorded windows, Cartesian states and SMA comparisons.
- [results/backlog/](results/backlog/): six short data-link runs, static prior-outflow
  probes, exhaustive matching comparisons, per-node and whole-network accounting,
  hand calculations, delivery and switching metrics, and step-refinement results.

The runner uses `Test` assertions, exits unsuccessfully on failure, and writes the
report only after all checks pass. Standard wide simulation output is disabled;
the saved CSVs deliberately contain the states and evidence needed for this task.
Schedule CSVs come from the actual simulation's `link_schedule` history.

## Experiments

| Fixture | Purpose |
|---|---|
| Forty-spacecraft + no-link control | Original 30-helper / 10-target configuration, shortened to 0.25 target orbits; animation off; every candidate checked at accepted steps |
| Empty producer | Demonstrate a known continuous delivery rate from an empty producing store |
| Production competition | Demonstrate that positive production does not guarantee selection; record switching, delivery and remaining backlog |
| Busy relay | Demonstrate previous feasible outflow credit, preference for a relay over a direct route, and explicitly permitted local-production overflow |
| Direct-route control | Same stores/rates/policy, with relay input deliberately excluded; compare destination delivery and loss |
| Refined competition and relay | Repeat with half the maximum step to expose scheduler-resolution sensitivity |
| Three static outflow probes | Distinguish no prior link, a feasible prior link and an out-of-range prior link |

The original data cases have finite stores at every endpoint, including the destination.
All candidate terminal graphs are forests. In the forty-spacecraft run, the full
allowed graph is cyclic, but its feasible terminal components are checked to be
disjoint edges at every accepted step. These runs therefore do not exercise the
optional cyclic matching backend.

The independent backlog calculations use
`min(rate, stored_tx/horizon + production_tx, room_rx/horizon + previous_feasible_outflow_rx)`.
The replay audit reconstructs store flags using the package's internal
`DataStorage._start_data_state!` helper on **copies**, so the simulation is not
modified. That narrowly scoped helper is version-sensitive; provenance records
the tested checkout. The expected scores, exhaustive matching enumeration,
analytic rates and balance equations do not call the engine's scoring or flow
implementation.

For focused reruns during development:

```sh
julia --startup-file=no --project=. -e \
  'include("kuang_later_tasks/run.jl"); KuangLaterTasks.run_constellation()'
julia --startup-file=no --project=. -e \
  'include("kuang_later_tasks/run.jl"); KuangLaterTasks.run_backlog_suite()'
```

Focused reruns do not refresh the top-level report/status; finish with the full
runner before sharing results.

## Scope and interpretation

- Positions and velocities are inertial Cartesian SI values. Attitude is disabled.
- Momentum links are ideal, range-only links. Data links check spherical planet
  occultation; their fixed rates and virtual terminals isolate scheduling.
- Matching changes and contact windows are sampled at accepted solver steps, not
  continuously located geometry events.
- Internal transfers cancel in the network balance. Because the destination is
  also a store, subtracting both destination storage and delivered bits would
  double-count delivery.
- Overflow is permitted only for the busy relay's own production, and its amount
  is checked against the recorded mode durations. A receiver staying below
  capacity is not proof of no loss elsewhere.
- The routing comparator is not an end-to-end optimizer, and the engine's
  aggregate counters do not identify which source's packets arrived.
- Give Nihal the accounting and hand-comparison tables for review. Ground-station
  geometry and a diagnostic viewer are now covered by the separate follow-up above.
  The optional 3D renderer, genuinely cyclic feasible graphs, mounted-beam physics
  and resumed laser history remain separate work.

Reference task materials:
[interlinks](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/dafb5f10de8c7497efe3e898e9e2630a458b7a7c/docs/src/user/interlinks.md),
[forty-spacecraft example](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/dafb5f10de8c7497efe3e898e9e2630a458b7a7c/examples/targeted_40sat_sma_interlinks.jl).
