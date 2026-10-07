# MuJoCo proximity scenes: subpackage architecture (Stage 0 and 1)

Developer note for `packages/SpaceAGORAMuJoCo`. Design basis: the 2026-10-06 MuJoCo extension design memo
(coupling option c, chief frame, Encke gravity). This note records what Stage 1 built and the Stage 0
dependency decision.

## Packaging

The scene runner is a separate package in `packages/SpaceAGORAMuJoCo`, following `packages/SpaceAGORAHYPR`:
its own `Project.toml` (depends on `SpaceAGORA`, compat `0.2`), its own lazy `Artifacts.toml`, `src/` and
`test/`. The root `Project.toml` and `Artifacts.toml` are untouched, and nothing was added under `src/` or
`ext/`, so no ownership boundary moves and no contract gate scans the new code (the gates that list
`packages/SpaceAGORAHYPR/src` as a root do not include it). Putting `ext/` and this package under the
topology contract is Stage 2 work, together with the core hooks that integration needs.

## Stage 0: library choice

`MuJoCo_jll` (General registry) has only 2.3.7+0 and 3.1.6+0 (2024-06-21), and its Yggdrasil recipe is pinned
to 3.1.6, while upstream is at 3.15.0. The decision was a current version close to Basilisk's 3.11, so the
package pins the official 3.11.0 Linux x86_64 release tarball as a lazy artifact. Pkg checks the download
against the sha256 in `Artifacts.toml` (which equals the published `.sha256` and the GitHub asset digest) and
the unpacked tree against its git tree hash. `eq_active` is `mjtBool*` in `mjData` (`mjdata.h:182`).
Platforms other than Linux x86_64 raise an unsupported-platform error. License is Apache-2.0.

## Layout

| File | Role |
|---|---|
| `src/binding.jl` | The only file that calls libmujoco: 21 C entry points resolved with `dlsym`, `MjModel`/`MjData` owners with finalizers, field access at byte offsets from the 3.11.0 headers (`offsetof`). `api()` refuses any `mj_version()` other than 3011000. |
| `src/proximity_scene.jl` | `ProximityScene`, state types and the runner. |
| `test/binding_tests.jl` | Checks every offset against an `mj_*` readback or an MJCF-fixed value. |
| `test/scene_tests.jl` | Free-body agreement with ordinary spacecraft, determinism, state round trip, validation, a_fb, leak loop. |

## Runner

Public names: `ProximityScene`, `SceneBodyState`, `SceneState`, `body_state_from_initial_condition`,
`scene_reset!`, `scene_step!`, `scene_state`, `scene_set_state!`, `scene_time`, `scene_chief`,
`scene_body_names`, `scene_body_state`, `scene_ctrl`. Every scene owns its `mjModel` copy and `mjData` and
writes no global state; the step is deterministic and `scene_reset!(scene[, states])`, `scene_state`,
`scene_set_state!` and `copy(scene)` are the hooks a vectorized (RL) runner builds on.

One step, with `t = t0 + n dt` from an integer counter: `mj_step1` (fresh kinematics) -> per-body wrench
`m_i dg_i + F_ext_i - m_i a_fb` from those kinematics, at the COM in world (ECI-parallel) axes, into
`xfrc_applied` -> `mj_step2` -> chief RK4 under SpaceAGORA's gravity plus `a_fb` -> `n += 1`. The chief is the
scene's mass-weighted COM at reset; `opt.gravity` is zero; the integrator is `implicitfast` by default and
`dt` is required. Wrenches are held for one `dt`, so the scene is first order in `dt` (measured below).

SpaceAGORA internals used: see the Stage 1 report; all are reached through qualified names
(`SimulationModel.*`, `SimulationEngine.sample_planet_frame_with_lpi`).

## Measured agreement

Two free boxes (20 kg and 35 kg, about 50 m apart, 51.6 degree orbit at r = 7000 km) against the same two
bodies as ordinary spacecraft under `run_simulation` (dp8, reltol 1e-12), one orbit (5850 s). Maximum
position difference of either body: point mass 1.36e-2 m at dt = 0.1 s, 6.8e-3 at 0.05, 3.4e-3 at 0.025,
1.70e-3 at 0.0125; J2 identical to three digits. Halving ratios 2.00, 2.00, 2.00. Velocity difference at
dt = 0.0125 s: 2.2e-6 m/s; relative (target minus chaser) position difference 2.7e-3 m. The test tolerances
are set from these numbers.
