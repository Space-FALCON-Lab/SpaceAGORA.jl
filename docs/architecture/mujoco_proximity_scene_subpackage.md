# MuJoCo proximity scenes: subpackage architecture (Stages 0 to 2)

Developer note for `packages/SpaceAGORAMuJoCo`. Design basis: the 2026-10-06 MuJoCo extension design memo
(coupling option c, chief frame, Encke gravity). This note records what Stage 1 built and the Stage 0
dependency decision.

## Packaging

The scene runner is a separate package in `packages/SpaceAGORAMuJoCo`, following `packages/SpaceAGORAHYPR`:
its own `Project.toml` (depends on `SpaceAGORA`, compat `0.2`), its own lazy `Artifacts.toml`, `src/` and
`test/`. The root `Project.toml` and `Artifacts.toml` are untouched, and nothing was added under `src/` or
`ext/`, so no ownership boundary moves and no contract gate scans the new code (the gates that list
`packages/SpaceAGORAHYPR/src` as a root do not include it). Putting `ext/` and this package under the
topology contract was done in Stage 2, together with the core hooks that integration needs (see "Stage 2" below).

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

## Stage 2: scenes inside `run_simulation`

Core knows nothing about MuJoCo. It has a generic hook for spacecraft that an external integrator owns for a
whole run (`SimulationModel.ExternalPropagation`, internal names, documented in its module docstring and in the
topology contract item 11): `external_spacecraft`, `external_step`, `external_preflight`, `external_prepare`,
`external_sync!`, `external_time`, `external_state`, `external_acceleration`. The owner is configured through
`SimulationConfiguration.external_propagators` (default `()`; runs without it are bit-identical to before).

`ProximitySceneDynamics(scene, 1 => "chaser", 2 => "target")` is the subpackage's owner. The scene is a template:
each run copies the native model, starts from the initial state of the paired spacecraft (`u0`), and never steps
the template, so one template can serve every sample of a threaded Monte Carlo campaign. `deepcopy` of a scene
copies the native model.

Engine behavior. Owned spacecraft stay in `u.sc` as shadow entries. Between syncs their right-hand side is the
chief acceleration (SpaceAGORA's gravity at the chief, displaced by the chief's velocity to the current time) and,
with `orientation_sim`, plain attitude kinematics at the last body rate. After every accepted step a
`DiscreteCallback`, installed before every other discrete callback, advances the scene by whole `dt` while
`n dt <= t` (an integer counter; the solver is never forced to take `dt`-sized steps) and overwrites the shadow
entries. The scene lags `t` by `tau < dt`, so the written state is the scene state advanced by `tau` with the same
model the shadow right-hand side uses (midpoint chief acceleration, constant body rate); at a scene-step boundary
the scene state is written unchanged. Guidance, navigation and control callbacks therefore read the scene state at
their ticks exactly, and saved fields on a tick are bit-identical to the standalone runner.

Attitude convention. SpaceAGORA stores `q = (x, y, z, w)` as the active body-to-inertial rotation (`rot(q)` is
inertial to body, kinematics `qdot = 1/2 q (x) (omega_body, 0)`); MuJoCo's `xquat` is `(w, x, y, z)`, Hamilton,
body to world, with the free joint's rate in the body frame and world axes parallel to ECI. The mapping is a
reordering of components (`sa_to_mujoco_quaternion`, `mujoco_to_sa_quaternion`) and the rate needs no
transformation. `scene_body_state` now returns the body-frame rate `omega` (and `omega_world`).

Refusals (clear `ArgumentError`s at setup): checkpointing and resume; `run_constellation_ensemble`; solver modes
other than `:tsit5`, `:auto_stiff`, `:rodas5p`, `:dp8`; the flat RHS route; articulated joints, compliant
attachments or a robot arm on an owned spacecraft; a control effector on an owned spacecraft (actuation is Stage 3;
reading state for guidance or navigation is allowed); guidance, navigation or control rates that are not
integer multiples of the scene `dt`; orbit-count termination and atmosphere-interface events (continuous events on
every spacecraft). Impact detection is not applied to shadow entries. Process-pool campaigns are not supported
(native state is not serializable).

Known limits. Owned bodies feel gravity only (thrusters, wheels, joint motors and environment loads are Stages 3
and 4). The adaptive solver's first step and error norm see the whole state vector, so an ordinary spacecraft's
results change at round-off level (about 1e-8 m) when other spacecraft, owned or not, share the run.
Measured with `mj_step1`/`mj_step2`, `integrator=:rk4` gave the same first-order attitude errors as `:implicitfast` (2.11 rad and 2.22 rad at dt = 0.04 on an intermediate-axis tumbling case, 0.07 rad with implicitfast on a stable spin), so it is not a higher-order option here; this is consistent with the split step not running RK4 stages, which was not verified in the MuJoCo source.
