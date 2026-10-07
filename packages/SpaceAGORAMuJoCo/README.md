# SpaceAGORAMuJoCo

Optional MuJoCo proximity scenes for SpaceAGORA. Stage 1 is the standalone scene runner; Stage 2
runs a scene inside `run_simulation` (shadow spacecraft entries, a sync callback, guidance and navigation
reading the scene state; actuation and sensors are Stage 3):

```julia
config = SM.SimConfig._with_configuration(config;
    external_propagators = (ProximitySceneDynamics(scene, 1 => "chaser", 2 => "target"),))
SpaceAGORA.run_simulation(config)   # spacecraft 1 and 2 are now carried by the scene
```
 The root `SpaceAGORA` project
does not depend on this package, and a plain `using SpaceAGORA` never touches MuJoCo.

```julia
using SpaceAGORA, SpaceAGORAMuJoCo
SM = SpaceAGORA.SimulationModel
scene = ProximityScene(; mjcf_xml = xml, dt = 0.02,                 # dt is required, no default
    planet = SM.make_no_gram_planet(:earth), gravity_effectors = (SM.InverseSquaredGravityModel(),),
    initial_states = [SceneBodyState("chaser", r1, v1), SceneBodyState("target", r2, v2)])
scene_step!(scene)
scene_body_state(scene, "chaser")     # absolute inertial (r, v, q, ω)
```

## Installation

Package sources do not propagate from a dependency's `Project.toml`, so use the checkout:

```sh
julia --project=packages/SpaceAGORAMuJoCo -e 'using Pkg; Pkg.instantiate()'   # [sources] points at ../..
julia --project=packages/SpaceAGORAMuJoCo packages/SpaceAGORAMuJoCo/test/runtests.jl
```

or add `SpaceAGORA` and `SpaceAGORAMuJoCo` as path packages to your own environment
(`julia packages/SpaceAGORAMuJoCo/scripts/setup_env.jl <dir>` builds one from the root manifest, as CI does). The first load downloads the MuJoCo artifact; the download is verified against the sha256 in
`Artifacts.toml` and the unpacked tree against its git tree hash.

## Dependency review

| | |
|---|---|
| Library | MuJoCo 3.11.0, official Google DeepMind release `mujoco-3.11.0-linux-x86_64.tar.gz` |
| License | Apache-2.0 (`LICENSE` and `THIRD_PARTY_NOTICES` ship in the tarball) |
| sha256 | `681a9353b296a06682e9e757d6d21ab7a0a5a5bdb218a1c05cb9a38536817c79` (equals the published `.sha256` file and the GitHub asset digest) |
| Size | 31.7 MB download, 85 MB unpacked, `libmujoco.so.3.11.0` 6.1 MB |
| Links against | libc, libm, libdl, libpthread, librt only |
| Platforms | Linux x86_64. Other platforms raise an unsupported-platform error: the Windows zip and macOS dmg were not verified or tested, and `MuJoCo_jll` stops at 3.1.6 (2024) |
| Why not `MuJoCo_jll` | General has only 2.3.7 and 3.1.6; the Yggdrasil recipe is pinned to 3.1.6 |

MuJoCo runs in-process: a native fault ends the Julia process. `mjtBool` replaces `mjtByte` for
`eq_active` in 3.11.0 (still in `mjData`).
