# MuJoCo proximity scenes: external companion

MuJoCo proximity scenes now belong to the independent
[SpaceAGORAMuJoCo.jl repository](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl).
This page preserves the old architecture-note link after extraction.

The companion owns the native binding and artifact declaration, scene runner,
`ProximitySceneDynamics` adapter, scene and engine tests, and
[Basilisk comparison harness](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl/tree/main/benchmarks/basilisk_mujoco_xval).
Its [architecture reference](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl/blob/main/docs/architecture.md)
contains the Stage 1/2 frame conventions, numerical method, measured agreement,
refusals and remaining limits. Extraction does not change those scientific
claims or provide new benchmark results. The coupled timing study from PR253
remains in SpaceAGORA under `benchmarks/studies/coupled_local_comparison/`.

Core retains the generic external-propagation lifecycle and shared dynamics.
Adapters connect through the supported
[`SpaceAGORA.ExternalPropagationServices` facade](../src/maintainer/external_propagation_services.md)
and follow the [common adapter standard](../src/maintainer/adapter_contract.md);
they do not reach through private engine modules or runtime storage layouts.
See [installation and environment](../src/user/installation_environment.md)
for explicit companion setup. Loading SpaceAGORA alone does not load the
companion or download MuJoCo.
