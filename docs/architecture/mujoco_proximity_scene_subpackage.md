# MuJoCo proximity scenes: external companion

MuJoCo proximity scenes belong to the independent
[SpaceAGORAMuJoCo.jl repository](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl);
core pins its source under `packages/robotics/`. This page remains so that
older links to the architecture note still resolve.

The companion owns the native binding and its artifact declaration, the scene
runner, the `ProximitySceneDynamics` adapter, the scene and engine tests, and
the [Basilisk comparison harness](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl/tree/main/benchmarks/basilisk_mujoco_xval).
Its [architecture reference](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl/blob/main/docs/architecture.md)
covers the frame conventions, the numerical method, the measured agreement, the
refusals and the known limits. The coupled comparison study lives in the private
SpaceAGORAParallelizationStudies repository, which records the core revision it
pins and what is needed to reproduce it.

Core keeps the generic external-propagation lifecycle and the shared dynamics.
Adapters connect through the supported
[`SpaceAGORA.ExternalPropagationServices` facade](../src/reference/external_propagation_services.md)
and follow the [common adapter standard](../src/reference/adapter_contract.md);
they do not reach into private engine modules or runtime storage layouts. See
[installation and environment](../src/user/installation_environment.md) for
setting up the companion. Loading SpaceAGORA alone neither loads the companion
nor downloads MuJoCo.
