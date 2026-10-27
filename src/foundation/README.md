# Shared foundation

Code that many parts of SpaceAGORA use and that belongs to none of them.

`core/` holds:

- `types/`: the abstract model types (`AbstractTypes`), the state and
  environment sample types effectors receive (`EffectorSampling`), and the
  older integer model codes (`LegacyModelCodes`);
- `numerics/`: geodesy, orbital geometry and quaternion utilities;
- `state/` and `interfaces/reference_system.jl`: reference-system
  configuration and frame operations;
- `interfaces/hook_deprecation.jl`: once-per-key warnings for deprecated names
  and pre-run configuration checks (`HookDeprecation`).

`assets/` finds and checks data files; see [its README](assets/README.md).

Some shared pieces live elsewhere: the composition of the `SimulationModel`
module in `../orchestration/simulation/composition/`, run buffers and other
concrete run state in `../orchestration/simulation/state/`, and the onboarding
configurations in `../orchestration/simulation/config/presets/`. These folders
are not package boundaries; the files are still loaded in one fixed include
order.

Checks that belong to one model stay with it: the atmosphere models validate
their own epochs, coordinates and payloads, geometry readers stay with the
vehicle structure code (asset loaders may call them), and command-line
formatting and argument handling stay in `../interfaces/cli/`.
