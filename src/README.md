# Source directory

The source is grouped into four folders by responsibility:

```text
src/
├── SpaceAGORA.jl
├── precompile_workload.jl
├── interfaces/
│   ├── cli/
│   ├── io/
│   └── analysis/
├── orchestration/
│   ├── simulation/
│   ├── mission/
│   └── parallel/
├── models/
│   ├── vehicle/
│   ├── environment/
│   ├── dynamics/
│   └── gnc/
└── foundation/
    ├── core/
    └── assets/
```

- `interfaces/` is how people and tools reach a run: the command-line tool
  (`cli/`), result and checkpoint files (`io/`), and analysis code such as the
  conformance checks, planner-comparison reports, telemetry verification and
  the viewer export (`analysis/`). Hooks that belong to one kind of model stay
  with that model, not here.
- `orchestration/` prepares and runs simulations (`simulation/`), defines
  missions and their segments (`mission/`), and decides how work is spread over
  threads and processes (`parallel/`).
- `models/` holds vehicle construction, the environment, the equations of
  motion and the GNC contracts and algorithms.
- `foundation/` holds shared abstract types, frame and numerical utilities, and
  asset lookup. Concrete run buffers and caches are in
  `orchestration/simulation/state/` instead.

The folders are not Julia modules and not strict dependency layers.
`SpaceAGORA.jl` is the package entry point, with `precompile_workload.jl`
beside it. `orchestration/simulation/composition/simulation_model.jl` builds the
`SimulationModel` module by including files from several of the folders, so a
file's folder does not tell you which module it belongs to; the `module` line
and the include that loads it do. Use the package API rather than including
source files from a script.

Inside `orchestration/simulation/`, `config/` holds scenario settings and
presets, `engine/` runs the simulation, and `callbacks/` uses the internal
operations that `callback_services.jl` declares and the engine implements.
Result file names come from `interfaces/io/paths/` and checkpoint files are
written by `interfaces/io/checkpoints/`, but the engine decides when to save and
when to resume.

Planning code is under `models/gnc/planning/`: RPO, robot-arm references, and
the geometry and tree operations both share. `models/gnc/companions/` holds
SpaceAGORA's side of the optional companions. The HYPR search and the EDG
numerical decisions live in the companion packages; SpaceAGORA keeps the
service modules they use and the simulation side, including the EDG model
types, lifecycle and actuation. Completed planner-comparison reports are
produced in `interfaces/analysis/planning/`.

Environment code supplies atmospheric, ephemeris, terrain and gravity-field
quantities. Dynamics code turns those quantities into forces and torques on a
spacecraft. Some modules span both folders; the READMEs in
`models/environment/` and `models/dynamics/effectors/` say which.

Outside `src/`:

- [`ext/`](../ext/) holds SpaceAGORA's Julia package extensions, which load
  when you also load PlotlyJS, PythonCall, SpaceAGORAEDG, SpaceAGORAVisualizer,
  GRAMSuite or Test.
- [`packages/`](../packages/README.md) pins the companion repositories. HYPR's
  SpaceAGORA adapter is an extension inside HYPR, active when both packages
  are loaded.
- [`examples/`](../examples/) holds runnable scenarios and
  [`scripts/`](../scripts/) operational tools; the browser viewer is
  [SpaceAGORAVisualizer.jl](https://github.com/Space-FALCON-Lab/SpaceAGORAVisualizer.jl).

[How a Simulation Fits Together](../docs/src/user/architecture.md) describes
the runtime flow, who owns what, and where new work belongs.
