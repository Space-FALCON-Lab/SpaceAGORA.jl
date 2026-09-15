---
id: analysis.example_support_run_and_report
label: run_and_report
kind: function
source:
  file: src/analysis/verification/telemetry_verification/example_support.jl
  symbol: run_and_report
  lines:
  - 177
  - 177
inputs:
- id: args
  type: SM.SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: result
  type: Any
  units: n/a
  description: Return value of `run_and_report`. Returns `saved_csv_path`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- analysis
origin: agent
---

# run_and_report

## Purpose
The common tail of every example script: run the simulation, report wall time, and hand back the CSV path if one was written.

## Design & Implementation
Applies `_example_smoke_args` first so smoke mode is honoured uniformly, times `run_simulation` with `@elapsed`, then looks for `simulation_results.csv` under the effective results directory. If results were enabled and the file exists it reads it with `CSV.read` into a `DataFrame` purely to print the row count and absolute path, and returns that path; otherwise it returns `nothing`. The computational time is printed in every case.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | SM.SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `run_and_report`. Returns `saved_csv_path`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/example_support.jl`

**Downstream**

- `callees` → [[analysis.example_support__example_smoke_args|_example_smoke_args]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:178-178`
- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:184-184`
- `callees` → [[simulation.run_simulation|run_simulation]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:179-179`
- `callees` → [[simx.engine_execution_run_simulation|run_simulation]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:179-179`
- `callees` → [[spaceagora.run_simulation|run_simulation]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:179-179`
<!-- vulcan:connections:end -->

## Limitations
Reading the whole CSV back just to count rows doubles the I/O for large runs; the filename is hard-coded, so a configuration using `generate_filenames` to produce a different name reports no saved CSV even though one was written.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/example_support.jl` line 177.
