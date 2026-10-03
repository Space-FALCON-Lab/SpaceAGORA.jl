# Scenario matrix debug utilities

These diagnostics summarize position RMSE from the maintained comparison runners
in `test/gmat_scenario_matrix.jl`. They are verification utilities, owned with that
matrix. They can run long propagations when invoked normally; they are not quick
installation checks and do not enforce the matrix testsets' acceptance limits.

## Invocation and output

From a configured repository environment:

```sh
julia --project=. scripts/j2_parity_debug.jl
julia --project=. scripts/tb_matrix_debug.jl
```

The J2 utility prints per-axis RMSE in kilometres, the combined-axis norm, and the
five worst J2 cases. The TB utility prints both comparison targets, missing pairs
and exactly equal third-body on/off pairs. It writes
`scripts/tb_matrix_rmse_m.csv` with the existing header, ordering, target labels,
full-precision metre values and empty values for incomplete axis triplets.
Including either launcher only defines its functions; call
`run_j2_parity_debug()` or `run_tb_matrix_debug()` explicitly after inclusion.
The TB function accepts `csv_path` to redirect its output.

`SPACEAGORA_GMAT_SCENARIOS` retains its historical name and accepts comma-separated
scenario names, for example `earth_j2_tbfalse,earth_j2_tbtrue`. The TB report still
lists the full body/gravity/third-body grid; unselected cases are reported as
missing. Existing solver and J2 controls remain matrix inputs, including
`SPACEAGORA_GMAT_PARITY_SOLVER`, `SPACEAGORA_DEBUG_COMPARE_J2`,
`SPACEAGORA_TELEMETRY_J2_SOURCE_DEFAULT`, and
`SPACEAGORA_TELEMETRY_J2_SOURCE_PLANET_SCENARIOS`. These debug invocations call the
runners explicitly, so `SPACEAGORA_SKIP_GMAT_MATRIX` does not disable them.

## Reference identities and prerequisites

| Diagnostic label | Maintained runner | Reference input | Existing scientific settings |
|---|---|---|---|
| J2 summary and legacy TB `GMAT` target | `_run_basilisk_scenario_matrix_result_once` | `data/telemetry/Basilisk_Examples_Full/*.feather` | `reference_target=:basilisk` |
| TB `STK` target | `_run_stk_scenario_matrix_result_once` | `data/telemetry/stk_results/*.csv` | `reference_target=:stk` |

The TB `GMAT` console/CSV label is retained for compatibility with historical
reports. It identifies the first row above; it does not establish GMAT provenance
for the Basilisk files. The asset registry identifies these as Basilisk runs.
The separate CYGNSS GMAT comparison uses `GMAT_Examples`, and neither debug
launcher runs that comparison. Reference paths, constants, scientific target
symbols and acceptance limits are unchanged by the definition-loading cleanup.

Both utilities need the configured Julia project, selected Basilisk reference
files, the four gravity CSVs (`EarthGGM05C.csv`, `Mars50c.csv`, `MGNP180U.csv`,
`LP165P.csv`), and the matrix's existing SPICE assets, including `de430.bsp` under
`data/GRAMSuite.jl/GRAM Suite 2.0/SPICE/spk/planets/`. Other kernels required by the
configured planet models must also be available. TB additionally needs the STK
CSVs for the same selected scenarios. It checks both sets before starting either
comparison. Missing directories, empty scenario discovery, missing selected
files, gravity CSVs or the planetary kernel produce an actionable input error.

See [`data/telemetry/PRIVATE_TELEMETRY.md`](../../data/telemetry/PRIVATE_TELEMETRY.md)
and [`docs/src/assets.md`](../src/assets.md) for asset access. Basilisk references
can be synced with `scripts/dev/fetch_private_telemetry.sh references`; provide
STK CSVs separately. Do not create substitute scientific inputs or commit fetched
reference data. Native GRAM execution and CYGNSS flight data are not prerequisites
for these orbit comparison diagnostics.

## Definition loading and bounded regression check

`scripts/scenario_matrix_debug_support.jl` loads the matrix in the isolated
`ScenarioMatrixDebugSupport` module with `SCENARIO_MATRIX_DEFINITIONS_ONLY=true`.
`scripts/tb_matrix_debug_defs_only.jl` remains an idempotent compatibility include
path for this module, not an independent source of matrix definitions. It no
longer creates a generated source file. Access matrix helpers through the module.

All matrix testsets, plots and the later example export are behind
`_run_scenario_matrix_testsets()`. Direct inclusion of the matrix retains its
existing automatic execution, skip flags, testset ordering and missing-input
behavior, including the Venus and CYGNSS entrypoints. Definitions-only import
suppresses the entire runner, including CYGNSS and the late export.

Run the bounded probe with:

```sh
julia --project=. test/probes/scenario_matrix_debug_probes.jl
```

The probe uses deterministic summaries and stub dispatch to check loading,
reporting and CSV output without a numerical campaign. The launcher runner and
input-check keyword arguments support those fixtures. These checks do not
establish trajectory equivalence, scientific acceptance or performance gains.

## Full-arc primary exports

`scripts/xval_fullarc.jl` records the effective planetary kernel and solver settings,
requested model/reference overrides, reference and gravity-file digests, and a
versioned completion record bound to its result, manifest, and series files.
The record includes loaded SPICE kernel digests in precedence order. Export
requires Python 3.11 or later, plus the existing pandas/pyarrow dependencies.
GMAT/STK `committed` runs reject nonempty `XVAL_*`, `SPACEAGORA_SPICE_*`,
`SPACEAGORA_TELEMETRY_*`, `SPACEAGORA_SOLVER_*`, and
`SPACEAGORA_GMAT_PARITY_SOLVER` overrides. `XVAL_SCENARIOS` may select a diagnostic subset;
primary export still requires all 24 unique cases for each reference target.
Use an explicit sensitivity variant name for frame, PCK, kernel, GM, field, or
reference-directory experiments. The exporter rejects legacy, interrupted,
modified, dirty-source, or sensitivity output as primary, and retains model/input
identity in the primary CSV. Legacy primary runs need to be rerun; a directory
name is not sufficient evidence. These records identify execution inputs and do
not independently validate the scientific provenance of the reference data.
