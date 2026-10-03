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

## Full-arc runs and export provenance

`scripts/xval_fullarc.jl <target> <outdir> [variant]` compares every retained
reference epoch. Each run must use a new output directory. `run_info.json`
starts as `running` and becomes `complete` only after all selected cases finish
and their inputs, source identity and environment controls have been checked.
It records the source commit and tree, dirty state, selected scenarios, effective
solver settings, scenario models, input SHA-256 identities, loaded kernels in
order, and manifest, series and results checksums. Keep this record with the
run outputs. A failed or partial run cannot replace a completed run in place.

The Python exporter validates those records before writing any deliverable.
Primary export requires both the exact 24-case GMAT and STK matrices, the same
clean source revision, successful solves, matching output hashes and committed
configuration. Renaming a sensitivity directory or labelling it `committed`
does not qualify it. Any nonempty scientific `XVAL_*` or `SPACEAGORA_*` control
makes the run diagnostic, including frame tables, GM, gravity files, planetary
kernels, PCK overrides, reference-directory relocation and parity solver mode.
Selection and the two normalization/deprecation warning controls are recorded
but permitted. Committed runs refuse scientific overrides before creating output.
New controls are rejected for primary export until reviewed.
The primary CSV preserves configuration and reference identities with the run
record hash; keep the record to recover the full effective settings.

Historical runs without this metadata are refused by the new exporter. Do not
construct a retrospective clean record from a directory name. Preserve those
runs as historical evidence and reconcile their actual inputs separately.
No existing trajectory acceptance threshold is changed by this gate.

Here, *primary* means a comparison under the committed configuration. The
reference trajectories' generation settings remain unverified: the available
older GMAT generator selects JGM2 for Earth and has no Luna central-body case,
while this matrix selects EGM96 and LP165P. STK tide conventions were inferred
from reduced residuals, and the cited rerun record has not been recovered.
Every export therefore carries `unverified_generation_settings`. Neither these
comparisons nor the reported millimetre residuals with substituted rotations
establish default simulator accuracy or that all residuals are frame errors.

### Diagnostic rotations and PCK policy

Normal full-arc runs use the package's rotation method unchanged. Only an
explicit `XVAL_FRAME_TABLE_DIR` at script load enables the diagnostic replacement.
Use a fresh process for each mode; changing the directory afterwards is refused.

`SPACEAGORA_SPICE_PCK_OVERRIDES` is an ordered, comma-separated list, with relative
paths resolved under the SPICE directory. Constructors load these after all
standard kernels and reapply them if another body's construction loads more
kernels. Repeated cached construction does not grow the kernel table. The
ordered path list is fixed at first construction, including an empty list.
Changing or removing it requires a fresh process, or `SPICE.kclear()` followed
by `SpaceAGORA.SimulationModel.Planets._reset_furnished_kernels!()`. Do not modify
kernel files or manually mutate the pool during a run.

Later binary PCKs take precedence over earlier binary PCKs; binary orientation
data take precedence over text PCK orientation regardless of loading order.
See the [NAIF PCK documentation](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/pck.html).

The scenario utility, run-record and synthetic CSPICE precedence probes run in
the normal probe-driver suite. The export rejection tests run in the CI planning
job and require no private telemetry or propagation campaign.
