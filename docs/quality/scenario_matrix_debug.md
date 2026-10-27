# Scenario matrix debug utilities

Two diagnostics summarize position RMSE from the comparison runners in
`test/gmat_scenario_matrix.jl`, and are maintained together with that matrix.
Run normally, they propagate the full matrix, which takes a long time. They are
not installation checks, and they do not apply the matrix testsets' acceptance
limits.

## Invocation and output

From the repository root, with the examples environment set up
(`julia --startup-file=no examples/setup.jl`):

```sh
julia --project=examples scripts/dev/diagnostics/j2_parity_debug.jl
julia --project=examples scripts/dev/diagnostics/tb_matrix_debug.jl
```

`j2_parity_debug.jl` prints the per-axis RMSE in kilometres for every J2
scenario, the combined-axis norm, and the five worst cases. `tb_matrix_debug.jl`
prints the combined-axis RMSE in metres for both comparison targets over the full
body, gravity and third-body grid, then the missing scenario/target pairs and any
pair whose third-body-off and third-body-on values are exactly equal. It writes
`output/diagnostics/tb_matrix_rmse_m.csv` with the header
`planet,gravity,target,third_body,rmse_m`, full-precision metre values, and an
empty value where a scenario lacks one of its three axes.

Run as scripts, both files call their function. Included from other code, they
only define it; call `run_j2_parity_debug()` or `run_tb_matrix_debug()`
yourself. `run_tb_matrix_debug(; csv_path=…)` writes the CSV elsewhere. Both
functions also take the runner and input-check functions as keyword arguments,
which the probe below uses.

The matrix reads these environment variables:

- `SPACEAGORA_GMAT_SCENARIOS`: comma-separated scenario names, for example
  `earth_j2_tbfalse,earth_j2_tbtrue`. The name is historical; it selects
  scenarios for every target. The TB report still lists the whole grid and
  reports unselected cases as missing.
- `SPACEAGORA_GMAT_PARITY_SOLVER`: switches to the tighter parity solver
  settings.
- `SPACEAGORA_DEBUG_COMPARE_J2`, `SPACEAGORA_TELEMETRY_J2_SOURCE_DEFAULT`
  (`file_c20`, the default, or `planet_j2`) and
  `SPACEAGORA_TELEMETRY_J2_SOURCE_PLANET_SCENARIOS`: the J2 source and comparison
  controls. The J2 utility prints their current values.

These utilities call the runners directly, so `SPACEAGORA_SKIP_GMAT_MATRIX`,
which only gates the matrix testsets, does not disable them.

## References and prerequisites

| Diagnostic label | Runner | Reference input | Target setting |
| --- | --- | --- | --- |
| J2 summary, and the TB `GMAT` target | `_run_basilisk_scenario_matrix_result_once` | `data/telemetry/Basilisk_Examples_Full/*.feather` | `reference_target=:basilisk` |
| TB `STK` target | `_run_stk_scenario_matrix_result_once` | `data/telemetry/stk_results/*.csv` | `reference_target=:stk` |

The TB utility's `GMAT` console and CSV label is kept for compatibility with
older reports. It names the first row of the table, whose references are
Basilisk runs (the asset registry `data/assets_manifest.toml` lists them as
`basilisk_parity_references`), not GMAT output. The separate CYGNSS GMAT
comparison uses `data/telemetry/GMAT_Examples`, and neither utility runs it.

Neither reference set is in the public checkout. The Basilisk references live in
the lab's private repository `Space-FALCON-Lab/spaceagora-verification-data`;
lab members sync them into the gitignored `data/telemetry/Basilisk_Examples_Full/`
with `scripts/dev/fetch_private_telemetry.sh references`, which covers the
Basilisk and GMAT parity sets. Supply the STK CSVs in `data/telemetry/stk_results/`
separately.
See [`data/telemetry/PRIVATE_TELEMETRY.md`](../../data/telemetry/PRIVATE_TELEMETRY.md)
and [`docs/src/assets.md`](../src/assets.md) for access and handling. Do not
create substitute scientific inputs or commit fetched reference data.

Both utilities also need the four gravity files that the checkout tracks under
`data/planets/gravity/` (`EarthGGM05C.csv`, `Mars50c.csv`, `MGNP180U.csv`,
`LP165P.csv`) and the planetary kernel `de430.bsp` under
`data/GRAMSuite.jl/GRAM Suite 2.0/SPICE/spk/planets/`. That path is inside the
`GRAMSuite.jl` submodule, so the submodule's SPICE folder must be present. Other
kernels required by the configured planet models must be available too. Native
GRAM execution and CYGNSS flight data are not needed.

`ScenarioMatrixDebugSupport.require_matrix_inputs` checks these inputs before any
propagation; the TB utility checks the STK set as well, before starting either
comparison. A missing reference directory, an empty scenario list, a missing
selected reference file, a missing gravity file or a missing planetary kernel
raises an `ArgumentError` that names what is missing and, for the Basilisk set,
the fetch command.

## How the utilities load the matrix

`scripts/dev/diagnostics/scenario_matrix_debug_support.jl` includes
`test/gmat_scenario_matrix.jl` into the module `ScenarioMatrixDebugSupport` with
`SCENARIO_MATRIX_DEFINITIONS_ONLY = true`. Use the matrix helpers through that
module. `scripts/dev/diagnostics/tb_matrix_debug_defs_only.jl` is a compatibility
include path that loads the same module once; it holds no definitions of its own
and generates no source.

The matrix keeps all its testsets, plots and the example export inside
`_run_scenario_matrix_testsets()`. Including the matrix directly still runs them
with the usual skip flags, order and missing-input handling. With
`SCENARIO_MATRIX_DEFINITIONS_ONLY = true`, nothing runs, CYGNSS testsets and the
export included.

The bounded probe checks loading, reporting and CSV output with synthetic
summaries and stub runners, without a numerical campaign:

```sh
julia --project=. test/probes/scenario_matrix_debug_probes.jl
```

It does not establish trajectory equivalence, scientific acceptance or anything
about run time.

## Full-arc runs and export provenance

```sh
julia --startup-file=no --project=examples scripts/xval/xval_fullarc.jl <target> <outdir> [variant]
```

`<target>` is `gmat`, `stk` or `basilisk` (the last needs `XVAL_BASILISK_DIR`);
`[variant]` defaults to `committed`. The script compares every reference epoch in
each file. Its run directory, `<outdir>/<target>_<variant>`, must not exist yet.
`run_info.json` starts with status `running` and becomes `complete` only after all
selected cases finish and the inputs, source identity and environment controls
have been checked again. The record holds the source commit and tree, whether the
tree was dirty, the selected scenarios, solver settings, scenario models, input
SHA-256 hashes, the loaded kernels in order, and the manifest, series and results
checksums. Keep it with the run outputs. A failed or partial run cannot replace a
completed one in place.

`scripts/xval/xval_fullarc_export.py <xval outdir> <destination root>` validates
these records before writing anything. A primary export needs both the
`gmat_committed` and `stk_committed` runs, each with the exact 24-case matrix, the
same clean source revision, successful solves, matching output hashes and the
committed configuration. Renaming a sensitivity run's directory or labelling it
`committed` does not qualify it. Any non-empty `XVAL_*` or `SPACEAGORA_*`
variable other than `XVAL_SCENARIOS`, `SPACEAGORA_WARN_NORMALIZE` and
`SPACEAGORA_WARN_DEPRECATED_CONFIG` makes a run diagnostic. That covers frame
tables, GM, gravity files, planetary kernels, PCK overrides, relocated reference
directories and the parity solver mode, and any new control is refused for
primary export by default. A committed `gmat` or `stk` run refuses such overrides
before it creates any output. The primary CSV carries the configuration and
reference identities together with the run record's hash; keep the record to
recover the full settings.

The exporter refuses older runs that lack this metadata. Do not build a clean
record after the fact from a directory name; keep such runs as historical
evidence and reconcile their inputs separately. This export check changes no
trajectory acceptance threshold.

*Primary* here means a comparison under the committed configuration. The
settings that generated the reference trajectories are not verified: the older
GMAT generator selects JGM2 for Earth and has no Luna central-body case, while
this matrix uses EGM96 and LP165P, and the STK tide conventions were inferred from
reduced residuals because the record of the STK reruns has not been recovered.
Every run record and export therefore carries
`reference_provenance = "unverified_generation_settings"`. These comparisons,
including runs with substituted rotations, do not establish the simulator's
default accuracy or that every residual is a frame error.

### Diagnostic rotations and PCK overrides

Normal full-arc runs use the package's rotation unchanged. Only setting
`XVAL_FRAME_TABLE_DIR` before the script loads replaces it with the tabulated
rotation. Use a fresh process for each mode; changing the directory afterwards is
refused.

`SPACEAGORA_SPICE_PCK_OVERRIDES` is an ordered, comma-separated list of kernels;
relative paths resolve under the SPICE directory. Planet constructors load them
after all standard kernels and load them again if another body's construction
adds kernels, and repeated cached construction does not grow the kernel table.
The list is fixed at the first planet construction, an empty list included.
Changing or removing it needs a fresh process, or `SPICE.kclear()` followed by
`SpaceAGORA.SimulationModel.Planets._reset_furnished_kernels!()`; otherwise the
next construction throws `ArgumentError`. Do not edit kernel files or change the
kernel pool by hand during a run.

Later binary PCKs take precedence over earlier ones, and binary orientation data
take precedence over text PCK data whatever the loading order; see the
[NAIF PCK documentation](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/pck.html).

The scenario-utility probe, the run-record probe (`test/probes/xval_provenance_probes.jl`)
and the synthetic CSPICE precedence probe (`test/probes/pck_override_probes.jl`)
run in the probe-driver suite, `test/suites/09_probe_drivers.jl`. The exporter's
refusal tests, `test/scripts/test_xval_fullarc_export.py`, run in the `plan` job of
`.github/workflows/julia-ci.yml`; none of these needs private telemetry or a
propagation campaign.
