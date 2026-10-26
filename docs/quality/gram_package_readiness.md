# GRAM package loading readiness

The supported model route loads `SpaceAGORA` as a Julia package, imports
`GRAMSuite`, and uses the `SpaceAGORAGRAMSuiteExt` extension with
`SpaceAGORA.SimulationModel` types; the extension owns the integration methods.
`test/smoke/gram_package_readiness.jl` checks that route against real, already
prepared package sources, in fresh Julia processes. It is opt-in: no CI job and
no test suite runs it. It complements the synthetic
[loading-policy tests](gram_loading_contract.md) and changes no loader. The
runtime benchmark, which includes its own copies of the model types, is outside
what it checks.

## Running it

You need a Julia version that SpaceAGORA's `Project.toml` allows (CI uses 1.12)
and a prepared `GRAMSuite.jl` wrapper whose source revision you have checked
yourself. [GRAMSuite setup](../src/user/gramsuite_setup.md) describes
`scripts/prepare_gram_loading_environment.py`, which builds such a directory
offline. From the repository root:

```sh
JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_PKG_OFFLINE=true \
  julia --startup-file=no --compiled-modules=existing \
  test/smoke/gram_package_readiness.jl \
  --gram-project /absolute/path/to/prepared/GRAMSuite.jl \
  --output /absolute/path/to/new/readiness-report
```

Options:

- `--gram-project` (required): the prepared wrapper.
- `--output` (required): a directory that must not exist yet, so a report from
  an earlier run cannot satisfy a new check.
- `--space-project`: another prepared SpaceAGORA checkout (default: this one).
- `--timeout-seconds`: the limit for each child process (default 300).

Before starting any case, the script refuses a Julia version outside
SpaceAGORA's or the wrapper's `julia` compat, and a wrapper version outside
SpaceAGORA's `GRAMSuite` compat. It uses the existing depot for the positive
cases and never activates, instantiates, resolves or installs packages. Each
child runs with `JULIA_PKG_OFFLINE=true` and `JULIA_PKG_PRECOMPILE_AUTO=0`.
Existing compiled images may be used, but no new ones are built, because
SpaceAGORA's precompile workload runs a short simulation.

The SpaceAGORA environment you select must not already make GRAM discoverable;
the fallback and absent-package cases depend on that. If it does, the run reports
unmet acceptance, not a skip. Prepare a separate environment rather than altering
a working installation to satisfy the check.

## Cases

Each case runs in a new Julia process with an explicit `JULIA_LOAD_PATH`.

| Case | Required observation |
| --- | --- |
| `discoverable_space_first` | Both packages resolve to the requested sources, and the extension attaches to the canonical model types. |
| `discoverable_gram_first` | Importing in the opposite order reaches the same package, type and extension identities. |
| `fallback_space_first` | GRAM is initially undiscoverable; prepending the prepared wrapper to the load path enables the integration. |
| `fallback_gram_first` | The same fallback works with the opposite import order. |
| `absent_gram` | SpaceAGORA loads, the GRAM import fails because GRAM is absent, and the extension stays absent. |
| `absent_dependency` | The wrapper is discoverable in an empty depot, and its import fails because its StaticArrays dependency is absent. |

In the discoverable cases the selected SpaceAGORA environment stays ahead of the
wrapper. The fallback cases exercise the prepend policy that the examples and the
parallelization benchmark use; the runtime benchmark's append policy is not
exercised. The script undoes its own load-path change. It does not test any
production loader's project restoration, retry or installation behavior.

The positive cases record package UUIDs, source paths and hashes, versions, the
parent modules of the model types, extension attachment, constructor method
ownership, the shared lock and ephemeris hooks, repeated imports and preserved
caller aliases. One extension function is called, a counter probe that reads two
atomic counters. The script records what a direct call from the older frame does
and requires the call through `invokelatest` to return zero; no model is
constructed. This does not show world-age correctness for native
constructors or for full benchmark call chains.

Each report also records the active project, load path, depot path and working
directory. It checks that the environment variables are unchanged without
writing them to the report; the only values it records are `JULIA_LOAD_PATH`,
`JULIA_DEPOT_PATH`, `JULIA_PKG_OFFLINE` and `JULIA_PKG_PRECOMPILE_AUTO`. StaticArrays' path and version show which prepared dependency was
picked up. A Git revision is reported only when it belongs to the package
directory itself; for a source-only snapshot, the source hashes identify the
files. A folder name or a path you pass in is not proof of a revision.

## Reading the result

The run passes only when every case meets its own expectation; `summary.toml` in
the output directory records each case. The two `absent_*` cases confirm failure
boundaries and do not count as successful loading. Unexpected import errors,
wrong identities, missing reports, failed children and timeouts are failures or
unmet acceptance, never skips.

An attached extension is not enough on its own. The positive cases also require
canonical type and constructor-method ownership and unchanged repeated imports,
the GRAM native-wrapper references to stay unset, and the extension's native
construction counters to stay at zero. Loading the packages may load native
dependencies such as SPICE, but the check constructs no GRAM model, queries no
atmosphere, starts no worker and runs no simulation.

The report always marks native and worker readiness as unverified. It does not
qualify GRAM binaries or data, library architecture, kernel order, model
construction, copying or serialization, reset and release, scientific outputs,
installation on a clean machine, or the runtime benchmark's own types. A pass
covers the Julia version and the package sources it observed; repeat it for each
other combination you want to qualify.

A hosted job that runs this check would need its own independently prepared,
pinned wrapper inputs and would have to keep its observations. Neither the
existing hosted checks nor a local run count as that job.
