# Real-package GRAM loading readiness

The supported model route loads `SpaceAGORA` as a Julia package, imports
`GRAMSuite`, and uses the canonical `SpaceAGORAGRAMSuiteExt` extension with
`SpaceAGORA.SimulationModel` types. The extension owns the integration methods.
The runtime benchmark's separately included types and local adapters are outside
this package-loading check.

`test/smoke/gram_package_readiness.jl` is an opt-in diagnostic for already prepared
package sources and dependencies. It complements the synthetic
[loading-policy contracts](gram_loading_contract.md). It does not change an
existing loader or consolidate the four consumer policies.

## Run with a prepared environment

Use a Julia version allowed by SpaceAGORA's `Project.toml` and a prepared wrapper
whose source revision you have independently verified. The current CI runtime
line is Julia 1.12; one successful invocation qualifies only the exact runtime
and observed package sources in that report. Repeat the check for any additional
runtime or wrapper combination being qualified.

From the repository root:

```sh
JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_PKG_OFFLINE=true \
  julia --startup-file=no --compiled-modules=existing \
  test/smoke/gram_package_readiness.jl \
  --gram-project /absolute/path/to/prepared/GRAMSuite.jl \
  --output /absolute/path/to/new/readiness-report
```

`--space-project` can select another prepared SpaceAGORA checkout.
`--timeout-seconds` bounds each child process. The output directory must be new,
so reports from an earlier invocation cannot satisfy a new check. The diagnostic
uses the existing depot for positive cases and does not activate, instantiate,
resolve or install packages. Existing compiled images may be used; new compiled
images are disabled because SpaceAGORA's precompile workload includes a short
simulation.

The selected SpaceAGORA environment must initially leave GRAM undiscoverable for
the fallback and absent-package cases. An environment that already exposes GRAM
cannot demonstrate those boundaries; that condition is unmet acceptance, not a
passing skip. Prepare the environment deliberately rather than modifying a
working installation to satisfy this test.

## Case matrix

Each case runs in a new Julia process with an explicit package search path.

| Case | Required observation |
| --- | --- |
| Discoverable wrapper, SpaceAGORA first | Both packages resolve to the requested sources; the extension attaches to the canonical model types. |
| Discoverable wrapper, GRAM first | The opposite import order reaches the same package, type and extension identities. |
| Explicit fallback, SpaceAGORA first | GRAM is initially undiscoverable; explicitly prepending the prepared wrapper enables the intended package integration. |
| Explicit fallback, GRAM first | The same explicit fallback works with the opposite import order. |
| GRAM absent | SpaceAGORA loads, GRAM import fails specifically because GRAM is absent, and the extension remains absent. |
| Wrapper dependency absent | The real wrapper is discoverable in an isolated empty depot, but import fails specifically because its StaticArrays dependency is absent. |

The discoverable cases keep the selected SpaceAGORA environment ahead of the
wrapper. The explicit fallback cases characterize the prepend policy used by
examples and the package benchmark. The diagnostic restores its own deliberate
search-path change. This is not proof of any production loader's project
restoration, retry or installation behavior; those policies remain separate.
The runtime benchmark's append order is not replaced by this diagnostic.

Positive cases inspect package UUIDs, source paths and hashes, versions, model
parent modules, extension attachment, constructor method ownership, shared-lock
and ephemeris hooks, repeated imports and preserved caller aliases. A harmless
extension counter probe records direct older-frame behavior and checks the
`invokelatest` path without constructing a model. It does not establish world-age
correctness for native constructors or complete benchmark call chains.

Project, search-path, depot and working-directory observations are retained.
Environment-variable preservation is checked without writing variable values to
the report. StaticArrays' actual path and version identify which prepared
dependency was selected. A Git revision is reported only when it belongs to the
package directory itself. Source hashes identify source-only snapshots; a folder
name or supplied path is not independent proof of a revision pin.

## Interpret the result

All cases must pass their distinct expectations for overall diagnostic success.
Negative cases validate failure boundaries; they do not count as successful
package loading. Unexpected import errors, incorrect identities, missing reports,
child failures and timeouts are failures or unmet acceptance, never passing skips.

An attached extension alone is insufficient: the positive cases also require
canonical type and constructor-method ownership and unchanged repeated imports.
The GRAM native-wrapper references must remain unused, and the extension's native
construction accounting must remain zero. Ordinary package initialization may
load native dependencies such as SPICE; this check promises no GRAM native model
construction, native atmosphere queries, worker startup or simulation.

The reports always leave native and worker readiness unverified. They do not
qualify GRAM binaries/data, library architecture, kernel order, model construction,
copy/serialization, reset/release behavior, scientific outputs, performance,
clean-machine installation or the runtime benchmark's independent types.

This diagnostic is deliberately outside the default no-GRAM test suite. A future
hosted job must supply independently prepared pinned wrapper inputs and retain
its actual observations. Existing hosted checks or a successful local invocation
must not be relabelled as that job having run.
