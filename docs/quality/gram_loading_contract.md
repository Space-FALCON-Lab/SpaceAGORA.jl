# GRAM loading policies

Four places in the repository load the `GRAMSuite` Julia wrapper, and they do it
differently on purpose. They share one small helper for package discovery and
keep everything else to themselves. The native-free test
`test/unit/environment/gram_loading_contract_tests.jl` runs their loader
definitions against a synthetic environment so that the differences stay visible
and cannot be merged by accident.

## The shared helper

`RuntimeServices._prepare_gramsuite_load_path!(vendor; prepend)`, in
`src/orchestration/simulation/runtime_services.jl`, does one thing: if
`Base.find_package("GRAMSuite")` finds nothing and the vendored directory
`vendor` exists, it adds that directory to `LOAD_PATH`, at the front when
`prepend=true` and at the end when `prepend=false`. The caller always chooses.
The helper does not remove duplicate entries. Each caller keeps its own
`import GRAMSuite`, in its own module and at its own point in startup, and owns
project activation, retries, aliases, worker transport, native resources and
scientific settings.

## The four consumers

| Consumer | Path and failure policy |
| --- | --- |
| `examples/common.jl` (`ensure_gramsuite_loaded!`, `setup_gram_example!`) | Use a discoverable package first; otherwise prepend the vendored `data/GRAMSuite.jl`. Wrap import errors in an `ErrorException`. For the two recognized missing-dependency messages, instantiate the vendored project once, reactivate the previous project (or `examples/` when there was none), and retry the import once. `setup_gram_example!` adds `InitialTime` and `GRAMAtmosphereModel` aliases to the requested module without replacing existing bindings. |
| `benchmarks/studies/performance_runtime_analysis/main.jl` (`_ensure_perf_gramsuite_loaded!`) | Append the vendored project when discovery fails. Import unconditionally when the entry point loads, and let a failure propagate; no installation or retry. This benchmark includes its own copy of the model modules from source and defines its own adapter methods. |
| `benchmarks/studies/parallelization_performance/cases.jl` (`ppc_ensure_gramsuite_loaded!`) | Use the package's model modules; prepend the vendored project when discovery fails. Let import errors propagate; no installation or retry. The file calls the loader only when the command line names a GRAM case. |
| Process workers (`_bootstrap_process_worker!` in `src/orchestration/parallel/process/worker_pool.jl`) | Prepend the vendored project on the worker and import `GRAMSuite` into the worker's `Main`. On failure, warn and continue to SPICE kernel setup, so campaigns that use no GRAM model still run. Worker project precedence, preloading and warm-up order are unchanged. |

The runtime benchmark appends so that the caller's project stays ahead of the
vendored environment; prepending there would change dependency resolution for
fresh workers. The examples and the parallelization benchmark prepend. Do not
replace one policy with another without checking what it changes.

The examples and the parallelization benchmark skip loading when
`isdefined(SpaceAGORA, :GRAMSuite)` is true, but their `import` binds `GRAMSuite`
in the calling module, not in `SpaceAGORA`. The test keeps that distinction: a
successful import in the caller does not satisfy the guard, and an existing
`SpaceAGORA.GRAMSuite` binding skips loading even when the caller has none. This
records present behavior; it is not an endorsement of the guard as a check that
the extension is ready. The test also calls the loaders repeatedly while
discovery still fails after the synthetic import, and the loaders add the vendored
path again each time.

## What the tests run

`test/unit/environment/gram_loading_contract_tests.jl` parses seven definitions
from their source files: `_prepare_gramsuite_load_path!`,
`_instantiate_vendored_gramsuite!`, `_missing_dependency_error`,
`ensure_gramsuite_loaded!`, `setup_gram_example!`,
`_ensure_perf_gramsuite_loaded!` and `ppc_ensure_gramsuite_loaded!`. Their
branches, exception handling, path operations, retry sequence, `try`/`finally`
restoration and alias installation run as written. Only package discovery,
directory checks, the package import, `Pkg` acquisition and active-project
lookups are redirected, to a private load path and a fake package manager that
record the order of operations. The cases cover both recognized dependency-error
messages, unrelated errors, a failed retry, a failed instantiation and the case
with no previous project. Sentinel projects make prepend and append
distinguishable. A call or macro the test does not expect makes the extraction
fail, so a new side effect cannot slip in unnoticed. The shared helper runs in a
synthetic `RuntimeServices` module rather than being replaced by a test version.

`test/unit/parallel/process_gram_loading_tests.jl` does the same for
`_bootstrap_process_worker!` and the shared helper, with synthetic transport,
import and kernel steps. It checks an import failure and that the warning comes
before kernel setup.

The synthetic import creates or keeps a binding in the test's caller module and
can throw a supplied error. It does not run Julia's package loader or extension
mechanism. A passing synthetic import therefore says nothing about package
resolution, precompilation, extension attachment, dependency pins, native
availability, world age or worker startup. The tests do not run entry points,
case selection, distributed startup, simulations or native libraries, and they
leave the real active project, load path, command-line arguments and loaded
packages unchanged.

## Supported route

The supported model route loads `SpaceAGORA` as a package, imports `GRAMSuite`,
and uses the `SpaceAGORAGRAMSuiteExt` extension with `SpaceAGORA.SimulationModel`
types. A successful import, an attached extension and working methods are
separate facts, and none of them alone shows that native resources are ready. The
caller's prepared environment keeps precedence; the vendored wrapper remains an
explicit fallback when discovery fails.

The opt-in diagnostic described in [GRAM package readiness](gram_package_readiness.md)
checks real package loading in fresh processes. The synthetic tests above remain
the regression tests for the loading policies. Native GRAM and workers each need
their own checks.

`SpaceAGORAGRAMSuiteExt` owns the package-type constructors, native allocation
and release accounting, construction recipes, copying and serialization, locking,
and the ephemeris and cache hooks. The runtime benchmark's adapter targets its
own copies of the types and does not provide that lifecycle, even where its
forwarding methods read the same. Agree on a supported route with the benchmark
and native owners before changing these adapters or fallback policies.
