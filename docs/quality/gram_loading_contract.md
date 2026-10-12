# GRAM loading policies and shared discovery

GRAM consumers currently use different loading policies. The native-free fixture
`test/unit/environment/gram_loading_contract_tests.jl` executes selected production
loader definitions against a synthetic environment to keep these differences
observable while adapter ownership and migration remain open.

`RuntimeServices._prepare_gramsuite_load_path!(vendor; prepend)` shares only
package discovery and conditional fallback-path insertion. The caller explicitly
selects prepend or append; the helper preserves short-circuit discovery and
existing duplicate entries. Every import remains in its original caller module
and at its original point in startup. Project activation, retries, aliases,
worker transport, native resources and scientific settings remain caller-owned.

| Consumer | Path and failure policy |
| --- | --- |
| `examples/common.jl` | Use a discoverable package first; otherwise prepend an existing vendored project. Wrap import errors. For recognized missing-dependency errors, instantiate the vendored project once, restore the previous project (or repository project), and retry once. `setup_gram_example!` installs aliases in its requested module without replacing existing bindings. |
| `benchmarks/studies/performance_runtime_analysis/main.jl` | Append an existing vendored project when discovery fails. Import unconditionally at entrypoint load and propagate failure without loader-level installation or retry. The entrypoint raw-includes its own model modules and owns separate adapter methods. |
| `benchmarks/studies/parallelization_performance/cases.jl` | Use the package model modules; prepend an existing vendored project when discovery fails. Propagate import errors without installation or retry. The entrypoint invokes this loader only for selected GRAM cases. |
| Core process workers | Prepend an existing vendored project inside the existing worker callback. Import separately into worker `Main`; warn on failure and continue to SPICE setup for non-GRAM campaigns. Worker project precedence, preload and warmup ordering are unchanged. |

The runtime benchmark's append order preserves the caller project ahead of the
vendored environment. `benchmarks/studies/parallelization_performance/HANDOFF.md`
records that prepending caused fresh workers to resolve dependencies against the
vendored project. The examples and package benchmark currently prepend instead.
These policies must not be unified by a blind replacement.

The examples and package benchmark guard on `isdefined(SpaceAGORA, :GRAMSuite)`,
while their import binds `GRAMSuite` in the including caller module. The fixture
preserves that distinction: a successful caller import does not satisfy the
package-binding guard, and an existing package binding skips loading even if the
caller has no binding. This characterizes present behavior; it does not endorse
the guard as an extension-readiness check. Repeated calls are also exercised with
discovery still absent after synthetic import: the current loaders insert the
vendored path again rather than deduplicating it.

## What the fixture proves

The fixture parses the seven selected discovery/loader/setup definitions from their current
source files. It preserves their branches, exception handling, path operations,
retry sequence, `try/finally` restoration and alias installation. Only package
discovery, directory presence, package import, package-manager acquisition and
active-project observation are redirected. A private load path and fake package
manager record the resulting operation order. Both recognized dependency-error
messages, unrelated errors, failed retry, failed instantiation and no prior
project are covered. Sentinel projects make prepend versus append observable.

The synthetic import creates or preserves a binding in the actual fixture caller
module and can throw supplied failures. It does not run Julia's package loader or
extension mechanism. Unexpected macros and calls fail fixture extraction rather
than adding new effects implicitly. The shared helper itself is executed in the
synthetic runtime-services module, rather than replaced with a test implementation.
A separate `test/unit/parallel/process_gram_loading_tests.jl` fixture executes
the production bootstrap and shared helper with synthetic transport, import and
kernel boundaries, including import failure and warning-before-kernel continuation.

## Supported direction and next acceptance boundary

The supported model route is package-loaded `SpaceAGORA` with imported
`GRAMSuite` and the canonical `SpaceAGORAGRAMSuiteExt` extension. Use
`SpaceAGORA.SimulationModel` types. A successful import, an attached extension
and usable methods are separate facts; none alone qualifies native resources.
The caller's prepared environment keeps precedence. A vendored wrapper remains
an explicit fallback when discovery fails.

The shared discovery mechanism does not change installation/retry policy, search
order, the current package-binding guard, worker setup, SPICE ordering or the
runtime benchmark's independent types.
The real-package readiness diagnostic described in
[GRAM package readiness](gram_package_readiness.md) adds a separate fresh-process
acceptance boundary. The synthetic fixture above remains the loading-policy
regression test. Native and worker qualification still require their own checks.

## Limits of the synthetic fixtures

Successful synthetic import does not prove package resolution, precompilation,
extension attachment, supported dependency pins, native availability, world-age
behavior or worker startup. The fixture does not execute complete entrypoints,
case selection, distributed bootstrapping, simulations or native libraries. The
real active project, load path, command arguments and package-load state remain
unchanged by it.

Core process workers have another deliberate policy: a failed GRAM import warns
and continues before kernel setup, so non-GRAM campaigns can proceed. The runtime
benchmark also prepares worker projects separately. The bootstrap fixture makes this warning/continuation order observable; it does
not replace actual fresh-worker qualification.

The canonical `SpaceAGORAGRAMSuiteExt` owns package-type constructors, native
allocation and release accounting, construction recipes, copy/serialization,
locking and ephemeris/cache hooks. The legacy runtime adapter targets independent
types and does not provide the same lifecycle contract. Text-identical forwarding
methods do not establish equivalent ownership. Agree on a supported route with
the benchmark and native owners before changing these adapters or fallback
policies. This fixture closes no native, scientific or performance acceptance.
