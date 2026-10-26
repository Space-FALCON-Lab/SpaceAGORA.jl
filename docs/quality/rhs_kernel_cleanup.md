# Shared per-spacecraft dynamics kernels

`src/orchestration/simulation/engine/dynamics_rhs.jl` has one private
per-spacecraft operation for each physics partition of the right-hand side:
`_spacecraft_full_rhs!`, `_spacecraft_slow_rhs!`, `_spacecraft_implicit_rhs!`
(the implicit atmosphere part), `_spacecraft_explicit_rhs!` (the explicit
remainder) and `_spacecraft_fast_control_rhs!`. The serial loop and the
satellite batch loop of each partition call the same operation, so the two
routes cannot drift apart. The partitions stay separate because their physics
and the derivatives they assign differ.

Choosing the execution plan, prefilling the shared samples, sizing and
synchronizing the batches and applying interlinks stay in the outer loops. An
early return skips only the current spacecraft. The bounds-checking context of
the caller passes through to the operation (`Base.@propagate_inbounds`), so the
`@inbounds` route gains no checks. No public name or module ownership changes.

## How the kernels are checked

`test/unit/simulation/rhs_kernel_parity_tests.jl` builds four spacecraft inside
an exponential atmosphere, one of them inactive, and evaluates every partition
with attitude integration off and on, serially and in satellite batches. It
requires the two routes to agree to `1e-13`, the inactive spacecraft's
derivative to be zero, and the implicit-atmosphere and fast-control partitions
to leave the heat loads (and, with attitude, the quaternion) unchanged. With a
single thread the batch route does not actually batch, and the test skips only
its check that batching is enabled.

Allocation checks must use the callback environment snapshot that normal
simulation setup captures. A parameter fixture built by hand without that
snapshot parses the environment on every call instead, so its allocation
figures do not represent a real run.

When the operations were extracted from the old loop bodies, the derivatives,
allocation counts and run times of the original and extracted code were
compared on the threaded full, slow and fast-control routes; the measurements
are kept in the lab's records. They support unchanged derivatives and no
consistent change in run time on those routes. They do not establish a speedup,
zero allocations, or performance acceptance for whole campaigns.
