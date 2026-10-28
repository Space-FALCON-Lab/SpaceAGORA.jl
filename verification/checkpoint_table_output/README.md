# Checkpoint table cadence repair

A checkpoint changes solve bounds, not the requested output cadence. On main
 a664365499c08cd84b1421047cb2208bf2b46e5b, an ordinary1.04s run with0.1s output
and0.37s checkpoints omits seven required epochs(0.4…1.0s), both when segmented
and after interrupted resume with a valid saved prefix. The direct run has all
12required rows. before-observed.toml and before.log reproduce11passes/2failed
missing-grid assertions. No KS package or implementation is involved.

The table saver now constructs canonical epoch k as
`last(range(origin;step=rate,length=k+1))`, independently for each requested
index. Actual timestamp comparisons correct integer-bound guesses. Only this
segment's points are queued; explicit starts/ends remain, with exact-only
deduplication. Cached problem tspan is not authoritative for later segments;
the existing shared segment-end clock is used. New public SavingCallback
instances reset each segment's counter/heap. A typed stable handle lets the
existing pre-event output hook resolve the current inner saver after reinit,
preserving incoming-trajectory samples before state/RHS-changing events.

This policy also handles fractional-cadence mission extension, where a range
constructed against a changing final time can change the bits of an earlier
sample. Normal uninterrupted table saving uses the same canonical policy. This
intentionally corrects rare floating endpoint representations; a nominal0.3s
sample and an explicit nextfloat(0.3) endpoint are both retained. It does not
snap nearby distinct timestamps to a tolerance. Nonfinite/backward/overflowing
or unresolved Float64 clock domains fail clearly rather than queue duplicates.

After repair the unchanged ordinary reproduction passes13/13 and returns14
unique rows in segmented/resumed cases: all12nominal rows plus0.37/0.74
checkpoint endpoints. Fresh child processes verify cold resume and fractional
mission extension(17checks). Valid persisted-prefix behavior follows current
upstream; missing/incompatible-prefix behavior and long streaming are covered
by upstream regression tests, without changing that storage policy.

Selected regression coverage is242existing checkpoint/controller assertions,
266canonical/affine/cached/extension assertions,19event-output assertions,
42storage assertions and115streaming assertions. The archived combined log
contains681assertions with the earlier16-event subset; event-output-final.log
contains the final19-event superset. These are684distinct selected assertions,
not the sum of overlapping reruns. Independent review records additional
source-specific checks and preserved intermediate failures separately.

From the normally instantiated repository, use new output directories:

```sh
julia --startup-file=no --project=. verification/checkpoint_table_output/reproduce.jl /absolute/new/ordinary-reproduction
julia --startup-file=no --project=. verification/checkpoint_table_output/fresh_process.jl /absolute/new/fresh-process
```

The new unit file is registered alongside checkpoint_periodic_tick_tests.jl.
Run the selected event/storage/long-streaming files from the repository's normal
test environment. Independent source/evidence hashes are in independent-review.json.

The claim is limited to the tested standard results-table path and normal ODE
routes. The separate TrajectoryRecorder cadence/duplicate-boundary policy and
specialized split/multirate output paths are not repaired or qualified here.
Existing recorder pre-event behavior is retained. No propagator,force model,
controller or MPC method is replaced, and no KS integration is included.
