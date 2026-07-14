The legacy numbered-suite runtime harness that used to live here has been
fully migrated to `test/unit/` and `test/contracts/` (see
`test/TEST_RESTRUCTURE_PLAN.md`). `runtests.jl` here is now just the
`test/helpers/bootstrap.jl` trigger point in the default `test/runtests.jl`
chain.

This directory remains the intended home for genuine end-to-end integration
tests, in these sub-buckets:

- `simulation/`
- `persistence/`
- `cli/`
- `examples/` (already has real content: `runtests.jl` runs the examples regression suite)
- `telemetry/`

The other sub-buckets are currently empty placeholders — add real
integration tests there as they're written, following the same pattern as
`examples/`.
