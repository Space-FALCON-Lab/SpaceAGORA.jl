# Runs test/smoke/ci_examples_suite_smoke.jl on the examples assigned to this CI
# shard (SPACEAGORA_CI_SHARD_ITEMS) and reports which ones it ran. See
# ci_shard_hooks.jl.
include(joinpath(@__DIR__, "ci_shard_hooks.jl"))
CIShard.run_example_shard()
