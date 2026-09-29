# Entry point for the unit-suite subprocess on a sharded CI run: suite 09's unit
# driver starts this instead of test/unit/runtests.jl, and it includes only the
# unit files assigned to this shard. See ci_shard_hooks.jl.
include(joinpath(@__DIR__, "ci_shard_hooks.jl"))
CIShard.run_unit_shard()
