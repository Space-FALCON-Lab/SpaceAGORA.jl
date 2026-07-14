# The test/suites/NN_*.jl content that used to run from here has been fully
# migrated into test/unit/ and test/contracts/ (see test/TEST_RESTRUCTURE_PLAN.md,
# Phase 2) -- test/suites/ has been deleted. This file is kept as the
# bootstrap trigger point in test/runtests.jl's default chain, and as the
# home for future genuine end-to-end integration tests (examples, CLI,
# persistence, telemetry -- see the subdirectory READMEs here).
include(joinpath(@__DIR__, "..", "helpers", "bootstrap.jl"))
