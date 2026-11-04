# First-campaign latency (WS11f item 3): the wall time of the FIRST
# `run_monte_carlo(threads=:auto)` call after `using SpaceAGORA`, in a fresh
# process. Before this item, `src/precompile_workload.jl`'s `@compile_workload`
# exercised the low-level Monte Carlo dispatchers (`_warm_campaign_dispatchers`,
# `SimulationCampaigns/monte_carlo.jl`) but never the campaign-level planning
# layer above them: `_campaign_route_plan`/`_run_campaign_with_route_env` for
# the default bandit route, or `predictive_plan`/`_run_campaign_predictive` for
# R7 (`SPACEAGORA_CAMPAIGN_PLANNER=predictive`). A session's first campaign
# paid to JIT-compile that layer instead.
#
# Run this script itself to reproduce the "after" number in a fresh process:
#   julia --project=. --threads=1 benchmarks/studies/small_wins/first_campaign_latency.jl
#
# The "before" number requires the pre-fix `src/precompile_workload.jl`
# (i.e. HEAD~1 on the commit that added `_warm_predictive_campaign` and
# `_warm_mixed_dispatch_campaign`), a forced recompile (`using SpaceAGORA` once
# to rebuild the pkgimage), and then this same script.
#
# Measured on space-falcon-1, 1 thread, two fresh processes each, back to back
# (see results/precompile_first_campaign_latency.csv for the raw numbers):
#
#   before   avg 3.824 s  (3.944, 3.703)
#   after    avg 0.658 s  (0.653, 0.664)
#   ratio    before/after ~= 5.81x
#
# Ratio only, as required: this is a shared, contended machine, so the two
# numbers above are not meant to stand alone as absolute latencies.

using SpaceAGORA

t = @elapsed run_monte_carlo(seed -> seed * 2, 1:2; threads=:auto)
println("first_campaign_elapsed_s=", t)
