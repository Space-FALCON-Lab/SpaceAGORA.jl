# Energy depletion guidance (EDG) and KS aerobraking MPC tests. As a numbered suite they run in one CI
# shard per mode; included from the harness preamble they ran in every shard.
# First in the suite list, so they still run before every other suite.
include(joinpath(REPO_ROOT, "test", "gnc", "aerobraking", "energy_depletion_gnc_tests.jl"))
include(joinpath(REPO_ROOT, "test", "gnc", "aerobraking", "edg_pass_heat_tests.jl"))
include(joinpath(REPO_ROOT, "test", "gnc", "aerobraking", "edg_heat_observer_tests.jl"))
include(joinpath(REPO_ROOT, "test", "gnc", "aerobraking", "edg_journal_contract_tests.jl"))
include(joinpath(REPO_ROOT, "test", "gnc", "aerobraking", "edg_replan_callback_tests.jl"))
include(joinpath(REPO_ROOT, "test", "gnc", "aerobraking", "energy_depletion_lifecycle_tests.jl"))
include(joinpath(REPO_ROOT, "test", "gnc", "aerobraking", "energy_depletion_ownership_tests.jl"))
include(joinpath(REPO_ROOT, "test", "gnc", "aerobraking", "edg_guidance_family_tests.jl"))
include(joinpath(REPO_ROOT, "test", "unit", "dynamics", "ks_dynamics_tests.jl"))
include(joinpath(REPO_ROOT, "test", "gnc", "aerobraking", "mpc_campaign_tests.jl"))
