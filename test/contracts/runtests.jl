# Full contract superset used for local runs: the PR gates plus the gates that
# are not wired into CI.
include(joinpath(@__DIR__, "pr_runtests.jl"))
include(joinpath(@__DIR__, "..", "gates", "ci_gnc_typed_command_boundary_gate.jl"))
include(joinpath(@__DIR__, "..", "gates", "ci_project_compat_gate.jl"))
include(joinpath(@__DIR__, "..", "gates", "ci_runtime_legacy_surface_gate.jl"))
include(joinpath(@__DIR__, "..", "gates", "ci_translational_ownership_gate.jl"))
