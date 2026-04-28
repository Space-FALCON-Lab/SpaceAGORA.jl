include(joinpath(@__DIR__, "pr_runtests.jl"))
include(joinpath(@__DIR__, "nightly_runtests.jl"))

# Full contract superset used for local runs.
include(joinpath(@__DIR__, "..", "ci_final_clean_contract_gate.jl"))
include(joinpath(@__DIR__, "..", "ci_gnc_typed_command_boundary_gate.jl"))
include(joinpath(@__DIR__, "..", "ci_io_surface_nonempty_gate.jl"))
include(joinpath(@__DIR__, "..", "ci_project_compat_gate.jl"))
include(joinpath(@__DIR__, "..", "ci_runtime_legacy_surface_gate.jl"))
include(joinpath(@__DIR__, "..", "ci_thin_entry_files_gate.jl"))
include(joinpath(@__DIR__, "..", "ci_translational_ownership_gate.jl"))
