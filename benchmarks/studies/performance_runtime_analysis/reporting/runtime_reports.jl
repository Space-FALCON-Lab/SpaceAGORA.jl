# Reporting and rollout decisions, directly importable without simulation setup.
# Runtime supplies report provenance and the two simulation callbacks explicitly.
module RuntimeReports
using CSV, DataFrames, Dates, Statistics, LinearAlgebra
export ProfileSpec, summarize_results, summarize_per_orbit_results, summarize_entry_duration_results, summarize_density_backend_breakdown, selected_cases
include("profile.jl")
include("settings.jl")
include("formatting.jl")
include("case_selection.jl")
include("summaries.jl")
include("rollout_gates.jl")
include("report_writing.jl")
end
