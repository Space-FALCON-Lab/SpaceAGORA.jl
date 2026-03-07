using SpaceAGORA

const REPO_ROOT = dirname(dirname(@__FILE__))

include(joinpath(REPO_ROOT, "docs", "public_api_symbols.jl"))
using .PublicAPISymbols

documented = Set(spec.symbol for spec in public_api_specs(SpaceAGORA))
exported = Set(Symbol.(names(SpaceAGORA)))
delete!(exported, :SpaceAGORA)

missing_from_docs = sort!(collect(setdiff(exported, documented)))
missing_from_exports = sort!(collect(setdiff(documented, exported)))

isempty(missing_from_docs) || error("Exported public symbols missing from docs/public_api_symbols.jl: $(missing_from_docs)")
isempty(missing_from_exports) || error("Documented public symbols missing from SpaceAGORA exports: $(missing_from_exports)")

for internal_name in (:SimulationModel, :SimulationEngine, :ParallelProfiles, :TelemetryVerification, :SpaceAGORACLI)
    internal_name in exported && error("Internal module $(internal_name) must not be exported from SpaceAGORA")
end

println("public_api_surface_gate_ok")
