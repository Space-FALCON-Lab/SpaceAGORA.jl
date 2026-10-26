using SpaceAGORA

const REPO_ROOT = dirname(dirname(dirname(@__FILE__)))

include(joinpath(REPO_ROOT, "docs", "public_api_symbols.jl"))
using .PublicAPISymbols

specs = public_api_specs(SpaceAGORA)
# The public Conformance module exports its own names; its section lists exactly those.
documented = Set(spec.symbol for spec in specs if spec.owner !== :Conformance)
exported = Set(Symbol.(names(SpaceAGORA)))
delete!(exported, :SpaceAGORA)
conformance_documented = Set(spec.symbol for spec in specs if spec.owner === :Conformance)
conformance_exported = setdiff(Set(Symbol.(names(SpaceAGORA.Conformance))), (:Conformance,))
conformance_documented == conformance_exported || error("docs/public_api_symbols.jl and the exports of SpaceAGORA.Conformance differ: " *
    "undocumented $(sort!(collect(setdiff(conformance_exported, conformance_documented)))), " *
    "not exported $(sort!(collect(setdiff(conformance_documented, conformance_exported))))")

missing_from_docs = sort!(collect(setdiff(exported, documented)))
missing_from_exports = sort!(collect(setdiff(documented, exported)))

isempty(missing_from_docs) || error("Exported public symbols missing from docs/public_api_symbols.jl: $(missing_from_docs)")
isempty(missing_from_exports) || error("Documented public symbols missing from SpaceAGORA exports: $(missing_from_exports)")

for internal_name in (:RuntimeServices, :SimulationModel, :SimulationEngine, :ParallelProfiles, :TelemetryVerification, :SpaceAGORACLI)
    internal_name in exported && error("Internal module $(internal_name) must not be exported from SpaceAGORA")
end

println("public_api_surface_gate_ok")
