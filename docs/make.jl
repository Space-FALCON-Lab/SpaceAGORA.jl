using Documenter
using SpaceAGORA
using Dates

include(joinpath(@__DIR__, "public_api_symbols.jl"))
using .PublicAPISymbols

const REPO_URL = "https://github.com/Space-FALCON-Lab/SpaceAGORA.jl"
const DOCS_SRC = joinpath(@__DIR__, "src")
const GENERATED_SRC = joinpath(DOCS_SRC, "generated")
const GENERATED_API_PAGE = joinpath(GENERATED_SRC, "public_api.md")
const GENERATED_CONTRACTS_ROOT = joinpath(GENERATED_SRC, "contracts")
const BUILD_DIR = joinpath(@__DIR__, "build")
const PUBLIC_API_REPORT = joinpath(BUILD_DIR, "undocumented_public_exports.txt")
const CONTRACT_PAGES = [
    ("architecture", "canonical_topology_contract.md"),
    ("architecture", "final_clean_contract.md"),
    ("architecture", "gnc_aerobraking_boundary_contract.md"),
    ("architecture", "src_canonical_owner_audit.md"),
    ("architecture", "src_completeness_contract.md"),
    ("architecture", "src_restructure_migration_plan.md"),
    ("quality", "api_naming_contract.md"),
    ("quality", "verification_contract.md"),
    ("quality", "performance_runtime_analysis.md"),
]

function _write_generated_public_api_page()::String
    mkpath(GENERATED_SRC)
    open(GENERATED_API_PAGE, "w") do io
        write(io, render_public_api_markdown(SpaceAGORA))
    end
    return GENERATED_API_PAGE
end

function _mirror_contract_docs()
    for (group, filename) in CONTRACT_PAGES
        src_path = joinpath(@__DIR__, group, filename)
        dest_dir = joinpath(GENERATED_CONTRACTS_ROOT, group)
        dest_path = joinpath(dest_dir, filename)
        mkpath(dest_dir)
        write(dest_path, read(src_path, String))
    end
    return nothing
end

function _write_undocumented_public_api_report()::Vector
    mkpath(BUILD_DIR)
    missing = undocumented_public_api_specs(SpaceAGORA)
    open(PUBLIC_API_REPORT, "w") do io
        println(io, "SpaceAGORA undocumented public exports report")
        println(io, "generated_utc=$(Dates.now(Dates.UTC))")
        println(io)
        if isempty(missing)
            println(io, "All documented.")
        else
            for spec in missing
                println(io, "$(spec.section): $(nameof(spec.owner_module)).$(spec.symbol)")
            end
        end
    end
    if isempty(missing)
        @info "All generated public API symbols have docstrings."
    else
        for spec in missing
            @warn "Undocumented public export" section=spec.section owner_module_name=String(nameof(spec.owner_module)) symbol=String(spec.symbol)
        end
    end
    return missing
end

_write_generated_public_api_page()
_mirror_contract_docs()
missing_public_api = undocumented_public_api_specs(SpaceAGORA)

DocMeta.setdocmeta!(SpaceAGORA, :DocTestSetup, :(using SpaceAGORA); recursive=true)

makedocs(
    sitename = "SpaceAGORA.jl",
    modules = [
        SpaceAGORA,
        SpaceAGORA.SimulationEngine,
        SpaceAGORA.ParallelProfiles,
        SpaceAGORA.TelemetryVerification,
    ],
    source = DOCS_SRC,
    build = BUILD_DIR,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://space-falcon-lab.github.io/SpaceAGORA.jl",
        assets = String[],
    ),
    pages = [
        "Overview" => "index.md",
        "Getting Started" => "getting_started.md",
        "Public API" => "generated/public_api.md",
        "Architecture & Quality Contracts" => Any[
            "contracts.md",
            "generated/contracts/architecture/canonical_topology_contract.md",
            "generated/contracts/architecture/final_clean_contract.md",
            "generated/contracts/architecture/gnc_aerobraking_boundary_contract.md",
            "generated/contracts/architecture/src_canonical_owner_audit.md",
            "generated/contracts/architecture/src_completeness_contract.md",
            "generated/contracts/architecture/src_restructure_migration_plan.md",
            "generated/contracts/quality/api_naming_contract.md",
            "generated/contracts/quality/verification_contract.md",
            "generated/contracts/quality/performance_runtime_analysis.md",
        ],
    ],
    doctest = false,
    linkcheck = true,
    checkdocs = :none,
)

_write_undocumented_public_api_report()

if !isempty(missing_public_api)
    @warn "Public API documentation is incomplete. See undocumented_public_exports.txt for the current warning list." report=PUBLIC_API_REPORT count=length(missing_public_api)
end

deploydocs(
    repo = "github.com/Space-FALCON-Lab/SpaceAGORA.jl.git",
    devbranch = "main",
    versions = nothing,
    push_preview = false,
)
