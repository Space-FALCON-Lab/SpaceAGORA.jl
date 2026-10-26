const REPO_ROOT = dirname(dirname(dirname(@__FILE__)))

pr_template_path = joinpath(REPO_ROOT, ".github", "pull_request_template.md")
doc_policy_path = joinpath(REPO_ROOT, "docs", "contributing", "documentation_policy.md")
docs_make_path = joinpath(REPO_ROOT, "docs", "make.jl")
docs_index_path = joinpath(REPO_ROOT, "docs", "src", "index.md")
contributor_index_path = joinpath(REPO_ROOT, "docs", "contributing", "index.md")
contributing_path = joinpath(REPO_ROOT, "CONTRIBUTING.md")
user_page_paths = [
    joinpath(REPO_ROOT, "docs", "src", "user", "quickstart.md"),
    joinpath(REPO_ROOT, "docs", "src", "user", "installation_environment.md"),
    joinpath(REPO_ROOT, "docs", "src", "user", "first_script.md"),
    joinpath(REPO_ROOT, "docs", "src", "user", "verification_study.md"),
    joinpath(REPO_ROOT, "docs", "src", "user", "recipes.md"),
    joinpath(REPO_ROOT, "docs", "src", "user", "concepts.md"),
]

for required in (pr_template_path, doc_policy_path, docs_make_path, docs_index_path, contributor_index_path, contributing_path, user_page_paths...)
    isfile(required) || error("Missing documentation process file: $(required)")
end

pr_template = read(pr_template_path, String)
doc_policy = read(doc_policy_path, String)
docs_make = read(docs_make_path, String)
docs_index = read(docs_index_path, String)

for marker in (
    "If root exports changed, I updated `docs/public_api_symbols.jl`",
    "If assets changed, I updated `data/assets_manifest.toml` and `docs/src/assets.md`",
    "If new extension points were added, I updated `docs/src/extensibility.md`",
    "If a hook was added, renamed or deprecated, I updated `docs/extension_hooks.jl`",
    "If new CLI commands or flags were added, I updated `docs/src/cli.md`",
    "`julia --startup-file=no docs/setup.jl`, then `julia --project=docs docs/make.jl`",
    "`julia --project=. test/gates/ci_public_api_surface_gate.jl`",
    "`julia --project=. test/gates/ci_hpc_extensibility_docs_gate.jl`",
)
    occursin(marker, pr_template) || error("PR template missing docs drift/release marker: $(marker)")
end

for marker in (
    "Generated API docs",
    "Contract docs",
    "Operational docs",
    "docs/public_api_symbols.jl",
    "data/assets_manifest.toml",
    # The docs environment's manifest is generated, so the build starts with its setup.
    "Before merge to `main`, the documentation release gate is:\n\n```text\n" *
        "julia --startup-file=no docs/setup.jl\njulia --project=docs docs/make.jl\n",
)
    occursin(marker, doc_policy) || error("documentation_policy.md missing marker: $(marker)")
end

for marker in (
    "\"Home\" => \"index.md\"",
    "\"Get Started\" => Any[",
    "\"Quickstart\" => \"user/quickstart.md\"",
    "\"Installation & Environment\" => \"user/installation_environment.md\"",
    "\"Build a Simulation\" => Any[",
    "\"Your First Script\" => \"user/first_script.md\"",
    "\"Assets & Modes\" => \"assets.md\"",
    "\"Help\" => Any[",
    "\"Troubleshooting\" => \"user/troubleshooting.md\"",
    "\"Extend SpaceAGORA\" => Any[",
    "\"Extension Hook Reference\" => \"generated/extension_hook_reference.md\"",
    "\"API Reference\" => Any[",
    "\"HYPR Service Contract\" => \"reference/hypr_services.md\"",
    "\"RPO Planner Contract\" => \"reference/rpo_planner_contract.md\"",
    "\"External Propagation Services\" => \"reference/external_propagation_services.md\"",
    "\"EDG Service Contract\" => \"reference/edg_adapter_boundary.md\"",
    "\"Optional Adapter Contract\" => \"reference/adapter_contract.md\"",
    "\"Public API and Compatibility\" => \"public_api_policy.md\"",
    "\"Extension Conventions\" => \"reference/extension_conventions.md\"",
    "\"Migrating Extension Hooks\" => \"reference/migrating_extension_hooks.md\"",
    "\"Writing a GNC Algorithm\" => \"user/gnc_algorithm.md\"",
    "\"Building Your Own Research Repository\" => \"user/research_repository.md\"",
)
    occursin(marker, docs_make) || error("docs/make.jl missing documentation IA marker: $(marker)")
end

# Contributor requirements stay enforced in the repository, outside the public site.
for path in (
    joinpath(REPO_ROOT, "docs", "src", "maintainer"),
    joinpath(REPO_ROOT, "docs", "src", "documentation_policy.md"),
    joinpath(REPO_ROOT, "docs", "src", "contracts.md"),
)
    ispath(path) && error("Contributor material must stay outside published docs/src: $(path)")
end
for marker in ("Maintainer Guide", "maintainer/index.md", "_mirror_contract_docs")
    occursin(marker, docs_make) && error("docs/make.jl republishes contributor material: $(marker)")
end
for relative in (
    "docs/contributing/index.md",
    "docs/contributing/architecture.md",
    "docs/contributing/configuration_ownership.md",
    "docs/contributing/documentation_policy.md",
    "docs/contributing/contracts.md",
    "docs/contributing/adapter_maintenance.md",
    "docs/contributing/gnc_ownership.md",
    "docs/contributing/gnc_algorithm_families.md",
    "docs/contributing/contributor_kit_integrations.md",
    "docs/architecture/canonical_topology_contract.md",
    "docs/architecture/gnc_aerobraking_boundary_contract.md",
    "docs/architecture/src_canonical_owner_audit.md",
    "docs/architecture/src_completeness_contract.md",
    "docs/quality/api_naming_contract.md",
    "docs/quality/verification_contract.md",
    "docs/quality/performance_runtime_analysis.md",
)
    isfile(joinpath(REPO_ROOT, relative)) || error("Missing retained contributor contract: $(relative)")
end
occursin("docs/contributing/index.md", read(contributing_path, String)) ||
    error("CONTRIBUTING.md must link the repository contributor guide.")
occursin("## Contributor checklist", read(contributor_index_path, String)) ||
    error("Contributor guide must retain the review checklist.")

!occursin("Markdown.parse(read(readme_path, String))", docs_index) ||
    error("docs/src/index.md still mirrors README content instead of using a hand-written docs landing page.")

println("docs_contract_gate_ok")
