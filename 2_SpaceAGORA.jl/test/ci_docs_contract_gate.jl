const REPO_ROOT = dirname(dirname(@__FILE__))

pr_template_path = joinpath(REPO_ROOT, ".github", "pull_request_template.md")
doc_policy_path = joinpath(REPO_ROOT, "docs", "src", "documentation_policy.md")
docs_make_path = joinpath(REPO_ROOT, "docs", "make.jl")
docs_index_path = joinpath(REPO_ROOT, "docs", "src", "index.md")
maintainer_index_path = joinpath(REPO_ROOT, "docs", "src", "maintainer", "index.md")
user_page_paths = [
    joinpath(REPO_ROOT, "docs", "src", "user", "quickstart.md"),
    joinpath(REPO_ROOT, "docs", "src", "user", "installation_environment.md"),
    joinpath(REPO_ROOT, "docs", "src", "user", "first_simulation.md"),
    joinpath(REPO_ROOT, "docs", "src", "user", "verification_study.md"),
    joinpath(REPO_ROOT, "docs", "src", "user", "recipes.md"),
    joinpath(REPO_ROOT, "docs", "src", "user", "concepts.md"),
]

for required in (pr_template_path, doc_policy_path, docs_make_path, docs_index_path, maintainer_index_path, user_page_paths...)
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
    "If new CLI commands or flags were added, I updated `docs/src/cli.md`",
    "`julia --project=. docs/make.jl`",
    "`julia --project=. test/ci_public_api_surface_gate.jl`",
    "`julia --project=. test/ci_hpc_extensibility_docs_gate.jl`",
)
    occursin(marker, pr_template) || error("PR template missing docs drift/release marker: $(marker)")
end

for marker in (
    "Generated API docs",
    "Contract docs",
    "Operational docs",
    "docs/public_api_symbols.jl",
    "data/assets_manifest.toml",
    "Before merge to `main`, the documentation release gate is:",
)
    occursin(marker, doc_policy) || error("documentation_policy.md missing marker: $(marker)")
end

for marker in (
    "\"User Guide\" => Any[",
    "\"Start Here\" => \"getting_started.md\"",
    "\"Quickstart\" => \"user/quickstart.md\"",
    "\"Installation & Environment\" => \"user/installation_environment.md\"",
    "\"Assets & Modes\" => \"assets.md\"",
    "\"Reference\" => Any[",
    "\"Maintainer Guide\" => Any[",
    "\"Maintainer Overview\" => \"maintainer/index.md\"",
)
    occursin(marker, docs_make) || error("docs/make.jl missing documentation IA marker: $(marker)")
end

occursin("\"Documentation Policy\" => \"documentation_policy.md\"", docs_make) ||
    error("docs/make.jl is missing the Documentation Policy page entry.")

!occursin("Markdown.parse(read(readme_path, String))", docs_index) ||
    error("docs/src/index.md still mirrors README content instead of using a hand-written docs landing page.")

println("docs_contract_gate_ok")
