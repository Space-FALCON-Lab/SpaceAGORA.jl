const REPO_ROOT = dirname(dirname(@__FILE__))

pr_template_path = joinpath(REPO_ROOT, ".github", "pull_request_template.md")
doc_policy_path = joinpath(REPO_ROOT, "docs", "src", "documentation_policy.md")
docs_make_path = joinpath(REPO_ROOT, "docs", "make.jl")

for required in (pr_template_path, doc_policy_path, docs_make_path)
    isfile(required) || error("Missing documentation process file: $(required)")
end

pr_template = read(pr_template_path, String)
doc_policy = read(doc_policy_path, String)
docs_make = read(docs_make_path, String)

for marker in (
    "If root exports changed, I updated `docs/public_api_symbols.jl`",
    "If assets changed, I updated `data/assets_manifest.toml` and `docs/src/assets.md`",
    "If new extension points were added, I updated `docs/src/extensibility.md`",
    "If new CLI commands or flags were added, I updated `docs/src/cli.md`",
    "`julia --project=docs docs/make.jl`",
    "`julia --project=.AGORA test/ci_public_api_surface_gate.jl`",
    "`julia --project=.AGORA test/ci_hpc_extensibility_docs_gate.jl`",
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

occursin("\"Documentation Policy\" => \"documentation_policy.md\"", docs_make) ||
    error("docs/make.jl is missing the Documentation Policy page entry.")

println("docs_contract_gate_ok")
