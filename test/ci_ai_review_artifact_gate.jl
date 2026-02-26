const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

const REQUIRED_SECTIONS = [
    "Scope",
    "Changed Files",
    "Findings",
    "P1 Assessment",
    "Tests Added/Updated",
    "Residual Risk",
]

@inline normalize_heading(s::AbstractString) = strip(replace(s, r"\s+" => " "))

function parse_pr_number()
    if haskey(ENV, "SPACEAGORA_PR_NUMBER")
        pr = tryparse(Int, ENV["SPACEAGORA_PR_NUMBER"])
        pr === nothing && error("Invalid SPACEAGORA_PR_NUMBER=$(ENV["SPACEAGORA_PR_NUMBER"])")
        return pr
    end

    ref = get(ENV, "GITHUB_REF", "")
    m = match(r"refs/pull/([0-9]+)/", ref)
    m === nothing && error("Unable to determine PR number. Set SPACEAGORA_PR_NUMBER or run in pull_request context.")
    return parse(Int, m.captures[1])
end

function parse_sections(text::String)
    sections = Dict{String, Vector{String}}()
    current = nothing
    for raw in split(text, '\n'; keepempty=true)
        m = match(r"^\s*##\s+(.+?)\s*$", raw)
        if m !== nothing
            current = normalize_heading(m.captures[1])
            sections[current] = String[]
            continue
        end
        current === nothing && continue
        push!(sections[current], raw)
    end
    return sections
end

function parse_changed_files_from_artifact(section_lines::Vector{String})
    files = Set{String}()
    for line in section_lines
        for m in eachmatch(r"src/[A-Za-z0-9_\-./]+\.jl", line)
            push!(files, normpath(m.match))
        end
    end
    return files
end

function git_changed_files(base_ref::Union{Nothing, String})
    function run_diff(cmd::Cmd)
        out = IOBuffer()
        err = IOBuffer()
        proc = run(pipeline(ignorestatus(cmd), stdout=out, stderr=err))
        if !success(proc)
            return nothing
        end
        text = String(take!(out))
        return filter(!isempty, strip.(split(text, '\n'; keepempty=false)))
    end

    candidates = Cmd[]
    if base_ref !== nothing && !isempty(base_ref)
        push!(candidates, Cmd(["git", "diff", "--name-only", "--diff-filter=AMR", "origin/$base_ref...HEAD"]))
    end
    push!(candidates, Cmd(["git", "diff", "--name-only", "--diff-filter=AMR", "HEAD~1...HEAD"]))

    changed = nothing
    for cmd in candidates
        changed = run_diff(cmd)
        changed === nothing || break
    end

    changed === nothing && error("Unable to compute changed files via git diff.")

    src_jl = Set{String}()
    for path in changed
        rel = replace(path, "\\" => "/")
        if startswith(rel, "src/") && endswith(rel, ".jl")
            push!(src_jl, normpath(rel))
        end
    end

    return src_jl
end

function main()
    pr_number = parse_pr_number()
    artifact_path = joinpath(REPO_ROOT, "test", "ai_reviews", "PR_$(pr_number).md")
    isfile(artifact_path) || error("Missing AI review artifact: $(relpath(artifact_path, REPO_ROOT))")

    artifact_text = read(artifact_path, String)
    sections = parse_sections(artifact_text)
    section_keys = Set(keys(sections))

    missing_sections = [name for name in REQUIRED_SECTIONS if !(name in section_keys)]
    isempty(missing_sections) || error("AI review artifact missing required sections: $(join(missing_sections, ", "))")

    for name in REQUIRED_SECTIONS
        body = strip(join(sections[name], "\n"))
        isempty(body) && error("AI review artifact section '$name' is empty.")
    end

    listed_files = parse_changed_files_from_artifact(sections["Changed Files"])
    base_ref = get(ENV, "SPACEAGORA_BASE_REF", get(ENV, "GITHUB_BASE_REF", nothing))
    changed_src_files = git_changed_files(base_ref)

    missing_file_coverage = sort(collect(setdiff(changed_src_files, listed_files)))
    if !isempty(missing_file_coverage)
        println("changed_src_files=$(length(changed_src_files)) listed_files=$(length(listed_files))")
        println("missing_changed_files_in_artifact:")
        for path in missing_file_coverage
            println("  - $path")
        end
        error("AI review artifact does not cover all changed src/*.jl files.")
    end

    println("ai_review_artifact_gate_ok pr=$(pr_number) changed_src_files=$(length(changed_src_files))")
end

main()
