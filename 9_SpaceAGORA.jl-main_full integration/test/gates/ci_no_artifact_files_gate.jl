const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function _tracked_paths(cmd::Cmd)::Vector{String}
    output = read(cmd, String)
    lines = split(chomp(output), '\n')
    return [line for line in lines if !isempty(strip(line))]
end

artifact_paths = _tracked_paths(`sh -lc "cd '$REPO_ROOT' && git ls-files '*.cov' '.DS_Store' 'output/**' 'docs/build/**' 'docs/site/**' 'docs/src/generated/**'"`)
!isempty(artifact_paths) && error("Tracked generated/artifact files must be removed:\n$(join(artifact_paths, '\n'))")

text_paths = _tracked_paths(`sh -lc "cd '$REPO_ROOT' && git ls-files '*.md' '*.csv' '*.json' '*.toml' '*.txt' '*.html' '*.env'"`)
abs_path_patterns = (
    r"/Users/[^/[:space:]]+",
    r"/home/[^/[:space:]]+",
    r"[A-Za-z]:\\\\Users\\\\[^\\\\[:space:]]+"
)
violations = String[]
for relpath in text_paths
    fullpath = joinpath(REPO_ROOT, relpath)
    isfile(fullpath) || continue
    for (line_no, raw_line) in enumerate(eachline(fullpath))
        for pattern in abs_path_patterns
            if occursin(pattern, raw_line)
                push!(violations, "$relpath:$line_no contains committed absolute-path/user-specific content")
                break
            end
        end
    end
end

!isempty(violations) && error("Committed text artifacts must not contain absolute paths or usernames:\n$(join(violations, '\n'))")

println("no_artifact_files_gate_ok")
