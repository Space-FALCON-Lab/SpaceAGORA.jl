# Filesystem lookup only; importing this owner starts no benchmark.
module BenchmarkArtifacts
export _latest_artifact_path, _latest_artifact_path_optional
function _latest_artifact_path_optional(outdir::String, prefix::String, profile_name::String, suffix::String)::Union{Nothing, String}
    pattern = "$(prefix)_$(profile_name)_"
    candidates = String[]
    for entry in readdir(outdir)
        if startswith(entry, pattern) && endswith(entry, suffix)
            push!(candidates, joinpath(outdir, entry))
        end
    end
    isempty(candidates) && return nothing
    sort!(candidates; by=path -> stat(path).mtime)
    return last(candidates)
end

function _latest_artifact_path(outdir::String, prefix::String, profile_name::String, suffix::String)::String
    path = _latest_artifact_path_optional(outdir, prefix, profile_name, suffix)
    path === nothing && error("No artifact found in '$outdir' for prefix='$prefix' profile='$(profile_name)' suffix='$suffix'")
    return path
end
end
