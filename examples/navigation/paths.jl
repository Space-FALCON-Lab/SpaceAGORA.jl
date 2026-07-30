module NavigationPaths

export REPO_ROOT, NAVIGATION_OUTPUT_ROOT
export navigation_output_path, transient_navigation_path
export env_path, env_override_path, resolve_repo_path
export stored_output_path, resolve_output_path

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const NAVIGATION_OUTPUT_ROOT = joinpath(REPO_ROOT, "output", "navigation")

"""
    resolve_repo_path(path)

Resolve user-provided relative paths from the repository root. Absolute paths
remain valid overrides.
"""
function resolve_repo_path(path::AbstractString)::String
    value = strip(String(path))
    isempty(value) && return REPO_ROOT
    return isabspath(value) ? normpath(value) : normpath(joinpath(REPO_ROOT, value))
end

navigation_output_path(parts::AbstractString...) =
    normpath(joinpath(NAVIGATION_OUTPUT_ROOT, parts...))

transient_navigation_path(parts::AbstractString...) =
    normpath(joinpath(tempdir(), "spaceagora_navigation", parts...))

"""
    env_path(name, default_parts...)

Read a path override from `ENV[name]`. Relative overrides are resolved from the
repository root; otherwise the path defaults to
`output/navigation/default_parts...`.
"""
function env_path(name::AbstractString, default_parts::AbstractString...)::String
    value = strip(get(ENV, String(name), ""))
    return isempty(value) ?
        navigation_output_path(default_parts...) :
        resolve_repo_path(value)
end

"""
    env_override_path(name, default_path)

Read an optional path override while preserving a default path that may already
depend on another configured output root.
"""
function env_override_path(
    name::AbstractString,
    default_path::AbstractString
)::String
    value = strip(get(ENV, String(name), ""))
    return isempty(value) ? normpath(default_path) : resolve_repo_path(value)
end

"""
    stored_output_path(path, output_root)

Return a portable path for status and manifest files. Paths inside an output
root are stored relative to that root.
"""
function stored_output_path(
    path::AbstractString,
    output_root::AbstractString
)::String
    candidate = normpath(path)
    relative = relpath(candidate, normpath(output_root))
    parts = splitpath(relative)
    return !isempty(parts) && first(parts) == ".." ? candidate : relative
end

"""
    resolve_output_path(path, output_root)

Resolve paths read from status files. Absolute paths from older campaigns are
accepted for backward compatibility. If an old absolute or repository-relative
path no longer exists, its suffix below the former `output_navigation_*` or
`output_observer_*` root is relocated below `output_root`.
"""
function _legacy_output_suffix(path::AbstractString)::Union{Nothing, String}
    parts = splitpath(normpath(String(path)))
    root_index = findfirst(
        part -> startswith(part, "output_navigation_") ||
                startswith(part, "output_observer_"),
        parts
    )
    root_index === nothing && return nothing
    suffix_parts = parts[(root_index + 1):end]
    return isempty(suffix_parts) ? "" : joinpath(suffix_parts...)
end

function resolve_output_path(path::AbstractString, output_root::AbstractString)::String
    value = normpath(String(path))
    if isabspath(value) && ispath(value)
        return value
    end

    legacy_suffix = _legacy_output_suffix(value)
    legacy_suffix !== nothing &&
        return normpath(joinpath(output_root, legacy_suffix))

    return isabspath(value) ? value : normpath(joinpath(output_root, value))
end

end # module NavigationPaths
