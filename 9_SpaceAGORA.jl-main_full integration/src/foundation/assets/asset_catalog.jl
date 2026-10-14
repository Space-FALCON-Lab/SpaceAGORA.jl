"""
Repository asset manifest and readiness operations. This internal service
returns plain field records so command front ends retain their own public types.
It performs no downloads and does not interpret atmosphere or geometry schemas.
"""
module AssetCatalog

using TOML

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const ASSET_STATES = (:usable, :missing, :empty, :pointer_only, :invalid)

function load_asset_manifest(; repo_root::String=REPO_ROOT,
    manifest_path::String=joinpath(repo_root, "data", "assets_manifest.toml"))
    raw = TOML.parsefile(manifest_path)
    entries = NamedTuple[]
    for entry in get(raw, "asset", Any[])
        push!(entries, (
            name=String(entry["name"]),
            scope=String(entry["scope"]),
            relative_path=String(get(entry, "relative_path", ".")),
            kind=String(get(entry, "kind", "directory")),
            required=Bool(get(entry, "required", false)),
            licensing=String(get(entry, "licensing", "unspecified")),
            detail=String(get(entry, "detail", "")),
        ))
    end
    return entries
end

@inline function manifest_entry_path(repo_root::String, entry)
    entry.relative_path in ("", ".") && return repo_root
    return joinpath(repo_root, splitpath(entry.relative_path)...)
end

# A Git LFS pointer is a short text stub left in place of a payload that was not
# downloaded (for example after a clone with GIT_LFS_SKIP_SMUDGE=1). The preset
# verifier in the environment models recognises the same header for its own
# payloads; this check only classifies, it never fetches.
const _LFS_POINTER_HEADER = "version https://git-lfs.github.com/spec/v1"
const _LFS_POINTER_MAX_BYTES = 1024

# Git worktrees and initialized submodules use a .git file instead of a
# directory. Finding that marker needs no Git process or network access.
function _asset_git_root(path::String)
    candidate = isdir(path) ? path : dirname(path)
    while true
        marker = joinpath(candidate, ".git")
        (isdir(marker) || isfile(marker)) && return candidate
        parent = dirname(candidate)
        parent == candidate && return nothing
        candidate = parent
    end
end

_asset_shell_quote(value::String) = "'" * replace(value, "'" => "'\"'\"'") * "'"

function _lfs_recovery_hint(owner::String, include_pattern::String; windows::Bool=Sys.iswindows())
    if windows
        # Windows users may run cmd.exe, PowerShell or Git Bash. Show the
        # owning checkout and Git pattern as data instead of guessing a shell.
        pattern = replace(include_pattern, '\\' => '/')
        return "retrieve the payload through the documented asset route; owning Git checkout: $(repr(owner)); Git LFS include pattern: $(repr(pattern))"
    end
    command = "git -C $(_asset_shell_quote(owner)) lfs pull --include=$(_asset_shell_quote(include_pattern))"
    return "fetch the payload with `$(command)` or the documented asset route"
end

function _lfs_fetch_hint(path::String; asset_directory::Union{Nothing,String}=nothing)
    payload = realpath(path)
    owner = _asset_git_root(payload)
    isnothing(owner) && return "retrieve the payload through the documented asset route (no Git checkout found for this path)"
    include_pattern = relpath(payload, owner)
    if !isnothing(asset_directory) && _asset_git_root(asset_directory) == owner
        # Collapse ordinary directory pointers into one hint, but preserve the
        # actual owner and exact path for file symlinks outside this asset root.
        relative_parts = splitpath(relpath(payload, asset_directory))
        if first(relative_parts) != ".."
            directory_pattern = relpath(asset_directory, owner)
            include_pattern = directory_pattern == "." ? "*" : joinpath(directory_pattern, "*")
        end
    end
    return _lfs_recovery_hint(owner, include_pattern)
end

function _is_lfs_pointer(path::String)::Bool
    size = filesize(path)
    0 < size <= _LFS_POINTER_MAX_BYTES || return false
    header = open(io -> read(io, min(ncodeunits(_LFS_POINTER_HEADER), size)), path)
    return String(header) == _LFS_POINTER_HEADER
end

function _file_readiness(path::String)
    isfile(path) || return (:missing, "no file at this path")
    filesize(path) == 0 && return (:invalid, "the file is empty (0 bytes)")
    _is_lfs_pointer(path) && return (:pointer_only, "the file is an undownloaded Git LFS pointer; " * _lfs_fetch_hint(path))
    return (:usable, "")
end

# Hidden files and directories never count as payload. Symlinks to files
# count as whatever they point at; the walk does not follow symlinked
# directories.
function _directory_readiness(path::String)
    isdir(path) || return (:missing, "no directory at this path")
    payloads = 0
    pointers = 0
    empties = 0
    fetch_hints = String[]
    canonical_directory = realpath(path)
    for (root, directories, files) in walkdir(path)
        filter!(directory -> !startswith(directory, "."), directories)
        for file in files
            startswith(file, ".") && continue
            full = joinpath(root, file)
            isfile(full) || continue
            if filesize(full) == 0
                empties += 1
            elseif _is_lfs_pointer(full)
                pointers += 1
                hint = _lfs_fetch_hint(full; asset_directory=canonical_directory)
                hint in fetch_hints || push!(fetch_hints, hint)
            else
                payloads += 1
            end
        end
    end
    total = payloads + pointers + empties
    total == 0 && return (:empty, "the directory contains no files")
    if pointers > 0
        detail = payloads > 0 ? " alongside $(payloads) downloaded payload file$(payloads == 1 ? "" : "s")" : ""
        return (:pointer_only, "$(pointers) of $(total) files are undownloaded Git LFS pointers$(detail); " * join(fetch_hints, "; "))
    end
    payloads == 0 && return (:invalid, "all $(total) files are empty (0 bytes)")
    return (:usable, "")
end

"""
    manifest_entry_readiness(repo_root, entry) -> (state, reason)

Classify a manifest entry without fetching anything. Directory and file
entries distinguish a missing path, an empty directory, undownloaded Git LFS
pointers, empty payloads and usable payloads. Built-in entries are always
usable. The check confirms that payload bytes exist; it does not validate their
contents, which the named-preset `assets check --preset` route does for grids.
"""
function manifest_entry_readiness(repo_root::String, entry)
    path = manifest_entry_path(repo_root, entry)
    entry.kind == "builtin" && return (:usable, "")
    entry.kind == "directory" && return _directory_readiness(path)
    entry.kind == "file" && return _file_readiness(path)
    throw(ArgumentError("Unsupported asset manifest kind '$(entry.kind)' for $(entry.name)."))
end

@inline manifest_entry_available(repo_root::String, entry) =
    first(manifest_entry_readiness(repo_root, entry)) === :usable

function check_assets(; repo_root::String=REPO_ROOT,
    manifest_path::String=joinpath(repo_root, "data", "assets_manifest.toml"))
    items = NamedTuple[]
    for entry in load_asset_manifest(; repo_root, manifest_path)
        detail = isempty(entry.detail) ? "licensing=$(entry.licensing)" : "$(entry.detail) [licensing=$(entry.licensing)]"
        state, reason = manifest_entry_readiness(repo_root, entry)
        push!(items, (
            name=entry.name,
            scope=entry.scope,
            path=manifest_entry_path(repo_root, entry),
            available=state === :usable,
            required=entry.required,
            detail=detail,
            state=state,
            reason=reason,
        ))
    end
    return (; repo_root, items)
end

end # module AssetCatalog
