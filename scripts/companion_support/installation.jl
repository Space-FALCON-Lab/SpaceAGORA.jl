# Canonical source. Generate companion copies with scripts/sync_companion_support.jl.
# This module needs only standard libraries and works before core is installed.
module CompanionInstaller
using Pkg, TOML

function canonical_destination(directory)
    path = abspath(directory)
    tail = String[]
    while !ispath(path)
        islink(path) && throw(ArgumentError("Environment path contains a dangling symlink: $path"))
        pushfirst!(tail, basename(path))
        parent = dirname(path)
        parent == path && throw(ArgumentError("Cannot resolve environment destination"))
        path = parent
    end
    isdir(path) || throw(ArgumentError("Environment destination must be a directory"))
    return joinpath(realpath(path), tail...)
end
inside(path, root) = first(splitpath(relpath(path, root))) != ".."

# Called both while planning and immediately before setup, so a destination that
# changed after preflight cannot silently become an unmanaged installation target.
function validate_destination(destination, package_root, core, marker_name)
    destination = canonical_destination(destination)
    inside(destination, realpath(package_root)) && throw(ArgumentError("Choose an environment outside the companion checkout"))
    core.path === nothing || !inside(destination, realpath(core.path)) ||
        throw(ArgumentError("Choose an environment directory separate from the core checkout"))
    if isdir(destination) && !isempty(readdir(destination))
        marker = joinpath(destination, marker_name)
        islink(marker) && throw(ArgumentError("Refusing a symlinked environment ownership record: $marker"))
        isfile(marker) || throw(ArgumentError("Refusing to alter an existing unmanaged environment: $destination"))
        saved = TOML.parsefile(marker)
        identity = core.path === nothing ? core.revision : core.path
        get(saved, "companion", "") == package_root && get(saved, "core", "") == identity ||
            throw(ArgumentError("Existing environment belongs to a different source pair"))
    end
    return destination
end

function installation_plan(package_root::AbstractString, directory::AbstractString;
        pin_file::AbstractString, marker_name::AbstractString,
        core_path::AbstractString=get(ENV, "SPACEAGORA_CORE_PATH", ""),
        core_revision::AbstractString=get(ENV, "SPACEAGORA_CORE_REV", ""))
    project = TOML.parsefile(joinpath(package_root, "Project.toml"))
    destination = canonical_destination(directory)
    isempty(core_path) || isempty(core_revision) || throw(ArgumentError("Set only one of SPACEAGORA_CORE_PATH and SPACEAGORA_CORE_REV"))
    pin = TOML.parsefile(joinpath(package_root, pin_file))
    get(pin, "url", "") == "https://github.com/Space-FALCON-Lab/SpaceAGORA.jl.git" ||
        throw(ArgumentError("Unexpected core source repository"))
    # A standalone companion may integrate with core through an optional extension.
    # Validate both declarations without making optional packages direct dependencies.
    core_uuids = [section["SpaceAGORA"] for section in
        (get(project, "deps", Dict()), get(project, "weakdeps", Dict()))
        if haskey(section, "SpaceAGORA")]
    !isempty(core_uuids) && all(==("afbfb69f-5c0b-4832-b760-43725dff8540"), core_uuids) ||
        throw(ArgumentError("Declare the SpaceAGORA UUID consistently in deps or weakdeps"))
    core_uuid = first(core_uuids)
    get(pin, "uuid", "") == core_uuid ||
        throw(ArgumentError("Core pin UUID does not match the companion dependency"))
    if isempty(core_path)
        revision = lowercase(isempty(core_revision) ? get(pin, "rev", "") : String(core_revision))
        occursin(r"^[0-9a-fA-F]{40}$", revision) || throw(ArgumentError("Record a full core commit SHA in the companion pin file before standalone installation; development can explicitly set SPACEAGORA_CORE_PATH"))
        core = (; path=nothing, url=pin["url"], revision)
    else
        path = abspath(core_path)
        isfile(joinpath(path, "Project.toml")) || throw(ArgumentError("SPACEAGORA_CORE_PATH must contain a Project.toml"))
        path = realpath(path)
        core_project = TOML.parsefile(joinpath(path, "Project.toml"))
        get(core_project, "name", "") == "SpaceAGORA" &&
            get(core_project, "uuid", "") == core_uuid ||
            throw(ArgumentError("SPACEAGORA_CORE_PATH is not the SpaceAGORA package"))
        core = (; path, url=nothing, revision=nothing)
    end
    validate_destination(destination, package_root, core, marker_name)
    # Direct test and benchmark imports must resolve from the independent project,
    # with JULIA_LOAD_PATH=@:@stdlib and without an ambient developer environment.
    dependencies = merge(get(project, "deps", Dict{String,Any}()), get(project, "extras", Dict{String,Any}()))
    delete!(dependencies, "SpaceAGORA")
    return (; destination, core, dependencies, companion=package_root)
end

# Use a working-tree checkout rather than Pkg's bare LibGit2 extraction, which
# cannot materialize this core revision's uninitialized submodule entries.
function source_git(args...)
    command = Cmd(["git", "-c", "submodule.recurse=false", "-c", "filter.lfs.process=",
        "-c", "filter.lfs.smudge=", "-c", "filter.lfs.clean=", "-c", "filter.lfs.required=false",
        String.(args)...])
    return addenv(command, "GIT_LFS_SKIP_SMUDGE" => "1", "GIT_DIR" => nothing,
        "GIT_WORK_TREE" => nothing, "GIT_INDEX_FILE" => nothing, "GIT_COMMON_DIR" => nothing,
        "GIT_OBJECT_DIRECTORY" => nothing, "GIT_ALTERNATE_OBJECT_DIRECTORIES" => nothing)
end

function verify_core_checkout(path, core)
    islink(path) && throw(ArgumentError("Refusing a symlinked managed core checkout: $path"))
    gitdir = joinpath(path, ".git")
    isdir(gitdir) && !islink(gitdir) || throw(ArgumentError("Not a setup-owned core checkout: $path"))
    marker = joinpath(gitdir, "spaceagora-companion-source.toml")
    isfile(marker) || throw(ArgumentError("Missing setup ownership record in $path"))
    owned = TOML.parsefile(marker)
    get(owned, "url", nothing) == core.url && get(owned, "revision", nothing) == core.revision ||
        throw(ArgumentError("Managed core source identity does not match the requested pin: $path"))
    readchomp(source_git("-C", path, "remote", "get-url", "origin")) == core.url ||
        throw(ArgumentError("Managed core origin does not match the requested source: $path"))
    head = readchomp(source_git("-C", path, "rev-parse", "HEAD"))
    head == core.revision || throw(ArgumentError("Managed core HEAD $head does not match $(core.revision): $path"))
    clean = readchomp(source_git("-C", path, "status", "--porcelain=v1", "--untracked-files=all", "--ignore-submodules=none"))
    isempty(clean) || throw(ArgumentError("Managed core checkout has local changes; refusing to overwrite them: $path"))
    modules = joinpath(gitdir, "modules")
    !isdir(modules) || isempty(readdir(modules)) ||
        throw(ArgumentError("Managed core checkout must not have initialized submodules: $path"))
    project = TOML.parsefile(joinpath(path, "Project.toml"))
    project["name"] == "SpaceAGORA" && project["uuid"] == "afbfb69f-5c0b-4832-b760-43725dff8540" ||
        throw(ArgumentError("Pinned checkout is not the SpaceAGORA package: $path"))
    return path
end

function acquire_core(plan)
    plan.core.path === nothing || return plan.core.path  # the explicit core checkout is never mutated
    sources = joinpath(plan.destination, "sources")
    islink(sources) && throw(ArgumentError("Refusing a symlinked managed sources directory: $sources"))
    mkpath(sources)
    path = joinpath(sources, "SpaceAGORA-" * plan.core.revision)
    if ispath(path) || islink(path)
        return verify_core_checkout(path, plan.core)
    end
    # Build under a temporary sibling, then publish only a verified checkout.
    # A failed fetch or validation removes only this newly created temporary tree.
    mktempdir(sources; prefix=".core-") do temporary
        run(source_git("init", "--quiet", temporary))
        run(source_git("-C", temporary, "remote", "add", "origin", plan.core.url))
        run(source_git("-C", temporary, "fetch", "--quiet", "--no-tags", "--no-recurse-submodules", "--depth=1", "origin", plan.core.revision))
        run(source_git("-C", temporary, "checkout", "--quiet", "--no-recurse-submodules", "--detach", plan.core.revision))
        open(joinpath(temporary, ".git", "spaceagora-companion-source.toml"), "w") do io
            TOML.print(io, Dict("url" => plan.core.url, "revision" => plan.core.revision))
        end
        verify_core_checkout(temporary, plan.core)
        mv(temporary, path)
    end
    return path
end

function setup(plan; marker_name::AbstractString, package_name::AbstractString, precompile::Bool=true)
    destination = validate_destination(plan.destination, plan.companion, plan.core, marker_name)
    destination == plan.destination || throw(ArgumentError("Environment destination changed after planning"))
    mkpath(destination)
    # Record ownership before acquisition, so an interrupted first fetch is safely
    # retryable. Validation above refuses any existing unmanaged content.
    open(joinpath(destination, marker_name), "w") do io
        TOML.print(io, Dict("companion" => plan.companion,
            "core" => (plan.core.path === nothing ? plan.core.revision : plan.core.path)))
    end
    core_path = acquire_core(plan)  # verify source before any Pkg environment mutation
    previous = Base.active_project()
    try
        withenv("JULIA_PKG_PRECOMPILE_AUTO" => "0") do
            Pkg.activate(plan.destination)
            # Explicit core installation is intentional: dependency [sources]
            # propagation is not assumed, including on Julia 1.12.
            Pkg.develop(Pkg.PackageSpec(path=core_path))
            Pkg.develop(Pkg.PackageSpec(path=plan.companion); preserve=Pkg.PRESERVE_ALL)
            direct = [Pkg.PackageSpec(name=name, uuid=uuid) for (name, uuid) in sort!(collect(plan.dependencies))]
            Pkg.add(direct; preserve=Pkg.PRESERVE_ALL)
            Pkg.instantiate()
            precompile && Pkg.precompile(["SpaceAGORA", package_name]; strict=true)
        end
    finally
        previous === nothing ? Pkg.activate() : Pkg.activate(dirname(previous))
    end
    return plan.destination
end


end # module CompanionInstaller
