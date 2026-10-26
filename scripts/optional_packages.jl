"""Locate and install optional packages from the domain folders without changing core."""
module OptionalPackages

using TOML
const ROOT = normpath(joinpath(@__DIR__, ".."))
const NAMES = ("hypr", "mujoco", "edg")
const DOMAINS = Dict("hypr" => "robotics", "mujoco" => "robotics", "edg" => "gnc")
const HELP = """Usage: julia --startup-file=no scripts/optional_packages.jl status
       julia --startup-file=no scripts/optional_packages.jl checkout <hypr|mujoco|edg>
       julia --startup-file=no scripts/optional_packages.jl check <hypr|mujoco|edg>
       julia --startup-file=no scripts/optional_packages.jl setup <hypr|mujoco|edg> <environment-directory>
Use a Git clone of SpaceAGORA; source ZIPs omit linked repository contents.
Setup environments must be outside the source checkout."""

# Ambient Git variables must not redirect checks or updates to another repository.
function git(root, args...)
    cmd = Cmd(["git", "-c", "submodule.recurse=false", "-C", root, String.(args)...])
    return addenv(cmd, "GIT_LFS_SKIP_SMUDGE" => "1", "GIT_OPTIONAL_LOCKS" => "0", "GIT_DIR" => nothing,
        "GIT_WORK_TREE" => nothing, "GIT_INDEX_FILE" => nothing, "GIT_COMMON_DIR" => nothing,
        "GIT_OBJECT_DIRECTORY" => nothing, "GIT_ALTERNATE_OBJECT_DIRECTORIES" => nothing)
end
gitread(root, args...) = readchomp(pipeline(git(root, args...), stderr=devnull))
fail(message) = throw(ArgumentError(message))
canonical_url(url) = endswith(url, ".git") ? chop(url; tail=4) : url
same_origin(a, b) = canonical_url(a) == canonical_url(b)

function source_root(root)
    path = realpath(root)
    ispath(joinpath(path, ".git")) ||
        fail("Use a Git clone of SpaceAGORA. Source ZIPs do not contain optional package pins or checkouts.")
    realpath(gitread(path, "rev-parse", "--show-toplevel")) == path ||
        fail("The selected directory must be the SpaceAGORA Git checkout root.")
    return path
end

function no_symlinks(root, relative)
    current = root
    for part in splitpath(relative)
        current = joinpath(current, part)
        islink(current) && fail("Refusing a symlinked optional package path: $current")
    end
end

"""Read the index gitlink pin and check the accompanying installation metadata."""
function spec(name; root=ROOT)
    name in NAMES || fail("Unknown optional package '$name'. Choose hypr, mujoco or edg.")
    root = source_root(root)
    item = TOML.parsefile(joinpath(root, "packages", "OptionalPackages.toml"))[name]
    relative = item["path"]
    domain = joinpath("packages", DOMAINS[name])
    length(splitpath(relative)) == 3 && startswith(relative, domain * "/") &&
        normpath(relative) == relative && !isabspath(relative) ||
        fail("The $name optional package path must stay directly under $domain.")
    no_symlinks(root, relative)
    fields = split(gitread(root, "ls-files", "--stage", "--", relative))
    length(fields) == 4 && fields[1] == "160000" && fields[3] == "0" &&
        fields[4] == relative && occursin(r"^[0-9a-f]{40}$", fields[2]) ||
        fail("Missing or conflicted gitlink for $relative. Installation requires an immutable published package pin; use a complete, qualified Git clone.")
    revision = fields[2]
    for (key, expected) in (("path", relative), ("url", item["url"]))
        actual = gitread(root, "config", "--file", ".gitmodules", "--get", "submodule.$(item["name"]).$key")
        actual == expected || fail("The $name .gitmodules $key does not match its package metadata.")
    end
    if name == "hypr"
        pin = TOML.parsefile(joinpath(root, "packages", "HYPRSource.toml"))
        all(get(pin, key, nothing) == item[key] for key in ("url", "uuid", "version")) &&
            pin["rev"] == revision || fail("HYPRSource.toml and the HYPR gitlink must pin the same package revision.")
    end
    return (; name=item["name"], path=relative, uuid=item["uuid"], version=item["version"],
        url=item["url"], revision)
end

function check_source(root, package)
    path = joinpath(root, package.path)
    no_symlinks(root, package.path)
    isdir(path) && ispath(joinpath(path, ".git")) && !islink(joinpath(path, ".git")) ||
        fail("Optional package is not initialized at $path. Run checkout first.")
    realpath(gitread(path, "rev-parse", "--show-toplevel")) == realpath(path) ||
        fail("Optional package must have its own Git checkout at $path.")
    same_origin(gitread(path, "remote", "get-url", "origin"), package.url) ||
        fail("Optional package origin differs from the declared repository at $path.")
    head = gitread(path, "rev-parse", "HEAD")
    head == package.revision ||
        fail("Optional package HEAD $head differs from pin $(package.revision). Preserve your work and reconcile it manually.")
    clean = gitread(path, "status", "--porcelain=v1", "--untracked-files=all", "--ignore-submodules=none")
    isempty(clean) || fail("Optional package has local changes at $path. Setup will not overwrite them.")
    project = TOML.parsefile(joinpath(path, "Project.toml"))
    all(get(project, key, nothing) == getproperty(package, Symbol(key)) for key in ("name", "uuid", "version")) ||
        fail("Optional package name, UUID or version does not match the declared package at $path.")
    return path
end

"""Verify source identity and cleanliness; this does not certify runtime compatibility."""
function check(name; root=ROOT)
    package = spec(name; root)
    return check_source(realpath(root), package)
end

"""Report source availability without initialization, network requests or file writes."""
function status(name; root=ROOT)
    try
        package = spec(name; root)
        path = joinpath(root, package.path)
        if !ispath(path) || (isdir(path) && isempty(readdir(path)))
            return (; state=:absent, message="$(package.path) @ $(package.revision) (not initialized)")
        end
        check_source(realpath(root), package)
        return (; state=:ready, message="$(package.path) @ $(package.revision)")
    catch error
        return (; state=:invalid, message=sprint(showerror, error))
    end
end

"""Initialize only the selected pinned submodule, never reset an existing checkout."""
function checkout(name; root=ROOT)
    package = spec(name; root)
    root = realpath(root)
    path = joinpath(root, package.path)
    if ispath(path) && !(isdir(path) && isempty(readdir(path)))
        return check_source(root, package)
    end
    # A configured URL override could fetch a different repository. Refuse it.
    command = git(root, "config", "--get", "submodule.$(package.name).url")
    if success(pipeline(command, stdout=devnull, stderr=devnull))
        same_origin(gitread(root, "config", "--get", "submodule.$(package.name).url"), package.url) ||
            fail("The configured submodule URL differs from the declared $name repository.")
    end
    run(git(root, "submodule", "update", "--init", "--checkout", "--depth=1", "--", package.path))
    return check_source(root, package)
end

# Resolve existing ancestors too, so a not-yet-created environment under a symlink
# cannot bypass the source-directory refusal.
function canonical_destination(directory)
    path = abspath(directory)
    tail = String[]
    while !ispath(path)
        islink(path) && fail("Environment path contains a dangling symlink: $path")
        pushfirst!(tail, basename(path))
        parent = dirname(path)
        parent == path && fail("Cannot resolve environment destination.")
        path = parent
    end
    isdir(path) || fail("Environment destination must be a directory.")
    return joinpath(realpath(path), tail...)
end

function destination(directory; root=ROOT)
    root = realpath(root)
    path = canonical_destination(directory)
    relative = relpath(path, root)
    (first(splitpath(relative)) == "..") ||
        fail("Choose an environment directory outside the SpaceAGORA source checkout.")
    return path
end

"""Build the existing installer's command after all source/destination preflight."""
function installation_command(name, environment; root=ROOT)
    root = source_root(root)
    env = destination(environment; root)
    path = check(name; root)
    if name == "hypr"
        isfile(joinpath(root, "Manifest.toml")) || fail("HYPR setup requires the core Manifest.toml.")
        driver = "include(ARGS[1]); HYPRInstallation.setup(ARGS[2]; root=ARGS[3], local_path=ARGS[4])"
        script = joinpath(root, "scripts", "setup_hypr.jl")
        return addenv(`$(Base.julia_cmd()) --startup-file=no --project=@stdlib -e $driver $script $env $root $path`,
            "JULIA_LOAD_PATH" => "@:@stdlib")
    end
    script = joinpath(path, "scripts", "setup_env.jl")
    isfile(script) || fail("Missing $(name) companion installer.")
    return addenv(`$(Base.julia_cmd()) --startup-file=no --project=@stdlib $script $env`,
        "SPACEAGORA_CORE_PATH" => root, "SPACEAGORA_CORE_REV" => "", "JULIA_LOAD_PATH" => "@:@stdlib")
end

function setup(name, environment; root=ROOT)
    destination(environment; root) # fail before fetching a package
    name == "hypr" && !isfile(joinpath(root, "Manifest.toml")) &&
        fail("HYPR setup requires the core Manifest.toml.")
    checkout(name; root)
    run(installation_command(name, environment; root))
    return canonical_destination(environment)
end

function main(args=ARGS; root=ROOT)
    if args == ["status"]
        for name in NAMES
            result = status(name; root)
            println(name, ": ", result.state, " | ", result.message)
        end
        return nothing
    elseif length(args) == 2 && args[1] in ("checkout", "check")
        path = args[1] == "checkout" ? checkout(args[2]; root) : check(args[2]; root)
        println("Verified source: ", path)
    elseif length(args) == 3 && args[1] == "setup"
        println("Environment: ", setup(args[2], args[3]; root))
    else
        fail(HELP)
    end
end

end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && OptionalPackages.main()
