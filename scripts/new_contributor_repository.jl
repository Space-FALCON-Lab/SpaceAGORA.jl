"""Create a contributor repository from a starter in `templates/starters/`.

    julia --startup-file=no scripts/new_contributor_repository.jl <project|package|companion> <directory> --name <Name> [--holder <copyright holder>]

- `project`: a research study on the public API (committed Manifest, scripts, results kept out).
- `package`: an algorithm package for an existing family (an RPO reference planner).
- `companion`: a companion that serves a core service facade (the external-propagation boundary).

The new repository pins this checkout's commit: its `SpaceAGORASource.toml` and its
workflow call name `HEAD`, its `[compat]` admits this core's release series, and a
companion requires the facade's current contract version. It gets a fresh UUID; a
project also gets the commit's Manifest as its committed Manifest. A package or a
companion gets the MIT license, for this year and `--holder` (the copyright holder);
without `--holder` its placeholder `Your name or organization` stays for the author to
replace. The starter's files
come from the commit itself, so an untracked or edited file in the checkout never
reaches the new repository. The command refuses a checkout with tracked changes (the
pinned commit would not be the code the starter was made from), an existing non-empty
directory, a directory inside this checkout, and a name that is not a package name or
that contains the starters' placeholder `StarterName`.
"""
module NewContributorRepository

using Dates, TOML, UUIDs

const ROOT = normpath(joinpath(@__DIR__, ".."))
const KINDS = ("project", "package", "companion")
const STARTERS = joinpath("templates", "starters")
# The placeholders of the starter templates. Each keeps its file valid Julia, TOML
# or YAML, so the templates parse and the repository's gates can scan them.
const NAME_TOKEN = "StarterName"
const UUID_TOKEN = "00000000-0000-0000-0000-000000000000"
const REV_TOKEN = "0"^40
const COMPAT_TOKEN = "SpaceAGORA = \"0.0\""
const PIN_VERSION_TOKEN = "version = \"0.0.0\""
const CONTRACT_TOKEN = "const CORE_CONTRACT_VERSION = v\"0.0.0\""
# The LICENSE line: the year and the copyright holder.
const HOLDER_PLACEHOLDER = "Your name or organization"
const LICENSE_TOKEN = "Copyright (c) 0000 $(HOLDER_PLACEHOLDER)"
const FACADE_SOURCE = joinpath("src", "orchestration", "simulation", "external_propagation_services.jl")
const HELP = "Usage: julia --startup-file=no scripts/new_contributor_repository.jl <$(join(KINDS, '|'))> <directory> --name <Name> [--holder <copyright holder>]"

fail(message) = throw(ArgumentError(message))

function git(root, args...)
    command = `git -c maintenance.auto=false -c gc.auto=0 -C $root $args`
    return readchomp(command)
end

"""The checkout's commit, version and facade contract version; refuses tracked changes unless `require_clean=false`."""
function core_identity(root=ROOT; require_clean::Bool=true)
    rev = git(root, "rev-parse", "HEAD")
    occursin(r"^[0-9a-f]{40}$", rev) || fail("$root: HEAD is not a commit ($rev).")
    if require_clean
        changes = git(root, "status", "--porcelain=v1", "--untracked-files=no", "--ignore-submodules=all")
        isempty(changes) || fail("$root has tracked changes, so its HEAD $rev is not the code a starter would be " *
            "made from. Commit or set them aside first:\n$changes")
    end
    project = TOML.parsefile(joinpath(root, "Project.toml"))
    get(project, "name", "") == "SpaceAGORA" || fail("$root is not the SpaceAGORA package.")
    version = VersionNumber(project["version"])
    source = read(joinpath(root, FACADE_SOURCE), String)
    contract = match(r"const CONTRACT_VERSION = v\"([^\"]+)\"", source)
    contract === nothing && fail("$(FACADE_SOURCE): no CONTRACT_VERSION found.")
    return (; root, rev, uuid=project["uuid"], version, contract=VersionNumber(contract.captures[1]))
end

"""The release series a `[compat]` entry admits: `0.minor` before 1.0, the major version after it."""
compat_series(version::VersionNumber) = version.major == 0 ? "0.$(version.minor)" : "$(version.major)"

"""Refuse a name that is not a package name, or one the starters reserve."""
function check_name(name::AbstractString)
    occursin(r"^[A-Z][A-Za-z0-9]*$", name) ||
        fail("--name must be a package name: a capital letter, then letters and digits (got \"$name\").")
    name == "SpaceAGORA" && fail("--name $name is reserved; choose the repository's own name.")
    # The placeholder must not survive in a generated file, so it cannot be part of a name.
    occursin(NAME_TOKEN, name) &&
        fail("--name must not contain \"$(NAME_TOKEN)\", the starters' placeholder for the name (got \"$name\").")
    endswith(name, "jl") && fail("--name is the package name without .jl (got \"$name\").")
    return String(name)
end

inside(path, root) = first(splitpath(relpath(path, root))) != ".."

function check_target(directory, root)
    target = abspath(directory)
    inside(target, realpath(root)) && fail("$target is inside the SpaceAGORA checkout; create the repository elsewhere.")
    isdir(target) && !isempty(readdir(target)) && fail("$target exists and is not empty.")
    ispath(target) && !isdir(target) && fail("$target exists and is not a directory.")
    return target
end

# A file's text at `HEAD`, as committed.
committed_text(root, path) = read(`git -c maintenance.auto=false -c gc.auto=0 -C $root show HEAD:$path`, String)

"""
    template_files(root, kind; committed=true) -> Vector{Pair{String, String}}

The `kind` starter's files, relative to its folder, with their text: from the commit
`HEAD` (`committed=true`), or from the working tree, untracked files included, for
work on the templates themselves (`committed=false`).
"""
function template_files(root, kind; committed::Bool=true)
    folder = joinpath(STARTERS, kind)
    if committed
        listed = git(root, "ls-tree", "-r", "-z", "--name-only", "HEAD", "--", folder)
        paths = filter(!isempty, split(listed, '\0'))
        isempty(paths) && fail("$(folder) is not in the commit HEAD.")
        return [relpath(path, folder) => committed_text(root, path) for path in paths]
    end
    template = joinpath(root, folder)
    isdir(template) || fail("$template is missing.")
    return [relpath(joinpath(dir, file), template) => read(joinpath(dir, file), String)
            for (dir, _, files) in walkdir(template) for file in files]
end

"""The text of a starter file with every placeholder filled in."""
function fill_in(text::AbstractString, values)
    for (token, value) in values
        text = replace(text, token => value)
    end
    return text
end

"""Refuse a copyright holder that is empty or spans lines."""
function check_holder(holder::AbstractString)
    stripped = strip(holder)
    isempty(stripped) && fail("--holder must name the copyright holder (got an empty name).")
    occursin('\n', stripped) && fail("--holder must be a single line (got \"$holder\").")
    return String(stripped)
end

"""
    generate(kind, directory; name, holder=nothing, year=this year, root=ROOT, require_clean=true) -> directory

Write the `kind` starter, filled in for `name` and this checkout, into `directory`. A
starter's `LICENSE` names `year` and `holder`, or the placeholder
`Your name or organization` when `holder` is `nothing`.
With `require_clean=false` (tests and work on the templates) the checkout may have
changes, and the starter comes from the working tree.
"""
function generate(kind::AbstractString, directory::AbstractString; name::AbstractString,
        holder::Union{Nothing, AbstractString}=nothing, year::Integer=Dates.year(Dates.today()), root=ROOT,
        require_clean::Bool=true)
    kind in KINDS || fail("Unknown kind \"$kind\". Choose from $(join(KINDS, ", ")).")
    name = check_name(name)
    holder = holder === nothing ? HOLDER_PLACEHOLDER : check_holder(holder)
    core = core_identity(root; require_clean)
    target = check_target(directory, root)
    files = template_files(root, kind; committed=require_clean)
    values = (
        NAME_TOKEN => name,
        UUID_TOKEN => string(uuid4()),
        REV_TOKEN => core.rev,
        COMPAT_TOKEN => "SpaceAGORA = \"$(compat_series(core.version))\"",
        PIN_VERSION_TOKEN => "version = \"$(core.version)\"",
        CONTRACT_TOKEN => "const CORE_CONTRACT_VERSION = v\"$(core.contract)\"",
        LICENSE_TOKEN => "Copyright (c) $(year) $(holder)",
    )
    mkpath(target)
    for (path, text) in files
        destination = joinpath(target, fill_in(path, (NAME_TOKEN => name,)))
        mkpath(dirname(destination))
        write(destination, fill_in(text, values))
    end
    # A study commits the exact versions behind its results; it starts from the core's.
    if kind == "project"
        manifest = require_clean ? committed_text(root, "Manifest.toml") : read(joinpath(root, "Manifest.toml"), String)
        write(joinpath(target, "Manifest.toml"), manifest)
    end
    for (folder, _, files) in walkdir(target), file in files
        text = read(joinpath(folder, file), String)
        for token in (NAME_TOKEN, UUID_TOKEN, REV_TOKEN)
            occursin(token, text) && error("$(joinpath(folder, file)) still contains the placeholder $token.")
        end
    end
    return target
end

function main(args=ARGS)
    positional, name, holder = String[], nothing, nothing
    k = 1
    while k <= length(args)
        if args[k] == "--name" && k < length(args)
            name = args[k+1]
            k += 2
        elseif args[k] == "--holder" && k < length(args)
            holder = args[k+1]
            k += 2
        elseif startswith(args[k], "--")
            fail("Unknown option $(args[k]).\n$HELP")
        else
            push!(positional, args[k])
            k += 1
        end
    end
    length(positional) == 2 && name !== nothing || fail(HELP)
    kind, directory = positional
    target = generate(kind, directory; name, holder)
    core = core_identity(ROOT; require_clean=false)
    println("Created the $kind starter $name in $target, pinned to SpaceAGORA $(core.rev).")
    if holder === nothing && isfile(joinpath(target, "LICENSE"))
        println("Its LICENSE (MIT) names \"$(HOLDER_PLACEHOLDER)\": put the copyright holder there, or pass --holder.")
    end
    println("Next: read its README.md, then `git init` there and push it to its own repository.")
    return target
end

end # module NewContributorRepository

if abspath(PROGRAM_FILE) == @__FILE__
    try
        NewContributorRepository.main()
    catch error
        error isa ArgumentError || rethrow()
        println(stderr, "ERROR: ", error.msg)
        exit(1)
    end
end
