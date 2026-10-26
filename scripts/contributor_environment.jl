"""Build the separate environment in which a contributor repository is tested against this core.

The repository records its core commit in `SpaceAGORASource.toml`, beside its
`Project.toml`, not in `[sources]`: Pkg cannot install SpaceAGORA from a url and rev
(the core has submodules), and `Pkg.test` would read such an entry instead of the
environment's core. The environment develops this checkout by path, so the tested
core is exactly the checked-out commit, and the pin file must name that commit. The
reusable workflow `.github/workflows/contributor-ci.yml` calls this script; it also
works from a local clone of the pinned commit, after which `Pkg.test` runs in the
environment as usual. `repin`, run from a checkout of another commit, moves a
repository's pin to that commit.

    # SpaceAGORASource.toml
    url = "https://github.com/Space-FALCON-Lab/SpaceAGORA.jl.git"
    rev = "<40-character commit SHA>"
    uuid = "afbfb69f-5c0b-4832-b760-43725dff8540"
    version = "0.6.1"          # optional: the core version at rev
    companions = ["mujoco"]    # optional: hypr, mujoco, edg, at this commit's pins
"""
module ContributorEnvironment

using Pkg, TOML
include(joinpath(@__DIR__, "optional_packages.jl"))
using .OptionalPackages: OptionalPackages, gitread, fail, same_origin

const ROOT = normpath(joinpath(@__DIR__, ".."))
const CORE_URL = "https://github.com/Space-FALCON-Lab/SpaceAGORA.jl.git"
const KINDS = ("project", "package", "companion")
const MANIFEST_MODES = ("core", "committed", "resolve")
const SCANNED_KINDS = ("project", "package")
const INTERNAL_MODULES = (:SimulationModel, :SimulationEngine)
const RECORD = "ContributorEnvironment.toml"
const PIN_FILE = "SpaceAGORASource.toml"
const PIN_KEYS = ("url", "rev", "uuid", "version", "companions", "path")
const HELP = """Usage: julia --startup-file=no scripts/contributor_environment.jl <project|package|companion> <repository> <environment> --manifest <core|committed|resolve> [--rev <sha>]
       julia --startup-file=no scripts/contributor_environment.jl check <project|package|companion> <repository> [--rev <sha>]
       julia --startup-file=no scripts/contributor_environment.jl scan <project|package|companion> <repository>
       julia --startup-file=no scripts/contributor_environment.jl summary <report-directory> [<environment>]
       julia --startup-file=no scripts/contributor_environment.jl repin <project|package|companion> <repository> [--rev <sha>]
The environment directory must be new and outside both the repository and this checkout.
--rev defaults to this checkout's HEAD; the repository's $PIN_FILE must pin SpaceAGORA to it.
repin points the repository's pin file, workflow calls and [compat] at this checkout's commit.
--core-from-caller (core CI only) also accepts a pin file whose path is this checkout."""

full_sha(rev) = rev isa AbstractString && occursin(r"^[0-9a-f]{40}$", rev)
inside(path, root) = first(splitpath(relpath(path, root))) != ".."

function project_file(directory)
    for name in Base.project_names
        path = joinpath(directory, name)
        isfile(path) && return path
    end
    return nothing
end

function manifest_file(directory)
    versioned = joinpath(directory, "Manifest-v$(VERSION.major).$(VERSION.minor).toml")
    isfile(versioned) && return versioned
    plain = joinpath(directory, "Manifest.toml")
    return isfile(plain) ? plain : nothing
end

"""The core checkout at `rev`, with no tracked changes: the commit that is tested."""
function verified_core(root, rev)
    path = OptionalPackages.source_root(root)
    head = gitread(path, "rev-parse", "HEAD")
    rev = rev === nothing ? head : String(rev)
    full_sha(rev) || fail("The SpaceAGORA revision must be a full 40-character lowercase commit SHA, not '$rev'.")
    head == rev || fail("The SpaceAGORA checkout at $path is at $head, not the requested revision $rev.")
    # Submodules are not part of the tested core; optional_packages.jl verifies companions.
    changes = gitread(path, "status", "--porcelain=v1", "--untracked-files=no", "--ignore-submodules=all")
    isempty(changes) || fail("The SpaceAGORA checkout at $path has tracked changes, so it is not revision $rev:\n$changes")
    project = TOML.parsefile(joinpath(path, "Project.toml"))
    get(project, "name", "") == "SpaceAGORA" || fail("$path is not the SpaceAGORA package.")
    return (; path, rev, uuid=project["uuid"], version=get(project, "version", ""))
end

refuse_core_source(project, file) = haskey(get(project, "sources", Dict{String,Any}()), "SpaceAGORA") &&
    fail("$file: remove the SpaceAGORA entry from [sources]. Pkg cannot install SpaceAGORA from a url and rev " *
         "(the core has submodules), and Pkg.test would use that entry instead of the tested core. " *
         "The commit goes in $PIN_FILE.")

"""The repository's pin file, checked against the tested core: `(; description, companions)`."""
function core_pin(project, file, repository, core; core_from_caller=false)
    refuse_core_source(project, file)
    pinfile = joinpath(repository, PIN_FILE)
    isfile(pinfile) || fail("$pinfile is missing. It records the SpaceAGORA commit this repository is tested against:\n" *
        "url = \"$CORE_URL\"\nrev = \"$(core.rev)\"\nuuid = \"$(core.uuid)\"")
    pin = TOML.parsefile(pinfile)
    unknown = setdiff(keys(pin), PIN_KEYS)
    isempty(unknown) || fail("$pinfile: unknown key(s) $(join(sort!(collect(unknown)), ", ")); allowed: $(join(PIN_KEYS, ", ")).")
    get(pin, "uuid", nothing) == core.uuid || fail("$pinfile: uuid must be SpaceAGORA's \"$(core.uuid)\".")
    version = get(pin, "version", nothing)
    version === nothing || version == core.version ||
        fail("$pinfile: version $(repr(version)) is not $(core.version), the SpaceAGORA version at $(core.rev).")
    companions = get(pin, "companions", String[])
    companions isa AbstractVector && all(c -> c isa AbstractString, companions) ||
        fail("$pinfile: companions must be a list of names, e.g. companions = [\"mujoco\"].")
    if haskey(pin, "path")
        core_from_caller || fail("$pinfile: a path pin is for core CI only; record the commit with url and rev.")
        (haskey(pin, "url") || haskey(pin, "rev")) && fail("$pinfile: give either path or url and rev, not both.")
        target = normpath(joinpath(repository, pin["path"]))
        isdir(target) && samefile(target, core.path) ||
            fail("$pinfile: path $(pin["path"]) is not the core under test ($(core.path)).")
        return (; description="path $(pin["path"])", companions)
    end
    url = get(pin, "url", nothing)
    rev = get(pin, "rev", nothing)
    url isa AbstractString && same_origin(url, CORE_URL) || fail("$pinfile: url must be $CORE_URL, not $(repr(url)).")
    rev isa AbstractString || fail("$pinfile: rev is missing; a branch is not a pin. Use rev = \"$(core.rev)\".")
    full_sha(rev) || fail("$pinfile: rev \"$rev\" is not a full 40-character lowercase commit SHA.")
    rev == core.rev || fail("$pinfile: pins SpaceAGORA $rev, but the tested revision is $(core.rev). " *
        "Use the same commit in $PIN_FILE and in the workflow call, so the recorded pair is the tested pair.")
    return (; description="rev $rev", companions)
end

"""Name and UUID of the repository's own package, or `nothing` for a project without a module."""
function package_identity(kind, project, file, repository)
    name, uuid = get(project, "name", nothing), get(project, "uuid", nothing)
    has_module = name isa AbstractString && isfile(joinpath(repository, "src", name * ".jl"))
    if kind != "project"
        name isa AbstractString && uuid isa AbstractString || fail("$file: a $kind needs a name and a uuid.")
        has_module || fail("$(joinpath(repository, "src", string(name, ".jl"))): the $kind's module file is missing.")
    end
    name == "SpaceAGORA" && fail("$file: this is the SpaceAGORA package itself, not a contributor repository.")
    return has_module && uuid isa AbstractString ? (; name=String(name), uuid=String(uuid)) : nothing
end

"""The repository's `test/Project.toml` as `(data, file)`; empty data when there is none."""
function test_project(repository)
    file = project_file(joinpath(repository, "test"))
    file === nothing && return (Dict{String,Any}(), nothing)
    data = TOML.parsefile(file)
    refuse_core_source(data, file)
    return (data, file)
end

"""Dependencies that `test/runtests.jl` imports: `[targets] test`, and a `test/Project.toml`'s `[deps]`."""
function test_dependencies(project, file, tests)
    deps = Dict{String,Any}()
    for name in get(get(project, "targets", Dict{String,Any}()), "test", String[])
        uuid = get(get(project, "extras", Dict{String,Any}()), name,
                   get(get(project, "weakdeps", Dict{String,Any}()), name, nothing))
        uuid === nothing && fail("$file: test target $name is not in [extras] or [weakdeps].")
        deps[name] = uuid
    end
    merge!(deps, get(tests, "deps", Dict{String,Any}()))
    delete!(deps, get(project, "name", ""))
    return deps
end

"""A compat entry for the versions both bounds admit; refused when they admit none."""
function compat_intersection(name, a, b, file)
    a == b && return a
    spec = intersect(Pkg.Versions.semver_spec(a), Pkg.Versions.semver_spec(b))
    isempty(spec) && fail("$file: [compat] $name = \"$b\" does not overlap the bound \"$a\" of the repository's Project.toml.")
    bound(v) = v.n == 0 ? "0" : join(Int.(v.t[1:v.n]), ".")
    # Each range as "lower - upper": a single version would read as a caret bound.
    return join((r.upper.n == 0 ? ">= $(bound(r.lower))" : "$(bound(r.lower)) - $(bound(r.upper))" for r in spec.ranges), ", ")
end

"""Bounds for the environment's direct dependencies: the repository's and its test project's, intersected."""
function environment_compat(project, tests, testfile, deps)
    compat = Dict{String,Any}()
    for (data, file) in ((project, nothing), (tests, testfile))
        for (name, bound) in get(data, "compat", Dict{String,Any}())
            name == "julia" || haskey(deps, name) || continue
            compat[name] = haskey(compat, name) ? compat_intersection(name, compat[name], bound, file) : bound
        end
    end
    return compat
end

# A [sources] entry with a path written relative to `directory`, made absolute.
rerooted(entry, directory) = haskey(entry, "path") ? merge(entry, Dict("path" => normpath(joinpath(directory, entry["path"])))) : entry

companion_names(list::AbstractVector) = companion_names(join(list, ' '))
function companion_names(list::AbstractString)
    names = unique(filter!(!isempty, split(list, r"[\s,]+")))
    for name in names
        name in OptionalPackages.NAMES ||
            fail("Unknown companion '$name'. Choose from $(join(OptionalPackages.NAMES, ", ")).")
    end
    return String.(names)
end

"""Package names that each companion installs from its recorded source pin."""
function companion_package_names(name; root)
    registry = TOML.parsefile(joinpath(root, "packages", "OptionalPackages.toml"))
    return [registry[name]["name"]]
end

"""A dependency on a companion needs the companion listed in the pin file, at the core's recorded pin."""
function check_companion_dependencies(project, file, companions; root, tests=Dict{String,Any}())
    used = merge(get(project, "deps", Dict{String,Any}()), get(project, "extras", Dict{String,Any}()),
                 get(tests, "deps", Dict{String,Any}()))
    sources = merge(get(tests, "sources", Dict{String,Any}()), get(project, "sources", Dict{String,Any}()))
    for name in OptionalPackages.NAMES, dependency in companion_package_names(name; root)
        haskey(used, dependency) || continue
        name in companions || fail("$file: $dependency is a dependency; list \"$name\" in $PIN_FILE companions.")
        rev = get(get(sources, dependency, Dict{String,Any}()), "rev", nothing)
        rev === nothing && continue
        pinned = OptionalPackages.spec(name; root).revision
        rev == pinned || fail("$file: [sources] pins $dependency $rev, but this core pins $pinned.")
    end
end

"""Check the repository against the tested core: its pin file, package identity and companions."""
function check_repository(kind, repository; rev=nothing, root=ROOT, core_from_caller=false)
    kind in KINDS || fail("Unknown kind '$kind'. Choose from $(join(KINDS, ", ")).")
    core = verified_core(root, rev)
    isdir(repository) || fail("Repository directory not found: $repository")
    repository = realpath(repository)
    file = project_file(repository)
    file === nothing && fail("$repository has no Project.toml.")
    project = TOML.parsefile(file)
    package = package_identity(kind, project, file, repository)
    pin = core_pin(project, file, repository, core; core_from_caller)
    tests, testfile = test_project(repository)
    names = companion_names(pin.companions)
    check_companion_dependencies(project, file, names; root=core.path, tests)
    direct = get(get(project, "deps", Dict{String,Any}()), "SpaceAGORA", nothing)
    weak = get(get(project, "weakdeps", Dict{String,Any}()), "SpaceAGORA", nothing)
    if kind == "companion"
        # Only the separate integration environment makes an optional core
        # adapter's weak dependency mandatory; the package stays standalone.
        (direct == core.uuid || weak == core.uuid) ||
            fail("$file: [deps] or [weakdeps] must list SpaceAGORA = \"$(core.uuid)\".")
    else
        direct == core.uuid || fail("$file: [deps] must list SpaceAGORA = \"$(core.uuid)\".")
    end
    # A correct dependency must not hide a conflicting declaration in another
    # table or the test project when we construct the integration environment.
    for (data, source) in ((project, file), (tests, testfile)), table in ("deps", "weakdeps", "extras")
        declared = get(get(data, table, Dict{String,Any}()), "SpaceAGORA", core.uuid)
        declared == core.uuid || fail("$source: SpaceAGORA in [$table] has a conflicting UUID; expected $(core.uuid).")
    end
    return (; kind, core, repository, file, project, tests, testfile, package, pin=pin.description, companions=names)
end

"""Check every input and decide the environment's contents; nothing is written or fetched."""
function plan(kind, repository, environment; manifest="core", rev=nothing, root=ROOT, core_from_caller=false)
    manifest in MANIFEST_MODES || fail("Unknown manifest mode '$manifest'. Choose from $(join(MANIFEST_MODES, ", ")).")
    (; core, repository, file, project, tests, testfile, package, pin, companions) =
        check_repository(kind, repository; rev, root, core_from_caller)
    destination = OptionalPackages.canonical_destination(environment)
    inside(destination, core.path) && fail("Choose an environment directory outside the SpaceAGORA checkout.")
    inside(destination, repository) && fail("Choose an environment directory outside the repository.")
    isdir(destination) && !isempty(readdir(destination)) && fail("Choose a new or empty environment directory: $destination")
    seed = if manifest == "core"
        joinpath(core.path, "Manifest.toml")
    elseif manifest == "committed"
        committed = manifest_file(repository)
        committed === nothing && fail("$repository has no committed Manifest.toml for manifest mode 'committed'.")
        committed
    else
        nothing
    end
    deps = merge(get(project, "deps", Dict{String,Any}()), test_dependencies(project, file, tests))
    deps["SpaceAGORA"] = core.uuid
    sources = Dict{String,Any}("SpaceAGORA" => Dict("path" => core.path))
    developed = [core.path]
    if package !== nothing
        deps[package.name] = package.uuid
        sources[package.name] = Dict("path" => repository)
        push!(developed, repository)
    end
    # Every direct dependency keeps its source: listed in the environment, it no longer
    # takes one from the developed package. The test project's entry (relative to test/)
    # wins for a package, as under Pkg.test; a project's own entry wins, since the
    # project is the environment its tests run in.
    own = Dict(name => rerooted(entry, repository) for (name, entry) in get(project, "sources", Dict{String,Any}()))
    test = Dict(name => rerooted(entry, joinpath(repository, "test")) for (name, entry) in get(tests, "sources", Dict{String,Any}()))
    for (name, entry) in (kind == "project" ? merge(test, own) : merge(own, test))
        haskey(deps, name) && !haskey(sources, name) && (sources[name] = entry)
    end
    compat = environment_compat(project, tests, testfile, deps)
    return (; kind, repository, file, package, pin, core, destination, manifest, seed, companions,
              project=Dict{String,Any}("deps" => deps, "sources" => sources, "compat" => compat), developed)
end

"""Fetch the selected companions at the core's gitlinks and add them to the plan."""
function add_companions(p)
    linked = Dict{String,Any}()
    for name in p.companions
        package = OptionalPackages.spec(name; root=p.core.path)
        paths = [OptionalPackages.checkout(name; root=p.core.path)]
        for path in paths
            identity = TOML.parsefile(joinpath(path, "Project.toml"))
            p.project["deps"][identity["name"]] = identity["uuid"]
            p.project["sources"][identity["name"]] = Dict("path" => path)
            push!(p.developed, path)
        end
        linked[name] = Dict("package" => package.name, "revision" => package.revision, "path" => first(paths))
    end
    return linked
end

function write_toml(path, data)
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
end

"""Write, seed and instantiate the environment, then record the tested pair."""
function build(p)
    linked = add_companions(p)
    mkpath(p.destination)
    write_toml(joinpath(p.destination, "Project.toml"), p.project)
    p.seed === nothing || cp(p.seed, joinpath(p.destination, "Manifest.toml"))
    preserve = p.seed === nothing ? Pkg.PRESERVE_NONE : Pkg.PRESERVE_ALL
    previous = Base.active_project()
    try
        withenv("JULIA_PKG_PRECOMPILE_AUTO" => "0") do
            Pkg.activate(p.destination)
            try
                Pkg.develop([PackageSpec(path=path) for path in p.developed]; preserve)
            catch error
                error isa Pkg.Resolve.ResolverError && p.seed !== nothing &&
                    @error "The repository needs versions that the $(p.manifest) manifest does not pin. Manifest mode 'resolve' shows the newest compatible set."
                rethrow()
            end
            Pkg.instantiate()
        end
    finally
        previous === nothing ? Pkg.activate() : Pkg.activate(dirname(previous))
    end
    write_toml(joinpath(p.destination, RECORD), Dict(
        "kind" => p.kind, "repository" => p.repository, "spaceagora_rev" => p.core.rev,
        "spaceagora_pin" => p.pin, "core" => p.core.path, "manifest" => p.manifest,
        "package" => p.package === nothing ? "" : p.package.name, "companions" => linked))
    return p.destination
end

setup(kind, repository, environment; kw...) = build(plan(kind, repository, environment; kw...))

# --- moving the pin ------------------------------------------------------------------------------

const WORKFLOW_CALL = r"Space-FALCON-Lab/SpaceAGORA\.jl/\.github/workflows/contributor-ci\.yml@([0-9a-f]{40})"

# `text` with the string value of `key = "..."` replaced by `value`, on the first such line
# at the top level (`table=nothing`) or in `[table]`; `nothing` when there is no such line.
function replace_toml_string(text::AbstractString, key::AbstractString, value::AbstractString; table=nothing)
    lines = split(text, '\n')
    pattern = Regex("^(\\s*\\Q$(key)\\E\\s*=\\s*)\"[^\"]*\"(.*)\$")
    current = nothing
    for (k, line) in enumerate(lines)
        header = match(r"^\s*\[\s*([^\[\]]+?)\s*\]\s*(#.*)?$", line)
        if header !== nothing
            current = header.captures[1]
            continue
        end
        current == table || continue
        m = match(pattern, line)
        m === nothing && continue
        lines[k] = string(m.captures[1], '"', value, '"', m.captures[2])
        return join(lines, '\n')
    end
    return nothing
end

"""
    repin(kind, repository; rev=nothing, root=ROOT) -> Vector{String}

Point a contributor repository at this checkout's commit (`rev`, which must be its
`HEAD`; it defaults to `HEAD`): the `rev` of its pin file and, when the file has one,
its `version`; the ref after `@` and the `spaceagora-rev` of every call of
`contributor-ci.yml` in `.github/workflows/`; and the `SpaceAGORA` entry in `[compat]`
of its `Project.toml` and `test/Project.toml`, when that entry does not admit this
core's version, which then becomes this core's release series. Only those values
change; the rest of each file keeps its text. A contributor moving to a release runs
it from a checkout of the release, and core CI runs it on a checkout of a contributor
repository to test that repository against the commit under test. Returns the changed
files, relative to the repository, and changes nothing when it refuses (a path pin,
a pin file it cannot rewrite in place).
"""
function repin(kind, repository; rev=nothing, root=ROOT)
    kind in KINDS || fail("Unknown kind '$kind'. Choose from $(join(KINDS, ", ")).")
    core = verified_core(root, rev)
    isdir(repository) || fail("Repository directory not found: $repository")
    repository = realpath(repository)
    version = isempty(core.version) ? nothing : VersionNumber(core.version)
    series = version === nothing ? nothing : version.major == 0 ? "0.$(version.minor)" : "$(version.major)"
    # Every edit is decided before any file is written, so a refusal changes nothing.
    edits = Pair{String,String}[]
    edit!(path, text, updated) = (updated == text || push!(edits, path => updated); nothing)

    pinfile = joinpath(repository, PIN_FILE)
    isfile(pinfile) || fail("$pinfile is missing. It records the SpaceAGORA commit this repository is tested against.")
    text = read(pinfile, String)
    pin = TOML.parse(text)
    haskey(pin, "path") && fail("$pinfile: a path pin names a checkout, not a commit; there is no commit to move.")
    old = get(pin, "rev", nothing)
    full_sha(old) || fail("$pinfile: rev must be a full 40-character lowercase commit SHA, not $(repr(old)).")
    updated = replace_toml_string(text, "rev", core.rev)
    haskey(pin, "version") && updated !== nothing && (updated = replace_toml_string(updated, "version", core.version))
    (updated !== nothing && get(TOML.parse(updated), "rev", nothing) == core.rev &&
        get(TOML.parse(updated), "version", core.version) == core.version) ||
        fail("$pinfile: could not rewrite rev and version in place; set rev = \"$(core.rev)\" and version = \"$(core.version)\" by hand.")
    edit!(pinfile, text, updated)

    # Every commit the workflow calls name, and the old pin, become this commit.
    workflows = joinpath(repository, ".github", "workflows")
    for file in (isdir(workflows) ? sort!(readdir(workflows)) : String[])
        endswith(file, ".yml") || endswith(file, ".yaml") || continue
        path = joinpath(workflows, file)
        text = read(path, String)
        calls = [m.captures[1] for m in eachmatch(WORKFLOW_CALL, text)]
        isempty(calls) && continue
        edit!(path, text, foldl((t, sha) -> replace(t, sha => core.rev), unique!(push!(calls, old)); init=text))
    end

    # A [compat] bound that does not admit this core's version becomes its release series.
    for directory in (repository, joinpath(repository, "test"))
        file = project_file(directory)
        (series === nothing || file === nothing) && continue
        text = read(file, String)
        bound = get(get(TOML.parse(text), "compat", Dict{String,Any}()), "SpaceAGORA", nothing)
        (bound === nothing || version in Pkg.Versions.semver_spec(bound)) && continue
        updated = replace_toml_string(text, "SpaceAGORA", series; table="compat")
        updated === nothing && fail("$file: could not rewrite [compat] SpaceAGORA in place; set it to \"$series\" by hand.")
        edit!(file, text, updated)
    end

    for (path, updated) in edits
        write(path, updated)
    end
    return [relpath(path, repository) for (path, _) in edits]
end

# --- internal-name scan --------------------------------------------------------------------------

internal_name(name::Symbol) = name in INTERNAL_MODULES || startswith(String(name), "_")

"""`a.b.c` as `[:a, :b, :c]`, from a property access or an import path; `nothing` otherwise."""
function dotted(ex)
    ex isa Symbol && return Symbol[ex]
    ex isa Expr && ex.head === :. || return nothing
    if length(ex.args) == 2 && ex.args[2] isa QuoteNode && ex.args[2].value isa Symbol
        left = dotted(ex.args[1])
        return left === nothing ? nothing : push!(left, ex.args[2].value)
    end
    return all(arg -> arg isa Symbol, ex.args) ? collect(Symbol, ex.args) : nothing
end

"""The first internal segment of a path rooted at SpaceAGORA (or an alias of it)."""
function internal_path(chain, roots)
    (chain === nothing || length(chain) < 2 || !(first(chain) in roots)) && return nothing
    index = findfirst(internal_name, @view chain[2:end])
    return index === nothing ? nothing : join(chain[1:index+1], ".")
end

# Aliases of SpaceAGORA follow lexical scope. A module (or file) collects its top-level
# aliases: `import SpaceAGORA as S` and `const S = SpaceAGORA`. A function, `let`, loop,
# comprehension or `catch` scope drops each alias that its arguments, bindings or
# assignments shadow, and adds the aliases it assigns itself (`S = SpaceAGORA`).

"""Add the names that an argument or assignment target binds: `x`, `x::T`, `x=default`, `(a, b)`, `x...`."""
function bound_names!(names, ex)
    if ex isa Symbol
        push!(names, ex)
    elseif ex isa Expr
        if ex.head in (:tuple, :parameters, :block, :local)
            foreach(arg -> bound_names!(names, arg), ex.args)
        elseif ex.head in (:kw, :(=), :...) || (ex.head === :(::) && length(ex.args) == 2)
            bound_names!(names, ex.args[1])
        end
    end
    return names
end

is_signature(ex) = ex isa Expr && (ex.head === :call || (ex.head in (:where, :(::)) && is_signature(ex.args[1])))

# A `where` parameter: `T`, `T<:U`, `T>:L` or `L<:T<:U`. The parameter binds inside the
# function; its bounds are annotations of the signature.
function where_parameter!(parameters, annotations, ex)
    if ex isa Symbol
        push!(parameters, ex)
    elseif ex isa Expr && ex.head === :braces
        foreach(arg -> where_parameter!(parameters, annotations, arg), ex.args)
    elseif ex isa Expr && ex.head in (:<:, :>:) && length(ex.args) == 2
        bound_names!(parameters, ex.args[1])
        push!(annotations, ex.args[2])
    elseif ex isa Expr && ex.head === :comparison && length(ex.args) == 5
        bound_names!(parameters, ex.args[3])
        push!(annotations, ex.args[1], ex.args[5])
    else
        push!(annotations, ex)
    end
end

"""`(callee, arguments, body, annotations, parameters)` of a function, macro or anonymous
function, or `nothing`. The annotations are its return type and `where` bounds, which are
scanned with the callee in the enclosing scope; the `where` parameters bind inside it."""
function function_parts(ex)
    ex isa Expr && length(ex.args) == 2 || return nothing
    annotations, parameters = Any[], Set{Symbol}()
    if ex.head in (:function, :macro) || (ex.head === :(=) && is_signature(ex.args[1]))
        signature = ex.args[1]
        while signature isa Expr && signature.head in (:where, :(::)) && signature.args[1] isa Expr &&
                signature.args[1].head in (:call, :tuple, :where, :(::))
            if signature.head === :(::)
                push!(annotations, signature.args[2])
            else
                foreach(arg -> where_parameter!(parameters, annotations, arg), signature.args[2:end])
            end
            signature = signature.args[1]
        end
        signature isa Expr && signature.head === :call &&
            return (signature.args[1], signature.args[2:end], ex.args[2], annotations, parameters)
        signature isa Expr && signature.head === :tuple && return (nothing, signature.args, ex.args[2], annotations, parameters)
        return (nothing, Any[signature], ex.args[2], annotations, parameters)
    elseif ex.head === :->
        signature = ex.args[1]
        arguments = signature isa Expr && signature.head in (:tuple, :block) ?
            filter(arg -> !(arg isa LineNumberNode), signature.args) : Any[signature]
        return (nothing, arguments, ex.args[2], annotations, parameters)
    end
    return nothing
end

const SCOPE_HEADS = (:module, :let, :for, :while, :try, :generator)

"""Collect a module's top-level aliases of SpaceAGORA into `roots`, without entering local scopes."""
function module_aliases!(roots, ex)
    ex isa Expr || return roots
    (function_parts(ex) !== nothing || ex.head in SCOPE_HEADS) && return roots
    if ex.head === :as && dotted(ex.args[1]) == [:SpaceAGORA] && ex.args[2] isa Symbol
        push!(roots, ex.args[2])
    elseif ex.head === :(=) && length(ex.args) == 2 && ex.args[1] isa Symbol && ex.args[2] isa Symbol && ex.args[2] in roots
        push!(roots, ex.args[1])
    end
    foreach(arg -> module_aliases!(roots, arg), ex.args)
    return roots
end

"""Collect the names a scope body assigns (its locals), the aliases among them, and its `global` names."""
function assigned!(locals, aliases, globals, ex, roots)
    ex isa Expr || return
    (function_parts(ex) !== nothing || ex.head === :module) && return
    if ex.head === :(=) && length(ex.args) == 2
        bound_names!(locals, ex.args[1])
        ex.args[1] isa Symbol && ex.args[2] isa Symbol && (ex.args[2] in roots || ex.args[2] in aliases) &&
            push!(aliases, ex.args[1])
    elseif ex.head === :local
        bound_names!(locals, ex)
    elseif ex.head === :global
        foreach(arg -> bound_names!(globals, arg), ex.args)
    elseif ex.head === :tuple
        # Named-tuple fields, `(S = 1,)`, bind nothing.
        for arg in ex.args
            assigned!(locals, aliases, globals, arg isa Expr && arg.head === :(=) ? arg.args[2] : arg, roots)
        end
        return
    end
    foreach(arg -> assigned!(locals, aliases, globals, arg, roots), ex.args)
end

"""The aliases visible in a scope that binds `bound` and runs `body`, inside a scope that sees `roots`."""
function scope_roots(roots, bound, body)
    locals, aliases, globals = Set{Symbol}(bound), Set{Symbol}(), Set{Symbol}()
    assigned!(locals, aliases, globals, body, roots)
    return union(setdiff(roots, setdiff(locals, globals)), aliases)
end

internal!(found, line, name) = name === nothing || push!(found, (line, "$name is internal to SpaceAGORA"))

# Iteration specifications `x = iter` of a loop or generator: bind the targets and scan the
# iterables. In `for a in A, b in B`, B sees `a` (`sequential`); the comma-separated
# iterables of a comprehension are all evaluated outside it.
function iteration!(found, specs, line, roots, bound; sequential=false)
    for spec in specs
        spec isa LineNumberNode && continue
        if spec isa Expr && spec.head === :block
            line = iteration!(found, spec.args, line, roots, bound; sequential)
        elseif spec isa Expr && spec.head === :filter
            line = iteration!(found, spec.args[2:end], line, roots, bound; sequential)
        elseif spec isa Expr && spec.head === :(=)
            line = scan_expression!(found, spec.args[2], line, sequential ? setdiff(roots, bound) : roots)
            bound_names!(bound, spec.args[1])
        end
    end
    return line
end

function scan_expression!(found, ex, line, roots)
    ex isa LineNumberNode && return ex.line
    ex isa Expr || return line
    parts = function_parts(ex)
    if ex.head in (:error, :incomplete)
        push!(found, (line, "the file does not parse"))
        return line
    elseif ex.head === :module
        return scan_expression!(found, ex.args[end], line, module_aliases!(Set([:SpaceAGORA]), ex.args[end]))
    elseif parts !== nothing
        callee, arguments, body, annotations, parameters = parts
        callee === nothing || (line = scan_expression!(found, callee, line, roots))
        for annotation in annotations
            line = scan_expression!(found, annotation, line, roots)
        end
        inner = scope_roots(roots, bound_names!(copy(parameters), Expr(:tuple, arguments...)), body)
        for argument in arguments
            line = scan_expression!(found, argument, line, inner)
        end
        return scan_expression!(found, body, line, inner)
    elseif ex.head === :let
        bindings, body = ex.args
        bound, aliases = Set{Symbol}(), Set{Symbol}()
        for binding in (bindings isa Expr && bindings.head === :block ? bindings.args : Any[bindings])
            binding isa LineNumberNode && continue
            bound_names!(bound, binding)
            if binding isa Expr && binding.head === :(=)
                line = scan_expression!(found, binding.args[2], line, roots)
                binding.args[1] isa Symbol && binding.args[2] isa Symbol && binding.args[2] in roots && push!(aliases, binding.args[1])
            end
        end
        return scan_expression!(found, body, line, union(scope_roots(roots, bound, body), aliases))
    elseif ex.head in (:for, :while)
        bound = Set{Symbol}()
        ex.head === :for ? (line = iteration!(found, Any[ex.args[1]], line, roots, bound; sequential=true)) :
            (line = scan_expression!(found, ex.args[1], line, roots))
        return scan_expression!(found, ex.args[2], line, scope_roots(roots, bound, ex.args[2]))
    elseif ex.head === :generator
        bound = Set{Symbol}()
        line = iteration!(found, ex.args[2:end], line, roots, bound)
        inner = setdiff(roots, bound)
        for spec in ex.args[2:end]
            spec isa Expr && spec.head === :filter && (line = scan_expression!(found, spec.args[1], line, inner))
        end
        return scan_expression!(found, ex.args[1], line, inner)
    elseif ex.head === :try
        line = scan_expression!(found, ex.args[1], line, roots)
        if length(ex.args) >= 3
            variable = ex.args[2]
            bound = variable isa Symbol ? Set([variable]) : Set{Symbol}()
            line = scan_expression!(found, ex.args[3], line, scope_roots(roots, bound, ex.args[3]))
        end
        for rest in ex.args[4:end]
            line = scan_expression!(found, rest, line, roots)
        end
        return line
    elseif ex.head === :.
        name = internal_path(dotted(ex), roots)
        if name !== nothing
            internal!(found, line, name)
            return line
        end
    elseif ex.head in (:using, :import)
        for arg in ex.args
            arg = arg isa Expr && arg.head === :as ? arg.args[1] : arg
            if arg isa Expr && arg.head === :(:)
                base = dotted(arg.args[1])
                for item in arg.args[2:end]
                    item = item isa Expr && item.head === :as ? item.args[1] : item
                    chain = dotted(item)
                    (base === nothing || chain === nothing) && continue
                    internal!(found, line, internal_path(vcat(base, chain), roots))
                end
            else
                internal!(found, line, internal_path(dotted(arg), roots))
            end
        end
        return line
    elseif ex.head === :call && length(ex.args) == 3 && ex.args[1] in (:getproperty, :getfield) &&
            ex.args[3] isa QuoteNode && ex.args[3].value isa Symbol
        chain = dotted(ex.args[2])
        chain === nothing || internal!(found, line, internal_path(push!(chain, ex.args[3].value), roots))
    end
    for arg in ex.args
        line = scan_expression!(found, arg, line, roots)
    end
    return line
end

"""References to `SpaceAGORA.SimulationModel`, `SpaceAGORA.SimulationEngine` and underscore names."""
function scan_source(text::AbstractString, filename::AbstractString)
    ex = Meta.parseall(text; filename)
    found = Tuple{Int,String}[]
    scan_expression!(found, ex, 0, module_aliases!(Set([:SpaceAGORA]), ex))
    return found
end

function scan(kind, repository)
    kind in KINDS || fail("Unknown kind '$kind'. Choose from $(join(KINDS, ", ")).")
    kind in SCANNED_KINDS || return Tuple{String,Int,String}[]
    repository = realpath(repository)
    findings = Tuple{String,Int,String}[]
    for (directory, folders, files) in walkdir(repository)
        filter!(folder -> !startswith(folder, "."), folders)
        for file in sort!(filter(endswith(".jl"), files))
            path = joinpath(directory, file)
            for (line, message) in scan_source(read(path, String), path)
                push!(findings, (relpath(path, repository), line, message))
            end
        end
    end
    return findings
end

# --- CI summary ----------------------------------------------------------------------------------

cell(value) = replace(string(value), "|" => "\\|", "\n" => " ")

"""Markdown with one line per check id from the reports in `directory`."""
function summary(directory, environment=nothing)
    io = IOBuffer()
    println(io, "## Contributor conformance\n")
    println(io, "Interface conformance only; scientific validation is the author's.\n")
    record = environment === nothing ? nothing : joinpath(environment, RECORD)
    if record !== nothing && isfile(record)
        r = TOML.parsefile(record)
        package = isempty(r["package"]) ? "" : " (`$(r["package"])`)"
        println(io, "- Kind: $(r["kind"])$package; manifest: $(r["manifest"])")
        println(io, "- SpaceAGORA: `$(r["spaceagora_rev"])`, pinned in $PIN_FILE by $(r["spaceagora_pin"])")
        for (name, companion) in sort!(collect(get(r, "companions", Dict{String,Any}())); by=first)
            println(io, "- Companion $name: $(companion["package"]) `$(companion["revision"])`")
        end
        println(io)
    end
    files = isdir(directory) ? sort!(filter(endswith(".toml"), readdir(directory; join=true))) : String[]
    if isempty(files)
        println(io, "No conformance reports were written: the tests did not call the conformance checks.")
        return String(take!(io))
    end
    counts = Dict{String,Int}()
    println(io, "| Target | Check | Status | Since | Detail |")
    println(io, "|---|---|---|---|---|")
    for file in files
        report = TOML.parsefile(file)
        target = string(get(report, "target", basename(file)), " (", get(report, "kind", "?"), ")")
        for check in get(report, "checks", Any[])
            status = string(get(check, "status", "?"))
            counts[status] = get(counts, status, 0) + 1
            println(io, "| ", cell(target), " | ", cell(get(check, "id", "?")), " | ", cell(status), " | ",
                    cell(get(check, "since", "")), " | ", cell(get(check, "detail", "")), " |")
        end
    end
    println(io, "\n", join(("$(n) $(status)" for (status, n) in sort!(collect(counts); by=first)), ", "))
    return String(take!(io))
end

# --- command line --------------------------------------------------------------------------------

function options(args)
    positional, named = String[], Dict{String,String}()
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--core-from-caller"
            named["core-from-caller"] = "true"
        elseif startswith(arg, "--")
            key, value = occursin('=', arg) ? split(arg[3:end], '='; limit=2) :
                (arg[3:end], i < length(args) ? args[i += 1] : fail("Option $arg needs a value.\n$HELP"))
            key in ("manifest", "rev") || fail("Unknown option --$key.\n$HELP")
            named[key] = value
        else
            push!(positional, arg)
        end
        i += 1
    end
    return positional, named
end

function main(args=ARGS)
    positional, named = options(args)
    core_from_caller = get(named, "core-from-caller", "") == "true"
    if length(positional) == 3 && positional[1] == "scan" && isempty(named)
        findings = scan(positional[2], positional[3])
        for (file, line, message) in findings
            println(stderr, "$file:$line: $message")
        end
        isempty(findings) || exit(1)
        println(positional[2] in SCANNED_KINDS ?
            "Internal-name scan: no SimulationModel, SimulationEngine or underscore names in $(positional[3])" :
            "Internal-name scan: a $(positional[2]) is not scanned; it uses its facade's bindings")
    elseif length(positional) in (2, 3) && positional[1] == "summary" && isempty(named)
        print(summary(positional[2], get(positional, 3, nothing)))
    elseif length(positional) == 3 && positional[1] == "repin" && !haskey(named, "manifest") && !core_from_caller
        changed = repin(positional[2], positional[3]; rev=get(named, "rev", nothing))
        core = verified_core(ROOT, get(named, "rev", nothing))
        println(isempty(changed) ? "$(positional[3]) already pins SpaceAGORA $(core.rev)" :
            "Pinned $(positional[3]) to SpaceAGORA $(core.rev) ($(core.version)): $(join(changed, ", "))")
    elseif length(positional) == 3 && positional[1] == "check" && !haskey(named, "manifest")
        c = check_repository(positional[2], positional[3]; rev=get(named, "rev", nothing), core_from_caller)
        println("SpaceAGORA pin: ", c.pin, " in ", joinpath(c.repository, PIN_FILE), " matches ", c.core.path)
        isempty(c.companions) || println("Companions at this core's pins: ", join(c.companions, ", "))
    elseif length(positional) == 3 && haskey(named, "manifest")
        environment = setup(positional...; manifest=named["manifest"], rev=get(named, "rev", nothing), core_from_caller)
        println("Environment: ", environment)
    else
        fail(HELP)
    end
end

end

# A refusal prints its message and exits nonzero, without the stack trace of a crash;
# any other error keeps its trace.
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    try
        ContributorEnvironment.main()
    catch error
        error isa ArgumentError || rethrow()
        println(stderr, "ERROR: ", error.msg)
        exit(1)
    end
end
