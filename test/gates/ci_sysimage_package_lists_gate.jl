using TOML

# build/build_sysimage.jl names the packages it bakes into the sysimage by
# symbol, in two lists (FULL_PACKAGES and CI_PACKAGES). PackageCompiler resolves
# those names against the root Project.toml, so every name must be :SpaceAGORA
# (the package itself) or a key of the root [deps] table; a name that is
# neither asks for a package the project no longer depends on, which is how
# Quaternions and SatelliteToolbox stayed in both lists after their last use
# was removed. This gate parses the builder and the project file only; it does
# not load PackageCompiler or build anything.
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const SYSIMAGE_BUILDER_RELATIVE = "build/build_sysimage.jl"

function sysimage_package_list(builder::String, list::Symbol)::Vector{Symbol}
    parsed = Meta.parseall(read(builder, String); filename=builder)
    any(e -> e isa Expr && e.head === :error, parsed.args) &&
        error("$SYSIMAGE_BUILDER_RELATIVE does not parse")
    definitions = Expr[]
    for e in parsed.args
        e isa Expr && e.head === :const && length(e.args) == 1 || continue
        assignment = e.args[1]
        assignment isa Expr && assignment.head === :(=) && assignment.args[1] === list || continue
        push!(definitions, assignment)
    end
    length(definitions) == 1 ||
        error("$SYSIMAGE_BUILDER_RELATIVE must define `const $list` exactly once; found $(length(definitions))")
    value = definitions[1].args[2]
    # A literal vector of symbols is the contract; any other form (a computed
    # list, a broadcast, a symbol spelled some other way) cannot be checked
    # statically, so it fails rather than passing unchecked.
    value isa Expr && value.head === :vect && all(a -> a isa QuoteNode && a.value isa Symbol, value.args) ||
        error("$SYSIMAGE_BUILDER_RELATIVE: `const $list` must be a literal vector of symbols")
    return Symbol[a.value for a in value.args]
end

function sysimage_unknown_packages(repo_root::String)::Vector{Pair{Symbol,Vector{String}}}
    project = TOML.parsefile(joinpath(repo_root, "Project.toml"))
    known = Set{String}(keys(get(project, "deps", Dict{String,Any}())))
    push!(known, "SpaceAGORA")
    builder = joinpath(repo_root, "build", "build_sysimage.jl")
    unknown = Pair{Symbol,Vector{String}}[]
    for list in (:FULL_PACKAGES, :CI_PACKAGES)
        names = String.(sysimage_package_list(builder, list))
        bad = sort(unique(filter(n -> !(n in known), names)))
        isempty(bad) || push!(unknown, list => bad)
    end
    return unknown
end

let unknown = sysimage_unknown_packages(REPO_ROOT)
    if !isempty(unknown)
        detail = join(("$list: $(join(bad, ", "))" for (list, bad) in unknown), "; ")
        error("Sysimage package lists in $SYSIMAGE_BUILDER_RELATIVE name packages that are not " *
              ":SpaceAGORA or a Project.toml [deps] entry: $detail")
    end
end

println("sysimage_package_lists_gate_ok")
