# Service facades: the shape the optional adapter standard fixes for every
# versioned boundary module (docs/src/reference/adapter_contract.md).

"""
    check_service_facade(facade::Module; version::VersionNumber) -> ConformanceReport

Check a versioned service facade against the optional adapter standard. One check,
`facade_shape`, verifies that:

- the facade defines `CONTRACT_VERSION`, `CompatibilityError`, `require_version`,
  `IMPLEMENTED_FUNCTIONS` and `CONSUMED_SERVICES`;
- `CONTRACT_VERSION` equals `version`, the contract version the caller reviewed;
- `require_version(version)` returns `nothing`, and every neighbouring version (one
  patch, minor or major release up or down, and `v"0.0.0"`) throws the facade's
  `CompatibilityError`: the check is exact;
- the two inventories are collections of names or name paths (tuples of names), each
  listed once, none in both; every implemented name resolves to a function and every
  consumed one to a value.

A facade that breaks any of these gets a failing `facade_shape` with the problems in
its detail; the check does not throw on a malformed facade.
"""
function check_service_facade(facade::Module; version::VersionNumber)
    problems = _facade_problems(facade, version)
    check = _check(:facade_shape, FIRST_RELEASE, isempty(problems), join(problems, "\n"))
    return ConformanceReport(string(facade), :service_facade, [check])
end

"""
    test_service_facade(facade::Module; version::VersionNumber)

Run [`check_service_facade`](@ref) and record its report with [`test_report`](@ref).
Needs the `Test` standard library: `using Test` loads the method.
"""
function test_service_facade end

const _FACADE_NAMES = (:CONTRACT_VERSION, :CompatibilityError, :require_version, :IMPLEMENTED_FUNCTIONS, :CONSUMED_SERVICES)

function _version_neighbours(v::VersionNumber)
    candidates = [v"0.0.0", VersionNumber(v.major, v.minor, v.patch + 1), VersionNumber(v.major, v.minor + 1, 0),
        VersionNumber(v.major + 1, 0, 0)]
    v.patch > 0 && push!(candidates, VersionNumber(v.major, v.minor, v.patch - 1))
    v.minor > 0 && push!(candidates, VersionNumber(v.major, v.minor - 1, 0))
    v.major > 0 && push!(candidates, VersionNumber(v.major - 1, 0, 0))
    return unique!(filter!(!=(v), candidates))
end

_binding_path(name::Symbol) = (name,)
_binding_path(path::Tuple{Vararg{Symbol}}) = path
_binding_path(entry) = nothing

function _resolve_binding(facade::Module, path)
    try
        return Some(foldl(getproperty, path; init=facade))
    catch
        return nothing
    end
end

# The name paths an inventory lists. A value that is not a collection, or entries
# that are not names or tuples of names, become problems instead of errors.
function _inventory_paths!(problems::Vector{String}, name::Symbol, inventory)
    entries = try
        collect(Any, inventory)
    catch err
        push!(problems, "$(name) is $(repr(inventory)), not a collection of names or name paths " *
            "($(first(split(_summary(err), '\n'))))")
        return Tuple{Vararg{Symbol}}[]
    end
    bad = filter(entry -> _binding_path(entry) === nothing || isempty(_binding_path(entry)), entries)
    isempty(bad) || push!(problems, "$(name) entries must be names or non-empty tuples of names: $(join(repr.(bad), ", "))")
    return [_binding_path(entry) for entry in entries if !(entry in bad)]
end

function _facade_problems(facade::Module, version::VersionNumber)
    missing_names = [name for name in _FACADE_NAMES if !isdefined(facade, name)]
    isempty(missing_names) || return ["$(facade) does not define $(join(missing_names, ", "))"]
    problems = String[]
    declared = facade.CONTRACT_VERSION
    declared == version || push!(problems, "CONTRACT_VERSION is $(repr(declared)), not the reviewed $(version)")
    E = facade.CompatibilityError
    if !(E isa Type && E <: Exception)
        push!(problems, "CompatibilityError is $(repr(E)), not an Exception type")
        E = Union{}
    end
    accepted = try
        facade.require_version(version)
    catch err
        push!(problems, "require_version($(version)) threw $(_summary(err))")
        nothing
    end
    accepted === nothing || push!(problems, "require_version($(version)) returned $(repr(accepted)), not nothing")
    for neighbour in _version_neighbours(version)
        outcome = try
            facade.require_version(neighbour)
            :accepted
        catch err
            err
        end
        if outcome === :accepted
            push!(problems, "require_version($(neighbour)) accepted a version other than $(version); it must be exact")
        elseif !(outcome isa E)
            push!(problems, "require_version($(neighbour)) threw $(typeof(outcome)), not CompatibilityError")
        end
    end
    paths = (implemented=_inventory_paths!(problems, :IMPLEMENTED_FUNCTIONS, facade.IMPLEMENTED_FUNCTIONS),
        consumed=_inventory_paths!(problems, :CONSUMED_SERVICES, facade.CONSUMED_SERVICES))
    for (label, list) in pairs(paths)
        repeated = unique(filter(path -> count(==(path), list) > 1, list))
        isempty(repeated) || push!(problems, "$(label) entries listed more than once: $(join(join.(repeated, "."), ", "))")
    end
    shared = intersect(paths.implemented, paths.consumed)
    isempty(shared) || push!(problems, "listed as both implemented and consumed: $(join(join.(shared, "."), ", "))")
    for path in paths.implemented
        value = _resolve_binding(facade, path)
        value !== nothing && something(value) isa Function && continue
        push!(problems, "implemented $(join(path, ".")) " * (value === nothing ? "does not resolve" : "is not a function"))
    end
    for path in paths.consumed
        value = _resolve_binding(facade, path)
        value !== nothing && something(value) !== nothing && continue
        push!(problems, "consumed $(join(path, ".")) " * (value === nothing ? "does not resolve" : "is nothing"))
    end
    return problems
end
