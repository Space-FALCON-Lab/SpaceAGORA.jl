module EDGOwnershipChecks
# These checks cover the typed route as well as the legacy boundary gate.
const OWNER = "src/models/gnc/guidance/aerobraking/typed_edg/"
# The two decisions SpaceAGORA forwards to the companion, the names its facade
# lists in IMPLEMENTED_FUNCTIONS. Every other EDG kernel is the companion's alone.
const KERNELS = ("control_decision!", "guidance_decision!")
const EXPECTED = Dict(name => OWNER * "algorithms.jl" for name in KERNELS)
EXPECTED["_edg_control_sat_state"] = OWNER * "services.jl"
EXPECTED["_edg_control_pos_vel_mass"] = OWNER * "services.jl"
EXPECTED["_edg_environment_state"] = OWNER * "services.jl"
EXPECTED["_edg_in_drag_passage"] = OWNER * "services.jl"
EXPECTED["_edg_ephemeris_time"] = OWNER * "services.jl"
EXPECTED["_edg_planet_frame_lpi"] = OWNER * "services.jl"
EXPECTED["_edg_targeting_prediction_environment"] = OWNER * "services.jl"
EXPECTED["_edg_sample_prediction_atmosphere"] = OWNER * "services.jl"
const RETIRED_KERNEL_FILES = ("heat_rate.jl", "heat_load.jl", "structural_load.jl",
    "targeting.jl", "targeting_plan.jl", "pass_heat.jl", "guidance_decision.jl", "angle_decision.jl")
const SERVICE_OVERLOADS = Set(("_edg_environment_state", "_edg_in_drag_passage", "_edg_ephemeris_time",
    "_edg_planet_frame_lpi", "_edg_targeting_prediction_environment", "_edg_sample_prediction_atmosphere"))
const FACADE_FILE = "src/models/gnc/companions/edg/services.jl"
const FACADE_FORWARD = Dict(
    "_edg_environment_state" => :(_environment(p, u, t, spacecraft_index)),
    "_edg_in_drag_passage" => :(_in_passage(p, env)),
    "_edg_ephemeris_time" => :(_epoch(p, t)),
    "_edg_targeting_prediction_environment" => :(_trajectory(p, r, v, t)),
)

# Companion helper names. No GNC module binds one to EDGAlgorithms, imports one
# from it or reaches one through it.
_is_helper_name(name) = startswith(name, "_edg_") || startswith(name, "_energy_depletion_")

function source_map(root)
    result = Dict{String,String}()
    for (dir, _, files) in walkdir(joinpath(root, "src"))
        for file in files
            endswith(file, ".jl") || continue
            path = joinpath(dir, file)
            result[replace(relpath(path, root), '\\' => '/')] = read(path, String)
        end
    end
    return result
end


# Parse syntax without evaluating it. Parenthesized include arguments and
# indented, short-form, macro-wrapped or qualified methods remain visible.
function _walk(f, ex)
    ex isa Expr || return
    ex.head == :quote && return
    f(ex)
    foreach(arg -> _walk(f, arg), ex.args)
end

function _leaf_name(ex)
    ex isa Symbol && return String(ex)
    ex isa Expr && ex.head == :. && ex.args[end] isa QuoteNode &&
        return String(ex.args[end].value)
    return nothing
end

function _signature(ex)
    while ex isa Expr && ex.head in (:(::), :where)
        ex = ex.args[1]
    end
    return ex
end

function _definition(ex)
    ex.head in (:function, :(=)) || return nothing
    sig = _signature(ex.args[1])
    sig isa Expr && sig.head == :call || return nothing
    name = _leaf_name(sig.args[1])
    return name === nothing ? nothing : (name, sig.args[1], ex.args[2])
end

function _include_calls(source)
    calls = Expr[]
    _walk(Meta.parseall(source)) do ex
        ex.head == :call && _leaf_name(ex.args[1]) == "include" && push!(calls, ex)
    end
    return calls
end

function _contains_control_literal(ex)
    ex isa String && return occursin("control", ex)
    ex isa Expr || return false
    return any(_contains_control_literal, ex.args)
end

has_control_include(source) =
    any(call -> any(_contains_control_literal, call.args[2:end]), _include_calls(source))

_syntax(ex) = ex isa Expr ?
    Expr(ex.head, (_syntax(arg) for arg in ex.args if !(arg isa LineNumberNode))...) : ex

function aggregator_violations(source)
    errors = String[]
    isempty(_include_calls(source)) || push!(errors, "$(OWNER)algorithms.jl: numerical source includes are forbidden")
    definitions = Dict{String,Int}()
    _walk(Meta.parseall(source)) do ex
        # The extension checks the facade's IMPLEMENTED_FUNCTIONS; there is no second list.
        ex.head == :const && push!(errors, "$(OWNER)algorithms.jl: forwarded names are the facade's IMPLEMENTED_FUNCTIONS, not a list here")
        definition = _definition(ex)
        definition === nothing && return
        name, _, _ = definition
        expected = Meta.parse("$name(args...; kwargs...) = _invoke_kernel(:$name, args...; kwargs...)")
        name in KERNELS && _syntax(ex) == _syntax(expected) ||
            push!(errors, "$(OWNER)algorithms.jl: $name must be an exact package forwarder")
        definitions[name] = get(definitions, name, 0) + 1
    end
    all(get(definitions, name, 0) == 1 for name in KERNELS) && length(definitions) == length(KERNELS) ||
        push!(errors, "$(OWNER)algorithms.jl: expected exactly one forwarder for each of the two decisions")
    return errors
end

# The facade's IMPLEMENTED_FUNCTIONS, which the extension accepts, are the forwarded decisions.
function facade_violations(source)
    found = Any[]
    _walk(Meta.parseall(source)) do ex
        ex.head == :const && ex.args[1] isa Expr && ex.args[1].head == :(=) &&
            ex.args[1].args[1] == :IMPLEMENTED_FUNCTIONS && push!(found, _syntax(ex.args[1].args[2]))
    end
    expected = Expr(:tuple, (QuoteNode(Symbol(name)) for name in ("guidance_decision!", "control_decision!"))...)
    return found == [expected] ? String[] :
        ["$FACADE_FILE: IMPLEMENTED_FUNCTIONS must list exactly the two forwarded decisions"]
end

_mentions_algorithms(ex) = ex === :EDGAlgorithms || ex == QuoteNode(:EDGAlgorithms) ||
    (ex isa Expr && any(_mentions_algorithms, ex.args))

# Imports of helper names from EDGAlgorithms, and qualified uses or method
# definitions of a helper through it, in any file but the forwarder file.
function helper_route_violations(path, source)
    errors = String[]
    _walk(Meta.parseall(source)) do ex
        if ex.head in (:using, :import)
            for clause in ex.args
                if clause isa Expr && clause.head == :.
                    # `import ..EDGAlgorithms._edg_name`
                    at = findfirst(==(:EDGAlgorithms), clause.args)
                    at !== nothing && at < length(clause.args) && _is_helper_name(String(clause.args[at + 1])) &&
                        push!(errors, "$path: imports EDG helper $(clause.args[at + 1]) from EDGAlgorithms")
                end
                clause isa Expr && clause.head == :(:) || continue
                origin = clause.args[1]
                origin isa Expr && origin.head == :. && last(origin.args) == :EDGAlgorithms || continue
                for item in clause.args[2:end]
                    name = item isa Expr && item.head == :. ? String(last(item.args)) :
                        item isa Expr && item.head == :as ? String(last(item.args[1].args)) : nothing
                    name !== nothing && _is_helper_name(name) &&
                        push!(errors, "$path: imports EDG helper $name from EDGAlgorithms")
                end
            end
        elseif ex.head == :. && length(ex.args) == 2 && ex.args[2] isa QuoteNode && ex.args[2].value isa Symbol
            name = String(ex.args[2].value)
            _is_helper_name(name) && _mentions_algorithms(ex.args[1]) &&
                push!(errors, "$path: reaches EDG helper $name through EDGAlgorithms")
        end
    end
    return errors
end

# Runtime checks complement the source inventory: aliases share the same
# function object, while a qualified extension carries its defining module.
function runtime_violations(algorithms, services, gnc_modules)
    errors = String[]
    for (name, path) in EXPECTED
        symbol = Symbol(name)
        owner = endswith(path, "/services.jl") ? services : algorithms
        isdefined(owner, symbol) || (push!(errors, "$owner: missing $name"); continue)
        owned = getfield(owner, symbol)
        all(method -> method.module === owner, methods(owned)) ||
            push!(errors, "$owner: foreign method on $name")
        for mod in gnc_modules
            isdefined(mod, symbol) || continue
            getfield(mod, symbol) === owned || push!(errors, "$mod: separate EDG binding $name")
        end
    end
    for mod in gnc_modules, symbol in names(mod; all=true, imported=true, usings=true)
        _is_helper_name(String(symbol)) && isdefined(mod, symbol) || continue
        value = getfield(mod, symbol)
        value isa Function && (parentmodule(value) === algorithms ||
            any(method -> method.module === algorithms, methods(value))) &&
            push!(errors, "$mod: binds EDG helper $symbol to EDGAlgorithms")
    end
    return sort!(errors)
end

function violations(sources)
    errors = String[]
    counts = Dict(name => 0 for name in keys(EXPECTED))
    for (path, source) in sources
        # Strip line comments for dependency checks; function bodies remain parsed
        # as complete top-level blocks so multiline signatures are covered.
        code = join((first(split(line, '#'; limit=2)) for line in split(source, '\n')), '\n')
        # This relocated callback adapter was outside the numerical owner
        # before the move and retains its separate no-control rule below.
        if startswith(path, OWNER) && path != OWNER * "guidance_models.jl"
            for forbidden in (r"\bControlHooks\b", r"\b_control_module\s*\(",
                              r"\b(?:rotate_link|_apply_solar_panel_aoa!)\s*\(",
                              r"\b(?:update_control!|update_guidance!)\s*\(",
                              r"\.\s*(?:α|q|heat_loads)\s*(?:\.=|=(?!=))",
                              r"\b(?:AerobrakingEnergyDepletionControlModel|SolarPanelAngleOfAttackControlModel)\b")
                occursin(forbidden, code) && push!(errors, "$path: forbidden reverse dependency or mutation: $forbidden")
            end
        end
        if path == "src/models/gnc/guidance/aerobraking/typed_edg/guidance_models.jl"
            occursin(r"\b_control_module\s*\(|\bControlHooks\b", code) &&
                push!(errors, "$path: typed guidance depends on control")
        end
        _walk(Meta.parseall(source)) do ex
            definition = _definition(ex)
            definition === nothing && return
            name, callee, body = definition
            haskey(EXPECTED, name) || return
            if path == EXPECTED[name] && callee isa Symbol
                counts[name] += 1
            elseif path == OWNER * "adapter_queries.jl" && name in SERVICE_OVERLOADS && callee isa Symbol
                counts[name] += 1
            elseif path == FACADE_FILE && haskey(FACADE_FORWARD, name) && callee isa Symbol
                _syntax(body isa Expr && body.head == :block ? only(filter(x -> !(x isa LineNumberNode), body.args)) : body) == FACADE_FORWARD[name] || push!(errors, "$path: $name must forward its capability query unchanged")
            else
                push!(errors, "$path: duplicate EDG definition $name")
            end
        end
    end
    services = OWNER * "services.jl"
    if haskey(sources, services)
        _syntax.(_include_calls(sources[services])) == [_syntax(:(include(joinpath(@__DIR__, "adapter_queries.jl"))))] ||
            push!(errors, "$services: only the reviewed adapter_queries.jl include is allowed")
    end
    aggregator = OWNER * "algorithms.jl"
    haskey(sources, aggregator) ?
        append!(errors, aggregator_violations(sources[aggregator])) :
        push!(errors, "$aggregator: missing aggregator")
    haskey(sources, FACADE_FILE) && append!(errors, facade_violations(sources[FACADE_FILE]))
    for (path, source) in sources
        path == aggregator || append!(errors, helper_route_violations(path, source))
    end
    for (name, count) in counts
        expected_count = name in SERVICE_OVERLOADS ? 2 : 1
        count == expected_count || push!(errors, "$(EXPECTED[name]): expected $expected_count $name definitions, found $count")
    end
    for file in RETIRED_KERNEL_FILES
        haskey(sources, OWNER * file) && push!(errors, "$(OWNER)$file: numerical implementation belongs to SpaceAGORAEDG")
    end
    return sort!(errors)
end
end # module EDGOwnershipChecks
