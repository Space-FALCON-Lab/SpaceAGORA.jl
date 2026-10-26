const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# The extension hooks SpaceAGORA 0.5.0 removed (docs/src/reference/migrating_extension_hooks.md
# gives each one's replacement). This gate parses every Julia file under src/,
# ext/, templates/ and examples/ and fails on any method of a removed name,
# qualified or not and including the callable form `(f::typeof(hook))(...)`, and on
# any other use of one as an identifier: a declaration, a call, an import or export,
# a quoted symbol. It also fails on the removed `wrench` form that took the
# spacecraft index as a fifth positional argument: a method that accepts five
# positional arguments (defaulted and trailing `args...` parameters included) or a
# call that passes five (keywords, after a comma or a semicolon, are not counted).
# Comments and strings are not code and are not checked here, so a name built from
# a string (`Symbol("calcForceTorque")`, `Meta.parse`, `include_string`) is out of
# its reach.
const REMOVED_HOOKS = (:calcForceTorque, :calcControlForceTorque, :calcControlMassFlowRate, :calcReactionWheelTorque,
    :calcControlEffect!, :calcGuidanceEffect!, :calcNavigationEffect!,
    :planet_frame_lpi, :ephemerides_sun_direction_ii, :ephemerides_time_seconds)
const SCANNED_DIRS = ("src", "ext", "templates", "examples")

_last_name(name::Symbol) = name
function _last_name(ex::Expr)
    ex.head === :. && return _last_name(ex.args[end])
    # `(f::typeof(hook))(...)` and `(::typeof(hook))(...)`
    if ex.head === :(::)
        T = ex.args[end]
        T isa Expr && T.head === :call && T.args[1] === :typeof && length(T.args) == 2 && return _last_name(T.args[2])
    end
    return nothing
end
_last_name(q::QuoteNode) = q.value isa Symbol ? q.value : nothing
_last_name(::Any) = nothing

# The call expression of a method signature, through `where` and return types.
function _signature_call(sig)
    while sig isa Expr && sig.head in (:where, :(::))
        sig = sig.args[1]
    end
    return sig isa Expr && sig.head === :call ? sig : nothing
end

# Positional arguments. In a call, `f(a, k=1)` and `f(a; k=1)` both pass `k` as a
# keyword: the comma form is a `:kw` node, the other a `:parameters` block. In a
# method signature a `:kw` node is a positional parameter with a default value, and
# only the `:parameters` block holds keywords; a signature with defaults or a
# trailing `args...` defines a method for each arity in a range.
_is_keywords(arg) = arg isa Expr && arg.head === :parameters
_is_default(arg) = arg isa Expr && arg.head === :kw
_is_varargs(arg) = arg isa Expr && arg.head === :...
_call_positional_count(call::Expr) = count(arg -> !(_is_keywords(arg) || _is_default(arg)), call.args[2:end])
function _signature_takes(call::Expr, n::Int)
    positional = filter(!_is_keywords, call.args[2:end])
    required = count(arg -> !(_is_default(arg) || _is_varargs(arg)), positional)
    any(_is_varargs, positional) && return required <= n
    return required <= n <= length(positional)
end

# `Meta.parseall` does not throw on a syntax error: it leaves an `:error` or
# `:incomplete` node in the tree, which would hide the rest of the file.
_parse_failed(ex) = ex isa Expr && (ex.head in (:error, :incomplete) || any(_parse_failed, ex.args))

# Each use of a removed name in `ex`, as (kind, name), with the line it is on.
# `signature` is the call node of the method definition being walked, which is not
# itself a call.
function _removed_uses!(found, ex, line::Ref{Int}, signature=nothing)
    if ex isa LineNumberNode
        line[] = ex.line
    elseif ex isa Symbol
        ex in REMOVED_HOOKS && push!(found, (line[], "uses", ex))
    elseif ex isa QuoteNode
        ex.value isa Symbol && ex.value in REMOVED_HOOKS && push!(found, (line[], "uses", ex.value))
    elseif ex isa Expr
        if ex.head === :macrocall && length(ex.args) >= 2 && ex.args[2] isa LineNumberNode
            line[] = ex.args[2].line
        end
        if ex.head === :function || (ex.head === :(=) && _signature_call(ex.args[1]) !== nothing)
            call = _signature_call(ex.args[1])
            name = call === nothing ? nothing : _last_name(call.args[1])
            if name in REMOVED_HOOKS
                push!(found, (line[], "defines a method of", name))
                foreach(arg -> _removed_uses!(found, arg, line), ex.args[2:end])
                return found
            elseif name === :wrench && _signature_takes(call, 5)
                push!(found, (line[], "defines a five-argument method of", :wrench))
                foreach(arg -> _removed_uses!(found, arg, line), ex.args[2:end])
                return found
            elseif call === nothing && ex.head === :function && length(ex.args) == 1 && _last_name(ex.args[1]) in REMOVED_HOOKS
                push!(found, (line[], "declares", _last_name(ex.args[1])))
                return found
            end
            if call !== nothing
                _removed_uses!(found, ex.args[1], line, call)
                foreach(arg -> _removed_uses!(found, arg, line), ex.args[2:end])
                return found
            end
        end
        if ex.head === :call && ex !== signature && _last_name(ex.args[1]) === :wrench && _call_positional_count(ex) == 5
            push!(found, (line[], "calls the five-argument form of", :wrench))
        end
        foreach(arg -> _removed_uses!(found, arg, line, signature), ex.args)
    end
    return found
end

violations = String[]
for dir in SCANNED_DIRS
    base = joinpath(REPO_ROOT, dir)
    isdir(base) || continue
    for (root, _, files) in walkdir(base), file in sort(files)
        endswith(file, ".jl") || continue
        path = joinpath(root, file)
        rel = relpath(path, REPO_ROOT)
        parsed = try
            Meta.parseall(read(path, String); filename=rel)
        catch err
            push!(violations, "$(rel): does not parse ($(sprint(showerror, err)))")
            continue
        end
        if _parse_failed(parsed)
            push!(violations, "$(rel): does not parse, so it cannot be checked.")
            continue
        end
        for (line, kind, name) in unique(_removed_uses!(Tuple{Int, String, Symbol}[], parsed, Ref(0)))
            what = name === :wrench ? "wrench, the removed form with the spacecraft index as an argument; " *
                "read x.spacecraft_index instead" : "the removed hook $(name); use its replacement " *
                "(docs/src/reference/migrating_extension_hooks.md)"
            push!(violations, "$(rel):$(line): $(kind) $(what).")
        end
    end
end

# The contributor starters are what new repositories copy, so they do not name a
# removed hook at all: not in a call, a comment, a README or a workflow.
const STARTERS = joinpath("templates", "starters")
_named(text, hook) = occursin(Regex("(?<![A-Za-z0-9_])" * String(hook) * "(?![A-Za-z0-9_!])"), text)
for (root, _, files) in (isdir(joinpath(REPO_ROOT, STARTERS)) ? walkdir(joinpath(REPO_ROOT, STARTERS)) : ()), file in sort(files)
    path = joinpath(root, file)
    rel = relpath(path, REPO_ROOT)
    text = read(path, String)
    for hook in REMOVED_HOOKS
        _named(text, hook) && push!(violations,
            "$(rel): names the removed hook $(hook); the starters use only current hook names.")
    end
end

if !isempty(violations)
    println("removed_hook_names_violations:")
    foreach(v -> println("  - ", v), violations)
    error("Removed hook names found in src/, ext/, templates/ or examples/.")
end
println("removed_hook_names_gate_ok")
