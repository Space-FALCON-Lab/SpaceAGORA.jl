module RuntimeAnyHotpathChecks

_relative_path(path, root) = replace(relpath(path, root), '\\' => '/')

# Syntax fence for explicit Vector{Any} and Dict{Symbol, Any} types, including
# nested and Base/Core-qualified forms. Not type inference: aliases, Any[...]
# literals and arbitrary macro expansion are outside this contract. Parse only.
_without_lines(x::Expr) = Expr(x.head, (_without_lines(a) for a in x.args if !(a isa LineNumberNode))...)
_without_lines(x) = x

function _type_name(x)
    x isa Symbol && return x
    if x isa Expr && x.head === :. && length(x.args) == 2 && x.args[1] in (:Base, :Core)
        return x.args[2] isa QuoteNode ? x.args[2].value : x.args[2]
    end
    return nothing
end

function _forbidden_type(x::Expr)
    x.head === :curly || return nothing
    a = x.args
    length(a) == 2 && _type_name(a[1]) === :Vector && _type_name(a[2]) === :Any && return "Vector{Any}"
    length(a) == 3 && _type_name(a[1]) === :Dict && _type_name(a[2]) === :Symbol &&
        _type_name(a[3]) === :Any && return "Dict{Symbol, Any}"
    return nothing
end

function _owner_name(x)
    x isa Symbol && return string(x)
    if x isa Expr
        x.head in (:call, :curly, :where, :(::), :(<:)) && return _owner_name(first(x.args))
        x.head === :. && return join((_owner_name(a isa QuoteNode ? a.value : a) for a in x.args), ".")
    end
    return "<anonymous>"
end

struct BoundaryException
    path::String
    owners::Tuple
    statement::Expr
    reason::String
end

function _exceptions(runtime_types_path, callback_dir)
    event_path = joinpath(callback_dir, "event_callbacks.jl")
    entry(path, owners, statement, reason) = BoundaryException(
        normpath(path), owners, _without_lines(Meta.parse(statement)), reason)
    # ONE complete statement at its actual enclosing owner per entry. Preserve
    # the gate's occurrence inventory: missing and duplicated entries both fail.
    # Removing a boundary requires an explicit inventory update in the same review.
    return [
        entry(runtime_types_path, ((:module, "ConfigTypes"), (:struct, "GramDensityPerturbationState")),
            "walk_models::Vector{Any}", "Retained opt-in GRAM walk hook handles."),
        entry(runtime_types_path, ((:module, "ConfigTypes"), (:struct, "SharedBuffers")),
            "external_runtimes::Vector{Any} = Any[]", "Optional adapters return their own per-run runtime types."),
        entry(runtime_types_path, ((:module, "ConfigTypes"),),
            "const SaveData = Dict{Symbol, Any}", "Public output/persistence payload boundary."),
        entry(event_path, ((:struct, "TrajectoryRecorder"),),
            "storage::Dict{Symbol, Any}", "Built-in array dimensions differ; accessors assert concrete types."),
        entry(event_path, ((:struct, "TrajectoryRecorder"),),
            "fallback::Dict{Symbol, Vector{Any}}", "Custom SaveField getters may return arbitrary payloads."),
        entry(event_path, ((:function, "TrajectoryRecorder"),),
            "storage = Dict{Symbol, Any}()", "Constructor for the recorder array boundary."),
        entry(event_path, ((:function, "TrajectoryRecorder"),),
            "fallback = Dict{Symbol, Vector{Any}}()", "Constructor for the custom payload boundary."),
        entry(event_path, ((:function, "TrajectoryRecorder"),),
            "fallback[field.name] = Vector{Any}(undef, cap)", "Allocate custom SaveField payload storage."),
        entry(event_path, ((:function, "_grow_trajectory_recorder!"),),
            "grown = Vector{Any}(undef, new_cap)", "Grow the same custom SaveField payload storage."),
    ]
end

function _parse_failure(x)
    x isa Expr || return false
    x.head in (:error, :incomplete) && return true
    return any(_parse_failure, x.args)
end

function _scan_source!(violations, source, path, root, exceptions, uses)
    rel = _relative_path(path, root)
    tree = try
        Meta.parseall(source; filename=path)
    catch err
        push!(violations, "$rel:1 cannot parse source: $(sprint(showerror, err))")
        return
    end
    if _parse_failure(tree)
        push!(violations, "$rel:1 cannot parse source")
        return
    end
    function visit(x, owners, line, statement_position)
        x isa Expr || return
        for (index, allowed) in enumerate(exceptions)
            if statement_position && allowed.path == path && allowed.owners == owners &&
                    isequal(allowed.statement, _without_lines(x))
                uses[index] += 1
                uses[index] == 1 && return
                break
            end
        end
        forbidden = _forbidden_type(x)
        forbidden === nothing || push!(violations, "$rel:$line contains forbidden $forbidden")
        if x.head === :module || x.head === :struct
            child_owners = (owners..., (x.head, _owner_name(x.args[2])))
            visit(x.args[2], child_owners, line, false)
            visit(x.args[3], child_owners, line, statement_position)
            return
        elseif x.head === :function || x.head === :(=)
            signature = x.args[1]
            cursor = signature
            while cursor isa Expr && cursor.head in (:where, :(::))
                cursor = first(cursor.args)
            end
            if x.head === :function || (cursor isa Expr && cursor.head === :call)
                child_owners = (owners..., (:function, _owner_name(signature)))
                visit(signature, child_owners, line, false)
                visit(x.args[2], child_owners, line, statement_position)
                return
            end
        elseif x.head === :->
            owners = (owners..., (:function, "<anonymous>"))
        end
        # Exception eligibility follows statement/field positions, not arbitrary
        # descendants. In particular an assignment/call/return/quote containing
        # an otherwise allowed assignment does not inherit its exception. A
        # nested block must never turn eligibility back on after it was lost.
        for (index, a) in enumerate(x.args)
            if a isa LineNumberNode
                line = a.line
                continue
            end
            child_statement = false
            if statement_position
                if x.head in (:toplevel, :block)
                    child_statement = true
                elseif x.head in (:for, :while, :let)
                    child_statement = index == length(x.args)
                elseif x.head in (:if, :elseif)
                    child_statement = index > 1
                elseif x.head === :try
                    child_statement = a isa Expr && a.head === :block
                elseif x.head === :macrocall && index == length(x.args)
                    macro_name = x.args[1] isa GlobalRef ?
                        (x.args[1].mod in (Base, Core) ? x.args[1].name : nothing) :
                        _type_name(x.args[1])
                    # These existing declaration wrappers do not make a field
                    # or a documented declaration into an expression operand.
                    child_statement = macro_name === Symbol("@doc") ||
                        (macro_name === Symbol("@kwdef") && a isa Expr && a.head === :struct)
                end
            end
            # Strings and quoted symbols are data. Interpolations and quoted
            # Julia expressions are scanned, but quotes never grant exceptions.
            visit(a, owners, line, child_statement)
        end
        return
    end
    visit(tree, (), 1, true)
    return
end

function scan(root::String, hotpath_files::Vector{String}, callback_dir::String;
        runtime_types_path::String=first(hotpath_files))
    isdir(callback_dir) || error("Missing runtime callback scan root: $(_relative_path(callback_dir, root))")
    targets = copy(hotpath_files)
    for (dir, _, files) in walkdir(callback_dir), file in files
        endswith(file, ".jl") && push!(targets, joinpath(dir, file))
    end
    exceptions = _exceptions(runtime_types_path, callback_dir)
    uses = zeros(Int, length(exceptions))
    violations = String[]
    for path in sort!(unique!(normpath.(targets)))
        isfile(path) || error("Missing runtime Any hot-path target: $(_relative_path(path, root))")
        _scan_source!(violations, read(path, String), path, root, exceptions, uses)
    end
    for (entry, observed) in zip(exceptions, uses)
        observed == 1 && continue
        owner = join((name for (_, name) in entry.owners), ".")
        push!(violations,
            "$(_relative_path(entry.path, root)): heterogeneous boundary in $owner expected 1 occurrence(s), observed $observed")
    end
    return violations
end

end # module RuntimeAnyHotpathChecks
