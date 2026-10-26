# Checks how a companion's runtime code reaches SpaceAGORA. The same file serves
# the service-boundary tests of SpaceAGORAEDG, HYPR and SpaceAGORAMuJoCo; each
# states its own policy: the names it may take from a SpaceAGORA module
# (`allowed_member(module, name)`) and the SpaceAGORA modules it may step through.
# A permitted module gives access to its permitted names only. A private module
# that a permitted one happens to bind, such as a facade's `SimulationModel`, is
# private access.
#
# - `runtime_violations`: bindings of the companion's modules whose value comes
#   from SpaceAGORA and is not permitted, whether imported or defined there (a
#   local alias such as `const T = SomeAlias.SimulationModel.ODEParams` too).
# - `source_violations`: in module code and function bodies, every access chain
#   that starts at SpaceAGORA, at a name the source binds to a SpaceAGORA module, or
#   at a `getfield`/`getproperty`/`getglobal` lookup on one (also `Base.`- or
#   `Core.`-qualified), and steps through a name that is not permitted, also past a
#   type assertion inside the chain (`(Facade::Module).x`); such a lookup with a
#   computed name; and the SpaceAGORA module passed on as a value.
#   Strings and comments are not code, so documentation is not checked.
#   It follows these bindings of SpaceAGORA modules, at module level and in
#   function bodies: `X = <path>` and `X::T = <path>` (with `const`, `global` or
#   `local`), chained assignments (`a = b = <path>`), type-asserted right-hand
#   sides (`X = <path>::Module`), positional and keyword argument defaults
#   (`f(m = <path>)`, `f(; m = <path>)`), `import A.B` and `using A.B` (which bind
#   `B`), `import ... as X`, each item of `import A: x, y as z` and
#   `using A: x, y as z`, the exported modules of a plain `using A`, and relative
#   paths such as `..SpaceAGORA`. It keeps one table of names for all the files,
#   not one per scope, and reports a name bound to two different SpaceAGORA
#   modules. It also returns that table.
# - `unknown_module_bindings`: module-level bindings of the companion's modules
#   whose value is a SpaceAGORA module under a name the source table does not have:
#   a backstop for a binding form the source scan does not follow.
#
# Limits. The source check reads syntax only and follows the forms above, so it
# does not see:
# - a module value passed through a call or a data structure and used there:
#   `identity(Facade)`, tuple destructuring or indexing, lambda and comprehension
#   variables, `parentmodule(Facade.f)`, or a permitted module handed to a
#   reflective or higher-order call (`Core.eval(Facade, :x)`,
#   `Base.invokelatest(getproperty, Facade, :x)`, `getfield.(Ref(Facade), names)`);
# - a chain whose base is a block, a conditional or another compound expression
#   rather than a name, an access, a lookup or a type assertion:
#   `(begin Facade end).x`, `(c ? Facade : Other).x`, `(let m = Facade; m end).x`;
# - code made at run time: `Core.eval` or `include_string` of a string, `eval` of a
#   constructed expression, string macros, `Base.require` or `Base.root_module` of
#   a package id;
# - files it is not given: a file included from outside the scanned directories,
#   or one without the `.jl` suffix.
# Because the table is not per scope, a local that is not a module but has the name
# of a bound one is checked as that module, which reports rather than hides. The
# runtime check and the backstop see the module-level bindings of the package's
# modules only; at that level they catch the forms above that leave a binding.
module CoreAccess
export CorePolicy, runtime_violations, source_violations, unknown_module_bindings

struct CorePolicy
    core::Module
    allowed_member::Function
    modules::Vector{Module}
    values::Vector{Any}
end

function CorePolicy(core::Module, allowed_member::Function, modules)
    values = Any[]
    for m in modules, name in names(m; all=true, imported=true, usings=true)
        isdefined(m, name) && allowed_member(m, name) || continue
        value = getglobal(m, name)
        value isa Module || push!(values, value)
    end
    return CorePolicy(core, allowed_member, collect(Module, modules), values)
end

owner(value) = try
    value isa Module ? value : value isa Union{Function, Type} ? parentmodule(value) : parentmodule(typeof(value))
catch
    nothing
end
from_core(policy, value) = (o = owner(value); o !== nothing && Base.moduleroot(o) === policy.core)
allowed_module(policy, m) = any(x -> x === m, policy.modules)
permitted(policy, value) = value isa Module ? allowed_module(policy, value) :
    value isa Union{Function, Type} ? any(v -> v === value, policy.values) :
    any(v -> v isa Type && typeof(value) <: v, policy.values)

# The package's own modules, from `root` down.
function own_modules(root::Module)
    modules = Module[root]
    for M in modules, name in names(M; all=true, imported=true, usings=true)
        isdefined(M, name) || continue
        value = getglobal(M, name)
        value isa Module && value !== M && parentmodule(value) === M &&
            Base.moduleroot(value) === Base.moduleroot(root) && push!(modules, value)
    end
    return modules
end

"""Bindings of `root` and its submodules whose SpaceAGORA value is not permitted."""
function runtime_violations(policy::CorePolicy, root::Module)
    found = String[]
    own = own_modules(root)
    for M in own, name in names(M; all=true, imported=true, usings=true)
        isdefined(M, name) || continue
        value = getglobal(M, name)
        any(m -> m === value, own) && continue
        from_core(policy, value) && !permitted(policy, value) && push!(found, "$M.$name")
    end
    return found
end

"""Bindings of `root` and its submodules to a SpaceAGORA module that `known`, the source table, does not have."""
function unknown_module_bindings(policy::CorePolicy, root::Module, known::AbstractDict)
    found = String[]
    for M in own_modules(root), name in names(M; all=true, imported=true, usings=true)
        isdefined(M, name) || continue
        value = getglobal(M, name)
        value isa Module && Base.moduleroot(value) === policy.core || continue
        get(known, name, nothing) === value || push!(found, "$M.$name")
    end
    return found
end

"""Source accesses in `files` that are not permitted, the permitted `(module, name)` steps seen, and the names bound to SpaceAGORA modules."""
function source_violations(policy::CorePolicy, files)
    core = policy.core
    exprs = [Meta.parseall(read(file, String); filename=file) for file in files]
    aliases = Dict{Symbol, Module}(nameof(core) => core)
    found = String[]
    used = Set{Tuple{Module, Symbol}}()
    # A permitted step reads a name the policy allows from a module it allows; a
    # module read that way must be one it allows too.
    function step(value, name, record)
        value isa Module && Base.moduleroot(value) === core || return nothing
        member = allowed_module(policy, value) && policy.allowed_member(value, name) && isdefined(value, name) ?
            Some(getglobal(value, name)) : nothing
        if member === nothing || (something(member) isa Module && !allowed_module(policy, something(member)))
            record && push!(found, "$(value).$(name)")
            return nothing
        end
        record && push!(used, (value, name))
        return member
    end
    # What a path names, permitted or not: the source table follows private
    # modules too, so that every access through them is checked.
    lookup(value, name) = value isa Module && Base.moduleroot(value) === core && isdefined(value, name) ?
        Some(getglobal(value, name)) : nothing
    is_access(ex) = Meta.isexpr(ex, :.) && length(ex.args) == 2 && ex.args[2] isa QuoteNode
    # A lookup's function: `getfield`, or a `Base.`/`Core.`-qualified one.
    lookup_name(f) = f isa Symbol ? f :
        Meta.isexpr(f, :.) && length(f.args) == 2 && f.args[1] in (:Base, :Core) && f.args[2] isa QuoteNode ?
        f.args[2].value : nothing
    is_lookup(ex) = Meta.isexpr(ex, :call) && length(ex.args) >= 3 &&
        lookup_name(ex.args[1]) in (:getfield, :getproperty, :getglobal, :isdefined)
    is_reader(ex) = is_lookup(ex) && lookup_name(ex.args[1]) !== :isdefined
    # A type assertion `x::T` inside a chain stands for `x`.
    is_assert(ex) = Meta.isexpr(ex, :(::)) && length(ex.args) == 2
    function resolve(ex, record)
        if ex isa Symbol
            return haskey(aliases, ex) ? Some(aliases[ex]) : nothing
        elseif is_access(ex)
            base = resolve(ex.args[1], record)
            return base === nothing ? nothing : step(something(base), ex.args[2].value, record)
        elseif is_assert(ex)
            return resolve(ex.args[1], record)
        elseif is_reader(ex)
            base = resolve(ex.args[2], record)
            base === nothing && return nothing
            key = ex.args[3]
            if !(key isa QuoteNode && key.value isa Symbol)
                record && push!(found, "computed name: $(ex)")
                return nothing
            end
            return step(something(base), key.value, record)
        end
        return nothing
    end
    function value_of(ex)
        if ex isa Symbol
            return haskey(aliases, ex) ? Some(aliases[ex]) : nothing
        elseif is_assert(ex) || Meta.isexpr(ex, :(=))
            # `x::T` and the value of a chained `a = b = x`.
            return value_of(ex.args[Meta.isexpr(ex, :(=)) ? 2 : 1])
        elseif is_access(ex) || (is_reader(ex) && ex.args[3] isa QuoteNode)
            base = value_of(ex.args[is_access(ex) ? 1 : 2])
            name = is_access(ex) ? ex.args[2].value : ex.args[3].value
            return base === nothing || !(name isa Symbol) ? nothing : lookup(something(base), name)
        end
        return nothing
    end
    # A path's names, with the leading dots of a relative path dropped: the table
    # resolves `..SpaceAGORA` through the name `SpaceAGORA`.
    function import_path(ex)
        Meta.isexpr(ex, :.) && all(a -> a isa Symbol, ex.args) || return nothing
        start = findfirst(a -> a !== :., ex.args)
        return start === nothing ? nothing : ex.args[start:end]
    end
    function resolve_path(path, record)
        (isempty(path) || !haskey(aliases, path[1])) && return nothing
        value = Some(aliases[path[1]])
        for name in path[2:end]
            value = step(something(value), name, record)
            value === nothing && return nothing
        end
        return value
    end
    function path_value(path)
        (isempty(path) || !haskey(aliases, path[1])) && return nothing
        value = Some(aliases[path[1]])
        for name in path[2:end]
            value = lookup(something(value), name)
            value === nothing && return nothing
        end
        return value
    end
    binding_name(lhs) = lhs isa Symbol ? lhs :
        Meta.isexpr(lhs, :(::)) && length(lhs.args) == 2 && lhs.args[1] isa Symbol ? lhs.args[1] : nothing
    # One table for all scopes: a name bound to two different SpaceAGORA modules is
    # reported instead of being read as the first.
    ambiguous = Set{Symbol}()
    function bind!(name, value)
        name isa Symbol && value !== nothing && something(value) isa Module || return false
        if haskey(aliases, name)
            aliases[name] === something(value) || push!(ambiguous, name)
            return false
        end
        aliases[name] = something(value)
        return true
    end
    # Bindings, wherever they are, until no new one appears.
    function collect!(ex)
        ex isa Expr || return false
        changed = false
        if Meta.isexpr(ex, :(=)) || Meta.isexpr(ex, :kw)
            # Assignments and argument defaults.
            changed |= bind!(binding_name(ex.args[1]), value_of(ex.args[2]))
        elseif ex.head in (:import, :using)
            for arg in ex.args
                if Meta.isexpr(arg, :as)
                    path = import_path(arg.args[1])
                    path === nothing || (changed |= bind!(arg.args[2], path_value(path)))
                elseif Meta.isexpr(arg, :(:))
                    path = import_path(arg.args[1])
                    base = path === nothing ? nothing : path_value(path)
                    base === nothing && continue
                    for item in arg.args[2:end]
                        target, name = Meta.isexpr(item, :as) ? (item.args[1], item.args[2]) : (item, nothing)
                        member = import_path(target)
                        (member === nothing || length(member) != 1) && continue
                        changed |= bind!(something(name, only(member)), lookup(something(base), only(member)))
                    end
                else
                    path = import_path(arg)
                    path === nothing && continue
                    value = path_value(path)
                    changed |= bind!(last(path), value)
                    # A plain `using` also binds the module's exported names.
                    ex.head === :using && value !== nothing && something(value) isa Module || continue
                    M = something(value)
                    for n in names(M)
                        Base.isexported(M, n) && isdefined(M, n) || continue
                        changed |= bind!(n, Some(getglobal(M, n)))
                    end
                end
            end
        end
        for arg in ex.args
            changed |= collect!(arg)
        end
        return changed
    end
    while any(collect!, exprs)
    end
    base_of(ex) = is_access(ex) || is_assert(ex) ? base_of(ex.args[1]) : ex
    # The types asserted along a chain, down to a lookup's module argument: resolving
    # the chain skips them, so they are checked on their own.
    asserted_types(ex) = is_access(ex) ? asserted_types(ex.args[1]) :
        is_assert(ex) ? Any[ex.args[2]; asserted_types(ex.args[1])] :
        is_lookup(ex) ? asserted_types(ex.args[2]) : Any[]
    function check!(ex)
        if ex isa Symbol
            get(aliases, ex, nothing) === core && push!(found, "the $(nameof(core)) module used as a value ($ex)")
        elseif !(ex isa Expr)
            return
        elseif is_access(ex)
            # The whole chain, when it starts at a known module or at a lookup on one.
            base = base_of(ex)
            foreach(check!, asserted_types(ex))
            if base isa Symbol && haskey(aliases, base)
                resolve(ex, true)
            elseif is_reader(base) && resolve(base.args[2], false) !== nothing
                resolve(ex, true)
                foreach(check!, base.args[4:end])
            else
                check!(base)
            end
        elseif is_lookup(ex) && resolve(ex.args[2], false) !== nothing
            is_reader(ex) && resolve(ex, true)
            foreach(check!, asserted_types(ex))
            foreach(check!, ex.args[4:end])
        # A binding's right-hand side as one chain; a chained or type-asserted one goes
        # to the general case below, part by part.
        elseif Meta.isexpr(ex, :(=)) && binding_name(ex.args[1]) !== nothing && value_of(ex.args[2]) !== nothing &&
                !(Meta.isexpr(ex.args[2], :(=)) || Meta.isexpr(ex.args[2], :(::)))
            resolve(ex.args[2], true)
        elseif ex.head in (:import, :using)
            for arg in ex.args
                target = Meta.isexpr(arg, :as) ? arg.args[1] : arg
                if Meta.isexpr(target, :(:))
                    path = import_path(target.args[1])
                    value = path === nothing ? nothing : resolve_path(path, true)
                    value === nothing && continue
                    for item in target.args[2:end]
                        name = import_path(Meta.isexpr(item, :as) ? item.args[1] : item)
                        name === nothing || step(something(value), only(name), true)
                    end
                else
                    path = import_path(target)
                    path === nothing || resolve_path(path, true)
                end
            end
        else
            foreach(check!, ex.args)
        end
        return nothing
    end
    foreach(check!, exprs)
    foreach(name -> push!(found, "name bound to two SpaceAGORA modules: $name"), sort!(collect(ambiguous)))
    return found, used, aliases
end
end # module CoreAccess
