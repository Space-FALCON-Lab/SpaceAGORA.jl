using SpaceAGORA
using Markdown

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# Checks the extension hook registry, docs/extension_hooks.jl, against the loaded
# package and the docs tree:
# 1. every row is an exported function listed in docs/public_api_symbols.jl;
# 2. every exported function in the hook sections is a row or in NOT_HOOKS;
# 3. the signature cell is the first line of the docstring's signature block;
# 4. the signature's arity matches a method, and the fallback cell matches the
#    method table at each documented arity, counting as fallbacks only methods
#    that every model of the row's `base` type reaches;
# 5. the row is deprecated exactly when its docstring opens with "Deprecated"; a
#    deprecated row names a replacement row and is removed in REMOVAL_RELEASE;
# 6. the row's page exists under docs/src and names the hook;
# 7. current rows name the spacecraft index `spacecraft_index`, never `i`,
#    `sat_idx` or `idx`.
# KNOWN_GAPS lists the rows exempt from rules 3, 6 or 7 today; an entry that no
# longer applies fails the gate. A hook with no methods yet has no arity check.
# The gate reads method tables, not behaviour: what a fallback returns is for the
# hook's own tests.

# The registry and the public API list load into a module of their own, so that
# the public API gate's copy in the same process is left alone.
_hook_gate_inputs = Module(:ExtensionHookReferenceGateInputs)
Base.include(_hook_gate_inputs, joinpath(REPO_ROOT, "docs", "public_api_symbols.jl"))
Base.include(_hook_gate_inputs, joinpath(REPO_ROOT, "docs", "extension_hooks.jl"))

const _HOOK_SHAPES = ("T", "R", "Q", "D", "A", "H", "P", "L")
const _HOOK_GAP_RULES = (:signature_line, :index_word, :page)
const _HOOK_OLD_INDEX_WORDS = (:i, :sat_idx, :idx)
# The argument a hook dispatches on is the first one with one of these names.
# `link` is the shared object of a GNC link (`check_link`, `link_initialize!`);
# `producer` the model that writes one (`before_reference_control!`).
const _HOOK_DISPATCH_NAMES = (:model, :effector, :controller, :actuator, :sensor,
    :planner, :algorithm, :family, :binding, :result, :provider, :link, :producer)

violations = String[]

# The binding's first docstring, as Julia and the API pages show it first.
function _hook_first_docstring(owner::Module, symbol::Symbol)
    doc = Docs.doc(Docs.Binding(owner, symbol))
    while doc isa Markdown.MD && !isempty(doc.content) && first(doc.content) isa Markdown.MD
        doc = first(doc.content)
    end
    return doc
end

# The first line of the signature block that opens the first docstring, or
# `nothing` when the docstring does not open with one.
function _hook_signature_line(owner::Module, symbol::Symbol)
    doc = _hook_first_docstring(owner, symbol)
    isempty(doc.content) && return nothing
    block = first(doc.content)
    block isa Markdown.Code || return nothing
    return String(strip(first(split(block.code, '\n'))))
end

function _hook_first_paragraph(owner::Module, symbol::Symbol)
    for block in _hook_first_docstring(owner, symbol).content
        block isa Markdown.Paragraph && return strip(Markdown.plain(block))
    end
    return ""
end

# Whether a docstring paragraph announces a deprecation ("Deprecated", or
# "**Deprecated**").
_hook_says_deprecated(paragraph::AbstractString) = startswith(lstrip(paragraph, ['*', '_']), "Deprecated")

function _hook_argument_name(arg)
    arg isa Symbol && return arg
    arg isa Expr || return nothing
    arg.head === :(::) && return length(arg.args) == 2 ? _hook_argument_name(arg.args[1]) : nothing
    arg.head in (:kw, :...) && return _hook_argument_name(arg.args[1])
    return nothing
end

# The name, positional argument names and arity range of a signature cell.
# Optional trailing arguments are written `[, p]`; text after ` -> ` is the result.
function _hook_parse_signature(signature::AbstractString)
    call = strip(first(split(signature, " -> "; limit=2)))
    required = Meta.parse(replace(call, r"\[,[^\[\]]*\]" => ""))
    full = Meta.parse(replace(call, r"\[,\s*([^\[\]]*)\]" => s", \1"))
    for expr in (required, full)
        expr isa Expr && expr.head === :call || error("not a call: $(signature)")
    end
    positional(expr) = [a for a in expr.args[2:end] if !(a isa Expr && a.head === :parameters)]
    defaulted(args) = count(a -> a isa Expr && a.head === :kw, args)
    required_args, full_args = positional(required), positional(full)
    names = [_hook_argument_name(a) for a in full_args]
    return (name = full.args[1], names,
            min_arity = length(required_args) - defaulted(required_args),
            max_arity = length(full_args))
end

_hook_positional_count(m::Method) = m.nargs - 1

_hook_accepts(m::Method, arity::Int) =
    m.isva ? arity >= _hook_positional_count(m) - 1 : arity == _hook_positional_count(m)

function _hook_dispatch_type(m::Method, position::Int)
    params = Base.unwrap_unionall(m.sig).parameters
    position + 1 <= length(params) || return nothing
    T = params[position + 1]
    T isa Core.TypeofVararg && (T = T.T)
    T isa TypeVar && (T = T.ub)
    return T
end

# The row's contract type: `:Any`, or a public abstract type of SpaceAGORA.
function _hook_base(name::Symbol)
    name === :Any && return Any
    name in names(SpaceAGORA) && isdefined(SpaceAGORA, name) || return nothing
    T = getproperty(SpaceAGORA, name)
    return T isa Type && isabstracttype(Base.unwrap_unionall(T)) ? T : nothing
end

# A fallback is a method every model of the row's base type reaches: its dispatch
# argument is untyped or typed with the base or a supertype of it. A method typed
# with some other abstract type (an RPO planner type on a force hook, say) is not.
_hook_is_fallback(T, base) = T === Any || (T isa Type && base <: T)

_hook_fallback_phrase(name, base) = base === Any ?
    "no method with an untyped `$(name)`" :
    "no method with `$(name)` untyped or typed with $(nameof(base)) or a supertype"

_hook_mentions(text::AbstractString, symbol::Symbol) =
    occursin(Regex("(?<![A-Za-z0-9_])\\Q$(symbol)\\E(?![A-Za-z0-9_!])"), text)

_hook_unverified_arity = let hook_api = _hook_gate_inputs.PublicAPISymbols,
        hook_registry = _hook_gate_inputs.ExtensionHooks
    rows = hook_registry.extension_hook_rows()
    row_symbols = [row.symbol for row in rows]
    not_hook_symbols = [entry.symbol for entry in hook_registry.NOT_HOOKS]
    gaps = Set((gap.symbol, gap.rule) for gap in hook_registry.KNOWN_GAPS)
    has_gap(symbol, rule) = (symbol, rule) in gaps
    docs_src = joinpath(REPO_ROOT, "docs", "src")

    specs = hook_api.public_api_specs(SpaceAGORA)
    spec_by_symbol = Dict(spec.symbol => spec for spec in specs)
    section_titles = Set(section.title for section in hook_api.PUBLIC_API_SECTIONS)

    for symbol in unique(filter(s -> count(==(s), row_symbols) > 1, row_symbols))
        push!(violations, "$(symbol): more than one registry row.")
    end
    for symbol in intersect(row_symbols, not_hook_symbols)
        push!(violations, "$(symbol): both a registry row and in NOT_HOOKS.")
    end
    for title in hook_registry.HOOK_SECTIONS
        title in section_titles || push!(violations, "HOOK_SECTIONS names \"$(title)\", which is not a section of docs/public_api_symbols.jl.")
    end

    # Rule 2: every exported function in the hook sections is accounted for.
    for spec in specs
        spec.section in hook_registry.HOOK_SECTIONS || continue
        getproperty(spec.owner_module, spec.symbol) isa Function || continue
        spec.symbol in row_symbols || spec.symbol in not_hook_symbols ||
            push!(violations, "$(spec.symbol) (\"$(spec.section)\"): exported function with no registry row and no NOT_HOOKS reason.")
    end
    for entry in hook_registry.NOT_HOOKS
        spec = get(spec_by_symbol, entry.symbol, nothing)
        if spec === nothing || !(spec.section in hook_registry.HOOK_SECTIONS) ||
                !(getproperty(spec.owner_module, entry.symbol) isa Function)
            push!(violations, "NOT_HOOKS lists $(entry.symbol), which is not an exported function in a hook section.")
        end
        isempty(strip(entry.reason)) && push!(violations, "NOT_HOOKS entry $(entry.symbol) gives no reason.")
    end

    for gap in hook_registry.KNOWN_GAPS
        gap.rule in _HOOK_GAP_RULES || push!(violations, "KNOWN_GAPS entry $(gap.symbol): unknown rule $(repr(gap.rule)).")
        gap.symbol in row_symbols || push!(violations, "KNOWN_GAPS entry $(gap.symbol) is not a registry row.")
        isempty(strip(gap.reason)) && push!(violations, "KNOWN_GAPS entry $(gap.symbol) gives no reason.")
    end

    unverified_arity = Symbol[]

    for group in hook_registry.EXTENSION_HOOK_GROUPS, row in group.rows
        s = row.symbol
        isempty(strip(group.title)) && push!(violations, "$(s): its group has no title.")
        base_name = row.base === nothing ? (hasproperty(group, :base) ? group.base : nothing) : row.base
        base = base_name isa Symbol ? _hook_base(base_name) : nothing
        if base === nothing
            push!(violations, "$(s): base $(repr(base_name)) is not :Any or a public abstract type of SpaceAGORA (set `base` on the group or the row).")
            base = Any
        end
        row.shape in _HOOK_SHAPES || push!(violations, "$(s): shape $(repr(row.shape)) is not one of $(join(_HOOK_SHAPES, ", ")).")
        for field in (:called, :arguments)
            isempty(strip(getproperty(row, field))) && push!(violations, "$(s): empty $(field) cell.")
        end

        # Rule 1: exported and listed.
        spec = get(spec_by_symbol, s, nothing)
        if !(isdefined(SpaceAGORA, s) && Base.isexported(SpaceAGORA, s))
            push!(violations, "$(s): not exported from SpaceAGORA.")
            continue
        end
        if spec === nothing
            push!(violations, "$(s): not listed in docs/public_api_symbols.jl.")
            continue
        end
        spec.section in hook_registry.HOOK_SECTIONS ||
            push!(violations, "$(s): listed in \"$(spec.section)\", which is not in HOOK_SECTIONS.")
        f = getproperty(SpaceAGORA, s)
        if !(f isa Function)
            push!(violations, "$(s): not a function.")
            continue
        end

        parsed = try
            _hook_parse_signature(row.signature)
        catch err
            push!(violations, "$(s): signature cell does not parse ($(sprint(showerror, err))).")
            continue
        end
        parsed.name === s || push!(violations, "$(s): signature cell names $(parsed.name).")

        # Rule 3: the signature is the docstring's.
        docline = _hook_signature_line(SpaceAGORA, s)
        if has_gap(s, :signature_line)
            docline === nothing ||
                push!(violations, "$(s): the docstring now opens with a signature line; remove its :signature_line entry from KNOWN_GAPS and use that line.")
        elseif docline === nothing
            push!(violations, "$(s): the docstring does not open with a signature line.")
        elseif docline != row.signature
            push!(violations, "$(s): signature cell \"$(row.signature)\" differs from the docstring's first line \"$(docline)\".")
        end

        # Rule 4: arity and fallback against the method table. A fallback must take
        # each documented arity, so that a type with no method of its own, called as
        # documented, reaches it.
        ms = collect(methods(f))
        arities = parsed.min_arity:parsed.max_arity
        if isempty(ms)
            push!(unverified_arity, s)
        else
            for arity in arities
                any(m -> _hook_accepts(m, arity), ms) ||
                    push!(violations, "$(s): no method takes $(arity) positional arguments, as the signature cell does.")
            end
        end
        position = findfirst(in(_HOOK_DISPATCH_NAMES), parsed.names)
        if position === nothing
            push!(violations, "$(s): the signature names no dispatch argument ($(join(_HOOK_DISPATCH_NAMES, ", "))).")
        else
            fallbacks = filter(ms) do m
                T = _hook_dispatch_type(m, position)
                T !== nothing && _hook_is_fallback(T, base)
            end
            # Only fallbacks a documented call can reach count against `none`.
            documented = filter(m -> any(a -> _hook_accepts(m, a), arities), fallbacks)
            if row.fallback == "none"
                isempty(documented) ||
                    push!(violations, "$(s): fallback cell says none, but $(first(documented)) is a fallback method.")
            elseif startswith(row.fallback, "returns ") || startswith(row.fallback, "throws ")
                for arity in arities
                    any(m -> _hook_accepts(m, arity), fallbacks) ||
                        push!(violations, "$(s): fallback cell says \"$(row.fallback)\", but $(_hook_fallback_phrase(parsed.names[position], base)) takes $(arity) positional arguments.")
                end
            else
                push!(violations, "$(s): fallback cell must be none, \"returns ...\" or \"throws ...\".")
            end
        end

        # Rule 5: the status agrees with the docstring, both ways, and the
        # deprecation fields are consistent.
        paragraph = _hook_first_paragraph(SpaceAGORA, s)
        if row.status == "current"
            row.replacement === nothing && row.removed_in === nothing ||
                push!(violations, "$(s): a current row has no replacement or removed_in.")
            _hook_says_deprecated(paragraph) &&
                push!(violations, "$(s): the docstring's first paragraph opens with \"Deprecated\", but the row says current.")
        elseif row.status == "deprecated"
            replacement = row.replacement
            if replacement === nothing || !(replacement in row_symbols) || replacement === s
                push!(violations, "$(s): a deprecated row must name a replacement that is another registry row.")
            end
            row.removed_in == hook_registry.REMOVAL_RELEASE ||
                push!(violations, "$(s): removed_in must be $(hook_registry.REMOVAL_RELEASE).")
            _hook_says_deprecated(paragraph) ||
                push!(violations, "$(s): the docstring's first paragraph does not open with \"Deprecated\".")
            replacement === nothing || _hook_mentions(paragraph, replacement) ||
                push!(violations, "$(s): the docstring's first paragraph does not name the replacement $(replacement).")
            occursin(hook_registry.REMOVAL_RELEASE, Markdown.plain(_hook_first_docstring(SpaceAGORA, s))) ||
                push!(violations, "$(s): the docstring does not give the removal release $(hook_registry.REMOVAL_RELEASE).")
        else
            push!(violations, "$(s): status must be current or deprecated.")
        end

        # Rule 6: the page exists and names the hook.
        if row.page === nothing
            has_gap(s, :page) || push!(violations, "$(s): no page; name the page that explains it.")
        else
            has_gap(s, :page) &&
                push!(violations, "$(s): has a page now; remove its :page entry from KNOWN_GAPS.")
            path = joinpath(docs_src, row.page)
            if startswith(row.page, "generated/") || !endswith(row.page, ".md") || !isfile(path)
                push!(violations, "$(s): page docs/src/$(row.page) does not exist as a hand-written page.")
            elseif !_hook_mentions(read(path, String), s)
                push!(violations, "$(s): page docs/src/$(row.page) does not name the hook.")
            end
        end

        # Rule 7: the index word, in current rows.
        old_words = intersect(parsed.names, _HOOK_OLD_INDEX_WORDS)
        if row.status == "current"
            if has_gap(s, :index_word)
                isempty(old_words) &&
                    push!(violations, "$(s): the signature now uses spacecraft_index; remove its :index_word entry from KNOWN_GAPS.")
            elseif !isempty(old_words)
                push!(violations, "$(s): a current row names the spacecraft index $(join(old_words, ", ")); use spacecraft_index.")
            end
        elseif has_gap(s, :index_word)
            push!(violations, "$(s): deprecated rows need no :index_word entry in KNOWN_GAPS.")
        end
    end

    for facade in hook_registry.VERSIONED_FACADES
        s = facade.symbol
        value = isdefined(SpaceAGORA, s) ? getproperty(SpaceAGORA, s) : nothing
        if !(value isa Module && isdefined(value, :CONTRACT_VERSION) && value.CONTRACT_VERSION isa VersionNumber)
            push!(violations, "$(s): not a SpaceAGORA facade module with a CONTRACT_VERSION.")
        end
        path = joinpath(docs_src, facade.page)
        isfile(path) && _hook_mentions(read(path, String), s) ||
            push!(violations, "$(s): page docs/src/$(facade.page) does not exist or does not name the facade.")
    end
    unverified_arity
end

if !isempty(violations)
    println("extension_hook_reference_violations:")
    for v in violations
        println("  - ", v)
    end
    error("Extension hook reference gate failed")
end

isempty(_hook_unverified_arity) ||
    println("extension_hook_reference_no_methods_yet: ", join(_hook_unverified_arity, ", "), " (arity unchecked; fallback none)")
let registry = _hook_gate_inputs.ExtensionHooks
    println("extension_hook_reference_gate_ok rows=$(length(registry.extension_hook_rows())) ",
            "not_hooks=$(length(registry.NOT_HOOKS)) known_gaps=$(length(registry.KNOWN_GAPS)) ",
            "facades=$(length(registry.VERSIONED_FACADES))")
end
