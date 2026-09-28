const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SRC_ROOT = joinpath(REPO_ROOT, "src")

const CANONICAL_AGGREGATOR_FILES = Set([
    joinpath("src", "simulation", "engine", "simulation_engine.jl"),
    joinpath("src", "simulation", "callbacks", "callbacks.jl"),
    joinpath("src", "parallel", "routing", "parallel_profiles.jl"),
    joinpath("src", "parallel", "process", "parallel_process.jl"),
    joinpath("src", "core", "simulation_model.jl")
])

const STUB_PATTERN = r"\b[a-zA-Z0-9_]*_stub\(\)\s*=\s*nothing\b"
const FORBIDDEN_TOKENS = (
    "Compatibility wrapper: canonical path forwarding to legacy implementation.",
    "No-op ",
)

function _contains_behavior_ast(ex)::Bool
    ex isa Expr || return false

    if ex.head == :function || ex.head == :struct
        return true
    end

    # Short-form function definition: f(x) = ...
    if ex.head == :(=)
        lhs = ex.args[1]
        if lhs isa Expr && lhs.head == :call
            return true
        end
    end

    for arg in ex.args
        _contains_behavior_ast(arg) && return true
    end
    return false
end

function _aggregator_has_behavior(src::String)::Bool
    wrapped = "begin\n" * src * "\nend"
    parsed = try
        Meta.parse(wrapped)
    catch
        return true
    end
    return _contains_behavior_ast(parsed)
end

violations = String[]
for (root, _, files) in walkdir(SRC_ROOT)
    for file in files
        endswith(file, ".jl") || continue
        path = joinpath(root, file)
        rel = relpath(path, REPO_ROOT)
        src = read(path, String)

        occursin(STUB_PATTERN, src) &&
            push!(violations, "$rel: contains forbidden silent stub pattern '*_stub() = nothing'.")

        for token in FORBIDDEN_TOKENS
            occursin(token, src) || continue
            push!(violations, "$rel: contains forbidden token '$token'.")
        end

        if rel in CANONICAL_AGGREGATOR_FILES && !occursin("# Canonical aggregator: no behavior ownership.", src)
            push!(violations, "$rel: canonical aggregator must declare standard header comment.")
        end
        if rel in CANONICAL_AGGREGATOR_FILES && _aggregator_has_behavior(src)
            push!(violations, "$rel: canonical aggregator contains behavior (function/type definitions detected).")
        end
    end
end

if !isempty(violations)
    println("src_completeness_contract_violations:")
    for v in sort(unique(violations))
        println("  - $v")
    end
    error("Source completeness contract gate failed")
end

println("src_completeness_contract_gate_ok")
