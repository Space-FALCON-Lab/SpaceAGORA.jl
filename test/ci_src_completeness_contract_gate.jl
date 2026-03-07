const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const SRC_ROOT = joinpath(REPO_ROOT, "src")

const SCAFFOLD_INTERFACE_FILES = Set([
    joinpath("src", "vehicle", "resources", "resources.jl"),
    joinpath("src", "vehicle", "resources", "battery", "battery_model.jl"),
    joinpath("src", "vehicle", "resources", "solar_array", "solar_array_model.jl"),
    joinpath("src", "vehicle", "resources", "power_bus", "power_bus_model.jl"),
    joinpath("src", "vehicle", "resources", "loads", "load_model.jl"),
    joinpath("src", "mission", "constellation", "network", "network.jl"),
    joinpath("src", "gnc", "estimation", "estimation.jl"),
    joinpath("src", "vehicle", "actuators", "laser_terminal", "laser_terminal.jl")
])

const CANONICAL_AGGREGATOR_FILES = Set([
    joinpath("src", "simulation", "engine", "simulation_engine.jl"),
    joinpath("src", "simulation", "callbacks", "callbacks.jl"),
    joinpath("src", "parallel", "routing", "parallel_profiles.jl"),
    joinpath("src", "core", "simulation_model.jl")
])
const RETIRED_SOURCE_DIRS = Set([
    joinpath("src", "examples"),
    joinpath("src", "analysis", "plotting"),
    joinpath("src", "analysis", "reports"),
    joinpath("src", "parallel", "state"),
    joinpath("src", "parallel", "tuning"),
    joinpath("src", "vehicle", "sensors"),
])
const RETIRED_SOURCE_FILES = Set([
    joinpath("src", "core", "utils", "typed_example_utils.jl"),
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

        if rel in SCAFFOLD_INTERFACE_FILES && !occursin("Not implemented: ", src)
            push!(violations, "$rel: scaffold interface must throw explicit 'Not implemented:' errors.")
        end

        if rel in CANONICAL_AGGREGATOR_FILES && !occursin("# Canonical aggregator: no behavior ownership.", src)
            push!(violations, "$rel: canonical aggregator must declare standard header comment.")
        end
        if rel in CANONICAL_AGGREGATOR_FILES && _aggregator_has_behavior(src)
            push!(violations, "$rel: canonical aggregator contains behavior (function/type definitions detected).")
        end
    end
end

for rel in RETIRED_SOURCE_DIRS
    isdir(joinpath(REPO_ROOT, rel)) &&
        push!(violations, "$rel: retired source directory must not exist.")
end

for rel in RETIRED_SOURCE_FILES
    isfile(joinpath(REPO_ROOT, rel)) &&
        push!(violations, "$rel: retired source file must not exist.")
end

if !isempty(violations)
    println("src_completeness_contract_violations:")
    for v in sort(unique(violations))
        println("  - $v")
    end
    error("Source completeness contract gate failed")
end

println("src_completeness_contract_gate_ok")
