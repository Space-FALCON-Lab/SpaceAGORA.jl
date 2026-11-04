const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function _active_lines(path::String)
    lines = Tuple{Int, String}[]
    for (line_no, raw_line) in enumerate(eachline(path))
        active_line = strip(first(split(raw_line, '#', limit=2)))
        isempty(active_line) && continue
        push!(lines, (line_no, active_line))
    end
    return lines
end

const HOTPATH_FILES = String[
    joinpath(REPO_ROOT, "src", "orchestration", "simulation", "state", "runtime_types.jl"),
    joinpath(REPO_ROOT, "src", "orchestration", "simulation", "engine", "setup.jl"),
    joinpath(REPO_ROOT, "src", "models", "dynamics", "coupled", "aerodynamic_wrench_models.jl"),
    joinpath(REPO_ROOT, "src", "models", "dynamics", "coupled", "perturbations.jl"),
]

# Split source keeps exactly the original hot-path scan coverage. New children
# are covered recursively so later source splits cannot bypass the contract.
const SPLIT_HOTPATH_DIR = joinpath(REPO_ROOT, "src", "models", "dynamics", "coupled", "perturbations")
const MULTIBODY_POLICY_FILE = joinpath(REPO_ROOT, "src", "orchestration", "parallel", "policy", "multibody_routing.jl")

const CALLBACK_DIR = joinpath(REPO_ROOT, "src", "orchestration", "simulation", "callbacks")
const FORBIDDEN_VECTOR_ANY = r"\bVector\{Any\}"
const FORBIDDEN_DICT_SYMBOL_ANY = r"\bDict\{Symbol,\s*Any\}"
const ALLOWED_SAVEDATA_LINE = r"^const\s+SaveData\s*=\s*Dict\{Symbol,\s*Any\}$"

# These existing heterogeneous boundaries carry provider-owned native/runtime
# objects or arbitrary recorder field values. Keep their runtime types intact.
# Each exception is one exact active statement at one exact owner, not a file
# exemption. Missing entries and additional occurrences both fail the gate.
const EXISTING_HETEROGENEOUS_BOUNDARIES = (
    (path="src/orchestration/simulation/state/runtime_types.jl",
     statement="walk_models::Vector{Any}", expected_count=1,
     reason="Optional GRAM providers own the concrete per-spacecraft walk types."),
    (path="src/orchestration/simulation/state/runtime_types.jl",
     statement="external_runtimes::Vector{Any} = Any[]", expected_count=1,
     reason="Independent external propagators own the concrete runtime types."),
    (path="src/orchestration/simulation/callbacks/event_callbacks.jl",
     statement="storage::Dict{Symbol, Any}", expected_count=1,
     reason="Recorder fields can have different scalar, vector or structured storage types."),
    (path="src/orchestration/simulation/callbacks/event_callbacks.jl",
     statement="fallback::Dict{Symbol, Vector{Any}}", expected_count=1,
     reason="User-defined recorder fields use heterogeneous fallback value histories."),
    (path="src/orchestration/simulation/callbacks/event_callbacks.jl",
     statement="storage = Dict{Symbol, Any}()", expected_count=1,
     reason="Constructs the recorder's existing heterogeneous field storage."),
    (path="src/orchestration/simulation/callbacks/event_callbacks.jl",
     statement="fallback = Dict{Symbol, Vector{Any}}()", expected_count=1,
     reason="Constructs the recorder's existing fallback-history collection."),
    (path="src/orchestration/simulation/callbacks/event_callbacks.jl",
     statement="fallback[field.name] = Vector{Any}(undef, cap)", expected_count=1,
     reason="Allocates the existing user-defined recorder field history."),
    (path="src/orchestration/simulation/callbacks/event_callbacks.jl",
     statement="grown = Vector{Any}(undef, new_cap)", expected_count=1,
     reason="Grows the same fallback recorder history without changing its value contract."),
)
const HETEROGENEOUS_BOUNDARY_COUNTS = Dict(
    (entry.path, entry.statement) => entry.expected_count for entry in EXISTING_HETEROGENEOUS_BOUNDARIES
)
length(HETEROGENEOUS_BOUNDARY_COUNTS) == length(EXISTING_HETEROGENEOUS_BOUNDARIES) ||
    error("Duplicate heterogeneous boundary inventory entry")
baseline_seen = Dict(key => 0 for key in keys(HETEROGENEOUS_BOUNDARY_COUNTS))

violations = String[]
targets = [HOTPATH_FILES; MULTIBODY_POLICY_FILE]

for scan_root in (CALLBACK_DIR, SPLIT_HOTPATH_DIR)
    isdir(scan_root) || error("Missing runtime hot-path scan root: $(relpath(scan_root, REPO_ROOT))")
    for (root, _, files) in walkdir(scan_root)
        for file in files
            endswith(file, ".jl") || continue
            push!(targets, joinpath(root, file))
        end
    end
end

for path in sort(unique(targets))
    isfile(path) || error("Missing runtime Any hot-path target: $(relpath(path, REPO_ROOT))")
    rel = replace(relpath(path, REPO_ROOT), '\\' => '/')
    for (line_no, active_line) in _active_lines(path)
        vector_any = occursin(FORBIDDEN_VECTOR_ANY, active_line)
        dict_any = occursin(FORBIDDEN_DICT_SYMBOL_ANY, active_line) && !occursin(ALLOWED_SAVEDATA_LINE, active_line)
        (vector_any || dict_any) || continue
        key = (rel, active_line)
        if haskey(HETEROGENEOUS_BOUNDARY_COUNTS, key)
            baseline_seen[key] += 1
            continue
        end
        vector_any && push!(violations, "$rel:$line_no contains forbidden Vector{Any}")
        dict_any && push!(violations, "$rel:$line_no contains forbidden Dict{Symbol, Any}")
    end
end

for ((relative, statement), expected_count) in HETEROGENEOUS_BOUNDARY_COUNTS
    observed = baseline_seen[(relative, statement)]
    observed == expected_count || push!(violations,
        "$relative: heterogeneous boundary '$statement' expected $expected_count occurrence(s), observed $observed")
end

if !isempty(violations)
    println("runtime_any_hotpath_violations:")
    for violation in sort(unique(violations))
        println("  - $violation")
    end
    error("Runtime Any hot-path gate failed")
end

println("runtime_any_hotpath_gate_ok")
