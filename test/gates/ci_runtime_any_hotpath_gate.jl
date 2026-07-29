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
    joinpath(REPO_ROOT, "src", "core", "types", "runtime_types.jl"),
    joinpath(REPO_ROOT, "src", "simulation", "engine", "setup.jl"),
    joinpath(REPO_ROOT, "src", "dynamics", "coupled", "aerodynamic_wrench_models.jl"),
    joinpath(REPO_ROOT, "src", "dynamics", "coupled", "perturbations.jl"),
]

const CALLBACK_DIR = joinpath(REPO_ROOT, "src", "simulation", "callbacks")
const FORBIDDEN_VECTOR_ANY = r"\bVector\{Any\}\b"
const FORBIDDEN_DICT_SYMBOL_ANY = r"\bDict\{Symbol,\s*Any\}\b"
const ALLOWED_SAVEDATA_LINE = r"^const\s+SaveData\s*=\s*Dict\{Symbol,\s*Any\}$"

violations = String[]
targets = copy(HOTPATH_FILES)

if isdir(CALLBACK_DIR)
    for (root, _, files) in walkdir(CALLBACK_DIR)
        for file in files
            endswith(file, ".jl") || continue
            push!(targets, joinpath(root, file))
        end
    end
end

for path in sort(unique(targets))
    isfile(path) || error("Missing runtime Any hot-path target: $(relpath(path, REPO_ROOT))")
    rel = relpath(path, REPO_ROOT)
    for (line_no, active_line) in _active_lines(path)
        occursin(FORBIDDEN_VECTOR_ANY, active_line) &&
            push!(violations, "$rel:$line_no contains forbidden Vector{Any}")
        if occursin(FORBIDDEN_DICT_SYMBOL_ANY, active_line) && !occursin(ALLOWED_SAVEDATA_LINE, active_line)
            push!(violations, "$rel:$line_no contains forbidden Dict{Symbol, Any}")
        end
    end
end

if !isempty(violations)
    println("runtime_any_hotpath_violations:")
    for violation in sort(unique(violations))
        println("  - $violation")
    end
    error("Runtime Any hot-path gate failed")
end

println("runtime_any_hotpath_gate_ok")
