const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

const THIN_ENTRY_PATHS = [
    joinpath("src", "analysis", "verification", "telemetry_verification.jl"),
    joinpath("src", "parallel", "policy", "parallel_policy.jl"),
    joinpath("benchmarks", "studies", "performance_runtime_analysis.jl"),
    joinpath("benchmarks", "studies", "performance_smart_parallel_ladder.jl"),
    joinpath("benchmarks", "studies", "telemetry_hybrid_tuner.jl"),
    joinpath("benchmarks", "studies", "telemetry_odyssey_tuner.jl"),
]

function _noncomment_nonblank_lines(text::String)
    return [line for line in split(text, '\n') if begin
        stripped = strip(line)
        !isempty(stripped) && !startswith(stripped, '#')
    end]
end

for rel in THIN_ENTRY_PATHS
    path = joinpath(REPO_ROOT, rel)
    isfile(path) || error("Missing thin entry file: $rel")
    src = read(path, String)
    lines = _noncomment_nonblank_lines(src)
    length(lines) <= 80 || error("Thin entry file grew too large: $rel ($(length(lines)) non-comment lines)")

    occursin(r"^function\s"m, src) && error("Thin entry file must not own function definitions: $rel")
    occursin(r"^(Base\.@kwdef\s+)?mutable\s+struct\s"m, src) && error("Thin entry file must not own mutable struct definitions: $rel")
    occursin(r"^(Base\.@kwdef\s+)?struct\s"m, src) && error("Thin entry file must not own struct definitions: $rel")
end

println("thin_entry_files_gate_ok")
