const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const SRP_SENSITIVE_FILES = (
    joinpath("src", "gnc", "control", "aerobraking", "control_commands.jl"),
    joinpath("src", "gnc", "control", "aerobraking", "constraint_tracking.jl"),
    joinpath("src", "gnc", "guidance", "aerobraking", "t_edg", "eom_predictor.jl"),
    joinpath("src", "gnc", "guidance", "aerobraking", "t_edg", "trajectory_predictor.jl")
)

function _non_comment_lines(src::String)
    out = String[]
    for raw in split(src, '\n')
        line = strip(first(split(raw, '#', limit=2)))
        isempty(line) && continue
        push!(out, line)
    end
    return out
end

violations = String[]
for rel in SRP_SENSITIVE_FILES
    path = joinpath(REPO_ROOT, rel)
    isfile(path) || error("Missing SRP-sensitive file: $rel")

    lines = _non_comment_lines(read(path, String))
    srp_tokens = [line for line in lines if occursin(r"\bsrp_ii\b", line)]
    has_assignment = any(occursin(r"\bsrp_ii\s*=", line) for line in srp_tokens)
    has_usage = length(srp_tokens) > (has_assignment ? 1 : 0)

    if has_assignment && !has_usage
        push!(violations, "$rel: srp_ii is assigned but never used in active force accumulation.")
    end
end

if !isempty(violations)
    println("inert_force_term_violations:")
    for v in sort(unique(violations))
        println("  - $v")
    end
    error("Inert force-term gate failed")
end

println("no_inert_force_terms_gate_ok")
