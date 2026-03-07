const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

required_files = (
    joinpath("src", "gnc", "guidance", "aerobraking", "interfaces.jl"),
    joinpath("src", "gnc", "guidance", "aerobraking", "dispatcher.jl"),
    joinpath("src", "gnc", "guidance", "aerobraking", "e_edg", "e_edg_strategy.jl"),
    joinpath("src", "gnc", "guidance", "aerobraking", "t_edg", "t_edg_strategy.jl"),
    joinpath("src", "gnc", "control", "aerobraking", "tracking_executor.jl"),
    joinpath("src", "mission", "operations", "aerobraking_policy", "policy_types.jl"),
)
for rel in required_files
    isfile(joinpath(REPO_ROOT, rel)) || error("Missing aerobraking boundary file: $rel")
end

isdir(joinpath(REPO_ROOT, "src", "gnc", "guidance", "targeting_control")) &&
    error("Retired guidance targeting_control path still exists.")

violations = String[]

# Guidance must not include control source files.
guidance_root = joinpath(REPO_ROOT, "src", "gnc", "guidance", "aerobraking")
for (root, _, files) in walkdir(guidance_root)
    for file in files
        endswith(file, ".jl") || continue
        path = joinpath(root, file)
        rel = relpath(path, REPO_ROOT)
        src = read(path, String)
        occursin(r"include\([^\)]*control", src) && push!(violations, "$rel: guidance includes control source directly")
    end
end

for root_rel in (
    joinpath("src", "gnc", "guidance"),
    joinpath("src", "gnc", "control"),
)
    root = joinpath(REPO_ROOT, root_rel)
    for (dir, _, files) in walkdir(root)
        for file in files
            endswith(file, ".jl") || continue
            path = joinpath(dir, file)
            rel = relpath(path, REPO_ROOT)
            src = read(path, String)
            occursin("DynamicEffectors", src) && push!(violations, "$rel: GNC source still depends on DynamicEffectors directly")
        end
    end
end

# Control must not own strategy/planning solver definitions.
for rel in (
    joinpath("src", "gnc", "control", "aerobraking", "tracking_executor.jl"),
    joinpath("src", "gnc", "control", "aerobraking", "control_commands.jl"),
    joinpath("src", "gnc", "control", "aerobraking", "constraint_tracking.jl"),
)
    src = read(joinpath(REPO_ROOT, rel), String)
    for forbidden_fn in (
        "function switch_calculation_with_integration",
        "function switch_calculation(",
        "function second_time_switch_recalc_with_integration",
        "function second_time_switch_recalc(",
        "function security_mode(",
        "function target_planning(",
    )
        occursin(forbidden_fn, src) && push!(violations, "$rel: contains strategy planning function `$forbidden_fn`")
    end
end

if !isempty(violations)
    println("gnc_aerobraking_boundary_violations:")
    for v in sort(unique(violations))
        println("  - $v")
    end
    error("GNC aerobraking boundary gate failed")
end

println("gnc_aerobraking_boundary_gate_ok")
