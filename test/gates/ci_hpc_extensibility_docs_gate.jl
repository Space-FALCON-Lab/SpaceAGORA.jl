const REPO_ROOT = dirname(dirname(dirname(@__FILE__)))

required_paths = [
    joinpath(REPO_ROOT, "docs", "src", "public_api_policy.md"),
    joinpath(REPO_ROOT, "docs", "src", "distributed_hpc.md"),
    joinpath(REPO_ROOT, "docs", "src", "extensibility.md"),
    joinpath(REPO_ROOT, "templates", "force_torque_model_template.jl"),
    joinpath(REPO_ROOT, "templates", "density_model_template.jl"),
    joinpath(REPO_ROOT, "templates", "control_hook_template.jl"),
    joinpath(REPO_ROOT, "templates", "actuator_template.jl"),
    joinpath(REPO_ROOT, "templates", "gnc_planner_template.jl"),
    joinpath(REPO_ROOT, "docs", "src", "user", "gnc_algorithm.md"),
    [joinpath(REPO_ROOT, "examples", env, name) for env in ("gnc_algorithm_env", "walkthrough_env") for name in
        ("README.md", "Project.template.toml", "setup.jl", "smoke.jl")]...,
    joinpath(REPO_ROOT, "templates", "starters", "README.md"),
    [joinpath(REPO_ROOT, "templates", "starters", kind, rel) for kind in ("project", "package", "companion")
        for rel in ("README.md", "Project.toml", "SpaceAGORASource.toml", joinpath("src", "StarterName.jl"),
            joinpath("test", "runtests.jl"), joinpath(".github", "workflows", "ci.yml"))]...,
    joinpath(REPO_ROOT, "templates", "starters", "companion", "test", "incompatible_core.jl"),
    joinpath(REPO_ROOT, "templates", "starters", "companion", "test", "load_order.jl"),
    joinpath(REPO_ROOT, "scripts", "new_contributor_repository.jl"),
    joinpath(REPO_ROOT, "scripts", "hpc", "local_process_runtime_analysis.sh"),
    joinpath(REPO_ROOT, "scripts", "hpc", "slurm_process_runtime_analysis.sh"),
]

missing = filter(path -> !isfile(path), required_paths)
isempty(missing) || error("Missing Wave 3/Wave 5 documentation assets: $(missing)")

hpc_doc = read(joinpath(REPO_ROOT, "docs", "src", "distributed_hpc.md"), String)
ext_doc = read(joinpath(REPO_ROOT, "docs", "src", "extensibility.md"), String)
api_doc = read(joinpath(REPO_ROOT, "docs", "src", "public_api_policy.md"), String)

for marker in (
    "SPACEAGORA_PERF_PROCS",
    "SPACEAGORA_PERF_WORKER_PROJECT",
    "SPACEAGORA_PARALLEL_POLICY_STATE_PATH",
    "shared depot",
    "scripts/hpc/slurm_process_runtime_analysis.sh",
)
    occursin(marker, hpc_doc) || error("distributed_hpc.md missing marker: $(marker)")
end

for marker in (
    "wrench",
    "getDensity",
    "update_control!",
    "control_wrench",
    "reference/extension_conventions.md",
    "templates/force_torque_model_template.jl",
    "templates/control_hook_template.jl",
    "templates/gnc_planner_template.jl",
    "user/gnc_algorithm.md",
)
    occursin(marker, ext_doc) || error("extensibility.md missing marker: $(marker)")
end

for marker in (
    "root `SpaceAGORA` module",
    "Public API",
    "internal",
    "generated public API",
)
    occursin(marker, api_doc) || error("public_api_policy.md missing marker: $(marker)")
end

# Docs pages whose code CI runs. Each script below runs in CI from a separate
# project, through that project's `smoke.jl`. Each region of a script, between a
# `# snippet: <name>` line and its `# end snippet` (regions may nest; marker lines
# belong to no region), must appear on the script's page as one Julia code block,
# character for character, so the page shows the code CI runs.
const SNIPPET_SCRIPTS = (
    joinpath("examples", "gnc_algorithm_env", "feedback_controller.jl") => joinpath("docs", "src", "user", "gnc_feedback.md"),
    joinpath("examples", "gnc_algorithm_env", "estimator.jl") => joinpath("docs", "src", "user", "gnc_feedback.md"),
    joinpath("examples", "gnc_algorithm_env", "rpo_planner.jl") => joinpath("docs", "src", "user", "gnc_planning.md"),
    joinpath("examples", "gnc_algorithm_env", "burn_planner.jl") => joinpath("docs", "src", "user", "gnc_planning.md"),
    joinpath("examples", "gnc_algorithm_env", "joint_planner.jl") => joinpath("docs", "src", "user", "gnc_planning.md"),
    joinpath("examples", "walkthrough_env", "my_effector.jl") => joinpath("docs", "src", "user", "custom_effector.md"),
    joinpath("examples", "walkthrough_env", "my_actuator.jl") => joinpath("docs", "src", "user", "custom_actuator.md"),
)

function snippet_regions(path)
    regions = Pair{String, String}[]
    open_regions = Pair{String, Vector{String}}[]
    for (number, line) in enumerate(eachline(path))
        start = match(r"^# snippet: (\S+)$", line)
        if start !== nothing
            push!(open_regions, String(start[1]) => String[])
        elseif line == "# end snippet"
            isempty(open_regions) && error("$(path):$(number): `# end snippet` without a start")
            name, lines = pop!(open_regions)
            push!(regions, name => join(lines, "\n"))
        else
            foreach(region -> push!(last(region), line), open_regions)
        end
    end
    isempty(open_regions) || error("$(path): snippet $(first(last(open_regions))) has no `# end snippet`")
    return regions
end

function julia_code_blocks(path)
    blocks, lines, inside = String[], String[], false
    for line in eachline(path)
        if !inside && line == "```julia"
            inside, lines = true, String[]
        elseif inside && line == "```"
            push!(blocks, join(lines, "\n"))
            inside = false
        elseif inside
            push!(lines, line)
        end
    end
    return blocks
end

snippet_dirs = unique(dirname(first(pair)) for pair in SNIPPET_SCRIPTS)
for dir in snippet_dirs
    listed = Set(basename(first(pair)) for pair in SNIPPET_SCRIPTS if dirname(first(pair)) == dir)
    present = Set(filter(name -> endswith(name, ".jl") && !(name in ("setup.jl", "smoke.jl")),
        readdir(joinpath(REPO_ROOT, dir))))
    present == listed || error("$(dir) holds the scripts $(sort!(collect(present))), but the gate lists " *
        "$(sort!(collect(listed))); list every script with the page that shows it.")
end
for (script, page) in SNIPPET_SCRIPTS
    runner = read(joinpath(REPO_ROOT, dirname(script), "smoke.jl"), String)
    occursin("\"$(basename(script))\"", runner) || error("$(dirname(script))/smoke.jl does not run $(basename(script))")
    regions = snippet_regions(joinpath(REPO_ROOT, script))
    isempty(regions) && error("$(script) has no snippet shown on $(page)")
    blocks = Set(julia_code_blocks(joinpath(REPO_ROOT, page)))
    for (name, text) in regions
        text in blocks || error("$(page) has no Julia code block equal to snippet `$(name)` of $(script); " *
            "change the page and the script together.")
    end
end

println("hpc_extensibility_docs_gate_ok")
