using SpaceAGORA

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const EXAMPLES_DIR = joinpath(REPO_ROOT, "examples")
# The examples run in the examples environment (examples/setup.jl prepares it).
const PROJECT_PATH = joinpath(REPO_ROOT, "examples")

const HYPR_EXAMPLES = Set([
    "Earth_RPO_CubeSat_MPC.jl", "Earth_RPO_CubeSat_MPC_Batch.jl",
    "Earth_RPO_CubeSat_MPC_PlannerComparison.jl", "Earth_RPO_CubeSat_MPC_Replanning.jl",
    "Robot_Arm_Planner_Cloth_Demo.jl",
])
const EDG_EXAMPLES = Set([
    "AGORA_Odyssey_Energy_Depletion_Control_Test.jl", "odyssey_surrogate.jl",
])
function example_project(path)
    name = basename(path)
    if name in EDG_EXAMPLES
        # Reuse the CLI's explicit/active paired-project validation. The child
        # loads EDG before common.jl activates core, while LOAD_PATH stays isolated.
        return SpaceAGORA.SpaceAGORACLI._edg_project()
    end
    return name in HYPR_EXAMPLES ?
        joinpath(REPO_ROOT, "examples", "rpo_planner_env") : PROJECT_PATH
end
function example_command(path)
    project = example_project(path)
    return `$(Base.julia_cmd()) --startup-file=no --compiled-modules=existing --depwarn=error --project=$(project) $(path)`
end

function list_examples()
    # setup.jl prepares the examples environment; it is not an example.
    helper_files = Set(["common.jl", "setup.jl"])
    # Odyssey uses its pinned example environment, not this repository environment.
    # surrogate-user.yml runs its full native-free passage and verifies the outputs.
    dedicated_environment_examples = Set(["odyssey_surrogate.jl"])
    files = sort(filter(f -> endswith(f, ".jl"), readdir(EXAMPLES_DIR; join=true)))
    files = filter(f -> !(basename(f) in helper_files || basename(f) in dedicated_environment_examples), files)
    token = strip(get(ENV, "SPACEAGORA_EXAMPLE_FILTER", ""))
    if !isempty(token)
        files = filter(f -> occursin(token, basename(f)), files)
    end
    return files
end

function run_example(example_path::String)
    output = IOBuffer()
    ok = false
    has_nan = false
    text = ""

    # Isolate each example's transient outputs so smoke runs do not share cwd state.
    mktempdir() do tmp
        cmd = example_command(example_path)
        cmd = Cmd(cmd; dir=tmp)
        cmd = addenv(
            cmd,
            "JULIA_LOAD_PATH" => "@:@stdlib",
            "SPACEAGORA_EXAMPLE_SMOKE" => "1",
            "SPACEAGORA_EXAMPLE_SMOKE_RESULTS" => "1",
            "SPACEAGORA_EXAMPLE_SMOKE_MISSION_TIME" => "120.0",
            "SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0",
            "SPACEAGORA_WARN_NORMALIZE" => "0"
        )
        proc = run(pipeline(ignorestatus(cmd), stdout=output, stderr=output))
        text = String(take!(output))
        has_nan = occursin("init_NaN", text) || occursin("First function call produced NaNs", text)
        ok = success(proc)
    end

    return ok, has_nan, text
end

examples = list_examples()
println("Running $(length(examples)) examples in smoke mode...")

failures = String[]
for (idx, example) in enumerate(examples)
    name = basename(example)
    print("[$idx/$(length(examples))] $name ... ")
    ok, has_nan, text = run_example(example)
    if ok && !has_nan
        println("ok")
    else
        println("FAILED")
        push!(failures, name)
        println("----- begin $name output -----")
        print(text)
        println("----- end $name output -----")
    end
end

if !isempty(failures)
    error("Example smoke failures: $(join(failures, ", "))")
end

println("examples_suite_smoke_ok")
