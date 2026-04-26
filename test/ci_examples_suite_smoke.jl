const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const EXAMPLES_DIR = joinpath(REPO_ROOT, "examples")
const PROJECT_PATH = joinpath(REPO_ROOT, ".AGORA")

function list_examples()
    helper_files = Set(["common.jl", "aerobraking_mission_plot_utils.jl"])
    files = sort(filter(f -> endswith(f, ".jl"), readdir(EXAMPLES_DIR; join=true)))
    files = filter(f -> !(basename(f) in helper_files), files)
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
        cmd = `$(Base.julia_cmd()) --startup-file=no --compiled-modules=existing --depwarn=error --project=$(PROJECT_PATH) $(example_path)`
        cmd = Cmd(cmd; dir=tmp)
        cmd = addenv(
            cmd,
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
