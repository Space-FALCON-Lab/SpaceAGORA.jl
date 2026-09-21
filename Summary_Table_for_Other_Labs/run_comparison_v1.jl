include("run_comparison_v2.jl")

##################
# Worker imports #
##################
if length(ARGS) == 4 && ARGS[1] == "--worker"
    if ARGS[2] == "prototype"
        include("1_Kuang's Prototype Code/test16_feather.jl")
        if get(TOML.parsefile(ARGS[3]), "video", false)
            import GLMakie
            using GeometryBasics, FileIO
            include("1_Kuang's Prototype Code/functions/10_Animation_ver2.jl")
        end
    elseif ARGS[2] == "spaceagora"
        include("2_SpaceAGORA.jl/ORACLE/run_case2_laser_links.jl")
    else
        error("Unknown worker: $(ARGS[2])")
    end
end

########
# Main #
########
function main(arguments = ARGS)
    ### 1. Execute a worker request ###
    if length(arguments) == 4 && arguments[1] == "--worker"
        config = TOML.parsefile(arguments[3])
        settings = (; (Symbol(key)=>value for (key, value) in config["case"])...)
        result = arguments[2] == "prototype" ?
            run_prototype(settings, config["smoke"]; video = get(config, "video", false)) :
            run_spaceagora(settings, config["smoke"])
        open(arguments[4], "w") do io
            TOML.print(io, result)
        end
        return nothing
    end

    ### 2. Validate settings and resolve the common duration ###
    all(argument -> argument in ("--smoke", "--video"), arguments) ||
        error("Usage: julia run_comparison.jl [--smoke] [--video]")
    smoke = "--smoke" in arguments
    validate_case(TEST_CASE)
    settings = resolve_duration(TEST_CASE; smoke)

    ### 3. Run each model in its own Julia process ###
    runs = mktempdir() do directory
        config_path = joinpath(directory, "case.toml")
        open(config_path, "w") do io
            TOML.print(io, Dict(
                "case"  => Dict(string(key) => value for (key, value) in pairs(settings)),
                "smoke" => smoke,
                "video" => "--video" in arguments,
            ))
        end
        results = Dict{String,Any}[]
        for model in ("prototype", "spaceagora")
            result_path = joinpath(directory, model*".toml")
            println("Running $model for $(settings.duration_s) seconds...")
            run(`$(Base.julia_cmd()) --startup-file=no --project=$SPACEAGORA_ROOT $(@__FILE__) --worker $model $config_path $result_path`)
            push!(results, TOML.parsefile(result_path))
        end
        results
    end

    ### 4. Write the comparison report ###
    report_path = joinpath(WORKSPACE_ROOT, smoke ? "comparison_summary_smoke.md" : "comparison_summary.md")
    write_report(report_path, settings, runs; smoke)
    println("Report saved to $report_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end