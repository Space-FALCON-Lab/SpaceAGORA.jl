include("run_grid.jl")
using Dates, TOML

function run_jobs()
    paths = CampaignPaths
    config = ProductionBatch.StudyConfig
    job = paths.JOB_DIRECTORY
    mkpath(job)
    status_path = joinpath(job, "status.toml")
    status = Dict{String,Any}("status" => "running", "phase" => "startup",
        "pid" => getpid(), "started_utc" => string(now(UTC)), "total_cases" => 243,
        "laser_cases" => 162, "no_laser_cases" => 81,
        "laser_output" => paths.OUTPUT, "no_laser_output" => paths.NO_LASER_OUTPUT,
        "difference_output" => paths.DIFFERENCE_OUTPUT)
    sources = paths.campaign_sources(config.source_hashes)
    status["source_signature"] = config.source_signature(sources)
    function phase!(name)
        paths.campaign_sources(config.source_hashes) == sources ||
            error("Source code changed during the comparison job")
        status["phase"] = name
        status["updated_utc"] = string(now(UTC))
        config.atomic_toml(status_path, status)
        println("PHASE $name at $(now(UTC)) UTC")
        flush(stdout)
    end
    function pause_requested()
        isfile(joinpath(job, "STOP_AFTER_CASE"))
    end
    function plot_child(script, name)
        phase!(name)
        julia = Base.julia_cmd()
        ProductionBatch.run_child(
            `$julia --startup-file=no --threads=1 --gcthreads=1 $script`,
            joinpath(job, "plotting.log"), status_path, status)
        delete!(status, "child_pid")
    end
    try
        for (campaign, root, plot_script) in (
                ("range_along", paths.OUTPUT, joinpath(paths.CODE_ROOT, "plot_campaign.jl")),
                ("no_laser", paths.NO_LASER_OUTPUT, joinpath(paths.NO_LASER_CODE, "plot_results.jl")))
            selection = config.campaign_settings(campaign)
            phase!("$(campaign)_simulations")
            status["active_batch_status"] = joinpath(root, selection.batch_directory, "status.toml")
            config.atomic_toml(status_path, status)
            run_campaign(campaign, root)
            batch_status = TOML.parsefile(status["active_batch_status"])
            status["$(campaign)_verified_cases"] = batch_status["verified_cases"]
            if batch_status["status"] != "complete" || pause_requested()
                status["status"] = "paused"
                status["reason"] = pause_requested() ? "stop_after_case" : batch_status["status"]
                config.atomic_toml(status_path, status)
                return status
            end
            plot_child(plot_script, "$(campaign)_plots_and_encounter_evaluation")
        end
        plot_child(joinpath(paths.DIFFERENCE_CODE, "plot_results.jl"), "laser_minus_no_laser_plots")
        phase!("complete")
        status["status"] = "complete"
        status["finished_utc"] = string(now(UTC))
        config.atomic_toml(status_path, status)
        status
    catch err
        status["status"] = "failed"
        status["error"] = sprint(showerror, err)
        status["failed_utc"] = string(now(UTC))
        config.atomic_toml(status_path, status)
        @error "Comparison job failed; completed results were retained" exception=(err, catch_backtrace())
        rethrow()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) || error("Usage: run_jobs.jl")
    run_jobs()
end
