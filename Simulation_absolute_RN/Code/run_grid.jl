include("StudyConfig.jl")
using .StudyConfig
using Dates, TOML

function batch_source_hashes(hashes)
    # Standalone plotting scripts are not loaded by simulation or verification.
    filter(pair -> !(startswith(first(pair), "Code/plot_") &&
        endswith(first(pair), ".jl")), hashes)
end

function same_simulation_sources(stored, current)
    verification_and_launch = Set(("Code/VerifyOutput.jl", "Code/verify_output.jl",
        "Code/run_grid.jl", "Code/launch_study.sh"))
    filter(pair -> first(pair) ∉ verification_and_launch, batch_source_hashes(stored)) ==
        filter(pair -> first(pair) ∉ verification_and_launch, batch_source_hashes(current))
end

function run_child(command, log_path, status_path, status)
    open(log_path, "a") do io
        process = run(pipeline(command; stdout=io, stderr=io); wait=false)
        status["child_pid"] = getpid(process)
        atomic_toml(status_path, status)
        wait(process)
        success(process) || error("Child process failed; see $log_path")
    end
end

function append_progress(output_root, event, scenario; elapsed_s=nothing, detail="", campaign="absolute_rn")
    stamp = Dates.format(now(UTC), "yyyy-mm-ddTHH:MM:SS") * "Z"
    elapsed = elapsed_s === nothing ? "" : " | elapsed=$(round(elapsed_s; digits=2)) s"
    message = replace(detail, '\n' => ' ', '\r' => ' ')
    open(joinpath(output_root, campaign_settings(campaign).progress_file), "a") do io
        println(io, "$stamp | $event | $scenario$elapsed",
            isempty(message) ? "" : " | $message")
        flush(io)
    end
end

function tracked_child(command, log, status_path, status, output_root, scenario, phase; campaign="absolute_rn")
    started = time_ns()
    append_progress(output_root, "$(phase)_STARTED", scenario; campaign, detail="details=$log")
    try
        run_child(command, log, status_path, status)
    catch
        append_progress(output_root, "$(phase)_FAILED", scenario;
            campaign, elapsed_s=(time_ns() - started) / 1e9, detail="See $log")
        rethrow()
    end
    elapsed = (time_ns() - started) / 1e9
    append_progress(output_root, "$(phase)_FINISHED", scenario; campaign, elapsed_s=elapsed)
    elapsed
end

function case_disk_requirement(case)
    rows = ceil(Int, mission_duration(case) / ENCOUNTER_SAMPLE_S) + 1
    dense_spool = rows * (6 * (case.helpers + 1) + 1) * sizeof(Float64)
    raw_trajectory = TRAJECTORY_ROWS * length(trajectory_columns(case.helpers + 1)) * sizeof(Float64)
    5 * 1024^3 + dense_spool + raw_trajectory + 64 * 1024^2
end

function simulation_command(output_root, case)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --threads=1 --gcthreads=1 $(joinpath(@__DIR__, "run_case.jl")) $output_root $(case.altitude_km) $(case.inclination_deg) $(case.helpers) $(case.policy) $(case.orbits) $(case.power_W) $(case.helper_altitude_km)`
end

function run_grid(output_root; available_space=free_disk_bytes, campaign="absolute_rn")
    output_root = abspath(output_root)
    selection = campaign_settings(campaign)
    batch = joinpath(output_root, selection.batch_directory)
    mkpath(joinpath(batch, "logs"))
    lock = joinpath(batch, "lock")
    if isdir(lock)
        owner_path = joinpath(lock, "owner.toml")
        isfile(owner_path) || error("Lock has no owner information: $lock")
        owner = TOML.parsefile(owner_path)
        ccall(:kill, Cint, (Cint, Cint), owner["pid"], 0) == 0 &&
            error("A batch process is already running with PID $(owner["pid"])")
        mv(lock, lock * ".stale." * Dates.format(now(UTC), "yyyymmddTHHMMSS"))
    end
    mkdir(lock)
    atomic_toml(joinpath(lock, "owner.toml"), Dict("pid" => getpid(), "started_utc" => string(now(UTC))))
    status_path = joinpath(batch, "status.toml")
    cases = grid_cases(campaign)
    hashes = source_hashes()
    signature = source_signature(hashes)
    status = Dict{String,Any}("status" => "running", "pid" => getpid(),
        "campaign" => campaign, "helper_counts" => collect(selection.helper_counts),
        "helper_altitudes_km" => collect(selection.helper_altitudes_km),
        "target_altitudes_km" => collect(selection.target_altitudes_km),
        "total_cases" => length(cases), "verified_cases" => 0, "reused_cases" => 0,
        "source_signature" => signature, "started_utc" => string(now(UTC)))
    atomic_toml(status_path, status)
    append_progress(output_root, "BATCH_STARTED", "$(length(cases)) scenarios";
        campaign, detail="campaign=$campaign; pid=$(getpid()); all timestamps UTC; elapsed time is wall-clock seconds")
    try
        for (index, case) in enumerate(cases)
            if isfile(joinpath(batch, "STOP_AFTER_CASE"))
                status["status"] = "stopped_after_case"
                append_progress(output_root, "BATCH_PAUSED", "stop-after-case requested"; campaign)
                break
            end
            batch_source_hashes(source_hashes()) == batch_source_hashes(hashes) ||
                error("Production code changed during the batch; stopping to avoid mixed versions")
            directory = case_directory(output_root, case)
            available, required = available_space(output_root), case_disk_requirement(case)
            if !isdir(directory) && available < required
                status["status"] = "paused_low_disk"
                status["next_case"] = directory
                status["available_disk_bytes"] = available
                status["required_free_bytes"] = required
                @warn "Batch paused before the next case; free space and resume to continue" available required
                append_progress(output_root, "BATCH_PAUSED_LOW_DISK", relpath(directory, output_root);
                    campaign, detail="available_bytes=$available required_bytes=$required")
                break
            end
            id = case_signature(case)
            status["case_index"] = index
            status["active_case"] = directory
            status["updated_utc"] = string(now(UTC))
            atomic_toml(status_path, status)
            log = joinpath(batch, "logs", id * ".log")
            info_path = joinpath(directory, "run_info.toml")
            scenario = "[$index/$(length(cases))] " * relpath(directory, output_root)
            scenario_started = time_ns()
            if isdir(directory)
                isfile(info_path) || error("Incomplete output directory: $directory")
                info = TOML.parsefile(info_path)
                info["status"] in ("complete", "verified") ||
                    error("Case needs investigation before retry: $directory")
                info["settings"] == case_settings(case) && info["run_id"] == id &&
                    same_simulation_sources(info["source_hashes"], hashes) &&
                    info["record_events"] || error("Existing case settings/code do not match: $directory")
                status["reused_cases"] += 1
                append_progress(output_root, "SIMULATION_REUSED", scenario;
                    campaign, elapsed_s=info["wall_elapsed_s"],
                    detail="previous simulation elapsed time; verifying saved files again")
            else
                command = simulation_command(output_root, case)
                tracked_child(command, log, status_path, status, output_root, scenario, "SIMULATION"; campaign)
            end
            julia = Base.julia_cmd()
            verify = `$julia --startup-file=no --threads=1 --gcthreads=1 $(joinpath(@__DIR__, "verify_output.jl")) $directory`
            tracked_child(verify, log, status_path, status, output_root, scenario, "VERIFICATION"; campaign)
            status["verified_cases"] += 1
            println("VERIFIED $(status["verified_cases"])/$(length(cases)): $directory")
            flush(stdout)
            append_progress(output_root, "SCENARIO_VERIFIED", scenario;
                campaign, elapsed_s=(time_ns() - scenario_started) / 1e9,
                detail="simulation plus verification, or re-verification of a reused case")
            atomic_toml(status_path, status)
        end
        status["verified_cases"] == length(cases) && (status["status"] = "complete")
        status["finished_utc"] = string(now(UTC))
        atomic_toml(status_path, status)
        append_progress(output_root, "BATCH_ENDED", status["status"];
            campaign, detail="verified=$(status["verified_cases"])/$(length(cases))")
    catch err
        status["status"] = "failed"
        status["error"] = sprint(showerror, err)
        status["failed_utc"] = string(now(UTC))
        atomic_toml(status_path, status)
        append_progress(output_root, "BATCH_FAILED", get(status, "active_case", "startup");
            campaign, detail=status["error"])
        @error "Study batch stopped; completed outputs were retained" exception=(err, catch_backtrace())
        rethrow()
    finally
        rm(joinpath(lock, "owner.toml"))
        rm(lock)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 2 || error("Usage: julia run_grid.jl [OUTPUT_ROOT] [absolute_rn]")
    run_grid(isempty(ARGS) ? DEFAULT_OUTPUT : ARGS[1];
        campaign=length(ARGS) == 2 ? ARGS[2] : "absolute_rn")
end
