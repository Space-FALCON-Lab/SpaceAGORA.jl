include(joinpath(@__DIR__, "..", "..", "Simulation_absolute_RN", "Code", "StudyConfig.jl"))
using .StudyConfig
using Dates, TOML

const RESULTS = joinpath(@__DIR__, "results")
const BATCH = joinpath(DEFAULT_OUTPUT, "_batch_absolute_RN")
const SESSION = "oracle-study-absolute-rn"
alive(pid) = ccall(:kill, Cint, (Cint, Cint), pid, 0) == 0

function main()
    preflight = TOML.parsefile(joinpath(RESULTS, "preflight.toml"))
    preflight["passed"] || error("Preflight did not pass")
    first_directory = case_directory(DEFAULT_OUTPUT, first(grid_cases()))
    deadline = time() + 600
    while time() < deadline
        status_path = joinpath(BATCH, "status.toml")
        if isfile(status_path)
            status = TOML.parsefile(status_path)
            status["status"] == "running" || error("Campaign stopped: $(status["status"])")
            status["total_cases"] == 270 || error("Wrong production case count")
            status["source_signature"] == preflight["source_signature"] ||
                error("Production source differs from validated source")
            alive(status["pid"]) || error("Batch process is no longer alive")
            if status["verified_cases"] >= 1 && haskey(status, "child_pid") && alive(status["child_pid"])
                info = TOML.parsefile(joinpath(first_directory, "run_info.toml"))
                verification = TOML.parsefile(joinpath(first_directory, "verification.toml"))
                progress = TOML.parsefile(joinpath(first_directory, "progress.toml"))
                info["status"] == "verified" && verification["passed"] && verification["rows"] == 1001 ||
                    error("First production case was not verified")
                info["settings"]["schema_version"] == 3 &&
                    info["settings"]["initial_target_periods"] == 100 ||
                    error("First case has incorrect schema/duration")
                info["julia_threads"] == info["blas_threads"] == 1 ||
                    error("Production simulation thread limits were not applied")
                progress["fraction"] == 1.0 || error("First case progress did not reach its endpoint")
                text = read(joinpath(DEFAULT_OUTPUT, "simulation_progress_absolute_RN.log"), String)
                all(occursin(marker, text) for marker in
                    ("SIMULATION_STARTED", "SIMULATION_FINISHED", "VERIFICATION_FINISHED", "SCENARIO_VERIFIED")) ||
                    error("Campaign progress log is missing completed phase transitions")
                success(`tmux has-session -t $SESSION`) || error("Persistent tmux session is absent")
                report = Dict("passed" => true, "checked_utc" => string(now(UTC)),
                    "tmux_session" => SESSION, "campaign" => "absolute_rn",
                    "total_cases" => 270, "verified_cases_at_check" => status["verified_cases"],
                    "parent_pid" => status["pid"], "active_child_pid" => status["child_pid"],
                    "active_case" => status["active_case"], "first_case" => first_directory,
                    "first_case_rows" => verification["rows"], "schema_version" => 3,
                    "first_case_simulation_wall_s" => info["wall_elapsed_s"],
                    "source_signature" => preflight["source_signature"],
                    "output_root" => DEFAULT_OUTPUT,
                    "progress_log" => joinpath(DEFAULT_OUTPUT, "simulation_progress_absolute_RN.log"))
                atomic_toml(joinpath(RESULTS, "launch_confirmation.toml"), report)
                TOML.print(stdout, report; sorted=true)
                return
            end
        end
        sleep(5)
    end
    error("First production verification was not observed within ten minutes; campaign was not interrupted")
end

main()
