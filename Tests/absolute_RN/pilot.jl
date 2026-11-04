length(ARGS) == 9 || error("pilot.jl new|baseline LABEL ALT INC HELPERS POLICY ORBITS POWER diagnostic|benchmark|control")
kind, label, altitude, inclination, helpers, policy, orbits, power, mode = ARGS
kind in ("new", "baseline") || error("Unknown pilot implementation")
mode in ("diagnostic", "benchmark", "control") || error("Unknown pilot mode")
occursin(r"^[A-Za-z0-9_-]+$", label) || error("Pilot labels must not contain paths")
kind == "baseline" && mode != "benchmark" && error("Use benchmark mode with the frozen old library")
const RUNNER_ROOT = kind == "new" ? normpath(joinpath(@__DIR__, "..", "..", "Simulation_absolute_RN")) :
    joinpath(@__DIR__, "results", "source_snapshot_pre_absolute_RN", "project", "Simulation")
include(joinpath(RUNNER_ROOT, "Code", "StudyRunner.jl"))
using .OracleStudy
using .OracleStudy.StudyConfig
using Arrow, DataFrames, LinearAlgebra, StaticArrays, TOML
include("diagnostic_probe.jl")

function main()
    c = StudyCase(parse(Float64, altitude), parse(Float64, inclination), parse(Int, helpers),
        Symbol(policy); orbits=parse(Float64, orbits), power_W=parse(Float64, power))
    root = joinpath(@__DIR__, "PilotOutput", label)
    report_path = joinpath(@__DIR__, "results", label * ".toml")
    (ispath(root) || ispath(report_path)) && error("Refusing to overwrite pilot $label")
    BLAS.set_num_threads(1)
    probe = mode == "diagnostic" ? diagnostic_probe(c) : nothing
    timing = timing_probe(c)
    observer = function (integrator, start, pairs)
        timing.observe(integrator, start, pairs)
        probe === nothing || probe.observe(integrator, start, pairs)
        nothing
    end
    started = time()
    directory = run_case(c; output_root=root, record_events=mode != "control",
        observer, save_schedule=true)
    wall = time() - started
    report = probe === nothing ? Dict{String,Any}() : probe.report()
    merge!(report, Dict("kind" => kind, "label" => label, "mode" => mode,
        "output_directory" => directory, "pilot_wall_s" => wall,
        "timing" => timing.report(), "settings" => case_settings(c)))
    if probe !== nothing
        reference_path = joinpath(root, "independent_reference.feather")
        Arrow.write(reference_path, probe.reference; file=true, compress=:zstd)
        report["reference_file"] = reference_path
    end
    info = TOML.parsefile(joinpath(directory, "run_info.toml"))
    for key in ("wall_elapsed_s", "peak_rss_bytes", "output_bytes", "temporary_dense_state_bytes",
        "solver_accepted_steps", "solver_rejected_steps", "schedule_sha256", "source_signature")
        report[key] = info[key]
    end
    atomic_toml(report_path, report)
    TOML.print(stdout, report; sorted=true)
    probe === nothing || report["passed"] || error("Independent diagnostic validation failed: $report_path")
end

main()
