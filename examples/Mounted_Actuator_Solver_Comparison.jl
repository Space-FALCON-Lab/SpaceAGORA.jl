# Reproduce docs/src/tutorials/solver_comparison.md. No external assets.
include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "support", "mounted_actuator_case.jl"))
using CSV, DataFrames, TOML, SHA

# Same case and physical clocks in every run, including the refined reference.
const DURATION_S = 120.0
const SAMPLE_S = 0.25
const SAMPLE_TIMES = collect(0.0:SAMPLE_S:DURATION_S)
const OUTPUT_ROOT = joinpath(pwd(), "output", "solver_comparison")

function comparison_tolerances(relative, absolute, cap)
    IntegrationTolerances(reltol_orbit=relative, abstol_orbit=absolute,
        reltol_quaternion=relative, abstol_quaternion=absolute,
        reltol_angular_rate=relative, abstol_angular_rate=absolute,
        dt_max_orbit=cap)
end

function run_case(; mode, cap, substeps=8, relative=1e-8, absolute=1e-10)
    config = mounted_actuator_configuration(duration_s=DURATION_S, sample_s=SAMPLE_S,
        solver=SolverConfig(solver_mode=mode, parallel=false,
            multirate_slow_dt_s=cap, multirate_fast_substeps=substeps),
        tolerances=comparison_tolerances(relative, absolute, cap))
    # The metadata route cannot be combined with return_results=true. Sample the
    # returned solution on a common grid instead of comparing adaptive step indices.
    elapsed = @elapsed result = run_simulation(config;
        return_solution=true, return_solver_metadata=true)
    result.retcode == "Success" || error("$mode did not finish successfully: $(result.retcode)")
    result.solution.t[end] == DURATION_S || error("$mode stopped before the common horizon")
    samples = [result.solution(t).sc[1] for t in SAMPLE_TIMES]
    all(s -> all(isfinite, s), samples) || error("$mode returned non-finite samples")
    return (; samples, result, elapsed)
end

# q and -q encode the same attitude. atan avoids acos losing precision near zero.
function attitude_error(a, b)
    qa, qb = normalize(SVector{4,Float64}(a)), normalize(SVector{4,Float64}(b))
    dot(qa, qb) < 0 && (qb = -qb)
    4 * atan(norm(qa - qb), norm(qa + qb))
end

function errors(candidate, reference)
    rows = map(eachindex(SAMPLE_TIMES)) do i
        a, b = candidate.samples[i], reference.samples[i]
        (; time_s=SAMPLE_TIMES[i], position_m=norm(a.pos-b.pos),
            velocity_m_s=norm(a.vel-b.vel), attitude_rad=attitude_error(a.q,b.q),
            body_rate_rad_s=norm(a.ω-b.ω), wheel_momentum_n_m_s=norm(a.h_wheels-b.h_wheels))
    end
    DataFrame(rows)
end

function comparison()
    mkpath(OUTPUT_ROOT)
    output = mktempdir(OUTPUT_ROOT; prefix="attempt-", cleanup=false)
    println("This attempt writes to ", output)
    # These controls remove timing/calibration variability; the explicit
    # SolverConfig pins the numerical mode independently of shell overrides.
    controls = ("SPACEAGORA_RHS_CALIBRATE"=>"off", "SPACEAGORA_RHS_IDENTIFY"=>"0",
        "SPACEAGORA_RHS_EXECUTION_MODE"=>"serial", "SPACEAGORA_VISUALIZATION"=>"0",
        "SPACEAGORA_SOLVER_SAVE_EVERYSTEP"=>"true", "SPACEAGORA_SOLVER_SAVE_ON"=>"true",
        "SPACEAGORA_SOLVER_SAVE_START"=>"true", "SPACEAGORA_SOLVER_SAVE_END"=>"true")
    withenv(controls...) do
        metadata = Dict("julia_version"=>string(VERSION), "threads"=>Threads.nthreads(),
            "duration_s"=>DURATION_S, "sample_s"=>SAMPLE_S, "control_period_s"=>1.0,
            "revision"=>strip(read(`git -C $REPO_ROOT rev-parse HEAD`, String)),
            "working_tree_changes"=>read(`git -C $REPO_ROOT status --porcelain`, String),
            "root_manifest_sha256"=>bytes2hex(sha256(read(joinpath(REPO_ROOT,"Manifest.toml")))),
            "examples_manifest_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"Manifest.toml")))),
            "candidate_relative_tolerance"=>1e-8, "candidate_absolute_tolerance"=>1e-10,
            "reference_relative_tolerance"=>1e-13, "reference_absolute_tolerance"=>1e-14,
            "refined_relative_tolerance"=>1e-14, "refined_absolute_tolerance"=>1e-15,
            "controlled_environment"=>Dict(controls),
            "solver_environment"=>Dict(k=>v for (k,v) in ENV if startswith(k,"SPACEAGORA_SOLVER_") || k in ("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS", "SPACEAGORA_MULTIRATE_SLOW_DT_S")),
            "status"=>"running")
        open(joinpath(output,"settings.toml"),"w") do io
            TOML.print(io, metadata)
        end
        try
            reference = run_case(mode=:dp8, cap=0.03125, relative=1e-13, absolute=1e-14)
            refined = run_case(mode=:dp8, cap=0.015625, relative=1e-14, absolute=1e-15)
            reference_errors = errors(reference, refined)
            CSV.write(joinpath(output, "reference_refinement.csv"), reference_errors)
            metrics = names(reference_errors)[2:end]
            reference_max = Dict(key => maximum(reference_errors[!, key]) for key in metrics)
            cases = [
                (; name="default_tsit5", mode=:tsit5, cap=1.0, substeps=8),
                [(; name="substeps_$(n)", mode=:multirate_mri, cap=1.0, substeps=n) for n in (1,2,4,8)]...,
                [(; name="outer_$(cap)", mode=:multirate_mri, cap, substeps=8) for cap in (0.5,0.25)]...,
            ]
            summary = NamedTuple[]
            traces = Dict{String,Any}("reference"=>repr(reference.result.solver_trace),
                "refined_reference"=>repr(refined.result.solver_trace))
            for case in cases
                println("Comparing ", case.name)
                run = run_case(; case.mode, case.cap, case.substeps)
                table = errors(run, refined)
                CSV.write(joinpath(output, case.name * ".csv"), table)
                traces[case.name] = repr(run.result.solver_trace)
                for metric in metrics
                    peak = maximum(table[!, metric])
                    # A reporting aid, not an accepted mission tolerance. If the
                    # reference changes by >10% of this difference, do not rank it.
                    resolved = reference_max[metric] <= 0.1 * peak && peak > 0
                    push!(summary, (; case=case.name, mode=string(case.mode), outer_cap_s=case.cap,
                        fast_substeps=case.substeps, metric, max_error=peak,
                        rms_error=sqrt(sum(abs2, table[!,metric])/nrow(table)),
                        reference_change=reference_max[metric], reference_resolved=resolved,
                        retcode=run.result.retcode, elapsed_s=run.elapsed))
                end
            end
            CSV.write(joinpath(output, "summary.csv"), DataFrame(summary))
            metadata["solver_traces"] = traces
            metadata["status"] = "complete"
            open(joinpath(output,"settings.toml"),"w") do io
                TOML.print(io, metadata)
            end
            show(stdout, MIME"text/plain"(), DataFrame(summary)); println()
            println("Comparison written to ", output)
        catch err
            metadata["status"] = "failed"
            metadata["error"] = sprint(showerror, err)
            open(joinpath(output,"settings.toml"),"w") do io
                TOML.print(io, metadata)
            end
            rethrow()
        end
    end
end
comparison()
