@testset "Parallel Unit Tests" begin
    # The router branch's unit tests under test/unit/parallel/ are standalone
    # files (each does `using SpaceAGORA` and binds its own module aliases), so
    # they run as subprocesses exactly like the probe drivers in suite 09 and
    # inherit the coverage flag, which is what makes the cost-model, contention
    # and dispatch modules visible to the coverage gates. Three files are left
    # to the standalone `test/unit/runtests.jl` driver because they assert
    # machine-specific facts (`machine_topology_tests.jl` assumes physical
    # cores <= Sys.CPU_THREADS, false on Apple Silicon where Julia counts only
    # performance cores; `mc_route_tie_tests.jl` and `policy_v2_tests.jl`
    # assume the host's full thread and process budget).
    unit_files = [
        "contention_inputs_tests.jl",
        "cost_machine_calibration_tests.jl",
        "cost_robust_timing_tests.jl",
        "cost_work_counts_tests.jl",
        "lock_width_cap_tests.jl",
        "mixed_dispatch_tests.jl",
        "native_lock_stats_tests.jl",
        "outer_route_persistence_tests.jl",
        "outer_split_budget_tests.jl",
        "rhs_batch_budget_tests.jl",
        "streaming_trial_tests.jl",
        "usl_tests.jl",
    ]
    coverage_flags = Base.JLOptions().code_coverage == 0 ? String[] : ["--code-coverage=user"]
    for unit in unit_files
        unit_script = joinpath(REPO_ROOT, "test", "unit", "parallel", unit)
        cmd = Cmd([
            Base.julia_cmd().exec...,
            "--startup-file=no",
            "--depwarn=error",
            "--project=$(REPO_ROOT)",
            coverage_flags...,
            # Four threads, not the probe drivers' two: the batch-budget test
            # asserts the inner budget up to three workers, which needs the
            # threads to exist. The CI runners have four.
            "--threads=4",
            unit_script,
        ])
        cmd = addenv(
            cmd,
            "SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0",
            "SPACEAGORA_WARN_NORMALIZE" => "0"
        )

        output = IOBuffer()
        proc = run(pipeline(ignorestatus(cmd), stdout=output, stderr=output))
        text = String(take!(output))
        if !success(proc)
            println("----- begin $(unit) output -----")
            println(text)
            println("----- end $(unit) output -----")
        end
        @test success(proc)
    end
end
