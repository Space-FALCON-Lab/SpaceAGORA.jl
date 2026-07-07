@testset "Coverage Debt Probe Drivers" begin
    # Each probe file is standalone (bootstraps its own SimulationModel include),
    # so it runs as a subprocess like the Coverage Threaded Probe Driver above.
    # Unlike that driver these run in every suite pass — they are regular tests
    # for the July 2026 coverage-debt register files; the subprocess inherits
    # the coverage flag so the coverage gates see their line data.
    probe_files = [
        "coverage_debt_reference_system_probes.jl",
        "coverage_debt_parallel_policy_probes.jl",
        "coverage_debt_density_selection_probes.jl",
        "coverage_debt_rpo_planning_probes.jl",
        "coverage_debt_rpo_nav_probes.jl",
        "coverage_debt_thruster_guidance_probes.jl",
        "coverage_debt_hypr_search_probes.jl",
        "coverage_debt_robot_arm_hypr_probes.jl",
        "coverage_debt_allocator_assets_probes.jl",
    ]
    coverage_flags = Base.JLOptions().code_coverage == 0 ? String[] : ["--code-coverage=user"]
    for probe in probe_files
        probe_script = joinpath(REPO_ROOT, "test", probe)
        cmd = Cmd([
            Base.julia_cmd().exec...,
            "--startup-file=no",
            "--depwarn=error",
            "--project=$(REPO_ROOT)",
            coverage_flags...,
            "--threads=2",
            probe_script,
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
            println("----- begin $(probe) output -----")
            println(text)
            println("----- end $(probe) output -----")
        end
        @test success(proc)
    end
end
