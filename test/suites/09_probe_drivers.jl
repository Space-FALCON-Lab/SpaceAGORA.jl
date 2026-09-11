@testset "Probe Drivers" begin
    # Each probe file is standalone (it does `using SpaceAGORA` and binds its own
    # module aliases, so its top-level helpers never collide with this scope), so
    # it runs as a subprocess like the Coverage Threaded Probe Driver above.
    # Unlike that driver these run in every suite pass — they are regular tests
    # for the probe files under test/probes/ (the paid-off July 2026 coverage-debt register); the subprocess inherits
    # the coverage flag so the coverage gates see their line data.
    probe_files = [
        "reference_system_probes.jl",
        "parallel_policy_probes.jl",
        "density_selection_probes.jl",
        "ei_partition_drag_probes.jl",
        "fm_incidence_probes.jl",
        "eddy_damping_probes.jl",
        "aero_torque_probes.jl",
        "rpo_planning_probes.jl",
        "rpo_nav_probes.jl",
        "thruster_guidance_probes.jl",
        "hypr_search_probes.jl",
        "robot_arm_hypr_probes.jl",
        "allocator_assets_probes.jl",
        "process_pool_probes.jl",
        "campaign_process_route_probes.jl",
        "state_anchor_probes.jl",
        "flat_route_parity_probes.jl",
        "kinematics_probes.jl",
    ]
    # This suite is 46% of the default entrypoint's wall clock, and every probe is
    # an independent subprocess, so it is the one worth sharding across CI jobs.
    # `SPACEAGORA_PROBE_SHARD="i/n"` runs probes i, i+n, i+2n, ... of the list.
    # Unset runs all of them. Round-robin rather than contiguous blocks, which is
    # a cheap heuristic and not a good one: measured at 2 shards it splits 9
    # probes each into 4m29s and 8m09s, because the probes differ in cost by far
    # more than index parity captures. The workflow compensates by giving the
    # lighter shard the unit-tree driver (about 2m05s) rather than by balancing
    # here, so nobody has to maintain a per-probe cost table. If the shard count
    # changes, re-measure rather than assuming this still evens out.
    probe_files = let raw = strip(get(ENV, "SPACEAGORA_PROBE_SHARD", ""))
        if isempty(raw)
            probe_files
        else
            parts = split(raw, "/")
            length(parts) == 2 || error("SPACEAGORA_PROBE_SHARD must look like \"1/3\", got \"$(raw)\"")
            idx = tryparse(Int, parts[1]); count = tryparse(Int, parts[2])
            (idx === nothing || count === nothing || count < 1 || idx < 1 || idx > count) &&
                error("SPACEAGORA_PROBE_SHARD=\"$(raw)\" is not a valid i/n with 1 <= i <= n")
            selected = [probe_files[i] for i in idx:count:length(probe_files)]
            println("probe shard $(idx)/$(count): $(length(selected)) of $(length(probe_files)) probes")
            selected
        end
    end
    coverage_flags = Base.JLOptions().code_coverage == 0 ? String[] : ["--code-coverage=user"]
    for probe in probe_files
        probe_script = joinpath(REPO_ROOT, "test", "probes", probe)
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

# The unit tree is one subprocess and cannot be split, so when suite 09 is
# sharded it must run in exactly one shard. `SPACEAGORA_SKIP_UNIT_DRIVER=1` lets
# the other shards leave it alone; the shard gate checks exactly one omits it.
const _RUN_UNIT_DRIVER = get(ENV, "SPACEAGORA_SKIP_UNIT_DRIVER", "0") != "1"

@testset "Standalone Unit Suite Driver" begin
    if !_RUN_UNIT_DRIVER
        @test true   # claimed by another shard
    else
    # test/unit/**/*_tests.jl are standalone (`using SpaceAGORA`) function-level
    # tests for the parallel routing, cost and calibration layers: machine
    # topology, the cost hierarchy, robust timing, the streaming paired trial,
    # the R6 route rules, RHS re-verification. No CI job ran test/unit/runtests.jl,
    # so the code they exercise reached the coverage gate with no line data at
    # all -- parallel/cost/machine_calibration.jl measured 0%. Run the whole
    # unit tree here as one subprocess: one Julia start-up, the coverage flag
    # forwarded, and four threads where the machine has them because several
    # testsets (rounds ties, the split race) skip below that.
    unit_script = joinpath(REPO_ROOT, "test", "unit", "runtests.jl")
    coverage_flags = Base.JLOptions().code_coverage == 0 ? String[] : ["--code-coverage=user"]
    unit_threads = clamp(Sys.CPU_THREADS, 2, 4)
    cmd = Cmd([
        Base.julia_cmd().exec...,
        "--startup-file=no",
        "--project=$(REPO_ROOT)",
        coverage_flags...,
        "--threads=$(unit_threads)",
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
        println("----- begin unit suite output -----")
        println(text)
        println("----- end unit suite output -----")
    end
    @test success(proc)
    end
end
