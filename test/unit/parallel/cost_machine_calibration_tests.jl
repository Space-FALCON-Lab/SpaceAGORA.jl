using Test
using SpaceAGORA

const PC = SpaceAGORA.SimulationModel.ParallelCost

@testset "Machine fingerprint" begin
    a = PC.machine_fingerprint()
    @test a == PC.machine_fingerprint()          # stable within a session
    @test length(a) == 16                        # 8 bytes hex

    # An explicit machine label must produce a distinct identity, so two
    # machines that happen to share a CPU model do not silently share constants.
    b = withenv("SPACEAGORA_PERF_MACHINE_LABEL" => "some-other-box") do
        PC.machine_fingerprint()
    end
    @test a != b
end

function _synthetic_constants()
    return PC.MachineConstants(
        usl_alpha_base = 0.05,
        usl_beta_alloc = 0.01,
        fingerprint = "testfingerprint",
        schema_version = PC.CALIBRATION_SCHEMA_VERSION,
        measured_threads = 8,
    )
end

const _KEPT_CONSTANT_KEYS = Set(["schema_version", "fingerprint", "measured_threads",
                                 "usl_alpha_base", "usl_beta_alloc"])

@testset "Constants round-trip through TOML, and only the kept keys are written" begin
    mktempdir() do dir
        path = joinpath(dir, "constants.toml")
        mc = _synthetic_constants()
        PC.save_machine_constants(mc, path)
        @test isfile(path)
        @test Set(keys(PC.TOML.parsefile(path))) == _KEPT_CONSTANT_KEYS

        back = PC.load_machine_constants(path)
        @test back !== nothing
        for field in fieldnames(PC.MachineConstants)
            @test getfield(back, field) == getfield(mc, field)
        end

        # Built with the defaults, a value saves a file the loader reads: the
        # schema defaults to the current one.
        @test PC.MachineConstants().schema_version == PC.CALIBRATION_SCHEMA_VERSION
        built = joinpath(dir, "defaults.toml")
        PC.save_machine_constants(PC.MachineConstants(usl_alpha_base = 0.1, usl_beta_alloc = 0.2), built)
        @test PC.load_machine_constants(built) == PC.MachineConstants(usl_alpha_base = 0.1, usl_beta_alloc = 0.2)

        # A campaign table written by the trace extraction survives a re-save,
        # and is the only other table the file holds.
        write(path, read(path, String) * "\n[campaign]\npool_cold_start_s = 1.5\n")
        PC.save_machine_constants(mc, path)
        @test Set(keys(PC.TOML.parsefile(path))) == union(_KEPT_CONSTANT_KEYS, Set(["campaign"]))
        @test PC.TOML.parsefile(path)["campaign"]["pool_cold_start_s"] == 1.5
    end
end

@testset "Stale or damaged constants are declined, not reinterpreted" begin
    mktempdir() do dir
        @test PC.load_machine_constants(joinpath(dir, "absent.toml")) === nothing

        garbage = joinpath(dir, "garbage.toml")
        write(garbage, "this is not toml {{{")
        @test PC.load_machine_constants(garbage) === nothing

        # A file from a different schema version must be discarded rather than
        # read: term semantics changed by definition when the version moved, so
        # its numbers do not mean what the current readers would take them to.
        stale = joinpath(dir, "stale.toml")
        PC.save_machine_constants(_synthetic_constants(), stale)
        body = replace(read(stale, String),
                       "schema_version = $(PC.CALIBRATION_SCHEMA_VERSION)" => "schema_version = 999")
        write(stale, body)
        @test PC.load_machine_constants(stale) === nothing

        # A current-schema file without one of the two terms the planner reads,
        # or without the measuring thread count, is not current either, so it is
        # measured again rather than read with a made-up value.
        for key in ("usl_alpha_base", "usl_beta_alloc", "measured_threads")
            partial = joinpath(dir, "without_$(key).toml")
            PC.save_machine_constants(_synthetic_constants(), partial)
            write(partial, join(filter(l -> !startswith(l, key), readlines(partial)), "\n"))
            @test PC.load_machine_constants(partial) === nothing
        end
    end
end

@testset "calibrate_machine measures the USL terms the planner reads" begin
    mc = PC.calibrate_machine(k = 1)
    # Non-negative but not bounded above: with two threads the fit has one rung
    # above serial, and a noisy reading can put alpha past 1 (the speedup it
    # predicts is then floored at 1).
    @test isfinite(mc.usl_alpha_base) && mc.usl_alpha_base >= 0.0
    @test isfinite(mc.usl_beta_alloc) && mc.usl_beta_alloc >= 0.0
    @test mc.measured_threads == Threads.nthreads()
    @test mc.schema_version == PC.CALIBRATION_SCHEMA_VERSION
    @test mc.fingerprint == PC.machine_fingerprint()
end

@testset "the calibration keeps only what something reads" begin
    @test fieldnames(PC.MachineConstants) ==
        (:usl_alpha_base, :usl_beta_alloc, :fingerprint, :schema_version, :measured_threads)
    # Removed with the constants only the RHS-plan predictor read: the rate
    # curves, the staleness canary and its reference kernels, and the kernels of
    # the ladders that measured those constants.
    for name in (:RateCurve, :rate_at, :constants_are_current, :reference_kernel_ns,
                 :reference_memory_kernel_ns, :calib_touch_kernel!, :calib_scalar_kernel!,
                 :calib_stream_kernel!, :_calibrate_coeff_touch, :_calibrate_simd_lane,
                 :_calibrate_scalar_item, :_calibrate_dispatch, :_calibrate_polyester_dispatch,
                 :_calibrate_parallel_speedup, :_llc_bytes, :_least_squares_line)
        @test !isdefined(PC, name)
    end
    for name in (:calib_batch_kernel!, :calib_table, :calib_alloc_kernel!, :_calibrate_usl,
                 :_speedup_ladder, :timing_noise_floor)
        @test isdefined(PC, name)
    end
end

@testset "the unwired RHS-plan predictor is gone; its live neighbours stay" begin
    SE = SpaceAGORA.SimulationEngine
    # Removed in 0.6.0: the analytic predictor, its work counts and its
    # per-effector probe. Nothing in the package called them.
    for name in (:WorkCounts, :effector_cost_terms, :constellation_work_counts,
                 :model_in_cost_domain, :flat_queue_node_effector, :EffectorProbe,
                 :validate_declaration, :probe_ns_map, :PlanCandidate, :PlanPrediction,
                 :predict_plan_ns, :plan_candidates, :select_plan, :LOCK_CONVOY_ABSTAIN_RATIO,
                 :workload_usl_parameters, :_cost_batchable_effector,
                 :_cost_harmonics_prepass_effector, :_cost_aero_prepass_effector,
                 :_alf_recurrence_iterations)
        @test !isdefined(PC, name)
    end
    @test !isdefined(SE, :probe_effector_costs!)
    @test !isdefined(SE, :warn_on_declaration_mismatch)
    @test !hasfield(PC.ContentionInputs, :workspace_bytes_per_sat)

    # Kept: what the native-lock width cap, the predictive planner and the RHS
    # calibration read.
    for name in (:ContentionInputs, :alloc_bytes_per_second, :lock_duty_cycle,
                 :lock_wait_hold_ratio, :MachineConstants, :calibrate_machine,
                 :ensure_machine_constants!, :load_machine_constants,
                 :save_machine_constants, :usl_speedup, :usl_peak_workers, :_fit_usl,
                 :timed_min, :StreamingPairedTrial)
        @test isdefined(PC, name)
    end
    @test isdefined(SE, :probe_contention_inputs!)
end
