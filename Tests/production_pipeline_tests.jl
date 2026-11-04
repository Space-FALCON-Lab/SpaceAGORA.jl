include(joinpath(@__DIR__, "..", "Simulation", "Code", "StudyRunner.jl"))
using .OracleStudy
using .OracleStudy.StudyConfig
using .OracleStudy.EncounterTimingRecorder
using Arrow, DataFrames, LinearAlgebra, StaticArrays, Test, TOML

const SE = OracleStudy.SimulationEngine
const SM = OracleStudy.SM
const TEST_ROOT = joinpath(@__DIR__, "results", "production_pipeline_compact")

module BatchCompatibility
include(joinpath(@__DIR__, "..", "Simulation", "Code", "run_grid.jl"))
end

function diagnostic_probe(c)
    count = c.helpers + 1
    totals = zeros(5 + count)
    max_rtn = 0.0
    max_signed = 0.0
    max_work = 0.0
    max_local_quadrature = 0.0
    steps = 0
    times = trajectory_times(c)
    next_sample = 2
    reference = DataFrame(time=[0.0], target_dv_R_mps=[0.0], target_dv_T_mps=[0.0],
        target_dv_N_mps=[0.0], target_dv_T_positive_mps=[0.0], target_dv_T_negative_mps=[0.0])
    for sc in 1:count
        reference[!, Symbol("sc$(sc)_laser_work_J")] = [0.0]
    end

    function rates(u, pairs)
        result = zeros(5 + count)
        isempty(pairs) && return result
        helper, target = only(pairs)
        r = SVector{3,Float64}(u.sc[target].pos)
        v = SVector{3,Float64}(u.sc[target].vel)
        displacement = r - SVector{3,Float64}(u.sc[helper].pos)
        force = (c.power_W * B / LIGHT_SPEED) * displacement / norm(displacement)
        radial = r / norm(r)
        angular_momentum = cross(r, v)
        normal = angular_momentum / norm(angular_momentum)
        transverse = cross(normal, radial)
        result[1:3] .= (dot(force, radial), dot(force, transverse), dot(force, normal)) ./ MASS_KG
        result[4] = max(result[2], 0.0)
        result[5] = min(result[2], 0.0)
        result[5 + target] = dot(force, v)
        result[5 + helper] = -dot(force, u.sc[helper].vel)
        result
    end

    function observe(integrator, start, pairs)
        stop = Float64(integrator.t)
        dt = stop - start
        function integrate_to(end_time)
            isempty(pairs) && return zeros(length(totals))
            width = end_time - start
            values = [rates(integrator(start + width * fraction), pairs) for fraction in 0.0:0.125:1.0]
            fine = width / 24 * (values[1] + 4values[2] + 2values[3] +
                4values[4] + 2values[5] + 4values[6] + 2values[7] + 4values[8] + values[9])
            coarse = width / 12 * (values[1] + 4values[3] + 2values[5] + 4values[7] + values[9])
            max_local_quadrature = max(max_local_quadrature, norm(fine[1:3] - coarse[1:3]))
            fine
        end
        while next_sample <= length(times) && times[next_sample] <= stop
            sample_time = times[next_sample]
            values = totals + integrate_to(sample_time)
            push!(reference, (sample_time, values...))
            next_sample += 1
        end
        totals .+= integrate_to(stop)
        target = integrator.u.sc[count]
        max_rtn = max(max_rtn, norm(target.laser_rtn - totals[1:3]))
        max_signed = max(max_signed, abs(target.laser_dv_T_positive - totals[4]),
            abs(target.laser_dv_T_negative - totals[5]))
        for sc in 1:count
            max_work = max(max_work, abs(integrator.u.sc[sc].laser_work - totals[5 + sc]))
        end
        steps += 1
        nothing
    end

    function report()
        Dict{String,Any}("steps" => steps, "independent_RTN_mps" => totals[1:3],
            "independent_T_positive_mps" => totals[4], "independent_T_negative_mps" => totals[5],
            "independent_work_J" => totals[6:end], "max_RTN_error_mps" => max_rtn,
            "max_signed_T_error_mps" => max_signed, "max_work_error_J" => max_work,
            "max_local_RTN_quadrature_difference_mps" => max_local_quadrature,
            "RTN_tolerance_mps" => 1e-8, "signed_T_tolerance_mps" => 1e-7,
            "work_tolerance_J" => 0.1 + 1e-8maximum(abs, totals[6:end]))
    end
    (; observe, report, reference)
end

function fixtures()
    @testset "Study grid and configuration" begin
        cases = grid_cases()
        @test length(cases) == 210
        @test length(unique(case_signature.(cases))) == 210
        @test length(unique(case_directory.(Ref(DEFAULT_OUTPUT), cases))) == 210
        @test all(c.orbits == 100.0 for c in cases)
        @test_throws ArgumentError StudyCase(1050, 5, 50, :gve_sma)
        @test_throws ArgumentError StudyCase(1050, 1, 20, :gve_sma)
        @test_throws ArgumentError StudyCase(1050, 1, 50, :naive)
        @test_throws ArgumentError StudyCase(1050, 1, 50, :gve_sma; orbits=NaN)
        @test_throws ArgumentError StudyCase(1050, 1, 50, :gve_sma; power_W=-1)
        for c in cases[1:2:end]
            args = study_configuration(c)
            initial = SE.build_initial_conditions(args)
            @test OracleStudy.check_initial_geometry(initial) > 1e-3
            @test isapprox(args.dynamics_model.spacecraft[end].initial_condition.ν, deg2rad(1.0))
            @test args.scheduling_policy_model.target_idx == c.helpers + 1
        end
    end
    @testset "Solver-integrated diagnostics and force" begin
        for power in (0.0, 10_000.0)
            c = StudyCase(1050, 1, 1, :gve_sma; power_W=power)
            args = study_configuration(c)
            state = SE.build_initial_conditions(args)
            state.sc[1].mass = 2MASS_KG
            args.interlink_model.graph.active .= true
            derivative = zero(state)
            SE._apply_interlink_rhs!(derivative, state, (args=args, is_active=trues(2)))
            for i in 1:2
                sc, du = state.sc[i], derivative.sc[i]
                force = SM.laser_force_on_spacecraft(args.interlink_model, state, i)
                radial = normalize(sc.pos)
                normal = normalize(cross(sc.pos, sc.vel))
                transverse = cross(normal, radial)
                @test isapprox(norm(force), B * power / LIGHT_SPEED; atol=1e-15)
                @test du.vel == du.laser_dv == force / sc.mass
                @test isapprox(du.laser_rtn, [dot(du.vel, radial), dot(du.vel, transverse),
                    dot(du.vel, normal)]; atol=1e-18)
                @test du.laser_work == dot(force, sc.vel)
                @test du.laser_dv_T_positive >= 0 && du.laser_dv_T_negative <= 0
                @test du.laser_dv_T_positive + du.laser_dv_T_negative == du.laser_rtn[2]
                @test iszero(sc.laser_rtn) && sc.laser_work == 0
            end
            @test norm(state.sc[1].mass * derivative.sc[1].vel +
                state.sc[2].mass * derivative.sc[2].vel) < 1e-15
            @test Set((:laser_rtn, :laser_work, :laser_dv_T_positive, :laser_dv_T_negative)) ⊆
                Set(field.name for field in SM.default_save_fields(args))
            args.interlink_model.graph.active .= false
            SE._apply_interlink_rhs!(derivative, state, (args=args, is_active=trues(2)))
            @test all(iszero(du.laser_rtn) && du.laser_dv_T_positive == 0 &&
                du.laser_dv_T_negative == 0 && du.laser_work == 0 for du in derivative.sc)
            legacy = SM.SimConfig._with_configuration(args; interlink_diagnostics=false)
            @test !hasproperty(SE.build_initial_conditions(legacy).sc[1], :laser_rtn)
            invalid = SM.SimConfig._with_configuration(args; interlink_model=nothing)
            @test_throws ArgumentError SE._validate_interlinks!(invalid, :tsit5)
            state.sc[2].pos .= state.sc[1].pos
            @test_throws ArgumentError OracleStudy.check_initial_geometry(state)
        end
    end
    @testset "Multiple crossings and range audit" begin
        state(t) = [0.0 16t^2-16t+3; 0.0 0.0; 0.0 0.0;
            0.0 32t-16; 0.0 0.0; 0.0 0.0]
        recorder = Recorder(state(0), 0.5)
        observe_firing!(recorder, [(1, 2)], state(0), 0.0; initial=true)
        observe_geometry!(recorder, state(0), state(1), 0.0, 1.0, state)
        finish!(recorder, state(1), 1.0)
        @test nrow(recorder.geometry) == 2
        segments = range_segments(recorder.geometry, recorder.laser; run_id="multiple-crossings")
        @test nrow(segments) == 5
        @test sum(segments.duration_s) == 1.0
        @test segments.within_range == [false, true, false, true, false]
        @test segments.end_time_s[1:end-1] == segments.start_time_s[2:end]
        outside = [0.0 2.0; zeros(4, 2); 0.0 0.0]
        recorder = Recorder(outside, 1.0)
        observe_firing!(recorder, [(1, 2)], outside, 0.0; initial=true)
        finish!(recorder, outside, 1.0)
        segments = range_segments(recorder.geometry, recorder.laser; run_id="outside")
        @test nrow(segments) == 1 && !only(segments.within_range)
        empty = Recorder(outside, 1.0)
        finish!(empty, outside, 1.0)
        @test isempty(range_segments(empty.geometry, empty.laser; run_id="empty"))
        tiny(t) = [0.0 1e-4*(-2+4t); 0.0 0.0; 0.0 0.0;
            0.0 4e-4; 0.0 0.0; 0.0 0.0]
        small = Recorder(tiny(0), 1e-4)
        observe_geometry!(small, tiny(0), tiny(1), 0.0, 1.0, tiny)
        finish!(small, tiny(1), 1.0)
        @test nrow(small.geometry) == 1
        @test only(small.geometry.duration_s) ≈ 0.5
        final_entry(t) = [0.0 2-t; 0.0 0.0; 0.0 0.0;
            0.0 -1.0; 0.0 0.0; 0.0 0.0]
        ending = Recorder(final_entry(0), 1.0)
        observe_geometry!(ending, final_entry(0), final_entry(1), 0.0, 1.0, final_entry)
        finish!(ending, final_entry(1), 1.0)
        @test isempty(ending.geometry)
    end
    @testset "Compact trajectory and dense encounter products" begin
        case = StudyCase(1050, 1, 300, :gve_sma)
        @test length(trajectory_times(case)) == 1001
        @test first(trajectory_times(case)) == 0
        @test last(trajectory_times(case)) == mission_duration(case)
        state(t) = [0.0 -2+4t; 0.0 0.2; 0.0 0.0; 0.0 4.0; 0.0 0.0; 0.0 0.0]
        recorder = Recorder(state(0), 1.0)
        observe_geometry!(recorder, state(0), state(1), 0.0, 1.0, state)
        finish!(recorder, state(1), 1.0)
        mktemp() do path, io
            for t in 0.0:0.1:1.0
                OracleStudy.write_dense_sample(io, t, state(t))
            end
            flush(io)
            products = OracleStudy.build_encounter_products!(recorder.geometry, path, 2, "fixture", 1.0; cadence=0.1)
            @test nrow(products) == 7
            @test only(recorder.geometry.sample_count) == 7
            @test only(recorder.geometry.min_relative_speed_mps) == 4.0
            @test only(recorder.geometry.max_relative_speed_mps) == 4.0
            @test only(recorder.geometry.min_separation_m) ≈ 0.2
            @test only(recorder.geometry.max_separation_m) ≈ 1.0
            @test only(recorder.geometry.min_range_rate_mps) < 0
            @test only(recorder.geometry.max_range_rate_mps) > 0
            @test all(products.is_all_pair_longest) && all(products.is_target_helper_longest)
            @test first(products.time) == only(recorder.geometry.start_time_s)
            @test last(products.time) == only(recorder.geometry.end_time_s)
        end
        empty = Recorder(state(0), 1.0)
        products = OracleStudy.build_encounter_products!(empty.geometry, "unused", 2, "empty", 1.0)
        @test isempty(products)
        @test eltype(products.time) == Float64
        @test eltype(empty.geometry.min_range_rate_mps) == Float64
    end
    mkpath(TEST_ROOT)
    atomic_toml(joinpath(TEST_ROOT, "fixtures.toml"), Dict("passed" => true))
end

function pilot(args)
    7 <= length(args) <= 8 || error(
        "pilot ROOT ALT INC HELPERS POLICY ORBITS MODE(recorded|control|zero) [HELPER_KM=1000]")
    root, h, i, n, policy, orbits, mode = args[1:7]
    mode in ("recorded", "control", "zero") || error("Unknown pilot mode")
    helper_altitude = length(args) == 8 ? parse(Float64, args[8]) : 1000.0
    c = StudyCase(parse(Float64, h), parse(Float64, i), parse(Int, n), Symbol(policy);
        orbits=parse(Float64, orbits), power_W=mode == "zero" ? 0.0 : 10_000.0,
        helper_altitude_km=helper_altitude)
    probe = diagnostic_probe(c)
    directory = run_case(c; output_root=root, record_events=mode != "control",
        observer=probe.observe, save_schedule=true)
    report = probe.report()
    report["passed"] = report["steps"] > 0 &&
        report["max_RTN_error_mps"] <= report["RTN_tolerance_mps"] &&
        report["max_signed_T_error_mps"] <= report["signed_T_tolerance_mps"] &&
        report["max_work_error_J"] <= report["work_tolerance_J"]
    report["output_directory"] = directory
    evidence = joinpath(TEST_ROOT, "diagnostics")
    mkpath(evidence)
    label = case_signature(c) * "_" * mode * "_" * basename(normpath(root))
    reference_path = joinpath(evidence, label * ".feather")
    Arrow.write(reference_path, probe.reference; file=true, compress=:zstd)
    report["reference_file"] = reference_path
    atomic_toml(joinpath(evidence, label * ".toml"), report)
    info_path = joinpath(directory, "run_info.toml")
    info = TOML.parsefile(info_path)
    info["independent_diagnostics_passed"] = report["passed"]
    report["passed"] || (info["status"] = "validation_failed")
    atomic_toml(info_path, info)
    TOML.print(stdout, report; sorted=true)
    @testset "Independent production RTN/work validation" begin
        @test report["steps"] > 0
        @test report["max_RTN_error_mps"] <= report["RTN_tolerance_mps"]
        @test report["max_signed_T_error_mps"] <= report["signed_T_tolerance_mps"]
        @test report["max_work_error_J"] <= report["work_tolerance_J"]
    end

end

function verify_reference(report_path)
    report = TOML.parsefile(report_path)
    reference = Arrow.Table(report["reference_file"])
    actual = Arrow.Table(joinpath(report["output_directory"], "trajectory.feather"))
    @testset "Fresh-file RTN/work against independent dense quadrature" begin
        @test report["passed"]
        @test length(actual.time) == TRAJECTORY_ROWS
        @test actual.time == reference.time
        for name in (:target_dv_R_mps, :target_dv_T_mps, :target_dv_N_mps)
            @test maximum(abs.(actual[name] - reference[name])) <= report["RTN_tolerance_mps"]
        end
        for name in (:target_dv_T_positive_mps, :target_dv_T_negative_mps)
            @test maximum(abs.(actual[name] - reference[name])) <= report["signed_T_tolerance_mps"]
        end
        for name in propertynames(reference)
            endswith(string(name), "_laser_work_J") || continue
            @test maximum(abs.(actual[name] - reference[name])) <= report["work_tolerance_J"]
        end
    end
    report["fresh_file_reference_passed"] = true
    atomic_toml(report_path, report)
end

function compare(recorded, control; check_schedule_file=true)
    @testset "Production recorder noninterference" begin
        a, b = TOML.parsefile(joinpath(recorded, "run_info.toml")), TOML.parsefile(joinpath(control, "run_info.toml"))
        @test a["run_id"] == b["run_id"]
        @test BatchCompatibility.same_simulation_sources(a["source_hashes"], b["source_hashes"])
        @test a["schedule_sha256"] == b["schedule_sha256"]
        @test a["solver_accepted_steps"] == b["solver_accepted_steps"]
        @test a["solver_rejected_steps"] == b["solver_rejected_steps"]
        left = Arrow.Table(joinpath(recorded, "trajectory.feather"))
        right = Arrow.Table(joinpath(control, "trajectory.feather"))
        @test propertynames(left) == propertynames(right)
        @test all(left[name] == right[name] for name in propertynames(left))
        if check_schedule_file
            @test DataFrame(Arrow.Table(joinpath(recorded, "schedule.feather"))) ==
                DataFrame(Arrow.Table(joinpath(control, "schedule.feather")))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("Use fixtures, pilot, or compare")
    if ARGS[1] == "fixtures"
        fixtures()
    elseif ARGS[1] == "pilot"
        pilot(ARGS[2:end])
    elseif ARGS[1] in ("compare", "compare-plain")
        length(ARGS) == 3 || error("compare RECORDED_DIRECTORY CONTROL_DIRECTORY")
        compare(ARGS[2], ARGS[3]; check_schedule_file=ARGS[1] == "compare")
    elseif ARGS[1] == "verify-reference"
        length(ARGS) == 2 || error("verify-reference REPORT_TOML")
        verify_reference(ARGS[2])
    else
        error("Unknown test command")
    end
end
