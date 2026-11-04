include(joinpath(@__DIR__, "..", "..", "Simulation_absolute_RN", "Code", "StudyRunner.jl"))
using .OracleStudy
using .OracleStudy.StudyConfig
using Arrow, ComponentArrays, DataFrames, LinearAlgebra, StaticArrays, Test, TOML

include(joinpath(@__DIR__, "..", "..", "Simulation_absolute_RN", "Code", "VerifyOutput.jl"))
module BatchControl
include(joinpath(@__DIR__, "..", "..", "Simulation_absolute_RN", "Code", "run_grid.jl"))
end

const SE = OracleStudy.SimulationEngine
const SM = OracleStudy.SM
const RESULTS = joinpath(@__DIR__, "results")
const COMPONENTS = ("R", "T", "N")
const BRANCH_FIELDS = Tuple(Symbol("laser_dv_$(component)_$sign")
    for component in COMPONENTS for sign in ("positive", "negative"))
include("diagnostic_probe.jl")
mkpath(joinpath(RESULTS, "tmp"))
ENV["TMPDIR"] = joinpath(RESULTS, "tmp")
atomic_toml(joinpath(RESULTS, "fixtures.toml"), Dict("passed" => false, "status" => "running"))

@testset "Exact 270-case grid and preserved configuration" begin
    cases = grid_cases()
    expected = Set((Float64(h), Float64(i), n, p)
        for h in (900, 950, 1000, 1050, 1100), i in (0, 0.5, 1),
        n in (1, 100, 200, 300, 400, 500, 600, 700, 800), p in (:gve_sma, :gve_inc))
    @test length(cases) == length(unique(case_signature.(cases))) == 270
    @test Set((c.altitude_km, c.inclination_deg, c.helpers, c.policy) for c in cases) == expected
    @test all(c.orbits == 100 && c.power_W == 10_000 && c.helper_altitude_km == 1000 for c in cases)
    @test length(unique(case_directory.(Ref(DEFAULT_OUTPUT), cases))) == 270
    @test basename(dirname(DEFAULT_OUTPUT)) == "Simulation_absolute_RN"
    @test campaign_settings("absolute_rn").progress_file == "simulation_progress_absolute_RN.log"
    @test_throws ArgumentError campaign_settings("helpers700")
    @test_throws ArgumentError campaign_settings("main")
    for h in (850, 1150)
        @test_throws ArgumentError StudyCase(h, 1, 100, :gve_sma)
    end
    for n in (50, 150, 250, 701)
        @test_throws ArgumentError StudyCase(1000, 1, n, :gve_sma)
    end
    @test_throws ArgumentError StudyCase(1000, 1, 100, :gve_sma; helper_altitude_km=900)
    @test_throws ArgumentError StudyCase(1000, 1, 100, :gve_sma; orbits=NaN)
    @test_throws ArgumentError StudyCase(1000, 1, 100, :gve_sma; power_W=-1)
    for c in cases
        settings = case_settings(c)
        times = trajectory_times(c)
        @test settings["schema_version"] == 3
        @test settings["solver"] == "Tsit5"
        @test settings["dt_max_s"] == settings["encounter_sample_interval_s"] == 10
        @test settings["reltol"] == settings["abstol"] == 1e-12
        @test length(times) == 1001
        @test first(times) == 0 && last(times) == mission_duration(c)
        @test times[801] ≈ 80 * 2pi * sqrt((EARTH_RADIUS_M + 1000c.altitude_km)^3 / MU) atol=1e-9 rtol=0
        @test length(trajectory_columns(c.helpers + 1)) == 10(c.helpers + 1) + 10
    end
    for n in (1, 100, 800), policy in (:gve_sma, :gve_inc)
        c = StudyCase(950, 1, n, policy)
        args = study_configuration(c)
        state = SE.build_initial_conditions(args)
        @test args.scheduling_policy_model.target_idx == n + 1
        @test args.interlink_diagnostics
        @test length(state.sc) == n + 1
        @test OracleStudy.check_initial_geometry(state) > 1e-3
        for k in 1:n
            expected_r = (EARTH_RADIUS_M + 1_000_000) *
                SVector(cosd(360(k - 1) / n), sind(360(k - 1) / n), 0)
            @test state.sc[k].pos ≈ expected_r atol=1e-8 rtol=0
            @test iszero(state.sc[k].laser_rtn) && state.sc[k].laser_work == 0
            @test all(iszero(getproperty(state.sc[k], field)) for field in BRANCH_FIELDS)
        end
        target_r = (EARTH_RADIUS_M + 950_000) *
            SVector(cosd(1), sind(1) * cosd(1), sind(1)^2)
        @test state.sc[n + 1].pos ≈ target_r atol=1e-8 rtol=0
    end
end

@testset "Independent absolute quadrature across sign changes" begin
    function rates(t)
        net = [t - 0.37, 0.5, 0.4 - t^2]
        vcat(net, max.(net, 0), min.(net, 0), abs.(net))
    end
    result = split_integral(rates, 0.0, 1.0)
    net = [0.13, 0.5, 0.4 - 1 / 3]
    positive = [(1 - 0.37)^2 / 2, 0.5, 2 / 3 * 0.4^1.5]
    negative = net - positive
    @test result.fine[1:3] ≈ net atol=1e-14 rtol=0
    @test result.fine[4:6] ≈ positive atol=1e-14 rtol=0
    @test result.fine[7:9] ≈ negative atol=1e-14 rtol=0
    @test result.fine[10:12] ≈ positive - negative atol=1e-14 rtol=0
    @test maximum(result.difference) <= 1e-14
end

@testset "Laser-only RHS, all branch signs, zero resets and disabled callers" begin
    for power in (0.0, 10_000.0), sr in (-1, 1), st in (-1, 1), sn in (-1, 1)
        c = StudyCase(1000, 1, 1, :gve_sma; power_W=power)
        args = study_configuration(c)
        state = SE.build_initial_conditions(args)
        radial = SVector(1.0, 0.0, 0.0)
        transverse = SVector(0.0, cosd(1), sind(1))
        normal = cross(radial, transverse)
        state.sc[2].pos .= (EARTH_RADIUS_M + 1_000_000) * radial
        state.sc[2].vel .= 7350 * transverse
        state.sc[1].pos .= state.sc[2].pos - 1000(sr * radial + st * transverse + sn * normal)
        state.sc[1].mass = 2MASS_KG
        args.interlink_model.graph.active .= true
        derivative = zero(state)
        gravity = SVector(-7.3, 0.2, -0.1)
        for du in derivative.sc
            du.vel .= gravity
        end
        SE._apply_interlink_rhs!(derivative, state, (args=args, is_active=trues(2)))
        for index in 1:2
            sc, du = state.sc[index], derivative.sc[index]
            displacement = SVector{3}(sc.pos - state.sc[3 - index].pos)
            force = (B * power / LIGHT_SPEED) * normalize(displacement)
            rhat = normalize(SVector{3}(sc.pos))
            nhat = normalize(cross(SVector{3}(sc.pos), SVector{3}(sc.vel)))
            that = cross(nhat, rhat)
            projected = SVector(dot(force / sc.mass, rhat), dot(force / sc.mass, that),
                dot(force / sc.mass, nhat))
            @test du.vel ≈ gravity + force / sc.mass atol=1e-15 rtol=0
            @test du.laser_dv ≈ force / sc.mass atol=1e-18 rtol=0
            @test du.laser_rtn ≈ projected atol=1e-18 rtol=0
            @test du.laser_work ≈ dot(force, sc.vel) atol=1e-12 rtol=0
            for (axis, component) in enumerate(COMPONENTS)
                positive = getproperty(du, Symbol("laser_dv_$(component)_positive"))
                negative = getproperty(du, Symbol("laser_dv_$(component)_negative"))
                @test positive >= 0 && negative <= 0
                @test positive + negative == du.laser_rtn[axis]
                @test positive - negative == abs(du.laser_rtn[axis])
                @test positive ≈ max(projected[axis], 0) atol=1e-18 rtol=0
                @test negative ≈ min(projected[axis], 0) atol=1e-18 rtol=0
            end
        end
        legacy_args = SM.SimConfig._with_configuration(args; interlink_diagnostics=false)
        legacy = SE.build_initial_conditions(legacy_args)
        for index in 1:2
            for field in (:pos, :vel)
                getproperty(legacy.sc[index], field) .= getproperty(state.sc[index], field)
            end
            legacy.sc[index].mass = state.sc[index].mass
        end
        legacy_derivative = zero(legacy)
        for du in legacy_derivative.sc
            du.vel .= gravity
        end
        SE._apply_interlink_rhs!(legacy_derivative, legacy, (args=legacy_args, is_active=trues(2)))
        @test !hasproperty(legacy.sc[1], :laser_rtn)
        @test all(!hasproperty(legacy.sc[1], field) for field in BRANCH_FIELDS)
        for index in 1:2, field in (:vel, :laser_dv, :laser_delta_sma)
            @test getproperty(legacy_derivative.sc[index], field) == getproperty(derivative.sc[index], field)
        end
        @test SE._validate_interlink_checkpoint_state(args, state) === nothing
        @test SE._validate_interlink_checkpoint_state(legacy_args, legacy) === nothing
        @test_throws ArgumentError SE._validate_interlink_checkpoint_state(args, legacy)
        @test_throws ArgumentError SE._validate_interlink_checkpoint_state(legacy_args, state)
        args.interlink_model.graph.active .= false
        SE._apply_interlink_rhs!(derivative, state, (args=args, is_active=trues(2)))
        @test all(iszero(du.laser_rtn) && du.laser_work == 0 &&
            all(iszero(getproperty(du, field)) for field in BRANCH_FIELDS) for du in derivative.sc)
    end
    args = study_configuration(StudyCase(1000, 1, 1, :gve_sma))
    state = SE.build_initial_conditions(args)
    old_fields = filter(name -> name ∉ (:laser_dv_R_positive, :laser_dv_R_negative,
        :laser_dv_N_positive, :laser_dv_N_negative), propertynames(state.sc[1]))
    old_state = ComponentVector(sc=[NamedTuple{Tuple(old_fields)}(
        Tuple(getproperty(sc, name) for name in old_fields)) for sc in state.sc])
    @test_throws ArgumentError SE._validate_interlink_checkpoint_state(args, old_state)
    positional = SM.SimulationConfiguration(
        (getfield(args, name) for name in fieldnames(typeof(args))[1:end-1])...)
    @test !positional.interlink_diagnostics
    no_interlinks = SM.SimulationConfiguration(
        (getfield(args, name) for name in fieldnames(typeof(args))[1:12])...)
    @test no_interlinks.interlink_model === nothing && !no_interlinks.interlink_diagnostics
    invalid = SM.SimConfig._with_configuration(args; interlink_model=nothing)
    @test_throws ArgumentError SE._validate_interlinks!(invalid, :tsit5)
end

@testset "Target indexing, fields, units, branch guards and schema rejection" begin
    c = StudyCase(1000, 1, 100, :gve_sma)
    args = study_configuration(c)
    state = SE.build_initial_conditions(args)
    for (index, sc) in enumerate(state.sc), (axis, component) in enumerate(COMPONENTS)
        sc.laser_rtn[axis] = axis * index
        setproperty!(sc, Symbol("laser_dv_$(component)_positive"), 2.0axis * index)
        setproperty!(sc, Symbol("laser_dv_$(component)_negative"), -1.0axis * index)
    end
    fields = SM.default_save_fields(args)
    for field in BRANCH_FIELDS
        selected = only(filter(saved -> saved.name == field, fields))
        @test selected.per_satellite
        @test selected.getter(state, 0.0, nothing) ==
            [getproperty(sc, field) for sc in state.sc]
    end
    disabled = SM.default_save_fields(SM.SimConfig._with_configuration(args; interlink_diagnostics=false))
    @test all(field.name ∉ BRANCH_FIELDS && field.name != :laser_rtn for field in disabled)
    mktempdir(RESULTS) do root
        path = joinpath(root, "target_index.feather")
        writer = OracleStudy.TrajectoryWriter(path, 101, "target-index"; chunk_rows=1)
        OracleStudy.save_sample!(writer, 0.0, state)
        OracleStudy.save_sample!(writer, 1.0, state)
        OracleStudy.close_trajectory!(writer)
        table = Arrow.Table(path)
        @test collect(propertynames(table)) == trajectory_columns(101)
        @test table.time == [0, 1]
        for (axis, component) in enumerate(COMPONENTS)
            @test table[Symbol("target_dv_$(component)_mps")] == fill(101axis, 2)
            @test table[Symbol("target_dv_$(component)_positive_mps")] == fill(202axis, 2)
            @test table[Symbol("target_dv_$(component)_negative_mps")] == fill(-101axis, 2)
        end
        metadata = Arrow.getmetadata(table)
        @test metadata["target_id"] == "101" && metadata["schema_version"] == "3"
        @test occursin("RTN:m/s", metadata["units"])
        @test metadata["rtn_branch_definition"] == RTN_BRANCH_DEFINITION
        atomic_toml(joinpath(root, "run_info.toml"), Dict("status" => "verified",
            "record_events" => true, "settings" => Dict("schema_version" => 2)))
        @test_throws ErrorException StudyVerification.verify_output(root)
        directory = case_directory(root, c)
        mkpath(directory)
        @test_throws ErrorException run_case(c; output_root=root)
    end
    for component in COMPONENTS
        @test StudyVerification.check_signed_totals(1.0, -0.5, 0.5, 0.9, -0.4, component) == 0
        @test_throws ErrorException StudyVerification.check_signed_totals(1.0, -0.5, 0.6, 0.9, -0.4, component)
        @test_throws ErrorException StudyVerification.check_signed_totals(-0.1, -0.5, -0.6, 0.0, -0.4, component)
        @test_throws ErrorException StudyVerification.check_signed_totals(1.0, 0.1, 1.1, 0.9, 0.0, component)
        @test_throws ErrorException StudyVerification.check_signed_totals(0.8, -0.5, 0.3, 0.9, -0.4, component)
    end
end

@testset "Isolated batch stop, disk, locks, source guards and safe reuse" begin
    mktempdir(RESULTS) do root
        batch = joinpath(root, "_batch_absolute_RN")
        mkpath(batch)
        stop = joinpath(batch, "STOP_AFTER_CASE")
        touch(stop)
        BatchControl.run_grid(root)
        status = TOML.parsefile(joinpath(batch, "status.toml"))
        @test status["status"] == "stopped_after_case"
        @test status["total_cases"] == 270 && status["verified_cases"] == 0
        @test !ispath(joinpath(root, "paper_plot_mode")) && !ispath(joinpath(batch, "lock"))
        rm(stop)
        BatchControl.run_grid(root; available_space=_ -> 0)
        status = TOML.parsefile(joinpath(batch, "status.toml"))
        @test status["status"] == "paused_low_disk" && !haskey(status, "child_pid")
        @test status["required_free_bytes"] > 5 * 1024^3
        @test isfile(joinpath(root, "simulation_progress_absolute_RN.log"))
        @test !isfile(joinpath(root, "simulation_progress.log"))
        lock = joinpath(batch, "lock")
        mkdir(lock)
        atomic_toml(joinpath(lock, "owner.toml"), Dict("pid" => getpid()))
        @test_throws ErrorException BatchControl.run_grid(root)
        rm(joinpath(lock, "owner.toml"))
        rm(lock)
        hashes = BatchControl.source_hashes()
        changed = copy(hashes)
        changed["Code/plot_fixture.jl"] = "postprocessing-only"
        @test BatchControl.same_simulation_sources(hashes, changed)
        @test BatchControl.batch_source_hashes(hashes) == BatchControl.batch_source_hashes(changed)
        for name in ("Code/StudyConfig.jl", "Code/StudyRunner.jl", "Code/TrajectoryWriter.jl",
            "../9_SpaceAGORA.jl-main_full integration/src/orchestration/simulation/engine/dynamics_rhs.jl")
            changed = copy(hashes)
            changed[name] = "changed-runtime"
            @test !BatchControl.same_simulation_sources(hashes, changed)
        end
        changed = copy(hashes)
        changed["Code/VerifyOutput.jl"] = "new-verifier"
        @test BatchControl.same_simulation_sources(hashes, changed)
        @test BatchControl.batch_source_hashes(hashes) != BatchControl.batch_source_hashes(changed)
        case = first(BatchControl.grid_cases())
        directory = BatchControl.case_directory(root, case)
        mkpath(directory)
        settings = BatchControl.case_settings(case)
        settings["schema_version"] = 2
        atomic_toml(joinpath(directory, "run_info.toml"), Dict("status" => "verified",
            "run_id" => BatchControl.case_signature(case), "settings" => settings,
            "source_hashes" => hashes, "record_events" => true))
        @test_throws ErrorException BatchControl.run_grid(root)
        status = TOML.parsefile(joinpath(batch, "status.toml"))
        @test status["status"] == "failed" && !haskey(status, "child_pid")
        @test !ispath(joinpath(batch, "lock"))
        command = BatchControl.simulation_command(root, case)
        @test "--threads=1" in command.exec && "--gcthreads=1" in command.exec
    end
end

@testset "Existing default-disabled interlink and scheduler regression suites" begin
    include(joinpath(LIBRARY_ROOT, "test", "unit", "simulation", "interlink_tests.jl"))
    include(joinpath(LIBRARY_ROOT, "test", "unit", "simulation", "inclination_scheduler_tests.jl"))
end

@testset "Enabled R/T/N checkpoint round trip" begin
    mktempdir(RESULTS) do directory
        base = InterlinkTests.gve_sma_configuration(duration=20.0)
        enabled = SM.SimConfig._with_configuration(base; interlink_diagnostics=true)
        checkpointed = InterlinkTests.with_settings(enabled; results=false,
            results_directory=directory, checkpoint_enabled=true, checkpoint_interval_s=10.0)
        first_run = OracleStudy.run_simulation(checkpointed; return_results=true, return_solution=true)
        resumed_args = InterlinkTests.with_settings(SM.SimConfig._with_configuration(checkpointed;
            mission_configuration=InterlinkTests.gve_sma_configuration(duration=30.0).mission_configuration);
            resume_from_checkpoint=true)
        resumed = OracleStudy.run_simulation(resumed_args; return_results=true, return_solution=true)
        for index in eachindex(first_run.solution.u[end].sc), field in BRANCH_FIELDS
            @test getproperty(resumed.solution.prob.u0.sc[index], field) ==
                getproperty(first_run.solution.u[end].sc[index], field)
        end
        full_args = SM.SimConfig._with_configuration(
            InterlinkTests.gve_sma_configuration(duration=30.0); interlink_diagnostics=true)
        full = OracleStudy.run_simulation(full_args; return_solution=true)
        @test resumed.solution.u[end] ≈ full.u[end] atol=1e-5 rtol=1e-9
        for sc in resumed.solution.u[end].sc, (axis, component) in enumerate(COMPONENTS)
            @test getproperty(sc, Symbol("laser_dv_$(component)_positive")) +
                getproperty(sc, Symbol("laser_dv_$(component)_negative")) ≈ sc.laser_rtn[axis] atol=1e-12
        end
    end
end

atomic_toml(joinpath(RESULTS, "fixtures.toml"), Dict("passed" => true,
    "schema_version" => SCHEMA_VERSION, "case_count" => length(grid_cases()),
    "source_signature" => source_signature(source_hashes())))
