module SolverConfigValidationTests
# SolverConfig refuses solver names and counts the engine cannot run, and the
# environment adapter accepts exactly the same names, so a typo can no longer
# reach the final Tsit5 branch of the solver policy under another name.
using Test
using SpaceAGORA
const SE = SpaceAGORA.SimulationEngine
const SC = SpaceAGORA.SimulationModel.SimConfig

function argument_error(f)
    try
        f()
    catch e
        e isa ArgumentError && return e
        rethrow()
    end
    return nothing
end

env_getter(env) = (name, default) -> String(get(env, name, default))

@testset "SolverConfig validation" begin
    @testset "accepted settings construct" begin
        @test SolverConfig().solver_mode === :tsit5
        for mode in SC._SOLVER_MODES
            @test SolverConfig(solver_mode=mode).solver_mode === mode
        end
        for name in SC._SPLIT_IMEX_SOLVERS
            @test SolverConfig(split_imex_solver=name).split_imex_solver === name
        end
        for name in SC._MULTIRATE_SOLVERS
            cfg = SolverConfig(multirate_slow_solver=name, multirate_fast_solver=name)
            @test (cfg.multirate_slow_solver, cfg.multirate_fast_solver) === (name, name)
        end
        cfg = SolverConfig(maxiters=1, multirate_fast_substeps=1, symplectic_dt_s=0.5,
            gravity_backbone_dt_s=2, multirate_slow_dt_s=1.5)
        @test cfg.maxiters == 1
        @test cfg.gravity_backbone_dt_s === 2.0
        @test SolverConfig(maxiters=nothing, gravity_backbone_dt_s=nothing).maxiters === nothing
    end

    @testset "unknown typed solver_mode is refused" begin
        for typo in (:tsit, :Tsit5, :rk4, :auto, :split)
            err = argument_error(() -> SolverConfig(solver_mode=typo))
            @test err isa ArgumentError
            @test occursin("SolverConfig.solver_mode=$(repr(typo))", err.msg)
            for mode in SC._SOLVER_MODES
                @test occursin(repr(mode), err.msg)
            end
        end
    end

    @testset ":symplectic points to :gravity_backbone_split" begin
        err = argument_error(() -> SolverConfig(solver_mode=:symplectic))
        @test err isa ArgumentError
        @test occursin(":symplectic", err.msg)
        @test occursin(":gravity_backbone_split", err.msg)
        @test !(:symplectic in SC._SOLVER_MODES)
        for spelling in ("symplectic", "KahanLi8", " verlet ")
            env_err = argument_error(() -> SE._parse_solver_mode_sym(spelling))
            @test env_err isa ArgumentError
            @test occursin(":gravity_backbone_split", env_err.msg)
        end
    end

    @testset "unknown typed solver names are refused" begin
        for (field, bad) in ((:split_imex_solver, :kencarp5), (:split_imex_solver, :tsit5),
                             (:multirate_slow_solver, :vern9), (:multirate_fast_solver, :kencarp47))
            err = argument_error(() -> SolverConfig(; field => bad))
            @test err isa ArgumentError
            @test occursin("SolverConfig.$(field)=$(repr(bad))", err.msg)
        end
    end

    @testset "counts and steps must be positive" begin
        for (field, bad) in ((:maxiters, 0), (:maxiters, -5), (:multirate_fast_substeps, 0),
                             (:multirate_fast_substeps, -1), (:gravity_backbone_dt_s, 0.0),
                             (:multirate_slow_dt_s, -2.0), (:symplectic_dt_s, 0.0),
                             (:gravity_backbone_dt_s, NaN))
            err = argument_error(() -> SolverConfig(; field => bad))
            @test err isa ArgumentError
            @test occursin("SolverConfig.$(field) must be positive", err.msg)
        end
    end

    @testset "environment parser shares the accepted names" begin
        for mode in SC._SOLVER_MODES
            @test SE._parse_solver_mode_sym(String(mode)) === mode
        end
        for name in SC._SPLIT_IMEX_SOLVERS
            @test SE._parse_split_imex_solver_sym(uppercase(String(name))) === name
        end
        for name in SC._MULTIRATE_SOLVERS
            @test SE._parse_multirate_solver_sym(String(name), "SPACEAGORA_MULTIRATE_FAST_SOLVER") === name
        end
        @test all(in(SC._SOLVER_MODES), values(SE._SOLVER_MODE_ALIASES))
        @test all(in(SC._SPLIT_IMEX_SOLVERS), values(SE._SPLIT_IMEX_SOLVER_ALIASES))
        @test all(in(SC._MULTIRATE_SOLVERS), values(SE._MULTIRATE_SOLVER_ALIASES))
        @test SE._parse_solver_mode_sym("") === :tsit5
        @test SE._parse_solver_mode_sym("Gravity-Backbone") === :gravity_backbone_split

        err = argument_error(() -> SE._parse_solver_mode_sym("tsit"))
        @test occursin("SPACEAGORA_SOLVER_MODE='tsit'", err.msg)
        @test occursin(join(SC._SOLVER_MODES, ", "), err.msg)
        @test argument_error(() -> SE._parse_split_imex_solver_sym("tsit5")) isa ArgumentError
        @test argument_error(() -> SE._parse_multirate_solver_sym("kencarp47", "X")) isa ArgumentError

        # Typed -> environment -> typed: every accepted mode survives the
        # engine-config override round trip unchanged.
        for mode in SC._SOLVER_MODES
            overrides = SE._engine_env_overrides(SE.SimulationEngineConfig(solver=SolverConfig(solver_mode=mode)))
            @test SE._solver_config_from_env(env_getter(overrides)).solver_mode === mode
        end
    end

    @testset "every accepted solver name has a solver" begin
        for name in SC._SPLIT_IMEX_SOLVERS
            @test SE._split_imex_solver_spec(SolverConfig(split_imex_solver=name)).label isa String
        end
        for name in SC._MULTIRATE_SOLVERS
            @test SE._multirate_solver_spec_from_sym(name, "multirate_fast_solver").label isa String
        end
    end

    @testset "strict environment parsing still validates every knob" begin
        @test SE._solver_config_from_env(env_getter(Dict("SPACEAGORA_SYMPLECTIC_DT_S" => "4.0"))).symplectic_dt_s == 4.0
        for (name, bad) in ("SPACEAGORA_SOLVER_MODE" => "tsit", "SPACEAGORA_SOLVER_MODE" => "symplectic",
                            "SPACEAGORA_SPLIT_IMEX_SOLVER" => "kencarp5", "SPACEAGORA_SOLVER_MAXITERS" => "0",
                            "SPACEAGORA_MULTIRATE_FAST_SUBSTEPS" => "-1", "SPACEAGORA_SYMPLECTIC_DT_S" => "0.0")
            @test_throws ArgumentError SE._solver_config_from_env(env_getter(Dict(name => bad)))
        end
    end

    @testset "non-strict environment fallback keeps its default and warns" begin
        cfg = @test_logs (:warn, r"SPACEAGORA_SOLVER_MODE='tsitt'.*:tsit5") SE.simulation_engine_config_from_env(
            Dict("SPACEAGORA_SOLVER_MODE" => "tsitt"))
        @test cfg.solver.solver_mode === :tsit5
        cfg = @test_logs (:warn, r"gravity_backbone_split") SE.simulation_engine_config_from_env(
            Dict("SPACEAGORA_SOLVER_MODE" => "symplectic"))
        @test cfg.solver.solver_mode === :tsit5
        cfg = @test_logs (:warn, r"SPACEAGORA_MULTIRATE_FAST_SUBSTEPS") SE.simulation_engine_config_from_env(
            Dict("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS" => "0"))
        @test cfg.solver.multirate_fast_substeps == 8
        cfg = @test_logs SE.simulation_engine_config_from_env(Dict("SPACEAGORA_SOLVER_MODE" => "rodas"))
        @test cfg.solver.solver_mode === :rodas5p
        @test_throws ArgumentError SE.simulation_engine_config_from_env(
            Dict("SPACEAGORA_SOLVER_MODE" => "tsitt"); solver_strict=true)
    end
end

end # module SolverConfigValidationTests
