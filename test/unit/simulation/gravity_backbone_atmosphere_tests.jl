module GravityBackboneAtmosphereTests
# solver_mode=:gravity_backbone_split integrates position and velocity under
# gravity alone, and no density or thermal callback can read that state, so a
# run with a density model set is refused before it starts, however the mode is
# chosen. Without the refusal, SolverConfig(solver_mode=:gravity_backbone_split)
# failed at the first step (the density callback read `integrator.u.sc` from the
# backbone's ArrayPartition state), and SPACEAGORA_SOLVER_MODE=gravity_backbone_split
# ran with the atmosphere ignored and reported zero heating. The backbone's
# admission rules are checked before any run setup, so a resumed run already at
# its end, which solves nothing, is refused too.
using Test
using Logging
using OrdinaryDiffEq: DiscreteCallback
using SpaceAGORA
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine

const PLANET = SM.make_no_gram_planet(:earth)
const BACKBONE = SolverConfig(solver_mode=:gravity_backbone_split)

exponential_atmosphere() = SM.ExponentialAtmosphereModel(1.0e-11, 4.0e5, 6.0e4; temperature_k=1000.0,
    valid_min_altitude_m=0.0, valid_max_altitude_m=2.0e6)

function backbone_config(density_model; solver_config=BACKBONE,
        effectors=(SM.InverseSquaredJ2GravityModel(),))
    root = SM.Link(root=true, m=500.0, ref_area=4.0)
    spacecraft = SM.SpacecraftModel(SM.Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0,
        SM.InitialCondition(a=PLANET.Rp_e + 4.0e5, e=0.0, i=53.0, Ω=0.0, ω=0.0, ν=0.0), 1)
    return SM.make_example_config(planet=PLANET, spacecraft=spacecraft, mission_time=60.0,
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=effectors, density_model=density_model,
        ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=false, keplerian=true,
        verbose=false, results=false, results_directory=mktempdir(), solver_config=solver_config)
end

# The atmosphere-without-drag warning is expected for these configurations.
function run_error(args; solver_mode_env, kwargs...)
    withenv("SPACEAGORA_SOLVER_MODE" => solver_mode_env) do
        with_logger(NullLogger()) do
            try
                run_simulation(args; kwargs...)
            catch err
                return err
            end
            return nothing
        end
    end
end

refuses_density(msg) = occursin("gravity_backbone_split requires density_model=NoAtmosphereModel()", msg) &&
    occursin("ExponentialAtmosphereModel", msg)

@testset "gravity backbone and a density model" begin
    @testset "the admission check names the density model" begin
        args = backbone_config(exponential_atmosphere())
        reason = SE._gravity_backbone_reject_reason(args)
        @test reason !== nothing && refuses_density(reason)
        @test !SE._gravity_backbone_eligible(args)
    end

    @testset "chosen by SolverConfig under another process default" begin
        err = run_error(backbone_config(exponential_atmosphere()); solver_mode_env="tsit5")
        @test err isa ArgumentError && refuses_density(err.msg)
    end

    @testset "chosen by SPACEAGORA_SOLVER_MODE" begin
        args = backbone_config(exponential_atmosphere(); solver_config=nothing)
        err = run_error(args; solver_mode_env="gravity_backbone_split")
        @test err isa ArgumentError && refuses_density(err.msg)
    end

    @testset "a vacuum backbone run is still admitted" begin
        args = backbone_config(SM.NoAtmosphereModel())
        @test SE._gravity_backbone_reject_reason(args) === nothing
        sol = withenv("SPACEAGORA_SOLVER_MODE" => "tsit5") do
            run_simulation(args; return_solution=true)
        end
        @test sol.t[end] == 60.0
    end

    @testset "refused before any run setup" begin
        # A callback factory is called during setup, before the solver runs.
        calls = Ref(0)
        spy = _ -> (calls[] += 1; DiscreteCallback((u, t, i) -> false, i -> nothing; save_positions=(false, false)))
        err = run_error(backbone_config(exponential_atmosphere()); solver_mode_env="tsit5",
            extra_callback_factories=(spy,))
        @test err isa ArgumentError && refuses_density(err.msg)
        @test calls[] == 0
        # The existing rules are checked at the same point.
        drag = backbone_config(SM.NoAtmosphereModel();
            effectors=(SM.InverseSquaredJ2GravityModel(), SM.AerodynamicCoefficientConstant()))
        err = run_error(drag; solver_mode_env="tsit5", extra_callback_factories=(spy,))
        @test err isa ArgumentError && occursin("does not support AerodynamicCoefficientConstant", err.msg)
        @test calls[] == 0
    end

    @testset "a resumed run already at its end" begin
        mktempdir() do dir
            settings(resume) = SM.SimulationSettings(results=false, verbose=false, generate_plots=false,
                normalize=false, results_directory=dir, checkpoint_enabled=true, checkpoint_interval_s=30.0,
                checkpoint_directory=joinpath(dir, "checkpoints"), resume_from_checkpoint=resume)
            with_settings(args, resume) = SM.SimConfig._with_configuration(args; simulation_settings=settings(resume))
            @test run_error(with_settings(backbone_config(SM.NoAtmosphereModel()), false); solver_mode_env="tsit5") === nothing
            # The checkpoint is at the mission's end, so a resumed run has nothing left to solve.
            @test SE._load_checkpoint(with_settings(backbone_config(SM.NoAtmosphereModel()), true)).t == 60.0
            @test run_error(with_settings(backbone_config(SM.NoAtmosphereModel()), true); solver_mode_env="tsit5") === nothing
            err = run_error(with_settings(backbone_config(exponential_atmosphere()), true); solver_mode_env="tsit5")
            @test err isa ArgumentError && refuses_density(err.msg)
        end
    end
end

end
