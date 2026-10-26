module CallbackSolverModeTests
# get_callbacks builds a run's callbacks for the solver mode the run solves with:
# the mode of its SolverConfig, or SPACEAGORA_SOLVER_MODE's when it has none, as
# run_simulation resolves it. It used to read the process setting only, so a run
# whose SolverConfig chose another mode than SPACEAGORA_SOLVER_MODE got the
# callbacks of the wrong mode: a gravity-backbone run chosen by SolverConfig got
# the first-order set, and a Tsit5 run chosen by SolverConfig under
# SPACEAGORA_SOLVER_MODE=gravity_backbone_split lost its planet-frame, GNC and
# result-saving callbacks.
using Test
using StaticArrays
using SpaceAGORA
const SM = SpaceAGORA.SimulationModel
const CB = SM.SimulationCallbacks
const CH = SM.ControlHooks

const PLANET = SM.make_no_gram_planet(:earth)
const BACKBONE = SolverConfig(solver_mode=:gravity_backbone_split)
const TSIT5 = SolverConfig(solver_mode=:tsit5)
const ENV_KEY = "SPACEAGORA_SOLVER_MODE"

mutable struct TickCounter <: SM.AbstractControlEffectorModel
    ticks::Int
end
CH.control_wrench(::TickCounter, x, p, t::Float64) =
    (SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0))
CH.update_control!(model::TickCounter, u, p, t::Float64, i::Int64) = (model.ticks += 1; nothing)

# One spacecraft under J2 gravity in vacuum for 60 s; the backbone admits it
# without controls.
function configuration(solver_config; controls=(), results=false)
    root = SM.Link(root=true, m=500.0, ref_area=4.0)
    spacecraft = SM.SpacecraftModel(SM.Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0,
        SM.InitialCondition(a=PLANET.Rp_e + 4.0e5, e=0.0, i=53.0, Ω=0.0, ω=0.0, ν=0.0), 1)
    args = SM.make_example_config(planet=PLANET, spacecraft=spacecraft, mission_time=60.0,
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(SM.InverseSquaredJ2GravityModel(),), density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=false, keplerian=true,
        verbose=false, results=results, results_directory=mktempdir(), solver_config=solver_config)
    isempty(controls) && return args
    return SM.SimConfig._with_configuration(args; control_model=SM.ControlModel(
        control_effectors=controls, control_rates=fill(10.0, length(controls))))
end

callbacks(args; kwargs...) = CB.get_callbacks(1, args.dynamics_model.dynamic_effectors, args; kwargs...)

# The types of the affects, in order, identify the callbacks a set holds.
signature(set) = (map(cb -> typeof(cb.affect!), set.continuous_callbacks),
    map(cb -> typeof(cb.affect!), set.discrete_callbacks))

const PLANET_FRAME = typeof(CB.update_planet_frame_callback().affect!)
updates_planet_frame(set) = any(cb -> typeof(cb.affect!) === PLANET_FRAME, set.discrete_callbacks)

# The set the default (no `solver_mode` keyword) builds, under `env_mode`.
default_set(args, env_mode; kwargs...) = withenv(() -> callbacks(args; kwargs...), ENV_KEY => env_mode)

@testset "callbacks follow the run's solver mode" begin
    @testset "the mode passed in decides" begin
        args = configuration(TSIT5; controls=(TickCounter(0),), results=true)
        backbone = callbacks(args; solver_mode=:gravity_backbone_split)
        first_order = callbacks(args; solver_mode=:tsit5)
        @test !updates_planet_frame(backbone)
        @test updates_planet_frame(first_order)
        # Control ticks and result saving are in the first-order set only.
        @test length(first_order.discrete_callbacks) > length(backbone.discrete_callbacks)
        @test signature(backbone) != signature(first_order)
        # The keyword wins over both the SolverConfig and the environment.
        @test signature(withenv(() -> callbacks(args; solver_mode=:gravity_backbone_split),
            ENV_KEY => "tsit5")) == signature(backbone)
    end

    @testset "chosen by SolverConfig only" begin
        for (solver_config, mode) in ((BACKBONE, :gravity_backbone_split), (TSIT5, :tsit5))
            args = configuration(solver_config)
            expected = signature(callbacks(args; solver_mode=mode))
            @test signature(default_set(args, nothing)) == expected
        end
        @test !updates_planet_frame(default_set(configuration(BACKBONE), nothing))
        @test updates_planet_frame(default_set(configuration(TSIT5), nothing))
    end

    @testset "chosen by SPACEAGORA_SOLVER_MODE only" begin
        args = configuration(nothing)
        for (env_mode, mode) in (("gravity_backbone_split", :gravity_backbone_split), ("tsit5", :tsit5))
            @test signature(default_set(args, env_mode)) == signature(callbacks(args; solver_mode=mode))
        end
        @test !updates_planet_frame(default_set(args, "gravity_backbone_split"))
        @test updates_planet_frame(default_set(args, "tsit5"))
    end

    @testset "SolverConfig and SPACEAGORA_SOLVER_MODE disagree: SolverConfig wins" begin
        backbone_args = configuration(BACKBONE; results=true)
        set = default_set(backbone_args, "tsit5"; record_saved_values=true)
        @test signature(set) == signature(callbacks(backbone_args; solver_mode=:gravity_backbone_split,
            record_saved_values=true))
        @test !updates_planet_frame(set)

        tsit5_args = configuration(TSIT5; controls=(TickCounter(0),), results=true)
        set = default_set(tsit5_args, "gravity_backbone_split"; record_saved_values=true)
        @test signature(set) == signature(callbacks(tsit5_args; solver_mode=:tsit5, record_saved_values=true))
        @test updates_planet_frame(set)
    end

    @testset "run_simulation installs the callbacks of the mode it solves with" begin
        # A Tsit5 run chosen by SolverConfig under a backbone environment ticks its
        # control and saves its results, as it does with the environment unset.
        function tsit5_run(env_mode)
            withenv(ENV_KEY => env_mode) do
                run_simulation(configuration(TSIT5; controls=(TickCounter(0),)); return_results=true)
            end
        end
        reference, disagreeing = tsit5_run(nothing), tsit5_run("gravity_backbone_split")
        ticks(res) = only(res.configuration.control_model.control_effectors).ticks
        @test ticks(reference) > 0
        @test ticks(disagreeing) == ticks(reference)
        @test size(reference.table, 1) > 1
        @test isequal(disagreeing.table, reference.table)

        # A backbone run chosen by SolverConfig under the Tsit5 default runs as
        # one chosen by SPACEAGORA_SOLVER_MODE.
        function backbone_run(solver_config, env_mode)
            withenv(ENV_KEY => env_mode) do
                run_simulation(configuration(solver_config); return_results=true, return_solution=true)
            end
        end
        by_config = backbone_run(BACKBONE, "tsit5")
        by_environment = backbone_run(nothing, "gravity_backbone_split")
        @test by_config.solution.t[end] == 60.0
        @test by_config.solution.t == by_environment.solution.t
        @test by_config.solution.u == by_environment.solution.u
        @test isequal(by_config.table, by_environment.table)
    end
end

end
