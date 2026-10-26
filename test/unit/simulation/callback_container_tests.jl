# The container a run's callbacks reach the solver in: a tuple up to a limit, and
# beyond it the discrete callbacks in a vector, in the same order, beside the
# continuous ones. A run gives the same bits either way: the same final state,
# saved states and saved times.
module CallbackContainerTests
using Test, StaticArrays, LinearAlgebra, SpaceAGORA
import DiffEqBase
using DiffEqCallbacks: PresetTimeCallback
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const ZERO3 = SVector(0.0, 0.0, 0.0)

# The laser-interlink example's configuration: the engine wraps the run's
# callbacks in the interlink scheduler after building them.
module InterlinkExample
    include(joinpath(@__DIR__, "..", "..", "..", "examples", "gve_sma_interlinks.jl"))
end

# Solve with the tuple limit set to `limit`, restoring it afterwards.
function with_limit(f, limit)
    previous = SE._DISCRETE_CALLBACK_TUPLE_LIMIT[]
    SE._DISCRETE_CALLBACK_TUPLE_LIMIT[] = limit
    try
        return f()
    finally
        SE._DISCRETE_CALLBACK_TUPLE_LIMIT[] = previous
    end
end

noop() = DiffEqBase.DiscreteCallback((u, t, i) -> false, i -> nothing; save_positions=(false, false))
never() = DiffEqBase.ContinuousCallback((u, t, i) -> 1.0, i -> nothing; save_positions=(false, false))

# A controller bound to one spacecraft that holds a small push along its velocity.
mutable struct Push <: SM.AbstractControlEffectorModel
    spacecraft_index::Int
    held::SVector{3, Float64}
end
SpaceAGORA.update_control!(m::Push, u, p, t::Float64, i::Int) =
    (m.held = 1.0e-3 * normalize(SVector{3, Float64}(u.sc[i].vel)); nothing)
SpaceAGORA.control_wrench(m::Push, x::StateSample, p, t::Float64) =
    x.spacecraft_index == m.spacecraft_index ? (m.held, ZERO3) : (ZERO3, ZERO3)
SpaceAGORA.spacecraft_scope(m::Push) = BoundSpacecraft((m.spacecraft_index,))
SpaceAGORA.bind_spacecraft(m::Push, k::Int) = Push(k, ZERO3)
SpaceAGORA.initialize_control!(m::Push, u, p, t) = (m.held = ZERO3; nothing)

const PLANET = SM.make_no_gram_planet(:earth)

# `n` spacecraft, each declaring a bound pusher, with J2 and, when `drag`, an
# exponential atmosphere and constant-coefficient drag.
function run_args(n; drag=false, mission_time=600.0)
    sats = map(1:n) do id
        root = SM.Link(root=true, m=100.0, ref_area=2.0)
        ic = SM.InitialCondition(ra=PLANET.Rp_e + (drag ? 2000e3 : 700e3), rp=PLANET.Rp_e + (drag ? 140e3 : 650e3),
            i=45.0, ω=0.0, Ω=30.0id, ν=drag ? 160.0 : 10.0)
        SM.SpacecraftModel(joints=SM.Joint[], links=[root], root=root, prop_mass=10.0, inertia_tensor=root.inertia,
            initial_condition=ic, id=id, control=SM.ControlModel((Push(0, ZERO3),), [7.0]))
    end
    density = drag ? SM.ExponentialAtmosphereModel(1.0e-9, 1.4e5, 7.0e3; temperature_k=800.0,
        valid_min_altitude_m=0.0, valid_max_altitude_m=2.0e6) : SM.NoAtmosphereModel()
    effectors = drag ? (SM.InverseSquaredJ2GravityModel(), SM.AerodynamicCoefficientConstant()) :
        (SM.InverseSquaredJ2GravityModel(),)
    env = SM.EnvironmentModel(planet=PLANET, EI=200.0, density_model=density, ephemerides_model=SM.SimpleEphemeridesModel(),
        thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=PLANET), topography=false, wind=false)
    return SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=mission_time, orientation_sim=false, num_steps_to_save=100),
        environment_model=env,
        dynamics_model=SM.DynamicsModel(SM.SpacecraftModel[sats...], effectors),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
end

# The container the solver holds the discrete callbacks in, seen from a callback.
function solver_container(args, limit)
    seen = Ref{Any}(nothing)
    probe = DiffEqBase.DiscreteCallback((u, t, i) -> (seen[] === nothing && (seen[] = (i.opts.callback.continuous_callbacks,
        i.opts.callback.discrete_callbacks)); false), i -> nothing; save_positions=(false, false))
    with_limit(limit) do
        run_simulation(args; extra_callbacks=(probe,))
    end
    return seen[]
end

deactivate(i, t_s) = PresetTimeCallback([t_s], integrator -> (integrator.p.is_active[i] = false); save_positions=(false, false))

same_solution(a, b) = a.t == b.t && length(a.u) == length(b.u) && all(isequal(collect(x), collect(y)) for (x, y) in zip(a.u, b.u))

@testset "Callback container" begin
    @testset "a tuple up to the limit, then the discrete callbacks in a vector" begin
        @test SE._DISCRETE_CALLBACK_TUPLE_LIMIT[] == 32
        continuous = never()
        discrete = [noop() for _ in 1:5]
        set = DiffEqBase.CallbackSet(continuous, discrete...)
        with_limit(5) do
            @test SE._solver_callback_set(set) === set                       # at the limit: unchanged
        end
        held = with_limit(4) do
            SE._solver_callback_set(set)
        end
        @test held.continuous_callbacks === set.continuous_callbacks && held.continuous_callbacks isa Tuple
        @test held.discrete_callbacks isa Vector{Any}
        @test length(held.discrete_callbacks) == 5 && all(held.discrete_callbacks .=== discrete)   # same order
        @test with_limit(() -> SE._solver_callback_set(held), 0) === held         # a vector stays as it is
    end

    @testset "the solver holds the run's discrete callbacks in a vector beyond the limit" begin
        args = run_args(3)
        continuous, discrete = solver_container(args, typemax(Int))
        @test continuous isa Tuple && discrete isa Tuple
        continuous, discrete = solver_container(args, 0)
        @test continuous isa Tuple && discrete isa Vector{Any}
        # A run with interlinks couples its spacecraft, so it has no block-diagonal
        # Jacobian prototype and the solver library holds both kinds in vectors
        # itself, whatever the limit.
        for limit in (typemax(Int), 0)
            continuous, discrete = solver_container(InterlinkExample.gve_sma_configuration(duration=20.0), limit)
            @test continuous isa Vector{Any} && discrete isa Vector{Any}
        end
    end

    @testset "a cached integrator is reused with the vector" begin
        cache = SE.SolverIntegratorCache()
        reused = with_limit(0) do
            run_simulation(run_args(3); solver_cache=cache, return_solution=true)
            held = cache.integrator
            @test held.opts.callback.discrete_callbacks isa Vector{Any}
            solution = run_simulation(run_args(3); solver_cache=cache, return_solution=true)
            @test cache.integrator === held
            solution
        end
        fresh = with_limit(typemax(Int)) do
            run_simulation(run_args(3); return_solution=true)
        end
        @test same_solution(reused, fresh)
    end

    @testset "a run gives the same bits with either container" begin
        for (label, args, extra) in (("bound controllers, a spacecraft deactivated", () -> run_args(4),
                                      () -> (deactivate(2, 200.0),)),
                                     ("drag and atmosphere events", () -> run_args(2; drag=true, mission_time=900.0), () -> ()),
                                     ("laser interlinks", () -> InterlinkExample.gve_sma_configuration(duration=60.0), () -> ()))
            tuple_run = with_limit(typemax(Int)) do
                run_simulation(args(); return_solution=true, extra_callbacks=extra())
            end
            vector_run = with_limit(0) do
                run_simulation(args(); return_solution=true, extra_callbacks=extra())
            end
            @test same_solution(tuple_run, vector_run)
            @test length(tuple_run.t) > 2
        end
    end
end

end
