module EffectorScopeTests

# Dynamic effectors acting on some spacecraft only (`scoped`): resolution and its
# refusals, the evaluation on every right-hand-side route, the bits against the
# workaround of a model that returns zero outside its spacecraft, the
# accelerometer, the aerodynamic diagnostics, the constellation ensemble and
# the checkpoint entry.

using Test
using LinearAlgebra
using StaticArrays
using RecursiveArrayTools: ArrayPartition
using Logging
using OrdinaryDiffEq: DiscreteCallback
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const Z3 = SVector(0.0, 0.0, 0.0)

# --- Models -------------------------------------------------------------------

# A small along-track push, with an explicit-kick form for the gravity backbone.
struct Push <: SpaceAGORA.AbstractForceTorqueModel
    newtons::Float64
end
SpaceAGORA.wrench(m::Push, x::StateSample, env::EnvironmentSample, t::Float64) = (m.newtons * normalize(x.vel_ii), Z3)
SpaceAGORA.gravity_backbone_kick_structure(::Push) = :velocity_kick_explicit
SpaceAGORA.gravity_backbone_kick_acceleration_ii(m::Push, x::StateSample, env::EnvironmentSample, t::Float64) =
    m.newtons * normalize(x.vel_ii) / x.mass_kg

# The workaround a scope replaces: the model selects its spacecraft from the sample.
struct GatedPush <: SpaceAGORA.AbstractForceTorqueModel
    newtons::Float64
    indices::Tuple{Vararg{Int}}
end
SpaceAGORA.wrench(m::GatedPush, x::StateSample, env::EnvironmentSample, t::Float64) =
    x.spacecraft_index in m.indices ? (m.newtons * normalize(x.vel_ii), Z3) : (Z3, Z3)

# An atmosphere that counts the density queries made of it.
mutable struct CountingAtmosphere <: SM.AbstractDensityModel
    queries::Int
end
function SM.EnvironmentModels.getDensity(m::CountingAtmosphere, h::Float64, lat::Float64, lon::Float64, t::Float64,
                                         wind::Bool, p)
    m.queries += 1
    return 1.0e-11, 800.0, SVector(0.0, 0.0, 0.0)
end

# A push that binds itself to a run position by its spacecraft_index field.
struct BoundPush <: SpaceAGORA.AbstractForceTorqueModel
    spacecraft_index::Int
    newtons::Float64
end
SpaceAGORA.wrench(m::BoundPush, x::StateSample, env::EnvironmentSample, t::Float64) =
    x.spacecraft_index == m.spacecraft_index ? (m.newtons * normalize(x.vel_ii), Z3) : (Z3, Z3)

# Central gravity of the user's own, so a scoped built-in leaves every spacecraft some gravity.
struct PointMass <: SpaceAGORA.AbstractForceTorqueModel end
SpaceAGORA.is_gravitational(::PointMass) = true
function SpaceAGORA.wrench(::PointMass, x::StateSample, env::EnvironmentSample, t::Float64)
    r = norm(x.pos_ii)
    return (-3.986004418e14 * x.mass_kg / r^3 * x.pos_ii, Z3)
end

# A uniform field of the user's own; the accelerometer does not feel gravity.
struct Uniform <: SpaceAGORA.AbstractForceTorqueModel end
SpaceAGORA.is_gravitational(::Uniform) = true
SpaceAGORA.wrench(::Uniform, x::StateSample, env::EnvironmentSample, t::Float64) = (x.mass_kg * SVector(0.0, 0.0, -0.5), Z3)
struct GatedUniform <: SpaceAGORA.AbstractForceTorqueModel
    indices::Tuple{Vararg{Int}}
end
SpaceAGORA.is_gravitational(::GatedUniform) = true
SpaceAGORA.wrench(m::GatedUniform, x::StateSample, env::EnvironmentSample, t::Float64) =
    x.spacecraft_index in m.indices ? (x.mass_kg * SVector(0.0, 0.0, -0.5), Z3) : (Z3, Z3)

# Probes recording the spacecraft each route evaluates them for.
struct IndexProbe <: SpaceAGORA.AbstractForceTorqueModel
    seen::Vector{Int}
    lock::ReentrantLock
    partition::Symbol
end
IndexProbe(partition::Symbol=:explicit) = IndexProbe(Int[], ReentrantLock(), partition)
record!(m, x) = lock(() -> push!(m.seen, x.spacecraft_index), m.lock)
SpaceAGORA.wrench(m::IndexProbe, x::StateSample, env::EnvironmentSample, t::Float64) = (record!(m, x); (Z3, Z3))
SpaceAGORA.solver_partition(m::IndexProbe) = m.partition
SpaceAGORA.gravity_backbone_kick_structure(::IndexProbe) = :velocity_kick_explicit
SpaceAGORA.gravity_backbone_kick_acceleration_ii(m::IndexProbe, x::StateSample, env::EnvironmentSample, t::Float64) =
    (record!(m, x); Z3)

struct GravityProbe <: SpaceAGORA.AbstractForceTorqueModel
    seen::Vector{Int}
    lock::ReentrantLock
end
GravityProbe() = GravityProbe(Int[], ReentrantLock())
SpaceAGORA.is_gravitational(::GravityProbe) = true
SpaceAGORA.gravity_backbone_structure(::GravityProbe) = :position_only_static_gravity
SpaceAGORA.gravity_backbone_acceleration_ii(m::GravityProbe, x::StateSample, env::EnvironmentSample, t::Float64) =
    (record!(m, x); Z3)
SpaceAGORA.wrench(::GravityProbe, x::StateSample, env::EnvironmentSample, t::Float64) = (Z3, Z3)

# An effector type that declares its own scope, and one whose declaration is not a scope.
struct OwnScope <: SpaceAGORA.AbstractForceTorqueModel end
SpaceAGORA.spacecraft_scope(::OwnScope) = BoundSpacecraft((2,))
SpaceAGORA.wrench(::OwnScope, x::StateSample, env::EnvironmentSample, t::Float64) = (SVector(1e-3, 0.0, 0.0), Z3)
struct BadScope <: SpaceAGORA.AbstractForceTorqueModel end
SpaceAGORA.spacecraft_scope(::BadScope) = :everyone
SpaceAGORA.wrench(::BadScope, x::StateSample, env::EnvironmentSample, t::Float64) = (Z3, Z3)

# A control model whose scope is a selection, which GNC does not take.
struct SelectingControl <: SpaceAGORA.AbstractControlEffectorModel end
SpaceAGORA.spacecraft_scope(::SelectingControl) = SpacecraftSelection(planes=1)

# An externally propagated spacecraft (only its ownership is read before the refusal).
struct OwnsSecond <: SM.ExternalPropagation.AbstractExternalPropagator end
SM.ExternalPropagation.external_spacecraft(::OwnsSecond) = [2]

# --- Configurations -----------------------------------------------------------

const PLANET = SM.make_no_gram_planet(:earth)
const T0 = SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0)
const DENSITY = SM.ExponentialAtmosphereModel(1.0e-11, 4.0e5, 6.0e4; temperature_k=1000.0,
    valid_min_altitude_m=0.0, valid_max_altitude_m=2.0e6)

function walker4(; orientation=false)
    factory = (plane, slot, index) -> begin
        root = SM.Link(root=true, m=500.0, ref_area=4.0)
        SM.SpacecraftModel(SM.Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0,
            SM.InitialCondition(a=PLANET.Rp_e + 4.0e5, e=0.0, i=53.0, Ω=0.0, ω=0.0, ν=0.0), index)
    end
    # Walker 4/2/1: positions 1-2 in plane 1, 3-4 in plane 2.
    return walker_constellation(factory; total=4, planes=2, phasing=1, altitude=4.0e5, planet=PLANET, i=53)
end

function configuration(effectors; solver_mode=:tsit5, density=SM.NoAtmosphereModel(), mission_time=300.0,
                       spacecraft=walker4(), results=false)
    base = SM.make_example_config(planet=PLANET, spacecraft=spacecraft[1], mission_time=mission_time, initial_time=T0,
        dynamic_effectors=effectors, density_model=density, ephemerides_model=SM.SimpleEphemeridesModel(),
        orientation_sim=false, keplerian=true, verbose=false, results=results, results_directory=mktempdir(),
        solver_config=SM.SolverConfig(solver_mode=solver_mode))
    return SM.SimConfig._with_configuration(base; dynamics_model=SM.DynamicsModel(spacecraft, effectors))
end

function solve(args; route="serial", kwargs...)
    sol = withenv("SPACEAGORA_RHS_EXECUTION_MODE" => route, "SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
        run_simulation(args; return_solution=true, kwargs...)
    end
    return sol isa Tuple ? last(sol) : sol
end
same_run(a, b) = isequal(a.t, b.t) && all(isequal(collect(x), collect(y)) for (x, y) in zip(a.u, b.u))

function params(args, mode)
    p = SM.ODEParams(n_sats=length(args.dynamics_model.spacecraft), args=args)
    p.is_active .= true
    p.shared_buffers.rhs_env_config[] = withenv(SE._snapshot_rhs_plan_env_config,
        "SPACEAGORA_RHS_EXECUTION_MODE" => mode)
    return p
end

const FLAT_PLAN = (
    mode=:flat_constellation_effector_queue, allotment=1, scheduler=:static,
    dominant_axis=:effector, policy_applied=false,
    effector_decision=(use_threads=false, allotment=1, mode=:off, policy_applied=false),
)

# One derivative evaluation of `args` on `route`.
function derivative(args, route)
    u = SE.build_initial_conditions(args)
    du = zero(u)
    if route in ("serial", "satellite")
        SE.spacecraft_dynamics!(du, u, params(args, route), 0.0)
    elseif route == "flat"
        SE._spacecraft_dynamics_flat_constellation_effector_queue!(du, u, params(args, "auto"), 0.0, FLAT_PLAN;
            rhs_kind=:full, partition=nothing)
    elseif route == "explicit_remainder"
        SE.spacecraft_dynamics_explicit_remainder!(du, u, params(args, "auto"), 0.0)
    elseif route == "multirate_slow"
        SE.spacecraft_dynamics_slow!(du, u, params(args, "auto"), 0.0)
    end
    return du
end

@testset "scoped and its resolution" begin
    push = Push(1e-3)
    s = scoped(push, BoundSpacecraft((1, 3)))
    @test s isa ScopedEffector && s.effector === push && s.scope == BoundSpacecraft((1, 3))
    @test SpaceAGORA.spacecraft_scope(s) == BoundSpacecraft((1, 3))
    args = configuration((InverseSquaredGravityModel(), s))
    r = SE._resolve_effector_scopes(args)
    @test r.dynamics_model.dynamic_effectors[2] === push            # the effector itself, not a copy
    @test r.dynamics_model.effector_scopes.positions == [nothing, [1, 3]]
    @test r.dynamics_model.effector_scopes.masks[2] == BitVector([1, 0, 1, 0])
    @test SE._resolve_effector_scopes(r) === r                       # idempotent
    # No scope: the configuration itself, with no table.
    plain = configuration((InverseSquaredGravityModel(), push))
    @test SE._resolve_effector_scopes(plain) === plain
    @test plain.dynamics_model.effector_scopes === nothing
    # A scope that names every spacecraft is the same as none.
    for all_of_them in (AllSpacecraft(), SpacecraftSelection(planes=1:2), BoundSpacecraft((1, 2, 3, 4)))
        full = SE._resolve_effector_scopes(configuration((InverseSquaredGravityModel(), scoped(push, all_of_them))))
        @test full.dynamics_model.effector_scopes === nothing
        @test full.dynamics_model.dynamic_effectors[2] === push
    end
    # A type's own declaration is honoured without scoped.
    own = SE._resolve_effector_scopes(configuration((InverseSquaredGravityModel(), OwnScope())))
    @test own.dynamics_model.effector_scopes.positions == [nothing, [2]]
    # GNC flattening keeps a resolved table (positions do not move).
    @test SE._flatten_spacecraft_gnc(r) === r
end

@testset "refusals" begin
    push = Push(1e-3)
    resolve(effs; kw...) = SE._resolve_effector_scopes(configuration(effs; kw...))
    @test_throws "scoped takes a scope" scoped(push, :some)
    @test_throws "does not take a CoupledSpacecraft" scoped(push, CoupledSpacecraft(1, (2,)))
    @test_throws "already scoped" scoped(scoped(push, BoundSpacecraft((1,))), BoundSpacecraft((1,)))
    @test_throws "already declares its scope" scoped(OwnScope(), BoundSpacecraft((1,)))
    @test_throws "has the scope :everyone" resolve((InverseSquaredGravityModel(), BadScope()))
    @test_throws "CoupledSpacecraft scope" resolve((InverseSquaredGravityModel(),
        ScopedEffector(push, CoupledSpacecraft(1, (2,)))))
    @test_throws "outside the run's positions 1 to 4" resolve((InverseSquaredGravityModel(), scoped(push, BoundSpacecraft((5,)))))
    @test_throws "id 9, which no spacecraft" resolve((InverseSquaredGravityModel(), scoped(push, SpacecraftSelection(ids=9))))
    twins = walker4(); twins[4].id = 1
    @test_throws "positions 1, 4 share" resolve((InverseSquaredGravityModel(), scoped(push, SpacecraftSelection(ids=1)));
        spacecraft=twins)
    @test_throws "plane 3, which no spacecraft" resolve((InverseSquaredGravityModel(), scoped(push, SpacecraftSelection(planes=3))))
    @test_throws "returned nothing" resolve((InverseSquaredGravityModel(),
        scoped(push, SpacecraftSelection(where=(i, s) -> nothing))))
    @test_throws "selects no spacecraft of the run" resolve((InverseSquaredGravityModel(),
        scoped(push, SpacecraftSelection(planes=1, ids=[3]))))
    @test_throws "without any gravitational effector" resolve((scoped(InverseSquaredGravityModel(), BoundSpacecraft((1,))),))
    @test resolve((PointMass(), scoped(InverseSquaredGravityModel(), BoundSpacecraft((1,))))).dynamics_model.effector_scopes !== nothing
    @test_throws "picks its spacecraft itself" resolve((InverseSquaredGravityModel(),
        scoped(SM.RobotArmReactionEffector(), BoundSpacecraft((1,)))))
    external = SM.SimConfig._with_configuration(configuration((InverseSquaredGravityModel(), scoped(push, BoundSpacecraft((2, 3)))));
        external_propagators=(OwnsSecond(),))
    @test_throws "propagated externally" SE._resolve_effector_scopes(external)
    # EDG's prediction reads the run's gravity models for every spacecraft.
    edg = SM.AerobrakingEnergyDepletionGuidanceModel(SM.AerobrakingEnergyDepletionConfig(),
        SM.AerobrakingEnergyDepletionState(num_sats=4))
    with_edg(effs) = SM.SimConfig._with_configuration(configuration(effs);
        guidance_model=SM.GuidanceModel(guidance_effectors=(edg,), guidance_rates=[1.0]))
    @test_throws "EDG guidance" SE._resolve_effector_scopes(with_edg((PointMass(),
        scoped(InverseSquaredJ2GravityModel(), BoundSpacecraft((1,))))))
    @test SE._resolve_effector_scopes(with_edg((InverseSquaredGravityModel(), scoped(push, BoundSpacecraft((1,)))))) isa
        SimulationConfiguration
    # An algorithm guidance of the EDG family predicts the same way; one of another family does not.
    with_binding(g) = SM.SimConfig._with_configuration(configuration((PointMass(),
            scoped(InverseSquaredJ2GravityModel(), BoundSpacecraft((1,)))));
        guidance_model=SM.GuidanceModel(guidance_effectors=(g,), guidance_rates=[1.0]))
    @test_throws "EDG guidance" SE._resolve_effector_scopes(with_binding(
        AlgorithmGuidance(algorithm=EDGGuidanceAlgorithm(), spacecraft_index=1)))
    @test SE._resolve_effector_scopes(with_binding(
        AlgorithmGuidance(algorithm=ApoapsisBurnPlanner(delta_v_mps=1.0), spacecraft_index=1))) isa SimulationConfiguration
    # A GNC model takes a SpacecraftSelection too (from 0.7.0), resolved when the run starts.
    gnc = SM.SimConfig._with_configuration(configuration((InverseSquaredGravityModel(),));
        control_model=SM.ControlModel(control_effectors=(SelectingControl(),), control_rates=[1.0]))
    @test SM.SimulationCallbacks.check_spacecraft_scopes(gnc) === nothing
    # A run set up by hand must resolve scoped entries first.
    unresolved = configuration((InverseSquaredGravityModel(), scoped(push, BoundSpacecraft((1,)))))
    u0 = SE.build_initial_conditions(unresolved)
    p = SM.ODEParams(n_sats=4, args=unresolved)
    problems = SE._missing_model_hooks(unresolved, SE._hook_call_types(u0, p))
    @test any(occursin("dynamic_effectors[2] is a scoped entry", line) for line in problems)
    # Checks that read effectors by type see the scoped effector itself.
    @test_throws "does not support solar-radiation-pressure" solve(configuration((InverseSquaredGravityModel(),
        scoped(SM.SolarRadiationPressureModel(1.8, 10.0), BoundSpacecraft((1,))))))
    # Conformance checks the effector itself.
    @test_throws "checks the effector itself" SpaceAGORA.Conformance.check_force_torque_model(scoped(push, BoundSpacecraft((1,))))
end

@testset "a scoped effector is evaluated only in its scope, on every route" begin
    for route in ("serial", "satellite", "flat", "explicit_remainder", "multirate_slow")
        probe = IndexProbe()
        args = SE._resolve_effector_scopes(configuration((InverseSquaredGravityModel(),
            scoped(probe, BoundSpacecraft((1, 3))))))
        derivative(args, route)
        @test sort!(unique(probe.seen)) == [1, 3]
        # Control: unscoped, the probe sees every spacecraft on that route.
        control = IndexProbe()
        derivative(configuration((InverseSquaredGravityModel(), control)), route)
        @test sort!(unique(control.seen)) == [1, 2, 3, 4]
    end
    # The per-effector thread paths (taken for thread-safe built-ins when the policy decides).
    threads = (use_threads=true, allotment=2, mode=:threads, policy_applied=false)
    probe, other = IndexProbe(), IndexProbe()
    args = SE._resolve_effector_scopes(configuration((InverseSquaredGravityModel(), scoped(probe, BoundSpacecraft((2,))), other)))
    u = SE.build_initial_conditions(args)
    p = params(args, "serial")
    for k in 1:4
        forces, torques = MVector(0.0, 0.0, 0.0), MVector(0.0, 0.0, 0.0)
        SE._accumulate_dynamic_effectors!(forces, torques, u.sc[k], p, k, 0.0, args.dynamics_model.dynamic_effectors, threads)
        SE._accumulate_dynamic_effectors_partitioned!(forces, torques, u.sc[k], p, k, 0.0,
            args.dynamics_model.dynamic_effectors, threads, :explicit)
    end
    @test unique(probe.seen) == [2]
    @test sort!(unique(other.seen)) == [1, 2, 3, 4]
    # The implicit side of split IMEX.
    probe = IndexProbe(:implicit)
    args = SE._resolve_effector_scopes(configuration((InverseSquaredGravityModel(), scoped(probe, SpacecraftSelection(planes=2)));
        density=DENSITY))
    u = SE.build_initial_conditions(args)
    SE.spacecraft_dynamics_implicit_atmosphere!(zero(u), u, params(args, "auto"), 0.0)
    @test isempty(probe.seen) || sort!(unique(probe.seen)) == [3, 4]
    # The gravity backbone's core and kicks.
    gravity, kick = GravityProbe(), IndexProbe()
    bargs = SE._resolve_effector_scopes(configuration((InverseSquaredGravityModel(), scoped(gravity, BoundSpacecraft((2, 4))),
        scoped(kick, BoundSpacecraft((1,)))); solver_mode=:gravity_backbone_split))
    bu = SE.build_initial_conditions(bargs)
    bp = params(bargs, "auto")
    q, dq = SE._gravity_backbone_initial_states(bu, bargs)
    SE.spacecraft_dynamics_gravity_backbone!(zero(dq), dq, q, bp, 0.0)
    @test sort!(unique(gravity.seen)) == [2, 4]
    SE._gravity_backbone_half_kick!(ArrayPartition(deepcopy(dq), deepcopy(q)), bp, 0.0, 1.0)
    @test unique(kick.seen) == [1]
end

@testset "each spacecraft gets the bits of its own effectors" begin
    # A batchable built-in scoped next to gravity of the user's own (the flat
    # route's batch kernels), a custom push scoped by plane.
    entries = (PointMass(), scoped(InverseSquaredGravityModel(), BoundSpacecraft((1, 3))),
               scoped(Push(1e-3), SpacecraftSelection(planes=2)), Push(2e-4))
    resolved = SE._resolve_effector_scopes(configuration(entries))
    in_scope = [[1, 2, 4], [1, 4], [1, 2, 3, 4], [1, 3, 4]]
    for route in ("serial", "satellite", "flat", "explicit_remainder", "multirate_slow")
        du = derivative(resolved, route)
        for k in 1:4
            own = Tuple(SE._bare_effector(entries[j]) for j in in_scope[k])
            reference = derivative(configuration(own), route)
            @test isequal(collect(du.sc[k]), collect(reference.sc[k]))
        end
    end
end

@testset "a scoped run gives the bits of the gated workaround" begin
    gated = (InverseSquaredJ2GravityModel(), SM.AerodynamicCoefficientConstant(), GatedPush(1e-3, (1, 3)))
    for (solver_mode, route) in ((:tsit5, "serial"), (:tsit5, "flat"), (:split_imex, "satellite"))
        reference = solve(configuration(gated; solver_mode=solver_mode, density=DENSITY); route=route)
        for scope in (BoundSpacecraft((1, 3)), SpacecraftSelection(ids=[1, 3]), SpacecraftSelection(slots=[(1, 1), (2, 1)]))
            effs = (InverseSquaredJ2GravityModel(), SM.AerodynamicCoefficientConstant(), scoped(Push(1e-3), scope))
            @test same_run(solve(configuration(effs; solver_mode=solver_mode, density=DENSITY); route=route), reference)
        end
        # Control: the gated run differs from the run with the push everywhere.
        unscoped = (InverseSquaredJ2GravityModel(), SM.AerodynamicCoefficientConstant(), Push(1e-3))
        everywhere = solve(configuration(unscoped; solver_mode=solver_mode, density=DENSITY); route=route)
        @test !same_run(everywhere, reference)
        covering = (InverseSquaredJ2GravityModel(), SM.AerodynamicCoefficientConstant(), scoped(Push(1e-3), SpacecraftSelection(planes=1:2)))
        @test same_run(solve(configuration(covering; solver_mode=solver_mode, density=DENSITY); route=route), everywhere)
    end
    # Routes and threads change timing only.
    effs = (InverseSquaredJ2GravityModel(), scoped(Push(1e-3), BoundSpacecraft((2,))))
    @test same_run(solve(configuration(effs); route="serial"), solve(configuration(effs); route="satellite"))
end

@testset "aerodynamic diagnostics outside the scope stay zero" begin
    effs = (InverseSquaredJ2GravityModel(), scoped(SM.AerodynamicCoefficientConstant(), SpacecraftSelection(planes=1)))
    res = withenv("SPACEAGORA_RHS_EXECUTION_MODE" => "serial") do
        run_simulation(configuration(effs; density=DENSITY); return_results=true)
    end
    drag(i) = [getproperty(row, Symbol("sc$(i)_drag_$(c)")) for row in eachrow(res.table) for c in 1:3]
    @test all(iszero, drag(3)) && all(iszero, drag(4))
    @test any(!iszero, drag(1)) && any(!iszero, drag(2))
end

@testset "the per-spacecraft atmosphere pre-sample samples only spacecraft an atmosphere consumer acts on" begin
    function presampled(drag)
        atmosphere = CountingAtmosphere(0)
        args = SE._resolve_effector_scopes(configuration((InverseSquaredJ2GravityModel(), drag); density=atmosphere))
        p = params(args, "auto")
        # An inactive spacecraft rules out the constellation batch, so the flat
        # route's pre-sample takes its per-spacecraft path.
        p.is_active[4] = false
        @test SE._uniform_light_density_model(p, 4) === nothing
        @test !SM.SimulationCallbacks._rhs_density_service_candidate(p, 4)
        SE._prefill_environment_samples!(p, 0.0, SE.build_initial_conditions(args).sc; atmosphere=true)
        return atmosphere.queries, copy(p.shared_buffers.densities)
    end
    queries, densities = presampled(scoped(SM.AerodynamicCoefficientConstant(), SpacecraftSelection(planes=1)))
    @test queries == 2
    @test densities[1] > 0 && densities[2] > 0 && densities[3] == 0
    # Control: unscoped, the pre-sample also samples spacecraft 3 (4 is inactive).
    queries, densities = presampled(SM.AerodynamicCoefficientConstant())
    @test queries == 3
    @test densities[3] > 0
end

@testset "the accelerometer subtracts only the gravity applied" begin
    function sensed(effs)
        sc = walker4()
        acc = Accelerometer(period_s=1.0)
        root = sc[2].root
        add_sensor!(sc[2], root, acc)
        base = configuration(effs; spacecraft=sc, mission_time=5.0)
        args = SM.SimConfig._with_configuration(base; mission_configuration=SM.MissionConfiguration(
            mission_type=SM.MissionTime, keplerian=true, number_of_orbits=1, mission_time=5.0, orientation_sim=true,
            num_steps_to_save=10))
        solve(args; isolate_state=false)
        return sensor_measurement(acc)
    end
    scoped_read = sensed((InverseSquaredGravityModel(), scoped(Uniform(), BoundSpacecraft((1,)))))
    @test norm(scoped_read) < 1e-7                                     # free fall: spacecraft 2 feels nothing
    @test isequal(scoped_read, sensed((InverseSquaredGravityModel(), GatedUniform((1,)))))
end

@testset "ensemble members keep the effectors whose scope holds them" begin
    entries = (InverseSquaredJ2GravityModel(), scoped(Push(1e-3), BoundSpacecraft((1, 3))))
    args = configuration(entries)
    members = run_constellation_ensemble(args; return_solution=true)
    for k in 1:4
        member = members.successful[k].value
        member = member isa Tuple ? last(member) : member
        own = k in (1, 3) ? (InverseSquaredJ2GravityModel(), Push(1e-3)) : (InverseSquaredJ2GravityModel(),)
        single = configuration(own; spacecraft=[walker4()[k]])
        @test same_run(member, solve(single))
    end
    # Control: member 2 without the push differs from member 2 with it.
    @test !same_run(solve(configuration((InverseSquaredJ2GravityModel(),); spacecraft=[walker4()[2]])),
                    solve(configuration((InverseSquaredJ2GravityModel(), Push(1e-3)); spacecraft=[walker4()[2]])))
    @test SpaceAGORA.SimulationCampaigns._ensemble_member_effectors(SE._resolve_effector_scopes(args), 2) ==
        (InverseSquaredJ2GravityModel(),)
    # The ensemble's own checks see the effectors themselves: one bound to a run
    # position by its own field is refused inside a scope as without one.
    @test_throws "binds it to member 2" run_constellation_ensemble(configuration((InverseSquaredJ2GravityModel(),
        scoped(BoundPush(2, 1e-3), BoundSpacecraft((2,))))))
    @test_throws "binds it to member 2" run_constellation_ensemble(configuration((InverseSquaredJ2GravityModel(),
        BoundPush(2, 1e-3))))
end

@testset "a checkpoint keeps the scopes, and a resume needs the same" begin
    struct Interrupted <: Exception end
    armed = Ref(false)
    interrupter = DiscreteCallback((u, t, integrator) -> armed[] && t > 150.0, integrator -> throw(Interrupted());
        save_positions=(false, false))
    function checkpointed(effs, directory; resume=false)
        base = configuration(effs)
        settings = SM.SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false,
            save_csv=false, checkpoint_enabled=true, checkpoint_interval_s=150.0,
            checkpoint_directory=String(directory), resume_from_checkpoint=resume)
        return SM.SimConfig._with_configuration(base; simulation_settings=settings)
    end
    run!(args) = with_logger(() -> solve(args; extra_callbacks=(interrupter,)), ConsoleLogger(stderr, Logging.Error))
    function interrupted(effs)
        directory = mktempdir()
        armed[] = true
        try
            run!(checkpointed(effs, directory))
        catch
        finally
            armed[] = false
        end
        return directory
    end
    entries = (InverseSquaredJ2GravityModel(), scoped(Push(1e-3), BoundSpacecraft((1, 3))))
    straight = run!(checkpointed(entries, mktempdir()))
    directory = interrupted(entries)
    checkpoint = SE._load_checkpoint(checkpointed(entries, directory))
    @test checkpoint.t == 150.0
    entry = checkpoint.runtime_state.effector_scopes
    @test entry.version == 1 && [e.positions for e in entry.effectors] == [nothing, [1, 3]]
    @test entry.effectors[2].kind.name === :Push
    resumed = run!(checkpointed(entries, directory; resume=true))
    @test same_run(resumed, straight)                                   # the second half, bit for bit
    # The same scope written another way resumes too (from a checkpoint not yet resumed to the end).
    again = interrupted(entries)
    @test same_run(run!(checkpointed((InverseSquaredJ2GravityModel(), scoped(Push(1e-3), SpacecraftSelection(ids=[1, 3]))),
        again; resume=true)), straight)
    directory = interrupted(entries)
    @test_throws "acts on run positions 1, 2 here and on run positions 1, 3" run!(checkpointed(
        (InverseSquaredJ2GravityModel(), scoped(Push(1e-3), BoundSpacecraft((1, 2)))), directory; resume=true))
    @test_throws "but this run has none" run!(checkpointed((InverseSquaredJ2GravityModel(), Push(1e-3)), directory; resume=true))
    plain = (InverseSquaredJ2GravityModel(), Push(1e-3))
    plain_directory = interrupted(plain)
    @test !hasproperty(SE._load_checkpoint(checkpointed(plain, plain_directory)).runtime_state, :effector_scopes)
    @test_throws "records no effector scopes" run!(checkpointed(entries, plain_directory; resume=true))
end

end # module EffectorScopeTests
