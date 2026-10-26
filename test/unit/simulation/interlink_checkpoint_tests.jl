module InterlinkCheckpointTests

using Test
using Serialization
using SpaceAGORA

const SE = SpaceAGORA.SimulationEngine

# A second link kind, for a checkpoint whose link kind changes.
struct OtherLink <: AbstractInterlink
    range_m::Float64
end
SpaceAGORA.link_range_m(link::OtherLink) = link.range_m

# Kinds whose type parameter can select their physics: same fields, different
# hooks for different parameters.
struct ModeLink{M} <: AbstractInterlink
    range_m::Float64
end
SpaceAGORA.link_range_m(link::ModeLink) = link.range_m
struct TypedLink{T} <: AbstractInterlink
    range_m::Float64
end
SpaceAGORA.link_range_m(link::TypedLink) = link.range_m

# Same-named types in two modules of one package (here the tests): a link kind
# in each, and a type in each to parameterize a link with.
module Left
import SpaceAGORA
struct Mode end
struct SameLink <: SpaceAGORA.AbstractInterlink
    range_m::Float64
end
SpaceAGORA.link_range_m(link::SameLink) = link.range_m
end
module Right
import SpaceAGORA
struct Mode end
struct SameLink <: SpaceAGORA.AbstractInterlink
    range_m::Float64
end
SpaceAGORA.link_range_m(link::SameLink) = link.range_m
end

# A kind with the laser link's name and fields, from another package (here the
# tests): its checkpoint identity differs from the built-in kind's.
module ForeignKinds
import SpaceAGORA
struct LaserMomentumLink <: SpaceAGORA.AbstractInterlink
    power_w::Float64
    magnification::Float64
    range_m::Float64
end
SpaceAGORA.link_range_m(link::LaserMomentumLink) = link.range_m
end

function checkpoint_fixture(directory)
    spacecraft = [(n_terminal=2,) for _ in 1:3]
    model = InterlinkModel(forbidden_pairs=[(3, 2)])
    register_link!(model, spacecraft, (1, 1), (2, 1),
        LaserMomentumLink(power_w=12.0, magnification=3.0, range_m=100.0))
    policy = SchedulingPolicyModel(SemimajorAxisRatePolicy(targets=[1, 2]))
    args = (interlink_model=model, scheduling_policy_model=policy,
        dynamics_model=(spacecraft=spacecraft,),
        simulation_settings=(checkpoint_directory=directory, results_directory=directory))
    state = (sc=[(pos=[Float64(i), 0.0, 0.0], interlink_dv=zeros(3), interlink_delta_sma=0.0) for i in 1:3],)
    return args, state
end

@testset "Interlink checkpoint identity and activity" begin
    mktempdir() do directory
        args, state = checkpoint_fixture(directory)
        flags = Bool[true, false, true]
        snapshot = SE._interlink_checkpoint_runtime_state(args, flags)
        flags[2] = true
        @test snapshot.interlinks.is_active == [true, false, true]
        @test snapshot.interlinks.version === 2
        @test snapshot.interlinks.configuration.candidates == [(endpoints=((1, 1), (2, 1)),
            link=(kind=(module_path=(:SpaceAGORA, :SimulationModel, :InterLinkModels), name=:LaserMomentumLink,
                parameters=()),
                fields=(power_w=12.0, magnification=3.0, range_m=100.0)))]
        foreign = SE._interlink_link_identity(ForeignKinds.LaserMomentumLink(12.0, 3.0, 100.0))
        @test foreign.kind.name === :LaserMomentumLink &&
            foreign.kind.module_path == (fullname(@__MODULE__)..., :ForeignKinds)
        @test foreign.fields == snapshot.interlinks.configuration.candidates[1].link.fields
        @test snapshot.interlinks.configuration.forbidden_pairs == [(2, 3)]
        @test snapshot.interlinks.configuration.terminal_counts == [2, 2, 2]
        SE._write_checkpoint!(args, 12.5, state, "tsit5"; is_active=snapshot.interlinks.is_active)
        checkpoint = SE._load_checkpoint(args)
        @test checkpoint.t == 12.5
        # The state keeps the interlink states' names before 0.6.0, with its values.
        @test propertynames(checkpoint.u.sc[1]) == (:pos, :laser_dv, :laser_delta_sma)
        @test checkpoint.u == (sc=[(pos=sc.pos, laser_dv=sc.interlink_dv, laser_delta_sma=sc.interlink_delta_sma)
            for sc in state.sc],)
        @test SE._resumed_checkpoint_state(args, checkpoint.u) == state
        # A state whose spacecraft name the interlink states differently is refused. A
        # run's own state cannot be built that way; a state given as named tuples can.
        uneven = (sc=[i == 1 ? sc : (pos=sc.pos, laser_dv=sc.interlink_dv, laser_delta_sma=sc.interlink_delta_sma)
            for (i, sc) in enumerate(state.sc)],)
        @test SE._interlink_state_naming(uneven) === :mixed
        @test_throws ArgumentError SE._check_interlink_checkpoint_state(args, (u=uneven,))
        # Every checkpoint also records the flags in its activity entry.
        @test checkpoint.runtime_state == merge(snapshot, (activity=(version=1, is_active=[true, false, true]),))
        restored = trues(3)
        @test isnothing(SE._restore_interlink_checkpoint_state!(args, checkpoint, restored))
        @test restored == [true, false, true]
        restored = trues(3)
        @test isnothing(SE._restore_activity_checkpoint_state!(args, checkpoint, restored))
        @test restored == [true, false, true]

        changes = (
            "topology" => a -> begin
                register_link!(a.interlink_model, a.dynamics_model.spacecraft, (1, 2), (3, 1), LaserMomentumLink())
                a
            end,
            "power" => a -> begin
                a.interlink_model.graph.edge_properties[1] = LaserMomentumLink(power_w=24.0, magnification=3.0, range_m=100.0)
                a
            end,
            "magnification" => a -> begin
                a.interlink_model.graph.edge_properties[1] = LaserMomentumLink(power_w=12.0, magnification=4.0, range_m=100.0)
                a
            end,
            "range" => a -> begin
                a.interlink_model.graph.edge_properties[1] = LaserMomentumLink(power_w=12.0, magnification=3.0, range_m=101.0)
                a
            end,
            "forbidden pairs" => a -> begin
                empty!(a.interlink_model.forbidden_pairs)
                a
            end,
            "link kind" => a -> begin
                a.interlink_model.graph.edge_properties[1] = OtherLink(100.0)
                a
            end,
            "same-named kind of another package" => a -> begin
                a.interlink_model.graph.edge_properties[1] = ForeignKinds.LaserMomentumLink(12.0, 3.0, 100.0)
                a
            end,
            "policy" => a -> merge(a, (scheduling_policy_model=SchedulingPolicyModel(SemimajorAxisRatePolicy(targets=:lower_index)),)),
            "targets" => a -> merge(a, (scheduling_policy_model=SchedulingPolicyModel(SemimajorAxisRatePolicy(targets=[1, 3])),)),
            "terminal counts" => a -> begin
                a.dynamics_model.spacecraft[3] = (n_terminal=3,)
                a
            end,
        )
        for (name, change) in changes
            @testset "$name" begin
                changed, _ = checkpoint_fixture(directory)
                destination = trues(3)
                @test_throws ArgumentError SE._restore_interlink_checkpoint_state!(change(changed), checkpoint, destination)
                @test all(destination) # Failed admission must not partly restore state.
            end
        end

        @test_throws ArgumentError SE._write_checkpoint!(args, 12.5, state, "tsit5")
        for invalid_runtime in (nothing,
            (interlinks=merge(snapshot.interlinks, (version=3,)),),
            (interlinks=merge(snapshot.interlinks, (is_active=Bool[true, false],)),),
            (interlinks=merge(snapshot.interlinks, (is_active=[1, 0, 1],)),))
            invalid = merge(checkpoint, (runtime_state=invalid_runtime,))
            @test_throws ArgumentError SE._restore_interlink_checkpoint_state!(args, invalid, trues(3))
        end
        # A version 1 snapshot (SpaceAGORA 0.4 and earlier: one laser link type
        # for the run, P, B and range per candidate) resumes the same run; with
        # another link type it is refused.
        v1_configuration = (link_type=:laser, vertices=[(1, 1), (2, 1)],
            candidates=[(endpoints=((1, 1), (2, 1)), P=12.0, B=3.0, range=100.0)],
            forbidden_pairs=[(2, 3)], policy=:gve_sma, target_idx=1, target_idxs=[1, 2], terminal_counts=[2, 2, 2])
        v1 = merge(checkpoint, (runtime_state=(interlinks=(version=1, configuration=v1_configuration,
            is_active=[true, false, true]),),))
        restored = trues(3)
        @test isnothing(SE._restore_interlink_checkpoint_state!(args, v1, restored))
        @test restored == [true, false, true]
        radio = merge(checkpoint, (runtime_state=(interlinks=(version=1,
            configuration=merge(v1_configuration, (link_type=:radio,)), is_active=[true, false, true]),),))
        @test_throws ArgumentError SE._restore_interlink_checkpoint_state!(args, radio, trues(3))

        # A kind's module and its concrete type parameters are part of its identity:
        # a checkpoint, read back from its file, resumes the same kind and refuses
        # another. The cases are value, type, nested type, Union and tuple
        # parameters; same-named link kinds of two modules of one package, and such
        # types in a parameter; value parameters equal as numbers but of
        # different types or bits, alone and in a tuple; and composite bits values
        # whose own isequal merges them (two NaN payloads in a ComplexF64).
        nan_payload = reinterpret(Float64, 0x7ff8000000000001)
        for (written, same, changed) in (
                (TypedLink{1}(100.0), TypedLink{1}(100.0), TypedLink{1.0}(100.0)),
                (TypedLink{1}(100.0), TypedLink{1}(100.0), TypedLink{true}(100.0)),
                (TypedLink{Int64(1)}(100.0), TypedLink{Int64(1)}(100.0), TypedLink{Int32(1)}(100.0)),
                (TypedLink{(1, :a)}(100.0), TypedLink{(1, :a)}(100.0), TypedLink{(1.0, :a)}(100.0)),
                (TypedLink{(1, :a)}(100.0), TypedLink{(1, :a)}(100.0), TypedLink{(Int32(1), :a)}(100.0)),
                (TypedLink{0.0}(100.0), TypedLink{0.0}(100.0), TypedLink{-0.0}(100.0)),
                (TypedLink{NaN}(100.0), TypedLink{NaN}(100.0), TypedLink{nan_payload}(100.0)),
                (TypedLink{complex(NaN, 0.0)}(100.0), TypedLink{complex(NaN, 0.0)}(100.0),
                    TypedLink{complex(nan_payload, 0.0)}(100.0)),
                (TypedLink{complex(1.0, 2.0)}(100.0), TypedLink{complex(1.0, 2.0)}(100.0),
                    TypedLink{complex(1.0, -2.0)}(100.0)),
                (TypedLink{(complex(NaN, 0.0), :a)}(100.0), TypedLink{(complex(NaN, 0.0), :a)}(100.0),
                    TypedLink{(complex(nan_payload, 0.0), :a)}(100.0)),
                (Left.SameLink(100.0), Left.SameLink(100.0), Right.SameLink(100.0)),
                (TypedLink{Left.Mode}(100.0), TypedLink{Left.Mode}(100.0), TypedLink{Right.Mode}(100.0)),
                (TypedLink{Vector{Left.Mode}}(100.0), TypedLink{Vector{Left.Mode}}(100.0), TypedLink{Vector{Right.Mode}}(100.0)),
                (ModeLink{1}(100.0), ModeLink{1}(100.0), ModeLink{-1}(100.0)),
                (TypedLink{Float64}(100.0), TypedLink{Float64}(100.0), TypedLink{Float32}(100.0)),
                (TypedLink{Vector{Int}}(100.0), TypedLink{Vector{Int}}(100.0), TypedLink{Vector{Int32}}(100.0)),
                (TypedLink{Union{Int, Nothing}}(100.0), TypedLink{Union{Int, Nothing}}(100.0), TypedLink{Union{Int, Missing}}(100.0)),
                (TypedLink{(1, :a)}(100.0), TypedLink{(1, :a)}(100.0), TypedLink{(1, :b)}(100.0)))
            written_args, written_state = checkpoint_fixture(directory)
            written_args.interlink_model.graph.edge_properties[1] = written
            SE._write_checkpoint!(written_args, 12.5, written_state, "tsit5"; is_active=[true, false, true])
            from_file = SE._load_checkpoint(written_args)
            positive, _ = checkpoint_fixture(directory)
            positive.interlink_model.graph.edge_properties[1] = same
            restored = trues(3)
            @test isnothing(SE._restore_interlink_checkpoint_state!(positive, from_file, restored))
            @test restored == [true, false, true]
            negative, _ = checkpoint_fixture(directory)
            negative.interlink_model.graph.edge_properties[1] = changed
            untouched = trues(3)
            @test_throws ArgumentError SE._restore_interlink_checkpoint_state!(negative, from_file, untouched)
            @test all(untouched)
        end
        int64 = (module_path=(:Core,), name=:Int64, parameters=())
        @test SE._interlink_link_identity(ModeLink{-1}(1.0)).kind ==
            (module_path=fullname(@__MODULE__), name=:ModeLink, parameters=((type=int64, value=-1),))
        @test SE._interlink_link_identity(TypedLink{Vector{Int}}(1.0)).kind.parameters ==
            ((module_path=(:Core,), name=:Array, parameters=(int64, (type=int64, value=1))),)
        float64 = (module_path=(:Core,), name=:Float64, parameters=())
        @test SE._interlink_link_identity(TypedLink{-0.0}(1.0)).kind.parameters ==
            ((type=float64, value=-0.0, bits=0x8000000000000000),)
        @test isequal(complex(NaN, 0.0), complex(nan_payload, 0.0))          # why fields, not the value
        @test SE._interlink_link_identity(TypedLink{complex(1.0, -0.0)}(1.0)).kind.parameters ==
            ((type=(module_path=(:Base,), name=:Complex, parameters=(float64,)),
              fields=((type=float64, value=1.0, bits=0x3ff0000000000000),
                      (type=float64, value=-0.0, bits=0x8000000000000000))),)
        @test SE._interlink_link_identity(TypedLink{Left.Mode}(1.0)).kind.parameters ==
            ((module_path=(fullname(@__MODULE__)..., :Left), name=:Mode, parameters=()),)
        SE._write_checkpoint!(args, 12.5, state, "tsit5"; is_active=snapshot.interlinks.is_active)

        legacy_args = merge(args, (interlink_model=nothing,))
        @test_throws ArgumentError SE._restore_interlink_checkpoint_state!(legacy_args, checkpoint, trues(3))

        # An actual old payload has no runtime_state member. The generic reader
        # still loads it, and ordinary non-interlink restoration is unchanged.
        paths = SE._checkpoint_paths(args)
        legacy_state = (sc=[(pos=copy(spacecraft.pos),) for spacecraft in state.sc],)
        open(paths.data, "w") do io
            serialize(io, (schema_version="1", t=12.5, u=legacy_state, solver_mode="tsit5"))
        end
        legacy = SE._load_checkpoint(legacy_args)
        @test get(legacy, :runtime_state, nothing) === nothing
        @test !haskey(legacy, :runtime_state)
        @test legacy.u == legacy_state
        @test_throws ArgumentError SE._restore_interlink_checkpoint_state!(args, legacy, trues(3))
        destination = Bool[true, false, true]
        @test isnothing(SE._restore_interlink_checkpoint_state!(legacy_args, legacy, destination))
        @test destination == [true, false, true]
        # An explicitly empty runtime member also retains the ordinary reader shape.
        open(paths.data, "w") do io
            serialize(io, (schema_version="1", t=12.5, u=legacy_state,
                solver_mode="tsit5", runtime_state=nothing))
        end
        @test SE._load_checkpoint(legacy_args) == legacy
        SE._write_checkpoint!(legacy_args, 13.0, legacy_state, "tsit5")
        ordinary = SE._load_checkpoint(legacy_args)
        @test get(ordinary, :runtime_state, nothing) === nothing
        @test keys(ordinary) == keys(legacy)
        @test keys(deserialize(paths.data)) == (:schema_version, :created_utc, :t, :solver_mode, :u)
    end
end

end # module InterlinkCheckpointTests
