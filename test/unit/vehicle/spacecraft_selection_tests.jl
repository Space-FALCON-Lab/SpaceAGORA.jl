module SpacecraftSelectionTests

# SpacecraftSelection and resolve_spacecraft_scope: the spacecraft of a run chosen
# by position, id, constellation plane or slot, or a predicate.

using Test
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const resolve = SM.resolve_spacecraft_scope

function bus_factory(plane, slot, index)
    bus = Link(root=true, m=100.0 + index)
    return SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia)
end

# Walker 6/3/1: planes 1-3 with slots 1-2; positions 1-6 and ids 1-6 in builder order.
walker() = walker_constellation(bus_factory; total=6, planes=3, phasing=1, altitude=550e3,
    planet=make_no_gram_planet(:earth), i=53)

plain(n) = [(bus = Link(root=true, m=100.0); SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia, id=i))
            for i in 1:n]

@testset "SpacecraftSelection construction" begin
    @test SpacecraftSelection(planes=1).planes == [1]
    @test SpacecraftSelection(positions=(1, 3)).positions == [1, 3]
    @test SpacecraftSelection(ids=101:103).ids == [101, 102, 103]
    @test SpacecraftSelection(slots=(2, 1)).slots == [ConstellationSlot(2, 1)]
    @test SpacecraftSelection(slots=[(1, 1), ConstellationSlot(2, 2)]).slots == [ConstellationSlot(1, 1), ConstellationSlot(2, 2)]
    @test SpacecraftSelection(planes=1) == SpacecraftSelection(planes=[1])
    @test SpacecraftSelection(planes=1) != SpacecraftSelection(planes=2)
    @test SpacecraftSelection(planes=1) isa SpacecraftScope
    @test_throws "at least one" SpacecraftSelection()
    @test_throws "numbered from 1" SpacecraftSelection(positions=0)
    @test_throws "numbered from 1" SpacecraftSelection(planes=[1, 0])
    @test_throws "integers" SpacecraftSelection(ids=["a"])
    @test_throws ArgumentError SpacecraftSelection(slots=(0, 1))
    @test_throws "ConstellationSlot" SpacecraftSelection(slots=[3])
    @test_throws "function" SpacecraftSelection(where=3)
end

@testset "resolution against a Walker layout" begin
    sc = walker()
    @test [s.constellation_slot for s in sc] == [ConstellationSlot(p, s) for p in 1:3 for s in 1:2]
    @test resolve(AllSpacecraft(), sc) == 1:6
    @test resolve(BoundSpacecraft((5, 1, 5)), sc) == [1, 5]
    @test resolve(SpacecraftSelection(positions=[4, 2]), sc) == [2, 4]
    @test resolve(SpacecraftSelection(ids=[2, 4]), sc) == [2, 4]
    @test resolve(SpacecraftSelection(planes=2), sc) == [3, 4]
    @test resolve(SpacecraftSelection(planes=[1, 3]), sc) == [1, 2, 5, 6]
    @test resolve(SpacecraftSelection(slots=[(3, 2)]), sc) == [6]
    @test resolve(SpacecraftSelection(slots=ConstellationSlot(1, 1)), sc) == [1]
    @test resolve(SpacecraftSelection(where=(i, s) -> isodd(i)), sc) == [1, 3, 5]
    @test resolve(SpacecraftSelection(where=(i, s) -> s.dry_mass > 103.5), sc) == [4, 5, 6]
    # Every criterion given must hold.
    @test resolve(SpacecraftSelection(planes=1:2, where=(i, s) -> isodd(i)), sc) == [1, 3]
    @test resolve(SpacecraftSelection(planes=1, ids=[3]), sc) == Int[]
end

@testset "resolution refusals" begin
    sc = walker()
    @test_throws "outside the run's positions 1 to 6" resolve(BoundSpacecraft((7,)), sc; owner="probe")
    @test_throws "outside the run's positions 1 to 6" resolve(SpacecraftSelection(positions=[2, 9]), sc)
    @test_throws "id 99, which no spacecraft" resolve(SpacecraftSelection(ids=99), sc)
    shared = plain(3)
    shared[3].id = 1
    @test_throws "positions 1, 3 share" resolve(SpacecraftSelection(ids=1), shared)
    @test resolve(SpacecraftSelection(ids=2), shared) == [2]
    @test_throws "no spacecraft of the run has a constellation_slot" resolve(SpacecraftSelection(planes=1), plain(2))
    @test_throws "plane 4, which no spacecraft" resolve(SpacecraftSelection(planes=4), sc)
    @test_throws "(plane 1, slot 3), which no" resolve(SpacecraftSelection(slots=(1, 3)), sc)
    @test_throws "threw for the spacecraft at run position 1" resolve(SpacecraftSelection(where=(i, s) -> error("no")), sc)
    @test_throws "returned 1" resolve(SpacecraftSelection(where=(i, s) -> 1), sc)
    @test_throws "CoupledSpacecraft" resolve(CoupledSpacecraft(1, (2,)), sc)
    # The owner names the model in each message.
    @test_throws "Dynamic effector 2 (Drag) names spacecraft 8" resolve(BoundSpacecraft((8,)), sc;
        owner="Dynamic effector 2 (Drag)")
end

end # module SpacecraftSelectionTests
