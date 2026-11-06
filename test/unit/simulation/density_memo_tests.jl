module DensityMemoTests
using Test, SpaceAGORA, StaticArrays
const SM = SpaceAGORA.SimulationModel
const CB = SM.SimulationCallbacks

@testset "density query memo" begin
    memos = Union{Nothing, SM.DensityQueryMemo}[nothing, nothing]
    model = Ref(1)   # identity key only
    calls = Ref(0)
    query(alt; sat=1, m=model, t=10.0, wind=false) = CB._density_memoized(memos, sat, m, alt, 0.1, 0.2, t, wind) do
        calls[] += 1
        (alt * 1e-12, 200.0 + calls[], SVector(1.0, 2.0, 3.0))
    end

    r1 = query(1.0e5)
    @test calls[] == 1
    @test query(1.0e5) === r1 && calls[] == 1                 # identical inputs: served from the slot
    @test query(nextfloat(1.0e5)) != r1 && calls[] == 2       # one ulp in altitude: miss
    @test query(nextfloat(1.0e5)) == query(nextfloat(1.0e5)) && calls[] == 2
    query(1.0e5; t=nextfloat(10.0)); @test calls[] == 3       # time differs: miss
    query(1.0e5; t=nextfloat(10.0), wind=true); @test calls[] == 4   # wind flag differs: miss
    query(1.0e5; t=nextfloat(10.0), wind=true, m=Ref(1)); @test calls[] == 5   # other model: miss
    query(1.0e5; sat=2); @test calls[] == 6                   # slots are per satellite
    query(1.0e5; sat=2); @test calls[] == 6
    query(2.0e5; sat=2); query(1.0e5; sat=1, t=nextfloat(10.0), wind=true, m=memos[1].model)
    @test calls[] == 7   # sat 2 missed once; sat 1 hit its stored record
end

@testset "density memo eligibility" begin
    planet = SM.make_no_gram_planet(:earth)
    p = (; shared_buffers=(; density_memos=Union{Nothing, SM.DensityQueryMemo}[nothing]))
    # Only a deterministic native GRAM mean is ever memoized.
    @test !CB._density_memo_eligible(p, SM.ExponentialAtmosphereModel(planet))
    @test !CB._density_memo_eligible(p, SM.NoAtmosphereModel())
end
end
