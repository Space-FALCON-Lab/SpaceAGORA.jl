module GramDensityPerturbationTests
using Test, SpaceAGORA, StaticArrays
const CB = SpaceAGORA.SimulationModel.SimulationCallbacks

@testset "Optional perturbation preserves ordinary atmosphere contexts" begin
    wind = SVector(12.0, -3.0, 0.0)
    for rho in (0.0, -0.0, 1.23456789e-9), slot in (false, true)
        buffers = slot ? (gram_density_perturbation=Ref(nothing),) : (;)
        p = (shared_buffers=buffers,)
        @test isequal(CB._apply_gram_density_perturbation(p, 1, 3.0, 100_000.0, rho, 210.0, wind),
                      (rho, 210.0, wind))
    end
end

@testset "Held and pass factors preserve temperature, wind and per-spacecraft state" begin
    wind = SVector(12.0, -3.0, 0.0)
    for mode in (:step, :pass)
        st = CB._new_gram_density_perturbation_state(mode, 2, 120_000.0, 1.0, 10.0, "")
        p = (shared_buffers=(gram_density_perturbation=Ref(st),),)
        st.held_r .= [2.0, 3.0]
        st.pass_active .= true
        st.pass_t0 .= 5.0
        st.pass_r[1] = [2.0, 4.0, 6.0]
        st.pass_r[2] = [3.0, 5.0, 7.0]
        before = deepcopy(st)
        for sat in 1:2
            factor = mode === :step ? sat + 1.0 : sat + 2.0
            @test CB._apply_gram_density_perturbation(p, sat, 5.5, 120_000.0, 0.5, 210.0, wind) ==
                  (0.5 * factor, 210.0, wind)
            @test CB._apply_gram_density_perturbation(p, sat, 5.5, 120_001.0, 0.5, 210.0, wind) ==
                  (0.5, 210.0, wind)
        end
        for field in fieldnames(typeof(st))
            @test isequal(getfield(st, field), getfield(before, field))
        end
    end
end
end
