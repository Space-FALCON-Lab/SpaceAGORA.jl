module GramDensityPerturbationTests
using Test, SpaceAGORA, StaticArrays
const CB = SpaceAGORA.SimulationModel.SimulationCallbacks
const EM = SpaceAGORA.SimulationModel.EnvironmentModels

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

# The distributed RHS density service, with only the pool evaluator replaced: it
# serves a mean density of 2 at every query. Eligibility, the fill and the buffer
# writes are the production functions. No worker or native library is involved.
struct ServiceProbeParams
    args::Any
    shared_buffers::Any
    is_active::Vector{Bool}
end

function CB._gram_process_pool_batch_eval!(
    rhos::AbstractVector{Float64}, Ts::AbstractVector{Float64},
    winds::AbstractVector{SVector{3, Float64}}, model,
    hs::AbstractVector{<:Real}, lats::AbstractVector{<:Real}, lons::AbstractVector{<:Real},
    t::Union{Float64, AbstractVector{<:Real}}, wind::Bool, p::ServiceProbeParams,
)::Bool
    fill!(rhos, 2.0)
    fill!(Ts, 200.0)
    fill!(winds, SVector(0.0, 0.0, 0.0))
    return true
end

function service_probe(st; n::Int=2, slot::Bool=true)
    model = EM.GRAMAtmosphereModel(nothing, ReentrantLock(), Dict{Symbol, Any}(:seed => 11))
    args = (environment_model=(planet=(T_ref=200.0,), EI=120.0, density_model=model, wind=false),
            mission_configuration=(keplerian=true,))
    zero = SVector(0.0, 0.0, 0.0)
    sb = (densities=zeros(n), temperatures=zeros(n), winds=fill(zero, n), density_sample_t=fill(NaN, n),
          density_models=EM.GRAMAtmosphereModel[])
    slot && (sb = merge(sb, (gram_density_perturbation=Ref{Any}(st),)))
    return ServiceProbeParams(args, sb, fill(true, n))
end

const SERVICE_ENV = ("SPACEAGORA_GRAM_PROCESS_POOL" => "on", "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing,
                     "SPACEAGORA_GRAM_TRACK_CACHE" => "off", "SPACEAGORA_VACUUM_GRAM_CACHE" => "off",
                     "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "off")

@testset "Distributed RHS density service applies the perturbation factor" begin
    wind = SVector(0.0, 0.0, 0.0)
    withenv(SERVICE_ENV...) do
        for mode in (:step, :pass)
            st = CB._new_gram_density_perturbation_state(mode, 2, 120_000.0, 1.0, 10.0, "")
            st.held_r .= [5.0, 3.0]
            st.pass_active .= true
            st.pass_t0 .= 0.0
            st.pass_r[1] = [5.0, 5.0, 5.0]
            st.pass_r[2] = [3.0, 3.0, 3.0]
            p = service_probe(st)
            @test CB._rhs_density_service_candidate(p, 2)
            # Spacecraft 1 at the entry interface itself, 2 just above it (the
            # keplerian run still asks GRAM there, and the factor is 1).
            alts = [120_000.0, 120_001.0]
            @test CB._rhs_density_service_fill!(p, 1.0, 2, alts, [0.0, 0.0], [0.0, 0.0])
            for i in 1:2
                scalar = CB._apply_gram_density_perturbation(p, i, 1.0, alts[i], 2.0, 200.0, wind)
                @test p.shared_buffers.densities[i] == scalar[1]
                @test p.shared_buffers.temperatures[i] == scalar[2]
                @test p.shared_buffers.density_sample_t[i] == 1.0
            end
            @test p.shared_buffers.densities == [10.0, 2.0]
            @test CB._rhs_density_service_fill!(p, 1.0, 2, [100_000.0, 100_000.0], [0.0, 0.0], [0.0, 0.0])
            @test p.shared_buffers.densities == [10.0, 6.0]
        end

        # The naive control reads the coordinator's own instance, which a served
        # query never updates: the service declines it.
        naive = CB._new_gram_density_perturbation_state(:naive_rhs, 2, 120_000.0, 1.0, 10.0, "")
        @test !CB._rhs_density_service_candidate(service_probe(naive), 2)

        # No mode installed, or no slot at all: the served mean, unchanged.
        for p in (service_probe(nothing), service_probe(nothing; slot=false))
            @test CB._rhs_density_service_candidate(p, 2)
            @test CB._rhs_density_service_fill!(p, 1.0, 2, [100_000.0, 100_000.0], [0.0, 0.0], [0.0, 0.0])
            @test p.shared_buffers.densities == [2.0, 2.0]
        end
    end
end
end
