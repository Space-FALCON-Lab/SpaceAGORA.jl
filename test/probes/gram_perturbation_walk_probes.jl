# Native probe: requires GRAM (Mars). Standalone:
#   julia --project=. test/probes/gram_perturbation_walk_probes.jl
#
# Checks the instance discipline the SPACEAGORA_GRAM_DENSITY_PERTURBATION modes
# rely on (src/simulation/callbacks/density_callbacks/gram_density_perturbation.jl):
#
#  1. A walk instance built with _gram_walk_clone reproduces the same r sequence
#     whether or not the MEAN instance is queried between its samples -- i.e.
#     right-hand-side calls on the mean instance cannot move the walk.
#  2. Interleaving those same extra queries on the walk instance ITSELF changes
#     the sequence: sharing one instance between the RHS and the walk would
#     corrupt it, which is why there are two.
#  3. A fresh clone replays the sequence bit for bit (seed carried by the recipe).
#  4. A different seed gives a different sequence; the mean density does not move.
#  5. With scales 0 the walk is silent (r == 1 exactly).
module GRAMPerturbationWalkProbes
using Test, SpaceAGORA
if Base.find_package("GRAMSuite") === nothing
    pushfirst!(LOAD_PATH, joinpath(dirname(dirname(pathof(SpaceAGORA))), "data", "GRAMSuite.jl"))
end
using GRAMSuite

const SM = SpaceAGORA.SimulationModel
const EM = SM.EnvironmentModels

mars_model(seed, scale) = EM.GRAMAtmosphereModel(
    planet_name="mars",
    initial_time=SM.InitialTime(year=2001, month=11, day=6, hour=19, minute=0, second=32.0),
    seed=seed,
    gram_perturbation_scales=(scale, scale, scale, scale),
    mars_map_year=2,
)

# A pass-like path: 0.5 s spacing, 130 -> 100 -> 130 km, ~4.5 km/s ground speed.
const N = 240
path_h(k) = 1e3 * (100.0 + 30.0 * abs(1.0 - 2.0 * (k - 1) / (N - 1)))
path_lat(k) = deg2rad(60.0 + 0.017 * (k - 1))
path_lon(k) = deg2rad(30.0 + 0.020 * (k - 1))
path_t(k) = 0.5 * (k - 1)
# "RHS stage" points near, but not at, the path point (as integrator stages are).
stage(k, j) = (path_h(k) + 50.0 * j, path_lat(k) + 1e-5 * j, path_lon(k) - 1e-5 * j, path_t(k) + 0.05 * j)

walk_seq(w; interleave=nothing) = map(1:N) do k
    if interleave !== nothing
        for j in 1:5
            EM._gram_point_density(interleave, stage(k, j)..., false)
        end
    end
    EM._gram_walk_sample(w, path_h(k), path_lat(k), path_lon(k), path_t(k))
end

@testset "GRAM perturbation walk instance isolation" begin
    @test Base.get_extension(SpaceAGORA, :SpaceAGORAGRAMSuiteExt) !== nothing
    withenv("SPACEAGORA_GRAM_STATIC_GRID" => "0") do
        mean_model = mars_model(11, 1.0)
        # Furnish the Mars kernels the GRAM ephemeris queries need, as a solve would.
        SpaceAGORA.TelemetryVerification._planet_from_name("mars")

        ref = walk_seq(EM._gram_walk_clone(mean_model))
        r_ref = [s[1] / s[2] for s in ref]
        @test all(isfinite, r_ref)
        @test any(r -> r != 1.0, r_ref)                    # the walk is live at scale 1

        # 1. RHS-style queries on the MEAN instance do not move the walk.
        isolated = walk_seq(EM._gram_walk_clone(mean_model); interleave=mean_model)
        @test isolated == ref

        # 2. The same queries on the walk instance itself corrupt it.
        shared_walk = EM._gram_walk_clone(mean_model)
        shared = walk_seq(shared_walk; interleave=shared_walk)
        @test shared != ref
        @test [s[2] for s in shared] == [s[2] for s in ref]  # mean unaffected, only the walk moved

        # 3. A fresh clone replays it bit for bit.
        @test walk_seq(EM._gram_walk_clone(mean_model)) == ref

        # 4. Another seed: different walk, identical mean.
        other = walk_seq(EM._gram_walk_clone(mars_model(22, 1.0)))
        @test [s[1] for s in other] != [s[1] for s in ref]
        @test [s[2] for s in other] == [s[2] for s in ref]

        # Mean density from the walk instance equals what the RHS path returns.
        rhs_mean = [EM._gram_point_density(mean_model, path_h(k), path_lat(k), path_lon(k), path_t(k), false)[1] for k in 1:N]
        @test rhs_mean == [s[2] for s in ref]

        # 5. Scales 0: the walk is silent.
        silent = walk_seq(EM._gram_walk_clone(mars_model(11, 0.0)))
        @test all(s -> s[1] == s[2], silent)

        r = r_ref
        μ = sum(r) / N
        σ = sqrt(sum((r .- μ) .^ 2) / (N - 1))
        println("probe: r mean=$(μ) std=$(σ) GRAM sigma median=$(sort([s[3] for s in ref])[N ÷ 2])")
    end
end
end
