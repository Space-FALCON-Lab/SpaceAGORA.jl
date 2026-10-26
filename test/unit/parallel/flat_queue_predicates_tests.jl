using Test
using SpaceAGORA

const FQ_SM = SpaceAGORA.SimulationModel
const FQ_SE = SpaceAGORA.SimulationEngine
const FQ_ER = FQ_SM.EffectorRouting

@testset "Flat-queue predicates: the engine routes with the shared owner" begin
    # The engine's pre-pass predicates are the EffectorRouting functions, not
    # copies, including the value-dependent gravity-gradient and per-link
    # exclusions.
    @test FQ_SE._batchable_effector === FQ_ER.batchable_effector
    @test FQ_SE._harmonics_prepass_effector === FQ_ER.harmonics_prepass_effector
    @test FQ_SE._aero_prepass_effector === FQ_ER.aero_prepass_effector

    E = FQ_SM.Earth()
    harm_file = joinpath(@__DIR__, "..", "..", "..", "data",
                         "planets", "gravity", "EarthGGM05C.csv")
    nbody = FQ_SM.NBodyGravityModel(body_names=("Sun",), planet=E, body_mus=(1.32712440041279419e20,))

    # (effector, batchable, harmonics pre-pass, aero pre-pass). gravity_gradient
    # takes the central terms off the batch, and per-link atmosphere keeps the fM
    # model on the queue, so both variants are pinned: a type-only predicate
    # would pass the default case and get these wrong.
    cases = Any[
        (FQ_SM.InverseSquaredGravityModel(), true, false, false),
        (FQ_SM.InverseSquaredJ2GravityModel(), true, false, false),
        (nbody, true, false, false),
        (FQ_SM.SolarRadiationPressureModel(1.8, 10.0), true, false, false),
        (FQ_SM.InverseSquaredGravityModel(gravity_gradient=true), false, false, false),
        (FQ_SM.InverseSquaredJ2GravityModel(gravity_gradient=true), false, false, false),
        (FQ_SM.AerodynamicCoefficientfM(), false, false, true),
        (FQ_SM.AerodynamicCoefficientfM(per_link_atmosphere=true), false, false, false),
        # Unknown and custom effectors stay on the queue.
        (nothing, false, false, false),
    ]
    if isfile(harm_file)
        push!(cases, (FQ_SM.GravitationalHarmonicsModel(4, 4, harm_file, E), false, true, false))
    end
    for (e, batch, harmonics, aero) in cases
        @test FQ_ER.batchable_effector(e) == batch
        @test FQ_ER.harmonics_prepass_effector(e) == harmonics
        @test FQ_ER.aero_prepass_effector(e) == aero
    end

    # The engine's queue-only count: the two gravity-gradient variants, per-link
    # aero and the unknown effector.
    effectors = Tuple(first.(cases))
    @test FQ_SE._count_flat_queue_only_effectors(effectors) == 4
    @test count(FQ_SE._flat_selection_mask(effectors, nothing)) == 4

    # A (gravity, aero) stack: both are pre-passed, so the queue builds no nodes.
    aero_stack = (FQ_SM.InverseSquaredJ2GravityModel(), FQ_SM.AerodynamicCoefficientfM())
    @test FQ_SE._count_flat_queue_only_effectors(aero_stack) == 0
    per_link_stack = (FQ_SM.InverseSquaredJ2GravityModel(),
                      FQ_SM.AerodynamicCoefficientfM(per_link_atmosphere=true))
    @test FQ_SE._count_flat_queue_only_effectors(per_link_stack) == 1
    @test FQ_SE._count_flat_queue_only_effectors((nothing, aero_stack...)) == 1

    # Under a solver partition the pre-passes do not apply: every effector of
    # the partition goes through the queue, the pre-passed ones included.
    for (e, _, _, _) in cases
        e === nothing && continue
        partition = FQ_SE._solver_partition_validated(e)
        @test FQ_SE._flat_selection_mask((e,), partition) == (true,)
    end
end
