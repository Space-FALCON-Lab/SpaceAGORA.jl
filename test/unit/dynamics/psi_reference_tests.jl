using Test
using SpaceAGORA

# PR-tier regression tests against the published plume-surface interaction
# reference cases. The numbers are not repeated here: the test loads the same
# frozen manifests the study runs from
# (`benchmarks/studies/psi_validation/manifests/`), so a reference value and its
# source can never drift out of step between the study and the test.
#
# Three of the study's gate-eligible cases are cheap enough for every PR (all
# three are closed-form evaluations of the model's public functions, with no
# simulation, no kernels and no file I/O beyond the manifests):
#
#   pfgt1_erosion_threshold_shear  the shear the model starts eroding at, against
#                                  NASA's subscale vacuum-chamber measurement
#   apollo_lm_surface_shear        the shear it lays on the ground, against the
#                                  Apollo LM CFD fit
#   surveyor3_ejecta_speed         the speed it throws grains at, against the
#                                  Surveyor III sandblasting estimates
#
# The tolerances are the manifests' own, and every one of them is wide and
# justified in the manifest's `tolerance.reason`: these are order-of-magnitude
# reference comparisons, not regression fixtures. A fourth block records the
# case the model currently fails (the Apollo 12 erosion rate) with a band loose
# enough that fixing the model passes it and only a collapse fails it.
module PSIValidationHarness
    using SpaceAGORA
    include(joinpath(@__DIR__, "..", "..", "..", "benchmarks", "studies", "psi_validation", "common.jl"))
end

const PSIH = PSIValidationHarness

_psi_case(id) = PSIH.load_case(joinpath(PSIH.MANIFEST_DIR, id * ".toml"))

"Scored (non-context) reference rows of a case."
_psi_scored(case) = [r for r in case.references if !Bool(get(r, "context_only", false))]

@testset "PSI reference cases" begin
    cfg = PlumeSurfaceConfig()

    @testset "every manifest is well formed and honest about gating" begin
        cases = PSIH.load_cases()
        @test length(cases) >= 6
        @test length(unique(c.id for c in cases)) == length(cases)
        for case in cases
            # A case may only gate if its reference is a published measurement of
            # a quantity the model actually produces.
            @test !case.gate_eligible || case.status == PSIH.GATEABLE_STATUS
            # Every case names its source completely enough to be checked.
            for key in ("authors", "title", "venue", "year", "locator")
                @test haskey(case.source, key)
                @test !isempty(string(case.source[key]))
            end
            @test !isempty(_psi_scored(case))
            if case.gate_eligible
                @test Float64(case.tolerance["low"]) > 0.0
                @test Float64(case.tolerance["high"]) > Float64(case.tolerance["low"])
                @test !isempty(String(case.tolerance["reason"]))
            end
        end
        # At least one case must be on record as unsourced or not modeled: the
        # ground-effect force and the ejection angle are both in that state, and
        # silently dropping them would hide what the model cannot defend.
        @test any(c -> c.status in ("unsourced", "not_modeled"), cases)
    end

    @testset "erosion threshold shear against the NASA subscale vacuum tests" begin
        # Stubbs and Mehta (AIAA SciTech 2026, NTRS 20250011216): the CFD surface
        # shear at the measured crater edge averages 0.25 Pa over mono-disperse
        # sand. The model's own threshold is read back as the peak wall shear at
        # its erosion onset height, so this does not depend on any config field.
        case = _psi_case("pfgt1_erosion_threshold_shear")
        thrust = Float64(case.conditions["thrust_n"])
        modelval = PSIH.model_quantity(cfg, "threshold_shear_pa", thrust, 0.0)
        @test modelval !== nothing
        @test modelval > 0.0
        for ref in _psi_scored(case)
            ratio = modelval / Float64(ref["value"])
            @test Float64(case.tolerance["low"]) <= ratio <= Float64(case.tolerance["high"])
        end
        # The threshold is a property of the soil in this model, not of the
        # engine: reading it back at four times the thrust must give the same
        # number, which is what makes the single-thrust comparison above valid.
        @test PSIH.model_quantity(cfg, "threshold_shear_pa", 4 * thrust, 0.0) ≈ modelval rtol = 1e-9
    end

    @testset "surface shear against the Apollo LM CFD fit" begin
        # Lane and Metzger (Acta Geophysica 63, 568-599, 2015), Eq. (14):
        # tau(h) = 6.21 exp(-0.123 h) Pa, area-averaged over the erosion radius.
        case = _psi_case("apollo_lm_surface_shear")
        thrust = Float64(case.conditions["thrust_n"])
        lo, hi = Float64(case.tolerance["low"]), Float64(case.tolerance["high"])
        for ref in _psi_scored(case)
            h = Float64(ref["height_m"])
            modelval = PSIH.model_quantity(cfg, "peak_shear_pa", thrust, h)
            @test modelval !== nothing
            ratio = modelval / Float64(ref["value"])
            @test lo <= ratio <= hi
        end
        # Both the model and the reference fall off with height; a model that
        # ever rose with height would be wrong for a reason no tolerance covers.
        shears = [PSIH.model_quantity(cfg, "peak_shear_pa", thrust, Float64(ref["height_m"]))
                  for ref in _psi_scored(case)]
        @test issorted(shears; rev=true)
    end

    @testset "ejecta speed against the Surveyor III sandblasting" begin
        # Four published estimates spanning 40 to 2000 m/s, compiled by Immer,
        # Lane, Metzger and Clements (Earth & Space 2008). Passing means only
        # that the model's characteristic speed lands inside that range.
        case = _psi_case("surveyor3_ejecta_speed")
        thrust = Float64(case.conditions["thrust_n"])
        for ref in _psi_scored(case)
            h = Float64(ref["height_m"])
            modelval = PSIH.model_quantity(cfg, "ejecta_speed_mps", thrust, h)
            @test modelval !== nothing
            ratio = modelval / Float64(ref["value"])
            @test Float64(case.tolerance["low"]) <= ratio <= Float64(case.tolerance["high"])
        end
    end

    @testset "the Apollo 12 erosion rate gap is recorded, not enforced" begin
        # This case fails the study's factor-of-three tolerance today: the model
        # returns 2 to 10 times less soil than Lane and Metzger measured from the
        # Apollo 12 descent film, and exactly zero above its own onset height
        # where the reference still measures erosion. The band here is loose
        # enough that closing the gap passes and only a collapse fails, so the
        # test guards the shape of the answer without freezing the defect in.
        case = _psi_case("apollo12_erosion_rate")
        thrust = Float64(case.conditions["thrust_n"])
        rows = _psi_scored(case)
        rates = [(Float64(r["height_m"]),
                  PSIH.model_quantity(cfg, "erosion_kg_s", thrust, Float64(r["height_m"])),
                  Float64(r["value"])) for r in rows]
        for (h, modelval, refval) in rates
            @test modelval !== nothing
            @test modelval >= 0.0
            h <= 30.0 || continue        # above the model's onset it is identically zero
            @test 0.02 <= modelval / refval <= 3.0
        end
        # Below the onset the model erodes something everywhere the reference does.
        @test all(m > 0.0 for (h, m, _) in rates if h <= 30.0)
    end
end
