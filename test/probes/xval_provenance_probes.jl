using Test, SpaceAGORA, JSON
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const EPHEMERIDES = parentmodule(first(methods(SpaceAGORA.SimulationModel.planet_frame_lpi)))
const ORIGINAL_FRAME_METHOD = first(methods(EPHEMERIDES._spice_planet_frame_lpi))
# Loading definitions must not perform a campaign or replace the package method.
withenv("XVAL_FRAME_TABLE_DIR" => nothing) do
    include(joinpath(ROOT, "scripts", "xval_fullarc.jl"))
end
const XP = XvalProvenance

@testset "Full-arc frame isolation" begin
    @test first(methods(EPHEMERIDES._spice_planet_frame_lpi)) === ORIGINAL_FRAME_METHOD
    @test _FRAME_TABLE[] === nothing
    withenv("XVAL_FRAME_TABLE_DIR" => "changed-after-load") do
        @test_throws ErrorException _install_frame_table!("earth_j2_tbfalse")
    end
    @test first(methods(EPHEMERIDES._spice_planet_frame_lpi)) === ORIGINAL_FRAME_METHOD
end

# Synthetic lifecycle cases isolate CI controls from the committed-run policy.
withenv((key => nothing for key in keys(XP.controls()))...) do
@testset "Full-arc provenance lifecycle" begin
    mktempdir() do root
        input = joinpath(root, "input.csv")
        write(input, "immutable input\n")
        identity = XP.file_identity(input)
        @test length(identity["sha256"]) == 64
        @test XP.verify_inputs([identity]) === nothing
        write(input, "changed input\n")
        @test_throws ErrorException XP.verify_inputs([identity])
        identity = XP.file_identity(input)
        rundir = joinpath(root, "run")
        record = XP.start_record(ROOT, rundir, "gmat", "committed", ["earth_j0_tbfalse"])
        @test JSON.parsefile(joinpath(rundir, "run_info.json"))["status"] == "running"
        @test_throws ErrorException XP.start_record(ROOT, rundir, "gmat", "committed", [])
        @test_throws ErrorException XP.finish_record!(record, ROOT, rundir)
        record["cases"]["earth_j0_tbfalse"] = Dict("inputs" => [identity], "kernels" => [identity])
        write(joinpath(rundir, "results.csv"), "fixture\n")
        withenv("XVAL_J0_GM" => "123") do
            @test_throws ErrorException XP.finish_record!(record, ROOT, rundir)
        end
        @test JSON.parsefile(joinpath(rundir, "run_info.json"))["status"] == "running"
        XP.finish_record!(record, ROOT, rundir)
        saved = JSON.parsefile(joinpath(rundir, "run_info.json"))
        @test saved["status"] == "complete"
        @test saved["reference_provenance"] == "unverified_generation_settings"
        @test saved["source"] == saved["source_end"]
        @test saved["results_sha256"] == XP.file_identity(joinpath(rundir, "results.csv"))["sha256"]
        write(input, "changed again\n")
        @test_throws ErrorException XP.finish_record!(record, ROOT, rundir)
    end
end

end # isolated lifecycle environment

@testset "Diagnostic frame replacement requires explicit opt-in" begin
    mktempdir() do root
        write(joinpath(root, "Earth_E_spice_to_gmat.csv"),
              "0,0,-1,0,1,0,0,0,0,1\n1,0,-1,0,1,0,0,0,0,1\n")
        code = """
            using Test, SpaceAGORA, SPICE
            eph = parentmodule(first(methods(SpaceAGORA.SimulationModel.planet_frame_lpi)))
            original = first(methods(eph._spice_planet_frame_lpi))
            include($(repr(joinpath(ROOT, "scripts", "xval_fullarc.jl"))))
            @test first(methods(eph._spice_planet_frame_lpi)) !== original
            _install_frame_table!("earth_j2_tbfalse")
            SPICE.pdpool("BODY399_POLE_RA", [0., 0., 0.])
            SPICE.pdpool("BODY399_POLE_DEC", [90., 0., 0.])
            SPICE.pdpool("BODY399_PM", [0., 1., 0.])
            expected = SPICE.pxform("J2000", "IAU_EARTH", 0.5) * [0. -1. 0.; 1. 0. 0.; 0. 0. 1.]
            @test eph._spice_planet_frame_lpi(SpaceAGORA.SimulationModel.Earth(), 0.5) ≈ expected
            withenv("XVAL_FRAME_TABLE_DIR" => nothing) do
                @test_throws ErrorException _install_frame_table!("earth_j2_tbfalse")
            end
            SPICE.kclear()
            """
        command = `$(Base.julia_cmd()) --startup-file=no --compiled-modules=existing --project=$ROOT -e $code`
        @test success(addenv(command, "XVAL_FRAME_TABLE_DIR" => root))
    end
end
