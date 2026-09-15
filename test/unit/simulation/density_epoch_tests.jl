using Test
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel

@testset "DensityModelEpoch" begin
    planet = make_no_gram_planet(:mars)
    it = SM.InitialTime(year=1993, month=5, day=26, hour=1, minute=36, second=57.0)
    for model in (SM.ExponentialAtmosphereModel(planet), SM.NoAtmosphereModel())
        @test SM.with_density_model_epoch(model, it) === model
    end
    # run_simulation leaves a configuration without an epoch-bearing density model untouched
    bus = SM.Link(root=true, m=100.0, dims=SM.MVector{3, Float64}(1.0, 1.0, 1.0), ref_area=1.0)
    sc = SM.SpacecraftModel(root=bus, initial_condition=SM.InitialCondition(ra=planet.Rp_e + 300e3, rp=planet.Rp_e + 300e3, i=10.0, ω=0.0, Ω=0.0, ν=0.0), id=1)
    args = SpaceAGORA.TelemetryVerification.make_example_config(planet=planet, spacecraft=sc, mission_time=10.0, initial_time=it,
        density_model=SM.ExponentialAtmosphereModel(planet), ephemerides_model=SM.SimpleEphemeridesModel(), verbose=false, results=false)
    @test SpaceAGORA.SimulationEngine._with_density_model_epoch(args) === args

    gram_ext = try
        @eval import GRAMSuite
        Base.get_extension(SpaceAGORA, :SpaceAGORAGRAMSuiteExt)
    catch
        nothing
    end
    if gram_ext === nothing || !isfile(joinpath(@__DIR__, "..", "..", "..", "data", "GRAMSuite.jl", "GRAM Suite 2.0", "Build", "lib", "libGRAM.so"))
        @info "GRAM not available; skipping the GRAM epoch alignment test"
    else
        model = SM.GRAMAtmosphereModel(planet_name="venus")
        @test model.constructor_kwargs[:planet_name] == "venus"
        @test Int(model.core.initial_time.year) == 2000
        aligned = SM.with_density_model_epoch(model, it)
        @test aligned !== model
        @test Int(aligned.core.initial_time.year) == 1993 && Int(aligned.core.initial_time.day) == 26
        @test aligned.constructor_kwargs[:planet_name] == "venus"
        @test SM.with_density_model_epoch(aligned, it) === aligned
        # an explicit epoch at construction needs no rebuild
        built = SM.GRAMAtmosphereModel(planet_name="venus", initial_time=it)
        @test SM.with_density_model_epoch(built, it) === built
        # the run configuration is rebuilt only when needed
        venus = SM.Venus("", joinpath(@__DIR__, "..", "..", "..", "data", "GRAMSuite.jl", "GRAM Suite 2.0", "SPICE"))
        sc_v = SM.SpacecraftModel(root=bus, initial_condition=SM.InitialCondition(ra=venus.Rp_e + 300e3, rp=venus.Rp_e + 300e3, i=10.0, ω=0.0, Ω=0.0, ν=0.0), id=1)
        args_v = SpaceAGORA.TelemetryVerification.make_example_config(planet=venus, spacecraft=sc_v, mission_time=10.0, initial_time=it, density_model=model, verbose=false, results=false)
        rebuilt = SpaceAGORA.SimulationEngine._with_density_model_epoch(args_v)
        @test rebuilt !== args_v
        @test Int(rebuilt.environment_model.density_model.core.initial_time.year) == 1993
        @test rebuilt.environment_model.EI == args_v.environment_model.EI
        @test SpaceAGORA.SimulationEngine._with_density_model_epoch(rebuilt) === rebuilt
    end
end
