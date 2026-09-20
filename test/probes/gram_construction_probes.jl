# Native probe: standalone execution requires GRAM. The threaded coverage
# driver includes it only inside its existing HAS_GRAMSUITE native gate.
module GRAMConstructionProbes
using Test, Serialization, SpaceAGORA
if Base.find_package("GRAMSuite") === nothing
    pushfirst!(LOAD_PATH, joinpath(dirname(dirname(pathof(SpaceAGORA))), "data", "GRAMSuite.jl"))
end
using GRAMSuite

const SM = SpaceAGORA.SimulationModel
const EM = SM.EnvironmentModels

function roundtrip(model)
    io = IOBuffer()
    serialize(io, model)
    seekstart(io)
    return deserialize(io)
end

mutable struct CallerEpoch
    year::Int
    month::Int
    day::Int
    hour::Int
    minute::Int
    second::Float64
end

@testset "GRAM construction settings survive copy and transfer" begin
    @test Base.get_extension(SpaceAGORA, :SpaceAGORAGRAMSuiteExt) !== nothing
    withenv("SPACEAGORA_GRAM_STATIC_GRID" => "0") do
        epoch = CallerEpoch(2020,1,2,3,4,5.5)
        model = EM.GRAMAtmosphereModel(
            planet_name="earth", initial_time=epoch,
            gram_perturbation_scales=(0.0,0.0,0.0,0.0),
            gram_min_relative_step_size=0.02)
        recipe = deepcopy(model.constructor_kwargs)
        @test !model.offline_surrogate_supported
        @test recipe[:gram_perturbation_scales] == (0.0,0.0,0.0,0.0)
        @test recipe[:gram_min_relative_step_size] == 0.02
        @test recipe[:gram_root_directory] == model.gram_root
        @test recipe[:gram_data_directory] == model.gram_data_root
        @test recipe[:spice_directory] == model.spice_root
        epoch.year = 2035
        @test Int(model.initial_time.year) == 2020
        @test Int(recipe[:initial_time].year) == 2020

        # Load the same public SPICE assets used by the native model. Compare
        # native mean-state evaluations as well as the retained recipe.
        SM.Earth("", model.spice_root)
        sites = ((150e3,0.0,0.0,0.0), (200e3,0.1,0.2,60.0))
        reference = [GRAMSuite.point_density_state(model.core, site..., false) for site in sites]
        @test all(value -> isfinite(value[1]) && value[1] > 0.0 && isfinite(value[2]), reference)
        for transform in (deepcopy, roundtrip)
            copied = transform(model)
            @test !copied.offline_surrogate_supported
            @test copied.offline_surrogate_unsupported_reason == model.offline_surrogate_unsupported_reason
            @test copied.constructor_kwargs == recipe
            @test copied.constructor_kwargs !== model.constructor_kwargs
            @test copied.instance_lock !== model.instance_lock
            @test copied.gram_atmosphere !== model.gram_atmosphere
            @test copied.initial_time == model.initial_time
            values = [GRAMSuite.point_density_state(copied.core, site..., false) for site in sites]
            @test values == reference
            copied.constructor_kwargs[:gram_min_relative_step_size] = 0.5
            @test model.constructor_kwargs == recipe
        end

        graph = deepcopy((first=model,second=model))
        @test graph.first === graph.second
        @test graph.first.instance_lock !== model.instance_lock
        surrogate = EM.GRAMAtmosphereModelSurrogate(model, "fixture.arrow", 1234.0)
        for transform in (deepcopy, roundtrip)
            copied = transform(surrogate)
            @test copied.surrogate_file == "fixture.arrow"
            @test copied.point_fallback_below_m == 1234.0
            @test copied.base_model.constructor_kwargs == recipe
            @test !copied.base_model.offline_surrogate_supported
            @test copied.base_model.instance_lock !== model.instance_lock
        end

        # Unknown raw-core construction retains its previous contract. This
        # does not claim recovery of its custom native construction options.
        raw = EM.GRAMAtmosphereModel(model.core)
        @test raw.constructor_kwargs === nothing
        supplied_lock = ReentrantLock()
        @test EM.GRAMAtmosphereModel(model.core,supplied_lock).instance_lock === supplied_lock
        for transform in (deepcopy, roundtrip)
            copied = transform(raw)
            @test copied.constructor_kwargs === nothing
            @test copied.instance_lock !== raw.instance_lock
        end
    end
end
end
