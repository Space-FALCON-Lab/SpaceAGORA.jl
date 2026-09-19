using Test, Serialization

function small_grid_payload()
    Dict{String,Any}("status"=>"ok","planet"=>"mars","type"=>"full_grid","format"=>"spaceagora_gram_static_grid_v1",
        "grid"=>Dict("alt_km"=>[100.,260.],"lat_deg"=>[40.,90.],"lon_deg"=>[0.,180.]),
        "fields"=>Dict(k=>fill(v,2,2,2) for (k,v) in zip(("density_kgm3","temperature_K","wind_ew_ms","wind_ns_ms","wind_up_ms"),(1e-8,180.,2.,3.,4.))),
        "source_note"=>"Synthetic software fixture")
end

@testset "Full-grid schema compatibility and legacy restriction" begin
    mktempdir() do dir
        path=joinpath(dir,"small.jls");payload=small_grid_payload()
        save()=serialize(path,payload)
        save();model=SpaceAGORA.GRAMGridAtmosphereModel(planet="mars",surrogate_file=path)
        @test model.core.metadata["type"]=="full_grid"
        @test_throws ArgumentError GRAMSuite._gram_load_offline_surrogate(path,"mars")
        for bad in (nothing,"unknown","spaceagora_gram_surrogate_trilinear_v1")
            bad===nothing ? delete!(payload,"format") : (payload["format"]=bad);save()
            @test_throws ArgumentError SpaceAGORA.GRAMGridAtmosphereModel(planet="mars",surrogate_file=path)
        end
        payload=small_grid_payload();payload["type"]="surrogate_trilinear"
        for bad in ("unknown","spaceagora_gram_static_grid_v1")
            payload["format"]=bad;save()
            @test_throws ArgumentError SpaceAGORA.GRAMGridAtmosphereModel(planet="mars",surrogate_file=path)
        end
        for format in (nothing,"spaceagora_gram_surrogate_trilinear_v1")
            format===nothing ? delete!(payload,"format") : (payload["format"]=format);save()
            @test SpaceAGORA.GRAMGridAtmosphereModel(planet="mars",surrogate_file=path).core.metadata["type"]=="surrogate_trilinear"
            @test GRAMSuite._gram_load_offline_surrogate(path,"mars").planet_name=="mars"
        end
        for mutate in (p->p["grid"]["lat_deg"][:]=[90.,40.], p->p["fields"]["density_kgm3"][1]=NaN,
                       p->p["fields"]["temperature_K"][1]=0.,p->p["fields"]["wind_ew_ms"][1]=Inf,
                       p->p["type"]="unknown",p->p["status"]="error",p->p["planet"]="earth")
            payload=small_grid_payload();mutate(payload);save()
            @test_throws ArgumentError SpaceAGORA.GRAMGridAtmosphereModel(planet="mars",surrogate_file=path)
        end
    end
end
