using Test, Serialization, StaticArrays, Libdl
const SM = SpaceAGORA.SimulationModel
const SV = SM.SceneVisualization
import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft

function grid_file(path; heights=[60.0,90.0], latitudes=[-90.0,0.0,90.0])
    longitudes=[0.0,90.0,180.0,270.0]
    keys=("density_kgm3","temperature_K","wind_ew_ms","wind_ns_ms","wind_up_ms")
    fields(h,a,b)=(1e-8*(2+h/100+a/180+0.2cosd(b)),200.0,0.0,0.0,0.0)
    serialize(path,Dict{String,Any}(
        "status"=>"ok","type"=>"surrogate_trilinear","planet"=>"mars",
        "grid"=>Dict("alt_km"=>heights,"lat_deg"=>latitudes,"lon_deg"=>longitudes),
        "fields"=>Dict(k=>[fields(h,a,b)[i] for h in heights,a in latitudes,b in longitudes]
            for (i,k) in enumerate(keys)),"provenance"=>"synthetic visualization regression"))
    return path
end
function config(model,dir; ei=100.0)
    planet=make_no_gram_planet(:mars)
    spacecraft=make_three_body_spacecraft(
        bus_dims=(2.0,2.0,2.5),panel_dims=(0.01,2.5,1.0),bus_mass=500.0,
        panel_mass_each=10.0,panel_offset_y=2.3,prop_mass=0.0,id=1,
        ic=SM.InitialCondition(ra=4_500e3,rp=3_800e3,i=30.0,ω=0.0,Ω=45.0,ν=0.0))
    return make_example_config(planet=planet,spacecraft=spacecraft,mission_time=2.0,
        initial_time=SM.InitialTime(year=2020,month=1,day=1),
        dynamic_effectors=(SM.InverseSquaredGravityModel(),),density_model=model,
        ephemerides_model=SM.SimpleEphemeridesModel(),orientation_sim=false,
        keplerian=true,EI_km=ei,verbose=false,results=true,results_directory=dir)
end

@testset "automatic scene sampling uses fixed grids without native GRAM" begin
    mktempdir() do dir
        model=SM.GRAMGridAtmosphereModel(planet="mars",surrogate_file=grid_file(joinpath(dir,"grid.jls")))
        args=with_visualization_scene(config(model,dir),true)
        before=deepcopy(model.core.surrogate.rho)
        spec=atmosphere_spec(args;profile_points=7,map_step_deg=30.0)
        @test spec.profile_altitude_m == collect(range(60e3,90e3;length=7))
        @test spec.profile_density_kg_m3 ≈ [1e-8*(2+h/1e5+0.2) for h in spec.profile_altitude_m]
        @test spec.map_altitude_m == 60e3
        @test length(spec.map_lat_deg)==6
        @test length(spec.map_lon_deg)==12
        @test length(spec.map_density_kg_m3)==72
        for (i,a) in enumerate(spec.map_lat_deg),(j,b) in enumerate(spec.map_lon_deg)
            @test spec.map_density_kg_m3[(i-1)*12+j] ==
                SM.getDensity(model,60e3,deg2rad(a),deg2rad(b),0.0,false)[1]
        end
        scene=build_visualization_scene(args;rotation_max_samples=4)
        @test length(scene.atmosphere.profile_altitude_m)==64
        @test !isempty(scene.atmosphere.map_density_kg_m3)
        path=SV.write_visualization_scene!(args)
        @test read_visualization_scene(path)==scene
        @test model.core.surrogate.rho==before
        @test GRAMSuite._GRAM_WRAPPER[]===nothing
        @test !any(p->occursin("libgram",lowercase(p)),Libdl.dllist())
        # Caller-specified map heights are never silently moved into the grid.
        @test_throws ArgumentError atmosphere_spec(args;map_altitude_m=50e3)
        @test_throws ArgumentError atmosphere_spec(args;map_altitude_m=91e3)
        @test atmosphere_spec(args;map_altitude_m=75e3).map_altitude_m==75e3
        # The default display height may be bounded, and the sidecar records it.
        @test atmosphere_spec(config(model,dir;ei=200.0)).map_altitude_m==90e3
        for bad in (0.0,-1.0,Inf,NaN,181.0)
            @test_throws ArgumentError atmosphere_spec(args;map_step_deg=bad)
        end
        @test_throws ArgumentError atmosphere_spec(args;map_altitude_m=NaN)
        # The display never enables the native/hybrid wrapper by broad type matching.
        raw=SM.GRAMAtmosphereModel(nothing)
        @test isempty(atmosphere_spec(config(raw,dir)).profile_density_kg_m3)
        unknown=SM.GRAMGridAtmosphereModel((unexpected_core=true,))
        @test isempty(atmosphere_spec(config(unknown,dir)).profile_density_kg_m3)
        regional=SM.GRAMGridAtmosphereModel(planet="mars",
            surrogate_file=grid_file(joinpath(dir,"regional.jls");latitudes=[-30.0,0.0,30.0]))
        spec=atmosphere_spec(config(regional,dir);map_step_deg=30.0)
        @test length(spec.profile_density_kg_m3)==64
        @test isempty(spec.map_density_kg_m3) # no stretching a regional map over the globe
        northern=SM.GRAMGridAtmosphereModel(planet="mars",
            surrogate_file=grid_file(joinpath(dir,"north.jls");latitudes=[10.0,20.0]))
        @test isempty(atmosphere_spec(config(northern,dir)).profile_density_kg_m3)
        high=SM.GRAMGridAtmosphereModel(planet="mars",
            surrogate_file=grid_file(joinpath(dir,"high.jls");heights=[500.0,600.0]))
        @test isempty(atmosphere_spec(config(high,dir)).profile_density_kg_m3)
    end
end
