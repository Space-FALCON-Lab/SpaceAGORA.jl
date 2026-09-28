module TerrainExportTests
using Test, SpaceAGORA, JSON, Base64, StaticArrays, DataFrames, Arrow
import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft
const SM=SpaceAGORA.SimulationModel
const SV=SM.SceneVisualization
const TM=SM.TerrainModels
const R=1_737_400.0
const TILE_BASE64="/9j/4AAQSkZJRgABAQAAAQABAAD/2wBDAAMCAgMCAgMDAwMEAwMEBQgFBQQEBQoHBwYIDAoMDAsKCwsNDhIQDQ4RDgsLEBYQERMUFRUVDA8XGBYUGBIUFRT/2wBDAQMEBAUEBQkFBQkUDQsNFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBT/wAARCAAIAAgDASIAAhEBAxEB/8QAHwAAAQUBAQEBAQEAAAAAAAAAAAECAwQFBgcICQoL/8QAtRAAAgEDAwIEAwUFBAQAAAF9AQIDAAQRBRIhMUEGE1FhByJxFDKBkaEII0KxwRVS0fAkM2JyggkKFhcYGRolJicoKSo0NTY3ODk6Q0RFRkdISUpTVFVWV1hZWmNkZWZnaGlqc3R1dnd4eXqDhIWGh4iJipKTlJWWl5iZmqKjpKWmp6ipqrKztLW2t7i5usLDxMXGx8jJytLT1NXW19jZ2uHi4+Tl5ufo6erx8vP09fb3+Pn6/8QAHwEAAwEBAQEBAQEBAQAAAAAAAAECAwQFBgcICQoL/8QAtREAAgECBAQDBAcFBAQAAQJ3AAECAxEEBSExBhJBUQdhcRMiMoEIFEKRobHBCSMzUvAVYnLRChYkNOEl8RcYGRomJygpKjU2Nzg5OkNERUZHSElKU1RVVldYWVpjZGVmZ2hpanN0dXZ3eHl6goOEhYaHiImKkpOUlZaXmJmaoqOkpaanqKmqsrO0tba3uLm6wsPExcbHyMnK0tPU1dbX2Nna4uPk5ebn6Onq8vP09fb3+Pn6/9oADAMBAAIRAxEAPwDg6KKK+aPpz//Z"
_write_json(path,value)=(write(path,JSON.json(value));path)
function _grid(dir,name,rows,cols; south=-2.0,north=3.0,west=350.0,east=371.0,offset=100.0,radius=R)
    height(lat,lon)=offset+2lat+3(lon-west)
    h=Float32[height(north-(r-0.5)*(north-south)/rows,west+(c-0.5)*(east-west)/cols) for r=1:rows,c=1:cols]
    meta=Dict("rows"=>rows,"cols"=>cols,"lat_min"=>south,"lat_max"=>north,"lon_min"=>west,"lon_max"=>east,"reference_radius_m"=>radius,"source"=>"synthetic plane","units"=>"m")
    _write_json(joinpath(dir,name*".json"),meta)
    bytes=UInt8[]
    for row=1:rows,col=1:cols
        bits=reinterpret(UInt32,h[row,col])
        append!(bytes,UInt8[(bits>>shift)&0xff for shift in (0,8,16,24)])
    end
    write(joinpath(dir,name*".f32"),bytes)
    return Dict("name"=>name,"reference_radius_m"=>radius),h
end
_site(dir,entries; tiles=nothing,lat=0.5,lon=360.5)=_write_json(joinpath(dir,"site.json"),Dict("site"=>Dict("lat_deg"=>lat,"lon_deg"=>lon,"name"=>"synthetic site"),"dem"=>entries,"tiles"=>tiles))
_decode(g)=collect(reinterpret(Float32,base64decode(g["heights"])))
function _tiles(dir)
    mkdir(joinpath(dir,"imagery"));write(joinpath(dir,"imagery","root.jpg"),base64decode(TILE_BASE64))
    return Dict{String,Any}("scheme"=>"quadtree","tile_px"=>8,"max_level"=>0,
        "root"=>Dict{String,Any}("lat_min"=>-2.0,"lat_max"=>3.0,"lon_min"=>350.0,"lon_max"=>371.0),
        "nodes"=>[Dict{String,Any}("level"=>0,"x"=>0,"y"=>0,"m_per_px"=>10.0,"file"=>"root.jpg")],
        "attribution"=>["Synthetic test image"],"source"=>"generated test pixels",
        "resolution"=>Dict("finest_m_per_px"=>10.0,"feature_scale_m"=>20.0))
end
function _scene(dir)
    craft=make_three_body_spacecraft(bus_dims=(1.0,1.0,1.0),panel_dims=(0.01,1.0,0.5),bus_mass=100.0,panel_mass_each=1.0,panel_offset_y=1.0,
        ic=SM.InitialCondition(ra=4.5e6,rp=3.8e6,i=0.0,ω=0.0,Ω=0.0,ν=0.0),prop_mass=0.0,id=1)
    args=make_example_config(planet=make_no_gram_planet(:mars),spacecraft=craft,mission_time=10.0,initial_time=SM.InitialTime(year=2020,month=1,day=1,hour=0,minute=0,second=0.0),dynamic_effectors=(SM.ConstantGravityModel(),),density_model=SM.NoAtmosphereModel(),ephemerides_model=SM.SimpleEphemeridesModel(),orientation_sim=false,keplerian=true,verbose=false,results=true,results_directory=dir)
    return build_visualization_scene(args;rotation_max_samples=2)
end

function _page_payload(path)
    html=read(path,String)
    start=findfirst("window.SPACEAGORA_VIEWER = ",html)
    stop=findnext(";\n</script>",html,last(start))
    return JSON.parse(html[last(start)+1:first(stop)-1])
end

function _combined_export(dir,scene,sitepath)
    # Two IDs distinct from their row indices exercise channel indexing and
    # model lookup together. One OBJ is shared, but transforms stay separate.
    g=only(scene.spacecraft)
    fleet=SV.VisualizationScene(scene.schema,scene.epoch_et_start_s,scene.epoch_utc,scene.planet,
        [SV.SpacecraftGeometry(id,"terrain craft $(id)",g.links,g.thrusters,g.facets,g.joints,
            g.bounding_radius_m,g.stl_path,g.arm) for id in (11,22)],
        scene.orientation_sim,"combined.feather",scene.link_pose_field,scene.link_pose_stride,scene.atmosphere)
    df=DataFrame(time=collect(0.0:4.0),
        sc1_pos_1=fill(R+1000,5),sc1_pos_2=zeros(5),sc1_pos_3=zeros(5),
        sc2_pos_1=fill(R+2000,5),sc2_pos_2=zeros(5),sc2_pos_3=zeros(5),
        sc1_metric=[1.0,2.0,NaN,4.0,5.0],sc2_metric=[10.0,20.0,30.0,40.0,50.0])
    obj=joinpath(dir,"shared.obj")
    write(obj,"v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n")
    options=(max_frames=3,channels=[(column="metric",label="Saved metric",unit="m",digits=1,log=false)],
        models=Dict(11=>obj,22=>obj),model_scale=Dict(11=>1.0,22=>2.0),
        model_rotation_deg=Dict(22=>(0.0,0.0,90.0)))
    plain=SV.viewer_payload(fleet,df;include_textures=false,options...)
    SV.write_visualization_scene(joinpath(dir,"combined_scene.json"),fleet)
    Arrow.write(joinpath(dir,"combined.feather"),df)
    page=export_visualization(joinpath(dir,"combined");textures=false,
        terrain=sitepath,terrain_max_grid=3,options...)
    combined=_page_payload(page)
    @test plain["terrain"]===nothing
    @test combined["terrain"]==SV.terrain_payload(sitepath;max_grid=3)
    @test combined["frames"]==plain["frames"]
    @test combined["models"]==plain["models"]
    @test combined["scene"]==plain["scene"]
    @test combined["frames"]["sats"]==2
    @test combined["frames"]["count"]==3
    @test collect(reinterpret(Float64,base64decode(combined["frames"]["t_s"])))==[0.0,2.0,4.0]
    channel=only(combined["frames"]["channels"])
    @test (channel["name"],channel["label"],channel["unit"],channel["digits"],channel["log"])==
        ("metric","Saved metric","m",1,false)
    @test isequal(collect(reinterpret(Float32,base64decode(channel["data"]))),
        Float32[1,10,NaN,30,5,50])
    @test haskey(combined["models"]["11"],"url")
    @test !haskey(combined["models"]["22"],"url")
    @test combined["models"]["22"]["url_from"]=="11"
    @test [combined["models"][string(id)]["scale"] for id in (11,22)]==[1.0,2.0]
    @test combined["models"]["11"]["rotation_deg"]==[0.0,0.0,0.0]
    @test combined["models"]["22"]["rotation_deg"]==[0.0,0.0,90.0]
    return combined
end

@testset "Terrain export contracts and cross-language fixtures" begin
    cases=Any[]
    mktempdir() do dir
        for (label,rows,cols,limit) in (("nonsquare",3,5,512),("odd_reduction",5,7,3),("single_row",2,7,3),("single_col",7,2,3),("singleton",5,7,1))
            sub=joinpath(dir,label);mkdir(sub)
            entry,source=_grid(sub,"grid",rows,cols);path=_site(sub,[entry])
            model,_=load_site_terrain(path);payload=SV.terrain_payload(path;max_grid=limit);g=only(payload["grids"])
            factor=max(1.0,max(rows,cols)/limit);r2=max(1,floor(Int,rows/factor));c2=max(1,floor(Int,cols/factor))
            @test (g["rows"],g["cols"])==(r2,c2)
            @test [g[k] for k in ("lat_min","lat_max","lon_min","lon_max")]==[-2.0,3.0,350.0,371.0]
            @test payload["reference_radius_m"]==R && payload["tiles"]===nothing
            @test payload["site"]["height_m"]==terrain_height(model,0.5,360.5)
            values=_decode(g);queries=Any[]
            for r=1:r2,c=1:c2
                lat=3.0-(r-0.5)*5/r2;lon=350.0+(c-0.5)*21/c2
                expected=Float32(terrain_height(model,lat,lon))
                @test values[(r-1)*c2+c]==expected
                push!(queries,Dict("lat"=>lat,"lon"=>lon,"expected"=>Float64(expected)))
                push!(queries,Dict("lat"=>lat,"lon"=>lon-360,"expected"=>Float64(expected)))
            end
            # Edges clamp to the exported edge samples, not shifted original locations.
            for (lat,lon,k) in ((3.0,350.0,1),(-2.0,371.0,r2*c2),(3.0,371.0,c2),(-2.0,350.0,(r2-1)*c2+1))
                push!(queries,Dict("lat"=>lat,"lon"=>lon,"expected"=>Float64(values[k])))
            end
            for (lat,lon) in ((4.0,360.0),(0.0,prevfloat(350.0)),(0.0,nextfloat(11.0)))
                push!(queries,Dict("lat"=>lat,"lon"=>lon,"expected"=>0.0))
            end
            if label=="nonsquare"
                for lon in (1e300,-1e300,floatmax(Float64))
                    push!(queries,Dict("lat"=>0.5,"lon"=>lon,"expected"=>terrain_height(model,0.5,lon)))
                end
            end
            push!(cases,Dict("name"=>label,"payload"=>payload,"queries"=>queries))
        end
        sub=joinpath(dir,"priority");mkdir(sub)
        first_entry,_=_grid(sub,"first",3,5;offset=1000.0)
        second_entry,_=_grid(sub,"second",7,9;offset=0.0)
        priority=SV.terrain_payload(_site(sub,[first_entry,second_entry]))
        @test [g["name"] for g in priority["grids"]]==["first","second"]
        model,_=load_site_terrain(joinpath(sub,"site.json"))
        push!(cases,Dict("name"=>"priority","payload"=>priority,"queries"=>[Dict("lat"=>0.5,"lon"=>360.5,"expected"=>terrain_height(model,0.5,360.5))]))
        sub=joinpath(dir,"zero_edge");mkdir(sub)
        entry,_=_grid(sub,"grid",2,2;south=0.0,north=2.0,west=0.0,east=20.0)
        zero=SV.terrain_payload(_site(sub,[entry];lat=1.0,lon=10.0));model,_=load_site_terrain(joinpath(sub,"site.json"))
        qs=[Dict("lat"=>1.0,"lon"=>lon,"expected"=>terrain_height(model,1.0,lon)) for lon in (-nextfloat(0.0),0.0,20.0,nextfloat(20.0),prevfloat(360.0),360.0)]
        push!(cases,Dict("name"=>"zero_edge","payload"=>zero,"queries"=>qs))

        sub=joinpath(dir,"imagery");mkdir(sub);entry,_=_grid(sub,"grid",5,7);index=_tiles(sub)
        indexpath=_write_json(joinpath(sub,"imagery","tiles.json"),index);sitepath=_site(sub,[entry];tiles="imagery/tiles.json")
        payload=SV.terrain_payload(sitepath;max_grid=3)
        @test payload["tiles"]["root"]["lon_min"]==350.0
        @test payload["tiles"]["attribution"]==["Synthetic test image"]
        @test base64decode(split(only(payload["tiles"]["nodes"])["url"],",";limit=2)[2])==base64decode(TILE_BASE64)
        @test payload["tiles"]["resolution"]["feature_scale_m"]==20.0
        push!(cases,Dict("name"=>"imagery","payload"=>payload,"queries"=>Any[]))
        for limit in (0,-1,true,1.5,2049,typemax(Int),"512")
            @test_throws ArgumentError SV.terrain_payload(sitepath;max_grid=limit)
        end
        @test_throws ArgumentError SV.terrain_tiles_payload(sub,"missing.json")
        @test_throws ArgumentError SV.terrain_tiles_payload(sub,indexpath)
        @test_throws ArgumentError SV.terrain_tiles_payload(sub,"../imagery/imagery/tiles.json")
        invalids=Any[]
        for (key,value) in (("lat_min",4.0),("lat_max",91.0),("lon_max",350.0),("lon_max",710.0),("lat_min",true),("lon_min",Inf))
            bad=deepcopy(index);bad["root"][key]=value;push!(invalids,bad)
        end
        for (key,value) in (("level",-1),("level",21),("level",true),("x",1),("y",-1),("x",0.5),("m_per_px",0.0),("m_per_px","1"),("file","missing.jpg"),("file","../../grid.f32"))
            bad=deepcopy(index);bad["nodes"][1][key]=value;push!(invalids,bad)
        end
        bad=deepcopy(index);push!(bad["nodes"],copy(bad["nodes"][1]));push!(invalids,bad)
        for (key,value) in (("max_level",1),("nodes",Any[]),("tile_px",9),("tile_px",true),("scheme","unknown"),("attribution","bad"),("source",9))
            bad=deepcopy(index);bad[key]=value;push!(invalids,bad)
        end
        for bad in invalids
            _write_json(indexpath,bad)
            @test_throws Exception SV.terrain_payload(sitepath)
        end
        _write_json(indexpath,index)
        # A declared file cannot escape through a symlink, even under a safe relative name.
        mktempdir() do outside
            write(joinpath(outside,"tile.jpg"),base64decode(TILE_BASE64))
            symlink(joinpath(outside,"tile.jpg"),joinpath(sub,"imagery","escape.jpg"))
            bad=deepcopy(index);bad["nodes"][1]["file"]="escape.jpg";_write_json(indexpath,bad)
            @test_throws ArgumentError SV.terrain_payload(sitepath)
        end
        _write_json(indexpath,index)
        scene=_scene(sub);df=DataFrame(time=[0.0,10.0],sc1_pos_1=[R+1000,R+500],sc1_pos_2=[0.,0.],sc1_pos_3=[0.,0.])
        raw=SV.scene_dict(scene);raw["results"]["feather"]="terrain.feather";scene=SV._scene_from_dict(raw)
        SV.write_visualization_scene(joinpath(sub,"terrain_scene.json"),scene);Arrow.write(joinpath(sub,"terrain.feather"),df)
        page=export_visualization(joinpath(sub,"terrain");textures=false,terrain=sitepath,terrain_max_grid=3,title="Synthetic terrain acceptance")
        html=read(page,String);a=findfirst("window.SPACEAGORA_VIEWER = ",html);z=findnext(";\n</script>",html,last(a));exported=JSON.parse(html[last(a)+1:first(z)-1])
        @test exported["terrain"]==payload
        @test exported["scene"]["planet"]["equatorial_radius_m"]==scene.planet.equatorial_radius_m
        @test exported["terrain"]["reference_radius_m"] != scene.planet.equatorial_radius_m
        @test occursin("MIT License",html) && occursin("mp4-muxer",html)
        @test !occursin("__PAYLOAD__",html)
        combined=_combined_export(sub,scene,sitepath)
        output=get(ENV,"SPACEAGORA_TERRAIN_TEST_FIXTURE_DIR","")
        if !isempty(output)
            mkpath(output);_write_json(joinpath(output,"fixtures.json"),cases)
            _write_json(joinpath(output,"combined-export.json"),combined)
            # The browser artifact is explicitly synthetic and uses a matching sphere.
            browser_raw=SV.scene_dict(scene)
            browser_raw["planet"]["name"]="Synthetic Moon terrain"
            browser_raw["planet"]["texture"]="moon"
            browser_raw["planet"]["equatorial_radius_m"]=R
            browser_raw["planet"]["polar_radius_m"]=R
            SV.write_visualization_scene(joinpath(sub,"terrain_scene.json"),SV._scene_from_dict(browser_raw))
            export_visualization(joinpath(sub,"terrain");out=joinpath(output,"terrain-browser-fixture.html"),textures=false,terrain=sitepath,terrain_max_grid=3,title="Synthetic terrain acceptance fixture")
        end
    end
end
end
