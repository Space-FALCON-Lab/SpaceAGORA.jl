using Test, Serialization, SHA, StaticArrays, Libdl, CSV

# These optional inputs are a retained Odyssey software regression, not a
# redistributable fixture or a new comparison with native GRAM.
const INPUTS = Dict(name => abspath(expanduser(ENV["SPACEAGORA_TEST_GRID_" * name]))
                    for name in ("FILE", "POINTS", "REFERENCE"))
const INPUT_SHA = Dict(name => lowercase(ENV["SPACEAGORA_TEST_GRID_" * name * "_SHA256"])
                      for name in keys(INPUTS))
const FINE_FILE = INPUTS["FILE"]
const FINE_SHA = INPUT_SHA["FILE"]
const POINT_FILE = INPUTS["POINTS"]
const GRID_CSV = INPUTS["REFERENCE"]
state_tuple(x) = (x[1], x[2], Tuple(x[3])...)
function read_controlled_csv(path)
    table = CSV.File(path; types=String)
    [Dict(String(k) => String(coalesce(getproperty(row,k), "")) for k in propertynames(row)) for row in table]
end
function grid_digests(s)
    Dict(String(k) => bytes2hex(sha256(reinterpret(UInt8, vec(getproperty(s,k)))))
         for k in (:rho, :T, :wind_e, :wind_n, :wind_u))
end

function test_retained_odyssey_grid()
    SM=SpaceAGORA.SimulationModel;EM=SM.EnvironmentModels;CB=SM.SimulationCallbacks
    source_before=bytes2hex(open(sha256,FINE_FILE))
    @test source_before==FINE_SHA
    points=read_controlled_csv(POINT_FILE)
    @test length(points)==654
    @test length(unique(p["id"] for p in points))==654
    @test all(p["expected_spatial_domain"]=="inside" for p in points)
    retained_rows=read_controlled_csv(GRID_CSV)
    retained=Dict(r["id"]=>r for r in retained_rows if r["method"]=="grid_fine" && r["group"]=="trajectory_with_events")
    @test Set(keys(retained))==Set(p["id"] for p in points)
    model=SpaceAGORA.GRAMGridAtmosphereModel(planet="mars",surrogate_file=FINE_FILE,expected_sha256=FINE_SHA)
    @test model isa SM.AbstractDensityModel
    @test CB.density_model_threadsafe(model)
    @test !CB._is_gram_density_model(model)
    @test model.core.source_sha256==FINE_SHA
    @test model.core.metadata["type"]=="full_grid"
    @test model.core.metadata["format"]=="spaceagora_gram_static_grid_v1"
    @test size(model.core.surrogate.rho)==(161,51,144)
    payload=open(deserialize,FINE_FILE)
    @test model.core.metadata==Dict{String,Any}(String(k)=>v for (k,v) in payload if k ∉ ("grid","fields"))
    @test model.core.metadata["generation_config"]["is_planetocentric"]===false
    @test model.core.metadata["generation_config"]["start_time_frame"]==0
    axes=payload["grid"];fields=payload["fields"]
    hs=[parse(Float64,p["height_km"])*1000 for p in points]
    aa=[deg2rad(parse(Float64,p["geodetic_lat_deg"])) for p in points]
    bb=[deg2rad(parse(Float64,p["east_lon_deg"])) for p in points]
    times=[parse(Float64,p["elapsed_from_reference_et_s"]) for p in points]
    expected=[Tuple(parse(Float64,retained[p["id"]][f]) for f in ("rho","T","wind_e","wind_n","wind_u")) for p in points]
    before_fields=grid_digests(model.core.surrogate)
    scalar=[state_tuple(SpaceAGORA.getDensity(model,hs[i],aa[i],bb[i],times[i],true)) for i in eachindex(points)]
    @test isequal(scalar,expected)
    seven=[state_tuple(SpaceAGORA.getDensity(model,hs[i],aa[i],bb[i],times[i],true,nothing)) for i in eachindex(points)]
    @test isequal(seven,expected)
    for time_arg in (0.0,times), wind_flag in (false,true)
        rho=zeros(654);temperature=zeros(654);wind=fill(SVector(0.,0.,0.),654)
        @test SpaceAGORA.getDensityBatch!(rho,temperature,wind,model,hs,aa,bb,time_arg,wind_flag,nothing)===nothing
        @test isequal([(rho[i],temperature[i],Tuple(wind[i])...) for i in eachindex(rho)],expected)
    end
    @test all(isequal(state_tuple(SpaceAGORA.getDensity(model,hs[i],aa[i],bb[i],-1e8,false)),expected[i]) for i in eachindex(points))
    for (h,a) in ((99999.,aa[1]),(260001.,aa[1]),(hs[1],deg2rad(39.99)),(hs[1],deg2rad(90.01)))
        @test_throws DomainError SpaceAGORA.getDensity(model,h,a,bb[1],0.,true)
    end
    for b in (bb[1]-2pi,bb[1]+2pi)
        actual=state_tuple(SpaceAGORA.getDensity(model,hs[1],aa[1],b,0.,true))
        @test all(isapprox.(actual,expected[1];rtol=2e-13,atol=1e-13))
    end
    @test_throws DomainError SpaceAGORA.getDensity(model,hs[1],aa[1],bb[1],NaN,true)
    @test_throws DomainError SpaceAGORA.getDensity(model,hs[1],NaN,bb[1],0.,true)
    # Compare query arithmetic with the accepted interpolated rows. A nominal
    # degree node can incur roundoff during periodic radian interpolation, so
    # exact equality to raw storage is a separate loading check below.
    # Longitude labels at the pole do not define one unique physical ENU vector.
    for longitude_deg in (0.,90.,180.)
        longitude=deg2rad(longitude_deg)
        q=state_tuple(SpaceAGORA.getDensity(model,180e3,pi/2,longitude,0.,true))
        matches=filter(retained_rows) do r
            r["method"]=="grid_fine" && r["group"]=="poles" &&
                parse(Float64,r["height_km"])==180.0 &&
                parse(Float64,r["geodetic_lat_deg"])==90.0 &&
                parse(Float64,r["east_lon_deg"])==longitude_deg
        end
        @test length(matches)==1
        pole=only(matches)
        @test pole["status"]=="finite" && pole["expected_spatial_domain"]=="inside"
        expected_pole=Tuple(parse(Float64,pole[k]) for k in ("rho","T","wind_e","wind_n","wind_u"))
        @test isequal(q,expected_pole)
        indices=(findfirst(==(180.0),axes["alt_km"]), findfirst(==(90.0),axes["lat_deg"]),
                 findfirst(==(longitude_deg),axes["lon_deg"]))
        stored=Tuple(fields[k][indices...] for k in ("density_kgm3","temperature_K","wind_ew_ms","wind_ns_ms","wind_up_ms"))
        loaded=Tuple(getproperty(model.core.surrogate,k)[indices...] for k in (:rho,:T,:wind_e,:wind_n,:wind_u))
        @test isequal(loaded,stored)
        @test q[3]!=0 || q[4]!=0
    end
    thread_results=Vector{NTuple{5,Float64}}(undef,8*654);thread_ids=zeros(Int,length(thread_results))
    Threads.@threads :static for j in eachindex(thread_results)
        i=mod1(j,654)
        thread_results[j]=state_tuple(SpaceAGORA.getDensity(model,hs[i],aa[i],bb[i],times[i],isodd(j)))
        thread_ids[j]=Threads.threadid()
    end
    @test Threads.nthreads()>=4
    @test length(unique(thread_ids))>=2
    @test all(isequal(thread_results[j],expected[mod1(j,654)]) for j in eachindex(thread_results))
    @test grid_digests(model.core.surrogate)==before_fields
    copied=deepcopy(model)
    @test copied.core.surrogate.rho!==model.core.surrogate.rho
    copied.core.surrogate.rho[1]*=2
    @test grid_digests(model.core.surrogate)==before_fields
    copied=nothing
    mktempdir() do dir
        local_file=joinpath(dir,"fine_unmodified.jls");cp(FINE_FILE,local_file)
        transported=SpaceAGORA.GRAMGridAtmosphereModel(planet="mars",surrogate_file=local_file,expected_sha256=FINE_SHA)
        io=IOBuffer();serialize(io,transported);seekstart(io)
        restored=deserialize(io);rm(local_file)
        @test !isfile(local_file)
        @test restored.core.source_sha256==FINE_SHA
        @test restored.core.surrogate.rho!==transported.core.surrogate.rho
        @test all(isequal(state_tuple(SpaceAGORA.getDensity(restored,hs[i],aa[i],bb[i],times[i],true)),expected[i]) for i in eachindex(points))
    end
    @test bytes2hex(open(sha256,FINE_FILE))==source_before
    @test GRAMSuite._GRAM_WRAPPER[]===nothing && GRAMSuite._GRAM_WRAPPER_FILE[]==""
    @test !any(p->occursin("libgram",lowercase(p)),Libdl.dllist())
    for name in keys(INPUTS)
        @test bytes2hex(open(sha256,INPUTS[name])) == INPUT_SHA[name]
    end
end

@testset "Retained Odyssey input identities" begin
    for name in keys(INPUTS)
        @test occursin(r"^[0-9a-f]{64}$",INPUT_SHA[name])
        @test bytes2hex(open(sha256,INPUTS[name])) == INPUT_SHA[name]
    end
end
@testset "Retained Odyssey fine-grid normal-import regression" begin
    test_retained_odyssey_grid()
end
