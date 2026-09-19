using Test, Serialization, SHA, StaticArrays, LinearAlgebra
println("Adapter test source SHA256: ",bytes2hex(sha256(read(@__FILE__))))
println("Julia threads available: ",Threads.nthreads())
flush(stdout)

const EM = SpaceAGORA.SimulationModel.EnvironmentModels
const CB = SpaceAGORA.SimulationModel.SimulationCallbacks
const EN = SpaceAGORA.SimulationEngine
const AT = SpaceAGORA.SimulationModel.AbstractTypes
const ES = SpaceAGORA.SimulationModel.EffectorSampling

# This synthetic spatial field is affine within every non-seam cell. Its wide
# altitude range deliberately distinguishes the adapter from native-specific
# clamps, polynomial fallbacks and the old 2000 km vacuum rule.
fields(h, a, b) = (1e-8*(1+h/1e6+a/100+b/1000), 180+h/1e5+a/10+b/100, 10.0, 20.0, 3.0)
function fixture(path)
    hs, lats, lons = [-0.1, 0.0, 200.0, 2500.0], [-30.0, 0.0, 30.0], [0.0, 80.0, 210.0, 300.0]
    keys = ("density_kgm3", "temperature_K", "wind_ew_ms", "wind_ns_ms", "wind_up_ms")
    payload = Dict{String,Any}(
        "status"=>"ok", "type"=>"surrogate_trilinear", "planet"=>"mars",
        "grid"=>Dict("alt_km"=>hs, "lat_deg"=>lats, "lon_deg"=>lons),
        "fields"=>Dict(k=>[fields(h*1000,a,b)[i] for h in hs, a in lats, b in lons] for (i,k) in enumerate(keys)),
        "provenance"=>"Synthetic interface test only, not validated Mars atmosphere",
        "epoch"=>"synthetic frozen epoch", "datum"=>"synthetic interface fixture",
    )
    serialize(path, payload)
    return path
end
function check_state(state, expected)
    values = (state[1],state[2],Tuple(state[3])...)
    @test length(values) == 5
    for (value, want) in zip(values,expected)
        @test isapprox(value,want; rtol=2e-13,atol=1e-22)
    end
end
struct PoisonParams end
Base.getproperty(::PoisonParams, name::Symbol) = error("Unexpected parameter access: $name")

function test_params(model; freeze=false, vacuum=false)
    cfg = withenv("SPACEAGORA_DENSITY_FREEZE_PER_STEP"=>(freeze ? "1" : "0"),
                  "SPACEAGORA_VACUUM_GRAM_CACHE"=>(vacuum ? "1" : "0"),
                  "SPACEAGORA_GRAM_TRACK_CACHE"=>"0",
                  "SPACEAGORA_GRAM_RUNTIME_STATS"=>"0") do
        CB._snapshot_callback_env_config()
    end
    return (
        args=(environment_model=(density_model=model,planet=(T_ref=190.0,Rp_e=3396190.0,Rp_p=3396190.0,ω=SVector(0.0,0.0,0.0)),EI=100.0,
                ephemerides_model=SpaceAGORA.SimpleEphemeridesModel(prime_meridian_at_reference_rad=0.0)),
              dynamics_model=(spacecraft=fill(nothing,4),dynamic_effectors=()),
              mission_configuration=(keplerian=false,mission_time=100.0)),
        shared_buffers=(density_models=AT.AbstractDensityModel[],
            callback_env_config=Ref(cfg),
            gram_density_cache=Union{Nothing,CB.GramTrackCache}[nothing],
            densities=[999.0], temperatures=[999.0], winds=[SVector(999.0,999.0,999.0)],
            density_sample_t=[0.0], current_time=Ref(0.0),in_atmosphere=[true],
            et_start=Ref(0.0),planet_frame_ephemeris_cache=Ref(nothing),
            spice_runtime_counters=(planet_pxform_runtime_calls=Threads.Atomic{Int64}(0),)),
    )
end
function cartesian_state(h,a,b)
    lat,lon=deg2rad(a),deg2rad(b)
    radius=3396190.0+h
    return SVector(radius*cos(lat)*cos(lon),radius*cos(lat)*sin(lon),radius*sin(lat),100.0,200.0,300.0,100.0)
end
function frame(h,a,b)
    ES.PlanetFrameSample(SMatrix{3,3,Float64}(I),SVector(1.0,2.0,3.0),SVector(4.0,5.0,6.0),h,deg2rad(a),deg2rad(b))
end

@testset "SpaceAGORA native-free grid adapter" begin
    @test GRAMSuite._GRAM_WRAPPER[] === nothing
    mktempdir() do dir
        path = fixture(joinpath(dir,"mars_surrogate.jls"))
        digest = bytes2hex(sha256(read(path)))
        model = EM.GRAMGridAtmosphereModel(planet="mars",surrogate_file=path,expected_sha256=digest)
        @test model isa AT.AbstractDensityModel
        @test model.core isa GRAMSuite.GRAMGridAtmosphereModel
        @test model.core.source_sha256 == digest
        # The current-main epoch hook must preserve a deliberately frozen grid;
        # changing the simulation epoch does not regenerate its stored fields.
        @test SpaceAGORA.with_density_model_epoch(model,
            SpaceAGORA.SimulationModel.InitialTime(year=2001,month=11,day=7)) === model
        @testset "Direct scalar contract, units, frozen time and no native fallbacks" begin
            for h in (-75.0,125000.0,2400000.0), t in (-50.0,0.0,10000.0), wind in (false,true)
                expected = fields(h,20.0,40.0)
                check_state(EM.getDensity(model,h,deg2rad(20.0),deg2rad(40.0),t,wind),expected)
                check_state(EM.getDensity(model,h,deg2rad(20.0),deg2rad(40.0),t,wind,PoisonParams()),expected)
            end
            for h in (-101.0,2500001.0)
                @test_throws DomainError EM.getDensity(model,h,0.0,0.0,0.0,true,PoisonParams())
            end
            for a in (-30.01,30.01)
                @test_throws DomainError EM.getDensity(model,1000.0,deg2rad(a),0.0,0.0,true)
            end
            for i in 1:4, bad in (NaN,Inf,-Inf)
                q=[1000.0,0.0,0.0,0.0]; q[i]=bad
                @test_throws DomainError EM.getDensity(model,q...,true)
            end
            for b in (330.0,-30.0,690.0)
                expected = (fields(1000.0,20.0,300.0) .+ fields(1000.0,20.0,0.0)) ./ 2
                check_state(EM.getDensity(model,1000.0,deg2rad(20.0),deg2rad(b),0.0,true),expected)
            end
            vacuum = EM.GRAMGridAtmosphereModel(planet="mars",surrogate_file=path,above_grid=:vacuum,vacuum_temperature=175.0)
            check_state(EM.getDensity(vacuum,2500001.0,0.0,0.0,0.0,true),(0.0,175.0,0.0,0.0,0.0))
            @test_throws DomainError EM.getDensity(vacuum,2500001.0,deg2rad(31.0),0.0,0.0,true)
            @test !EM.density_vanishes_above_entry_interface(model)
        end
        @testset "Existing generic batch is scalar-equivalent" begin
            hs=[-75.0,125000.0,2400000.0]; aa=deg2rad.([-20.0,0.0,20.0]); bb=deg2rad.([40.0,140.0,240.0])
            rhos=zeros(3); temps=zeros(3); winds=fill(SVector(0.0,0.0,0.0),3)
            for times in (0.0,[-20.0,0.0,10000.0]), wind in (false,true)
                @test EM.getDensityBatch!(rhos,temps,winds,model,hs,aa,bb,times,wind,PoisonParams()) === nothing
                for i in eachindex(hs)
                    check_state((rhos[i],temps[i],winds[i]),fields(hs[i],rad2deg(aa[i]),rad2deg(bb[i])))
                end
            end
            @test_throws ArgumentError EM.getDensityBatch!(zeros(2),temps,winds,model,hs,aa,bb,0.0,true,nothing)
            @test_throws ArgumentError EM.getDensityBatch!(rhos,temps,winds,model,hs,aa,bb,[0.0],true,nothing)
            @test EM.getDensityBatch!(Float64[],Float64[],SVector{3,Float64}[],model,Float64[],Float64[],Float64[],0.0,true,nothing) === nothing
        end
        @testset "Satellite ownership and native exclusion" begin
            p=test_params(model)
            withenv("SPACEAGORA_GRAM_PER_SAT_INSTANCES"=>"1") do
                @test EN._initialize_density_model_instances!(p) === nothing
            end
            @test isempty(p.shared_buffers.density_models)
            for i in 1:4
                @test CB._density_model_for_sat(p,i) === model
            end
            @test CB._density_batch_model_for_callback(p,4) === model
            @test CB.density_model_threadsafe(model)
            @test !CB._is_gram_density_model(model)
            @test !CB._gram_track_trajectory_supported(model)
            @test CB._gram_isolated_pool_batch_model_for_callback(p,4) === nothing
        end
        @testset "Stage coordinates cannot be masked by buffers or optional approximations" begin
            x=SVector(1.0,2.0,3.0,4.0,5.0,6.0,100.0)
            for freeze in (false,true), vacuum in (false,true)
                p=test_params(model;freeze,vacuum)
                @test CB._callback_env_config(p).density_freeze_per_step == freeze
                @test CB._callback_env_config(p).vacuum_gram_cache_enabled == vacuum
                # A stale equal-time sample must neither be read nor overwritten
                # by the grid-specific read-only buffered consumer.
                for elapsed in (0.0,50.0)
                    s=EN.sample_buffered_atmosphere(cartesian_state(125000.0,20.0,40.0),p,1,elapsed)
                    check_state((s.rho_kg_m3,s.temperature_k,s.wind_pp),fields(125000.0,20.0,40.0))
                    @test_throws DomainError EN.sample_buffered_atmosphere(cartesian_state(125000.0,31.0,40.0),p,1,elapsed)
                    @test p.shared_buffers.densities == [999.0]
                end
                for write_buffers in (false,true)
                    s=EN._sample_atmosphere_from_planet_frame(x,frame(125000.0,20.0,40.0),p,1,0.0;write_buffers)
                    check_state((s.rho_kg_m3,s.temperature_k,s.wind_pp),fields(125000.0,20.0,40.0))
                    @test_throws DomainError EN._sample_atmosphere_from_planet_frame(x,frame(125000.0,31.0,40.0),p,1,0.0;write_buffers)
                    @test_throws DomainError EN._sample_atmosphere_from_planet_frame(x,frame(2500001.0,20.0,40.0),p,1,0.0;write_buffers)
                end
                @test p.shared_buffers.winds[1] == SVector(10.0,20.0,3.0)
                @test p.shared_buffers.density_sample_t[1] == 0.0
            end
            # The pre-existing freeze/buffer policy still applies to other models.
            plain=EM.ExponentialAtmosphereModel(1e-8,0.0,10000.0;temperature_k=180.0)
            p=test_params(plain;freeze=true,vacuum=false)
            @test EN._buffered_atmosphere_valid(p,1,50.0)
            s=EN._sample_atmosphere_from_planet_frame(x,frame(125000.0,20.0,40.0),p,1,50.0)
            @test s.rho_kg_m3 == 999.0
        end
        @testset "Stored ENU wind receives exactly one downstream rotation" begin
            lat,lon=deg2rad(20.0),deg2rad(40.0)
            rho,temp,wind=EM.getDensity(model,125000.0,lat,lon,0.0,true)
            angle=0.37
            lpi=@SMatrix [cos(angle) sin(angle) 0.0; -sin(angle) cos(angle) 0.0; 0.0 0.0 1.0]
            pos=SVector(2e6,1e6,3e6); vel=SVector(1000.0,2500.0,-300.0)
            pf=ES.PlanetFrameSample(lpi,pos,vel,125000.0,lat,lon)
            spacecraft=(links=((root=true,ref_area=2.0),),)
            x=ES.StateSample(lpi'*pos,lpi'*vel,100.0;spacecraft)
            env=ES.EnvironmentSample((γ=1.3,R=190.0);planet_frame=pf,atmosphere=ES.AtmosphereSample(rho,temp,wind))
            ae=SpaceAGORA.SimulationModel.DynamicEffectors.AerodynamicEffectors
            force,torque,drag,lift,crossforce=ae._aero_pure_wrench(:constant,x,env)
            east=SVector(-sin(lon),cos(lon),0.0)
            north=SVector(-sin(lat)*cos(lon),-sin(lat)*sin(lon),cos(lat))
            up=SVector(cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat))
            vrel=vel-(10east+20north+3up)
            expected=lpi'*(-0.5*rho*2.2*2.0*norm(vrel)*vrel)
            @test force ≈ expected rtol=2e-13
            @test drag ≈ expected rtol=2e-13
            @test torque == lift == crossforce == SVector(0.0,0.0,0.0)
            # Deliberately incorrect direct-planet-frame interpretation differs.
            wrong=vel-wind
            @test !isapprox(force,lpi'*(-0.5*rho*2.2*2.0*norm(wrong)*wrong);rtol=1e-5)
        end
        @testset "Independent copy, in-memory serialization and threaded read-only queries" begin
            duplicate=deepcopy(model)
            @test duplicate.core.surrogate.rho !== model.core.surrogate.rho
            @test duplicate.core.metadata !== model.core.metadata
            duplicate.core.surrogate.rho[1,1,1] *= 2
            duplicate.core.metadata["epoch"]="changed copy only"
            @test duplicate.core.surrogate.rho[1,1,1] != model.core.surrogate.rho[1,1,1]
            @test model.core.metadata["epoch"] == "synthetic frozen epoch"
            io=IOBuffer(); serialize(io,model); seekstart(io); restored=deserialize(io)
            @test restored.core.surrogate.rho !== model.core.surrogate.rho
            rm(path)
            check_state(EM.getDensity(restored,125000.0,deg2rad(20.0),deg2rad(40.0),0.0,true),fields(125000.0,20.0,40.0))
            qs=[(Float64(1000+i),deg2rad(Float64(mod(i,41)-20)),deg2rad(Float64(mod(i,360))),Float64(i)) for i in 1:10000]
            serial=[EM.getDensity(model,q...,true) for q in qs]
            threaded=similar(serial)
            Threads.@threads for i in eachindex(qs)
                threaded[i]=EM.getDensity(model,qs[i]...,true)
            end
            @test threaded == serial
            @test GRAMSuite._GRAM_WRAPPER[] === nothing
        end
    end
    @test GRAMSuite._GRAM_WRAPPER[] === nothing
end
