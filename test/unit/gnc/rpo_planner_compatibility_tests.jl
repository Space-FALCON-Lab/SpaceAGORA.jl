# Prospective compatibility fixtures for the later HYPR adapter. This packet
# does not install that adapter or modify the legacy runtime route.
module RPOPlannerCompatibilityTests
using Test, SpaceAGORA, Random, StaticArrays, SHA, JSON
const S = SpaceAGORA.SimulationModel
const G = S.GuidanceHooks
const P = SpaceAGORA.RPOPlannerInterfaces

fingerprint(a::AbstractArray{Float64}) = (shape=size(a), sha256=bytes2hex(sha256(reinterpret(UInt8,vec(a)))))
fingerprint(x::Float64) = string(reinterpret(UInt64,x);base=16,pad=16)
function fixture_geometry()
    return S.RPOReferenceGeometry(S.RPOStationGeometry(
        reshape([0.,0.,50.],3,1);keepout_radius_m=0.25);
        chaser=S.RPOCubeSatGeometry(dims_m=(0.1,0.1,0.3)))
end
function fixture(mode, accel, seed)
    geometry=fixture_geometry()
    config=S.RPOPSOConfig(hypr_mode=mode, n_waypoints=2,n_particles=8,n_iters=3,
        adaptive_enable=false, adaptive_n_waypoints_min=2,adaptive_n_waypoints_max=2,
        adaptive_n_particles_min=8,adaptive_n_particles_max=8,
        adaptive_n_iters_min=3,adaptive_n_iters_max=3,
        sample_ds_m=0.2,safe_distance_m=0.2,retime_dt_s=0.1,
        retime_a_max_mps2=0.1,retime_max_speed_mps=0.5,
        retime_accel_limit_enable=accel,iteration_runtime_limit_s=Inf,
        rrt_warmstart_enable=true)
    result=S.rpo_pso_plan_path(SVector(3.,0.,0.),SVector(5.,1.,0.),geometry,config;
        safe_distance_m=0.2,rng=MersenneTwister(seed))
    times,positions,velocities=G.rpo_reference_from_path(result.path,geometry,result.config;safe_distance_m=0.2)
    return result,times,positions,velocities
end

snapshots=[]
@testset "Legacy seeded paths and effective retiming fixtures" begin
    @test P !== G
    @test :AbstractRPOPlanner in names(SpaceAGORA) # Packet 3 makes the accepted contract public.
    @test SpaceAGORA.AbstractRPOPlanner === SpaceAGORA.RPOPlannerInterfaces.AbstractRPOPlanner
    @test !hasfield(S.RPOGuidanceModel,:planner)
    @test fieldnames(S.RPOPlanBuffer)==(:valid,:plan,:updated_at_s)
    for (mode,accel,seed) in ((:legacy,false,741),(:legacy,true,742),(:manuscript,true,743))
        first=fixture(mode,accel,seed);second=fixture(mode,accel,seed)
        a=first[1];b=second[1]
        @test length(a.cost_history)>1 # Actual iterations, not a zero-waypoint shortcut.
        @test !a.iteration_timed_out
        @test fingerprint(a.path)==fingerprint(b.path)
        @test fingerprint(a.cost)==fingerprint(b.cost)
        @test fingerprint(a.cost_history)==fingerprint(b.cost_history)
        for i in 2:4
            @test fingerprint(first[i])==fingerprint(second[i])
            @test all(isfinite,first[i])
        end
        @test first[2][1]==0
        @test first[3][:,1]≈[3.,0.,0.]
        @test first[3][:,end]≈[5.,1.,0.]
        # Validate the unmodified reference against physical request limits.
        # Tangential retiming limits do not guarantee total-vector acceleration.
        req=P.RPOPlanningRequest(request_id=seed,chaser_id=101,target_id=201,epoch="fixture",
            time_s=0.,x_rtn=(3.,0.,0.,0.,0.,0.),target_state_ii=(7e6,0.,0.,0.,7500.,0.),
            goal_rtn_m=(5.,1.,0.),geometry=fixture_geometry(),geometry_revision="fixture-v1",
            constraints=P.RPOPlanningConstraints(clearance_m=0.2,max_speed_mps=0.5,
                max_acceleration_mps2=accel ? 0.1 : nothing),
            reference_dt_s=0.1,preview_horizon_steps=2,valid_until_s=first[2][end]+1.)
        ref=P.RPOReference(t_ref_s=first[2],r_ref_rtn_m=first[3],v_ref_rtn_mps=first[4],
            origin_time_s=0.,valid_until_s=req.valid_until_s,chaser_id=101,target_id=201,
            geometry_revision="fixture-v1")
        out=P.RPOPlanningResult(request_id=seed,status=:candidate,termination=:completed,reference=ref)
        checked=P.validate_rpo_result(req,out;clearance_at=(point,g)->S.rpo_clearance_distance_to_station(SVector{3}(point),g))
        if seed==742
            @test !checked.accepted && checked.reason===:acceleration_limit
            @test checked.metrics.max_acceleration_mps2 > 0.1*(1.0 + 1e-8)
        else
            @test checked.accepted
            @test checked.metrics.max_declared_speed_mps <= 0.5*(1.0 + 128eps(Float64))
            @test checked.metrics.max_implied_speed_mps <= 0.5*(1.0 + 128eps(Float64))
        end
        push!(snapshots,(mode=mode,acceleration_limited=accel,seed=seed,
            path=fingerprint(a.path),cost=fingerprint(a.cost),history=fingerprint(a.cost_history),
            t=fingerprint(first[2]),r=fingerprint(first[3]),v=fingerprint(first[4]),
            effective_config=Dict(string(k)=>string(getfield(a.config,k)) for k in fieldnames(typeof(a.config)))))
    end
end
snapshot_path=get(ENV,"SPACEAGORA_RPO_CONTRACT_SNAPSHOT","")
isempty(snapshot_path) || write(snapshot_path,JSON.json(snapshots))
end
