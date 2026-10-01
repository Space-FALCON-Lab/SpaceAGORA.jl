# Prospective compatibility fixtures for the later HYPR adapter. This packet
# does not install that adapter or modify the legacy runtime route.
module RPOPlannerCompatibilityTests
using Test, SpaceAGORA, Random, StaticArrays, SHA, JSON
const S = SpaceAGORA.SimulationModel
const G = S.GuidanceHooks
const P = SpaceAGORA.RPOPlannerInterfaces

fingerprint(a::AbstractArray{Float64}) = (shape=size(a), sha256=bytes2hex(sha256(reinterpret(UInt8,vec(a)))))
fingerprint(x::Float64) = string(reinterpret(UInt64,x);base=16,pad=16)
function fixture(mode, accel, seed)
    geometry=S.RPOReferenceGeometry(S.RPOStationGeometry(
        reshape([0.,0.,50.],3,1);keepout_radius_m=0.25);
        chaser=S.RPOCubeSatGeometry(dims_m=(0.1,0.1,0.3)))
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
    @test !(:AbstractRPOPlanner in names(SpaceAGORA))
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
        push!(snapshots,(mode=mode,acceleration_limited=accel,seed=seed,
            path=fingerprint(a.path),cost=fingerprint(a.cost),history=fingerprint(a.cost_history),
            t=fingerprint(first[2]),r=fingerprint(first[3]),v=fingerprint(first[4]),
            effective_config=Dict(string(k)=>string(getfield(a.config,k)) for k in fieldnames(typeof(a.config)))))
    end
end
snapshot_path=get(ENV,"SPACEAGORA_RPO_CONTRACT_SNAPSHOT","")
isempty(snapshot_path) || write(snapshot_path,JSON.json(snapshots))
end
