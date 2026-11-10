module ConstellationMPCTests
using Test, LinearAlgebra, SpaceAGORA
const CS=SpaceAGORA.ConstellationScheduling

@testset "Constellation laser MPC" begin
    c=walker27_constellation()
    times=collect(range(0.,.06c.period_s;length=19))
    @test_throws ArgumentError run_laser_mpc(c,times;horizon=4,execute=6)

    # The planner's per-interval terminal gain is what a single firing of that
    # link does in the nonlinear dynamics; the schedule is only as good as this.
    plan=plan_laser_schedule(c,times)
    e=findfirst(plan.available[:,1]);u=zeros(Int,size(plan.command));u[e,1]=1
    fired=laser_element_history(c,[last(times)],replay_laser_schedule(c,times,u))-
        laser_element_history(c,[last(times)],plan.reference)
    @test fired[1] ≈ plan.gains[1,e,1] rtol=1e-2

    # Executing the whole horizon at once is exactly the one-shot plan.
    whole=run_laser_mpc(c,times;horizon=18,execute=18)
    @test whole.command==plan.command
    @test length(whole.segments)==1

    r=run_laser_mpc(c,times;horizon=12,execute=6)
    @test [s.first_interval for s in r.segments]==[1,7,13]
    # MPC applies the leading commands of a plan from the current state, not a substitute.
    @test r.command[:,1:6]==plan_laser_schedule(c,times[1:13]).command[:,1:6]
    @test all(sum(r.command;dims=1).<=1)
    # Each segment's linear prediction is checked against the nonlinear replay of
    # the same commands; over short segments they must agree.
    for s in r.segments
        @test s.predicted[1] ≈ s.actual[1] rtol=1e-3
    end
    # Restarting each segment from the replayed state must not change the
    # trajectory: the chained replay ends where one continuous replay of the
    # executed schedule ends.
    y=replay_laser_schedule(c,times,r.command)(last(times))
    for s in axes(r.final_cartesian,2)
        x=ks_state_to_cartesian(view(y,CS.slice(s)))
        @test r.final_cartesian[1:3,s] ≈ x.position_ii_m atol=1e-3
        @test r.final_cartesian[4:6,s] ≈ x.velocity_ii_m atol=1e-6
    end
end
end
