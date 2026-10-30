module TerminalGainInterlinkEngineTests
using Test, SpaceAGORA, LinearAlgebra
const CS=SpaceAGORA.ConstellationScheduling
const KS=SpaceAGORA.SimulationModel.DynamicsKS

@testset "Fixed-interval schedule through the native interlink engine" begin
    case=walker27_constellation()
    # Deliberate emitter switch followed by coast, independent of optimization.
    times=[0.,10.,20.,30.]
    command=zeros(Int,26,3);command[1,1]=1;command[2,2]=1
    policy=TerminalGainPolicy(times,case.keys,command)
    args=constellation_replay_configuration(case,policy)
    native=run_simulation(args;return_solution=true,
        extra_callback_factories=(terminal_gain_schedule_callback,))
    replay=replay_laser_schedule(case,times,command)
    @test native.t[end]==last(times)
    history=native.prob.p.args.interlink_model.history
    @test all(t->any(h->h.time==t,history),times)
    for h in history
        expected=h.time<10 ? [case.keys[1]] : h.time<20 ? [case.keys[2]] : InterLinkKey[]
        @test h.active==expected
    end
    for s in 1:27
        x=KS.ks_state_to_cartesian(replay(last(times))[CS.slice(s)])
        @test norm(native.u[end].sc[s].pos-x.position_ii_m)<1e-3
        @test norm(native.u[end].sc[s].vel-x.velocity_ii_m)<1e-6
    end
    # Native state isolation must leave the planner's graph and configuration idle.
    @test !any(case.interlinks.graph.active)
    @test isempty(args.interlink_model.history)
end
end
