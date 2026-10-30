module ConstellationSchedulingTests
using Test, LinearAlgebra, SpaceAGORA
const CS=SpaceAGORA.ConstellationScheduling
const KS=SpaceAGORA.SimulationModel.DynamicsKS
const OE=SpaceAGORA.OrbitalElementSensitivity

@testset "Constellation scheduling integration" begin
    c=walker27_constellation();p=c.config.planet;y=CS.initial_ks(c,0.)[1:10]
    x=KS.ks_state_to_cartesian(y);r=x.position_ii_m;v=x.velocity_ii_m
    @test y[9] ≈ -(dot(v,v)/2-p.μ/norm(r))
    m=KS.ks_physical_tangent(y,p)
    @test m.H*m.G ≈ Matrix{Float64}(I,6,6) atol=1e-9
    I3=Matrix{Float64}(I,3,3)
    direct=[zeros(3,3) I3;p.μ*(3r*r'/norm(r)^5-I3/norm(r)^3)+KS.ks_j2_acceleration_jacobian_si(r,p) zeros(3,3)]
    @test m.A_cartesian ≈ direct atol=1e-12 rtol=1e-10
    @test m.B_cartesian ≈ [zeros(3,3);I3] atol=1e-12
    # Check the lift's energy row and all five output derivatives numerically.
    C=OE.orbital_element_jacobian(r,v,p.μ);z=vcat(r,v)
    for j in 1:6
        step=j<=3 ? .1 : 1e-4
        d=zeros(6);d[j]=step
        function values(z)
            yy=KS.cartesian_to_ks_state(z[1:3],z[4:6],p)
            CS.elements(yy,c.config)
        end
        numerical=(values(z+d)-values(z-d))/(2step)
        @test C[:,j] ≈ numerical rtol=2e-4 atol=2e-6
    end
    @test signed_link_schedule([-1. 2. 0.;-2. 3. -1.],trues(2,3))==[0 0 0;0 1 0]
    @test signed_link_schedule(ones(2,1),trues(2,1))==reshape([1,0],2,1)
    @test signed_link_schedule(ones(2,1),falses(2,1))==zeros(Int,2,1)
    state=CS.initial_ks(c,0.);g=laser_geometry(state,c)
    @test all(x->x.available,g)
    off=zeros(26,1);on=copy(off);on[1,1]=1
    foff=similar(state);fon=similar(state)
    CS.constellation_rhs!(foff,state,(c,[0.,1.],off,deepcopy(c.interlinks)),0.)
    CS.constellation_rhs!(fon,state,(c,[0.,1.],on,deepcopy(c.interlinks)),0.)
    @test fon[10:10:end] ≈ ones(27)
    target_change=m.H*(fon[1:10]-foff[1:10])
    emitter_map=KS.ks_physical_tangent(state[11:20],p).H
    emitter_change=emitter_map*(fon[11:20]-foff[11:20])
    @test target_change[4:6] ≈ c.config.reflection_coefficient*c.config.power_W/
        (c.config.mass_kg*299792458.)*g[1].direction atol=1e-12
    @test target_change+emitter_change ≈ zeros(6) atol=1e-12
    state[11:20].=KS.cartesian_to_ks_state(-r,-v,p)
    @test !laser_geometry(state,walker27_constellation(LaserSchedulingConfig(max_range_m=2e7)))[1].available
    F,V=CS.zoh([0. 1.;0. 0.],reshape([0.,1.],2,1),2.)
    @test F≈[1. 2.;0. 1.]
    @test V≈reshape([2.,2.],2,1)
    times=collect(range(0.,.02c.period_s;length=6))
    result=plan_laser_schedule(c,times)
    @test all(sum(result.command;dims=1).<=1)
    @test all(result.command.<=result.available)
    @test sum(result.scores.*result.command)≈result.predicted[1,end]
    @test sum(result.scores.*result.command)≈sum(result.scores.*signed_link_schedule(result.scores,result.available))
    @test result.policy isa AbstractSchedulingPolicy
    @test all(link->link isa LaserBeamLink,c.interlinks.graph.edge_properties)
    @test c.keys[1]==((1,1),(2,1)) # native undirected key, not emitter/receiver order
    @test all(count(==(1),column)<=1 for column in eachcol(result.command))
    native=CS.interlink_state(CS.initial_ks(c,0.),c)
    graph=deepcopy(c.interlinks)
    scheduler=SchedulingPolicyModel(result.policy)
    schedule_interlinks!(graph,scheduler,c.spacecraft,native,p.μ;context=CS.link_context(c,0.))
    @test graph.graph.active==Bool.(result.command[:,1])
    selected=findfirst(!iszero,result.command[:,1])
    push!(graph.forbidden_pairs,minmax(c.edges[selected]...))
    schedule_interlinks!(graph,scheduler,c.spacecraft,native,p.μ;context=CS.link_context(c,0.))
    @test !any(graph.graph.active) # do not substitute another helper during replay
    empty!(graph.forbidden_pairs)
    schedule_interlinks!(graph,scheduler,c.spacecraft,native,p.μ;context=CS.link_context(c,-1.))
    @test !any(graph.graph.active)
    # Zero or negative authority must select no link; no fallback after the horizon.
    for score in (0.,-1.)
        scheduler=SchedulingPolicyModel(TerminalGainPolicy([0.,1.],c.keys,fill(score,26,1)))
        schedule_interlinks!(graph,scheduler,c.spacecraft,native,p.μ;context=CS.link_context(c,0.))
        @test !any(graph.graph.active)
    end
    schedule_interlinks!(graph,SchedulingPolicyModel(result.policy),c.spacecraft,native,p.μ;
        context=CS.link_context(c,last(times)))
    @test !any(graph.graph.active)
    @test_throws DimensionMismatch TerminalGainPolicy([0.,1.],c.keys,zeros(25,1))
    @test_throws ArgumentError plan_laser_schedule(c,times;objective=:unknown)
    @test_throws ArgumentError replay_laser_schedule(c,times,fill(.5,26,5))
end
end # module
