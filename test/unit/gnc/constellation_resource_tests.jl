module ConstellationResourceTests
# Battery/thermal/data-constrained terminal-gain scheduling. Coefficients come
# from the hardware mounted on the spacecraft (Battery, SolarArray, PowerLoad,
# LaserTerminal, RadioTerminal, DataStore, DataSource, ThermalNode, DataLink);
# the scheduling checks are ported from Constellation-Controllability's
# test/test_resource_constrained_scheduling.jl.
using Test, LinearAlgebra, SpaceAGORA, JuMP, HiGHS
const CS=SpaceAGORA.ConstellationScheduling
const SM=SpaceAGORA.SimulationModel

# The target and two helpers of the Walker-27 resource constellation.
function three(;kw...)
    c=walker27_resource_constellation()
    sub=merge(c,(spacecraft=c.spacecraft[1:3],edges=c.edges[1:2],data_keys=c.data_keys[1:2]))
    sub,resource_model(sub;terminal_bits=[.1e9,.15e9,.15e9],kw...)
end
node(;kw...)=ThermalNode(;merge((heat_capacity_j_k=150e3,conductance_w_k=150.,initial_k=293.,minimum_k=270.,
    maximum_k=310.,optical_absorptivity=.05,optical_reflectivity=.95),kw)...)

@testset "Constellation resource scheduling" begin
    sub,p=three();times=[0.,100.,200.,300.];edges=sub.edges
    env=(sun=ones(3,3),sink_K=fill(290.,3,3),ground=ones(Int,3))

    @testset "ThermalNode part" begin
        sc=walker27_constellation().spacecraft[1]
        @test !has_thermal_node(sc) && thermal_node(sc)===nothing
        n=add_thermal!(sc,sc.root,node())
        @test thermal_node(sc)===n && has_thermal_node(sc)
        @test_throws ArgumentError add_thermal!(sc,sc.root,node())
        @test_throws ArgumentError node(optical_absorptivity=.1)          # absorbed + reflected > incident
        @test_throws ArgumentError node(initial_k=320.)
        @test_throws ArgumentError node(conductance_w_k=0.)
    end

    @testset "Coefficients come from the mounted hardware" begin
        target,helper=sub.spacecraft[1],sub.spacecraft[2]
        b=SM.battery(helper);laser=only(e.actuator for e in mounted_actuators(helper) if e.actuator isa LaserTerminal)
        radio=only(e.actuator for e in mounted_actuators(target) if e.actuator isa RadioTerminal)
        @test p.capacity_kWh[2]==b.capacity_wh/1e3 && p.initial_kWh[2]==SM.initial_battery_wh(helper)/1e3
        @test p.reserve_kWh[2]≈.2b.capacity_wh/1e3
        @test p.bus_kW≈[SM.bus_load_w(m) for m in sub.spacecraft]./1e3
        @test p.bus_kW[1]≈.35  # the target's receiver electronics are a constant PowerLoad
        @test p.solar_kW[2]≈1.5  # Sun-tracking array at its converter limit
        @test p.laser_kW[1]≈laser.max_power_w/laser.efficiency/1e3
        @test p.radio_kW≈radio.max_power_w/radio.efficiency/1e3 && p.downlink_Gbps≈radio.max_rate_bps/1e9
        @test p.generation_Gbps[2]≈3e-6 && p.buffer_Gbit[1]≈15 && p.initial_Gbit[2]≈.6
        @test p.link_Gbps≈[.01,.01]
        # Emitter waste heat is the draw less the optical output, plus the absorbed return.
        @test p.laser_heat_emitter_kW[1]≈10/.35-10+.05*.95*10
        @test p.laser_heat_receiver_kW[1]≈.05*10
        # Changing a part changes the model.
        helper.root.batteries[1]=Battery(capacity_wh=6000.,initial_soc=.5,charge_efficiency=.9)
        q=resource_model(sub;terminal_bits=[.1e9,.15e9,.15e9])
        @test q.capacity_kWh[2]≈6 && q.initial_kWh[2]≈3 && q.charge_efficiency[2]≈.9
        empty!(helper.root.batteries)
        @test_throws ArgumentError resource_model(sub;terminal_bits=[.1e9,.15e9,.15e9])
        mismatch=walker27_resource_constellation(LaserSchedulingConfig(power_W=5000.))
        mismatch.spacecraft[2].root.laser_terminals  # mounted at the configured power: consistent
        @test resource_model(mismatch).laser_kW[1]≈5/.35
        @test_throws ArgumentError resource_model(merge(mismatch,(config=LaserSchedulingConfig(),)))
        @test_throws ArgumentError three(reserve_soc=1.5)
        @test_throws ArgumentError three(return_capture_fraction=-.1)
    end

    @testset "Supplied laser-only warm start" begin
        initial=[1 0 0;0 0 0];mask=trues(2,3)
        seed=CS.resource_seed(mask,edges,times,env,p;battery=false,thermal=true,data=false,initial_u=initial)
        @test seed.u==initial
        @test all(iszero,seed.v) && all(iszero,seed.flow) && all(iszero,seed.downlink)
        @test seed.states.temperature≈resource_replay(initial,zeros(3),zeros(2,3),zeros(3),edges,times,env,p).temperature
        result=optimize_resources(ones(2,3),mask,edges,times,env,p;battery=false,thermal=true,data=false,initial_u=initial)
        @test result.feasible && result.objective>=sum(initial)-1e-7
        for invalid in (fill(.5,2,3),ones(Int,2,3),zeros(1,3))
            @test_throws ArgumentError CS.resource_seed(mask,edges,times,env,p;
                battery=false,thermal=true,data=false,initial_u=invalid)
        end
        @test_throws ArgumentError CS.resource_seed(falses(2,3),edges,times,env,p;
            battery=false,thermal=true,data=false,initial_u=initial)
        @test_throws ArgumentError CS.resource_seed(mask,edges,times,env,p;
            battery=false,thermal=true,data=true,initial_u=initial)
        @test_throws ArgumentError optimize_resources(ones(2,3),mask,edges,times,env,merge(p,(maximum_K=fill(294.,3),));
            battery=false,thermal=true,data=false,initial_u=initial)
    end

    @testset "Verified resource warm starts" begin
        seed=CS.resource_seed(trues(2,3),edges,times,env,p;battery=true,thermal=true,data=true)
        @test seed!==nothing
        @test all(seed.states.queue[:,end].<=p.terminal_Gbit.+1e-8)
        @test all(sum(seed.u;dims=1).<=1)
        @test all(seed.states.energy.>=p.reserve_kWh.-1e-8)
        @test all(seed.states.temperature.<=p.maximum_K.+1e-8)
        for (mask,environment,parameters) in (
            (falses(2,3),env,p),
            (trues(2,3),merge(env,(ground=zeros(Int,3),)),p),
            (trues(2,3),env,merge(p,(maximum_K=fill(294.,3),))))
            @test CS.resource_seed(mask,edges,times,environment,parameters;battery=true,thermal=true,data=true)===nothing
        end
    end

    @testset "First-order resource scheduling" begin
        scores=[3. 3 3;2 2 2];mask=trues(2,3)
        free=optimize_resources(scores,mask,edges,times,env,p;battery=false,thermal=false,data=false)
        @test free.feasible && free.objective≈9
        @test all(free.u[1,:].==1) && all(sum(free.u;dims=1).<=1)
        # A battery just above its reserve, without sunlight, cannot power any laser.
        pb=merge(p,(initial_kWh=fill(1.5,3),reserve_kWh=fill(1.,3),solar_kW=zeros(3)))
        bat=optimize_resources(scores,mask,edges,times,env,pb;thermal=false,data=false)
        @test bat.feasible && sum(bat.u)==0
        pt=merge(p,(maximum_K=fill(305.,3),))
        hot=optimize_resources(scores,mask,edges,times,env,pt;battery=false,data=false)
        @test hot.feasible && sum(hot.u[1,:])<3
        @test maximum(hot.temperature)<=305+1e-6
        # Data that must be delivered forces links on even when they cost orbital gain.
        pd=merge(p,(initial_Gbit=[0.,1.,1.],generation_Gbps=zeros(3),
            terminal_Gbit=[0.,0.,0.],link_Gbps=[.01,.01],downlink_Gbps=.02))
        dat=optimize_resources([-3. -3 -3;-2 -2 -2],mask,edges,times,env,pd;battery=false,thermal=false)
        @test dat.feasible && sum(dat.u)>=2 && dat.objective<0
        @test 100sum(dat.downlink)≈2 atol=1e-6
        @test dat.downlink[1]≈0 atol=1e-8  # store-and-forward: nothing received yet
        @test all(100dat.flow[:,1].<=pd.initial_Gbit[2:3].+1e-6)
        @test maximum(abs,vec(sum(dat.queue;dims=1)).+[0.;cumsum(100dat.downlink)].-2)<1e-6
        impossible=optimize_resources(scores,falses(2,3),edges,times,env,pd;battery=false,thermal=false)
        @test !impossible.feasible && occursin("INFEASIBLE",impossible.status)
        allres=optimize_resources(scores,mask,edges,times,env,p;time_limit=30)
        @test allres.feasible
        replay=resource_replay(allres.u,allres.v,allres.flow,allres.downlink,edges,times,env,p;solar=allres.solar)
        @test replay.energy≈allres.energy atol=1e-7
        @test replay.temperature≈allres.temperature atol=1e-7
        @test replay.queue≈allres.queue atol=1e-7
        coast=resource_replay(zeros(Int,2,3),zeros(3),zeros(2,3),zeros(3),edges,times,env,p)
        equilibrium=290+p.bus_kW[2]/p.conductance_kW_K[2]
        expected=equilibrium+(p.initial_K[2]-equilibrium)*exp(-300p.conductance_kW_K[2]/p.heat_capacity_kJ_K[2])
        @test coast.temperature[2,end]≈expected
    end

    @testset "Battery efficiencies" begin
        # Discharging at efficiency 1/2 costs twice the energy.
        one_block=[0.,100.];env1=(sun=zeros(3,1),sink_K=fill(290.,3,1),ground=zeros(Int,1))
        lossy=merge(p,(discharge_efficiency=fill(.5,3),))
        u=reshape([1,0],2,1);z=zeros(Int,1);f=zeros(2,1);d=zeros(1)
        ideal=resource_replay(u,z,f,d,edges,one_block,env1,p);half=resource_replay(u,z,f,d,edges,one_block,env1,lossy)
        @test p.initial_kWh[2]-half.energy[2,2]≈2(p.initial_kWh[2]-ideal.energy[2,2])
        r=optimize_resources(ones(2,1),trues(2,1),edges,one_block,env1,lossy;thermal=false,data=false)
        @test r.feasible && sum(r.u)==1
        # Below unit efficiency the plan may waste energy (charging and discharging in one
        # block), never create it: the battery, collecting what fits, holds at least as much.
        physical=resource_replay(r.u,r.v,r.flow,r.downlink,edges,one_block,env1,lossy).energy
        @test all(physical.>=r.energy.-1e-7) && all(lossy.reserve_kWh.<=physical.<=lossy.capacity_kWh)
        # A charge limit below the surplus caps what the battery stores.
        sunny=(sun=ones(3,1),sink_K=fill(290.,3,1),ground=zeros(Int,1))
        capped=merge(p,(initial_kWh=fill(1.,3),max_charge_kW=fill(.5,3)))
        c=optimize_resources(zeros(2,1),trues(2,1),edges,one_block,sunny,capped;thermal=false,data=false)
        @test maximum(c.energy[:,2].-c.energy[:,1])<=.5*100/3600+1e-9
    end

    @testset "Reflected optical heating" begin
        times=[0.,100.];env=(sun=ones(3,1),sink_K=fill(290.,3,1),ground=zeros(Int,1))
        u=reshape([1,0],2,1);v=zeros(Int,1);flow=zeros(2,1);downlink=zeros(1)
        _,legacy=three(return_capture_fraction=0.)
        old=resource_replay(u,v,flow,downlink,edges,times,env,legacy)
        new=resource_replay(u,v,flow,downlink,edges,times,env,p)
        phi=exp(-.15*100/150)
        @test new.temperature[2,2]-old.temperature[2,2]≈(1-phi)*.475/.15
        @test new.temperature[[1,3],:]==old.temperature[[1,3],:]
        @test new.energy==old.energy && new.queue==old.queue
        _,half=three(return_capture_fraction=.5)
        @test resource_replay(u,v,flow,downlink,edges,times,env,half).temperature[2,2]-old.temperature[2,2]≈
            (new.temperature[2,2]-old.temperature[2,2])/2
        # No return heat without a reflecting target or an absorbing emitter.
        for (who,change) in ((1,(optical_reflectivity=0.,)),(2,(optical_absorptivity=0.,)))
            s,_=three();s.spacecraft[who].root.thermal_nodes[1]=node(;change...)
            q=resource_model(s;terminal_bits=[.1e9,.15e9,.15e9])
            @test resource_replay(u,v,flow,downlink,edges,times,env,q).temperature[2,:]≈old.temperature[2,:]
        end
        # A limit between the old and reflected one-block temperatures: only return heat forbids the link.
        bound=(maximum_K=fill((new.temperature[2,2]+old.temperature[2,2])/2,3),)
        mask=reshape([true,false],2,1);scores=ones(2,1)
        without=optimize_resources(scores,mask,edges,times,env,merge(legacy,bound);battery=false,data=false)
        with=optimize_resources(scores,mask,edges,times,env,merge(p,bound);battery=false,data=false)
        @test without.feasible && with.feasible
        @test sum(without.u)==1 && sum(with.u)==0
    end

    @testset "Walker-27 scheduler integration" begin
        c=walker27_resource_constellation();fine=collect(range(0.,.25c.period_s;length=61))
        @test_throws ArgumentError plan_resource_schedule(c,fine[1:60];stride=7)
        free=plan_resource_schedule(c,fine;battery=false,thermal=false,data=false)
        @test free.feasible && size(free.u)==(26,12) && length(free.times)==13
        # Without resource coupling the MILP must reproduce the independent
        # best-positive-link choice on the same block gains; otherwise it is not
        # the original scheduler's problem.
        @test free.u==signed_link_schedule(free.scores,free.available)
        @test free.predicted[1]≈sum(free.scores.*free.u) rtol=1e-10
        @test free.command==repeat(free.u;inner=(1,5))
        env=free.env
        @test all(in((0,1)),env.sun) && all(in((0,1)),env.ground)
        @test all((env.sink_K.==290).|((env.sink_K.==280).&(env.sun.==0)))  # a fully lit block is hot
        # The five-orbit delivery requirement cannot be met in a quarter orbit with one
        # ground-contact block, and must say so.
        @test !plan_resource_schedule(c,fine;reference=free.reference,battery=false,thermal=false,time_limit=30).feasible
        short=resource_model(c;terminal_bits=vcat(15e9,fill(.61e9,26)))
        res=plan_resource_schedule(c,fine;reference=free.reference,parameters=short,time_limit=60)
        @test res.feasible
        p=res.parameters
        replay=resource_replay(res.u,res.v,res.flow,res.downlink,c.edges,res.times,env,p;solar=res.solar)
        @test replay.energy≈res.energy atol=1e-6
        @test replay.temperature≈res.temperature atol=1e-6
        @test replay.queue≈res.queue atol=1e-6
        @test all(p.reserve_kWh.-1e-6 .<=replay.energy.<=p.capacity_kWh.+1e-6)
        @test all(p.minimum_K.-1e-6 .<=replay.temperature.<=p.maximum_K.+1e-6)
        @test all(-1e-6 .<=replay.queue.<=p.buffer_Gbit.+1e-6) && all(replay.queue[:,end].<=p.terminal_Gbit.+1e-6)
        @test all(res.u.<=res.available) && all(sum(res.u;dims=1).<=1)
        # Resources can only remove orbital authority from the free problem.
        @test res.objective<=free.objective+1e-9*abs(free.objective)
        delivered=dot(res.downlink,diff(res.times))
        @test sum(replay.queue[:,end])+delivered≈sum(p.initial_Gbit)+last(res.times)*sum(p.generation_Gbps) atol=1e-6
        # Mounting the hardware leaves the laser force model alone: same plan as the bare constellation.
        bare=plan_laser_schedule(walker27_constellation(),fine;reference=free.reference)
        @test bare.command==plan_laser_schedule(c,fine;reference=free.reference).command
    end
end
end
