module CheckpointOutputGridTests
# The run-clock table saver: canonical output epochs indexed from the run origin,
# checked against an affine trajectory whose exact solution is known.
using SpaceAGORA, Test, LinearAlgebra
const S=SpaceAGORA.SimulationModel
const SC=S.SimulationCallbacks
const SE=SpaceAGORA.SimulationEngine
const SEG=SC._segment_output_times
truth(t)=[1.0 + 2t,5.0 - 3t]
rhs!(dy,y,p,t)=(dy[1]=2.;dy[2]=-3.;nothing)
function clockparams(origin,stop;horizon=1.04)
    return (shared_buffers=(periodic_clock_origin_s=Ref(origin),periodic_clock_horizon_s=Ref(horizon),
        periodic_clock_pending_s=Ref(NaN),solve_segment_end_time=Ref(stop)),)
end
@testset "Mission grid exact boundary and floating cadence" begin
    for (origin,start,stop,rate) in ((0.,0.,1.,.1),(0.,.37,.81,.1),
            (0.,prevfloat(.3),.5,.1),(0.,nextfloat(.3),.5,.1),
            (0.,0.,nextfloat(.3),.1),(2.5,2.65,3.17,.2),(0.,.37,.37,.1))
        points=SEG(origin,start,stop,rate;horizon=max(origin+1.04,stop))
        @test first(points)==start && last(points)==stop
        @test all(x->start<=x<=stop,points)
        @test all(diff(points).>0.)
    end
    @test .3 in SEG(0.,0.,nextfloat(.3),.1;horizon=1.04)
    @test prevfloat(.3) in SEG(0.,prevfloat(.3),.5,.1;horizon=1.04)
    @test .3 in SEG(0.,prevfloat(.3),.5,.1;horizon=1.04)
    @test !(.3 in SEG(0.,nextfloat(.3),.5,.1;horizon=1.04))
    @test SEG(0.,.37,.81,.1)==[.37,.4,.5,.6,.7,.8,.81]
    @test SEG(2.5,2.65,3.17,.2)==[2.65,2.7,2.9,3.1,3.17]
    @test_throws ArgumentError SEG(0.,0.,1.,.1;horizon=.9)
    @test_throws ArgumentError SEG(0.,0.,1.,.1;horizon=NaN)
    @test SEG(0.,0.,nextfloat(.3),.1;horizon=nextfloat(.3))==[0.,.1,.2,.3,nextfloat(.3)]
    @test SEG(0.,0.,.3,.1;horizon=.3)==[0.,.1,.2,.3]
    @test SEG(2.5,2.5,2.8,.1)==[2.5,2.6,2.7,2.8]
    @test SEG(0.,1e15,1e15+4.,1.)==[1e15,1e15+1.,1e15+2.,1e15+3.,1e15+4.]
    @test SEG(0.,1e16,1e16+8.,2.)==[1e16,1e16+2.,1e16+4.,1e16+6.,1e16+8.]
    @test_throws ArgumentError SEG(1.,1e16,1e16+8.,2.)
    @test_throws ArgumentError SEG(1e16,1e16,1e16+8.,1.)
    @test_throws ArgumentError SEG(0.,0.,1e300,1e-300)
    @test_throws ArgumentError SEG(-floatmax(Float64),0.,floatmax(Float64),floatmax(Float64))
    tiny=nextfloat(0.)
    @test SEG(0.,0.,2tiny,tiny)==[0.,tiny,2tiny]
    for rate in (0.,-1.,NaN,Inf)
        @test_throws ArgumentError SEG(0.,0.,1.,rate)
    end
    for triple in ((NaN,0.,1.),(0.,NaN,1.),(0.,0.,Inf),(1.,0.,2.),(0.,2.,1.))
        @test_throws ArgumentError SEG(triple...,.1)
    end
end
function make_saver(params,active_start,rate)
    values=SC.SavedValues(Float64,Vector{Float64})
    observations=Tuple{Float64,Float64,Float64}[]
    function saver(u,t,int)
        @test t>=active_start[]
        if params!==nothing
            @test t<=params.shared_buffers.solve_segment_end_time[]
        end
        @test norm(u-truth(t))<=2e-11
        push!(observations,(t,int.tprev,int.t))
        return copy(u)
    end
    return SC._run_clock_saving_callback(saver,values,rate),values,observations
end
function run_fresh(start,stop;origin=0.,ordinary=false)
    p=ordinary ? nothing : clockparams(origin,stop)
    cb,values,observations=make_saver(p,Ref(start),.1)
    prob=SE.ODEProblem(rhs!,truth(start),(start,stop),p)
    sol=SE.solve(prob,SE.DP8();callback=cb,reltol=1e-11,abstol=1e-12,
        dtmax=.17,save_everystep=false,dense=false)
    @test sol.t[end]==stop
    return copy(values.t),deepcopy(values.saveval),observations
end
@testset "Normal ordinary affine trajectory and fresh resumed clock" begin
    t,y,_=run_fresh(0.,1.04)
    @test t==vcat(collect(0.:.1:1.),1.04)
    t,y,_=run_fresh(.37,1.04)
    @test t==[.37,.4,.5,.6,.7,.8,.9,1.,1.04]
    @test all(t.>=.37)
    t,y,_=run_fresh(prevfloat(.3),.75)
    @test t==[prevfloat(.3),.3,.4,.5,.6,.7,.75]
    ordinary,_,_=run_fresh(.37,1.04;ordinary=true)
    @test ordinary==vcat(collect(range(.37,1.04;step=.1)),1.04)
end
@testset "Cached reinit preserves original grid and real segment endpoints" begin
    p=clockparams(0.,.37);active_start=Ref(0.)
    cb,values,observations=make_saver(p,active_start,.1)
    prob=SE.ODEProblem(rhs!,truth(0.),(0.,.37),p)
    int=SE.DiffEqBase.init(prob,SE.DP8();callback=cb,reltol=1e-11,abstol=1e-12,
        dtmax=.17,save_everystep=false,dense=false)
    combined=Float64[]
    for (j,(start,stop)) in enumerate(((0.,.37),(.37,.81),(.81,1.04)))
        if j>1
            empty!(values.t);empty!(values.saveval)
            active_start[]=start;p.shared_buffers.solve_segment_end_time[]=stop
            SE.SciMLBase.reinit!(int,truth(start);t0=start,tf=stop,erase_sol=true,
                reset_dt=true,reinit_callbacks=true)
        end
        SE.DiffEqBase.solve!(int)
        @test int.t==stop
        @test first(values.t)==start && last(values.t)==stop
        @test values.t==SEG(0.,start,stop,.1;horizon=1.04)
        append!(combined,j==1 ? values.t : values.t[2:end])
    end
    @test all(diff(combined).>0.)
    @test combined==[0.,.1,.2,.3,.37,.4,.5,.6,.7,.8,.81,.9,1.,1.04]
    @test all(x->x in combined,collect(0.:.1:1.))
end
@testset "Canonical output survives extended horizon with cached problem bounds" begin
    short=nextfloat(.3)
    p=clockparams(0.,short;horizon=short);active_start=Ref(0.)
    cb,values,obs=make_saver(p,active_start,.1)
    int=SE.DiffEqBase.init(SE.ODEProblem(rhs!,truth(0.),(0.,short),p),SE.DP8();
        callback=cb,reltol=1e-11,abstol=1e-12,dtmax=.17,save_everystep=false,dense=false)
    SE.DiffEqBase.solve!(int)
    @test values.t==[0.,.1,.2,.3,short]
    retained=copy(values.t)
    active_start[]=short;p.shared_buffers.solve_segment_end_time[]=1.04
    p.shared_buffers.periodic_clock_horizon_s[]=1.04
    empty!(values.t);empty!(values.saveval)
    SE.SciMLBase.reinit!(int,truth(short);t0=short,tf=1.04,erase_sol=true,
        reset_dt=true,reinit_callbacks=true)
    SE.DiffEqBase.solve!(int)
    @test values.t==[short,.4,.5,.6,.7,.8,.9,1.,1.04]
    append!(retained,values.t[2:end])
    @test retained==[0.,.1,.2,.3,short,.4,.5,.6,.7,.8,.9,1.,1.04]
    @test all(diff(retained).>0.)
    @test all(t->t in retained,[0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.])
end
end # module
