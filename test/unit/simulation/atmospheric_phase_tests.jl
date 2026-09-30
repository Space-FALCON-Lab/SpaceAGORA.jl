module AtmosphericPhaseTests
using Test, SpaceAGORA, OrdinaryDiffEq, LinearAlgebra
const SM = SpaceAGORA.SimulationModel
const SC = SM.SimConfig
const SE = SpaceAGORA.SimulationEngine
const TV = SpaceAGORA.TelemetryVerification
const CB = SM.SimulationCallbacks
const planet = SM.make_no_gram_planet(:earth)
member(id; ν, rp=110_000.0, ra=400_000.0) = begin
    body = SM.Link(root=true, m=400.0, ref_area=4.0)
    ic = SM.InitialCondition(ra=planet.Rp_e + ra, rp=planet.Rp_e + rp, i=40.0, ω=0.0, Ω=0.0, ν=ν)
    SM.SpacecraftModel(SM.Joint[], [body], body, true, 400.0, 0.0, body.inertia, 0, 0, ic, id)
end
function cfg(members; T, tol=nothing, checkpoint=nothing, dir=tempdir(), resume=false)
    eff = (SM.InverseSquaredJ2GravityModel(), SM.AerodynamicCoefficientfM())
    args = TV.make_example_config(planet=planet, spacecraft=members[1], mission_time=T,
        initial_time=SM.InitialTime(year=2021, month=3, day=4, hour=5), dynamic_effectors=eff,
        density_model=SM.ExponentialAtmosphereModel(planet), ephemerides_model=SM.SimpleEphemeridesModel(),
        EI_km=300.0, results=false, verbose=false, results_directory=dir, solver_config=SM.SolverConfig(solver_mode=:tsit5))
    mc = SC.MissionConfiguration(mission_type=SM.MissionTime, keplerian=false, number_of_orbits=1, mission_time=T,
        orientation_sim=false, num_steps_to_save=1000, data_rate=5.0)
    args = SC._with_configuration(args; mission_configuration=mc, dynamics_model=SM.DynamicsModel(members, eff))
    tol === nothing || (args = SC._with_configuration(args; integration_tolerances=tol))
    if checkpoint !== nothing || resume
        ss = args.simulation_settings
        f = NamedTuple{fieldnames(typeof(ss))}(Tuple(getfield(ss, k) for k in fieldnames(typeof(ss))))
        args = SC._with_configuration(args; simulation_settings=SM.SimulationSettings(; merge(f, (results=true,
            checkpoint_enabled=checkpoint !== nothing, checkpoint_interval_s=checkpoint === nothing ? 200.0 : checkpoint, resume_from_checkpoint=resume, checkpoint_directory=joinpath(dir, "ckpt")))...))
    end
    args
end

const distinct = SM.IntegrationTolerances(reltol_orbit=1e-6, abstol_orbit=1e-8,
    reltol_atmosphere=1e-9, abstol_atmosphere=1e-11, dt_max_orbit=20.0, dt_max_atmosphere=0.2)
maxstep(sol) = maximum(diff(sol.t))
function trace_callback(records)
    record(i) = push!(records, (t=Float64(i.t), cap=i.opts.dtmax,
        rtol=i.opts.reltol.sc[1].pos[1], atol=i.opts.abstol.sc[1].pos[1],
        inside=copy(i.p.shared_buffers.in_atmosphere), active=copy(i.p.is_active)))
    DiscreteCallback((u,t,i)->true, record; save_positions=(false,false),
        initialize=(c,u,t,i)->record(i))
end
function check_phases(records)
    @test !isempty(records)
    @test all(records) do r
        inside=any(r.inside .& r.active)
        r.cap == (inside ? 0.2 : 20.0) &&
        r.rtol == (inside ? 1e-9 : 1e-6) && r.atol == (inside ? 1e-11 : 1e-8)
    end
end

@testset "Atmospheric phase lifecycle" begin
@testset "Atmospheric settings from initial and active state" begin
    cache=SE.SolverIntegratorCache(); records=NamedTuple[]
    a=cfg([member(1; ν=-35.0)]; T=10.0,tol=distinct)
    inside=run_simulation(a;solver_cache=cache,return_solution=true,extra_callbacks=(trace_callback(records),))
    @test first(records).t == 0.0
    check_phases(records)
    @test maxstep(inside) <= 0.2 + 1e-12
    @test cache.integrator.opts.reltol.sc[1].mass == distinct.reltol_mass
    @test cache.integrator.opts.abstol.sc[1].mass == distinct.abstol_mass

    # The atmospheric cap is effective from startup, including when changed
    # without touching the orbital cap. Default phase tolerances also apply.
    fields=NamedTuple{fieldnames(typeof(distinct))}(Tuple(getfield(distinct,k) for k in fieldnames(typeof(distinct))))
    finer=SM.IntegrationTolerances(;merge(fields,(dt_max_atmosphere=0.1,))...)
    finer_records=NamedTuple[]
    finer_sol=run_simulation(cfg([member(1;ν=-35.0)];T=10.0,tol=finer);
        return_solution=true,extra_callbacks=(trace_callback(finer_records),))
    @test maxstep(finer_sol) <= 0.1 + 1e-12
    @test collect(finer_sol.u[end]) != collect(inside.u[end])
    defaults=SM.IntegrationTolerances(); default_records=NamedTuple[]
    run_simulation(cfg([member(1;ν=-35.0)];T=1.0,tol=defaults);
        extra_callbacks=(trace_callback(default_records),))
    @test first(default_records).rtol == defaults.reltol_atmosphere
    @test first(default_records).atol == defaults.abstol_atmosphere

    # Entry and exit update both cap and componentwise tolerances.
    records=NamedTuple[]
    crossings=run_simulation(cfg([member(1;ν=-130.0,rp=200_000.0)];T=4000.0,tol=distinct);
        return_solution=true,extra_callbacks=(trace_callback(records),))
    check_phases(records)
    @test any(r->r.inside[1],records)
    @test !first(records).inside[1] && !last(records).inside[1]

    # The first member remains inside after the second member exits.
    records=NamedTuple[]
    mixed=run_simulation(cfg([member(1;ν=-130.0),member(2;ν=60.0)];T=1000.0,tol=distinct);
        return_solution=true,extra_callbacks=(trace_callback(records),))
    check_phases(records)
    @test first(records).inside == [false,true]
    @test any(r->r.inside == [true,false],records)
    @test maxstep(mixed) <= 0.2 + 1e-10

    records=NamedTuple[]
    same=run_simulation(cfg([member(1;ν=-130.0,rp=200_000.0),member(2;ν=-130.0,rp=200_000.0)];T=4000.0,tol=distinct);
        return_solution=true,extra_callbacks=(trace_callback(records),))
    check_phases(records)
    @test any(r->all(r.inside),records)
    @test all(!,last(records).inside)
    @test all(r->r.inside[1] == r.inside[2],records)

    # A real impact event removes the only atmospheric member from the aggregate.
    records=NamedTuple[]
    impacted=run_simulation(cfg([member(1;ν=-35.0),member(2;ν=0.0,rp=600_000.0,ra=600_000.0)];
        T=1500.0,tol=distinct);return_solution=true,extra_callbacks=(trace_callback(records),))
    check_phases(records)
    @test !last(records).active[1] && last(records).active[2]
    @test last(records).cap == 20.0
end

@testset "Integrator reuse resets phase before selecting initial step" begin
    cache=SE.SolverIntegratorCache()
    run_simulation(cfg([member(1;ν=-130.0)];T=500.0,tol=distinct);solver_cache=cache,return_solution=true)
    old=cache.integrator
    @test old.opts.dtmax == 0.2
    orbit()=cfg([member(1;ν=0.0,rp=600_000.0,ra=600_000.0)];T=500.0,tol=distinct)
    reused=run_simulation(orbit();solver_cache=cache,return_solution=true)
    fresh=run_simulation(orbit();return_solution=true)
    @test cache.integrator === old
    @test cache.integrator.opts.dtmax == 20.0
    @test reused.t == fresh.t
    @test reused.u == fresh.u

    # Reuse must refresh closures that capture the current caller's output.
    traced_cache=SE.SolverIntegratorCache(); before=NamedTuple[]; after=NamedTuple[]
    run_simulation(orbit();solver_cache=traced_cache,return_solution=true,
        extra_callbacks=(trace_callback(before),))
    traced_old=traced_cache.integrator; old_count=length(before)
    run_simulation(orbit();solver_cache=traced_cache,return_solution=true,
        extra_callbacks=(trace_callback(after),))
    @test traced_cache.integrator === traced_old
    @test length(before) == old_count
    @test !isempty(after) && first(after).t == 0.0
    # An event-layout change requires a fresh callback cache.
    additional=DiscreteCallback((u,t,i)->false, i->nothing)
    run_simulation(orbit();solver_cache=traced_cache,return_solution=true,
        extra_callbacks=(trace_callback(after),additional))
    @test traced_cache.integrator !== traced_old

    # A changed root condition with the same closure type must also rebuild,
    # because newer solver libraries retain conditions in bracketing caches.
    root_at(target, times)=ContinuousCallback((u,t,i)->t-target,
        i->push!(times,Float64(i.t));save_positions=(false,false))
    for n in (1,2)
        root_orbit()=cfg([member(i;ν=0.0,rp=600_000.0,ra=600_000.0) for i in 1:n];T=500.0,tol=distinct)
        root_cache=SE.SolverIntegratorCache(); early=Float64[]; later=Float64[]
        run_simulation(root_orbit();solver_cache=root_cache,return_solution=true,
            extra_callbacks=(root_at(200.0,early),))
        first_root_integrator=root_cache.integrator
        run_simulation(root_orbit();solver_cache=root_cache,return_solution=true,
            extra_callbacks=(root_at(300.0,later),))
        @test root_cache.integrator !== first_root_integrator
        @test only(early) ≈ 200.0
        @test only(later) ≈ 300.0
    end
end

@testset "Checkpoint segments and atmospheric resume" begin
    mktempdir() do dir
        records=NamedTuple[]
        a=cfg([member(1;ν=-130.0)];T=600.0,tol=distinct,checkpoint=200.0,dir=dir)
        final=run_simulation(a;return_solution=true,extra_callbacks=(trace_callback(records),))
        check_phases(records)
        @test first(final.t) == 400.0
        @test maxstep(final) <= 0.2 + 1e-10
        @test any(r->r.t == 400.0 && r.cap == 0.2,records)
        # Write a checkpoint inside the atmosphere but retain an orbital initial
        # condition, so a resume cannot silently reuse the original altitude flags.
        SE._write_checkpoint!(a,500.0,final(500.0),"tsit5")
        records=NamedTuple[]
        resumed=run_simulation(cfg([member(1;ν=-130.0)];T=600.0,tol=distinct,dir=dir,resume=true);
            return_solution=true,extra_callbacks=(trace_callback(records),))
        check_phases(records)
        @test first(records).t == 500.0
        @test first(records).inside == [true]
        @test maxstep(resumed) <= 0.2 + 1e-10
    end
end

@testset "Explicit step override and plain ODE compatibility" begin
    a=cfg([member(1;ν=-35.0)];T=2.0,tol=distinct)
    p=SM.ODEParams(n_sats=1,args=a);p.shared_buffers.in_atmosphere[1]=true
    prob=ODEProblem((du,u,p,t)->(du .= -u),[1.0],(0.0,2.0),p)
    cache=SE.SolverIntegratorCache();c=SE.SolverConfig(solver_mode=:tsit5)
    SE._solve_with_explicit_solver(prob,c,a,Tsit5(),1e-8,1e-8;dtmax_override=0.3,solver_cache=cache)
    @test cache.integrator.opts.dtmax == 0.3
    @test cache.integrator.opts.reltol == 1e-8
    @test cache.integrator.opts.abstol == 1e-8
    plain=ODEProblem((du,u,p,t)->(du .= -u),[1.0],(0.0,2.0))
    sol=SE._solve_with_explicit_solver(plain,c,a,Tsit5(),1e-8,1e-8)
    @test sol.u[end][1] ≈ exp(-2) rtol=1e-7
end
end
end
