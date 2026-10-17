# The burn-plan family and the generic algorithm binding (`AlgorithmGuidance`):
# trusted validation reasons, the planner, runs through the shared request driver,
# and the burn scheduler's recorded refusals.
module BurnPlanFamilyTests
using Test, Random, SpaceAGORA, SpaceAGORA.SimulationModel
using OrdinaryDiffEq: DiscreteCallback
const GA = SpaceAGORA.GNCAlgorithms

const PLANET = SimulationModel.make_no_gram_planet(:earth)

# A second family, written only against the general contract (the last testset):
# it shows that the binding carries nothing burn-specific, with a configured
# link, the default `:fail` policy and an algorithm that throws.
struct CounterFamily <: GA.AbstractAlgorithmFamily end
struct CountRequest <: GA.AbstractAlgorithmRequest
    request_id::Int
    time_s::Float64
    reason::Symbol
end
struct Counter <: GA.AbstractGNCAlgorithm
    fail_after::Int
end
GA.algorithm_family(::Counter) = CounterFamily()
GA.family_request(::CounterFamily, settings, u, p, t, i, reason, request_id) =
    CountRequest(request_id, t, reason)
GA.propose!(state, c::Counter, r::CountRequest, rng) = r.request_id > c.fail_after ? error("counter stopped") :
    GA.AlgorithmResult(request_id=r.request_id, status=:candidate, termination=:counted, product=r.request_id)
GA.validate_proposal(::CounterFamily, r::CountRequest, result, trusted) =
    GA.ValidationResult(result.product == r.request_id, :checked)
GA.install!(link::Base.RefValue{Int}, ::CounterFamily, product, r, t; validation, diagnostics) =
    (link[] = product; nothing)

function burn_sat(; prop_mass=50.0, id=7)
    root = Link(root=true, m=100.0, ref_area=2.0)
    ic = InitialCondition(ra=PLANET.Rp_e + 2000e3, rp=PLANET.Rp_e + 400e3, i=30.0, ω=0.0, Ω=0.0, ν=10.0)
    return SpacecraftModel(joints=Joint[], links=[root], root=root, prop_mass=prop_mass,
        inertia_tensor=root.inertia, initial_condition=ic, id=id)
end

burn_scheduler(; Δv=0.0, start=-1.0, stop=-1.0) = BaseThrusterModel(thrust=[50.0], direction=[0.0], Δv=[Δv],
    start_burn_time=[start], stop_burn_time=[stop], Isp=[300.0])

function burn_config(; scheduler=burn_scheduler(), guidance=(), rates=fill(1e6, length(guidance)),
                     prop_mass=50.0, mission_time=7000.0, settings=SimulationSettings(results=false,
                     verbose=false, generate_plots=false, normalize=false),
                     spacecraft=SpacecraftModel[burn_sat(prop_mass=prop_mass)])
    env = EnvironmentModel(planet=PLANET, EI=120.0, density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=PLANET), topography=false, wind=false)
    return SimulationConfiguration(simulation_settings=settings,
        mission_configuration=MissionConfiguration(mission_type=MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=mission_time, orientation_sim=false, num_steps_to_save=1000),
        environment_model=env,
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel(Tuple(guidance), Float64.(rates)),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel(control_effectors=(scheduler,), control_rates=[1.0]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
end

burn_binding(; delta_v_mps=10.0, spacecraft_index=1, events=[AlgorithmEvent(0.0)], kwargs...) =
    AlgorithmGuidance(; algorithm=ApoapsisBurnPlanner(delta_v_mps=delta_v_mps), spacecraft_index=spacecraft_index,
        events=events, kwargs...)

ran_scheduler(solution) = only(solution.prob.p.args.control_model.control_effectors)
same_states(a, b) = a.t == b.t && length(a.u) == length(b.u) && all(collect(x) == collect(y) for (x, y) in zip(a.u, b.u))

# A request at apoapsis-planning conditions, and the plan the planner proposes for it.
function request(; kwargs...)
    defaults = (request_id=1, reason=:forced, spacecraft_index=1, time_s=100.0,
        position_m=(PLANET.Rp_e + 400e3, 0.0, 0.0), velocity_mps=(0.0, 8200.0, 0.0), mass_kg=150.0,
        available_propellant_kg=50.0, thrust_n=50.0, isp_s=300.0, mu_m3ps2=PLANET.μ,
        planet_radius_m=PLANET.Rp_e, entry_interface_altitude_m=120e3)
    return BurnPlanRequest(; merge(defaults, (; kwargs...))...)
end
planned(r) = GA.propose!(nothing, ApoapsisBurnPlanner(delta_v_mps=10.0), r, MersenneTwister(1))
with_plan(result; kwargs...) = GA.AlgorithmResult(request_id=result.request_id, status=:candidate,
    termination=result.termination, product=PropulsiveBurnPlan(; merge(plan_fields(result.product), (; kwargs...))...))
plan_fields(p::PropulsiveBurnPlan) = NamedTuple{fieldnames(PropulsiveBurnPlan)}(Tuple(getfield(p, f) for f in fieldnames(PropulsiveBurnPlan)))
reason(r, result) = validate_burn_plan(r, result).reason

@testset "Burn-plan family declarations" begin
    planner = ApoapsisBurnPlanner(delta_v_mps=10.0)
    family = algorithm_family(planner)
    @test family === BurnPlanFamily() && planner isa AbstractBurnPlanner <: AbstractGNCAlgorithm
    @test on_rejection(family) === :keep_previous
    @test family_contract_version(family) == v"1.0.0"
    @test algorithm_result_type(family) === AlgorithmResult
    @test !consumes_rng(planner)
    @test_throws ArgumentError ApoapsisBurnPlanner(delta_v_mps=0.0)
    @test_throws ArgumentError ApoapsisBurnPlanner(delta_v_mps=NaN)
    @test_throws ArgumentError ApoapsisBurnPlanner(delta_v_mps=1.0, direction_rad=Inf)
    @test_throws ArgumentError request(request_id=0)
    @test_throws ArgumentError request(position_m=(1.0, 2.0))
    # The fallbacks, for a family that defines none of the three.
    @test_throws ArgumentError family_request(RPOReferenceFamily(), nothing, nothing, nothing, 0.0, 1, :forced, 1)
    @test_throws ArgumentError family_link(RPOReferenceFamily(), nothing, nothing, nothing, 1)
    @test family_link(RPOReferenceFamily(), :target, nothing, nothing, 1) === :target
    @test preflight_family(RPOReferenceFamily(), nothing, nothing, nothing, 1) === nothing
    @test preflight_family(RPOReferenceFamily(), nothing, nothing, nothing, 1; isolate_state=false) === nothing
end

@testset "Apoapsis planner proposes the scheduler's window" begin
    r = request()
    rng = MersenneTwister(2026); witness = copy(rng)
    result = GA.propose!(nothing, ApoapsisBurnPlanner(delta_v_mps=10.0), r, rng)
    @test rand(rng, UInt64) == rand(witness, UInt64)   # draws nothing
    @test result.status === :candidate && result.termination === :planned
    plan = result.product
    sizing = SimulationModel.ControlHooks._burn_plan_sizing(150.0, 10.0, 0.0, 50.0, 300.0, nothing)
    @test plan.valid && plan.delta_v_mps == 10.0 && plan.direction_rad == 0.0
    @test (plan.thrust_n, plan.isp_s) == (50.0, 300.0)
    @test plan.propellant_required_kg == sizing.propellant_required_kg
    @test plan.commanded_impulse_n_s == sizing.commanded_impulse_n_s
    @test r.time_s < plan.start_burn_s < plan.stop_burn_s
    @test validate_burn_plan(r, result).accepted
    @test planned(request(position_m=(PLANET.Rp_e + 50e3, 0.0, 0.0))).termination === :inside_atmosphere
    @test planned(request(velocity_mps=(-100.0, 7000.0, 0.0))).termination === :past_apoapsis
    @test planned(request(mass_kg=NaN)).termination === :invalid_mass
    @test planned(request(velocity_mps=(0.0, 20000.0, 0.0))).termination === :invalid_orbit
    @test planned(request(thrust_n=0.0)).termination === :invalid_thrust
end

@testset "A burn planner from root names only" begin
    for name in (:BaseThrusterModel, :PropulsiveBurnPlan, :burn_plan_sizing, :BurnPlanRequest, :AlgorithmResult)
        @test Base.isexported(SpaceAGORA, name)
    end
    @test SpaceAGORA.BaseThrusterModel === SimulationModel.ThrusterModels.BaseThrusterModel
    @test SpaceAGORA.PropulsiveBurnPlan === SimulationModel.PropulsiveBurnPlan
    r = request()
    sized = SpaceAGORA.burn_plan_sizing(r.mass_kg, 5.0, 0.0, r.thrust_n, r.isp_s; available_propellant_kg=r.available_propellant_kg)
    @test sized == SimulationModel.ControlHooks._burn_plan_sizing(150.0, 5.0, 0.0, 50.0, 300.0, 50.0)
    @test sized.reason === :accepted
    start = r.time_s + 60.0
    plan = SpaceAGORA.PropulsiveBurnPlan(valid=true, delta_v_mps=5.0, direction_rad=0.0, source_orbit=-1,
        start_burn_s=start, stop_burn_s=start + sized.burn_duration_s, thrust_n=r.thrust_n, isp_s=r.isp_s,
        commanded_impulse_n_s=sized.commanded_impulse_n_s, propellant_required_kg=sized.propellant_required_kg)
    result = SpaceAGORA.AlgorithmResult(request_id=r.request_id, status=:candidate, termination=:planned, product=plan)
    @test validate_burn_plan(r, result).accepted
    @test SpaceAGORA.burn_plan_sizing(150.0, 5.0, 0.0, 50.0, 300.0; available_propellant_kg=0.01).reason ===
        :insufficient_propellant
    @test isnan(SpaceAGORA.burn_plan_sizing(150.0, -1.0, 0.0, 50.0, 300.0).propellant_required_kg)
end

@testset "Trusted validation reasons" begin
    r = request()
    good = planned(r)
    plan = good.product
    @test validate_burn_plan(r, good).reason === :accepted
    @test keys(validate_burn_plan(r, good).metrics) == (:propellant_required_kg, :commanded_impulse_n_s, :burn_duration_s)
    stale = AlgorithmResult(request_id=2, status=:candidate, termination=:planned, product=plan)
    mismatch = validate_burn_plan(r, stale)
    @test mismatch.reason === :request_mismatch && !mismatch.accepted
    @test mismatch.metrics == (result_request_id=2, request_id=1)
    @test reason(r, AlgorithmResult(request_id=2, status=:infeasible, termination=:past_apoapsis)) === :request_mismatch
    @test reason(request(request_id=2), stale) === :accepted
    @test reason(r, AlgorithmResult(request_id=1, status=:infeasible, termination=:past_apoapsis)) === :infeasible
    @test reason(r, AlgorithmResult(request_id=1, status=:unsupported, termination=:x)) === :unsupported
    @test reason(r, AlgorithmResult(request_id=1, status=:candidate, termination=:x, product=(start=1.0,))) === :invalid_product
    @test reason(r, with_plan(good; valid=false)) === :invalid_product
    @test reason(request(installed_window_s=(50.0, 150.0)), good) === :burn_in_progress
    @test reason(request(installed_window_s=(150.0, 160.0)), good) === :accepted   # a later burn may be replaced
    @test reason(request(mass_kg=NaN), good) === :invalid_mass
    @test reason(r, with_plan(good; delta_v_mps=0.0)) === :invalid_delta_v
    @test reason(r, with_plan(good; direction_rad=NaN)) === :invalid_direction
    @test reason(request(thrust_n=0.0), good) === :invalid_thrust
    @test reason(request(isp_s=-1.0), good) === :invalid_isp
    @test reason(r, with_plan(good; delta_v_mps=1e-300)) === :invalid_mass_fraction
    insufficient = validate_burn_plan(request(available_propellant_kg=0.1), good)
    @test insufficient.reason === :insufficient_propellant && !insufficient.accepted
    @test insufficient.metrics.propellant_required_kg == plan.propellant_required_kg
    @test reason(request(available_propellant_kg=nothing), good) === :accepted
    @test reason(r, with_plan(good; thrust_n=49.0)) === :thrust_mismatch
    @test reason(r, with_plan(good; isp_s=301.0)) === :isp_mismatch
    @test reason(r, with_plan(good; propellant_required_kg=plan.propellant_required_kg * (1 + 1e-6))) === :propellant_mismatch
    @test reason(r, with_plan(good; commanded_impulse_n_s=plan.commanded_impulse_n_s * (1 - 1e-6))) === :impulse_mismatch
    @test reason(r, with_plan(good; stop_burn_s=plan.start_burn_s)) === :invalid_window
    @test reason(r, with_plan(good; start_burn_s=NaN)) === :invalid_window
    shift = plan.start_burn_s - 50.0
    @test reason(r, with_plan(good; start_burn_s=50.0, stop_burn_s=plan.stop_burn_s - shift)) === :window_before_request
    @test reason(r, with_plan(good; stop_burn_s=plan.stop_burn_s + 1.0)) === :duration_mismatch
    # Validation reads copies only: the request and the proposal are unchanged.
    r_copy, good_copy = deepcopy(r), deepcopy(good)
    validate_burn_plan(r, good)
    @test plan_fields(good.product) == plan_fields(good_copy.product) && r.installed_window_s === r_copy.installed_window_s
    @test burn_plan_is_current(r, plan, plan.stop_burn_s - 1.0)
    @test !burn_plan_is_current(r, plan, plan.stop_burn_s)
    @test !burn_plan_is_current(r, PropulsiveBurnPlan(), 0.0)
end

@testset "Family through the generic binding equals the scheduler's own burn" begin
    # The scheduler plans the burn itself (nonzero Δv); the binding plans the same
    # burn through the request driver and the scheduler only executes it.
    own = run_simulation(burn_config(scheduler=burn_scheduler(Δv=10.0)); return_solution=true)
    binding = burn_binding()
    viafamily = run_simulation(burn_config(guidance=(binding,)); return_solution=true)
    @test same_states(own, viafamily)
    @test ran_scheduler(own).burn_windows == ran_scheduler(viafamily).burn_windows
    @test isfinite(ran_scheduler(viafamily).burn_windows[1][1])
    @test viafamily.u[1].sc[1].mass - viafamily.u[end].sc[1].mass > 0.5
    report = algorithm_run_report(viafamily).bindings[1]
    @test [r.event for r in report.records] == [:request, :installed]
    @test report.request_id == 1 && report.product isa PropulsiveBurnPlan
    @test report.records[1].request isa BurnPlanRequest && report.records[1].request.time_s == 0.0
    @test report.product.start_burn_s == ran_scheduler(viafamily).burn_windows[1][1]
    @test report.stream_words == gnc_stream_words("AlgorithmGuidance/v1", 0, :BurnPlanFamily, 7)
    @test isempty(burn_plan_rejections(ran_scheduler(viafamily), 1))
    @test binding.runtime === nothing   # isolate_state=true: the configured binding never runs
end

@testset "A plan installed at a later tick stops the solver at its window" begin
    steps = Float64[]
    spy = DiscreteCallback((u, t, integrator) -> (push!(steps, t); false), integrator -> nothing;
        save_positions=(false, false))
    solution = run_simulation(burn_config(guidance=(burn_binding(events=[AlgorithmEvent(500.0)]),), rates=[10.0]);
        return_solution=true, extra_callbacks=(spy,))
    report = algorithm_run_report(solution).bindings[1]
    plan = report.product
    @test [(r.event, r.time_s) for r in report.records] == [(:request, 500.0), (:installed, 500.0)]
    @test ran_scheduler(solution).burn_windows == [(plan.start_burn_s, plan.stop_burn_s, plan.thrust_n)]
    @test plan.start_burn_s in steps && plan.stop_burn_s in steps
    burned = solution.u[1].sc[1].mass - solution.u[end].sc[1].mass
    @test isapprox(burned, plan.propellant_required_kg; rtol=1e-9)
end

# A planner whose burn starts at once: the scheduler must take it up at the tick
# the binding installs it, which needs the binding's tick before control's.
struct StartNowPlanner <: AbstractBurnPlanner end
function GA.propose!(state, ::StartNowPlanner, r::BurnPlanRequest, rng)
    sized = SimulationModel.ControlHooks._burn_plan_sizing(r.mass_kg, 5.0, 0.0, r.thrust_n, r.isp_s, nothing)
    plan = PropulsiveBurnPlan(valid=true, delta_v_mps=5.0, direction_rad=0.0, source_orbit=-1,
        start_burn_s=r.time_s, stop_burn_s=r.time_s + sized.burn_duration_s, thrust_n=r.thrust_n, isp_s=r.isp_s,
        commanded_impulse_n_s=sized.commanded_impulse_n_s, propellant_required_kg=sized.propellant_required_kg)
    return AlgorithmResult(request_id=r.request_id, status=:candidate, termination=:planned, product=plan)
end

@testset "A burn installed at a tick starts at that tick" begin
    trace = tempname() * ".csv"
    solution = withenv("SPACEAGORA_MANEUVER_TRACE_CSV" => trace) do
        run_simulation(burn_config(guidance=(AlgorithmGuidance(algorithm=StartNowPlanner(), spacecraft_index=1,
            events=[AlgorithmEvent(500.0)]),), rates=[10.0], mission_time=700.0); return_solution=true)
    end
    plan = algorithm_run_report(solution).bindings[1].product
    @test plan.start_burn_s == 500.0
    events = [split(line, ",") for line in readlines(trace)[2:end]]
    @test [(e[1], parse(Float64, e[2])) for e in events] == [("burn_start", 500.0), ("burn_end", plan.stop_burn_s)]
    @test isapprox(solution.u[1].sc[1].mass - solution.u[end].sc[1].mass, plan.propellant_required_kg; rtol=1e-9)
    rm(trace; force=true)
end

@testset "One binding per spacecraft, each with its own stream" begin
    two = BaseThrusterModel(thrust=[50.0, 50.0], direction=[0.0, 0.0], Δv=[0.0, 0.0],
        start_burn_time=[-1.0, -1.0], stop_burn_time=[-1.0, -1.0], Isp=[300.0, 300.0])
    template = burn_binding()
    bindings = (template, SpaceAGORA.bind_spacecraft(template, 2))
    solution = run_simulation(burn_config(scheduler=two, guidance=bindings,
        spacecraft=SpacecraftModel[burn_sat(id=7), burn_sat(id=8)]); return_solution=true)
    reports = algorithm_run_report(solution).bindings
    @test [r.spacecraft_index for r in reports] == [1, 2]
    @test all(r -> [x.event for x in r.records] == [:request, :installed], reports)
    @test reports[1].stream_words != reports[2].stream_words
    @test reports[2].stream_words == gnc_stream_words("AlgorithmGuidance/v1", 0, :BurnPlanFamily, 8)
    windows = ran_scheduler(solution).burn_windows
    @test [w[1] for w in windows] == [r.product.start_burn_s for r in reports]
    @test all(i -> solution.u[1].sc[i].mass - solution.u[end].sc[i].mass > 0.5, 1:2)
end

# A planner that returns its first result again for every later request, as a
# stale or cached result would arrive.
struct CachedResultPlanner <: AbstractBurnPlanner end
GA.initialize_algorithm(::CachedResultPlanner, run_context) = Ref{Any}(nothing)
function GA.propose!(state, ::CachedResultPlanner, r::BurnPlanRequest, rng)
    state[] === nothing && (state[] = GA.propose!(nothing, ApoapsisBurnPlanner(delta_v_mps=10.0), r, rng))
    return state[]
end

@testset "Negative control: a result for another request is refused" begin
    events = [AlgorithmEvent(0.0), AlgorithmEvent(500.0; token=2)]
    cached = run_simulation(burn_config(guidance=(AlgorithmGuidance(algorithm=CachedResultPlanner(),
        spacecraft_index=1, events=events),), rates=[10.0]); return_solution=true)
    once = run_simulation(burn_config(guidance=(AlgorithmGuidance(algorithm=CachedResultPlanner(),
        spacecraft_index=1, events=events[1:1]),), rates=[10.0]); return_solution=true)
    report = algorithm_run_report(cached).bindings[1]
    @test [r.event for r in report.records] == [:request, :installed, :request, :rejected]
    rejection = report.records[4]
    @test rejection.rejection === :request_mismatch && rejection.request_id == 2
    @test rejection.validation.metrics == (result_request_id=1, request_id=2)
    # The installed plan is the first request's, and the run equals one without the
    # second request.
    @test report.installed_request.request_id == 1
    @test report.product == algorithm_run_report(once).bindings[1].product
    @test same_states(cached, once)
end

@testset "Negative control: a refused plan keeps the coast" begin
    coast = run_simulation(burn_config(prop_mass=0.05); return_solution=true)
    refused = run_simulation(burn_config(guidance=(burn_binding(),), prop_mass=0.05); return_solution=true)
    @test same_states(coast, refused)
    @test refused.u[end].sc[1].mass == refused.u[1].sc[1].mass
    report = algorithm_run_report(refused).bindings[1]
    @test [r.event for r in report.records] == [:request, :rejected]
    rejection = report.records[2]
    @test rejection.rejection === :insufficient_propellant && rejection.reason === :forced
    @test rejection.spacecraft_index == 1 && rejection.termination === :planned
    @test !rejection.validation.accepted && isapprox(rejection.validation.metrics.available_propellant_kg, 0.05; atol=1e-12)
    @test report.product === nothing
    @test isequal(ran_scheduler(refused).burn_windows, [(NaN, NaN, 0.0)])
    # A second, affordable request installs; the earlier refusal stays on record.
    two = run_simulation(burn_config(guidance=(AlgorithmGuidance(algorithm=ApoapsisBurnPlanner(delta_v_mps=0.01),
        spacecraft_index=1, events=[AlgorithmEvent(0.0)]),), prop_mass=0.05); return_solution=true)
    @test [r.event for r in algorithm_run_report(two).bindings[1].records] == [:request, :installed]
end

@testset "The scheduler records the burns it refuses" begin
    refused = run_simulation(burn_config(scheduler=burn_scheduler(Δv=10.0), prop_mass=0.05); return_solution=true)
    records = burn_plan_rejections(ran_scheduler(refused), 1)
    @test length(records) == 1
    @test records[1] == (time_s=0.0, rejection=:insufficient_propellant, delta_v_mps=10.0, direction_rad=0.0, source_orbit=-1)
    @test refused.u[end].sc[1].mass == refused.u[1].sc[1].mass
    @test_throws ArgumentError burn_plan_rejections(ran_scheduler(refused), 2)
    # No record for a zero Δv, and every run starts with none.
    @test isempty(burn_plan_rejections(ran_scheduler(run_simulation(burn_config(); return_solution=true)), 1))
    reused = burn_scheduler(Δv=10.0)
    run_simulation(burn_config(scheduler=reused, prop_mass=0.05); isolate_state=false)
    run_simulation(burn_config(scheduler=reused, prop_mass=0.05); isolate_state=false)
    @test length(burn_plan_rejections(reused, 1)) == 1
end

@testset "Configuration checks before the run" begin
    @test_throws ArgumentError run_simulation(burn_config(guidance=(burn_binding(), burn_binding())))
    @test_throws ArgumentError run_simulation(burn_config(scheduler=burn_scheduler(Δv=1.0), guidance=(burn_binding(),)))
    @test_throws ArgumentError run_simulation(burn_config(scheduler=burn_scheduler(start=10.0, stop=20.0), guidance=(burn_binding(),)))
    @test_throws ArgumentError run_simulation(burn_config(guidance=(burn_binding(spacecraft_index=2),)))
    @test_throws ArgumentError run_simulation(burn_config(guidance=(burn_binding(settings=(margin=1.0,)),)))
    @test_throws ArgumentError run_simulation(burn_config(guidance=(burn_binding(replan_interval_s=0.0),)))
    @test_throws ArgumentError run_simulation(burn_config(guidance=(burn_binding(events=[AlgorithmEvent(0.0), AlgorithmEvent(5.0)]),)))
    thruster_guidance = ApoapsisTargetPeriapsisRaiseGuidanceModel(target_apoapsis_radius_m=PLANET.Rp_e + 3000e3,
        target_periapsis_altitude_m=500e3)
    @test_throws ArgumentError run_simulation(burn_config(guidance=(burn_binding(), thruster_guidance), rates=[1e6, 1e6]))
    checkpointing = SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false,
        checkpoint_enabled=true)
    @test_throws ArgumentError run_simulation(burn_config(guidance=(burn_binding(),), settings=checkpointing))
    # A configured link is refused when the run starts, with the binding's error.
    err = try
        run_simulation(burn_config(guidance=(burn_binding(link=Ref(0)),)))
        nothing
    catch e
        e
    end
    @test err isa AlgorithmGuidanceError && err.reason === :initialization_exception && err.cause isa ArgumentError
    @test err isa GNCAlgorithmError
end

@testset "The binding runs any family" begin
    link = Ref(0)
    counter = AlgorithmGuidance(algorithm=Counter(10), spacecraft_index=1, link=link, seed=3,
        events=[AlgorithmEvent(0.0; reason=:start), AlgorithmEvent(55.0; reason=:middle, token=2)],
        replan_interval_s=100.0)
    solution = run_simulation(burn_config(guidance=(counter,), rates=[10.0], mission_time=320.0); return_solution=true)
    report = algorithm_run_report(solution).bindings[1]
    requests = [r.request for r in report.records if r.event === :request]
    @test [(q.time_s, q.reason) for q in requests] == [(0.0, :start), (60.0, :middle), (160.0, :interval), (260.0, :interval)]
    @test report.request_id == 4 && report.product == 4
    @test report.stream_words == gnc_stream_words("AlgorithmGuidance/v1", 3, :CounterFamily, 7)
    @test link[] == 0   # isolate_state=true: the configured link object was copied
    shared = Ref(0)
    run_simulation(burn_config(guidance=(AlgorithmGuidance(algorithm=Counter(10), spacecraft_index=1, link=shared,
        events=[AlgorithmEvent(0.0)]),)); isolate_state=false)
    @test shared[] == 1
    # Algorithm code that throws fails the binding's request with its own error.
    failing = AlgorithmGuidance(algorithm=Counter(1), spacecraft_index=1, link=Ref(0),
        events=[AlgorithmEvent(0.0)], replan_interval_s=10.0)
    err = try
        run_simulation(burn_config(guidance=(failing,), rates=[10.0], mission_time=100.0))
        nothing
    catch e
        e
    end
    @test err isa AlgorithmGuidanceError && err.reason === :algorithm_exception && err.request_id == 2
    @test err.cause isa ErrorException && last(err.report.records).event === :failure
    # Two bindings on one link are refused; a bound copy has fresh run state.
    @test_throws ArgumentError run_simulation(burn_config(guidance=(failing, failing), rates=[10.0, 10.0]))
    bound = SpaceAGORA.bind_spacecraft(counter, 3)
    @test bound.spacecraft_index == 3 && bound.algorithm === counter.algorithm && bound.events === counter.events
    @test bound.runtime === nothing && bound.target === nothing
    @test SimulationModel.SimulationCallbacks._guidance_precedes_control(counter)
end

end # module BurnPlanFamilyTests
