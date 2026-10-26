# The EDG guidance family of the GNC algorithm contract: its declarations, its
# trusted validation, and its runs through the generic binding (`AlgorithmGuidance`),
# compared with the EDG guidance model's.
module EDGGuidanceFamilyTests
using Test, SpaceAGORA, StaticArrays, LinearAlgebra, Random
import .._edg_test_context, .._ODYSSEY_RA_INITIAL_M
const S = SpaceAGORA.SimulationModel
const SM = S
const GA = SpaceAGORA.GNCAlgorithms
const L = SpaceAGORA.SimulationLifecycle
const F = SpaceAGORA.EDGGuidanceFamilies
const CB = S.SimulationCallbacks
const SC = S.SimConfig
include(joinpath(@__DIR__, "..", "..", "fixtures", "edg_guidance_conformance.jl"))

const FAMILY = EDGGuidanceFamily()
const ALGORITHM = EDGGuidanceAlgorithm()
state_values(s) = Tuple(deepcopy(getfield(s, f)) for f in fieldnames(typeof(s)))
replace_fields(d::EDGGuidanceDecision; kw...) =
    EDGGuidanceDecision((get(kw, f, getfield(d, f)) for f in fieldnames(EDGGuidanceDecision))...)
candidate(request, decision) = AlgorithmResult(request_id=request.request_id, status=:candidate,
    termination=:decided, product=decision)

# Records the EDG state after the EDG control's tick at the same rate.
struct FamilyProbe <: S.AbstractTypes.AbstractControlEffectorModel
    state::SpaceAGORA.AerobrakingEnergyDepletionState
    records::Vector{NamedTuple}
end
function S.update_control!(probe::FamilyProbe, u, p, t::Float64, i::Int64)
    st = probe.state
    push!(probe.records, (time=t, count=st.energy_bracketing_count[i], last=st.last_switch_solve_t[i],
        switch=st.targeting_switch_s[i], diagnostic=st.targeting_diagnostics[i],
        alpha=p.args.dynamics_model.spacecraft[i].links[2].α))
    return nothing
end

# An algorithm of the family that supports maximum depletion only.
struct MaxOnlyEDG <: SpaceAGORA.AbstractGNCAlgorithm end
SpaceAGORA.algorithm_family(::MaxOnlyEDG) = EDGGuidanceFamily()
SpaceAGORA.algorithm_capabilities(::MaxOnlyEDG) = EDGGuidanceCapabilities(guidance_modes=(:max_energy_depletion,))

# An algorithm of the family whose decisions select a mode the configuration lacks.
struct UnconfiguredModeEDG <: SpaceAGORA.AbstractGNCAlgorithm end
SpaceAGORA.algorithm_family(::UnconfiguredModeEDG) = EDGGuidanceFamily()
SpaceAGORA.algorithm_capabilities(::UnconfiguredModeEDG) = EDGGuidanceCapabilities(guidance_modes=(:max_energy_depletion,))
SpaceAGORA.consumes_rng(::UnconfiguredModeEDG) = false
SpaceAGORA.propose!(state, ::UnconfiguredModeEDG, r::EDGGuidanceRequest, rng) = candidate(r,
    replace_fields(EDGGuidanceDecision(r.state, r.spacecraft_index); selected_mode=:targeting))

@testset "EDG guidance family declarations" begin
    @test SpaceAGORA.algorithm_family(ALGORITHM) == FAMILY
    @test SpaceAGORA.default_capabilities(FAMILY).guidance_modes == ()
    @test SpaceAGORA.algorithm_capabilities(ALGORITHM).guidance_modes == (:max_energy_depletion, :targeting)
    @test !SpaceAGORA.consumes_rng(ALGORITHM)
    @test SpaceAGORA.require_algorithm_support(ALGORITHM) === nothing
    @test SpaceAGORA.on_rejection(FAMILY) === :fail
    @test SpaceAGORA.family_contract_version(FAMILY) == SpaceAGORA.EDGServices.CONTRACT_VERSION == v"3.0.0"
    @test SpaceAGORA.algorithm_result_type(FAMILY) === AlgorithmResult
    # The decision holds exactly the state's guidance fields, one value per field.
    @test fieldnames(EDGGuidanceDecision) == (:spacecraft_index, F.DECISION_FIELDS...)
    @test all(f -> f in fieldnames(SpaceAGORA.AerobrakingEnergyDepletionState), F.DECISION_FIELDS)
end

@testset "EDG guidance requests are copies and proposals change only the decision" begin
    for modes in ((:max_energy_depletion,), (:targeting, :max_energy_depletion))
        ctx = _edg_test_context(guidance_modes=modes, density_reference_kg_m3=1e-7)
        settings = EDGGuidanceSettings(config=ctx.config, state=ctx.state)
        request = SpaceAGORA.family_request(FAMILY, settings, ctx.u, ctx.p, 0.0, 1, :interval, 3)
        @test request.request_id == 3 && request.spacecraft_index == 1 && request.time_s == 0.0
        @test request.config === ctx.config
        @test request.state !== ctx.state && isequal(state_values(request.state), state_values(ctx.state))
        @test propertynames(request.services) == ()
        @test_throws ArgumentError request.services.params
        @test deepcopy(request).services === request.services
        observed = request.observation.sc[1].pos
        ctx.u.sc[1].pos[3] += 1.0
        @test request.observation.sc[1].pos == observed
        ctx.u.sc[1].pos[3] -= 1.0

        before = state_values(request.state)
        rng = MersenneTwister(7)
        witness = copy(rng)
        result = SpaceAGORA.propose!(nothing, ALGORITHM, request, rng)
        @test rand(rng, UInt64) == rand(witness, UInt64)
        @test isequal(state_values(request.state), before)
        @test result.status === :candidate && result.request_id == 3
        # The decision is what the provider writes into the state it is given.
        reference = deepcopy(ctx.state)
        S.EDGAlgorithms.guidance_decision!(reference, ctx.config, ctx.u, ctx.p, 0.0, 1)
        @test isequal(result.product, EDGGuidanceDecision(reference, 1))
        @test !F._changed_outside(ctx.state, reference, 1)
        if :targeting in modes
            @test result.product.energy_bracketing_count == 1
            @test result.product.targeting_diagnostics.status == :no_admissible_family
        else
            @test result.product.selected_mode === :max_energy_depletion
        end
        accepted = validate_edg_guidance_decision(request, result)
        @test accepted.accepted && accepted.reason === :accepted
        # Installation writes the decision fields of the spacecraft and nothing else.
        SpaceAGORA.install!(ctx.state, FAMILY, result.product, request, 0.0; validation=accepted,
            diagnostics=NamedTuple())
        @test isequal(state_values(ctx.state), state_values(reference))
        @test SpaceAGORA.is_current(FAMILY, request, result.product, 0.0)
        @test !SpaceAGORA.is_current(FAMILY, request, result.product, -1.0)
    end
    # A provider that writes outside the decision is caught.
    a = SpaceAGORA.AerobrakingEnergyDepletionState(num_sats=2)
    b = deepcopy(a); b.selected_mode[1] = :targeting; b.energy_bracketing_count[1] = 4
    @test !F._changed_outside(a, b, 1)
    c = deepcopy(a); c.last_alpha_rad[1] = 0.2
    @test F._changed_outside(a, c, 1)
    d = deepcopy(a); d.selected_mode[2] = :targeting
    @test F._changed_outside(a, d, 1)
    e = deepcopy(a); e.control_clock_period_s = 1.0
    @test F._changed_outside(a, e, 1)
end

@testset "EDG guidance validation reasons" begin
    ctx = _edg_test_context(guidance_modes=(:max_energy_depletion,))
    request = SpaceAGORA.family_request(FAMILY, EDGGuidanceSettings(config=ctx.config, state=ctx.state),
        ctx.u, ctx.p, 0.0, 1, :interval, 1)
    request.state.energy_bracketing_count[1] = 2
    decision = replace_fields(EDGGuidanceDecision(request.state, 1); selected_mode=:max_energy_depletion)
    reason(result) = validate_edg_guidance_decision(request, result).reason
    @test reason(candidate(request, decision)) === :accepted
    @test reason(AlgorithmResult(request_id=2, status=:candidate, termination=:decided, product=decision)) ===
        :request_mismatch
    @test reason(AlgorithmResult(request_id=1, status=:infeasible, termination=:none)) === :infeasible
    @test reason(candidate(request, 1.0)) === :invalid_product
    @test reason(candidate(request, replace_fields(decision; spacecraft_index=2))) === :spacecraft_mismatch
    @test reason(candidate(request, replace_fields(decision; selected_mode=:coast))) === :unknown_mode
    @test reason(candidate(request, replace_fields(decision; selected_mode=:targeting))) === :mode_not_configured
    @test reason(candidate(request, replace_fields(decision; selected_mode=:inactive))) === :accepted
    @test reason(candidate(request, replace_fields(decision; selected_mode=:safe_low_drag))) === :accepted
    @test reason(candidate(request, replace_fields(decision; energy_bracketing_count=1))) === :count_decreased
    # NaN is no switch time (Inf is the no-switch sentinel), and an active targeting
    # plan is the targeting mode with a switch time.
    @test reason(candidate(request, replace_fields(decision; targeting_switch_s=NaN))) === :invalid_switch_time
    @test reason(candidate(request, replace_fields(decision; targeting_active=true))) === :inconsistent_targeting
    targeting = _edg_test_context(guidance_modes=(:targeting, :max_energy_depletion))
    plan_request = SpaceAGORA.family_request(FAMILY, EDGGuidanceSettings(config=targeting.config, state=targeting.state),
        targeting.u, targeting.p, 0.0, 1, :interval, 1)
    plan = replace_fields(EDGGuidanceDecision(plan_request.state, 1); selected_mode=:targeting, targeting_active=true,
        targeting_switch_s=42.0)
    plan_reason(d) = validate_edg_guidance_decision(plan_request, candidate(plan_request, d)).reason
    @test plan_reason(plan) === :accepted
    @test plan_reason(replace_fields(plan; targeting_switch_s=NaN)) === :invalid_switch_time
    @test plan_reason(replace_fields(plan; targeting_switch_s=Inf)) === :invalid_switch_time
    @test plan_reason(replace_fields(plan; targeting_active=false)) === :inconsistent_targeting
    @test plan_reason(replace_fields(plan; selected_mode=:max_energy_depletion)) === :inconsistent_targeting
    @test plan_reason(replace_fields(plan; selected_mode=:max_energy_depletion, targeting_active=false)) === :accepted
end

@testset "EDG guidance capabilities and configuration checks" begin
    ctx = _edg_test_context(guidance_modes=(:targeting, :max_energy_depletion))
    settings = EDGGuidanceSettings(config=ctx.config, state=ctx.state)
    request = SpaceAGORA.family_request(FAMILY, settings, ctx.u, ctx.p, 0.0, 1, :interval, 1)
    @test SpaceAGORA.validate_capabilities(FAMILY, ALGORITHM, request).accepted
    refused = SpaceAGORA.validate_capabilities(FAMILY, MaxOnlyEDG(), request)
    @test !refused.accepted && refused.reason === :unsupported_guidance_mode && refused.metrics.guidance_modes == (:targeting,)
    @test SpaceAGORA.validate_capabilities(FAMILY, ALGORITHM, request; checkpoint=true).reason === :restart_unsupported
    @test SpaceAGORA.family_link(FAMILY, nothing, settings, ctx.p, 1) === ctx.state
    @test_throws ArgumentError SpaceAGORA.family_link(FAMILY, ctx.state, settings, ctx.p, 1)
    @test SpaceAGORA.preflight_family(FAMILY, ALGORITHM, settings, ctx.args, 1) === nothing
    @test_throws ArgumentError SpaceAGORA.preflight_family(FAMILY, ALGORITHM, ctx.config, ctx.args, 1)
    two = (dynamics_model=(spacecraft=[ctx.spacecraft, ctx.spacecraft],),)
    @test_throws ArgumentError SpaceAGORA.preflight_family(FAMILY, ALGORITHM,
        EDGGuidanceSettings(config=ctx.config, state=SpaceAGORA.AerobrakingEnergyDepletionState(num_sats=2)), two, 1)
    @test_throws ArgumentError SpaceAGORA.preflight_family(FAMILY, ALGORITHM,
        EDGGuidanceSettings(config=ctx.config, state=SpaceAGORA.AerobrakingEnergyDepletionState(num_sats=2)), ctx.args, 1)

    # In a run the binding is the EDG guidance: it pairs with the EDG control by state
    # and configuration, and it cannot run beside the EDG guidance model.
    args = edg_guidance_configuration()
    binding = only(args.guidance_model.guidance_effectors)
    control = only(args.control_model.control_effectors)
    @test L.requires_atmosphere_events(binding)
    @test SpaceAGORA.gnc_links(binding) == (SpaceAGORA.Produces(control.state),)
    @test L.preflight_control(control, args; isolate_state=true) === nothing
    apart = SC._with_configuration(args; control_model=S.ControlModel((S.AerobrakingEnergyDepletionControlModel(
        control.config, SpaceAGORA.AerobrakingEnergyDepletionState(num_sats=1)),), [1 / 3]))
    @test_throws ArgumentError L.preflight_control(only(apart.control_model.control_effectors), apart; isolate_state=true)
    both = SC._with_configuration(args; guidance_model=S.GuidanceModel(
        (binding, S.AerobrakingEnergyDepletionGuidanceModel(control.config, control.state)), [1 / 3, 1 / 3]))
    @test_throws ArgumentError L.preflight_control(control, both; isolate_state=true)
    # Targeting control needs its paired guidance; the binding is that guidance.
    targeting = edg_guidance_configuration(guidance_modes=(:targeting, :max_energy_depletion))
    @test L.preflight_control(only(targeting.control_model.control_effectors), targeting; isolate_state=true) === nothing
    alone = SC._with_configuration(targeting; guidance_model=S.GuidanceModel((), Float64[]))
    @test_throws ArgumentError L.preflight_control(only(alone.control_model.control_effectors), alone; isolate_state=true)
    # Another control model beside the EDG control, before or after it and sharing
    # its supertype, is neither paired with the binding nor counted.
    panel = S.SolarPanelAngleOfAttackControlModel()
    unpaired = only(apart.control_model.control_effectors)
    for before in (true, false)
        mixed = SC._with_configuration(args; control_model=S.ControlModel(
            before ? (panel, control) : (control, panel), [1 / 3, 1 / 3]))
        @test L.preflight_control(control, mixed; isolate_state=true) === nothing
        mixed_apart = SC._with_configuration(args; control_model=S.ControlModel(
            before ? (panel, unpaired) : (unpaired, panel), [1 / 3, 1 / 3]))
        err = try
            L.preflight_control(unpaired, mixed_apart; isolate_state=true)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("must share the same mutable state", sprint(showerror, err))
    end
end

@testset "EDG guidance through AlgorithmGuidance runs as the EDG guidance model" begin
    model = run_simulation(edg_guidance_configuration(guidance=:model, duration_s=2.0);
        return_results=true, return_solution=true)
    bound = run_simulation(edg_guidance_configuration(guidance=:binding, duration_s=2.0);
        return_results=true, return_solution=true)
    @test model.solution.t == bound.solution.t
    @test model.solution.u == bound.solution.u
    model_state = only(model.configuration.control_model.control_effectors).state
    bound_state = only(bound.configuration.control_model.control_effectors).state
    @test bound_state === only(bound.configuration.guidance_model.guidance_effectors).settings.state
    @test isequal(state_values(model_state), state_values(bound_state))
    @test only(bound_state.selected_mode) === :max_energy_depletion
    @test isfinite(only(bound_state.last_alpha_rad))
    # One request at the run's start and one at each of the five guidance ticks
    # before the 2 s end (every 1/3 s), each installed.
    report = only(algorithm_run_report(bound.solution).bindings)
    events = [r.event for r in report.records]
    requests = count(==(:request), events)
    @test requests == report.request_id == 6
    @test count(==(:installed), events) == requests
    @test !(:failure in events)
    @test all(r.product isa EDGGuidanceDecision for r in report.records if r.event === :installed)
end

@testset "EDG targeting through AlgorithmGuidance replans as the EDG guidance model" begin
    # The assembled production callbacks of the replan-owner test, once with the EDG
    # guidance model and once with the binding, on identical copies of one state.
    c = _edg_test_context(guidance_modes=(:targeting, :max_energy_depletion), density_reference_kg_m3=1e-7,
        heat_load_limit_j_cm2=30.0, gravity_models=(S.InverseSquaredJ2GravityModel(),))
    c.args.control_model.control_rates[1] = 1 / 3
    c.state.control_clock_period_s = 1 / 3
    planet = c.args.environment_model.planet
    rp = planet.Rp_e + 100e3; ra = _ODYSSEY_RA_INITIAL_M; a = (rp + ra) / 2
    e = (ra - rp) / (ra + rp); sl = a * (1 - e^2); radius = planet.Rp_e + 150e3
    nu = acos((sl / radius - 1) / e)
    r = SVector(radius * cos(nu), radius * sin(nu), 0.0)
    v = sqrt(planet.μ / sl) * SVector(-sin(nu), e + cos(nu), 0.0)
    # The replan-owner test's reachable target, as SpaceAGORAEDG's kernels give it
    # (recorded on v0.6.0; SpaceAGORAEDG's tests check the kernels still do).
    target = 9.750451931596396e6
    fields = (; (k => getfield(c.config, k) for k in fieldnames(typeof(c.config)))...)
    config = S.AerobrakingEnergyDepletionConfig(; merge(fields,
        (target_apoapsis_radius_m=target, switch_recompute_interval_s=2.0))...)
    c.state.heat_load_entry_j_cm2[1] = [100.0, 6.0, 6.0]
    function fixture_rhs!(du, u, p, t)
        fill!(du, 0.0)
        du.sc[1].heat_loads[2] = 0.05
        du.sc[1].heat_loads[3] = 0.05
    end
    function run(kind)
        state = deepcopy(c.state)
        guide = kind === :model ? S.AerobrakingEnergyDepletionGuidanceModel(config, state) :
            AlgorithmGuidance(; algorithm=ALGORITHM, spacecraft_index=1,
                settings=EDGGuidanceSettings(config=config, state=state), replan_interval_s=1 / 3)
        control = S.AerobrakingEnergyDepletionControlModel(config, state)
        probe = FamilyProbe(state, NamedTuple[])
        afields = (; (k => getfield(c.args, k) for k in fieldnames(typeof(c.args)))...)
        args = S.SimulationConfiguration(; merge(afields, (
            simulation_settings=S.SimulationSettings(results=false, verbose=false, generate_plots=false, save_csv=false),
            guidance_model=S.GuidanceModel((guide,), [1 / 3]),
            control_model=S.ControlModel((control, probe), [1 / 3, 1 / 3])))...)
        params = S.ODEParams(n_sats=1, args=args)
        params.shared_buffers.et_start[] = c.p.shared_buffers.et_start[]
        u = deepcopy(c.u); u.sc[1].pos .= r; u.sc[1].vel .= v
        u.sc[1].heat_loads .= [100.0, 8.0, 9.0]
        # The guidance decides once at the start: the model by a direct call, the
        # binding by its initialization request.
        kind === :model ? S.update_guidance!(guide, u, params, 0.0, 1) : L.initialize_guidance!(guide, u, params, 0.0)
        callbacks = CB.get_callbacks(1, args.dynamics_model.dynamic_effectors, args)
        prob = CB.ODEProblem(fixture_rhs!, u, (0.0, 4.1), params)
        sol = CB.solve(prob, CB.Tsit5(); callback=callbacks, dtmax=0.1, abstol=1e-10, reltol=1e-10,
            save_everystep=false)
        return (state=state, log=probe.records, sol=sol, guide=guide)
    end
    model, bound = run(:model), run(:binding)
    @test length(model.log) == 12
    @test isequal(model.log, bound.log)
    @test model.sol.u[end] == bound.sol.u[end]
    @test isequal(state_values(model.state), state_values(bound.state))
    @test model.state.energy_bracketing_count[1] == 3
    @test only(model.state.targeting_diagnostics).mode === :targeting
    records = bound.guide.runtime.records
    @test count(r -> r.event === :installed, records) == count(r -> r.event === :request, records) == 13
end

@testset "A targeting binding that starts inside the drag passage plans at once" begin
    # The binding's first request comes before the solve starts; the run reports the
    # configured control clock when it starts, so that plan has its clock.
    c = _edg_test_context(guidance_modes=(:targeting, :max_energy_depletion), density_reference_kg_m3=1e-7)
    @test c.u.sc[1].pos[1] - c.args.environment_model.planet.Rp_e < 1e3 * c.args.environment_model.EI
    state = SpaceAGORA.AerobrakingEnergyDepletionState(num_sats=1)
    binding = AlgorithmGuidance(; algorithm=ALGORITHM, spacecraft_index=1,
        settings=EDGGuidanceSettings(config=c.config, state=state), replan_interval_s=1.0)
    args = SC._with_configuration(c.args;
        guidance_model=S.GuidanceModel((binding,), [1.0]),
        control_model=S.ControlModel((S.AerobrakingEnergyDepletionControlModel(c.config, state),), [1.0]),
        mission_configuration=S.MissionConfiguration(mission_type=S.MissionTime, mission_time=0.5, orientation_sim=false),
        simulation_settings=S.SimulationSettings(results=false, verbose=false, generate_plots=false),
        solver_config=S.SolverConfig(solver_mode=:tsit5))
    result = run_simulation(args; return_results=true, return_solution=true)
    report = only(algorithm_run_report(result.solution).bindings)
    first_install = first(r for r in report.records if r.event === :installed)
    @test first_install.time_s == 0.0
    @test first_install.product.energy_bracketing_count == 1
    @test first_install.product.targeting_diagnostics.status in (:no_admissible_family, :exact_sample, :quantized_bracket)
    run_state = only(result.configuration.control_model.control_effectors).state
    @test run_state.control_clock_period_s == 1.0 && run_state.control_clock_origin_s == 0.0
    @test !(:failure in [r.event for r in report.records])
end

@testset "A refused EDG guidance decision stops the run" begin
    configuration = edg_guidance_configuration(algorithm=UnconfiguredModeEDG())
    state = only(configuration.control_model.control_effectors).state
    # Any algorithm of the family run with EDG settings is the run's EDG guidance.
    @test S.GuidanceHooks._edg_guidance_pairing(only(configuration.guidance_model.guidance_effectors)).state === state
    failure = try
        run_simulation(configuration; isolate_state=false)
        nothing
    catch err
        err
    end
    @test failure isa AlgorithmGuidanceError
    @test failure.reason === :mode_not_configured
    @test failure.request_id == 1
    # The refused decision was never installed.
    @test only(state.selected_mode) === :inactive
end

end # module EDGGuidanceFamilyTests
