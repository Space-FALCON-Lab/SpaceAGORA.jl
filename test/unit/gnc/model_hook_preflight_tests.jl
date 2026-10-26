using SpaceAGORA, Test
using StaticArrays
using ComponentArrays: ComponentVector

const HP_SM = SpaceAGORA.SimulationModel
const HP_SE = SpaceAGORA.SimulationEngine
const HP_ZERO3 = SVector{3, Float64}(0.0, 0.0, 0.0)

# One model per accepted or refused case of the pre-run hook check.
struct HpGuidanceNoHook <: SpaceAGORA.AbstractGuidanceModel end
struct HpGuidanceWrongTime <: SpaceAGORA.AbstractGuidanceModel end
SpaceAGORA.update_guidance!(::HpGuidanceWrongTime, u, p, t::Int, spacecraft_index::Int) = nothing
mutable struct HpGuidanceCounter <: SpaceAGORA.AbstractGuidanceModel
    calls::Int
    first_t::Float64
end
HpGuidanceCounter() = HpGuidanceCounter(0, NaN)
function SpaceAGORA.update_guidance!(m::HpGuidanceCounter, u, p, t::Float64, spacecraft_index::Int)
    m.calls == 0 && (m.first_t = t)
    m.calls += 1
    return nothing
end

# A state argument typed narrower than the ComponentVector the engine passes is
# refused; one typed as that ComponentVector is accepted.
struct HpGuidanceVectorState <: SpaceAGORA.AbstractGuidanceModel end
SpaceAGORA.update_guidance!(::HpGuidanceVectorState, u::Vector{Float64}, p, t::Float64, spacecraft_index::Int) = nothing
struct HpGuidanceComponentState <: SpaceAGORA.AbstractGuidanceModel end
SpaceAGORA.update_guidance!(::HpGuidanceComponentState, u::ComponentVector, p::HP_SM.ODEParams, t::Float64, spacecraft_index::Int) = nothing

struct HpNavigationNoHook <: SpaceAGORA.AbstractNavigationModel end
mutable struct HpNavigationCounter <: SpaceAGORA.AbstractNavigationModel
    calls::Int
    first_t::Float64
end
HpNavigationCounter() = HpNavigationCounter(0, NaN)
function SpaceAGORA.update_navigation!(m::HpNavigationCounter, u, p, t::Float64, spacecraft_index::Int)
    m.calls == 0 && (m.first_t = t)
    m.calls += 1
    return nothing
end

# A control tick but no wrench hook: refused in :ideal actuation, accepted in :mounted.
struct HpControlTickOnly <: SpaceAGORA.AbstractControlEffectorModel end
SpaceAGORA.update_control!(::HpControlTickOnly, u, p, t::Float64, i::Int) = nothing
struct HpMountedControl <: SpaceAGORA.AbstractControlEffectorModel end
SpaceAGORA.update_control!(::HpMountedControl, u, p, t::Float64, i::Int) = nothing
SpaceAGORA.actuation_mode(::HpMountedControl) = :mounted

struct HpForceNoHook <: SpaceAGORA.AbstractForceTorqueModel end
# Counts the dynamics evaluations its wrench serves.
mutable struct HpForceCounter <: SpaceAGORA.AbstractForceTorqueModel
    calls::Int
end
function SpaceAGORA.wrench(m::HpForceCounter, x::SpaceAGORA.StateSample, env::SpaceAGORA.EnvironmentSample, t::Float64)
    m.calls += 1
    return HP_ZERO3, HP_ZERO3
end

# Backbone declarations: a gravity core with its acceleration hook, one without,
# and a kick without its kick hook; none has wrench. Then a
# model that opts into the per-link kernel without the base-body hooks.
struct HpBackboneCore <: SpaceAGORA.AbstractForceTorqueModel end
SpaceAGORA.gravity_backbone_structure(::HpBackboneCore) = :position_only_static_gravity
SpaceAGORA.gravity_backbone_acceleration_ii(::HpBackboneCore, x::SpaceAGORA.StateSample, env::SpaceAGORA.EnvironmentSample, t::Float64) = HP_ZERO3
struct HpBackboneOnly <: SpaceAGORA.AbstractForceTorqueModel end
SpaceAGORA.gravity_backbone_structure(::HpBackboneOnly) = :position_only_static_gravity
struct HpKickNoAccel <: SpaceAGORA.AbstractForceTorqueModel end
SpaceAGORA.gravity_backbone_kick_structure(::HpKickNoAccel) = :velocity_kick_explicit
struct HpLinkOnly <: SpaceAGORA.AbstractForceTorqueModel end
HP_SM.link_wrench_capable(::HpLinkOnly) = true
# A spacecraft the check reads as articulated (one with a moving joint).
struct HpArticulatedSpacecraft end
HP_SM.articulated_has_moving_joints(::HpArticulatedSpacecraft) = true

struct HpDensitySixArgs <: SpaceAGORA.AbstractDensityModel end
SpaceAGORA.getDensity(::HpDensitySixArgs, h::Float64, lat::Float64, lon::Float64, t::Float64, wind::Bool) =
    (0.0, 200.0, HP_ZERO3)
struct HpDensitySevenArgs <: SpaceAGORA.AbstractDensityModel end
SpaceAGORA.getDensity(::HpDensitySevenArgs, h::Float64, lat::Float64, lon::Float64, t::Float64, wind::Bool, p) =
    (0.0, 200.0, HP_ZERO3)

include(joinpath(@__DIR__, "..", "..", "..", "templates", "force_torque_model_template.jl"))

# The parts of a configuration the check reads.
hp_args(; guidance=(), navigation=(), control=(), effectors=(), density=SpaceAGORA.NoAtmosphereModel(),
        spacecraft=(), live_pose_loads=false) = (
    guidance_model=HP_SM.GuidanceModel(guidance, fill(1.0, length(guidance))),
    navigation_model=HP_SM.NavigationModel(navigation, fill(1.0, length(navigation))),
    control_model=HP_SM.ControlModel(control, fill(1.0, length(control))),
    dynamics_model=(dynamic_effectors=effectors, spacecraft=spacecraft),
    simulation_settings=(articulated_live_pose_loads=live_pose_loads,),
    environment_model=(density_model=density, ephemerides_model=SpaceAGORA.SimpleEphemeridesModel(),
        planet=HP_SM.make_no_gram_planet(:earth)),
    initial_time=HP_SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
)

# Every concrete type defined in SpaceAGORA's modules that subtypes `T`.
const hp_package_types = SpaceAGORA.Conformance.package_types

function hp_run_config(; guidance=(), navigation=(), effectors=(), solver_mode=:tsit5)
    planet = HP_SM.make_no_gram_planet(:earth)
    root = HP_SM.Link(root=true, m=100.0, ref_area=2.0)
    ic = HP_SM.InitialCondition(ra=planet.Rp_e + 700e3, rp=planet.Rp_e + 650e3, i=45.0, ω=0.0, Ω=0.0, ν=10.0)
    sc = HP_SM.SpacecraftModel(joints=HP_SM.Joint[], links=[root], root=root, prop_mass=0.0,
        inertia_tensor=root.inertia, initial_condition=ic, id=1)
    env = HP_SM.EnvironmentModel(planet=planet, EI=120.0, density_model=SpaceAGORA.NoAtmosphereModel(),
        ephemerides_model=SpaceAGORA.SimpleEphemeridesModel(),
        thermal_model=HP_SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false, wind=false)
    return HP_SM.SimulationConfiguration(
        simulation_settings=HP_SM.SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=HP_SM.MissionConfiguration(mission_type=HP_SM.MissionTime, keplerian=true,
            number_of_orbits=1, mission_time=60.0, orientation_sim=false, num_steps_to_save=100),
        environment_model=env,
        dynamics_model=HP_SM.DynamicsModel(HP_SM.SpacecraftModel[sc], (HP_SM.InverseSquaredJ2GravityModel(), effectors...)),
        solver_config=HP_SM.SolverConfig(; solver_mode),
        guidance_model=HP_SM.GuidanceModel(guidance, fill(10.0, length(guidance))),
        navigation_model=HP_SM.NavigationModel(navigation, fill(15.0, length(navigation))),
        control_model=HP_SM.ControlModel((), Float64[]),
        initial_time=HP_SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
end

# The concrete state and parameter types of a real run, which the check uses.
const HP_FIXTURE = hp_run_config()
const HP_U0 = HP_SE.build_initial_conditions(HP_FIXTURE)
const HP_PARAMS = HP_SM.ODEParams(n_sats=1, args=HP_FIXTURE)
const HP_CALL = HP_SE._hook_call_types(HP_U0, HP_PARAMS)
const HP_P = HP_CALL.params
hp_sampled(f, T) = hasmethod(f, Tuple{T, SpaceAGORA.StateSample, SpaceAGORA.EnvironmentSample, Float64})
# Whether `T` has its own method of a backbone declaration hook, not the default.
hp_declares(f, T) = which(f, Tuple{T}).sig !== Tuple{typeof(f), Any}

@testset "Model hook declarations and pre-run check" begin
    @testset "hooks are public with generic docstrings" begin
        for name in (:AbstractNavigationModel, :update_guidance!, :update_navigation!, :update_control!,
                     :getHeatRate, :touchdown_spec, :control_wrench, :control_propellant_rate, :control_wheel_torque)
            @test name in names(SpaceAGORA)
        end
        @test SpaceAGORA.control_wrench === HP_SM.ControlHooks.control_wrench
        @test SpaceAGORA.control_propellant_rate === HP_SM.ControlHooks.control_propellant_rate
        @test SpaceAGORA.control_wheel_torque === HP_SM.ControlHooks.control_wheel_torque
        @test isabstracttype(SpaceAGORA.AbstractNavigationModel)
        @test SpaceAGORA.AbstractNavigationModel === HP_SM.AbstractTypes.AbstractNavigationModel
        @test SpaceAGORA.update_guidance! === HP_SM.GuidanceHooks.update_guidance!
        @test SpaceAGORA.update_navigation! === HP_SM.NavigationHooks.update_navigation!
        @test SpaceAGORA.update_control! === HP_SM.ControlHooks.update_control!
        @test SpaceAGORA.touchdown_spec === HP_SM.ControlHooks.touchdown_spec
        @test SpaceAGORA.getHeatRate === HP_SM.VehicleThermalModels.getHeatRate

        # The docstring each public hook shows is the generic one, not a method's.
        root_doc(name) = string(Docs.meta(SpaceAGORA)[Docs.Binding(SpaceAGORA, name)].docs[Union{}].text...)
        @test startswith(root_doc(:update_guidance!), "    update_guidance!(model, u, p, t, spacecraft_index)")
        @test startswith(root_doc(:update_navigation!), "    update_navigation!(model, u, p, t, spacecraft_index)")
        @test startswith(root_doc(:update_control!), "    update_control!(model, u, p, t, spacecraft_index)")
        @test !occursin("update_control!(model::BaseThrusterModel", root_doc(:update_control!))
        for name in (:control_wrench, :control_propellant_rate, :control_wheel_torque)
            @test startswith(root_doc(name), "    $(name)(model, x::StateSample, p, t::Float64)")
            @test !occursin("::BaseThrusterModel", root_doc(name))
        end
        # No catch-all control_wrench or tick: the pre-run check refuses a model
        # without one. The propellant rate and the wheel torque default, for any
        # model, to none used and no wheel.
        @test !hasmethod(SpaceAGORA.control_wrench, Tuple{Any, Any, Any, Any})
        @test SpaceAGORA.control_propellant_rate(HpControlTickOnly(), nothing, nothing, 0.0) === 0.0
        @test SpaceAGORA.control_wheel_torque(HpControlTickOnly(), nothing, nothing, 0.0) === nothing
        for f in (SpaceAGORA.control_propellant_rate, SpaceAGORA.control_wheel_torque)
            @test which(f, Tuple{Any, Any, Any, Any}).sig === Tuple{typeof(f), Any, Any, Any, Any}
        end
        for f in (SpaceAGORA.update_guidance!, SpaceAGORA.update_navigation!, SpaceAGORA.update_control!)
            @test !hasmethod(f, Tuple{Any, Any, Any, Any, Any})
        end
        @test startswith(root_doc(:getHeatRate), "    getHeatRate(model, S, T, ρ, v, α)")
        @test startswith(root_doc(:touchdown_spec), "    touchdown_spec(effector, spacecraft_index::Int)")
        # The ephemerides hooks, model first.
        em = HP_SM.EphemeridesModels
        for (name, arity) in ((:ephemerides_frame_lpi, 3), (:ephemerides_sun_ii, 3), (:ephemerides_epoch_seconds, 2),
                              (:ephemerides_requires_spice, 1), (:ephemerides_cache_key, 1))
            @test name in names(SpaceAGORA) && getfield(SpaceAGORA, name) === getfield(em, name)
            @test startswith(root_doc(name), "    $(name)(model")
            @test !hasmethod(getfield(SpaceAGORA, name), Tuple{ntuple(_ -> Any, arity)...})
        end
        # No catch-all navigation method: an effector without one is still a
        # MethodError when called.
        @test_throws MethodError SpaceAGORA.update_navigation!(HpNavigationNoHook(), Float64[], nothing, 0.0, 1)
    end

    @testset "the hook names removed in 0.5.0 are gone" begin
        owners = (SpaceAGORA, HP_SM, HP_SM.DynamicEffectors, HP_SM.ControlHooks, HP_SM.GuidanceHooks,
            HP_SM.NavigationHooks, HP_SM.EphemeridesModels)
        for name in (:calcForceTorque, :calcControlForceTorque, :calcControlMassFlowRate, :calcReactionWheelTorque,
                     :calcControlEffect!, :calcGuidanceEffect!, :calcNavigationEffect!,
                     :planet_frame_lpi, :ephemerides_sun_direction_ii, :ephemerides_time_seconds)
            @test !any(owner -> isdefined(owner, name), owners)
        end
        # A method qualified with the package fails where it is defined.
        scope = Module(:HpRemovedHooks)
        Core.eval(scope, :(using SpaceAGORA))
        @test_throws UndefVarError Core.eval(scope, :(SpaceAGORA.calcForceTorque(model::Int, x, p, i) = nothing))
        @test_throws UndefVarError Core.eval(scope, :(SpaceAGORA.calcControlEffect!(model::Int, u, p, t, i) = nothing))
        # The plume model's wrench takes the index from the sample only.
        plume = Tuple{SpaceAGORA.PlumeSurfaceInteractionModel, SpaceAGORA.StateSample, SpaceAGORA.EnvironmentSample, Float64}
        @test hasmethod(SpaceAGORA.wrench, plume)
        @test !hasmethod(SpaceAGORA.wrench, Tuple{plume.parameters..., Int})
    end

    @testset "the new hook names collide with no package often loaded beside SpaceAGORA" begin
        # A script that loads these packages with SpaceAGORA must see the new
        # names unqualified: a name two of them export would be ambiguous.
        new_names = (:control_wrench, :control_propellant_rate, :control_wheel_torque,
            :update_control!, :update_guidance!, :update_navigation!,
            :ephemerides_frame_lpi, :ephemerides_sun_ii, :ephemerides_epoch_seconds,
            :ephemerides_requires_spice, :ephemerides_cache_key)
        scope = Module(:HpNameCollisions)
        Core.eval(scope, :(using SpaceAGORA, Plots, OrdinaryDiffEq, ComponentArrays, StaticArrays, DataFrames))
        for pkg in (:Plots, :OrdinaryDiffEq, :ComponentArrays, :StaticArrays, :DataFrames)
            exported = names(getfield(scope, pkg))
            @test isempty(filter(in(exported), new_names))
        end
        for name in new_names
            @test isdefined(scope, name) && getfield(scope, name) === getfield(SpaceAGORA, name)
        end
    end

    @testset "the check uses the run's concrete argument types" begin
        @test HP_CALL.state === typeof(HP_U0) && HP_CALL.state <: ComponentVector
        @test HP_P === typeof(HP_PARAMS)
    end

    @testset "built-in models pass the check" begin
        # Guidance, control and force/torque models: the conformance tests
        # (test/unit/conformance/builtin_conformance_tests.jl) run this check and
        # the sample signatures on an instance of each.
        forces = hp_package_types(SpaceAGORA.AbstractForceTorqueModel)
        # Every built-in backbone declaration has its acceleration hook.
        cores = filter(T -> hp_declares(SpaceAGORA.gravity_backbone_structure, T), forces)
        kicks = filter(T -> hp_declares(SpaceAGORA.gravity_backbone_kick_structure, T), forces)
        @test !isempty(cores) && !isempty(kicks)
        for T in cores
            @test hp_sampled(SpaceAGORA.gravity_backbone_acceleration_ii, T)
        end
        for T in kicks
            @test hp_sampled(SpaceAGORA.gravity_backbone_kick_acceleration_ii, T)
        end
        # The GRAM-backed density models get getDensity from the GRAMSuite extension.
        em = HP_SM.EnvironmentModels
        gram = (em.GRAMAtmosphereModel, em.GRAMAtmosphereModelSurrogate, em.GRAMGridAtmosphereModel, em.GRAMNearSurfaceAtmosphereModel)
        gram_loaded = Base.get_extension(SpaceAGORA, :SpaceAGORAGRAMSuiteExt) !== nothing
        densities = hp_package_types(SpaceAGORA.AbstractDensityModel)
        @test length(densities) >= 10
        for T in densities
            ok = hasmethod(SpaceAGORA.getDensity, Tuple{T, Float64, Float64, Float64, Float64, Bool, HP_P})
            @test ok == (gram_loaded || !any(G -> Base.unwrap_unionall(T) <: G, gram))
        end
    end

    @testset "missing hooks are named" begin
        args = hp_args(
            guidance=(HpGuidanceNoHook(), HpGuidanceWrongTime(), HpGuidanceCounter(), HpGuidanceVectorState(), HpGuidanceComponentState()),
            navigation=(HpNavigationNoHook(), HpNavigationCounter()),
            control=(HpControlTickOnly(), HpMountedControl()),
            effectors=(HpForceNoHook(), HpForceCounter(0), HP_SM.InverseSquaredJ2GravityModel()),
            density=HpDensitySixArgs(),
        )
        lines = HP_SE._missing_model_hooks(args, HP_CALL)
        has(parts...) = any(line -> all(part -> occursin(part, line), parts), lines)
        @test length(lines) == 7
        @test has("guidance_model.guidance_effectors[1] (HpGuidanceNoHook)", "update_guidance!(model, u, p, t::Float64, spacecraft_index::Int)")
        @test has("guidance_model.guidance_effectors[2] (HpGuidanceWrongTime)", "update_guidance!", "none accepts the arguments this run passes")
        @test has("guidance_model.guidance_effectors[4] (HpGuidanceVectorState)", "update_guidance!", "none accepts the arguments this run passes")
        @test !has("(HpGuidanceNoHook)", "none accepts")
        @test has("navigation_model.navigation_effectors[1] (HpNavigationNoHook)", "update_navigation!")
        @test has("control_model.control_effectors[1] (HpControlTickOnly)",
            "has no method control_wrench(model, x::StateSample, p, t::Float64)", "actuation_mode :ideal")
        @test has("dynamics_model.dynamic_effectors[1] (HpForceNoHook)", "has no method wrench(model, x::StateSample")
        @test has("environment_model.density_model (HpDensitySixArgs)", "getDensity(", "wind::Bool, p)")
        # A :mounted effector needs no control_wrench.
        for accepted in ("HpGuidanceCounter", "HpGuidanceComponentState", "HpNavigationCounter", "HpMountedControl", "HpForceCounter", "InverseSquaredJ2GravityModel")
            @test !has(accepted)
        end

        err = try
            HP_SE._validate_model_hooks!(args, HP_U0, HP_PARAMS)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test all(line -> occursin(line, err.msg), lines)

        @test HP_SE._validate_model_hooks!(hp_args(
            guidance=(HpGuidanceCounter(),), navigation=(HpNavigationCounter(),), control=(HpMountedControl(),),
            effectors=(HpForceCounter(0),), density=HpDensitySevenArgs(),
        ), HP_U0, HP_PARAMS) === nothing

        # Outside gravity_backbone_split the backbone hooks do not count. A per-link
        # model needs no `wrench` only in a run that applies live-pose loads to
        # articulated spacecraft, where link_wrench serves it.
        special = hp_args(effectors=(HpBackboneCore(), HpBackboneOnly(), HpKickNoAccel(), HpLinkOnly()))
        plain = HP_SE._missing_model_hooks(special, HP_CALL)
        @test length(plain) == 4
        @test all(k -> any(l -> occursin("dynamic_effectors[$k] (", l) && occursin("has no method wrench(", l), plain), 1:4)
        @test any(l -> occursin("(HpLinkOnly)", l) && occursin("articulated_live_pose_loads", l), plain)
        @test !any(l -> occursin("(HpBackboneOnly)", l) && occursin("articulated_live_pose_loads", l), plain)
        for (spacecraft, live) in (((HpArticulatedSpacecraft(),), false), ((), true))
            @test length(HP_SE._missing_model_hooks(hp_args(effectors=special.dynamics_model.dynamic_effectors,
                spacecraft=spacecraft, live_pose_loads=live), HP_CALL)) == 4
        end
        live = hp_args(effectors=special.dynamics_model.dynamic_effectors, spacecraft=(HpArticulatedSpacecraft(),),
            live_pose_loads=true)
        @test length(HP_SE._missing_model_hooks(live, HP_CALL)) == 3
        @test !any(l -> occursin("HpLinkOnly", l), HP_SE._missing_model_hooks(live, HP_CALL))
        # In gravity_backbone_split a declaration needs its acceleration hook.
        split = HP_SE._missing_model_hooks(special, HP_CALL, :gravity_backbone_split)
        @test length(split) == 2
        @test any(l -> occursin("dynamic_effectors[2] (HpBackboneOnly) has no method gravity_backbone_acceleration_ii(", l), split)
        @test any(l -> occursin("dynamic_effectors[3] (HpKickNoAccel) has no method gravity_backbone_kick_acceleration_ii(", l), split)
        @test !any(l -> occursin("HpBackboneCore", l) || occursin("HpLinkOnly", l), split)

        # A model in the wrong container is pointed to its own.
        @test occursin("Its type belongs in control_model.", only(HP_SE._missing_model_hooks(hp_args(guidance=(HpControlTickOnly(),)), HP_CALL)))
        # The energy depletion controller's panel actuator cannot run on its own.
        aoa_line = only(HP_SE._missing_model_hooks(hp_args(control=(SpaceAGORA.SolarPanelAngleOfAttackControlModel(),)), HP_CALL))
        @test occursin("(SolarPanelAngleOfAttackControlModel) has no method update_control!", aoa_line)
    end

    @testset "a run with a missing hook stops before integration" begin
        force = HpForceCounter(0)
        navigation = HpNavigationCounter()
        args = hp_run_config(guidance=(HpGuidanceNoHook(),), navigation=(navigation,), effectors=(force,))
        err = try
            run_simulation(args; isolate_state=false)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("guidance_model.guidance_effectors[1] (HpGuidanceNoHook) has no method update_guidance!", err.msg)
        # Nothing was integrated: no dynamics evaluation and no navigation tick.
        @test force.calls == 0
        @test navigation.calls == 0

        # A guidance method typed for a Vector state is refused before the run too:
        # the run passes a ComponentVector.
        err = try
            run_simulation(hp_run_config(guidance=(HpGuidanceVectorState(),), effectors=(force,)); isolate_state=false)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("(HpGuidanceVectorState) has no method update_guidance!", err.msg)
        @test force.calls == 0

        # A backbone gravity core without its acceleration hook, in gravity_backbone_split.
        err = try
            run_simulation(hp_run_config(effectors=(HpBackboneOnly(),), solver_mode=:gravity_backbone_split); isolate_state=false)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("(HpBackboneOnly) has no method gravity_backbone_acceleration_ii(", err.msg)

        # With the hook present the same run integrates.
        guidance = HpGuidanceCounter()
        run_simulation(hp_run_config(guidance=(guidance,), navigation=(navigation,), effectors=(force,)); isolate_state=false)
        @test force.calls > 0
        @test guidance.calls > 0 && navigation.calls > 0
        # The first tick comes one period after the start, as the docstrings say.
        @test guidance.first_t == 10.0
        @test navigation.first_t == 15.0
    end

    @testset "the force template passes the check" begin
        template = ExampleForceTorqueModelTemplate.ExampleForceTorqueModel(HP_ZERO3, HP_ZERO3)
        @test isempty(HP_SE._missing_model_hooks(hp_args(effectors=(template,)), HP_CALL))
    end

    @testset "the shared deprecation helper honours --depwarn=error" begin
        helper = joinpath(pkgdir(SpaceAGORA), "src", "foundation", "core", "interfaces", "hook_deprecation.jl")
        script = "include($(repr(helper))); try HookDeprecation._hook_deprecation((Int, :hook), \"msg\"); " *
                 "print(\"no error\") catch err; print(err isa ErrorException && err.msg == \"msg\") end"
        @test read(`$(Base.julia_cmd()) --startup-file=no --depwarn=error -e $script`, String) == "true"
        @test read(`$(Base.julia_cmd()) --startup-file=no --depwarn=no -e $script`, String) == "no error"
    end
end
