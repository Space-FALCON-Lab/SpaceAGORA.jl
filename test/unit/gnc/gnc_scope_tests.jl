# Spacecraft scope, inactive spacecraft and held commands (GNCRuntime): the
# spacecraft the tick loops call a model for, the release that replaces the tick of
# an inactive spacecraft, the check before the run, and the built-in declarations.
module GNCScopeTests
using Test, Logging, StaticArrays, SpaceAGORA
using DiffEqCallbacks: PresetTimeCallback
const SM = SpaceAGORA.SimulationModel
const L = SpaceAGORA.SimulationLifecycle
const ZERO3 = SVector(0.0, 0.0, 0.0)

# Models that count, per spacecraft, their ticks and the releases in their place.
mutable struct Counts
    ticks::Vector{Int}
    released::Vector{Int}
end
Counts(n) = Counts(zeros(Int, n), zeros(Int, n))

struct ScopedControl <: SM.AbstractControlEffectorModel
    scope::Any
    skip::Any
    counts::Counts
end
ScopedControl(scope, n; skip=false) = ScopedControl(scope, skip, Counts(n))
SpaceAGORA.update_control!(m::ScopedControl, u, p, t, i) = (m.counts.ticks[i] += 1; nothing)
SpaceAGORA.control_wrench(::ScopedControl, x::StateSample, p, t) = (ZERO3, ZERO3)

struct ScopedGuidance <: AbstractGuidanceModel
    scope::Any
    skip::Any
    counts::Counts
end
ScopedGuidance(scope, n; skip=false) = ScopedGuidance(scope, skip, Counts(n))
SpaceAGORA.update_guidance!(m::ScopedGuidance, u, p, t::Float64, i::Int) = (m.counts.ticks[i] += 1; nothing)

struct ScopedNavigation <: AbstractNavigationModel
    scope::Any
    skip::Any
    counts::Counts
end
ScopedNavigation(scope, n; skip=false) = ScopedNavigation(scope, skip, Counts(n))
SpaceAGORA.update_navigation!(m::ScopedNavigation, u, p, t::Float64, i::Int) = (m.counts.ticks[i] += 1; nothing)

const Scoped = Union{ScopedControl, ScopedGuidance, ScopedNavigation}
SpaceAGORA.spacecraft_scope(m::Scoped) = m.scope
SpaceAGORA.skip_inactive(m::Scoped) = m.skip
SpaceAGORA.release_command!(m::Scoped, t, i) = (m.counts.released[i] += 1; nothing)

# A model with no declaration of its own.
struct PlainControl <: SM.AbstractControlEffectorModel end

# A model with the default declarations that counts its ticks per spacecraft.
struct CountingControl <: SM.AbstractControlEffectorModel
    ticks::Vector{Int}
end
SpaceAGORA.update_control!(m::CountingControl, u, p, t, i) = (m.ticks[i] += 1; nothing)
SpaceAGORA.control_wrench(::CountingControl, x::StateSample, p, t) = (ZERO3, ZERO3)

function constellation_args(; n=3, control=(), guidance=(), navigation=(), mission_time=5.0)
    planet = SM.make_no_gram_planet(:earth)
    sats = map(1:n) do id
        root = SM.Link(root=true, m=100.0, ref_area=2.0)
        ic = SM.InitialCondition(ra=planet.Rp_e + 700e3, rp=planet.Rp_e + 650e3, i=45.0, ω=0.0, Ω=30.0id, ν=10.0)
        SM.SpacecraftModel(joints=SM.Joint[], links=[root], root=root, prop_mass=0.0, inertia_tensor=root.inertia,
            initial_condition=ic, id=id)
    end
    env = SM.EnvironmentModel(planet=planet, EI=120.0, density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(),
        thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet), topography=false, wind=false)
    return SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=mission_time, orientation_sim=false, num_steps_to_save=100),
        environment_model=env,
        dynamics_model=SM.DynamicsModel(SM.SpacecraftModel[sats...], (SM.InverseSquaredJ2GravityModel(),)),
        guidance_model=SM.GuidanceModel(Tuple(guidance), fill(1.0, length(guidance))),
        navigation_model=SM.NavigationModel(Tuple(navigation), fill(1.0, length(navigation))),
        control_model=SM.ControlModel(control_effectors=Tuple(control), control_rates=fill(1.0, length(control))),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
end

# The message of the `ArgumentError` that `f()` refuses with, or `nothing` when it
# returns.
function refusal(f)
    try
        f()
        return nothing
    catch e
        e isa ArgumentError || rethrow()
        return sprint(showerror, e)
    end
end
scope_refusal(args) = refusal(() -> SM.SimulationCallbacks.check_spacecraft_scopes(args))

# Deactivate spacecraft `i` at time `t_s`, as an impact or touchdown event does.
deactivate(i, t_s) = PresetTimeCallback([t_s], integrator -> (integrator.p.is_active[i] = false); save_positions=(false, false))

@testset "GNC spacecraft scope" begin
    @testset "the defaults" begin
        @test spacecraft_scope(:any_model) === AllSpacecraft()
        @test skip_inactive(:any_model) === true
        @test held_command(:any_model, 1) === nothing
        @test isnan(held_command_time(:any_model, 1))
        @test release_command!(:any_model, 0.0, 1) === nothing
        @test spacecraft_scope === SpaceAGORA.GNCRuntime.spacecraft_scope
        @test BoundSpacecraft((Int8(2), 3)) == BoundSpacecraft((2, 3))
        @test CoupledSpacecraft(UInt(1), (Int8(2),)) == CoupledSpacecraft(1, (2,))
        @test BoundSpacecraft(()) isa SpacecraftScope
    end

    @testset "the tick loops call a model only for the spacecraft in its scope" begin
        all_c = ScopedControl(AllSpacecraft(), 3)
        bound_c = ScopedControl(BoundSpacecraft((2,)), 3)
        coupled_c = ScopedControl(CoupledSpacecraft(3, (1,)), 3)
        bound_g = ScopedGuidance(BoundSpacecraft((1, 3)), 3)
        bound_n = ScopedNavigation(BoundSpacecraft((3,)), 3)
        run_simulation(constellation_args(control=(all_c, bound_c, coupled_c), guidance=(bound_g,),
            navigation=(bound_n,)); isolate_state=false)
        k = all_c.counts.ticks[1]
        @test k > 0 && all_c.counts.ticks == [k, k, k]
        @test bound_c.counts.ticks == [0, k, 0]
        @test coupled_c.counts.ticks == [0, 0, k]          # the primary only
        @test bound_g.counts.ticks == [k, 0, k]
        @test bound_n.counts.ticks == [0, 0, k]
        @test all(m -> iszero(m.counts.released), (all_c, bound_c, coupled_c, bound_g, bound_n))
    end

    @testset "an inactive spacecraft's command is released instead of ticked" begin
        skipping = (ScopedControl(AllSpacecraft(), 2; skip=true), ScopedGuidance(AllSpacecraft(), 2; skip=true),
                    ScopedNavigation(AllSpacecraft(), 2; skip=true))
        ticking = ScopedControl(AllSpacecraft(), 2)
        default = CountingControl(zeros(Int, 2))
        args = constellation_args(n=2, control=(skipping[1], ticking, default), guidance=(skipping[2],),
            navigation=(skipping[3],))
        run_simulation(args; isolate_state=false, extra_callbacks=(deactivate(2, 2.5),))
        for m in skipping
            k = m.counts.ticks[1]
            @test m.counts.released[1] == 0
            @test 0 < m.counts.released[2] < k && m.counts.ticks[2] + m.counts.released[2] == k
        end
        # The default, skip_inactive = true: no tick for the inactive spacecraft.
        @test default.ticks == skipping[1].counts.ticks
        # A model that declares skip_inactive = false keeps ticking.
        @test ticking.counts.ticks[2] == ticking.counts.ticks[1] && iszero(ticking.counts.released)
        # A coupled model follows its primary only.
        coupled = ScopedControl(CoupledSpacecraft(1, (2,)), 2; skip=true)
        run_simulation(constellation_args(n=2, control=(coupled,)); isolate_state=false,
            extra_callbacks=(deactivate(2, 2.5),))
        @test coupled.counts.ticks[1] > 0 && coupled.counts.ticks[2] == 0 && iszero(coupled.counts.released)
    end

    @testset "the engine's step for one spacecraft" begin
        p = (is_active=[true, false],)
        tick!(m, u, p, t, i) = (m.counts.ticks[i] += 1; nothing)
        step!(m, i) = SM.SimulationCallbacks._gnc_tick!(tick!, m, nothing, p, 0.0, i)
        out = ScopedControl(BoundSpacecraft((1,)), 2; skip=true)
        @test !step!(out, 2) && iszero(out.counts.ticks) && iszero(out.counts.released)
        released = ScopedControl(AllSpacecraft(), 2; skip=true)
        @test step!(released, 2) && released.counts.released == [0, 1] && iszero(released.counts.ticks)
        ticked = ScopedControl(AllSpacecraft(), 2)
        @test step!(ticked, 2) && ticked.counts.ticks == [0, 1] && iszero(ticked.counts.released)
    end

    @testset "scopes are checked before the run" begin
        @test scope_refusal(constellation_args(control=(ScopedControl(BoundSpacecraft((3,)), 3),))) === nothing
        refused = scope_refusal(constellation_args(control=(ScopedControl(BoundSpacecraft((2, 4, 5)), 3),)))
        @test occursin("ScopedControl names spacecraft 4, 5", refused) && occursin("positions 1 to 3", refused)
        refused = scope_refusal(constellation_args(guidance=(ScopedGuidance(CoupledSpacecraft(1, (0,)), 3),)))
        @test occursin("ScopedGuidance names spacecraft 0", refused)
        # A refusal, not a warning, whatever --depwarn says: nothing is logged.
        logger = Test.TestLogger(min_level=Logging.Warn)
        refused = with_logger(logger) do
            scope_refusal(constellation_args(navigation=(ScopedNavigation(BoundSpacecraft((9,)), 3),)))
        end
        @test occursin("ScopedNavigation names spacecraft 9", refused) && isempty(logger.logs)
        # Malformed declarations are refused.
        for bad in (ScopedControl(:chaser, 3), ScopedControl(AllSpacecraft(), 3; skip=:yes))
            err = try
                SM.SimulationCallbacks.check_spacecraft_scopes(constellation_args(control=(bad,)))
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
        end
        # The engine runs the check before every run.
        refused = refusal(() -> run_simulation(constellation_args(n=2,
            control=(ScopedControl(BoundSpacecraft((3,)), 2),))))
        @test occursin("names spacecraft 3", refused) && occursin("positions 1 to 2", refused)
    end

    @testset "the built-in models declare their scope and held command" begin
        mpc = SM.RPOMPCControlModel(chaser_index=2, target_index=1)
        @test spacecraft_scope(mpc) == CoupledSpacecraft(2, (1,))
        mpc.held = SM.RPOHeldActuation(force_ii=SVector(1.0, 2.0, 3.0))
        @test held_command(mpc, 2) == (force_ii=SVector(1.0, 2.0, 3.0), torque_body=ZERO3,
            thruster_forces_n=zero(SVector{6, Float64}), rw_torque_body=ZERO3)
        @test held_command(mpc, 1) === nothing && isnan(held_command_time(mpc, 2))

        arm = SM.RobotArmControlEffector(spacecraft_index=2)
        @test spacecraft_scope(arm) == BoundSpacecraft((2,))
        arm.held.joint_torque_nm = [0.5, -0.5]
        held = held_command(arm, 2)
        @test held.joint_torque_nm == [0.5, -0.5] && held.joint_torque_nm !== arm.held.joint_torque_nm
        @test held_command(arm, 1) === nothing

        mm = SM.MagneticMomentumManagerModel(spacecraft_index=1)
        @test spacecraft_scope(mm) == BoundSpacecraft((1,))
        mm.held_torque_body = SVector(0.0, 0.0, 1e-3)
        mm.last_update_s = 4.0
        @test held_command(mm, 1) == (torque_body=SVector(0.0, 0.0, 1e-3), dipole_am2=ZERO3)
        @test held_command_time(mm, 1) == 4.0 && isnan(held_command_time(mm, 2))

        thruster = SM.BaseThrusterModel(thrust=[0.0, 2.0], direction=[0.0, 0.0], Δv=[0.0, 0.0],
            start_burn_time=[0.0, 0.0], stop_burn_time=[0.0, 0.0], Isp=[300.0, 300.0])
        @test spacecraft_scope(thruster) === AllSpacecraft()
        thruster.burn_windows[2] = (10.0, 20.0, 2.0)
        thruster.burning[2] = true
        @test held_command(thruster, 2) == (window_start_s=10.0, window_stop_s=20.0, thrust_n=2.0, burning=true)
        none = held_command(thruster, 1)
        @test isnan(none.window_start_s) && isnan(none.window_stop_s) && !none.burning
        @test held_command(thruster, 3) === nothing

        edg_state = SM.AerobrakingEnergyDepletionState(num_sats=2)
        edg = SM.AerobrakingEnergyDepletionControlModel(SM.AerobrakingEnergyDepletionConfig(), edg_state)
        @test spacecraft_scope(edg) === AllSpacecraft()
        edg_state.last_alpha_rad[2] = 0.25
        @test held_command(edg, 2) == 0.25 && held_command(edg, 3) === nothing

        binding = AlgorithmGuidance(algorithm=ApoapsisBurnPlanner(delta_v_mps=1.0), spacecraft_index=2)
        @test spacecraft_scope(binding) == BoundSpacecraft((2,))

        rpo = make_rpo_configuration(planner=DirectRPOPlanner(), mission_time_s=0.2)
        planner, control = only(rpo.guidance_model.guidance_effectors), only(rpo.control_model.control_effectors)
        @test spacecraft_scope(planner) == spacecraft_scope(control) == CoupledSpacecraft(1, (2,))
        # No built-in model ticks for an inactive spacecraft.
        edg_guidance = SM.AerobrakingEnergyDepletionGuidanceModel(SM.AerobrakingEnergyDepletionConfig(), edg_state)
        maneuver = SM.AerobrakingCampaignPropulsiveManeuverGuidanceModel(maneuver_orbit_number=[1], maneuver_Δv=[1.0])
        raise = SM.ApoapsisTargetPeriapsisRaiseGuidanceModel(target_apoapsis_radius_m=7.0e6)
        @test all(skip_inactive, (mpc, arm, mm, thruster, edg, SM.SolarPanelAngleOfAttackControlModel(), binding,
            planner, control, edg_guidance, maneuver, raise))
        @test scope_refusal(rpo) === nothing
    end

    @testset "the binding refusal names the scope" begin
        message(effector) = try
            bind_spacecraft(effector, 1)
            ""
        catch e
            e isa ArgumentError ? sprint(showerror, e) : rethrow()
        end
        @test occursin("couples several spacecraft", message(SM.RPOMPCControlModel()))
        @test occursin("acts on every spacecraft", message(PlainControl()))
        @test occursin("bound to spacecraft 2", message(ScopedControl(BoundSpacecraft((2,)), 2)))
    end
end
end
