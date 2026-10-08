# Bounded regression fixtures for the Gen-EDG repair candidate.
# Execute exact production definitions in isolated modules. Dependencies are deterministic
# fixtures, so this does not qualify the complete simulation, shooting model or native images.
module EDGE2cRepairTests
using Test, StaticArrays, LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const HEAT = "src/gnc/control/heat_load_control.jl"
const CONTROL = "src/gnc/control/targeting_control.jl"
const GUIDANCE = "src/gnc/guidance/target_energy_bracketing.jl"

function load_definitions!(target, path, names)
    # The override allows the same tests to discriminate the pinned parent revision.
    revision = get(ENV, "EDG_REPAIR_SOURCE_REVISION", "")
    source = isempty(revision) ? read(joinpath(ROOT, path), String) :
        read(`git -C $ROOT show $revision:$path`, String)
    definitions = Meta.parseall(source)
    for name in names
        found = Any[]
        for item in definitions.args
            item isa Expr && item.head == :function || continue
            signature = item.args[1]
            signature isa Expr && signature.head == :(::) && (signature = signature.args[1])
            signature isa Expr && signature.head == :call && signature.args[1] == name && push!(found, item)
        end
        length(found) == 1 || error("Expected exactly one production definition of $name in $path")
        Core.eval(target, only(found))
    end
end

module HeatFixture
using Test, StaticArrays, LinearAlgebra, Roots
struct ODEParams
    args
end
const TRACK = Ref((time=[0.0, 30.0, 60.0], rates=[1.0, 1.0, 1.0]))
const EXPECTED_PROFILE = Ref(0.0)
_edg_predict_mass(spacecraft, mass) = mass
_edg_closed_form_heat_load_trajectory(args...) = TRACK[]
function _edg_profile_heat_rates(config, p, track, profile; heat_rate_control)
    @test !heat_rate_control
    @test profile == fill(EXPECTED_PROFILE[], length(track.time))
    return copy(track.rates)
end
const ROOT_K = Ref(2.0)
const RESIDUAL_MODE = Ref(:crossing)
const K_SEEN = Float64[]
_edg_heat_load_coefficients(args...) = nothing
function _edg_heat_load_profile_for_k(config, p, spacecraft, pos, vel, mass, t, env, coeffs,
        k, links, heat_rate_control, structural_control, shooting_guess)
    push!(K_SEEN, k)
    residual = RESIDUAL_MODE[] == :crossing ? k - ROOT_K[] :
        RESIDUAL_MODE[] == :below ? -10.0 : 10.0
    # Integral equals limit + residual; the exact production quadrature is exercised.
    track = (time=[0.0, 1.0], rates=fill(config.heat_load_limit_j_cm2 + residual, 2))
    return track, [k, k], [k, k], nothing
end
# Encoding the solved k in the window exposes the actual root returned by Brent.
_edg_low_alpha_switch_window(config, t, track, profile) = (t + profile[1], t + profile[1] + 1.0)
end
load_definitions!(HeatFixture, HEAT, [:_edg_profile_heat_load, :_edg_heat_load_security_required])

# Root fixtures need variable profiles, independently from the all-minimum security fixture.
module RootFixture
using Test, StaticArrays, LinearAlgebra, Roots
using ..HeatFixture: ODEParams, ROOT_K, RESIDUAL_MODE, K_SEEN, _edg_predict_mass,
    _edg_heat_load_coefficients, _edg_heat_load_profile_for_k, _edg_low_alpha_switch_window
_edg_profile_heat_rates(config, p, track, profile; heat_rate_control) = copy(track.rates)
end
load_definitions!(RootFixture, HEAT, [:_edg_profile_heat_load, :_edg_solve_heat_load_switches])

module ScheduleFixture
using StaticArrays, LinearAlgebra
struct ODEParams
    args
end
struct AerobrakingEnergyDepletionControlModel
    config
    state
end
const CALLS = Ref(0)
_edg_in_drag_passage(p, env) = true
function _edg_recompute_second_heat_load_switch(config, p, spacecraft, pos, vel, mass, env, load, t, switches; kwargs...)
    CALLS[] += 1
    return (switches[1], switches[2] + 1.0)
end
end
load_definitions!(ScheduleFixture, CONTROL, [:_edg_recompute_switches!])

module VacuumFixture
using StaticArrays, LinearAlgebra
struct ODEParams
    args
end
const FORCE = Ref(:central)
function _edg_prediction_gravity_acceleration(p, r, v, mass, t)
    FORCE[] == :central && return -p.args.environment_model.planet.μ * r / norm(r)^3
    FORCE[] == :turn_at_cap && return SVector(-1.0, 0.0, 0.0)
    return SVector(0.0, 0.0, 0.0)
end
function _edg_orbit_metrics_from_rv(r, v, mass, planet)
    energy = dot(v, v) / 2 - planet.μ / norm(r)
    if !(isfinite(energy) && energy < 0)
        return (energy=energy, periapsis=NaN, apoapsis=Inf)
    end
    a = -planet.μ / (2energy)
    eccentricity = norm(cross(v, cross(r, v)) / planet.μ - r / norm(r))
    return (energy=energy, periapsis=a*(1-eccentricity), apoapsis=a*(1+eccentricity))
end
_edg_ephemeris_time(p, t) = t
r_intor_p!(r, v, planet, t, ephemerides) = (r, v)
rtolatlong(r, planet) = (norm(r) - planet.Rp, 0.0, 0.0)
end
load_definitions!(VacuumFixture, CONTROL, [:_edg_vacuum_drag_passage_exit, :_edg_vacuum_apoapsis_correction])

module GuidanceFixture
using StaticArrays
struct ODEParams
    args
end
struct AerobrakingEnergyDepletionGuidanceModel
    config
    state
end
const EXIT_OK = Ref(true)
const APO_OK = Ref(true)
const APO_CALLS = Ref(0)
_control_module() = @__MODULE__
_edg_environment_state(args...) = nothing
_edg_in_drag_passage(args...) = true
_edg_sat_state(u, i) = u
_edg_pos_vel_mass(u) = (u.pos, u.vel, u.mass)
_edg_max_heat_load_for_links(args...) = 0.0
_edg_pass_heat_load_for_links(args...) = 0.0
_edg_targeting_bracket_outcomes(args...; kwargs...) =
    ((energy_jkg=-20.0, periapsis_radius_m=2.0), (energy_jkg=-10.0, periapsis_radius_m=2.0))
_edg_vacuum_drag_passage_exit(p, pos, vel, mass, t) =
    (position=pos, velocity=vel, propagation_time_s=2000.0, event_reached=EXIT_OK[])
function _edg_vacuum_apoapsis_correction(args...)
    APO_CALLS[] += 1
    return (periapsis_radius_m=2.0, energy_change_jkg=0.0,
        propagation_time_s=20000.0, event_reached=APO_OK[])
end
_edg_corrected_target_energy_from_apoapsis(args...) = -15.0
end
load_definitions!(GuidanceFixture, GUIDANCE, [:_edg_run_target_energy_bracketing!])

const POS = SVector(1.0, 0.0, 0.0)
const VEL = SVector(1.0, 1.0, 0.0)

@testset verbose=true "Gen-EDG E2c regressions" begin
@testset "E2c repairs: security remainder" begin
    F = HeatFixture
    config = (min_alpha_rad=0.0, heat_load_limit_j_cm2=100.0)
    for temperature in (1.0, 150.0, 300.0), t in (0.0, 5000.0)
        p = F.ODEParams((environment_model=(planet=(T=temperature,),),))
        F.TRACK[] = (time=[0.0, 30.0, 60.0], rates=[1.0, 1.0, 1.0])
        @test F._edg_heat_load_security_required(config,p,nothing,POS,VEL,1.0,nothing,50.0,t) == (true,t+60.0)
        @test F._edg_heat_load_security_required(config,p,nothing,POS,VEL,1.0,nothing,40.0,t) == (false,t+60.0)
        # Heat concentrated near the current state must not be discarded.
        F.TRACK[] = (time=[0.0, 30.0, 60.0], rates=[3.0, 0.0, 0.0])
        @test F._edg_heat_load_security_required(config,p,nothing,POS,VEL,1.0,nothing,50.0,t)[1]
        F.TRACK[] = (time=[0.0, 30.0, 60.0], rates=zeros(3))
        @test !F._edg_heat_load_security_required(config,p,nothing,POS,VEL,1.0,nothing,99.0,t)[1]
    end
end

@testset "E2c repairs: verified Brent bracket" begin
    F = RootFixture
    for (planet,solver,root_k) in (("mars",:closed_form,2.0), ("venus",:closed_form,60.0),
            ("earth",:closed_form,5.0), ("mars",:tpbvp_integration,0.06),
            ("venus",:tpbvp_integration,0.4))
        config = (heat_load_limit_j_cm2=100.0, heat_load_switch_solver=solver, controlled_panel_links=())
        p = F.ODEParams((environment_model=(planet=(name=planet,),),))
        F.ROOT_K[]=root_k; F.RESIDUAL_MODE[]=:crossing; empty!(F.K_SEEN)
        result = F._edg_solve_heat_load_switches(config,p,nothing,POS,VEL,1.0,nothing,0.0,20.0;
            heat_rate_control=false,structural_control=false)
        @test result[1] ≈ 20.0 + root_k atol=1e-5
        @test result[2] - result[1] ≈ 1.0
        @test last(F.K_SEEN) ≈ root_k atol=1e-5
        for (mode,expected) in ((:below,(20.0,20.0)),
                (:above,solver == :closed_form ? (20.0,20.5) : (20.0,1020.0)))
            F.RESIDUAL_MODE[]=mode
            @test F._edg_solve_heat_load_switches(config,p,nothing,POS,VEL,1.0,nothing,0.0,20.0;
                heat_rate_control=false,structural_control=false) == expected
        end
    end
end

@testset "E2c repairs: closed-form cadence and numerical caching" begin
    F = ScheduleFixture
    # Cases straddle all existing strict cadence thresholds and include expired short windows.
    cases = [(60.0,11.0,true), (60.0,10.0,false), (40.0,3.1,true),
        (40.0,3.0,false), (2.0,2.0,true), (2.0,0.8,false),
        (0.0,0.0,false), (-1.0,4.0,false)]
    for solver in (:closed_form,:tpbvp_integration), enabled in (false,true),
            security in (false,true), (remaining,elapsed_since,closed_due) in cases
        t=60.0
        config=(heat_load_switch_solver=solver, second_switch_reevaluation=enabled,
            heat_load_security_mode=false, heat_load_limit_j_cm2=100.0, max_energy_submodes=(:heat_load,))
        state=(selected_mode=[:max_energy_depletion], targeting_active=[false],
            heat_load_drag_passage_active=[true], heat_load_switch_solved=[true],
            heat_load_switches_s=[(remaining < 0 ? t-2.0 : t-10.0,t+remaining)],
            heat_load_entry_time_s=[0.0],heat_load_last_reevaluation_s=[t-elapsed_since],
            heat_load_security_active=[security],heat_load_previous_j_cm2=[0.0])
        before=state.heat_load_switches_s[1]; F.CALLS[]=0
        model=F.AerobrakingEnergyDepletionControlModel(config,state)
        F._edg_recompute_switches!(model,F.ODEParams(nothing),nothing,nothing,POS,VEL,1.0,10.0,t,1)
        expected=solver == :closed_form && enabled && !security && closed_due
        @test F.CALLS[] == Int(expected)
        @test state.heat_load_switches_s[1] == (before[1],before[2]+Int(expected))
        @test state.heat_load_last_reevaluation_s[1] == (expected ? t : t-elapsed_since)
    end
end

@testset "E2c repairs: vacuum events and caps" begin
    F=VacuumFixture
    params(mu,radius,exit) = F.ODEParams((environment_model=(planet=(μ=mu,Rp=radius),EI=exit/1000,ephemerides_model=nothing),))
    F.FORCE[]=:zero
    p=params(1.0,1e6,2010.0)
    exit=F._edg_vacuum_drag_passage_exit(p,SVector(1e6+10,0.0,0.0),SVector(1.0,0.0,0.0),1.0,0.0)
    @test exit.event_reached
    @test exit.propagation_time_s == 2000.0 # Event exactly on cap is success.
    capped=F._edg_vacuum_drag_passage_exit(params(1.0,1e6,2011.0),SVector(1e6+10,0.0,0.0),SVector(1.0,0.0,0.0),1.0,0.0)
    @test !capped.event_reached
    @test capped.propagation_time_s == 2000.0
    descending=F._edg_vacuum_drag_passage_exit(params(1.0,1e6,0.0),SVector(2e6,0.0,0.0),SVector(-1.0,0.0,0.0),1.0,0.0)
    @test !descending.event_reached # High altitude alone is insufficient.
    unbound=F._edg_vacuum_apoapsis_correction(p,SVector(1e9,0.0,0.0),SVector(1.0,0.0,0.0),1.0,0.0)
    @test !unbound.event_reached
    @test unbound.propagation_time_s == 20000.0
    F.FORCE[]=:turn_at_cap
    turn=F._edg_vacuum_apoapsis_correction(p,SVector(1e9,0.0,0.0),SVector(20000.0,0.0,0.0),1.0,0.0)
    @test turn.event_reached
    @test turn.propagation_time_s == 20000.0
    F.FORCE[]=:central
    mu=3.986004418e14; rp=7e6; ra=9e6; a=(rp+ra)/2
    speed=sqrt(mu*(2/rp-1/a)); half_period=pi*sqrt(a^3/mu)
    orbit=F._edg_vacuum_apoapsis_correction(params(mu,6e6,100000.0),SVector(rp,0.0,0.0),SVector(0.0,speed,0.0),1.0,0.0)
    @test orbit.event_reached
    @test abs(orbit.propagation_time_s-half_period) <= 1.0
    @test orbit.apoapsis_radius_m ≈ ra rtol=1e-9
    @test abs(orbit.energy_change_jkg) < 0.01
end

@testset "E2c repairs: guidance refuses uncertified results" begin
    F=GuidanceFixture
    for exit_ok in (false,true), apo_ok in (false,true)
        F.EXIT_OK[]=exit_ok; F.APO_OK[]=apo_ok; F.APO_CALLS[]=0
        state=(energy_bracketing_evaluated=[false],energy_bracketing_count=[0],
            target_energy_jkg=[0.0],target_perturbation_energy_change_jkg=[0.0],
            bracket_min_energy_jkg=[0.0],bracket_max_energy_jkg=[0.0],
            targeting_active=[false],safe_low_drag=[false],selected_mode=[:unset])
        prior=deepcopy(state)
        config=(controlled_panel_links=(),max_energy_submodes=(),target_apoapsis_radius_m=10.0,guidance_modes=(:targeting,))
        model=F.AerobrakingEnergyDepletionGuidanceModel(config,state)
        p=F.ODEParams((dynamics_model=(spacecraft=[nothing],),environment_model=(planet=nothing,)))
        u=(pos=POS,vel=VEL,mass=1.0)
        if exit_ok && apo_ok
            F._edg_run_target_energy_bracketing!(model,u,p,10.0,1)
            @test state.energy_bracketing_evaluated == [true]
            @test state.energy_bracketing_count == [1]
            @test state.target_energy_jkg == [-15.0]
            @test state.selected_mode == [:targeting]
            F._edg_run_target_energy_bracketing!(model,u,p,11.0,1)
            @test F.APO_CALLS[] == 1 # Still cached after success.
        else
            err = try
                F._edg_run_target_energy_bracketing!(model,u,p,10.0,1)
                nothing
            catch e
                e
            end
            @test err isa ErrorException
            @test occursin(exit_ok ? "did not reach apoapsis" : "did not reach outbound EI", sprint(showerror,err))
            @test state == prior
            @test F.APO_CALLS[] == Int(exit_ok)
        end
    end
end
end # outer testset
end # module
