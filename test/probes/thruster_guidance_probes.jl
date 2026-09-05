using Test
using Dates
using StaticArrays
using ComponentArrays
using LinearAlgebra

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
using .SimulationModel
include(joinpath(REPO_ROOT, "src", "core", "interfaces", "reference_system.jl"))

if !isdefined(@__MODULE__, :SimulationEngine)
    include(joinpath(REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
end
if !isdefined(@__MODULE__, :build_initial_conditions)
    const build_initial_conditions = SimulationEngine.build_initial_conditions
end

const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
const EARTH = Earth("", SPICE_PATH)

const GH = SimulationModel.GuidanceHooks
const CT = SimulationModel.CommandTypes

function make_probe_spacecraft(;
    ra_alt_m::Float64=500e3,
    rp_alt_m::Float64=100e3,
    ν_deg::Float64=160.0
)
    root = Link(root=true, m=500.0, ref_area=12.0)
    ic = InitialCondition(
        ra=EARTH.Rp_e + ra_alt_m,
        rp=EARTH.Rp_e + rp_alt_m,
        i=0.0,
        ω=0.0,
        Ω=0.0,
        ν=ν_deg
    )
    return SpacecraftModel(
        Joint[],
        [root],
        root,
        true,
        root.m,
        0.0,
        root.inertia,
        0,
        0,
        ic,
        1
    )
end

function make_probe_config(;
    n_sc::Int=1,
    guidance_effecters::Tuple=(),
    guidance_rates::Vector{Float64}=Float64[]
)
    environment_model = EnvironmentModel(
        planet=EARTH,
        EI=120.0,
        density_model=NoAtmosphereModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
        topography=false,
        wind=false
    )
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=60.0,
            orientation_sim=false,
            num_steps_to_save=100
        ),
        environment_model=environment_model,
        dynamics_model=DynamicsModel(
            [make_probe_spacecraft() for _ in 1:n_sc],
            (InverseSquaredGravityModel(),)
        ),
        guidance_model=GuidanceModel(guidance_effecters, guidance_rates),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances()
    )
end

# Build an equatorial (i = 0, ω = 0, Ω = 0) inertial state at true anomaly ν
# from apoapsis/periapsis radii, independently of the guidance helper code.
function equatorial_state(ra_m::Float64, rp_m::Float64, ν::Float64; μ::Float64=EARTH.μ)
    a = 0.5 * (ra_m + rp_m)
    e = (ra_m - rp_m) / (ra_m + rp_m)
    p_slr = a * (1.0 - e^2)
    r = p_slr / (1.0 + e * cos(ν))
    pos = SVector{3, Float64}(r * cos(ν), r * sin(ν), 0.0)
    vel = sqrt(μ / p_slr) .* SVector{3, Float64}(-sin(ν), e + cos(ν), 0.0)
    return pos, vel
end

function set_state!(u, i::Int, pos::SVector{3, Float64}, vel::SVector{3, Float64})
    u.sc[i].pos .= pos
    u.sc[i].vel .= vel
    return nothing
end

@testset "Thruster Guidance Coverage Probes" begin
    @testset "_wrap_2pi_guidance" begin
        @test GH._wrap_2pi_guidance(0.0) == 0.0
        @test GH._wrap_2pi_guidance(Float64(pi)) ≈ Float64(pi)
        @test GH._wrap_2pi_guidance(3.0 * pi) ≈ Float64(pi)
        @test GH._wrap_2pi_guidance(2pi) ≈ 0.0 atol=1e-15
        # Julia's mod already returns a value in [0, 2pi) for negative input,
        # so the wrap must land in [0, 2pi) here too.
        @test GH._wrap_2pi_guidance(-pi / 2) ≈ 3pi / 2
        @test 0.0 <= GH._wrap_2pi_guidance(-7.5) < 2pi
    end

    @testset "_ensure_apo_target_state!" begin
        model = ApoapsisTargetPeriapsisRaiseGuidanceModel(target_apoapsis_radius_m=1.0e7)
        @test isempty(model.command_state)
        GH._ensure_apo_target_state!(model, 3)
        @test model.command_state == zeros(Int64, 3)
        # Growing to a smaller index is a no-op.
        GH._ensure_apo_target_state!(model, 1)
        @test length(model.command_state) == 3
    end

    @testset "_osculating_elements_and_periapsis_direction" begin
        μ = EARTH.μ
        ra_m = EARTH.Rp_e + 500e3
        rp_m = EARTH.Rp_e + 100e3
        a_ref = 0.5 * (ra_m + rp_m)
        e_ref = (ra_m - rp_m) / (ra_m + rp_m)

        # Periapsis state: ν = 0, periapsis direction = +x.
        pos_p, vel_p = equatorial_state(ra_m, rp_m, 0.0)
        el_p = GH._osculating_elements_and_periapsis_direction(pos_p, vel_p, EARTH)
        @test el_p !== nothing
        @test el_p.a ≈ a_ref rtol=1e-10
        @test el_p.e ≈ e_ref rtol=1e-8
        @test el_p.true_anomaly ≈ 0.0 atol=1e-6
        @test el_p.periapsis_direction ≈ SVector(1.0, 0.0, 0.0) atol=1e-8

        # Ascending leg (dot(pos, vel) > 0): ν recovered directly.
        ν_asc = deg2rad(70.0)
        pos_a, vel_a = equatorial_state(ra_m, rp_m, ν_asc)
        @test dot(pos_a, vel_a) > 0.0
        el_a = GH._osculating_elements_and_periapsis_direction(pos_a, vel_a, EARTH)
        @test el_a !== nothing
        @test el_a.true_anomaly ≈ ν_asc rtol=1e-8

        # Descending leg (dot(pos, vel) < 0): ν reflected to (π, 2π).
        ν_desc = deg2rad(250.0)
        pos_d, vel_d = equatorial_state(ra_m, rp_m, ν_desc)
        @test dot(pos_d, vel_d) < 0.0
        el_d = GH._osculating_elements_and_periapsis_direction(pos_d, vel_d, EARTH)
        @test el_d !== nothing
        @test el_d.true_anomaly ≈ ν_desc rtol=1e-8

        # Near-circular orbit: e <= 1e-10 branch pins ν = 0 and uses -pos/r.
        r_c = EARTH.Rp_e + 400e3
        v_c = sqrt(μ / r_c)
        pos_c = SVector(r_c, 0.0, 0.0)
        vel_c = SVector(0.0, v_c, 0.0)
        el_c = GH._osculating_elements_and_periapsis_direction(pos_c, vel_c, EARTH)
        @test el_c !== nothing
        @test el_c.e <= 1e-10
        @test el_c.true_anomaly == 0.0
        @test el_c.periapsis_direction ≈ SVector(-1.0, 0.0, 0.0) atol=1e-12

        # Degenerate inputs return nothing.
        @test GH._osculating_elements_and_periapsis_direction(
            SVector(0.0, 0.0, 0.0), vel_c, EARTH) === nothing            # rmag = 0
        @test GH._osculating_elements_and_periapsis_direction(
            pos_c, SVector(2v_c, 0.0, 0.0), EARTH) === nothing           # h = 0 (radial velocity)
        @test GH._osculating_elements_and_periapsis_direction(
            pos_c, SVector(0.0, 2v_c, 0.0), EARTH) === nothing           # hyperbolic (energy >= 0)
        @test GH._osculating_elements_and_periapsis_direction(
            pos_c, SVector(NaN, v_c, 0.0), EARTH) === nothing            # non-finite state

        # Numeric-consistency guard: a huge radius makes the specific energy
        # cancel down to a denormal, so a = -μ/(2ε) overflows to Inf even
        # though the earlier finiteness gates pass.
        r_huge = 4.0e294
        v_esc = sqrt(2.0 * μ / r_huge)
        overflow_state = nothing
        for k in 1:200
            v_k = prevfloat(v_esc, k)
            en_k = 0.5 * v_k^2 - μ / r_huge
            if en_k < 0.0 && !isfinite(-μ / (2.0 * en_k))
                overflow_state = (SVector(r_huge, 0.0, 0.0), SVector(0.0, v_k, 0.0))
                break
            end
        end
        @test overflow_state !== nothing
        if overflow_state !== nothing
            @test GH._osculating_elements_and_periapsis_direction(
                overflow_state[1], overflow_state[2], EARTH) === nothing
        end
    end

    @testset "oblate altitude helpers" begin
        # Surface radius: equator -> Rp_e, pole -> Rp_p.
        @test GH._oblate_surface_radius(SVector(1.0, 0.0, 0.0), EARTH) ≈ EARTH.Rp_e rtol=1e-12
        @test GH._oblate_surface_radius(SVector(0.0, 0.0, 1.0), EARTH) ≈ EARTH.Rp_p rtol=1e-12
        u45 = SVector(sqrt(0.5), 0.0, sqrt(0.5))
        r45 = GH._oblate_surface_radius(u45, EARTH)
        @test EARTH.Rp_p < r45 < EARTH.Rp_e

        # Geodetic altitude from radius: exact along equator and pole.
        @test GH._oblate_altitude_from_radius(EARTH.Rp_e + 250e3, SVector(1.0, 0.0, 0.0), EARTH) ≈ 250e3 atol=1e-3
        @test GH._oblate_altitude_from_radius(EARTH.Rp_p + 250e3, SVector(0.0, 0.0, 1.0), EARTH) ≈ 250e3 atol=1e-3
        # Bowring closed form: a point ON the ellipsoid at 45 deg latitude
        # must report ~0 m of geodetic altitude, not just at the equator and
        # pole where the e2*N*sin(lat) term vanishes.
        surf45_alt = GH._oblate_altitude_from_radius(r45, u45, EARTH)
        @test abs(surf45_alt) < 1.0

        # Radius-for-altitude inversion round-trips through the altitude model.
        @test isnan(GH._radius_for_oblate_altitude(-1.0, SVector(1.0, 0.0, 0.0), EARTH))
        r_eq = GH._radius_for_oblate_altitude(300e3, SVector(1.0, 0.0, 0.0), EARTH)
        @test r_eq ≈ EARTH.Rp_e + 300e3 atol=1.0
        r_mid = GH._radius_for_oblate_altitude(200e3, u45, EARTH)
        @test GH._oblate_altitude_from_radius(r_mid, u45, EARTH) ≈ 200e3 atol=1.0
        # Zero target altitude converges onto the surface.
        r_zero = GH._radius_for_oblate_altitude(0.0, u45, EARTH)
        @test abs(GH._oblate_altitude_from_radius(r_zero, u45, EARTH)) < 1.0

        # Strongly oblate synthetic body: the initial bracket upper bound can
        # sit below the target geodetic altitude, forcing the hi-expansion
        # loop to run before bisection. The helpers only touch Rp_e/Rp_p, so
        # a NamedTuple stands in for the planet.
        squashed = (Rp_e=10.0, Rp_p=1.0)
        u_low_lat = SVector(cosd(10.0), 0.0, sind(10.0))
        target_alt = 10.0
        lo0 = GH._oblate_surface_radius(u_low_lat, squashed)
        hi0 = lo0 + target_alt + abs(squashed.Rp_e - squashed.Rp_p) + 1.0
        @test GH._oblate_altitude_from_radius(hi0, u_low_lat, squashed) < target_alt
        r_squashed = GH._radius_for_oblate_altitude(target_alt, u_low_lat, squashed)
        @test r_squashed > hi0
        @test GH._oblate_altitude_from_radius(r_squashed, u_low_lat, squashed) ≈ target_alt atol=1e-8
    end

    @testset "AerobrakingCampaignPropulsiveManeuverGuidanceModel calcGuidanceEffect!" begin
        campaign = AerobrakingCampaignPropulsiveManeuverGuidanceModel(
            maneuver_orbit_number=[2, 3, 4],
            maneuver_Δv=[5.0, -3.0, 0.0]
        )
        args = make_probe_config(n_sc=2, guidance_effecters=(campaign,), guidance_rates=[1.0])
        p = ODEParams(n_sats=2, args=args)
        u = build_initial_conditions(args)

        # Out-of-range spacecraft indices are ignored.
        @test calcGuidanceEffect!(campaign, u, p, 0.0, 0) === nothing
        @test calcGuidanceEffect!(campaign, u, p, 0.0, 3) === nothing
        @test all(cmd -> !cmd.valid && cmd.delta_v_mps == 0.0, p.shared_buffers.maneuver_commands)

        # Orbit with no scheduled maneuver: valid command with zero delta-v.
        p.orbit_counter[1] = 1
        calcGuidanceEffect!(campaign, u, p, 0.0, 1)
        cmd_idle = p.shared_buffers.maneuver_commands[1]
        @test cmd_idle.valid
        @test cmd_idle.delta_v_mps == 0.0
        @test cmd_idle.direction_rad == 0.0
        @test cmd_idle.source_orbit == 1

        # Positive scheduled delta-v: prograde command (direction 0).
        p.orbit_counter[1] = 2
        calcGuidanceEffect!(campaign, u, p, 10.0, 1)
        cmd_pro = p.shared_buffers.maneuver_commands[1]
        @test cmd_pro.valid
        @test cmd_pro.delta_v_mps == 5.0
        @test cmd_pro.direction_rad == 0.0
        @test cmd_pro.source_orbit == 2

        # Negative scheduled delta-v: magnitude with retrograde direction (π).
        p.orbit_counter[2] = 3
        calcGuidanceEffect!(campaign, u, p, 10.0, 2)
        cmd_retro = p.shared_buffers.maneuver_commands[2]
        @test cmd_retro.valid
        @test cmd_retro.delta_v_mps == 3.0
        @test cmd_retro.direction_rad ≈ Float64(π)
        @test cmd_retro.source_orbit == 3

        # Zero scheduled delta-v falls into the retrograde direction branch.
        p.orbit_counter[2] = 4
        calcGuidanceEffect!(campaign, u, p, 20.0, 2)
        cmd_zero = p.shared_buffers.maneuver_commands[2]
        @test cmd_zero.valid
        @test cmd_zero.delta_v_mps == 0.0
        @test cmd_zero.direction_rad ≈ Float64(π)
        @test cmd_zero.source_orbit == 4
    end

    @testset "ApoapsisTargetPeriapsisRaiseGuidanceModel calcGuidanceEffect!" begin
        μ = EARTH.μ
        ra_m = EARTH.Rp_e + 500e3
        rp_m = EARTH.Rp_e + 100e3
        args = make_probe_config(n_sc=1)
        eph = args.environment_model.ephemerides_model

        # ---- Happy path: pre-apoapsis inside window, apoapsis below target ----
        model = ApoapsisTargetPeriapsisRaiseGuidanceModel(
            target_apoapsis_radius_m=EARTH.Rp_e + 600e3,
            target_periapsis_altitude_m=200e3
        )
        p = ODEParams(n_sats=1, args=args)
        u = build_initial_conditions(args)
        ν_burn = deg2rad(160.0)
        pos, vel = equatorial_state(ra_m, rp_m, ν_burn)
        set_state!(u, 1, pos, vel)
        p.orbit_counter[1] = 7

        @test calcGuidanceEffect!(model, u, p, 0.0, 1) === nothing
        cmd = p.shared_buffers.maneuver_commands[1]
        @test cmd.valid
        @test cmd.direction_rad == 0.0
        @test cmd.source_orbit == 7
        @test model.command_state[1] == GH._APO_TARGET_COMMAND_ISSUED

        # Independent delta-v check: apoapsis-tangential burn onto the target
        # (r_apo, r_peri_target) ellipse. Periapsis direction is +x (equatorial
        # orbit built with ω = Ω = i = 0).
        a_cur = 0.5 * (ra_m + rp_m)
        l_pi = GH.planet_frame_lpi(EARTH, p.shared_buffers.et_start[], eph)
        u_pp = SVector{3, Float64}(l_pi * SVector(1.0, 0.0, 0.0))
        u_pp = u_pp / norm(u_pp)
        rp_target = GH._radius_for_oblate_altitude(200e3, u_pp, EARTH)
        a_target = 0.5 * (ra_m + rp_target)
        dv_expected = sqrt(μ * (2.0 / ra_m - 1.0 / a_target)) - sqrt(μ * (2.0 / ra_m - 1.0 / a_cur))
        @test dv_expected > 0.0
        @test cmd.delta_v_mps ≈ dv_expected rtol=1e-6
        @test 5.0 < cmd.delta_v_mps < 100.0

        # ---- COMMAND_ISSUED with no valid burn plan: command left in place ----
        calcGuidanceEffect!(model, u, p, 1.0, 1)
        @test p.shared_buffers.maneuver_commands[1].valid
        @test model.command_state[1] == GH._APO_TARGET_COMMAND_ISSUED

        # ---- COMMAND_ISSUED with a valid burn plan: command retired ----
        p.shared_buffers.maneuver_burn_plans[1] = CT.PropulsiveBurnPlan(valid=true, delta_v_mps=cmd.delta_v_mps)
        calcGuidanceEffect!(model, u, p, 2.0, 1)
        @test !p.shared_buffers.maneuver_commands[1].valid
        @test model.command_state[1] == GH._APO_TARGET_DISCARDED

        # ---- DISCARDED state: command stays cleared, no re-issue ----
        p.shared_buffers.maneuver_commands[1] = CT.PropulsiveManeuverCommand(valid=true, delta_v_mps=99.0)
        calcGuidanceEffect!(model, u, p, 3.0, 1)
        @test !p.shared_buffers.maneuver_commands[1].valid
        @test model.command_state[1] == GH._APO_TARGET_DISCARDED

        # ---- Out-of-range index leaves state untouched ----
        model_oob = ApoapsisTargetPeriapsisRaiseGuidanceModel(target_apoapsis_radius_m=EARTH.Rp_e + 600e3)
        @test calcGuidanceEffect!(model_oob, u, p, 0.0, 0) === nothing
        @test calcGuidanceEffect!(model_oob, u, p, 0.0, 2) === nothing
        @test isempty(model_oob.command_state)

        # ---- Degenerate state: osculating-element extraction fails, no command ----
        model_bad = ApoapsisTargetPeriapsisRaiseGuidanceModel(target_apoapsis_radius_m=EARTH.Rp_e + 600e3)
        p_bad = ODEParams(n_sats=1, args=args)
        u_bad = build_initial_conditions(args)
        set_state!(u_bad, 1, pos, SVector(0.0, 0.0, 0.0))
        calcGuidanceEffect!(model_bad, u_bad, p_bad, 0.0, 1)
        @test !p_bad.shared_buffers.maneuver_commands[1].valid
        @test model_bad.command_state[1] == GH._APO_TARGET_IDLE

        # ---- Non-finite target apoapsis radius: rejected before windowing ----
        model_naninf = ApoapsisTargetPeriapsisRaiseGuidanceModel(target_apoapsis_radius_m=Inf)
        p_naninf = ODEParams(n_sats=1, args=args)
        u_naninf = build_initial_conditions(args)
        set_state!(u_naninf, 1, pos, vel)
        calcGuidanceEffect!(model_naninf, u_naninf, p_naninf, 0.0, 1)
        @test !p_naninf.shared_buffers.maneuver_commands[1].valid
        @test model_naninf.command_state[1] == GH._APO_TARGET_IDLE

        # ---- Apoapsis already above target (outside tolerance): no command ----
        model_high = ApoapsisTargetPeriapsisRaiseGuidanceModel(
            target_apoapsis_radius_m=EARTH.Rp_e + 300e3,
            apoapsis_tolerance_m=0.0
        )
        p_high = ODEParams(n_sats=1, args=args)
        u_high = build_initial_conditions(args)
        set_state!(u_high, 1, pos, vel)
        calcGuidanceEffect!(model_high, u_high, p_high, 0.0, 1)
        @test !p_high.shared_buffers.maneuver_commands[1].valid
        @test model_high.command_state[1] == GH._APO_TARGET_IDLE

        # ---- Same geometry accepted once the tolerance absorbs the excess ----
        model_tol = ApoapsisTargetPeriapsisRaiseGuidanceModel(
            target_apoapsis_radius_m=EARTH.Rp_e + 300e3,
            apoapsis_tolerance_m=250e3,
            target_periapsis_altitude_m=200e3
        )
        p_tol = ODEParams(n_sats=1, args=args)
        u_tol = build_initial_conditions(args)
        set_state!(u_tol, 1, pos, vel)
        calcGuidanceEffect!(model_tol, u_tol, p_tol, 0.0, 1)
        @test p_tol.shared_buffers.maneuver_commands[1].valid
        @test model_tol.command_state[1] == GH._APO_TARGET_COMMAND_ISSUED

        # ---- Pre-apoapsis but outside the anomaly window: no command ----
        model_early = ApoapsisTargetPeriapsisRaiseGuidanceModel(target_apoapsis_radius_m=EARTH.Rp_e + 600e3)
        p_early = ODEParams(n_sats=1, args=args)
        u_early = build_initial_conditions(args)
        pos_e, vel_e = equatorial_state(ra_m, rp_m, deg2rad(90.0))
        set_state!(u_early, 1, pos_e, vel_e)
        calcGuidanceEffect!(model_early, u_early, p_early, 0.0, 1)
        @test !p_early.shared_buffers.maneuver_commands[1].valid
        @test model_early.command_state[1] == GH._APO_TARGET_IDLE

        # ---- Post-apoapsis: no command even close to apoapsis ----
        model_late = ApoapsisTargetPeriapsisRaiseGuidanceModel(target_apoapsis_radius_m=EARTH.Rp_e + 600e3)
        p_late = ODEParams(n_sats=1, args=args)
        u_late = build_initial_conditions(args)
        pos_l, vel_l = equatorial_state(ra_m, rp_m, deg2rad(185.0))
        set_state!(u_late, 1, pos_l, vel_l)
        calcGuidanceEffect!(model_late, u_late, p_late, 0.0, 1)
        @test !p_late.shared_buffers.maneuver_commands[1].valid
        @test model_late.command_state[1] == GH._APO_TARGET_IDLE

        # ---- Target periapsis radius at/above current apoapsis: infeasible ----
        model_infeasible = ApoapsisTargetPeriapsisRaiseGuidanceModel(
            target_apoapsis_radius_m=EARTH.Rp_e + 600e3,
            target_periapsis_altitude_m=5_000e3
        )
        p_inf = ODEParams(n_sats=1, args=args)
        u_inf = build_initial_conditions(args)
        set_state!(u_inf, 1, pos, vel)
        calcGuidanceEffect!(model_infeasible, u_inf, p_inf, 0.0, 1)
        @test !p_inf.shared_buffers.maneuver_commands[1].valid
        @test model_infeasible.command_state[1] == GH._APO_TARGET_IDLE

        # ---- Non-positive delta-v (target periapsis below current): discard ----
        model_lower = ApoapsisTargetPeriapsisRaiseGuidanceModel(
            target_apoapsis_radius_m=EARTH.Rp_e + 600e3,
            target_periapsis_altitude_m=0.0
        )
        p_low = ODEParams(n_sats=1, args=args)
        u_low = build_initial_conditions(args)
        set_state!(u_low, 1, pos, vel)
        p_low.shared_buffers.maneuver_commands[1] = CT.PropulsiveManeuverCommand(valid=true, delta_v_mps=1.0)
        calcGuidanceEffect!(model_lower, u_low, p_low, 0.0, 1)
        @test !p_low.shared_buffers.maneuver_commands[1].valid
        @test model_lower.command_state[1] == GH._APO_TARGET_DISCARDED
        # Once discarded, later passes keep the command cleared.
        calcGuidanceEffect!(model_lower, u_low, p_low, 1.0, 1)
        @test !p_low.shared_buffers.maneuver_commands[1].valid
    end
end
