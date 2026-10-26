module ApolloDescentTests
using Test
using Logging
using SpaceAGORA
using SpaceAGORA.SimulationModel
using StaticArrays
using LinearAlgebra
import SpaceAGORA.TelemetryVerification: make_example_config

const SM = SpaceAGORA.SimulationModel
const GH = SM.GuidanceHooks
const CH = SM.ControlHooks

# A user's own guidance or control model that works through an ApolloDescentState.
struct CustomDescentGuidance <: SpaceAGORA.AbstractGuidanceModel
    state::ApolloDescentState
end
struct CustomDescentControl <: SpaceAGORA.AbstractControlEffectorModel
    state::ApolloDescentState
end

@testset "ApolloDescent" begin
    @testset "quadratic law" begin
        ph = DescentPhaseTargets(name=:test, r_T=SVector(0.0, 0.0, 30.0), v_T=SVector(0.0, 0.0, -1.0), a_T=SVector(0.0, 0.0, 0.0), t_go_initial_s=100.0, t_go_min_s=4.0)
        r = SVector(1000.0, 50.0, 500.0); v = SVector(-20.0, 0.0, -8.0); T = 100.0
        a = GH.descent_accel_command(r, v, ph, T)
        @test a ≈ ph.a_T - 6 * (ph.v_T + v) / T + 12 * (ph.r_T - r) / T^2
        # the cubic acceleration profile that starts at a_c reaches the targets at T
        # (integrate a(τ) = c0 + c1 τ + c2 τ^2 with the closed-form coefficients)
        A = ph.a_T * T - (ph.v_T - v); B = ph.a_T * T^2 / 2 - (ph.r_T - r - v * T)
        c2 = (24 * A * T - 36 * B) / T^4; c1 = (-30 * A + 48 * B / T) / T^2; c0 = ph.a_T - c1 * T - c2 * T^2
        @test c0 ≈ a
        rT = r + v * T + c0 * T^2 / 2 + c1 * T^3 / 6 + c2 * T^4 / 12
        vT = v + c0 * T + c1 * T^2 / 2 + c2 * T^3 / 3
        @test rT ≈ ph.r_T atol=1e-6
        @test vT ≈ ph.v_T atol=1e-9
        # time-to-go: a state exactly on a profile with zero jerk solves back to its T
        T_true = 80.0
        ph0 = DescentPhaseTargets(name=:t, r_T=SVector(0.0, 0.0, 0.0), v_T=SVector(-10.0, 0.0, 0.0), a_T=SVector(0.5, 0.0, 0.0), jerk_uprange=0.0, t_go_initial_s=80.0, t_go_min_s=1.0)
        # with zero jerk and zero snap: r(τ) = r_T + v_T τ + a_T τ²/2 for τ = -T
        r0 = ph0.r_T - ph0.v_T * T_true + ph0.a_T * T_true^2 / 2
        v0 = ph0.v_T - ph0.a_T * T_true
        @test GH.descent_time_to_go(r0, v0, ph0, 60.0) ≈ T_true atol=1e-3
        @test GH.descent_time_to_go(r0, v0, ph0, 100.0) ≈ T_true atol=1e-3
        cfg_thr = ApolloDescentConfig(reference_radius_m=1_737_400.0, site_lat_deg=0.0, site_lon_deg=0.0, braking=ph, approach=ph, max_thrust_n=1000.0)
        @test GH.descent_throttle(cfg_thr, 50.0) == 0.10
        @test GH.descent_throttle(cfg_thr, 400.0) == 0.40
        @test GH.descent_throttle(cfg_thr, 700.0) == 1.0
        b, ap = apollo11_descent_targets()
        @test b.name == :braking && ap.name == :approach
        @test b.altitude_switch_m == 2500.0
    end

    @testset "attitude command and error" begin
        up = SVector(0.0, 0.0, 1.0)
        # thrust straight up: body -z along +z, windows (+x) horizontal
        fwd = SVector(1.0, 0.0, 0.0)
        q = descent_attitude_command(up, up, fwd)
        A = SM.rot(q)
        @test SVector(A[3, 1], A[3, 2], A[3, 3]) ≈ -up atol=1e-12
        @test SVector(A[1, 1], A[1, 2], A[1, 3]) ≈ fwd atol=1e-12          # windows toward the site when hovering
        # thrust horizontal (retrograde at PDI): windows up
        dir = SVector(1.0, 0.0, 0.0)
        q2 = descent_attitude_command(dir, up, dir)
        A2 = SM.rot(q2)
        @test SVector(A2[3, 1], A2[3, 2], A2[3, 3]) ≈ -dir atol=1e-12
        @test SVector(A2[1, 1], A2[1, 2], A2[1, 3]) ≈ up atol=1e-12
        # the DCM round trip keeps the passive convention
        @test SM.rot(GH.quaternion_from_passive_dcm(A2)) ≈ A2 atol=1e-12
        @test SM.rot(GH.quaternion_from_passive_dcm(A)) ≈ A atol=1e-12
        # error vector: rotate the command by +10 deg about body y and check the turn back is -10 deg about y
        δ = deg2rad(10.0)
        q_rot = SVector(0.0, sin(δ / 2), 0.0, cos(δ / 2))         # rotation about body y
        q_actual = SM.quat_mult(q2, q_rot)                        # body rotated by δ relative to the command
        e = CH.attitude_error_vector(SVector{4, Float64}(q_actual), q2)
        @test e ≈ SVector(0.0, -δ, 0.0) atol=1e-9
        @test norm(CH.attitude_error_vector(q2, q2)) < 1e-12
        # Exact half turns have zero skew part in their rotation matrix; the
        # controller must still command a pi-radian turn rather than zero.
        identity = SVector(0.0, 0.0, 0.0, 1.0)
        half_turn = SVector(0.0, 1.0, 0.0, 0.0)
        e180 = CH.attitude_error_vector(half_turn, identity)
        @test norm(e180) ≈ pi atol=1e-12
        @test abs(e180[2]) ≈ pi atol=1e-12
        @test CH.attitude_error_vector(q_actual, -q2) ≈ e atol=1e-12
        @test_throws ArgumentError CH.attitude_error_vector(zero(identity), identity)
        @test_throws ArgumentError descent_attitude_command(zero(up), up)

    end

    @testset "descent datum, selection and constructor contracts" begin
        b, ap = apollo11_descent_targets()
        radius = SM.Earth().Rp_e - 2_000.0
        grid = SM.DEMGrid(Float32[3500 3500; 2500 2500], 20, 40, -5, 5;
                          reference_radius_m=radius)
        terrain = SM.DEMTerrainModel([grid]; reference_radius_m=radius)
        cfg = ApolloDescentConfig(reference_radius_m=radius, site_lat_deg=30.0,
            site_lon_deg=0.0, braking=b, approach=ap)
        state = ApolloDescentState(2)
        guidance = ApolloDescentGuidanceModel(cfg, state, terrain; spacecraft_indices=(2,))
        frame = descent_site_frame(cfg, terrain, radius)
        @test norm(frame.origin_p) ≈ radius + 3_000.0 atol=2e-9
        @test frame.reference_radius_m == radius
        @test guidance.spacecraft_indices == (2,)
        # Use the real guidance environment helper at a nonzero latitude on an
        # oblate planet. A geodetic lookup samples a different terrain height.
        earth = SM.Earth()
        direction = SVector(cosd(30.0), 0.0, sind(30.0))
        pos = (radius + 3_500.0) * direction
        ic = SM.CartesianInitialCondition(pos, cross(earth.ω, pos))
        root = SM.Link(root=true, m=500.0, ref_area=1.0)
        craft = SM.SpacecraftModel(links=[root], root=root, initial_condition=ic, id=97)
        ephem = SM.SimpleEphemeridesModel(prime_meridian_at_reference_rad=0.0)
        base = make_example_config(planet=earth, spacecraft=craft, mission_time=1.0,
            initial_time=SM.InitialTime(year=2000, month=1, day=1, hour=12),
            dynamic_effectors=(), density_model=SM.NoAtmosphereModel(), ephemerides_model=ephem,
            orientation_sim=false, keplerian=true, verbose=false, results=false)
        args = SM.SimConfig._with_configuration(base; dynamics_model=SM.DynamicsModel([deepcopy(craft),craft], ()))
        params = SM.ODEParams(n_sats=2, args=args)
        params.shared_buffers.et_start[] = 0.0
        u = SpaceAGORA.SimulationEngine.build_initial_conditions(args)
        env = GH.descent_environment(guidance, u, params, 0.0, 2)
        @test env.lat_deg ≈ 30.0 atol=1e-12
        @test env.radar_m ≈ 500.0 atol=2e-9
        @test norm(env.v_p) < 1e-10
        @test env.frame.reference_radius_m == radius

        control = ApolloDescentControlModel(ApolloDescentControlConfig(), cfg, state, terrain;
                                            spacecraft_indices=(2,))
        @test CH.touchdown_spec(control, 2).reference_radius_m == radius
        # The first guidance command must not bypass the engine slew limit.
        # Exercise actual control cycles with the real parameter/configuration
        # object; no propagation or synthetic force callback is substituted.
        sample = (q=SVector(0.0, 0.0, 0.0, 1.0), ω=SVector(0.0, 0.0, 0.0))
        CH.update_control!(control, sample, params, 0.0, 2)
        @test control.actuators.last_update_s[2] == 0.0
        @test control.actuators.thrust_n[2] == 0.0
        state.phase_start_s[1, 2] = 0.0
        state.thrust_cmd_n[2] = cfg.max_thrust_n
        CH.update_control!(control, sample, params, 0.05, 2)
        @test control.actuators.thrust_n[2] ≈ control.config.thrust_slew_n_s * 0.05
        previous = control.actuators.thrust_n[2]
        CH.update_control!(control, sample, params, 0.05, 2)
        @test control.actuators.thrust_n[2] == previous
        state.thrust_cmd_n[2] = 0.0
        CH.update_control!(control, sample, params, 0.10, 2)
        @test control.actuators.thrust_n[2] ≈ 0.0 atol=1e-10
        @test control.actuators.thrust_n[1] == 0.0

        @test_throws ArgumentError descent_site_frame(cfg, terrain, radius + 1.0)
        wrong = ApolloDescentConfig(reference_radius_m=radius + 1.0, site_lat_deg=30.0,
            site_lon_deg=0.0, braking=b, approach=ap)
        @test_throws ArgumentError ApolloDescentGuidanceModel(wrong, state, terrain)
        for ids in ((), (0,), (3,), (2,2), (1.5,), (true,))
            @test_throws ArgumentError ApolloDescentGuidanceModel(cfg, state, terrain; spacecraft_indices=ids)
        end
        @test_throws ArgumentError ApolloDescentControlState(0)
        @test_throws ArgumentError ApolloDescentControlModel(ApolloDescentControlConfig(dps_isp_s=0.0), cfg, state, terrain)
        malformed = deepcopy(state)
        pop!(malformed.thrust_cmd_n)
        @test_throws ArgumentError ApolloDescentGuidanceModel(cfg, malformed, terrain)
    end

    @testset "descent thruster diagnostics and actuator shutdown" begin
        # `d` is each jet's exhaust; a Thruster's direction is its thrust, the
        # opposite. Unit lever arms and deliberately non-unit directions give
        # independent +x/+y/+z torque oracles.
        jet(r, d, thrust) = SM.Thruster(location=MVector{3,Float64}(r),
            direction=MVector{3,Float64}((-1.0 .* d)...), max_thrust=thrust)
        engine = jet((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 100.0)
        root = SM.Link(root=true, m=500.0, thrusters=[
            jet((0.0, 1.0, 0.0), (0.0, 0.0, -2.0), 10.0), engine])
        pod = SM.Link(thrusters=[
            jet((0.0, 0.0, 1.0), (-3.0, 0.0, 0.0), 20.0),
            jet((1.0, 0.0, 0.0), (0.0, -4.0, 0.0), 30.0)])
        planet = SM.Moon()
        craft = SM.SpacecraftModel(links=[root, pod], root=root, id=97,
            inertia_tensor=SMatrix{3,3,Float64}(10, 0, 0, 0, 20, 0, 0, 0, 30),
            initial_condition=SM.CartesianInitialCondition(
                SVector(planet.Rp_e + 100.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)))
        layout = CH.descent_thruster_layout(craft)
        @test layout.engine == 2
        @test layout.engine_max_thrust_n == 100.0
        @test layout.jets == [1, 3, 4]
        @test layout.torque_arms_nm ≈ Diagonal([10.0, 20.0, 30.0])
        levels = fill(-1.0, 4)
        @test CH.descent_thruster_levels!(levels, layout, 50.0, SVector(5.0, 10.0, 15.0)) === levels
        @test levels ≈ fill(0.5, 4)
        @test layout.torque_arms_nm * levels[layout.jets] ≈ [5.0, 10.0, 15.0]
        @test_throws BoundsError CH.descent_thruster_levels!(zeros(2), layout,
            50.0, SVector(5.0, 10.0, 15.0))
        # Negative jet demand clips to zero; demands above capacity clip to one.
        CH.descent_thruster_levels!(levels, layout, 150.0, SVector(-10.0, 40.0, 0.0))
        @test levels ≈ [0.0, 1.0, 1.0, 0.0]
        CH.descent_thruster_levels!(levels, layout, NaN, SVector(NaN, 0.0, 0.0))
        @test all(iszero, levels)
        engine_only = CH.descent_thruster_layout((links=[SM.Link(root=true, thrusters=[engine])],))
        @test isempty(engine_only.jets)
        @test CH.descent_thruster_levels!([1.0], engine_only, -1.0, zero(SVector{3,Float64})) == [0.0]

        b, ap = apollo11_descent_targets()
        gcfg = ApolloDescentConfig(reference_radius_m=planet.Rp_e,
            site_lat_deg=0.0, site_lon_deg=0.0, braking=b, approach=ap)
        state = ApolloDescentState(1)
        control = ApolloDescentControlModel(ApolloDescentControlConfig(
            thrust_slew_n_s=30.0, rate_gain=1.0, rate_limit_rad_s=0.1,
            rate_bandwidth=2.0), gcfg, state, NoTerrainModel())
        args = make_example_config(planet=planet, spacecraft=craft, mission_time=1.0,
            initial_time=SM.InitialTime(year=2000, month=1, day=1, hour=12),
            dynamic_effectors=(), density_model=SM.NoAtmosphereModel(),
            ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=true,
            keplerian=false, verbose=false, results=false)
        params = SM.ODEParams(n_sats=1, args=args)
        sample = (q=SVector(0.0, 0.0, 0.0, 1.0), ω=zero(SVector{3,Float64}))
        CH.update_control!(control, sample, params, 0.0, 1)
        cached_layout = control.actuators.thruster_layout[1]
        saved_levels = CH.control_thruster_levels(control, 1)
        @test cached_layout.engine == 2
        @test saved_levels === control.actuators.thruster_level[1]
        @test saved_levels == zeros(4)
        state.phase_start_s[1, 1] = 0.0
        state.thrust_cmd_n[1] = 50.0
        state.attitude_cmd[1] = SVector(1.0, 0.0, 0.0, 0.0)
        CH.update_control!(control, sample, params, 1.0, 1)
        @test control.actuators.thruster_layout[1] === cached_layout
        @test CH.control_thruster_levels(control, 1) === saved_levels
        @test control.actuators.thrust_n[1] == 30.0
        # A half turn saturates the commanded rate at 0.1 rad/s, giving
        # 10 kg m² * 2 /s * 0.1 rad/s = 2 N m, below the RCS torque cap.
        @test control.actuators.torque_nm[1] ≈ SVector(2.0, 0.0, 0.0)
        @test saved_levels ≈ [0.2, 0.3, 0.0, 0.0]
        # The held command, read without the model's storage, and the tick that produced it.
        @test held_command(control, 1) == (thrust_n=30.0, torque_nm=control.actuators.torque_nm[1])
        @test held_command_time(control, 1) == 1.0
        @test held_command(control, 2) === nothing && isnan(held_command_time(control, 2))

        # Contact can precede the first guidance update that creates a site
        # frame. Store the planet-fixed velocity and clear every actuator.
        @test state.site_frame[1] === nothing
        v_p = SVector(1.0, 2.0, -3.0)
        CH.touchdown_spec(control, 1).on_touchdown(2.0,
            SVector(planet.Rp_e, 0.0, 0.0), v_p, 1)
        @test state.phase[1] == :landed
        @test state.touchdown_s[1] == 2.0
        @test state.touchdown_v_mps[1] == v_p
        @test isnan(state.touchdown_miss_m[1])
        @test state.thrust_cmd_n[1] == 0.0
        @test control.actuators.thrust_n[1] == 0.0
        @test iszero(control.actuators.torque_nm[1])
        @test all(iszero, saved_levels)
        # Once the spacecraft is inactive, the engine releases the command at each
        # control tick instead of ticking: engine and jets off.
        @test spacecraft_scope(control) == BoundSpacecraft((1,)) && skip_inactive(control)
        params.is_active[1] = false
        control.actuators.thrust_n[1] = 10.0
        control.actuators.torque_nm[1] = SVector(1.0, 0.0, 0.0)
        saved_levels .= 0.5
        @test SM.SimulationCallbacks._gnc_tick!(CH.update_control!, control, sample, params, 3.0, 1)
        @test CH.control_thruster_levels(control, 1) === saved_levels
        @test all(iszero, saved_levels)
        @test held_command(control, 1) == (thrust_n=0.0, torque_nm=zero(SVector{3, Float64}))
        @test held_command_time(control, 1) == 1.0       # a release is not a tick
    end

    @testset "closed-loop vertical descent on the mounted jets" begin
        # The same descent as the ideal-actuation testset below, flown by the
        # lunar module's engine and sixteen RCS jets in `:mounted` actuation:
        # the engine at the slewed thrust, the jets at their clipped
        # least-norm allocation of the attitude torque, propellant by each
        # jet's specific impulse.
        planet = SM.Moon()
        ephem = SM.SimpleEphemeridesModel()
        initial_time = SM.InitialTime(year=1969, month=7, day=20, hour=20, minute=10, second=0.0)
        lat, lon, h0 = 0.674, 23.473, 60.0
        et = SM.ephemerides_epoch_seconds(ephem, initial_time)
        l_pi = SM.ephemerides_frame_lpi(ephem, planet, et)
        upv = SVector(cosd(lat) * cosd(lon), cosd(lat) * sind(lon), sind(lat))
        r_p = (planet.Rp_e + h0) * upv
        r_i = SVector{3, Float64}(l_pi' * r_p)
        v_i = SVector{3, Float64}(l_pi' * cross(SVector{3, Float64}(planet.ω), r_p))
        up_i = normalize(r_i)
        tilt = SVector(0.0, sin(deg2rad(8.0) / 2), 0.0)
        fwd_i = normalize(SVector{3, Float64}(l_pi' * SVector(sind(lon), -cosd(lon), 0.0)))
        q0 = SM.quat_mult(descent_attitude_command(up_i, up_i, fwd_i), SVector(tilt[1], tilt[2], tilt[3], cos(deg2rad(8.0) / 2)))
        q0 = q0 / norm(q0)
        ic = SM.CartesianInitialCondition(r_i, v_i; q=SVector{4, Float64}(q0))
        thrust(exhaust) = MVector{3, Float64}((-1.0 .* exhaust)...)
        jets = [SM.Thruster(max_thrust=45_040.0, location=MVector{3, Float64}(0.0, 0.0, 1.5), direction=thrust((0.0, 0.0, 1.0)), Isp=311.0)]
        for sx in (-1.0, 1.0), sy in (-1.0, 1.0), exhaust in ((sx, 0.0, 0.0), (0.0, sy, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, -1.0))
            push!(jets, SM.Thruster(max_thrust=445.0, location=MVector{3, Float64}(1.65 * sx, 1.65 * sy, -0.3), direction=thrust(exhaust), Isp=290.0))
        end
        root = SM.Link(root=true, m=6_900.0, dims=MVector{3, Float64}(4.2, 4.2, 4.0), ref_area=1.0, q=MVector{4, Float64}(q0...), thrusters=jets)
        sc = SM.SpacecraftModel(links=[root], root=root, prop_mass=1_000.0, inertia_tensor=SMatrix{3, 3, Float64}(2.2e4, 0, 0, 0, 2.4e4, 0, 0, 0, 2.0e4), initial_condition=ic, id=97)
        ph = DescentPhaseTargets(name=:braking, r_T=SVector(0.0, 0.0, 0.0), v_T=SVector(0.0, 0.0, 0.0), a_T=SVector(0.0, 0.0, 0.0), t_go_initial_s=10.0, t_go_min_s=1.0, altitude_switch_m=Inf)
        gcfg = ApolloDescentConfig(reference_radius_m=planet.Rp_e, site_lat_deg=lat, site_lon_deg=lon, site_height_m=0.0, approach_azimuth_deg=270.0, braking=ph, approach=ph,
            ignition_trim_s=0.0, engine_cutoff_altitude_m=0.0)
        state = ApolloDescentState(1)
        guidance = ApolloDescentGuidanceModel(gcfg, state, NoTerrainModel())
        control = ApolloDescentControlModel(ApolloDescentControlConfig(touchdown_height_m=0.5), gcfg, state, NoTerrainModel(); actuation=:mounted)
        @test SpaceAGORA.actuation_mode(control) === :mounted
        @test_throws ArgumentError ApolloDescentControlModel(ApolloDescentControlConfig(), gcfg, state; actuation=:bogus)
        base = make_example_config(planet=planet, spacecraft=sc, mission_time=240.0, initial_time=initial_time,
            dynamic_effectors=(SM.InverseSquaredGravityModel(),), density_model=SM.NoAtmosphereModel(), ephemerides_model=ephem,
            orientation_sim=true, keplerian=false, EI_km=1.0, verbose=false, results=false)
        args = SM.SimConfig._with_configuration(base;
            guidance_model=SM.GuidanceModel(guidance_effectors=(guidance,), guidance_rates=[1.0]),
            control_model=SM.ControlModel(control_effectors=(control,), control_rates=[0.05]),
            integration_tolerances=SM.IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9,
                dt_max_orbit=0.25, reltol_atmosphere=1e-9, abstol_atmosphere=1e-9, dt_max_atmosphere=0.25),
            solver_config=SM.SolverConfig(solver_mode=:tsit5))
        result = run_simulation(args; isolate_state=false, return_solution=true)
        @test state.phase[1] == :landed
        @test 40.0 < state.touchdown_s[1] < 200.0
        v_td = state.touchdown_v_mps[1]
        @test -1.0 < v_td[3] < -0.1
        @test hypot(v_td[1], v_td[2]) < 0.5
        @test state.touchdown_miss_m[1] < 20.0
        @test control.actuators.attitude_error_rad[1] < deg2rad(3.0)
        @test result.u[end].sc[1].mass < sc.dry_mass + sc.prop_mass
        @test all(iszero, control.actuators.thruster_force_n[1])           # shut down at contact
        @test SpaceAGORA.actuator_levels(control, 1) == zeros(17)

        # The mounted jets, not the configured specific impulses, set the
        # propellant use: with very efficient jets a mounted run uses almost
        # none, while the ideal run of the same first seconds burns the DPS.
        function first_seconds(actuation, isp)
            craft = deepcopy(sc)
            foreach(t -> t.Isp = isp, craft.root.thrusters)
            st = ApolloDescentState(1)
            g = ApolloDescentGuidanceModel(gcfg, st, NoTerrainModel())
            c = ApolloDescentControlModel(ApolloDescentControlConfig(touchdown_height_m=0.5), gcfg, st, NoTerrainModel(); actuation=actuation)
            short = SM.SimConfig._with_configuration(args;
                dynamics_model=SM.DynamicsModel([craft], (SM.InverseSquaredGravityModel(),)),
                guidance_model=SM.GuidanceModel(guidance_effectors=(g,), guidance_rates=[1.0]),
                control_model=SM.ControlModel(control_effectors=(c,), control_rates=[0.05]),
                mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, mission_time=3.0,
                    number_of_orbits=1, keplerian=false, orientation_sim=true, num_steps_to_save=20))
            out = run_simulation(short; return_solution=true)
            return (craft.dry_mass + craft.prop_mass) - out.u[end].sc[1].mass
        end
        used_ideal = first_seconds(:ideal, 1.0e9)
        used_mounted = first_seconds(:mounted, 1.0e9)
        @test used_ideal > 1.0                       # kilograms, from the configured DPS specific impulse
        @test 0.0 <= used_mounted < 1.0e-4

        # The controller commands the thrusters only; every other mounted actuator
        # holds. A gimbal at zero on the descent engine leaves the run bit for bit,
        # and a reaction wheel, whose level vector was once too short, holds still.
        function mounted_seconds(craft; trim=false)
            st = ApolloDescentState(1)
            g = ApolloDescentGuidanceModel(gcfg, st, NoTerrainModel())
            c = ApolloDescentControlModel(ApolloDescentControlConfig(touchdown_height_m=0.5, engine_trim=trim), gcfg, st, NoTerrainModel(); actuation=:mounted)
            short = SM.SimConfig._with_configuration(args;
                dynamics_model=SM.DynamicsModel([craft], (SM.InverseSquaredGravityModel(),)),
                guidance_model=SM.GuidanceModel(guidance_effectors=(g,), guidance_rates=[1.0]),
                control_model=SM.ControlModel(control_effectors=(c,), control_rates=[0.05]),
                mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, mission_time=20.0,
                    number_of_orbits=1, keplerian=false, orientation_sim=true, num_steps_to_save=20))
            out = run_simulation(short; isolate_state=false, return_solution=true)
            return out.u[end].sc[1], c
        end
        x_plain, _ = mounted_seconds(deepcopy(sc))
        gimballed = deepcopy(sc)
        gimbal = add_actuator!(gimballed, gimballed.root, ThrusterGimbal(gimballed.root.thrusters[1]; axis=SVector(1.0, 0.0, 0.0),
            min_angle_rad=-0.1, max_angle_rad=0.1))
        x_gimbal, c_gimbal = mounted_seconds(gimballed)
        @test gimbal.angle_rad == 0.0
        @test length(SpaceAGORA.actuator_levels(c_gimbal, 1)) == 18 && SpaceAGORA.actuator_levels(c_gimbal, 1)[18] == 0.0
        @test reinterpret(UInt64, collect(x_gimbal)) == reinterpret(UInt64, collect(x_plain))
        tilted = deepcopy(sc)
        held = add_actuator!(tilted, tilted.root, ThrusterGimbal(tilted.root.thrusters[1]; axis=SVector(1.0, 0.0, 0.0), angle_rad=0.02))
        mounted_seconds(tilted)
        @test held.angle_rad == 0.02                                       # held at its own angle, not driven to zero

        # Engine trim: the descent engine off its axis, its gimbals turn the
        # thrust line through the centre of mass (the body origin), within the
        # gimbal's rate limit; with trim off they hold.
        @test_throws ArgumentError ApolloDescentControlModel(ApolloDescentControlConfig(engine_trim=true), gcfg,
            ApolloDescentState(1), NoTerrainModel(); actuation=:ideal)
        torque_of(engine) = norm(cross(SVector{3}(engine.location), SpaceAGORA.thrust_direction(engine)))
        offset = deepcopy(sc)
        offset.root.thrusters[1].location .+= SVector(0.05, 0.0, 0.0)
        pitch = add_actuator!(offset, offset.root, ThrusterGimbal(offset.root.thrusters[1]; axis=SVector(0.0, 1.0, 0.0),
            min_angle_rad=-0.1, max_angle_rad=0.1, max_rate_rad_s=0.01))
        nozzle = SVector{3}(offset.root.thrusters[1].location)
        mounted_seconds(offset; trim=true)
        @test pitch.angle_rad ≈ atan(nozzle[1], nozzle[3]) rtol=1e-12              # reached at the rate limit, then held
        @test torque_of(offset.root.thrusters[1]) < 1e-12
        untrimmed = deepcopy(offset)
        untrimmed.root.gimbals[1].angle_rad = 0.0
        mounted_seconds(untrimmed)
        @test untrimmed.root.gimbals[1].angle_rad == 0.0                         # trim off: the gimbal holds
        two = deepcopy(sc)
        two.root.thrusters[1].location .+= SVector(0.05, -0.03, 0.0)
        add_actuator!(two, two.root, ThrusterGimbal(two.root.thrusters[1]; axis=SVector(1.0, 0.0, 0.0)))
        add_actuator!(two, two.root, ThrusterGimbal(two.root.thrusters[1]; axis=SVector(0.0, 1.0, 0.0)))
        mounted_seconds(two; trim=true)
        @test torque_of(two.root.thrusters[1]) < 1e-9                             # two axes, solved in turn
        # One control cycle's solve from zero angles: the alternating solves converge within the cycle.
        once = deepcopy(sc)
        once.root.thrusters[1].location .+= SVector(0.05, -0.03, 0.0)
        add_actuator!(once, once.root, ThrusterGimbal(once.root.thrusters[1]; axis=SVector(1.0, 0.0, 0.0)))
        add_actuator!(once, once.root, ThrusterGimbal(once.root.thrusters[1]; axis=SVector(0.0, 1.0, 0.0)))
        levels = zeros(SM.mounted_actuator_count(once))
        SM.ControlHooks._trim_engine_gimbals!(levels, once, 1)
        SM.apply_gimbal_levels!(once, levels, 0.0)                                # as the engine applies a tick's levels
        @test torque_of(once.root.thrusters[1]) < 1e-12
        wheeled = deepcopy(sc)
        add_actuator!(wheeled, wheeled.root, ReactionWheel(spin_axis=SVector(0.0, 0.0, 1.0)))
        x_wheel, c_wheel = mounted_seconds(wheeled)
        @test length(SpaceAGORA.actuator_levels(c_wheel, 1)) == 18 && x_wheel.h_wheels[1] == 0.0
        @test collect(x_wheel.pos) ≈ collect(x_plain.pos) rtol=1e-12
    end

    @testset "closed-loop vertical descent with RCS attitude control" begin
        planet = SM.Moon()
        ephem = SM.SimpleEphemeridesModel()
        initial_time = SM.InitialTime(year=1969, month=7, day=20, hour=20, minute=10, second=0.0)
        # start 60 m above the site, at rest relative to the ground, engine axis 8 deg off vertical
        lat, lon, h0 = 0.674, 23.473, 60.0
        et = SM.ephemerides_epoch_seconds(ephem, initial_time)
        l_pi = SM.ephemerides_frame_lpi(ephem, planet, et)
        upv = SVector(cosd(lat) * cosd(lon), cosd(lat) * sind(lon), sind(lat))
        r_p = (planet.Rp_e + h0) * upv
        v_p = SVector(0.0, 0.0, 0.0)
        r_i = SVector{3, Float64}(l_pi' * r_p)
        v_i = SVector{3, Float64}(l_pi' * (v_p + cross(SVector{3, Float64}(planet.ω), r_p)))
        up_i = normalize(r_i)
        tilt = SVector(0.0, sin(deg2rad(8.0) / 2), 0.0)
        fwd_i = normalize(SVector{3, Float64}(l_pi' * SVector(sind(lon), -cosd(lon), 0.0)))   # flying west, as the approach azimuth says
        q0 = SM.quat_mult(descent_attitude_command(up_i, up_i, fwd_i), SVector(tilt[1], tilt[2], tilt[3], cos(deg2rad(8.0) / 2)))
        q0 = q0 / norm(q0)
        ic = SM.CartesianInitialCondition(r_i, v_i; q=SVector{4, Float64}(q0))
        root = SM.Link(root=true, m=6_900.0, dims=MVector{3, Float64}(4.2, 4.2, 4.0), ref_area=1.0, q=MVector{4, Float64}(q0...))
        sc = SM.SpacecraftModel(links=[root], root=root, prop_mass=1_000.0, inertia_tensor=SMatrix{3, 3, Float64}(2.2e4, 0, 0, 0, 2.4e4, 0, 0, 0, 2.0e4), initial_condition=ic, id=97)
        # both quadratic phases switch out at once, leaving the rate-of-descent phase
        ph = DescentPhaseTargets(name=:braking, r_T=SVector(0.0, 0.0, 0.0), v_T=SVector(0.0, 0.0, 0.0), a_T=SVector(0.0, 0.0, 0.0), t_go_initial_s=10.0, t_go_min_s=1.0, altitude_switch_m=Inf)
        gcfg = ApolloDescentConfig(reference_radius_m=planet.Rp_e, site_lat_deg=lat, site_lon_deg=lon, site_height_m=0.0, approach_azimuth_deg=270.0, braking=ph, approach=ph,
            ignition_trim_s=0.0, engine_cutoff_altitude_m=0.0)
        # Keep power available until the 0.5 m contact event. An engine cutoff
        # above contact would add a ballistic fall and change the speed bound.
        state = ApolloDescentState(2)
        guidance = ApolloDescentGuidanceModel(gcfg, state, NoTerrainModel(); spacecraft_indices=(2,))
        control = ApolloDescentControlModel(ApolloDescentControlConfig(touchdown_height_m=0.5), gcfg, state, NoTerrainModel(); spacecraft_indices=(2,))
        orbiter = deepcopy(sc)
        orbiter.id = 41
        orbit_r = planet.Rp_e + 100_000.0
        orbit_v = sqrt(planet.μ / orbit_r) * normalize(cross(SVector(0.0, 0.0, 1.0), up_i))
        orbiter.initial_condition = SM.CartesianInitialCondition(orbit_r * up_i, orbit_v;
            q=SVector{4,Float64}(q0))
        base = make_example_config(planet=planet, spacecraft=sc, mission_time=240.0, initial_time=initial_time,
            dynamic_effectors=(SM.InverseSquaredGravityModel(),), density_model=SM.NoAtmosphereModel(), ephemerides_model=ephem,
            orientation_sim=true, keplerian=false, EI_km=1.0, verbose=false, results=false)
        args = SM.SimConfig._with_configuration(base;
            dynamics_model=SM.DynamicsModel([orbiter, sc], (SM.InverseSquaredGravityModel(),)),
            guidance_model=SM.GuidanceModel(guidance_effectors=(guidance,), guidance_rates=[1.0]),
            control_model=SM.ControlModel(control_effectors=(control,), control_rates=[0.05]),
            integration_tolerances=SM.IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9,
                dt_max_orbit=0.25, reltol_atmosphere=1e-9, abstol_atmosphere=1e-9, dt_max_atmosphere=0.25),
            solver_config=SM.SolverConfig(solver_mode=:tsit5))
        # Isolation must preserve the alias shared by the paired models, while
        # leaving both input actuator and guidance state unchanged.
        short = SM.SimConfig._with_configuration(args; mission_configuration=
            SM.MissionConfiguration(mission_type=SM.MissionTime, mission_time=2.1,
                number_of_orbits=1, keplerian=false, orientation_sim=true, num_steps_to_save=20))
        isolated = run_simulation(short; return_solution=true)
        copied_guidance = only(isolated.prob.p.args.guidance_model.guidance_effectors)
        copied_control = only(isolated.prob.p.args.control_model.control_effectors)
        @test copied_guidance.state === copied_control.state
        @test copied_guidance.state !== state
        @test all(isnan, state.phase_start_s)
        @test all(iszero, control.actuators.thrust_n)
        @test isfinite(copied_guidance.state.phase_start_s[1, 2])
        @test isnan(copied_guidance.state.phase_start_s[1, 1])
        @test CH.touchdown_spec(control, 1) === nothing
        @test CH.control_thruster_levels(control, 1) === nothing
        # Guidance on a state of its own: the control would never receive a
        # command, so the run is refused before it starts.
        lone = ApolloDescentGuidanceModel(gcfg, ApolloDescentState(2), NoTerrainModel(); spacecraft_indices=(2,))
        split = SM.SimConfig._with_configuration(short;
            guidance_model=SM.GuidanceModel(guidance_effectors=(lone,), guidance_rates=[1.0]))
        unpaired = try
            run_simulation(split; return_solution=true)
            nothing
        catch e
            e
        end
        @test unpaired isa ArgumentError && occursin(
            r"not paired: spacecraft 2 have an ApolloDescentControlModel but no ApolloDescentGuidanceModel",
            sprint(showerror, unpaired))
        @test all(isnan, control.state.phase_start_s) && all(iszero, control.actuators.thrust_n)

        result = run_simulation(args; isolate_state=false, return_solution=true)
        @test state.phase[2] == :landed
        @test isfinite(state.touchdown_s[2])
        @test 40.0 < state.touchdown_s[2] < 200.0
        v_td = state.touchdown_v_mps[2]
        @test -0.8 < v_td[3] < -0.2
        @test hypot(v_td[1], v_td[2]) < 0.3
        @test state.touchdown_miss_m[2] < 15.0
        @test control.actuators.attitude_error_rad[2] < deg2rad(2.0)
        @test isnan(state.touchdown_s[1])
        @test isnan(state.phase_start_s[1, 1])
        @test result.prob.p.is_active == [true, false]
        @test result.t[end] == 240.0
        @test result.u[end].sc[1].mass ≈ orbiter.dry_mass + orbiter.prop_mass atol=1e-8
        @test result.u[end].sc[2].mass < sc.dry_mass + sc.prop_mass
        @test control.actuators.thrust_n[2] == 0
        @test control.actuators.torque_nm[2] == zero(SVector{3,Float64})
        x_end = SpaceAGORA.SimulationEngine._control_state_sample(result.u[end].sc[2], result.prob.p, 2)
        @test CH.control_wrench(control, x_end, result.prob.p, result.t[end]) ==
            (zero(SVector{3,Float64}), zero(SVector{3,Float64}))
        @test CH.control_propellant_rate(control, x_end, result.prob.p, result.t[end]) == 0.0

    end

    @testset "guidance and control pairing check" begin
        b, ap = apollo11_descent_targets()
        gcfg = ApolloDescentConfig(reference_radius_m=1_737_400.0, site_lat_deg=0.0, site_lon_deg=0.0, braking=b, approach=ap)
        guide(st, ids) = ApolloDescentGuidanceModel(gcfg, st; spacecraft_indices=ids)
        ctrl(st, ids) = ApolloDescentControlModel(ApolloDescentControlConfig(), gcfg, st; spacecraft_indices=ids)
        gnc(gs, cs) = (guidance_model=SM.GuidanceModel(guidance_effectors=Tuple(gs), guidance_rates=fill(1.0, length(gs))),
            control_model=SM.ControlModel(control_effectors=Tuple(cs), control_rates=fill(0.05, length(cs))))
        # run_simulation calls this once per control effector before the run.
        preflight(args) = foreach(c -> SpaceAGORA.SimulationLifecycle.preflight_control(c, args; isolate_state=true),
            args.control_model.control_effectors)
        state, other = ApolloDescentState(2), ApolloDescentState(2)
        unguided = r"not paired: spacecraft 2 have an ApolloDescentControlModel but no ApolloDescentGuidanceModel on the same ApolloDescentState, so their descent engine never receives a command\. Build both"
        uncontrolled = r"not paired: spacecraft 2 have an ApolloDescentGuidanceModel but no ApolloDescentControlModel on the same ApolloDescentState, so no control flies their guidance commands\. Build both"
        # A mismatch is refused at every preflight, not once per session as a
        # warning was, and nothing is logged.
        function check_mismatch(args, message)
            logger = Test.TestLogger(min_level=Logging.Warn)
            with_logger(logger) do
                for _ in 1:2
                    err = try
                        preflight(args)
                        nothing
                    catch e
                        e
                    end
                    @test err isa ArgumentError
                    @test occursin(message, sprint(showerror, err))
                end
            end
            @test isempty(logger.logs)
        end
        # Paired: one pair for both spacecraft, or one bound pair per spacecraft.
        paired = gnc((guide(state, (1, 2)),), (ctrl(state, (1, 2)),))
        @test_logs preflight(paired)
        @test_logs preflight(gnc((guide(state, (1,)), guide(state, (2,))), (ctrl(state, (2,)), ctrl(state, (1,)))))
        # A configuration with only one side is not checked.
        @test_logs preflight(gnc((), (ctrl(state, (1, 2)),)))
        @test_logs preflight(gnc((guide(state, (1, 2)),), ()))
        # Separate states (the unshared state of the design's link check): both
        # directions in one refusal, despite two control models.
        check_mismatch(gnc((guide(other, (1, 2)),), (ctrl(state, (1,)), ctrl(state, (2,)))),
            r"not paired: spacecraft 1, 2 have an ApolloDescentControlModel but no ApolloDescentGuidanceModel on the same ApolloDescentState, so their descent engine never receives a command; spacecraft 1, 2 have an ApolloDescentGuidanceModel but no ApolloDescentControlModel")
        # Same state, different spacecraft selections.
        missing_guidance = gnc((guide(state, (1,)),), (ctrl(state, (1,)), ctrl(state, (2,))))
        missing_control = gnc((guide(state, (1, 2)),), (ctrl(state, (1,)),))
        check_mismatch(missing_guidance, unguided)
        check_mismatch(missing_control, uncontrolled)
        # Each refusal is independent of the ones before it.
        for args in (missing_guidance, missing_control, gnc((guide(other, (1, 2)),), (ctrl(state, (1, 2)),)))
            @test_throws ArgumentError preflight(args)
        end
        # Rule A is broad: another model type exempts its direction even if it
        # has unrelated state. The opposite direction must still be checked.
        for custom_state in (state, other)
            @test_logs preflight(gnc((guide(state, (1,)), CustomDescentGuidance(custom_state)), (ctrl(state, (1, 2)),)))
            @test_logs preflight(gnc((guide(state, (1, 2)),), (ctrl(state, (1,)), CustomDescentControl(custom_state))))
            check_mismatch(gnc((guide(state, (1, 2)), CustomDescentGuidance(custom_state)), (ctrl(state, (1,)),)), uncontrolled)
            check_mismatch(gnc((guide(state, (1,)),), (ctrl(state, (1, 2)), CustomDescentControl(custom_state))), unguided)
        end

        # Under --depwarn=error the pre-run refusals are the same ArgumentErrors:
        # they no longer go through the strict-mode warning. A child process
        # checks the Apollo pairing through its preflight hook, in both
        # directions, and the engine's period, step and scope checks on the RPO
        # pilot.
        strict = raw"""
            using SpaceAGORA, Test, Logging
            const SM = SpaceAGORA.SimulationModel
            b, ap = SM.apollo11_descent_targets()
            cfg = SM.ApolloDescentConfig(reference_radius_m=1_737_400.0,
                site_lat_deg=0.0, site_lon_deg=0.0, braking=b, approach=ap)
            state = SM.ApolloDescentState(2)
            guide(ids) = SM.ApolloDescentGuidanceModel(cfg, state; spacecraft_indices=ids)
            ctrl(ids) = SM.ApolloDescentControlModel(SM.ApolloDescentControlConfig(), cfg, state; spacecraft_indices=ids)
            gnc(gs, cs) = (guidance_model=SM.GuidanceModel(guidance_effectors=(gs,), guidance_rates=[1.0]),
                control_model=SM.ControlModel(control_effectors=(cs,), control_rates=[0.05]))
            preflight(args) = SpaceAGORA.SimulationLifecycle.preflight_control(
                only(args.control_model.control_effectors), args; isolate_state=true)
            logger = Test.TestLogger(min_level=Logging.Warn, respect_maxlog=true)
            with_logger(logger) do
                preflight(gnc(guide((1, 2)), ctrl((1, 2))))
                for (args, expected) in (
                        (gnc(guide((1,)), ctrl((1, 2))), "but no ApolloDescentGuidanceModel"),
                        (gnc(guide((1, 2)), ctrl((1,))), "but no ApolloDescentControlModel"))
                    for _ in 1:2
                        err = try
                            preflight(args)
                            nothing
                        catch e
                            e
                        end
                        @test err isa ArgumentError
                        @test occursin(expected, sprint(showerror, err))
                        @test occursin("spacecraft 2", sprint(showerror, err))
                    end
                end
                checks = SM.SimulationCallbacks
                function rpo(change!)
                    args = make_rpo_configuration(planner=DirectRPOPlanner(), mission_time_s=0.2)
                    change!(only(args.control_model.control_effectors), args)
                    return args
                end
                for (args, check, expected) in (
                        (rpo((c, a) -> (a.control_model.control_rates[1] = 0.2)), checks.check_control_periods,
                         "declares a 0.1 s period"),
                        (rpo((c, a) -> (c.controller.dt_s = 0.2)), checks.check_control_periods,
                         "discretized with a 0.2 s step"),
                        (rpo((c, a) -> (c.target_index = 5)), checks.check_spacecraft_scopes,
                         "names spacecraft 5"))
                    err = try
                        check(args)
                        nothing
                    catch e
                        e
                    end
                    @test err isa ArgumentError
                    @test occursin(expected, sprint(showerror, err))
                end
            end
            @test isempty(logger.logs)
            print("strict_refusals_ok")
            """
        cmd = `$(Base.julia_cmd()) --startup-file=no --compiled-modules=existing --depwarn=error --project=$(pkgdir(SpaceAGORA)) -e $strict`
        @test read(cmd, String) == "strict_refusals_ok"
    end
end

end # module ApolloDescentTests
