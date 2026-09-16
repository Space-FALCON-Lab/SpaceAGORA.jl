using Test
using SpaceAGORA
using StaticArrays
using LinearAlgebra
import SpaceAGORA.TelemetryVerification: make_example_config

const SM = SpaceAGORA.SimulationModel
const GH = SM.GuidanceHooks
const CH = SM.ControlHooks

@testset "ApolloDescent" begin
    @testset "terrain grids" begin
        # 3 x 4 grid, north row first: heights rise 10 m per column eastward and 100 m per row southward
        h = [100.0 110.0 120.0 130.0; 200.0 210.0 220.0 230.0; 300.0 310.0 320.0 330.0]
        g = DEMGrid(h, 0.0, 3.0, 10.0, 14.0; name="test")
        model = DEMTerrainModel([g]; reference_radius_m=1_737_400.0, fallback_height_m=-5.0)
        # cell centers: row 1 at lat 2.5, column 1 at lon 10.5
        @test terrain_height(model, 2.5, 10.5) ≈ 100.0
        @test terrain_height(model, 0.5, 13.5) ≈ 330.0
        @test terrain_height(model, 1.5, 12.0) ≈ 215.0          # halfway between columns 2 and 3 on row 2
        @test terrain_height(model, 2.0, 10.5) ≈ 150.0          # halfway between rows 1 and 2
        @test terrain_height(model, 2.5, 10.5 + 360.0) ≈ 100.0  # longitudes compare modulo 360
        @test terrain_height(model, 2.5, 10.5 - 360.0) ≈ 100.0
        @test terrain_height(model, 5.0, 12.0) == -5.0          # outside: fallback
        @test terrain_radius(model, 2.5, 10.5) ≈ 1_737_500.0
        @test terrain_height(NoTerrainModel(), 1.0, 2.0) == 0.0
        # finest grid answers first
        fine = DEMGrid(fill(7.0f0, 2, 2), 1.0, 2.0, 11.0, 12.0; name="fine")
        @test terrain_height(DEMTerrainModel([fine, g]; reference_radius_m=1.0), 1.5, 11.5) == 7.0
        @test terrain_height(DEMTerrainModel([fine, g]; reference_radius_m=1.0), 2.5, 10.5) == 100.0
        @test_throws ArgumentError DEMGrid(fill(1.0, 1, 3), 0.0, 1.0, 0.0, 1.0)
        @test_throws ArgumentError DEMTerrainModel(DEMGrid[]; reference_radius_m=1.0)
        # round trip through the fetch script's file format
        dir = mktempdir()
        open(joinpath(dir, "dem_test.f32"), "w") do io
            write(io, Float32.(vec(permutedims(h))))   # row-major
        end
        open(joinpath(dir, "dem_test.json"), "w") do io
            write(io, """{"rows": 3, "cols": 4, "lat_min": 0.0, "lat_max": 3.0, "lon_min": 10.0, "lon_max": 14.0, "source": "unit", "reference_radius_m": 1737400.0}""")
        end
        g2 = load_dem_grid(joinpath(dir, "dem_test.json"))
        @test g2.heights == Float32.(h)
        open(joinpath(dir, "site.json"), "w") do io
            write(io, """{"site": {"lat_deg": 1.5, "lon_deg": 12.0, "name": "unit"}, "dem": [{"name": "dem_test", "reference_radius_m": 1737400.0}]}""")
        end
        model2, site = load_site_terrain(joinpath(dir, "site.json"))
        @test site.height_m ≈ 215.0
        @test site.name == "unit"
    end

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
        cfg_thr = ApolloDescentConfig(site_lat_deg=0.0, site_lon_deg=0.0, braking=ph, approach=ph, max_thrust_n=1000.0)
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
    end

    @testset "closed-loop vertical descent with RCS attitude control" begin
        planet = SM.Moon()
        ephem = SM.SimpleEphemeridesModel()
        initial_time = SM.InitialTime(year=1969, month=7, day=20, hour=20, minute=10, second=0.0)
        # start 60 m above the site, at rest relative to the ground, engine axis 8 deg off vertical
        lat, lon, h0 = 0.674, 23.473, 60.0
        et = SM.ephemerides_time_seconds(initial_time, ephem)
        l_pi = SM.planet_frame_lpi(planet, et, ephem)
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
        sc = SM.SpacecraftModel(links=[root], root=root, prop_mass=1_000.0, inertia_tensor=SMatrix{3, 3, Float64}(2.2e4, 0, 0, 0, 2.4e4, 0, 0, 0, 2.0e4), initial_condition=ic, id=1)
        # both quadratic phases switch out at once, leaving the rate-of-descent phase
        ph = DescentPhaseTargets(name=:braking, r_T=SVector(0.0, 0.0, 0.0), v_T=SVector(0.0, 0.0, 0.0), a_T=SVector(0.0, 0.0, 0.0), t_go_initial_s=10.0, t_go_min_s=1.0, altitude_switch_m=Inf)
        gcfg = ApolloDescentConfig(site_lat_deg=lat, site_lon_deg=lon, site_height_m=0.0, approach_azimuth_deg=270.0, braking=ph, approach=ph,
            ignition_trim_s=0.0, engine_cutoff_altitude_m=1.0)
        state = ApolloDescentState(1)
        guidance = ApolloDescentGuidanceModel(gcfg, state, NoTerrainModel())
        control = ApolloDescentControlModel(ApolloDescentControlConfig(touchdown_height_m=0.5), gcfg, state, NoTerrainModel())
        base = make_example_config(planet=planet, spacecraft=sc, mission_time=240.0, initial_time=initial_time,
            dynamic_effectors=(SM.InverseSquaredGravityModel(),), density_model=SM.NoAtmosphereModel(), ephemerides_model=ephem,
            orientation_sim=true, keplerian=false, EI_km=1.0, verbose=false, results=false)
        args = SM.SimulationConfiguration(file_paths=base.file_paths, simulation_settings=base.simulation_settings,
            mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=false, number_of_orbits=1, mission_time=240.0,
                orientation_sim=true, num_steps_to_save=500, data_rate=0.5),
            environment_model=base.environment_model, dynamics_model=base.dynamics_model,
            guidance_model=SM.GuidanceModel(guidance_effectors=(guidance,), guidance_rates=[1.0]),
            navigation_model=base.navigation_model,
            control_model=SM.ControlModel(control_effectors=(control,), control_rates=[0.05]),
            initial_time=initial_time,
            integration_tolerances=SM.IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=0.25, reltol_atmosphere=1e-9, abstol_atmosphere=1e-9, dt_max_atmosphere=0.25))
        # the effectors' shared state must be the live one: no isolating copy of the configuration
        result = run_simulation(args; isolate_state=false)
        @test state.phase[1] == :landed
        @test isfinite(state.touchdown_s[1])
        @test 40.0 < state.touchdown_s[1] < 200.0          # about 60 m at 1 m/s then 0.5 m/s
        v_td = state.touchdown_v_mps[1]
        @test -0.8 < v_td[3] < -0.2                          # descending gently at touchdown
        @test hypot(v_td[1], v_td[2]) < 0.3                  # horizontal drift nulled
        @test state.touchdown_miss_m[1] < 15.0               # landed near the site
        @test control.actuators.attitude_error_rad[1] < deg2rad(2.0)   # the RCS pulled the engine back to vertical
    end
end
