using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA
const SM = SpaceAGORA.SimulationModel
using .SM
const SE = SpaceAGORA.SimulationEngine
const CB = SM.SimulationCallbacks
const TV = SpaceAGORA.TelemetryVerification
const TH = SM.VehicleThermalModels
const EARTH = make_no_gram_planet(:earth)
const ATMO = ExponentialAtmosphereModel(1.0e-11, 300e3, 50e3;
    temperature_k=800.0, valid_max_altitude_m=1000e3)
const Q_ID = SVector{4, Float64}(0, 0, 0, 1)

# An angle-only thermal law gives an exact integrated reference at fixed attitude.
# Actual Maxwellian heating is checked separately at the same known angles.
struct IncidenceReferenceHeat <: SM.AbstractTypes.AbstractThermalModel end
TH.getHeatRate(::IncidenceReferenceHeat, S::Float64, T::Float64, rho::Float64,
    v::Float64, alpha::Float64) = 1.0 + sin(alpha)
struct LegacyIncidenceEffector <: SM.AbstractTypes.AbstractForceTorqueModel end

pitch_q(deg) = (0.0, sind(deg / 2), 0.0, cosd(deg / 2))
function thermal_spacecraft(; root_pitch=0.0, panel_pitch=20.0)
    ic = InitialCondition(6698e3, 0.0, 45.0, 0.0, 0.0, 0.0,
        Q_ID, SVector{3, Float64}(0, 0, 0))
    return TV.make_three_body_spacecraft(
        bus_dims=(0.5, 0.5, 0.3), panel_dims=(0.02, 0.9, 0.5),
        bus_mass=25.0, panel_mass_each=0.0, panel_offset_y=0.7, ic=ic,
        bus_attitude_q=pitch_q(root_pitch),
        panel_attitude_q_left=pitch_q(panel_pitch),
        panel_attitude_q_right=pitch_q(panel_pitch))
end

function thermal_config(sc, effectors; orientation=false, density=ATMO,
        thermal=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH), duration=600.0)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime, keplerian=false,
            number_of_orbits=1, mission_time=duration, orientation_sim=orientation, num_steps_to_save=200),
        environment_model=EnvironmentModel(planet=EARTH, EI=300.0, density_model=density,
            ephemerides_model=SimpleEphemeridesModel(), thermal_model=thermal, topography=false, wind=false),
        dynamics_model=DynamicsModel([sc], effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(year=2015, month=6, day=1),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-8, abstol_orbit=1e-8,
            dt_max_orbit=5.0, reltol_atmosphere=1e-8, abstol_atmosphere=1e-8, dt_max_atmosphere=1.0))
end

function thermal_fixture(effectors; sc=thermal_spacecraft(), orientation=false,
        density=ATMO, thermal=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH))
    args = thermal_config(sc, effectors; orientation, density, thermal)
    u = SE.build_initial_conditions(args)
    u.sc[1].pos .= SVector(6698e3, 0.0, 0.0)
    u.sc[1].vel .= SVector(1000.0, 7500.0, 500.0)
    p = SM.ODEParams(n_sats=1, args=args)
    p.shared_buffers.et_start[] = SM.ephemerides_time_seconds(args.initial_time, args.environment_model.ephemerides_model)
    p.shared_buffers.current_time[] = 0.0
    SE._initialize_save_cache_buffers!(p)
    SE._initialize_aero_workspace_buffers!(p)
    return (; args, sc, u, p)
end

sample_heat(fx; buffered=false) = copy(CB._compute_stage_heat_rates!(fx.p, fx.u.sc[1], 1, 0.0;
    use_buffered_density=buffered))
function expected_heat(fx, angles)
    pf = SE.sample_planet_frame(fx.u.sc[1], fx.p, 1, 0.0)
    atm = SE.sample_atmosphere(fx.u.sc[1], fx.p, 1, 0.0; write_buffers=false)
    # These fixtures have no wind. Angle expectations never call aero geometry.
    @assert all(iszero, atm.wind_pp)
    v = norm(pf.vel_pp)
    S = v / sqrt(2 * EARTH.R * atm.temperature_k)
    return [TH.getHeatRate(fx.args.environment_model.thermal_model, S, atm.temperature_k,
        atm.rho_kg_m3, v, Float64(a)) for a in angles]
end

@testset "thermal incidence follows current geometry, not a prior wrench" begin
    @testset "known fixed angles with poisoned stored incidence" begin
        for mode in (:max_drag, :attitude, :tumbling_average), scale in (nothing, 0.898)
            fm = AerodynamicCoefficientfM(fixed_attitude_incidence=mode)
            aero = isnothing(scale) ? fm : TV.ScaledAerodynamicCoefficientfM(fm, scale)
            fx = thermal_fixture((aero,); sc=thermal_spacecraft(root_pitch=15.0))
            for link in fx.sc.links
                link.α = -0.37
            end
            angles = mode === :attitude ? deg2rad.([75.0, 55.0, 55.0]) :
                mode === :max_drag ? deg2rad.([90.0, 70.0, 70.0]) : fill(pi / 2, 3)
            expected = expected_heat(fx, angles)
            first = sample_heat(fx) # No dynamic effector has run yet.
            @test first ≈ expected rtol=1e-12 atol=1e-16
            @test all(l -> l.α == -0.37, fx.sc.links)
            state = SE.build_state_sample(fx.u.sc[1], fx.sc, false)
            SE._evaluate_dynamic_effector(aero, fx.u.sc[1], state, fx.p, 1, 0.0)
            @test sample_heat(fx) ≈ expected rtol=1e-12 atol=1e-16
            @test sample_heat(fx) == first
            SE.sample_atmosphere(fx.u.sc[1], fx.p, 1, 0.0; write_buffers=true)
            @test sample_heat(fx; buffered=true) ≈ expected rtol=1e-12 atol=1e-16
        end
        for aero in (AerodynamicCoefficientConstant(), AerodynamicCoefficientNoBallisticFlight())
            fx = thermal_fixture((aero,))
            for link in fx.sc.links
                link.α = NaN
            end
            @test sample_heat(fx) ≈ expected_heat(fx, deg2rad.([90.0, 70.0, 70.0])) rtol=1e-12 atol=1e-16
            @test all(l -> isnan(l.α), fx.sc.links)
        end
    end

    @testset "propagated root attitude and root-relative panels" begin
        for wrapped in (false, true)
            fm = AerodynamicCoefficientfM(fixed_attitude_incidence=:tumbling_average)
            aero = wrapped ? TV.ScaledAerodynamicCoefficientfM(fm, 1.7) : fm
            fx = thermal_fixture((aero,); orientation=true)
            for link in fx.sc.links
                link.α = 0.12
            end
            results = Vector{Float64}[]
            for yaw in (0.0, 45.0)
                fx.u.sc[1].q .= (0.0, 0.0, sind(yaw / 2), cosd(yaw / 2))
                pf = SE.sample_planet_frame(fx.u.sc[1], fx.p, 1, 0.0)
                vi = pf.l_pi' * pf.vel_pp
                # Explicit passive rotations, independent of rot/_aero_link_angles.
                vr = SVector(cosd(yaw)*vi[1] + sind(yaw)*vi[2],
                    -sind(yaw)*vi[1] + cosd(yaw)*vi[2], vi[3])
                vp = SVector(cosd(20.0)*vr[1] - sind(20.0)*vr[3], vr[2],
                    sind(20.0)*vr[1] + cosd(20.0)*vr[3])
                angles = [atan(vr[1], vr[3]), atan(vp[1], vp[3]), atan(vp[1], vp[3])]
                actual = sample_heat(fx)
                @test actual ≈ expected_heat(fx, angles) rtol=1e-12 atol=1e-16
                @test all(l -> l.α == 0.12, fx.sc.links)
                push!(results, actual)
            end
            @test !isapprox(results[1], results[2]; rtol=1e-3)
            @test SVector{4, Float64}(fx.sc.root.q) == Q_ID
        end
    end

    @testset "control updates geometry without a force-order dependency" begin
        for mode in (:max_drag, :attitude), wrapped in (false, true), root_pitch in (0.0, 15.0)
            fm = AerodynamicCoefficientfM(fixed_attitude_incidence=mode)
            aero = wrapped ? TV.ScaledAerodynamicCoefficientfM(fm, 0.898) : fm
            fx = thermal_fixture((aero,); sc=thermal_spacecraft(; root_pitch))
            ctrl = SM.SolarPanelAngleOfAttackControlModel(controlled_panel_links=(2, 3))
            for alpha in (pi / 3, pi / 6)
                SM.ControlHooks._apply_solar_panel_aoa!(ctrl, fx.sc, alpha)
                root_alpha = mode === :attitude ? pi/2 - deg2rad(root_pitch) : pi/2
                panel_alpha = mode === :attitude ? alpha - deg2rad(root_pitch) : alpha
                @test sample_heat(fx) ≈ expected_heat(fx, [root_alpha, panel_alpha, panel_alpha]) rtol=1e-12 atol=1e-16
                @test fx.sc.links[2].α == alpha == fx.sc.links[3].α
            end
        end
    end

    @testset "legacy/custom fallback and vacuum" begin
        for effectors in ((), (LegacyIncidenceEffector(),), (InverseSquaredGravityModel(),))
            fx = thermal_fixture(effectors)
            angles = [0.3, 0.5, 0.7]
            for (link, alpha) in zip(fx.sc.links, angles)
                link.α = alpha
            end
            @test sample_heat(fx) ≈ expected_heat(fx, angles) rtol=1e-12 atol=1e-16
        end
        for density in (NoAtmosphereModel(), ExponentialAtmosphereModel(0.0, 300e3, 50e3; temperature_k=800.0))
            fx = thermal_fixture((AerodynamicCoefficientfM(),); density)
            fx.p.shared_buffers.heat_rates[1] .= 42.0
            @test all(iszero, sample_heat(fx))
        end
    end

    @testset "wind-relative flow in a rotated planet frame" begin
        fx = thermal_fixture((AerodynamicCoefficientfM(),); orientation=true)
        yaw = 45.0
        fx.u.sc[1].q .= (0.0, 0.0, sind(yaw/2), cosd(yaw/2))
        pf = SE.sample_planet_frame(fx.u.sc[1], fx.p, 1, 0.0)
        atm = SE.sample_atmosphere(fx.u.sc[1], fx.p, 1, 0.0; write_buffers=true)
        @test norm(pf.l_pi - SMatrix{3, 3, Float64, 9}(I)) > 0.1
        calm = sample_heat(fx; buffered=true)
        wind = SVector(700.0, 100.0, -200.0) # east, north, up
        fx.p.shared_buffers.winds[1] = wind
        phi, lam = pf.lat_rad, pf.lon_rad
        east = SVector(-sin(lam), cos(lam), 0.0)
        north = SVector(-sin(phi)*cos(lam), -sin(phi)*sin(lam), cos(phi))
        up = SVector(cos(phi)*cos(lam), cos(phi)*sin(lam), sin(phi))
        air_pp = pf.vel_pp - (wind[1]*east + wind[2]*north + wind[3]*up)
        vi = pf.l_pi' * air_pp
        vr = SVector(cosd(yaw)*vi[1] + sind(yaw)*vi[2],
            -sind(yaw)*vi[1] + cosd(yaw)*vi[2], vi[3])
        vp = SVector(cosd(20.0)*vr[1] - sind(20.0)*vr[3], vr[2],
            sind(20.0)*vr[1] + cosd(20.0)*vr[3])
        angles = [atan(vr[1], vr[3]), atan(vp[1], vp[3]), atan(vp[1], vp[3])]
        v = norm(air_pp)
        S = v / sqrt(2 * EARTH.R * atm.temperature_k)
        expected = [TH.getHeatRate(fx.args.environment_model.thermal_model, S,
            atm.temperature_k, atm.rho_kg_m3, v, alpha) for alpha in angles]
        @test sample_heat(fx; buffered=true) ≈ expected rtol=1e-12 atol=1e-16
        @test !isapprox(calm, expected; rtol=1e-3)
    end

    @testset "incidence ownership is explicit for multiple effectors" begin
        fm = AerodynamicCoefficientfM(fixed_attitude_incidence=:attitude)
        wrapped = TV.ScaledAerodynamicCoefficientfM(fm, 0.898)
        for effectors in ((fm, wrapped), (LegacyIncidenceEffector(), wrapped),
                (wrapped, LegacyIncidenceEffector()))
            fx = thermal_fixture(effectors; sc=thermal_spacecraft(root_pitch=15.0))
            @test sample_heat(fx) ≈ expected_heat(fx, deg2rad.([75.0, 55.0, 55.0])) rtol=1e-12 atol=1e-16
        end
        modes = (fm, AerodynamicCoefficientfM(fixed_attitude_incidence=:max_drag))
        @test_throws ArgumentError sample_heat(thermal_fixture(modes))
        # With propagated attitude these policy flags are ignored by aero and heat.
        mixed = thermal_fixture(modes; orientation=true)
        single = thermal_fixture((fm,); orientation=true)
        @test sample_heat(mixed) == sample_heat(single)
    end

    @testset "600-second canted-panel accumulated heat has an independent reference" begin
        results = []
        for wrapped in (false, true)
            sc = thermal_spacecraft(root_pitch=15.0)
            for link in sc.links
                link.α = 0.01
            end
            fm = AerodynamicCoefficientfM(fixed_attitude_incidence=:attitude)
            aero = wrapped ? TV.ScaledAerodynamicCoefficientfM(fm, 1.0) : fm
            args = thermal_config(sc, (InverseSquaredGravityModel(), aero); thermal=IncidenceReferenceHeat())
            sol = run_simulation(args; return_solution=true)
            @test string(sol.retcode) == "Success"
            @test sol.t[end] == 600.0
            expected = 600.0 .* (1.0 .+ sin.(deg2rad.([75.0, 55.0, 55.0])))
            @test collect(sol.u[end].sc[1].heat_loads) ≈ expected rtol=1e-10 atol=1e-8
            push!(results, sol.u[end])
        end
        @test results[1] ≈ results[2] rtol=1e-11 atol=1e-12
    end
end
