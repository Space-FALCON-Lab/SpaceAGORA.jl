using Test
using SpaceAGORA

const ET = SpaceAGORA.SimulationModel.DynamicEffectors.EjectaTransport

# Apollo lunar module descent engine, full throttle and the approach thrust near
# the end of the descent (about 7100 kg of lunar weight), and the Moon.
const EJECTA_FULL_THRUST_N = 45_040.0
const EJECTA_APPROACH_THRUST_N = 11_500.0
const EJECTA_MIN_THRUST_N = 4_500.0
const MOON_GRAVITY_M_S2 = 1.62
const MOON_RADIUS_M = 1.7374e6

# Closed-form solution of dv/dt = k (u - v)^2, v(0) = 0, in a uniform flow of
# speed u, with k = 3 rho C_D / (4 rho_p d) constant. Substituting w = u - v
# gives dw/dt = -k w^2, so w = u/(1 + k u t) and x = u t - ln(1 + k u t)/k.
_uniform_speed(u, k, t) = u * k * u * t / (1.0 + k * u * t)
_uniform_distance(u, k, t) = u * t - log1p(k * u * t) / k

function _time_to_distance(u, k, target)
    lo = 0.0
    hi = 1.0
    while _uniform_distance(u, k, hi) < target
        hi *= 2.0
        hi > 1.0e9 && error("uniform-flow distance never reaches the target")
    end
    for _ in 1:200
        mid = 0.5 * (lo + hi)
        if _uniform_distance(u, k, mid) < target
            lo = mid
        else
            hi = mid
        end
    end
    return 0.5 * (lo + hi)
end

@testset "EjectaTransport" begin
    cfg = EjectaTransportConfig()
    soil = EjectaSoil()
    field = EjectaReferenceGasField()

    @testset "soil and configuration defaults are the sourced lunar values" begin
        # Lunar Sourcebook ch. 9: specific gravity 3.1 (9.1.3), in situ bulk
        # density 1.5 g/cm^3 (9.1.4), average median grain size 70 um (9.1.1).
        @test soil.particle_density_kg_m3 == 3_100.0
        @test soil.bulk_density_kg_m3 == 1_500.0
        @test soil.median_particle_diameter_m == 70.0e-6
        # Immer, Lane, Metzger, Clements (Icarus 214, 2011): 1-3 degrees.
        @test cfg.ejection_angle_min_deg == 1.0
        @test cfg.ejection_angle_max_deg == 3.0
        # Metzger (2023, arXiv:2305.12234): about 3100 m/s for the Apollo LM.
        @test field.exhaust_velocity_mps == 3_100.0
    end

    @testset "the grain is not in continuum flow" begin
        # The regime claim the drag law rests on, checked at the radius where the
        # wall shear stress peaks, over the whole Apollo descent envelope.
        for (h, F) in ((2.0, EJECTA_FULL_THRUST_N), (10.0, EJECTA_FULL_THRUST_N),
                       (30.0, EJECTA_APPROACH_THRUST_N), (40.0, EJECTA_MIN_THRUST_N))
            R = max(h * tand(25.0), 0.75)
            gas = ejecta_gas_state(field, cfg, F, h, R / sqrt(2.0))
            _, _, kn = ET.ejecta_particle_flow_numbers(cfg, gas, gas.speed_mps, 70.0e-6)
            @test kn > 0.01                     # never continuum
            @test ejecta_flow_regime(kn) != ET.EJECTA_REGIME_CONTINUUM
        end
        # At the reference condition of the module header it is transitional.
        gas10 = ejecta_gas_state(field, cfg, EJECTA_FULL_THRUST_N, 10.0, 10.0 * tand(25.0) / sqrt(2.0))
        re, ma, kn = ET.ejecta_particle_flow_numbers(cfg, gas10, gas10.speed_mps, 70.0e-6)
        @test 1.0 < re < 100.0
        @test 1.0 < ma < 6.0
        @test 0.1 < kn <= 10.0
        @test ejecta_flow_regime(kn) == ET.EJECTA_REGIME_TRANSITIONAL
        @test ET.ejecta_regime_name(ejecta_flow_regime(kn)) === :transitional
    end

    @testset "regime boundaries follow Schaaf and Chambre" begin
        @test ejecta_flow_regime(0.005) == ET.EJECTA_REGIME_CONTINUUM
        @test ejecta_flow_regime(0.05) == ET.EJECTA_REGIME_SLIP
        @test ejecta_flow_regime(1.0) == ET.EJECTA_REGIME_TRANSITIONAL
        @test ejecta_flow_regime(50.0) == ET.EJECTA_REGIME_FREE_MOLECULAR
        @test ejecta_flow_regime(Inf) == ET.EJECTA_REGIME_FREE_MOLECULAR
    end

    @testset "Henderson's correlation reaches the free-molecular limit" begin
        # As Re -> 0 at supersonic Mach, Henderson's rarefied branch must become
        # the closed-form drag on a sphere in free molecular flow with diffuse
        # reflection (Schaaf and Chambre 1961).
        for ma in (2.0, 4.0, 6.0)
            thin = (pressure_pa=1.0, shear_pa=0.0, density_kg_m3=1.0e-12,
                    speed_mps=1.0, temperature_k=850.0, mach=ma)
            cd = ejecta_drag_coefficient(cfg, thin, 1.0, 70.0e-6)
            s = ma * sqrt(0.5 * cfg.gas_gamma)
            cd_fm = ET.ejecta_free_molecular_drag_coefficient(s, cfg.grain_temperature_k / 850.0)
            # Henderson's rarefied branch is a fit, not the kinetic-theory
            # result: it tracks the closed form to 5 percent from Ma = 2 and to
            # better than 1 percent by Ma = 6, where the hypersonic limit sets in.
            @test isapprox(cd, cd_fm; rtol=ma >= 6.0 ? 0.01 : 0.05)
        end
        # The hypersonic free-molecular limit of a sphere is 2, which is the
        # constant `plume_surface_interaction.jl` uses today.
        @test isapprox(ET.ejecta_free_molecular_drag_coefficient(1_000.0, 0.0), 2.0; atol=1.0e-5)
        @test ejecta_drag_coefficient(EjectaTransportConfig(drag_model=:constant), nothing, 1.0, 1.0) == 2.0
    end

    @testset "a grain in a uniform flow reaches the closed-form speed" begin
        # A constant drag coefficient makes the launch balance integrable, so the
        # integrator can be checked against the exact solution rather than
        # against itself.
        const_cfg = EjectaTransportConfig(drag_model=:constant, constant_drag_coefficient=2.0)
        u = 1_500.0
        rho = 2.0e-3
        gas = (pressure_pa=1.0, shear_pa=0.0, density_kg_m3=rho, speed_mps=u,
               temperature_k=850.0, mach=2.0)
        for d in (1.0e-5, 7.0e-5, 5.0e-4), len in (0.5, 5.0, 50.0)
            k = 3.0 * rho * const_cfg.constant_drag_coefficient /
                (4.0 * soil.particle_density_kg_m3 * d)
            t = _time_to_distance(u, k, len)
            expected = _uniform_speed(u, k, t)
            got = ejecta_launch_speed(gas, soil, d, len; config=const_cfg)
            @test got.reached_length
            @test isapprox(got.distance_m, len; rtol=1.0e-9)
            @test isapprox(got.speed_mps, expected; rtol=1.0e-6)
        end
        # The terminal speed of the balance is the gas speed itself: a grain can
        # never overtake the flow that is dragging it.
        far = ejecta_launch_speed(gas, soil, 1.0e-6, 1.0e5; config=const_cfg)
        @test far.speed_mps <= u
        @test far.speed_mps > 0.99 * u
    end

    @testset "zero drag gives the textbook ballistic range" begin
        vacuum(_r, _z) = ET.ejecta_vacuum_gas_state()
        for (v0, ang_deg) in ((100.0, 2.0), (350.0, 3.0), (25.0, 45.0))
            ang = deg2rad(ang_deg)
            traj = ejecta_trajectory(vacuum, soil, 70.0e-6, v0, ang, MOON_GRAVITY_M_S2;
                                     config=cfg, escape_speed_mps=Inf)
            @test !traj.escaped
            @test isapprox(traj.range_m, v0^2 * sin(2 * ang) / MOON_GRAVITY_M_S2; rtol=1.0e-9)
            @test isapprox(traj.flight_time_s, 2 * v0 * sin(ang) / MOON_GRAVITY_M_S2; rtol=1.0e-9)
            @test isapprox(traj.apex_m, (v0 * sin(ang))^2 / (2 * MOON_GRAVITY_M_S2); rtol=1.0e-3)
            @test isapprox(traj.impact_speed_mps, v0; rtol=1.0e-9)
        end
        # A grain starting beyond the escape speed never comes back.
        esc = ejecta_escape_speed(MOON_GRAVITY_M_S2, MOON_RADIUS_M)
        @test isapprox(esc, 2_373.0; atol=5.0)          # Immer et al. (2011): about 2373 m/s
        away = ejecta_trajectory(vacuum, soil, 70.0e-6, esc + 1.0, deg2rad(2.0),
                                 MOON_GRAVITY_M_S2; config=cfg, escape_speed_mps=esc)
        @test away.escaped
        @test !isfinite(away.range_m)
    end

    @testset "drag carries a grain further than the ballistic range" begin
        gas_source = ET._EjectaFlightGas(field, cfg, EJECTA_FULL_THRUST_N, 10.0)
        vacuum(_r, _z) = ET.ejecta_vacuum_gas_state()
        v0 = 200.0
        ang = deg2rad(2.0)
        with_drag = ejecta_trajectory(gas_source, soil, 70.0e-6, v0, ang, MOON_GRAVITY_M_S2;
                                      config=cfg, start_radius_m=3.0, escape_speed_mps=Inf)
        without = ejecta_trajectory(vacuum, soil, 70.0e-6, v0, ang, MOON_GRAVITY_M_S2;
                                    config=cfg, start_radius_m=3.0, escape_speed_mps=Inf)
        @test with_drag.range_m > without.range_m
        @test with_drag.impact_speed_mps > v0
    end

    @testset "smaller grains leave faster" begin
        # Metzger (2023): drag force to inertia scales as 1/d, so the fines
        # outrun the coarse fraction. This is the qualitative trend the model has
        # to reproduce; the magnitudes are checked against the Apollo bounds below.
        h = 10.0
        R = h * tand(25.0)
        gas = ejecta_gas_state(field, cfg, EJECTA_FULL_THRUST_N, h, R / sqrt(2.0))
        speeds = [ejecta_launch_speed(gas, soil, d, R; config=cfg).speed_mps
                  for d in (1.0e-6, 5.0e-6, 2.0e-5, 7.0e-5, 2.0e-4, 5.0e-4)]
        @test issorted(speeds; rev=true)
        @test all(s -> 0.0 < s < field.exhaust_velocity_mps, speeds)
    end

    @testset "the distribution's speeds sit inside the Apollo bounds" begin
        # The measured bounds are on the fast tail that pitted Surveyor 3, not on
        # the mean of the population, so the tolerances are deliberately wide.
        # Immer et al. (2011) sec. 1: lower bounds of 40 m/s (Nickle and Carroll
        # 1972), 70 m/s (Jaffe 1972) and 100 m/s (Cour-Palais et al. 1972), and
        # the pit-structure estimate they call most reliable, 300-2000 m/s
        # (Brownlee, Bucher et al. 1972). Metzger (2023): the finest dust
        # approaches the exhaust velocity, about 3100 m/s for the Apollo LM.
        dist = ejecta_distribution(field, cfg, soil, EJECTA_FULL_THRUST_N, 10.0)
        @test dist.sample_count > 0
        @test dist.total_weight > 0.0
        @test 40.0 < dist.mean_speed_mps < 2_000.0
        @test dist.median_speed_mps <= dist.mean_speed_mps * 2.0
        @test 300.0 < dist.max_speed_mps <= field.exhaust_velocity_mps
        # The ejection angle is an input, so all the histogram can do is stay
        # inside the range Immer et al. measured.
        @test cfg.ejection_angle_min_deg <= dist.mean_angle_deg <= cfg.ejection_angle_max_deg
        @test dist.angle_edges_deg[1] >= cfg.ejection_angle_min_deg - 1.0e-9
        @test dist.angle_edges_deg[end] <= cfg.ejection_angle_max_deg + 1.0e-9
        # Histograms are mass fractions of the weighted population.
        for frac in (dist.speed_fraction, dist.angle_fraction, dist.deposition_fraction)
            @test all(>=(0.0), frac)
            @test isapprox(sum(frac), 1.0; atol=1.0e-9)
        end
        @test length(dist.speed_edges_mps) == length(dist.speed_fraction) + 1
        @test dist.mean_deposition_radius_m > 0.0
        @test dist.p90_deposition_radius_m >= dist.mean_deposition_radius_m / 10.0
        @test dist.dominant_regime in (:slip, :transitional, :free_molecular)
        # The mean ejection angle is a first-class output in degrees above the
        # local horizontal at launch, the convention Immer et al. measured in:
        # their per-mission values are 2.6, 2.4, 8.1, 1.4 and 2.0 degrees.
        @test 1.4 <= dist.mean_angle_deg <= 8.1
    end

    @testset "the escape fraction is zero at low thrust and rises with it" begin
        # Nothing reaches 2373 m/s on a high, throttled approach; the fines do
        # once the engine is at full thrust a few meters off the ground.
        quiet = ejecta_distribution(field, cfg, soil, EJECTA_MIN_THRUST_N, 40.0)
        approach = ejecta_distribution(field, cfg, soil, EJECTA_APPROACH_THRUST_N, 30.0)
        mid = ejecta_distribution(field, cfg, soil, EJECTA_FULL_THRUST_N, 10.0)
        low = ejecta_distribution(field, cfg, soil, EJECTA_FULL_THRUST_N, 2.0)
        @test quiet.escape_fraction == 0.0
        # Not identically zero on the approach: a handful of micron grains at the
        # innermost radii already clear 2373 m/s at 11.5 kN and 30 m. Reported as
        # the model gives it (0.4 percent of the weighted mass), not rounded away.
        @test 0.0 < approach.escape_fraction < 0.01
        @test mid.escape_fraction > 0.05
        @test quiet.escape_fraction < approach.escape_fraction < mid.escape_fraction < low.escape_fraction
        @test all(d -> isapprox(d.escape_speed_mps, 2_373.0; atol=5.0),
                  (quiet, approach, mid, low))
        # Speed rises monotonically as the engine comes down on the surface.
        @test quiet.mean_speed_mps < approach.mean_speed_mps < mid.mean_speed_mps < low.mean_speed_mps
    end

    @testset "the per-grain functions do not allocate" begin
        function _hot()
            local_cfg = EjectaTransportConfig()
            local_soil = EjectaSoil()
            local_field = EjectaReferenceGasField()
            gas = ejecta_gas_state(local_field, local_cfg, EJECTA_FULL_THRUST_N, 10.0, 3.3)
            source = ET._EjectaFlightGas(local_field, local_cfg, EJECTA_FULL_THRUST_N, 10.0)
            ejecta_launch_speed(gas, local_soil, 70.0e-6, 3.3; config=local_cfg)
            ejecta_trajectory(source, local_soil, 70.0e-6, 100.0, 0.035, MOON_GRAVITY_M_S2;
                              config=local_cfg, start_radius_m=3.3, escape_speed_mps=2_373.0)
            ejecta_drag_coefficient(local_cfg, gas, 1_000.0, 70.0e-6)
            a_gas = @allocated ejecta_gas_state(local_field, local_cfg, EJECTA_FULL_THRUST_N, 10.0, 3.3)
            a_cd = @allocated ejecta_drag_coefficient(local_cfg, gas, 1_000.0, 70.0e-6)
            a_launch = @allocated ejecta_launch_speed(gas, local_soil, 70.0e-6, 3.3; config=local_cfg)
            a_flight = @allocated ejecta_trajectory(source, local_soil, 70.0e-6, 100.0, 0.035,
                                                    MOON_GRAVITY_M_S2; config=local_cfg,
                                                    start_radius_m=3.3, escape_speed_mps=2_373.0)
            return (a_gas, a_cd, a_launch, a_flight)
        end
        @test all(==(0), _hot())
    end

    @testset "degenerate inputs stay finite" begin
        dead = ET.ejecta_vacuum_gas_state()
        @test ejecta_launch_speed(dead, soil, 70.0e-6, 1.0; config=cfg).speed_mps == 0.0
        @test ejecta_launch_speed(dead, soil, 0.0, 1.0; config=cfg).speed_mps == 0.0
        gas = ejecta_gas_state(field, cfg, EJECTA_FULL_THRUST_N, 10.0, 3.3)
        @test ejecta_launch_speed(gas, soil, 70.0e-6, 0.0; config=cfg).speed_mps == 0.0
        vacuum(_r, _z) = dead
        @test ejecta_trajectory(vacuum, soil, 70.0e-6, 0.0, 0.035, MOON_GRAVITY_M_S2;
                                config=cfg).range_m == 0.0
        @test ejecta_gas_state(field, cfg, 0.0, 10.0, 1.0).density_kg_m3 == 0.0
        @test ejecta_gas_state(field, cfg, EJECTA_FULL_THRUST_N, 10.0, 1.0e6).density_kg_m3 >= 0.0
    end
end
