using Test
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const RE = SM.DynamicEffectors.RegolithErosion

const MOON_G = 1.625                       # m/s^2

# The threshold shear stress `PlumeSurfaceConfig` currently carries as a fitted
# constant. The point of this module is to derive it instead, so the two derived
# estimates are compared against it here rather than against each other alone.
const FITTED_THRESHOLD_SHEAR_PA = PlumeSurfaceConfig().threshold_shear_pa

# A surface gas state under an Apollo LM, built from the repository's Gaussian
# footprint at the approach thrust plus an ASSUMED exhaust temperature and molar
# mass -- the analytic plume field carries no density or temperature, so these
# two numbers are inputs to the test, not predictions of it. 500 K and
# 21.5 g/mol are the order for storable-hypergolic products expanded into
# vacuum; every assertion below that depends on them is written with a tolerance
# wide enough to survive a factor of two in either.
const ASSUMED_SURFACE_TEMPERATURE_K = 500.0
const ASSUMED_EXHAUST_MOLAR_MASS = 0.0215      # kg/mol
const UNIVERSAL_GAS_CONSTANT = 8.314462618     # J/(mol K)

const APPROACH_THRUST_N = 11_500.0             # Apollo 11 lunar weight near touchdown
const DPS_FULL_THRUST_N = 45_040.0             # LM descent engine at full throttle

"Gas state on the ground at radius `r_m` under a `thrust_n` plume standing `height_m` up."
function apollo_gas_state(thrust_n::Float64, height_m::Float64, r_m::Float64=0.0;
                          cfg::PlumeSurfaceConfig=PlumeSurfaceConfig(),
                          temperature_k::Float64=ASSUMED_SURFACE_TEMPERATURE_K)
    p0, R = plume_surface_footprint(cfg, thrust_n, height_m)
    x = r_m / R
    p = p0 * exp(-x * x)
    tau = cfg.friction_coefficient * p0 * 2.0 * x * exp(-x * x)
    rho = p * ASSUMED_EXHAUST_MOLAR_MASS / (UNIVERSAL_GAS_CONSTANT * temperature_k)
    gas = (pressure_pa=p, shear_pa=tau, density_kg_m3=rho, speed_mps=1000.0,
           temperature_k=temperature_k, mach=3.0)
    return gas, R
end

"The same, at the radius where the Gaussian shear law peaks (r = R / sqrt(2))."
function apollo_peak_shear_state(thrust_n::Float64, height_m::Float64; kwargs...)
    _, R = plume_surface_footprint(PlumeSurfaceConfig(), thrust_n, height_m)
    return apollo_gas_state(thrust_n, height_m, R / sqrt(2.0); kwargs...)
end

const DEAD_GAS = (pressure_pa=0.0, shear_pa=0.0, density_kg_m3=0.0, speed_mps=0.0,
                  temperature_k=0.0, mach=0.0)

@testset "RegolithErosion" begin

    soil = lunar_mare_regolith()

    @testset "lunar mare properties carry the Lunar Sourcebook values" begin
        # Chapter 9 of the Lunar Sourcebook (Carrier, Olhoeft, Mendell, 1991).
        @test soil.bulk_density_kg_m3 == 1_500.0          # Table 9.4, top 15 cm
        @test soil.particle_density_kg_m3 == 3_100.0      # section 9.1.3, specific gravity 3.1
        @test soil.median_diameter_m == 70.0e-6           # section 9.1.2, average of 40-130 um
        @test soil.cohesion_pa == 520.0                   # Table 9.12, 0-15 cm average
        @test soil.friction_angle_deg == 42.0             # Table 9.12, 0-15 cm average
        @test 1.0e-12 <= soil.permeability_m2 <= 7.0e-12  # section 9.1.8, Surveyor 5
        # the porosity must be consistent with the Sourcebook's own relation
        # n = 1 - rho / (G rho_w), which Table 9.5 applies to get 49 percent for
        # the denser top 30 cm
        @test soil.porosity ≈ 1.0 - soil.bulk_density_kg_m3 / soil.particle_density_kg_m3 atol = 0.01
        @test RE.lunar_mare_regolith === lunar_mare_regolith
        @test lunar_mare_regolith(cohesion_pa=1_600.0).cohesion_pa == 1_600.0

        # Metzger (2024a) section 4: <D> = 1.5 D84, D84 = 2.3 D50
        @test mean_lift_height_m(soil) ≈ 1.5 * 2.3 * 70.0e-6
        # Mohr-Coulomb tensile cutoff, 2 c cos(phi) / (1 + sin(phi))
        @test soil_tensile_strength_pa(soil) ≈ 463.0 rtol = 1e-3
        @test soil_tensile_strength_pa(soil) < soil.cohesion_pa
    end

    @testset "the viscous threshold derived from the soil vs. the fitted 0.15 Pa" begin
        tau_t = shields_threshold_shear_pa(soil, MOON_G)
        # Shao and Lu (2000) eq. 22 with the lunar defaults
        @test tau_t ≈ 0.0571 rtol = 1e-2
        # the fitted constant is 2.6 times higher: the headline comparison
        @test 2.0 < FITTED_THRESHOLD_SHEAR_PA / tau_t < 3.5
        # over Shao and Lu's full fitted range for gamma the derived value spans
        # 0.033 to 0.092 Pa, so the fitted constant sits above even the top of it
        lo = shields_threshold_shear_pa(lunar_mare_regolith(cohesion_parameter_kg_s2=1.65e-4), MOON_G)
        hi = shields_threshold_shear_pa(lunar_mare_regolith(cohesion_parameter_kg_s2=5.0e-4), MOON_G)
        @test lo ≈ 0.0333 rtol = 1e-2
        @test hi ≈ 0.0922 rtol = 1e-2
        @test lo < tau_t < hi < FITTED_THRESHOLD_SHEAR_PA

        # On the Moon the threshold is cohesion-dominated: dropping gravity to
        # zero barely changes it, while dropping the cohesion parameter to zero
        # collapses it by more than an order of magnitude.
        @test shields_threshold_shear_pa(soil, 0.0) / tau_t > 0.9
        cohesionless = shields_threshold_shear_pa(lunar_mare_regolith(cohesion_parameter_kg_s2=0.0), MOON_G)
        @test cohesionless / tau_t < 0.1
        # and on Earth, where g is six times larger, the weight term matters more
        @test shields_threshold_shear_pa(soil, 9.80665) > tau_t

        # monotone in the grain size on the cohesive branch
        @test shields_threshold_shear_pa(lunar_mare_regolith(median_diameter_m=40.0e-6), MOON_G) > tau_t
        @test shields_threshold_shear_pa(lunar_mare_regolith(median_diameter_m=130.0e-6), MOON_G) < tau_t
    end

    @testset "the energy-flux threshold is a second, independent estimate" begin
        gas, _ = apollo_gas_state(APPROACH_THRUST_N, 10.0)            # stagnation point
        # v_bar = sqrt(8 p / (pi rho)) must recover sqrt(8 R T / (pi M))
        @test mean_thermal_speed_mps(gas) ≈
              sqrt(8 * UNIVERSAL_GAS_CONSTANT * ASSUMED_SURFACE_TEMPERATURE_K /
                   (pi * ASSUMED_EXHAUST_MOLAR_MASS)) rtol = 1e-9
        tau_e = energy_flux_threshold_shear_pa(gas, soil)
        # Metzger's E_th = 0.123 J/(m^2 s), recast as a wall shear stress at the
        # stagnation state an Apollo plume lays down at 10 m, is 0.159 Pa --
        # six percent from the fitted 0.15 Pa, and reached from Apollo 16 dust
        # opacity rather than Apollo 11 crew reports. The band is a factor of
        # 1.5 because the gas density behind the number is assumed, not sourced.
        @test tau_e ≈ 0.158 rtol = 2e-2
        @test FITTED_THRESHOLD_SHEAR_PA / 1.5 < tau_e < 1.5 * FITTED_THRESHOLD_SHEAR_PA
        # at the radius where the shear peaks the footprint pressure is lower,
        # so the threshold falls with it, but stays within a factor of two
        peak, _ = apollo_peak_shear_state(APPROACH_THRUST_N, 10.0)
        @test FITTED_THRESHOLD_SHEAR_PA / 2 < energy_flux_threshold_shear_pa(peak, soil) < tau_e
        # it is not a soil constant. tau_t = sqrt(E_th rho vbar / 3) with
        # vbar = sqrt(8 p / (pi rho)) makes it proportional to (rho p)^(1/4).
        denser = (gas..., density_kg_m3=16 * gas.density_kg_m3)
        @test energy_flux_threshold_shear_pa(denser, soil) ≈ 2 * tau_e rtol = 1e-9
        harder = (gas..., pressure_pa=16 * gas.pressure_pa)
        @test energy_flux_threshold_shear_pa(harder, soil) ≈ 2 * tau_e rtol = 1e-9
        # and with no gas there is no threshold to cross
        @test !isfinite(energy_flux_threshold_shear_pa(DEAD_GAS, soil))
    end

    @testset "Sutherland viscosity and the pressure diffusion depth" begin
        @test gas_dynamic_viscosity_pa_s(500.0) ≈ 2.67e-5 rtol = 2e-2
        @test gas_dynamic_viscosity_pa_s(1500.0) ≈ 5.26e-5 rtol = 2e-2
        @test gas_dynamic_viscosity_pa_s(0.0) == 0.0
        @test gas_dynamic_viscosity_pa_s(NaN) == 0.0
        @test gas_dynamic_viscosity_pa_s(1500.0) > gas_dynamic_viscosity_pa_s(500.0)

        gas, _ = apollo_gas_state(APPROACH_THRUST_N, 10.0)
        mu = gas_dynamic_viscosity_pa_s(gas.temperature_k)
        d1 = pressure_diffusion_depth_m(gas, soil, 1.0, mu)
        # a centimeter-scale front after a second, per Darcy with the lunar
        # permeability: deep enough to matter, far shallower than the footprint
        @test 1.0e-3 < d1 < 5.0e-2
        # delta ~ sqrt(k p t / (mu n))
        @test pressure_diffusion_depth_m(gas, soil, 4.0, mu) ≈ 2 * d1 rtol = 1e-9
        @test pressure_diffusion_depth_m(gas, lunar_mare_regolith(permeability_m2=4 * soil.permeability_m2),
                                         1.0, mu) ≈ 2 * d1 rtol = 1e-9
        @test pressure_diffusion_depth_m(gas, soil, 0.0, mu) == 0.0
        @test pressure_diffusion_depth_m(DEAD_GAS, soil, 1.0, mu) == 0.0
    end

    @testset "bearing capacity against the Lunar Sourcebook's own anchor" begin
        # Sourcebook section 9.1.9: about 6000 kPa for a 1 m footing, and
        # 3000-11000 kPa for the Apollo 11 LM footpad, computed with Durgunoglu
        # and Mitchell's (1975) wedge-penetration factors on the deeper, denser
        # soil. This implementation uses the classical Prandtl/Reissner/Vesic
        # factors, which are several times smaller, so the test only asks that
        # it land within an order of magnitude and on the conservative side.
        deep = lunar_mare_regolith(bulk_density_kg_m3=1_660.0, cohesion_pa=1_600.0,
                                   friction_angle_deg=49.0)     # Table 9.12, 0-60 cm
        q = soil_bearing_capacity_pa(deep, MOON_G, 1.0)
        @test 3.0e5 < q < 1.1e7                                  # within an order of magnitude of 3-11 MPa
        @test q < 3.0e6                                          # and conservative, as documented

        # structural monotonicity
        @test soil_bearing_capacity_pa(soil, MOON_G, 2.0) > soil_bearing_capacity_pa(soil, MOON_G, 1.0)
        @test soil_bearing_capacity_pa(lunar_mare_regolith(cohesion_pa=2 * soil.cohesion_pa), MOON_G, 1.0) >
              soil_bearing_capacity_pa(soil, MOON_G, 1.0)
        @test soil_bearing_capacity_pa(lunar_mare_regolith(friction_angle_deg=49.0), MOON_G, 1.0) >
              soil_bearing_capacity_pa(soil, MOON_G, 1.0)
        # a soil with no friction angle has no bearing-capacity theory here
        @test !isfinite(soil_bearing_capacity_pa(lunar_mare_regolith(friction_angle_deg=0.0), MOON_G, 1.0))
        # the Moon's sixth of a gravity halves the width-dependent term relative
        # to Earth for a wide footing
        @test soil_bearing_capacity_pa(soil, 9.80665, 5.0) > soil_bearing_capacity_pa(soil, MOON_G, 5.0)
    end

    @testset "every regime is exactly zero below its onset" begin
        env = erosion_environment(footprint_radius_m=2.0)
        regimes = (ViscousErosionRoberts(), ViscousErosionEnergyFlux(),
                   DiffusionDrivenFlow(), BearingCapacityFailure())
        # no gas at all
        for regime in regimes
            @test erosion_rate(regime, DEAD_GAS, soil, MOON_G, env) == 0.0
            @test erosion_onset(regime, DEAD_GAS, soil, MOON_G, env) == false
        end
        # a non-finite gas state must not leak a NaN into the trajectory
        broken = (pressure_pa=NaN, shear_pa=NaN, density_kg_m3=NaN, speed_mps=NaN,
                  temperature_k=NaN, mach=NaN)
        for regime in regimes
            @test erosion_rate(regime, broken, soil, MOON_G, env) == 0.0
            @test erosion_onset(regime, broken, soil, MOON_G, env) == false
        end
        # high above the ground, where the footprint is wide and weak
        far, R_far = apollo_peak_shear_state(APPROACH_THRUST_N, 200.0)
        env_far = erosion_environment(footprint_radius_m=R_far)
        for regime in regimes
            @test erosion_rate(regime, far, soil, MOON_G, env_far) == 0.0
        end
        # the predicate and the rate always agree, at every height
        for h in (0.0, 1.0, 5.0, 20.0, 40.0, 80.0, 200.0), regime in regimes
            gas, R = apollo_peak_shear_state(APPROACH_THRUST_N, h)
            e = erosion_environment(footprint_radius_m=R)
            @test erosion_onset(regime, gas, soil, MOON_G, e) == (erosion_rate(regime, gas, soil, MOON_G, e) > 0.0)
        end
    end

    @testset "viscous erosion crosses exactly at the derived threshold" begin
        tau_t = shields_threshold_shear_pa(soil, MOON_G)
        base, _ = apollo_peak_shear_state(APPROACH_THRUST_N, 10.0)
        env = erosion_environment(footprint_radius_m=4.66)
        just_under = (base..., shear_pa=0.999 * tau_t)
        just_over = (base..., shear_pa=1.001 * tau_t)
        @test erosion_rate(ViscousErosionRoberts(), just_under, soil, MOON_G, env) == 0.0
        @test erosion_rate(ViscousErosionRoberts(), just_over, soil, MOON_G, env) > 0.0
        # and rises monotonically above it
        rates = [erosion_rate(ViscousErosionRoberts(), (base..., shear_pa=s), soil, MOON_G, env)
                 for s in (0.2, 0.5, 1.0, 2.0)]
        @test all(diff(rates) .> 0.0)

        # the energy-flux law crosses at its own threshold, E > E_th
        e_th = soil.erosion_energy_threshold_w_m2
        tau_e = energy_flux_threshold_shear_pa(base, soil)
        @test erosion_rate(ViscousErosionEnergyFlux(), (base..., shear_pa=0.999 * tau_e), soil, MOON_G, env) == 0.0
        @test erosion_rate(ViscousErosionEnergyFlux(), (base..., shear_pa=1.001 * tau_e), soil, MOON_G, env) > 0.0
        # exactly at the threshold the numerator vanishes
        @test RE._lift_height_energy_flux_w_m2((base..., shear_pa=tau_e), soil) ≈ e_th rtol = 1e-9
        # stiffer soil erodes more slowly: cohesive energy density in the
        # denominator of Metzger eq. 16
        stiff = lunar_mare_regolith(cohesive_energy_density_j_m3=10 * soil.cohesive_energy_density_j_m3)
        @test 0.0 < erosion_rate(ViscousErosionEnergyFlux(), base, stiff, MOON_G, env) <
                    erosion_rate(ViscousErosionEnergyFlux(), base, soil, MOON_G, env)
    end

    @testset "the onset heights the regimes predict, and their order" begin
        cfg = PlumeSurfaceConfig()
        function onset_height(regime; thrust_n=APPROACH_THRUST_N)
            for h in 200.0:-0.05:0.0
                gas, R = apollo_peak_shear_state(thrust_n, h)
                erosion_onset(regime, gas, soil, MOON_G, erosion_environment(footprint_radius_m=R)) && return h
            end
            return 0.0
        end

        h_roberts = onset_height(ViscousErosionRoberts())
        h_energy = onset_height(ViscousErosionEnergyFlux())
        h_ddf = onset_height(DiffusionDrivenFlow())
        h_bcf = onset_height(BearingCapacityFailure())

        # Viscous erosion is the regime that reaches highest; the deep-cratering
        # regimes are the ones that need the plume nearly on the ground. This is
        # the ordering Metzger (2024a, section 1) reports for lunar landings.
        @test h_roberts > h_ddf
        @test h_energy > h_ddf
        @test h_ddf >= h_bcf

        # Swapping the fitted 0.15 Pa for the derived 0.057 Pa raises the
        # Roberts onset from the 31 m the repository reports to about 50 m, by
        # h ~ 1 / sqrt(tau_t). Metzger (2024b, section 2.4) puts the observed
        # first dust in the Apollo 16 video at 31.5 m, so the derived threshold
        # overshoots the observation by about 60 percent on this shear law.
        @test h_roberts ≈ plume_erosion_onset_height(cfg, APPROACH_THRUST_N) *
                          sqrt(FITTED_THRESHOLD_SHEAR_PA / shields_threshold_shear_pa(soil, MOON_G)) rtol = 5e-3
        @test 45.0 < h_roberts < 55.0

        # Neither deep-cratering regime fires anywhere in an Apollo descent at
        # the approach throttle, which is what the literature reports.
        @test h_ddf == 0.0
        @test h_bcf == 0.0
        @test !erosion_onset(DiffusionDrivenFlow(), apollo_gas_state(APPROACH_THRUST_N, 0.0)[1], soil, MOON_G,
                             erosion_environment(footprint_radius_m=cfg.nozzle_exit_radius_m))

        # A harder throttle does reach diffusion-driven flow near the ground:
        # the full 45 kN descent engine crosses it from a couple of meters down.
        R_dps = plume_surface_footprint(cfg, DPS_FULL_THRUST_N, 0.0)[2]
        @test erosion_onset(DiffusionDrivenFlow(), apollo_gas_state(DPS_FULL_THRUST_N, 0.0)[1], soil, MOON_G,
                            erosion_environment(footprint_radius_m=R_dps))
        @test onset_height(DiffusionDrivenFlow(); thrust_n=DPS_FULL_THRUST_N) > 0.0
        @test onset_height(DiffusionDrivenFlow(); thrust_n=DPS_FULL_THRUST_N) < 5.0

        # Bearing capacity failure needs roughly ten times the Apollo descent
        # engine through the same nozzle; it is a big-lander regime.
        gas_big, R_big = apollo_gas_state(5.0e5, 0.0)
        @test erosion_onset(BearingCapacityFailure(), gas_big, soil, MOON_G,
                            erosion_environment(footprint_radius_m=R_big, bearing_width_m=2 * R_big))
        gas_mid, R_mid = apollo_gas_state(1.0e5, 0.0)
        @test !erosion_onset(BearingCapacityFailure(), gas_mid, soil, MOON_G,
                             erosion_environment(footprint_radius_m=R_mid, bearing_width_m=2 * R_mid))
    end

    @testset "diffusion-driven flow needs a lateral pressure gradient" begin
        gas, _ = apollo_gas_state(DPS_FULL_THRUST_N, 0.0)
        R = plume_surface_footprint(PlumeSurfaceConfig(), DPS_FULL_THRUST_N, 0.0)[2]
        tight = erosion_environment(footprint_radius_m=R)
        @test erosion_onset(DiffusionDrivenFlow(), gas, soil, MOON_G, tight)
        # spread the same surface pressure over a footprint a hundred times
        # wider and the uplift vanishes: without a lateral gradient the plume's
        # own pressure on the plug cancels the pore pressure under it
        broad = erosion_environment(footprint_radius_m=100 * R)
        @test !erosion_onset(DiffusionDrivenFlow(), gas, soil, MOON_G, broad)
        @test erosion_rate(DiffusionDrivenFlow(), gas, soil, MOON_G, broad) == 0.0
        # an impermeable soil cannot be lifted this way at all
        sealed = lunar_mare_regolith(permeability_m2=1.0e-18)
        @test !erosion_onset(DiffusionDrivenFlow(), gas, sealed, MOON_G, tight)
        # a stronger soil resists longer
        strong = lunar_mare_regolith(cohesion_pa=20 * soil.cohesion_pa)
        @test !erosion_onset(DiffusionDrivenFlow(), gas, strong, MOON_G, tight)
        # once it lifts, the flux is the mobilized layer over the loading time
        mu = gas_dynamic_viscosity_pa_s(gas.temperature_k)
        delta = pressure_diffusion_depth_m(gas, soil, tight.residence_time_s, mu)
        @test erosion_rate(DiffusionDrivenFlow(), gas, soil, MOON_G, tight) ≈
              soil.bulk_density_kg_m3 * delta / tight.residence_time_s rtol = 1e-9
    end

    @testset "bearing capacity failure flux follows the plug momentum balance" begin
        env = erosion_environment(footprint_radius_m=0.75, bearing_width_m=1.5)
        q = soil_bearing_capacity_pa(soil, MOON_G, 1.5)
        just_under = (DEAD_GAS..., pressure_pa=0.999 * q, temperature_k=500.0, density_kg_m3=1.0e-3)
        just_over = (DEAD_GAS..., pressure_pa=1.001 * q, temperature_k=500.0, density_kg_m3=1.0e-3)
        @test erosion_rate(BearingCapacityFailure(), just_under, soil, MOON_G, env) == 0.0
        @test erosion_rate(BearingCapacityFailure(), just_over, soil, MOON_G, env) > 0.0
        excess = 2.0e5
        loaded = (DEAD_GAS..., pressure_pa=q + excess, temperature_k=500.0, density_kg_m3=1.0e-3)
        @test erosion_rate(BearingCapacityFailure(), loaded, soil, MOON_G, env) ≈
              sqrt(2 * soil.bulk_density_kg_m3 * excess) rtol = 1e-9
        # NaN for the bearing width falls back to the footprint diameter
        fallback = erosion_environment(footprint_radius_m=0.75, bearing_width_m=NaN)
        @test erosion_rate(BearingCapacityFailure(), loaded, soil, MOON_G, fallback) ≈
              erosion_rate(BearingCapacityFailure(), loaded, soil, MOON_G, env) rtol = 1e-12
    end

    @testset "the dispatcher totals the rates and names the dominant regime" begin
        gas, R = apollo_peak_shear_state(APPROACH_THRUST_N, 10.0)
        env = erosion_environment(footprint_radius_m=R)
        regimes = default_erosion_regimes()
        @test regimes == (ViscousErosionEnergyFlux(), DiffusionDrivenFlow(), BearingCapacityFailure())
        @test !(ViscousErosionRoberts() in regimes)      # two viscous laws would double count

        out = regolith_erosion_rate(regimes, gas, soil, MOON_G, env)
        @test out.rate_kg_m2_s ≈ sum(erosion_rate(r, gas, soil, MOON_G, env) for r in regimes)
        @test out.dominant === ViscousErosion
        @test out.dominant_rate_kg_m2_s == erosion_rate(ViscousErosionEnergyFlux(), gas, soil, MOON_G, env)
        @test out.active_count == 1
        @test regolith_erosion_rate(gas, soil, MOON_G, env) == out       # default-regimes method

        # nothing active
        quiet = regolith_erosion_rate(regimes, DEAD_GAS, soil, MOON_G, env)
        @test quiet.rate_kg_m2_s == 0.0
        @test quiet.dominant === NoErosion
        @test quiet.dominant_rate_kg_m2_s == 0.0
        @test quiet.active_count == 0

        # A meganewton engine at contact drives all three regimes at once.
        # Which one comes out dominant is a genuine model output rather than a
        # foregone conclusion: on the synthesized gas state used here viscous
        # erosion still wins, because feeding the energy-flux law a density
        # proportional to the pressure inflates it (see the module docstring).
        # The test therefore checks the dispatcher's bookkeeping, not a ranking
        # the model has not earned.
        hard, R_hard = apollo_peak_shear_state(1.0e6, 0.0)
        env_hard = erosion_environment(footprint_radius_m=R_hard, bearing_width_m=2 * R_hard)
        heavy = regolith_erosion_rate(regimes, hard, soil, MOON_G, env_hard)
        @test heavy.active_count == 3
        @test erosion_onset(BearingCapacityFailure(), hard, soil, MOON_G, env_hard)
        @test erosion_onset(DiffusionDrivenFlow(), hard, soil, MOON_G, env_hard)
        @test heavy.dominant_rate_kg_m2_s ==
              maximum(erosion_rate(r, hard, soil, MOON_G, env_hard) for r in regimes)
        @test heavy.rate_kg_m2_s > heavy.dominant_rate_kg_m2_s      # more than one regime contributing

        # regime tags
        @test regime_kind(ViscousErosionRoberts()) === ViscousErosion
        @test regime_kind(ViscousErosionEnergyFlux()) === ViscousErosion
        @test regime_kind(DiffusionDrivenFlow()) === DiffusionDrivenFlowRegime
        @test regime_kind(BearingCapacityFailure()) === BearingCapacityFailureRegime
        @test ViscousErosionRoberts() isa AbstractErosionRegime
    end

    @testset "nothing allocates in the hot path" begin
        gas, R = apollo_peak_shear_state(APPROACH_THRUST_N, 10.0)
        env = erosion_environment(footprint_radius_m=R)
        regimes = default_erosion_regimes()
        # warm up, then measure
        regolith_erosion_rate(regimes, gas, soil, MOON_G, env)
        for regime in (ViscousErosionRoberts(), ViscousErosionEnergyFlux(),
                       DiffusionDrivenFlow(), BearingCapacityFailure())
            erosion_rate(regime, gas, soil, MOON_G, env)
            erosion_onset(regime, gas, soil, MOON_G, env)
            @test @allocated(erosion_rate(regime, gas, soil, MOON_G, env)) == 0
            @test @allocated(erosion_onset(regime, gas, soil, MOON_G, env)) == 0
        end
        @test @allocated(regolith_erosion_rate(regimes, gas, soil, MOON_G, env)) == 0
        @test isbitstype(typeof(regolith_erosion_rate(regimes, gas, soil, MOON_G, env)))
        @test isbitstype(RegolithProperties)
        @test isbitstype(ErosionEnvironment)
    end
end
