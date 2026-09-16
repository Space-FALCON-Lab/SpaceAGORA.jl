using Test
using SpaceAGORA

const PGF = SpaceAGORA.SimulationModel.DynamicEffectors.PlumeGasField

const PLUME_TABLE_PATH = joinpath(dirname(dirname(dirname(@__DIR__))), "data", "psi", "apollo_lmde.json")

# Apollo 11 near the end of the descent: about 7100 kg of lunar weight, and the
# 13.34 kN Morris 2012 uses for the published reference points.
const PGF_APPROACH_THRUST_N = 11_500.0
const PGF_MORRIS_THRUST_N = 13_340.0

"Integral of the surface pressure over the ground plane, by the trapezoidal rule."
function _plume_pressure_integral(field, cfg, thrust_n, height_m; span=8.0, n=4001)
    rs = range(0.0, span * height_m; length=n)
    total = 0.0
    prev = plume_gas_state(field, cfg, thrust_n, height_m, rs[1]).pressure_pa * rs[1]
    for i in 2:n
        cur = plume_gas_state(field, cfg, thrust_n, height_m, rs[i]).pressure_pa * rs[i]
        total += pi * (prev + cur) * (rs[i] - rs[i - 1])
        prev = cur
    end
    return total
end

@testset "PlumeGasField" begin
    analytic = PlumeAnalyticField()
    cfg_a = PlumeSurfaceConfig()
    @test cfg_a.field isa PlumeAnalyticField                 # the analytic field stays the default

    @testset "gas state shape and guards" begin
        st = plume_gas_state(analytic, cfg_a, PGF_APPROACH_THRUST_N, 10.0, 2.0)
        @test st isa PGF.PlumeGasState
        @test propertynames(st) == (:pressure_pa, :shear_pa, :density_kg_m3, :speed_mps, :temperature_k, :mach)
        @test all(isfinite, values(st))
        @test st.pressure_pa > 0.0 && st.shear_pa > 0.0 && st.density_kg_m3 > 0.0
        # the wall jet is at rest under the nozzle and supersonic out in the footprint
        @test plume_gas_state(analytic, cfg_a, PGF_APPROACH_THRUST_N, 10.0, 0.0).speed_mps == 0.0
        @test st.speed_mps > 0.0
        # a shock cannot take the gas past the plume's total temperature
        @test st.temperature_k <= PlumeAnalyticField().chamber_temperature_k
        for bad in ((0.0, 10.0, 2.0), (-1.0, 10.0, 2.0), (NaN, 10.0, 2.0),
                    (PGF_APPROACH_THRUST_N, -1.0, 2.0), (PGF_APPROACH_THRUST_N, 10.0, NaN))
            z = plume_gas_state(analytic, cfg_a, bad...)
            @test z.pressure_pa == 0.0 && z.shear_pa == 0.0 && z.density_kg_m3 == 0.0
        end
    end

    @testset "the analytic field is exactly the Gaussian footprint" begin
        for h in (3.0, 10.0, 40.0), x in (0.0, 0.5, 1.0, 2.0)
            R = max(h * tand(cfg_a.plume_half_angle_deg), cfg_a.nozzle_exit_radius_m)
            p0 = PGF_APPROACH_THRUST_N / (pi * R * R)
            st = plume_gas_state(analytic, cfg_a, PGF_APPROACH_THRUST_N, h, x * R)
            @test st.pressure_pa ≈ p0 * exp(-x^2) rtol = 1e-12
            @test st.shear_pa ≈ cfg_a.friction_coefficient * p0 * 2x * exp(-x^2) rtol = 1e-12
        end
        # its shear coefficient is the configured skin friction, so the dynamic
        # pressure the erosion closure uses is recoverable from the shear stress
        @test PGF.plume_field_shear_coefficient(analytic, cfg_a) == cfg_a.friction_coefficient
    end

    @testset "the surface pressure carries the thrust" begin
        # Classical Newtonian impingement turns all of the plume's axial momentum
        # on the plane, so the integral of the wall pressure is the thrust. The
        # Gaussian footprint is normalized to do the same by construction.
        for h in (5.0, 20.0, 50.0)
            ratio = _plume_pressure_integral(analytic, cfg_a, PGF_APPROACH_THRUST_N, h) / PGF_APPROACH_THRUST_N
            @test ratio ≈ 1.0 rtol = 1e-3
        end
    end

    @testset "shipped table" begin
        @test isfile(PLUME_TABLE_PATH)
        table = load_plume_field(PLUME_TABLE_PATH)
        cfg_t = PlumeSurfaceConfig(field=table)
        @test table isa PlumeFieldTable
        @test table.name == "apollo_lmde"
        @test issorted(table.height_over_diameter) && issorted(table.radius_over_height)
        @test first(table.radius_over_height) == 0.0
        @test table.surface_drag_coefficient == 0.2          # Roberts' rough-surface value

        @testset "momentum: the wall pressure integrates to the thrust" begin
            # The residual is the gap between the source flow's limiting velocity
            # and the engine's vacuum thrust, about one percent for this nozzle.
            for h in (5.0, 20.0, 50.0)
                ratio = _plume_pressure_integral(table, cfg_t, PGF_APPROACH_THRUST_N, h) / PGF_APPROACH_THRUST_N
                @test 0.97 < ratio < 1.03
            end
        end

        @testset "the table and the analytic field agree in magnitude" begin
            for h in (3.0, 5.0, 10.0, 20.0, 50.0)
                pa, Ra = plume_surface_footprint(cfg_a, PGF_APPROACH_THRUST_N, h)
                pt, Rt = plume_surface_footprint(cfg_t, PGF_APPROACH_THRUST_N, h)
                @test 0.5 < pt / pa < 2.0                    # stagnation pressure
                @test 0.5 < Rt / Ra < 2.0                    # footprint radius
                # the shear stress is where the two genuinely differ: Roberts'
                # drag coefficient on the wall jet against a skin-friction
                # coefficient on the static pressure
                qa = plume_quantities(cfg_a, PGF_APPROACH_THRUST_N, h)
                qt = plume_quantities(cfg_t, PGF_APPROACH_THRUST_N, h)
                @test 3.0 < qt.shear_pa / qa.shear_pa < 12.0
            end
        end

        @testset "reference point: Morris 2012 section 4.4.2" begin
            # LMDE hovering 5 m above the surface at 13.34 kN; that work's DSMC
            # peak laminar smooth-wall shear stress is 92 Pa. The shipped table
            # puts its point source at the exit plane, which the same work says
            # under-predicts the near-field stress, and it does: the tolerance
            # here records the gap rather than hiding it.
            peak = maximum(plume_gas_state(table, cfg_t, PGF_MORRIS_THRUST_N, 5.0, r).shear_pa
                           for r in range(0.0, 12.0; length=601))
            @test 20.0 < peak < 92.0
            @test peak / 92.0 > 0.3
        end

        @testset "radial profile accessors" begin
            h = 10.0
            # the wall shear vanishes under the nozzle, peaks off the axis and
            # decays outboard, so the area average is well below the peak
            @test plume_wall_shear(table, cfg_t, PGF_APPROACH_THRUST_N, h, 0.0) == 0.0
            profile = [plume_wall_shear(table, cfg_t, PGF_APPROACH_THRUST_N, h, r) for r in 0.0:0.25:25.0]
            peak, kpeak = findmax(profile)
            @test 1 < kpeak < length(profile)
            mean_in = plume_mean_shear(table, cfg_t, PGF_APPROACH_THRUST_N, h, 1.3 * h)
            @test 0.0 < mean_in < peak
            # the average over a wider disk is smaller, and it scales with thrust
            @test plume_mean_shear(table, cfg_t, PGF_APPROACH_THRUST_N, h, 2.2 * h) < mean_in
            @test plume_mean_shear(table, cfg_t, 2 * PGF_APPROACH_THRUST_N, h, 1.3 * h) ≈ 2 * mean_in rtol = 1e-9
            @test plume_mean_shear(table, cfg_t, PGF_APPROACH_THRUST_N, h, 0.0) == 0.0
            # the scour radius is where the shear last clears a threshold
            a0 = plume_scour_radius(table, cfg_t, PGF_APPROACH_THRUST_N, h, cfg_t.threshold_shear_pa)
            @test a0 > 0.0
            @test plume_wall_shear(table, cfg_t, PGF_APPROACH_THRUST_N, h, 0.95 * a0) > cfg_t.threshold_shear_pa
            @test plume_wall_shear(table, cfg_t, PGF_APPROACH_THRUST_N, h, 1.2 * a0) < cfg_t.threshold_shear_pa
            @test plume_scour_radius(table, cfg_t, PGF_APPROACH_THRUST_N, h, 1.0e9) == 0.0
            @test plume_scour_radius(table, cfg_t, PGF_APPROACH_THRUST_N, 0.0, 0.15) == 0.0
            # The eroding region is much wider than the pressure footprint,
            # because the wall jet's shear stress does not follow the static
            # pressure: the jet keeps a large fraction of the exhaust speed
            # outboard while its density thins. Lane and Metzger, "Estimation of
            # Apollo lunar dust transport using optical extinction measurements",
            # Acta Geophysica 63(2), 2015, put the eroding radius of the Apollo
            # 12 descent at 1.3 to 2.2 times the height.
            for hh in (5.0, 10.0, 15.0)
                scour = plume_scour_radius(table, cfg_t, PGF_APPROACH_THRUST_N, hh, 0.15) / hh
                footprint = plume_surface_footprint(cfg_t, PGF_APPROACH_THRUST_N, hh)[2] / hh
                @test 1.2 < scour < 2.3
                @test footprint < 0.6
                @test scour > 2 * footprint
            end
        end

        @testset "the gas state is not a pressure surrogate" begin
            # Density comes from the plume's own mass flux and temperature from
            # the surface shock, so neither tracks the surface pressure: the
            # temperature falls outboard with the normal velocity component
            # while the density falls far more slowly than the pressure.
            h = 10.0
            near = plume_gas_state(table, cfg_t, PGF_APPROACH_THRUST_N, h, 1.0)
            far = plume_gas_state(table, cfg_t, PGF_APPROACH_THRUST_N, h, 20.0)
            @test near.temperature_k > far.temperature_k
            @test near.pressure_pa / far.pressure_pa > 4 * (near.density_kg_m3 / far.density_kg_m3)
            @test far.mach > near.mach                       # the wall jet accelerates outward
            @test all(s -> s.temperature_k > 0.0 && s.density_kg_m3 > 0.0,
                      (plume_gas_state(table, cfg_t, PGF_APPROACH_THRUST_N, h, r) for r in 0.5:0.5:20.0))
        end

        @testset "thrust scaling" begin
            for h in (4.0, 25.0), r in (0.0, 3.0)
                a = plume_gas_state(table, cfg_t, 10_000.0, h, r)
                b = plume_gas_state(table, cfg_t, 30_000.0, h, r)
                @test b.pressure_pa ≈ 3 * a.pressure_pa rtol = 1e-12
                @test b.shear_pa ≈ 3 * a.shear_pa rtol = 1e-12
                @test b.density_kg_m3 ≈ 3 * a.density_kg_m3 rtol = 1e-12
                @test b.speed_mps ≈ a.speed_mps rtol = 1e-12       # velocity does not scale
                @test b.temperature_k ≈ a.temperature_k rtol = 1e-12
                @test b.mach ≈ a.mach rtol = 1e-12
            end
        end

        @testset "interpolation is continuous and monotone in height" begin
            hs = collect(2.0:0.25:120.0)
            p0s = [plume_surface_footprint(cfg_t, PGF_APPROACH_THRUST_N, h)[1] for h in hs]
            @test all(diff(p0s) .< 0.0)                      # the footprint pressure only falls
            radii = [plume_surface_footprint(cfg_t, PGF_APPROACH_THRUST_N, h)[2] for h in hs]
            @test all(diff(radii) .> 0.0)                    # and the footprint only widens
            # no interpolation seam: a small step in height or radius is a small
            # step in the answer
            for h in (2.5, 7.0, 33.0), r in (0.0, 1.0, 6.0)
                a = plume_gas_state(table, cfg_t, PGF_APPROACH_THRUST_N, h, r).pressure_pa
                b = plume_gas_state(table, cfg_t, PGF_APPROACH_THRUST_N, h + 1e-6, r + 1e-6).pressure_pa
                @test isapprox(a, b; rtol=1e-4)
            end
        end

        @testset "out-of-range queries clamp" begin
            far = last(table.radius_over_height)
            h = 10.0
            @test plume_gas_state(table, cfg_t, PGF_APPROACH_THRUST_N, h, far * h).pressure_pa ==
                  plume_gas_state(table, cfg_t, PGF_APPROACH_THRUST_N, h, 10 * far * h).pressure_pa
            # below the first height node the field freezes instead of diverging
            h_min = first(table.height_over_diameter) * table.exit_diameter_m
            @test plume_gas_state(table, cfg_t, PGF_APPROACH_THRUST_N, 0.0, 0.1) ==
                  plume_gas_state(table, cfg_t, PGF_APPROACH_THRUST_N, h_min, 0.1)
            @test isfinite(plume_surface_footprint(cfg_t, PGF_APPROACH_THRUST_N, 0.0)[1])
            # and above the last height node it clamps to that node's shape
            h_max = last(table.height_over_diameter) * table.exit_diameter_m
            big = plume_gas_state(table, cfg_t, PGF_APPROACH_THRUST_N, 4 * h_max, 0.0)
            @test big.pressure_pa > 0.0 && isfinite(big.pressure_pa)
            @test big.mach ≈ plume_gas_state(table, cfg_t, PGF_APPROACH_THRUST_N, h_max, 0.0).mach rtol = 1e-12
        end

        @testset "the effector runs on a table" begin
            # the ground-effect force does not read the field, so a table changes
            # the diagnostics and not the trajectory
            for h in (1.0, 5.0, 30.0)
                @test plume_ground_effect_force(cfg_t, PGF_APPROACH_THRUST_N, h) ==
                      plume_ground_effect_force(cfg_a, PGF_APPROACH_THRUST_N, h)
            end
            onset_t = plume_erosion_onset_height(cfg_t, PGF_APPROACH_THRUST_N)
            onset_a = plume_erosion_onset_height(cfg_a, PGF_APPROACH_THRUST_N)
            @test onset_t > onset_a                          # the table reports more shear
            @test plume_quantities(cfg_t, PGF_APPROACH_THRUST_N, onset_t * 1.05).erosion_kg_s == 0.0
            @test plume_quantities(cfg_t, PGF_APPROACH_THRUST_N, onset_t * 0.9).erosion_kg_s > 0.0
            @test plume_erosion_onset_height(cfg_t, 0.0) == 0.0
            q = plume_quantities(cfg_t, PGF_APPROACH_THRUST_N, 10.0)
            @test 0.0 < q.inner_m < q.outer_m
            @test q.ejecta_mps >= cfg_t.ejecta_speed_min_mps
        end
    end

    @testset "building a table" begin
        nozzle = PlumeNozzle()
        @test PGF.plume_limit_speed(nozzle) ≈ 3302.6 rtol = 1e-3
        theta_lim, theta_0, decay = PGF.plume_limit_angle(nozzle)
        @test 0.0 < theta_0 < theta_lim < pi                 # the core is inside the limiting streamline
        @test decay > 0.0
        # mass conservation of the angular distribution, the normalization the
        # whole field hangs on
        flux = PGF.plume_angular_mass_flux(nozzle, 45_040.0)
        total = 2 * pi * PGF._simpson(t -> flux(t) * sin(t), 0.0, Float64(pi), 4096)
        @test total ≈ 45_040.0 / (nozzle.specific_impulse_s * 9.80665) rtol = 1e-6
        # the distribution falls away from the axis
        @test flux(0.0) > flux(theta_0) > flux(theta_lim)

        small = build_plume_field_table(nozzle; height_nodes=6, radius_nodes=12)
        @test size(small.pressure_hat) == (6, 12)
        @test all(>=(0.0), small.pressure_hat)
        @test length(small.footprint_over_height) == 6
        @test_throws ArgumentError build_plume_field_table(nozzle; height_nodes=1)

        # moving the point source downstream, as Morris 2012 section 4.6 measured
        # for this nozzle, raises the near-field surface pressure
        shifted = build_plume_field_table(PlumeNozzle(virtual_source_offset_m=1.8); height_nodes=12, radius_nodes=24)
        cfg_s = PlumeSurfaceConfig(field=shifted)
        cfg_b = PlumeSurfaceConfig(field=build_plume_field_table(nozzle; height_nodes=12, radius_nodes=24))
        @test plume_surface_footprint(cfg_s, PGF_MORRIS_THRUST_N, 5.0)[1] >
              plume_surface_footprint(cfg_b, PGF_MORRIS_THRUST_N, 5.0)[1]

        # the published near-field denominator is off by default because it is a
        # density multiplier rather than a geometric change, so it does not
        # conserve the plume's momentum against a surface
        near = build_plume_field_table(PlumeNozzle(near_field_correction=true); height_nodes=12, radius_nodes=32)
        cfg_n = PlumeSurfaceConfig(field=near)
        cfg_p = PlumeSurfaceConfig(field=build_plume_field_table(nozzle; height_nodes=12, radius_nodes=32))
        @test plume_surface_footprint(cfg_n, PGF_MORRIS_THRUST_N, 5.0)[1] >
              1.5 * plume_surface_footprint(cfg_p, PGF_MORRIS_THRUST_N, 5.0)[1]
        @test _plume_pressure_integral(near, cfg_n, PGF_MORRIS_THRUST_N, 5.0) / PGF_MORRIS_THRUST_N > 1.5
        @test _plume_pressure_integral(cfg_p.field, cfg_p, PGF_MORRIS_THRUST_N, 5.0) / PGF_MORRIS_THRUST_N < 1.05

        mktempdir() do dir
            path = save_plume_field(joinpath(dir, "small.json"), small)
            back = load_plume_field(path)
            @test back.name == small.name
            @test back.height_over_diameter ≈ small.height_over_diameter rtol = 1e-5
            @test back.pressure_hat ≈ small.pressure_hat rtol = 1e-5
            @test back.mach ≈ small.mach rtol = 1e-5
            write(joinpath(dir, "junk.json"), "{\"schema\":\"nope\"}")
            @test_throws ArgumentError load_plume_field(joinpath(dir, "junk.json"))
            @test_throws ArgumentError load_plume_field(joinpath(dir, "missing.json"))
        end
    end
end
