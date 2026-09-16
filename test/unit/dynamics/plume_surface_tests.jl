using Test
using SpaceAGORA
using StaticArrays
using LinearAlgebra

const SM = SpaceAGORA.SimulationModel
const PSI = SM.DynamicEffectors.PlumeSurfaceInteraction
const ES = SM.EffectorSampling

# The descent control effector the plume model reads: the only thing it needs
# from it is the engine's actual thrust per spacecraft.
_plume_mock_control(thrusts::Vector{Float64}) = (actuators = (thrust_n = thrusts,),)

# Apollo 11 near the end of the descent: about 7100 kg of lunar weight.
const APPROACH_THRUST_N = 11_500.0

function _plume_samples(cfg_planet, height_m::Float64)
    Rp = Float64(cfg_planet.Rp_e)
    up = SVector{3, Float64}(1.0, 0.0, 0.0)
    pos_p = (Rp + height_m) * up
    q = descent_attitude_command(up, up, SVector{3, Float64}(0.0, 1.0, 0.0))
    x = ES.StateSample(pos_p, SVector{3, Float64}(0.0, 0.0, 0.0), 7_100.0; q_ib=q)
    pf = ES.PlanetFrameSample(SMatrix{3, 3, Float64}(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0), pos_p,
        SVector{3, Float64}(0.0, 0.0, 0.0), height_m, 0.0, 0.0)
    env = ES.EnvironmentSample(cfg_planet; planet_frame=pf)
    return x, env
end

@testset "PlumeSurfaceInteraction" begin
    cfg = PlumeSurfaceConfig()

    @testset "surface pressure footprint conserves the thrust" begin
        for h in (5.0, 20.0, 50.0)
            p0, R = plume_surface_footprint(cfg, APPROACH_THRUST_N, h)
            @test R ≈ max(h * tand(cfg.plume_half_angle_deg), cfg.nozzle_exit_radius_m)
            # ∫ p0 exp(-(r/R)^2) dA over the plane is p0 π R^2, the engine thrust
            @test p0 * pi * R^2 ≈ APPROACH_THRUST_N
        end
        # never narrower than the nozzle
        _, R_contact = plume_surface_footprint(cfg, APPROACH_THRUST_N, 0.0)
        @test R_contact == cfg.nozzle_exit_radius_m
    end

    @testset "shear stress and erosion rate at reference heights" begin
        # At 10 m under the Apollo approach thrust the published plume-erosion
        # work puts the surface pressure in the hundreds of pascals, the wall
        # shear in the pascal range, and the erosion rate at kilograms a second.
        q10 = plume_quantities(cfg, APPROACH_THRUST_N, 10.0)
        @test 10.0 < q10.pressure_pa < 1.0e4
        @test 0.1 < q10.shear_pa < 10.0
        @test 1.0 < q10.erosion_kg_s < 100.0
        @test 10.0 <= q10.ejecta_mps <= 200.0        # the order the Apollo films show
        @test 0.0 < q10.inner_m < q10.outer_m        # erosion happens in an annulus, not a disk
        # the shear peaks off the stagnation point, at r = R_p / sqrt(2)
        _, R = plume_surface_footprint(cfg, APPROACH_THRUST_N, 10.0)
        @test q10.inner_m < R / sqrt(2) < q10.outer_m
        # deeper in the descent the footprint tightens and the pressure rises
        q3 = plume_quantities(cfg, APPROACH_THRUST_N, 3.0)
        @test q3.pressure_pa > q10.pressure_pa
        @test q3.shear_pa > q10.shear_pa
        @test q3.outer_m < q10.outer_m
    end

    @testset "erosion onset height" begin
        # The default law derives its threshold from the soil, so the onset is a
        # prediction rather than a fit. It sits above the ~31 m at which the
        # Apollo crews first SAW dust, which is a visibility threshold and not
        # the onset of motion: Lane and Metzger still measure erosion at 36.6 m.
        onset = plume_erosion_onset_height(cfg, APPROACH_THRUST_N)
        @test 35.0 < onset < 80.0
        @test plume_quantities(cfg, APPROACH_THRUST_N, onset * 1.02).erosion_kg_s == 0.0
        @test plume_quantities(cfg, APPROACH_THRUST_N, onset * 1.02).ejecta_mps == 0.0
        @test plume_quantities(cfg, APPROACH_THRUST_N, onset * 0.9).erosion_kg_s > 0.0
        # a bigger engine reaches further up
        @test plume_erosion_onset_height(cfg, 4 * APPROACH_THRUST_N) > onset
        @test plume_erosion_onset_height(cfg, 0.0) == 0.0
    end

    @testset "the fitted Roberts law stays reachable and unchanged" begin
        fitted = PlumeSurfaceConfig(erosion_model=:roberts_fitted)
        onset = plume_erosion_onset_height(fitted, APPROACH_THRUST_N)
        @test 20.0 < onset < 40.0                    # Apollo 11: dust from about 30 m
        # the onset is exactly where the peak shear crosses the fitted threshold
        @test plume_quantities(fitted, APPROACH_THRUST_N, onset).shear_pa ≈ fitted.threshold_shear_pa rtol = 1e-9
        # the closed form scales as the square root of the thrust
        @test plume_erosion_onset_height(fitted, 4 * APPROACH_THRUST_N) ≈ 2 * onset rtol = 1e-9
        @test plume_quantities(fitted, APPROACH_THRUST_N, onset * 1.02).erosion_kg_s == 0.0
        @test plume_quantities(fitted, APPROACH_THRUST_N, onset * 0.9).erosion_kg_s > 0.0
        # the two laws report the same gas and the same geometry, only a
        # different mass flux, so switching laws must not move the footprint
        @test plume_surface_footprint(fitted, APPROACH_THRUST_N, 10.0) ==
              plume_surface_footprint(cfg, APPROACH_THRUST_N, 10.0)
        @test plume_quantities(fitted, APPROACH_THRUST_N, 10.0).shear_pa ==
              plume_quantities(cfg, APPROACH_THRUST_N, 10.0).shear_pa
        @test_throws ArgumentError PlumeSurfaceInteractionModel(_plume_mock_control([1.0]);
                                                                config=PlumeSurfaceConfig(erosion_model=:nonsense))
    end

    @testset "the dominant regime is reported" begin
        # On the Moon under an Apollo-class engine only viscous erosion fires;
        # the deep-cratering regimes need a far larger plume (see
        # regolith_erosion.jl). The dispatcher must say which one it was.
        @test plume_quantities(cfg, APPROACH_THRUST_N, 10.0).regime === ViscousErosion
        @test plume_regime_code(ViscousErosion) == 1.0
        @test plume_quantities(cfg, APPROACH_THRUST_N, 500.0).regime === NoErosion
        @test plume_regime_code(NoErosion) == 0.0
        # gravity is an argument, not a constant: lighter soil weight lowers the
        # threshold, so the same plume erodes at least as much
        weak = plume_quantities(cfg, APPROACH_THRUST_N, 30.0; gravity_m_s2=0.1).erosion_kg_s
        strong = plume_quantities(cfg, APPROACH_THRUST_N, 30.0; gravity_m_s2=9.81).erosion_kg_s
        @test weak >= strong
    end

    @testset "everything is zero with no engine and far above the ground" begin
        for q in (plume_quantities(cfg, 0.0, 5.0),
                  plume_quantities(cfg, APPROACH_THRUST_N, cfg.max_height_m + 1.0),
                  plume_quantities(cfg, APPROACH_THRUST_N, NaN))
            @test q.pressure_pa == 0.0
            @test q.shear_pa == 0.0
            @test q.erosion_kg_s == 0.0
            @test q.ejecta_mps == 0.0
            @test q.ground_effect_n == 0.0
        end
    end

    @testset "ground effect is monotone and vanishes at the cutoff" begin
        D = 2 * cfg.nozzle_exit_radius_m
        F = APPROACH_THRUST_N
        @test plume_ground_effect_force(cfg, F, 0.0) ≈ cfg.ground_effect_max_fraction * F
        heights = collect(0.0:0.05:(cfg.ground_effect_cutoff * D))
        forces = [plume_ground_effect_force(cfg, F, h) for h in heights]
        @test all(diff(forces) .< 0.0)                                   # strictly falling with height
        @test forces[end] ≈ 0.0 atol = 1e-12                             # continuous into the cutoff
        @test plume_ground_effect_force(cfg, F, cfg.ground_effect_cutoff * D) == 0.0
        @test plume_ground_effect_force(cfg, F, 10 * D) == 0.0
        @test plume_ground_effect_force(cfg, 0.0, 0.5) == 0.0
        # a modest fraction of the thrust, never more
        @test maximum(forces) < 0.05 * F
    end

    @testset "the effector's force runs along the engine axis" begin
        planet = SM.Moon()
        control = _plume_mock_control([APPROACH_THRUST_N])
        model = PlumeSurfaceInteractionModel(control, NoTerrainModel(); config=cfg)
        x, env = _plume_samples(planet, cfg.nozzle_offset_m + 1.0)
        force, torque = PSI._plume_wrench(model, x, env, 0.0, 1)
        @test torque == SVector{3, Float64}(0.0, 0.0, 0.0)
        up = SVector{3, Float64}(1.0, 0.0, 0.0)
        @test model.state.height_m[1] ≈ cfg.nozzle_offset_m + 1.0 rtol = 1e-9
        # the plume geometry is measured from the nozzle exit plane, the recorded
        # height from the vehicle reference point the trajectory is integrated at
        @test model.state.ground_effect_n[1] ≈ plume_ground_effect_force(cfg, APPROACH_THRUST_N, 1.0)
        @test force ≈ model.state.ground_effect_n[1] * up                # pushes the vehicle away from the ground
        @test norm(force) > 0.0
        # the engine axis points at the ground
        @test PSI.plume_engine_axis(x) ≈ -up atol = 1e-12
        # high above the surface nothing is left
        x_high, env_high = _plume_samples(planet, cfg.max_height_m + 10.0)
        f_high, _ = PSI._plume_wrench(model, x_high, env_high, 1.0, 1)
        @test f_high == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test model.state.erosion_kg_s[1] == 0.0
    end

    @testset "the cumulative eroded mass survives rejected steps" begin
        planet = SM.Moon()
        control = _plume_mock_control([APPROACH_THRUST_N])
        # The crater feedback is off here so the rate is genuinely constant from
        # one evaluation to the next: this test is about the time-integral
        # guard, and the feedback has its own test above.
        integral_cfg = PlumeSurfaceConfig(crater_height_feedback=false)
        model = PlumeSurfaceInteractionModel(control, NoTerrainModel(); config=integral_cfg)
        x, env = _plume_samples(planet, cfg.nozzle_offset_m + 5.0)
        # The effector weighs the soil against the local mu/r^2 of the planet it
        # is over, not the configuration's nominal gravity, so the expected rate
        # has to be evaluated at the same gravity.
        g_local = Float64(planet.μ) / Float64(planet.Rp_e)^2
        rate = plume_quantities(integral_cfg, APPROACH_THRUST_N, 5.0; gravity_m_s2=g_local).erosion_kg_s
        @test rate > 0.0
        PSI._plume_wrench(model, x, env, 0.0, 1)                # first evaluation only seeds the integral
        @test model.state.eroded_kg[1] == 0.0
        PSI._plume_wrench(model, x, env, 1.0, 1)
        @test model.state.eroded_kg[1] ≈ rate rtol = 1e-9
        after_one = model.state.eroded_kg[1]
        # a rejected step re-evaluates times already passed: nothing may be added
        for t in (0.25, 0.5, 1.0)
            PSI._plume_wrench(model, x, env, t, 1)
            @test model.state.eroded_kg[1] == after_one
        end
        PSI._plume_wrench(model, x, env, 3.0, 1)
        @test model.state.eroded_kg[1] ≈ 3 * rate rtol = 1e-9
        # the integral never runs backwards
        @test model.state.eroded_kg[1] > after_one
    end

    @testset "the crater deepens, feeds back, and keeps its shape" begin
        planet = SM.Moon()
        control = _plume_mock_control([APPROACH_THRUST_N])
        model = PlumeSurfaceInteractionModel(control, NoTerrainModel(); config=cfg)
        radii, depth = plume_crater_profile(model, 1)
        @test length(radii) == cfg.crater_bins
        @test length(depth) == cfg.crater_bins
        @test issorted(radii)                                   # logarithmic, ascending
        @test first(radii) ≈ cfg.crater_min_radius_m
        @test last(radii) ≈ cfg.crater_max_radius_m
        @test all(iszero, depth)                                # nothing eroded yet
        @test model.state.crater_depth_m[1] == 0.0
        @test model.state.crater_radius_m[1] == 0.0

        x, env = _plume_samples(planet, cfg.nozzle_offset_m + 4.0)
        PSI._plume_wrench(model, x, env, 0.0, 1)
        @test all(iszero, depth)                                # the first call only seeds
        PSI._plume_wrench(model, x, env, 30.0, 1)
        @test model.state.crater_depth_m[1] > 0.0
        @test maximum(depth) == model.state.crater_depth_m[1]
        # a crater, not a uniform scrape: the depth falls away from a peak that
        # sits inside the footprint, and vanishes far outside it
        _, peak = findmax(depth)
        @test peak > 1
        @test depth[peak] > depth[end]
        @test depth[end] < 0.01 * depth[peak]
        @test 0.0 < model.state.crater_radius_m[1] <= cfg.crater_max_radius_m
        # it only ever deepens
        before = copy(collect(depth))
        PSI._plume_wrench(model, x, env, 60.0, 1)
        @test all(collect(depth) .>= before)
        # ... and a replayed time adds nothing
        after = copy(collect(depth))
        PSI._plume_wrench(model, x, env, 45.0, 1)
        @test collect(depth) == after

        # the crater under the stagnation point moves the ground away from the
        # nozzle, so the height the field is queried at grows with it
        @test model.state.last_query_height_m[1] ≈ 4.0 + depth[1] rtol = 1e-9
        no_feedback = PlumeSurfaceConfig(crater_height_feedback=false)
        m2 = PlumeSurfaceInteractionModel(control, NoTerrainModel(); config=no_feedback)
        PSI._plume_wrench(m2, x, env, 0.0, 1)
        PSI._plume_wrench(m2, x, env, 30.0, 1)
        @test m2.state.last_query_height_m[1] ≈ 4.0 rtol = 1e-12
        @test m2.state.crater_depth_m[1] > 0.0
    end

    @testset "the ejecta diagnostic is off the right-hand side" begin
        planet = SM.Moon()
        control = _plume_mock_control([APPROACH_THRUST_N])
        model = PlumeSurfaceInteractionModel(control, NoTerrainModel(); config=cfg)
        x, env = _plume_samples(planet, cfg.nozzle_offset_m + 6.0)
        PSI._plume_wrench(model, x, env, 0.0, 1)
        # the right-hand side never populates the ejecta summary
        @test model.state.ejecta_angle_deg[1] == 0.0
        @test model.state.ejecta_range_m[1] == 0.0

        plume_refresh_ejecta!(model, 0.0, 1)
        angle = model.state.ejecta_angle_deg[1]
        # the angle is an input spread over the Apollo film range, so it must
        # land inside it (Immer et al. 2008: 1 to 3 degrees above the horizontal)
        @test cfg.ejecta.ejection_angle_min_deg <= angle <= cfg.ejecta.ejection_angle_max_deg
        @test model.state.ejecta_range_m[1] > 0.0
        @test 0.0 <= model.state.ejecta_escape_frac[1] <= 1.0
        # a second call at the same time recomputes nothing
        model.state.ejecta_angle_deg[1] = -1.0
        plume_refresh_ejecta!(model, 0.0, 1)
        @test model.state.ejecta_angle_deg[1] == -1.0
        plume_refresh_ejecta!(model, 1.0, 1)
        @test model.state.ejecta_angle_deg[1] ≈ angle rtol = 1e-9

        summary = plume_ejecta_summary(model, 1)
        @test summary !== nothing
        @test summary.mean_angle_deg ≈ angle rtol = 1e-12
        # the grain mass is weighted by a lognormal fitted to the soil's own D50
        # and D84/D50, so the median size carries most of it and the micron
        # fines, which drive the deposition radius, carry almost none
        w = ejecta_lognormal_mass_weights(cfg.ejecta_sizes, cfg.soil)
        @test sum(w) ≈ 1.0
        @test argmax(w) == argmin(abs.(cfg.ejecta_sizes .- cfg.soil.median_diameter_m))
        @test w[1] < 0.05

        # with the diagnostic off nothing is computed at all
        quiet = PlumeSurfaceInteractionModel(control, NoTerrainModel();
                                             config=PlumeSurfaceConfig(ejecta_diagnostic=false))
        PSI._plume_wrench(quiet, x, env, 0.0, 1)
        plume_refresh_ejecta!(quiet, 0.0, 1)
        @test quiet.state.ejecta_angle_deg[1] == 0.0
    end

    @testset "state sizing and construction" begin
        control = _plume_mock_control([1.0, 2.0, 3.0])
        model = PlumeSurfaceInteractionModel(control)
        @test length(model.state.height_m) == 3
        @test model.terrain isa NoTerrainModel
        @test PlumeSurfaceState(2).eroded_kg == zeros(2)
        @test_throws ArgumentError PlumeSurfaceState(0)
        @test_throws ArgumentError PlumeSurfaceInteractionModel((a = 1,))
    end
end
