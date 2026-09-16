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
        onset = plume_erosion_onset_height(cfg, APPROACH_THRUST_N)
        @test 20.0 < onset < 40.0                    # Apollo 11: dust from about 30 m
        @test plume_quantities(cfg, APPROACH_THRUST_N, onset * 1.02).erosion_kg_s == 0.0
        @test plume_quantities(cfg, APPROACH_THRUST_N, onset * 1.02).ejecta_mps == 0.0
        @test plume_quantities(cfg, APPROACH_THRUST_N, onset * 0.9).erosion_kg_s > 0.0
        # the onset is exactly where the peak shear crosses the threshold
        @test plume_quantities(cfg, APPROACH_THRUST_N, onset).shear_pa ≈ cfg.threshold_shear_pa rtol = 1e-9
        # a bigger engine reaches further up
        @test plume_erosion_onset_height(cfg, 4 * APPROACH_THRUST_N) ≈ 2 * onset rtol = 1e-9
        @test plume_erosion_onset_height(cfg, 0.0) == 0.0
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
        model = PlumeSurfaceInteractionModel(control, NoTerrainModel(); config=cfg)
        x, env = _plume_samples(planet, cfg.nozzle_offset_m + 5.0)
        rate = plume_quantities(cfg, APPROACH_THRUST_N, 5.0).erosion_kg_s
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
