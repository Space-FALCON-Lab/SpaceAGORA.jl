module MeshAeroTests
using Test
using LinearAlgebra
using StaticArrays
using SpecialFunctions
using JSON
using Random
using DataFrames
using Arrow
using SpaceAGORA
using SpaceAGORA.SimulationModel

import SpaceAGORA.TelemetryVerification: make_example_config

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const AE = SM.DynamicEffectors.AerodynamicEffectors
const ZERO3 = SVector{3, Float64}(0.0, 0.0, 0.0)
const E1 = SVector{3, Float64}(1.0, 0.0, 0.0)
const I3 = SMatrix{3, 3, Float64, 9}(I)

# Twelve counter-clockwise (outward) triangles of an axis-aligned box centerd on the origin.
function _box_triangles(lx, ly, lz)
    hx, hy, hz = lx / 2, ly / 2, lz / 2
    quads = (
        ((hx, -hy, -hz), (hx, hy, -hz), (hx, hy, hz), (hx, -hy, hz)),        # +x
        ((-hx, -hy, -hz), (-hx, -hy, hz), (-hx, hy, hz), (-hx, hy, -hz)),    # -x
        ((-hx, hy, -hz), (-hx, hy, hz), (hx, hy, hz), (hx, hy, -hz)),        # +y
        ((-hx, -hy, -hz), (hx, -hy, -hz), (hx, -hy, hz), (-hx, -hy, hz)),    # -y
        ((-hx, -hy, hz), (hx, -hy, hz), (hx, hy, hz), (-hx, hy, hz)),        # +z
        ((-hx, -hy, -hz), (-hx, hy, -hz), (hx, hy, -hz), (hx, -hy, -hz)),    # -z
    )
    tris = Float64[]
    for (a, b, c, d) in quads
        append!(tris, a..., b..., c..., a..., c..., d...)
    end
    return reshape(tris, 3, :)
end

function _uv_sphere(R; nlat=48, nlon=96)
    tris = Float64[]
    p(th, ph) = R .* (sin(th) * cos(ph), sin(th) * sin(ph), cos(th))
    for i in 0:(nlat - 1), j in 0:(nlon - 1)
        th0 = pi * i / nlat; th1 = pi * (i + 1) / nlat; ph0 = 2pi * j / nlon; ph1 = 2pi * (j + 1) / nlon
        a, b, c, d = p(th0, ph0), p(th1, ph0), p(th1, ph1), p(th0, ph1)
        i > 0 && append!(tris, a..., b..., d...)
        i < nlat - 1 && append!(tris, b..., c..., d...)
    end
    return reshape(tris, 3, :)
end

# One-sided square plate of side `w` with outward normal +x at x = x0.
function _plate_triangles(w, x0; y0=0.0, z0=0.0)
    h = w / 2
    a = (x0, y0 - h, z0 - h); b = (x0, y0 + h, z0 - h); c = (x0, y0 + h, z0 + h); d = (x0, y0 - h, z0 + h)
    return reshape(Float64[a..., b..., c..., a..., c..., d...], 3, :)
end

# Sphere drag coefficient in free-molecular flow, fully diffuse, Tw = T (Schaaf & Chambre).
_sphere_cd(s) = (2s^2 + 1) / (sqrt(pi) * s^3) * exp(-s^2) + (4s^4 + 4s^2 - 1) / (2s^4) * erf(s) + 2sqrt(pi) / (3s)

# Flow direction for the fM model's (alpha, beta) convention: body_frame_velocity
# = (cos(beta) sin(alpha), sin(beta), cos(beta) cos(alpha)), alpha = pi/2 is flow along +x.
_flow_from_angles(alpha, beta) = SVector{3, Float64}(cos(beta) * sin(alpha), sin(beta), cos(beta) * cos(alpha))

# Scalar-last quaternion (x, y, z, w) of a rotation by `angle` about `axis`.
function _quat_axis_angle(axis, angle)
    a = normalize(SVector{3, Float64}(axis))
    return SVector{4, Float64}(a[1] * sin(angle / 2), a[2] * sin(angle / 2), a[3] * sin(angle / 2), cos(angle / 2))
end

# Quaternion q with SM.rot(q) == R for a rotation matrix R whose rows are the
# body axes in the reference frame (the engine's inertial-to-body convention).
function _quat_from_rows(R::SMatrix{3, 3, Float64})
    w = sqrt(max(0.0, 1.0 + R[1, 1] + R[2, 2] + R[3, 3])) / 2
    q = SVector{4, Float64}((R[3, 2] - R[2, 3]) / (4w), (R[1, 3] - R[3, 1]) / (4w), (R[2, 1] - R[1, 2]) / (4w), w)
    isapprox(SM.rot(q), R; atol=1e-10) && return q
    qc = SVector{4, Float64}(-q[1], -q[2], -q[3], q[4])
    isapprox(SM.rot(qc), R; atol=1e-10) || error("quaternion convention mismatch in the test helper")
    return qc
end

# Direct wrench evaluation with a hand-built sample; the planet-fixed frame is
# taken as the inertial frame (l_pi = I) and the wind is zero.
function _samples(planet, sc, q_ib, pos_pp, vel_pp, rho, T)
    r = norm(pos_pp)
    pf = SM.PlanetFrameSample(I3, pos_pp, vel_pp, r - planet.Rp_e, asin(pos_pp[3] / r), atan(pos_pp[2], pos_pp[1]))
    at = SM.AtmosphereSample(rho, T, ZERO3)
    x = SM.StateSample(pos_pp, vel_pp, 400.0; q_ib=q_ib, ω_body=(q_ib === nothing ? nothing : ZERO3), spacecraft=sc)
    return x, SM.EnvironmentSample(planet; planet_frame=pf, atmosphere=at)
end

# ----------------------------------------------------------------------------
# Engine runs
# ----------------------------------------------------------------------------
function _record(args, fields)
    recorder = TrajectoryRecorder(args; save_fields=fields)
    metadata = withenv(
        "SPACEAGORA_RHS_EXECUTION_MODE" => "serial",
        "SPACEAGORA_RHS_CALIBRATE" => "off", "SPACEAGORA_RHS_IDENTIFY" => "0",
        "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => "0",
        "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => "0",
    ) do
        run_simulation(args; return_solver_metadata=true, visualization=false,
            extra_callbacks=(get_trajectory_recorder_callback(recorder),))
    end
    @test metadata.retcode == "Success"
    return collect(trajectory_times(recorder)), trajectory_save_data(recorder)
end
_cache_field(name, getter) = SaveField(name, (u, t, integrator) -> begin
    cache = getter(integrator.p.save_cache)
    [i <= length(cache) ? SVector{3, Float64}(cache[i]) : ZERO3 for i in eachindex(u.sc)]
end; per_satellite=true)
const STATE_FIELDS = SaveField[
    SaveField(:pos, (u, t, integrator) -> [SVector{3, Float64}(sc.pos) for sc in u.sc]; per_satellite=true),
    SaveField(:vel, (u, t, integrator) -> [SVector{3, Float64}(sc.vel) for sc in u.sc]; per_satellite=true),
    _cache_field(:drag, c -> c.drag_cache), _cache_field(:lift, c -> c.lift_cache), _cache_field(:cross, c -> c.cross_cache),
]
_trapz(t, y) = sum(0.5 * (y[k] + y[k + 1]) * (t[k + 1] - t[k]) for k in 1:(length(t) - 1))
function _mars_state(planet, alt_m, rhat)
    r_ii = (planet.Rp_e + alt_m) * rhat
    pole = SVector(cos(planet.δ) * cos(planet.α), cos(planet.δ) * sin(planet.α), sin(planet.δ))
    ω_ii = norm(planet.ω) * pole
    return r_ii, ω_ii
end
# Explicit configuration with a bounded solver step (dt_max 0.05 s): the drag/lift/cross
# caches hold the last RHS-stage value, so bounding the step keeps the recorded force
# within a fraction of a scale height of the saved state during a fast descent.
function _config(planet, sc, effectors, mission_s, orientation)
    return SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false, save_csv=false),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=mission_s, orientation_sim=orientation, num_steps_to_save=300, data_rate=0.02),
        environment_model=SM.EnvironmentModel(planet=planet, EI=250.0, density_model=SM.ExponentialAtmosphereModel(planet),
            thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false, ephemerides_model=SM.SimpleEphemeridesModel()),
        dynamics_model=SM.DynamicsModel([sc], effectors),
        guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=SM.ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=SM.IntegrationTolerances(reltol_orbit=1e-10, abstol_orbit=1e-10,
            reltol_quaternion=1e-11, abstol_quaternion=1e-12, reltol_angular_rate=1e-11, abstol_angular_rate=1e-12,
            dt_max_orbit=0.05),
        solver_config=SM.SolverConfig(solver_mode=:tsit5))
end

@testset "MeshAero" begin
    @testset "panels: normals, areas, normalization defaults" begin
        tris = _box_triangles(2.0, 3.0, 4.0)
        panels = mesh_aero_panels(tris)
        @test length(panels) == 12
        @test sum(panels.areas) ≈ 2 * (2 * 3 + 2 * 4 + 3 * 4)
        @test panels.reference_area_m2 ≈ sum(panels.areas) / 4
        @test panels.reference_length_m ≈ norm([2.0, 3.0, 4.0])
        # every normal points away from the origin (outward)
        for k in 1:12
            @test dot(panels.normals[:, k], panels.centroids[:, k]) > 0
        end
        flipped = mesh_aero_panels(tris[:, [1, 3, 2, 4, 6, 5, 7:end...]]; outward_normals=true)
        @test all(dot(flipped.normals[:, k], flipped.centroids[:, k]) > 0 for k in 1:12)
        @test_throws ArgumentError mesh_aero_panels(tris[1:2, :])
        @test_throws ArgumentError mesh_aero_panels(tris; reference_area_m2=0.0)
        custom = mesh_aero_panels(tris; reference_area_m2=12.0, reference_length_m=3.0, moment_reference_m=(0.1, 0.0, 0.0))
        @test custom.reference_area_m2 == 12.0 && custom.reference_length_m == 3.0 && custom.moment_reference_m[1] == 0.1
    end

    @testset "sphere matches the closed-form drag" begin
        panels = mesh_aero_panels(_uv_sphere(1.0); reference_area_m2=pi)
        for s in (3.0, 8.0, 15.0)
            cf, cm = panel_aero_coefficients(panels, SVector(1.0, 0.0, 0.0), s; shadowing=false)
            @test isapprox(-cf[1], _sphere_cd(s); rtol=5e-3)
            @test hypot(cf[2], cf[3]) < 1e-3 * abs(cf[1])
            @test norm(cm) < 1e-3 * abs(cf[1])
            # shadowing drops the limb's back facets, whose thermal flux matters most at low speed ratio
            cfs, _ = panel_aero_coefficients(panels, SVector(0.0, 0.0, 1.0), s; shadowing=true)
            @test isapprox(-cfs[3], _sphere_cd(s); rtol=(s < 5 ? 2e-2 : 6e-3))
        end
        # wall temperature enters through sqrt(Tw/T) only
        cf_a, cm_a, cf_b, cm_b = panel_aero_coefficients_split(panels, SVector(1.0, 0.0, 0.0), 8.0; shadowing=false)
        cf_hot, _ = panel_aero_coefficients(panels, SVector(1.0, 0.0, 0.0), 8.0; tw_ratio=4.0, shadowing=false)
        @test cf_hot ≈ cf_a + 2.0 * cf_b
        @test isapprox(-cf_b[1], 2sqrt(pi) / (3 * 8.0); rtol=5e-3)
        @test panel_projected_area(panels, SVector(0.0, 0.0, 1.0)) ≈ pi rtol = 2e-2
    end

    @testset "box reproduces the Hart closed forms" begin
        lx, ly, lz = 1.5, 2.5, 0.8
        sigma = 0.85
        link = SM.Link(root=true, dims=MVector{3, Float64}(lx, ly, lz), ref_area=ly * lz, reflection_coefficient=sigma)
        panels = mesh_aero_panels(_box_triangles(lx, ly, lz); reference_area_m2=ly * lz)
        T = 180.0
        # near-normal rather than exactly normal: at exact grazing the Hart forms count one
        # side face of each pair where the panel method keeps both (see panel_shadow_mask)
        for s in (4.0, 9.0), (alpha, beta) in ((pi / 2 - 0.02, 0.015), (pi / 2 - 0.35, 0.2), (1.1, -0.6), (2.4, 0.4))
            vhat = _flow_from_angles(alpha, beta)
            cf, _ = panel_aero_coefficients(panels, vhat, s; sigma_n=sigma, sigma_t=sigma, tw_ratio=1.0, shadowing=true, grid=512)
            cl, cd, cs, _, _, _ = AE.aerodynamic_coefficient_fM(link, T, s, alpha, beta, 0.0)
            @test isapprox(-dot(cf, vhat), cd; rtol=1e-9)
            lateral = sqrt(max(0.0, dot(cf, cf) - dot(cf, vhat)^2))
            @test isapprox(lateral, hypot(cl, cs); rtol=1e-8, atol=1e-12)
        end
    end

    @testset "articulations pose parts of a mesh" begin
        # two plates with normal +x; the far one (x = 5) turns 90 deg about z through its own center,
        # so it ends up spanning x (1 m) and z with its normal along y
        front = _plate_triangles(1.0, 0.0)
        far = _plate_triangles(1.0, 5.0)
        tris = hcat(front, far)
        posed = articulate_triangles(tris, [(region=(x_min=4.0,), axis=(0.0, 0.0, 1.0), angle_deg=90.0)])
        @test posed[:, 1:6] == tris[:, 1:6]                                  # untouched part
        moved = posed[:, 7:12]
        @test maximum(moved[1, :]) - minimum(moved[1, :]) ≈ 1.0 atol = 1e-9   # now spans x
        @test maximum(abs.(moved[2, :])) < 1e-9                                 # thin in y (rotated about its own center)
        @test maximum(moved[3, :]) - minimum(moved[3, :]) ≈ 1.0 atol = 1e-9
        @test_throws ArgumentError articulate_triangles(tris, [(region=(x_min=100.0,), axis=(0.0, 0.0, 1.0), angle_deg=90.0)])
        @test_throws ArgumentError articulate_triangles(tris, [(axis=(0.0, 0.0, 1.0), angle_deg=90.0)])
        payload = articulation_payload([(region=(x_min=4.0,), axis=(0.0, 0.0, 2.0), angle_deg=90.0)], tris)
        @test payload[1]["axis"] ≈ [0.0, 0.0, 1.0] && payload[1]["pivot"] ≈ [5.0, 0.0, 0.0] && payload[1]["region"]["min"] == [4.0, nothing, nothing]
        # the panel method sees the posed geometry: the far plate now lies edge-on to +x flow
        p_open = mesh_aero_panels(tris; reference_area_m2=1.0)
        p_posed = mesh_aero_panels(posed; reference_area_m2=1.0)
        cf_open, _ = panel_aero_coefficients(p_open, SVector(1.0, 0.0, 0.0), 8.0; shadowing=false)
        cf_posed, _ = panel_aero_coefficients(p_posed, SVector(1.0, 0.0, 0.0), 8.0; shadowing=false)
        @test -cf_posed[1] < 0.6 * -cf_open[1]
    end

    @testset "shadowing hides a plate behind another" begin
        front = _plate_triangles(1.0, 1.0)
        back = _plate_triangles(1.0, 0.0)
        one = mesh_aero_panels(front; reference_area_m2=1.0)
        both = mesh_aero_panels(hcat(front, back); reference_area_m2=1.0)
        v = SVector(1.0, 0.0, 0.0)
        cf1, _ = panel_aero_coefficients(one, v, 8.0)
        cf2, _ = panel_aero_coefficients(both, v, 8.0; shadowing=true)
        cf2_open, _ = panel_aero_coefficients(both, v, 8.0; shadowing=false)
        @test cf2 ≈ cf1
        @test cf2_open ≈ 2 * cf1
        @test panel_shadow_mask(both, v) == BitVector([false, false, true, true])
        # offset the rear plate sideways so it is exposed again
        apart = mesh_aero_panels(hcat(front, _plate_triangles(1.0, 0.0; y0=3.0)); reference_area_m2=1.0)
        cf3, _ = panel_aero_coefficients(apart, v, 8.0; shadowing=true)
        @test cf3 ≈ 2 * cf1
    end

    @testset "spherical-harmonic basis" begin
        Y = zeros(AE._sh_count(6))
        # orthonormality over a Fibonacci quadrature (approximate): Gram matrix near identity
        dirs = AE.fibonacci_directions(4000)
        G = zeros(length(Y), length(Y))
        for v in dirs
            AE.real_sh_basis!(Y, 6, v)
            G .+= Y * Y'
        end
        G .*= 4pi / length(dirs)
        @test maximum(abs.(G - I)) < 0.05
        @test AE._sh_index(0, 0) == 1 && AE._sh_index(2, -2) == 5 && AE._sh_index(2, 2) == 9
        AE.real_sh_basis!(Y, 6, SVector(0.0, 0.0, 1.0))
        @test Y[1] ≈ 1 / sqrt(4pi)
        @test all(isfinite, Y)
    end

    @testset "surrogate fit, evaluation and JSON round trip" begin
        panels = mesh_aero_panels(_box_triangles(2.0, 3.0, 1.0); reference_area_m2=3.0)
        sur = fit_mesh_aero_surrogate(panels; degree=8, poly_degree=2, n_directions=700, speed_ratios=(4.0, 8.0, 16.0), grid=128, holdout_directions=60)
        @test sur.degree == 8 && sur.poly_degree == 2
        @test size(sur.coeff_a) == (6, AE._basis_count(8, 2))
        @test sur.speed_ratio_min == 4.0 && sur.speed_ratio_max == 16.0
        fit = sur.metadata["fit"]
        @test fit["samples"] == 3 * 700
        hold = fit["holdout"]
        # a box has edges the harmonics smooth over; the drag-direction components stay within a few percent
        for c in ("CFx", "CFy", "CFz")
            @test hold[c]["max_error"] <= 0.12 * max(hold[c]["scale"], 1e-9) + 0.02
        end
        rng = MersenneTwister(3)
        worst = 0.0
        for _ in 1:40
            v = normalize(SVector(randn(rng), randn(rng), randn(rng)))
            s = 4.0 + 12.0 * rand(rng)
            cf_true, cm_true = panel_aero_coefficients(panels, v, s; grid=128)
            cf, cm = mesh_aero_coefficients(sur, v, s)
            worst = max(worst, norm(cf - cf_true) / max(norm(cf_true), 1e-9))
        end
        @test worst < 0.15
        # extrapolation clamps to the fitted range
        v = SVector(1.0, 0.0, 0.0)
        @test mesh_aero_coefficients(sur, v, 100.0) == mesh_aero_coefficients(sur, v, 16.0)
        @test mesh_aero_coefficients(sur, v, 8.0; tw_ratio=1.0) != mesh_aero_coefficients(sur, v, 8.0; tw_ratio=2.0)
        path = joinpath(mktempdir(), "box.json")
        @test write_mesh_aero_surrogate(path, sur) == path
        back = read_mesh_aero_surrogate(path)
        @test back.coeff_a == sur.coeff_a && back.coeff_b == sur.coeff_b
        @test back.reference_area_m2 == sur.reference_area_m2 && back.sigma_n == sur.sigma_n
        @test back.moment_reference_m == sur.moment_reference_m
        @test back.metadata["facets"] == 12
        @test mesh_aero_coefficients(back, v, 6.0) == mesh_aero_coefficients(sur, v, 6.0)
        @test_throws ArgumentError read_mesh_aero_surrogate(joinpath(@__DIR__, "mesh_aero_tests.jl"))
        @test_throws ArgumentError fit_mesh_aero_surrogate(panels; degree=8, n_directions=10)
        @test_throws ArgumentError AerodynamicCoefficientMeshSurrogate(Dict{Int, MeshAeroSurrogate}())
        @test_throws ArgumentError AerodynamicCoefficientMeshSurrogate(sur; wall_temperature_k=0.0)
        model = AerodynamicCoefficientMeshSurrogate(sur; wall_temperature_k=250.0)
        @test model.surrogates[1] === sur && model.wall_temperature_k == 250.0
        @test SM.environment_requirements(model).atmosphere
        @test SM.solver_partition(model) == :implicit
    end

    @testset "moment reference translation" begin
        tris = _box_triangles(2.0, 3.0, 1.0)
        r0 = SVector(0.4, -0.3, 0.2)
        L = 2.5
        p0 = mesh_aero_panels(tris; reference_area_m2=3.0, reference_length_m=L)
        pr = mesh_aero_panels(tris; reference_area_m2=3.0, reference_length_m=L, moment_reference_m=r0)
        # panel method: M about r0 = M about origin - r0 x F, exactly
        for v in (normalize(SVector(1.0, 0.4, -0.2)), normalize(SVector(-0.3, 1.0, 0.7)))
            cf0, cm0 = panel_aero_coefficients(p0, v, 9.0; grid=128)
            cfr, cmr = panel_aero_coefficients(pr, v, 9.0; grid=128)
            @test cfr ≈ cf0 rtol=1e-12
            @test cmr ≈ cm0 - cross(r0, cf0) / L rtol=1e-10 atol=1e-14
        end
        fitkw = (degree=8, poly_degree=2, n_directions=500, speed_ratios=(4.0, 8.0, 16.0), grid=96, holdout_directions=0)
        s0 = fit_mesh_aero_surrogate(p0; fitkw...)
        sr = fit_mesh_aero_surrogate(pr; fitkw...)
        planet = make_no_gram_planet(:mars)
        rho, T = 2.0e-9, 190.0
        pos = (planet.Rp_e + 130e3) * normalize(SVector(0.6, 0.7, 0.3))
        vel = 3.4e3 * normalize(cross(pos, SVector(0.1, -0.2, 1.0)))

        @testset "root link: torque about the root is independent of the fitted reference point" begin
            bus = SM.Link(root=true, m=400.0, dims=MVector{3, Float64}(2.0, 3.0, 1.0), ref_area=3.0)
            sc = SM.SpacecraftModel(root=bus, id=1)
            q_ib = _quat_axis_angle((0.3, 1.0, -0.4), 0.9)
            x, env = _samples(planet, sc, q_ib, pos, vel, rho, T)
            f0, t0 = SM.wrench(AerodynamicCoefficientMeshSurrogate(s0), x, env, 0.0)
            fr, tr = SM.wrench(AerodynamicCoefficientMeshSurrogate(sr), x, env, 0.0)
            @test norm(f0) > 0 && norm(t0) > 0
            @test fr ≈ f0 rtol=1e-9
            @test tr ≈ t0 rtol=1e-9
            # the translation term is not negligible here, so ignoring the reference point would fail this
            f_root = SM.rot(q_ib) * fr
            @test norm(cross(r0, f_root)) > 0.1 * norm(tr)
            # hand composition from the surrogate itself: M_O = M_P + P x F
            vhat_link = SM.rot(q_ib) * normalize(vel)
            s = norm(vel) / sqrt(2 * planet.R * T)
            cf, cm = mesh_aero_coefficients(sr, vhat_link, s; tw_ratio=300.0 / T)
            q_dyn = 0.5 * rho * dot(vel, vel)
            f_link = q_dyn * 3.0 * cf
            @test fr ≈ SM.rot(q_ib)' * f_link rtol=1e-12
            @test tr ≈ q_dyn * 3.0 * L * cm + cross(r0, f_link) rtol=1e-12
        end

        @testset "whole-vehicle form follows the actual root, not link one" begin
            # SpacecraftModel appends the root when links initially contains
            # only children. A constant asymmetric wrench makes the wrong
            # child's rotation and lever arm directly observable.
            coeff = reshape(sqrt(4pi) .* [-1.0, 0.2, 0.1, 0.05, 0.1, -0.07], 6, 1)
            constant = MeshAeroSurrogate(0, 0, coeff, zeros(6, 1),
                1.0, 1.0, ZERO3, 1.0, 1.0, 3.0, 20.0)
            bus = SM.Link(root=true, m=300.0)
            child = SM.Link(root=false, m=100.0,
                r=MVector{3, Float64}(1.5, -0.4, 0.2),
                q=MVector{4, Float64}(_quat_axis_angle((0.0, 0.0, 1.0), pi / 2)))
            sc = SM.SpacecraftModel(root=bus, links=[child], id=97)
            @test sc.links[1] === child && sc.links[2] === bus
            whole = AerodynamicCoefficientMeshSurrogate(constant)
            actual_root = AerodynamicCoefficientMeshSurrogate(Dict(2 => constant))
            first_link = AerodynamicCoefficientMeshSurrogate(Dict(1 => constant))
            for q_ib in (_quat_axis_angle((0.3, 1.0, -0.4), 0.9), nothing)
                x, env = _samples(planet, sc, q_ib, pos, vel, rho, T)
                force, torque = SM.wrench(whole, x, env, 0.0)
                expected_force, expected_torque = SM.wrench(actual_root, x, env, 0.0)
                child_force, child_torque = SM.wrench(first_link, x, env, 0.0)
                @test force ≈ expected_force rtol=1e-13
                @test torque ≈ expected_torque rtol=1e-13
                @test norm(force - child_force) > 0.1 * norm(force)
                if q_ib !== nothing
                    @test norm(torque) > 0
                    @test norm(torque - child_torque) > 0.1 * norm(torque)
                else
                    @test torque == ZERO3
                end
            end
        end

        @testset "child link: reference point, link attitude and root-frame offset compose" begin
            bus = SM.Link(root=true, m=300.0, dims=MVector{3, Float64}(1.0, 1.0, 1.0), ref_area=1.0)
            q_child = _quat_axis_angle((0.0, 0.0, 1.0), pi / 2)
            r_child = SVector(1.5, -0.4, 0.2)
            child = SM.Link(root=false, m=100.0, r=MVector{3, Float64}(r_child), q=MVector{4, Float64}(q_child), dims=MVector{3, Float64}(2.0, 3.0, 1.0), ref_area=3.0)
            sc = SM.SpacecraftModel(root=bus, links=[bus, child], id=1)
            @test sc.links[2] === child
            q_ib = _quat_axis_angle((1.0, 0.2, 0.5), -1.3)
            x, env = _samples(planet, sc, q_ib, pos, vel, rho, T)
            model = AerodynamicCoefficientMeshSurrogate(Dict(2 => sr); wall_temperature_k=300.0)
            f, tau = SM.wrench(model, x, env, 0.0)
            # expected from the same surrogate with the documented composition
            R_root = SM.rot(q_ib)'              # root -> inertial
            R_c = SM.rot(q_child)'              # link -> root
            vhat_link = (R_root * R_c)' * normalize(vel)
            s = norm(vel) / sqrt(2 * planet.R * T)
            cf, cm = mesh_aero_coefficients(sr, vhat_link, s; tw_ratio=300.0 / T)
            q_dyn = 0.5 * rho * dot(vel, vel)
            f_link = q_dyn * 3.0 * cf
            m_link = q_dyn * 3.0 * L * cm + cross(r0, f_link)
            f_root = R_c * f_link
            @test f ≈ R_root * f_root rtol=1e-12
            @test tau ≈ R_c * m_link + cross(r_child, f_root) rtol=1e-12
            @test norm(cross(r_child, f_root)) > 0.1 * norm(tau)
            # the same torque from a whole-vehicle surrogate on the root, independent of the child fit
            @test SM.wrench(AerodynamicCoefficientMeshSurrogate(Dict(1 => sr)), x, env, 0.0)[2] ≈ q_dyn * 3.0 * L * mesh_aero_coefficients(sr, R_root' * normalize(vel), s; tw_ratio=300.0 / T)[2] + cross(r0, q_dyn * 3.0 * mesh_aero_coefficients(sr, R_root' * normalize(vel), s; tw_ratio=300.0 / T)[1]) rtol=1e-12
        end

        @testset "child link composition agrees with the panel method on the composed mesh" begin
            # A coarse sphere on a child link versus the same sphere placed by hand in the root
            # frame and integrated directly about the root origin: the lever-arm torque and the
            # force must agree to the surrogate's fit accuracy.
            sphere = _uv_sphere(0.5; nlat=16, nlon=32)
            ps = mesh_aero_panels(sphere; reference_area_m2=pi * 0.25, reference_length_m=1.0, moment_reference_m=(0.1, 0.0, -0.1))
            ss = fit_mesh_aero_surrogate(ps; degree=6, poly_degree=2, n_directions=300, speed_ratios=(4.0, 8.0, 16.0), grid=64, holdout_directions=0)
            q_child = _quat_axis_angle((0.2, 1.0, 0.3), 0.7)
            r_child = SVector(2.0, 0.5, -1.0)
            bus = SM.Link(root=true, m=300.0, dims=MVector{3, Float64}(1.0, 1.0, 1.0), ref_area=1.0)
            child = SM.Link(root=false, m=50.0, r=MVector{3, Float64}(r_child), q=MVector{4, Float64}(q_child), dims=MVector{3, Float64}(1.0, 1.0, 1.0), ref_area=1.0)
            sc = SM.SpacecraftModel(root=bus, links=[bus, child], id=1)
            q_ib = _quat_axis_angle((0.0, 1.0, 0.0), 0.4)
            x, env = _samples(planet, sc, q_ib, pos, vel, rho, T)
            f, tau = SM.wrench(AerodynamicCoefficientMeshSurrogate(Dict(2 => ss)), x, env, 0.0)
            R_c = SM.rot(q_child)'
            composed = Matrix{Float64}(undef, 3, size(sphere, 2))   # a loop: splatting thousands of SVectors into hcat compiles for ever
            for i in 1:size(sphere, 2)
                composed[:, i] = R_c * SVector{3, Float64}(sphere[:, i]) + r_child
            end
            pc = mesh_aero_panels(composed; reference_area_m2=pi * 0.25, reference_length_m=1.0)   # moment about the root origin
            vhat_root = SM.rot(q_ib) * normalize(vel)
            s = norm(vel) / sqrt(2 * planet.R * T)
            cf, cm = panel_aero_coefficients(pc, vhat_root, s; tw_ratio=300.0 / T, grid=128)
            q_dyn = 0.5 * rho * dot(vel, vel)
            @test SM.rot(q_ib) * f ≈ q_dyn * pi * 0.25 * cf rtol=3e-2
            @test tau ≈ q_dyn * pi * 0.25 * 1.0 * cm rtol=3e-2
            @test norm(tau) > 0.5 * norm(r_child) * norm(f)
        end
    end

    @testset "radial flight keeps the aerodynamic force" begin
        planet = make_no_gram_planet(:mars)
        sphere = _uv_sphere(1.0; nlat=16, nlon=32)
        ps = mesh_aero_panels(sphere; reference_area_m2=pi)
        ss = fit_mesh_aero_surrogate(ps; degree=6, poly_degree=2, n_directions=300, speed_ratios=(4.0, 8.0, 16.0), grid=64, holdout_directions=0)
        model = AerodynamicCoefficientMeshSurrogate(ss; wall_temperature_k=200.0)
        bus = SM.Link(root=true, m=400.0, dims=MVector{3, Float64}(2.0, 2.0, 2.0), ref_area=pi)
        sc = SM.SpacecraftModel(root=bus, id=1)
        rho, T = 3.0e-9, 200.0
        rhat = normalize(SVector(0.3, -0.5, 0.8))
        pos = (planet.Rp_e + 110e3) * rhat
        speed = 2.5e3
        s = speed / sqrt(2 * planet.R * T)
        expected_mag = 0.5 * rho * speed^2 * pi * _sphere_cd(s)
        for q_ib in (nothing, _quat_axis_angle((0.1, 0.4, 1.0), 0.8))
            # straight down: the orbit normal vanishes and no lift basis exists
            x, env = _samples(planet, sc, q_ib, pos, -speed * rhat, rho, T)
            f, tau, drag, lift, crs = AE._mesh_aero_wrench_sampled(model, x, env)
            @test isapprox(norm(f), expected_mag; rtol=3e-2)
            @test dot(normalize(f), rhat) > 0.999
            # the surrogate of a meshed sphere is not perfectly isotropic: a fitted lateral
            # component of order 1e-5 of the drag separates the projection from the force
            @test norm(drag - f) < 1e-3 * norm(f)
            @test lift == ZERO3 && crs == ZERO3
            @test SM.wrench(model, x, env, 0.0)[1] == f
            # straight up, and radial with a transverse component below the basis tolerance
            @test norm(AE._mesh_aero_wrench_sampled(model, _samples(planet, sc, q_ib, pos, speed * rhat, rho, T)...)[1]) ≈ norm(f) rtol=5e-3
            tiny = -speed * rhat + 1e-12 * speed * normalize(cross(rhat, SVector(0.0, 0.0, 1.0)))
            ft, _, _, lt, ct = AE._mesh_aero_wrench_sampled(model, _samples(planet, sc, q_ib, pos, tiny, rho, T)...)
            @test ft ≈ f rtol=1e-9
            @test lt == ZERO3 && ct == ZERO3
            # a horizontal pass at the same speed: same magnitude, a defined basis
            vh = speed * normalize(cross(rhat, SVector(0.0, 0.0, 1.0)))
            fh, _, dh, lh, ch = AE._mesh_aero_wrench_sampled(model, _samples(planet, sc, q_ib, pos, vh, rho, T)...)
            @test isapprox(norm(fh), norm(f); rtol=5e-3)
            @test dh + lh + ch ≈ fh rtol=1e-9
            @test norm(lh) + norm(ch) < 5e-2 * norm(dh)
        end
        # the caching hook stores the diagnostics for the spacecraft index
        args = make_example_config(planet=planet, spacecraft=sc, mission_time=1.0,
            initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
            dynamic_effectors=(SM.InverseSquaredGravityModel(), model), density_model=SM.ExponentialAtmosphereModel(planet),
            ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=false, verbose=false, results=false)
        p = SM.ODEParams(n_sats=1, args=args)
        SE._initialize_save_cache_buffers!(p)
        x, env = _samples(planet, sc, nothing, pos, -speed * rhat, rho, T)
        f, tau = SM.wrench_caching!(model, x, env, 0.0, p, 1)
        @test norm(p.save_cache.drag_cache[1] - f) < 1e-3 * norm(f)
        @test p.save_cache.lift_cache[1] == ZERO3 && p.save_cache.cross_cache[1] == ZERO3
        @test tau == ZERO3   # fixed attitude: no torque, as for the other aerodynamic effectors
    end

    @testset "malformed inputs are rejected before evaluation" begin
        nb(L, P) = AE._basis_count(L, P)
        good(L=8, P=2) = (L, P, zeros(6, nb(L, P)), zeros(6, nb(L, P)), 1.0, 1.0, zeros(3), 1.0, 1.0, 3.0, 20.0)
        @test MeshAeroSurrogate(good()...) isa MeshAeroSurrogate
        @test MESH_AERO_MAX_DEGREE == 20
        @test MeshAeroSurrogate(good(20, 0)...).degree == 20
        # Addition overflow and multiplication-to-zero must not admit an
        # empty matrix into the evaluator. Never evaluate these invalid inputs.
        for huge_poly in (typemax(Int), typemax(Int) ÷ 2)
            err = try
                MeshAeroSurrogate(1, huge_poly, zeros(6, 0), zeros(6, 0),
                    1.0, 1.0, ZERO3, 1.0, 1.0, 3.0, 20.0)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError && occursin("basis size", err.msg)
        end
        # every constructor route validates: exact field types cannot bypass it
        @test_throws ArgumentError MeshAeroSurrogate(good(21, 2)...)
        @test_throws ArgumentError MeshAeroSurrogate(21, 2, zeros(6, nb(21, 2)), zeros(6, nb(21, 2)), 1.0, 1.0, SVector(0.0, 0.0, 0.0), 1.0, 1.0, 3.0, 20.0, Dict{String, Any}())
        @test length(methods(MeshAeroSurrogate)) == 2   # the validating constructor with and without metadata
        @test_throws ArgumentError MeshAeroSurrogate(-1, 2, zeros(6, 0), zeros(6, 0), 1.0, 1.0, zeros(3), 1.0, 1.0, 3.0, 20.0)
        @test_throws ArgumentError MeshAeroSurrogate(8, -1, zeros(6, 0), zeros(6, 0), 1.0, 1.0, zeros(3), 1.0, 1.0, 3.0, 20.0)
        @test_throws ArgumentError MeshAeroSurrogate(8, 2, zeros(6, nb(8, 2) - 1), zeros(6, nb(8, 2)), 1.0, 1.0, zeros(3), 1.0, 1.0, 3.0, 20.0)
        @test_throws ArgumentError MeshAeroSurrogate(8, 2, zeros(5, nb(8, 2)), zeros(6, nb(8, 2)), 1.0, 1.0, zeros(3), 1.0, 1.0, 3.0, 20.0)
        bad = zeros(6, nb(8, 2)); bad[3, 4] = NaN
        @test_throws ArgumentError MeshAeroSurrogate(8, 2, bad, zeros(6, nb(8, 2)), 1.0, 1.0, zeros(3), 1.0, 1.0, 3.0, 20.0)
        @test_throws ArgumentError MeshAeroSurrogate(8, 2, zeros(6, nb(8, 2)), zeros(6, nb(8, 2)), 0.0, 1.0, zeros(3), 1.0, 1.0, 3.0, 20.0)
        @test_throws ArgumentError MeshAeroSurrogate(8, 2, zeros(6, nb(8, 2)), zeros(6, nb(8, 2)), 1.0, Inf, zeros(3), 1.0, 1.0, 3.0, 20.0)
        @test_throws ArgumentError MeshAeroSurrogate(8, 2, zeros(6, nb(8, 2)), zeros(6, nb(8, 2)), 1.0, 1.0, zeros(2), 1.0, 1.0, 3.0, 20.0)
        @test_throws ArgumentError MeshAeroSurrogate(8, 2, zeros(6, nb(8, 2)), zeros(6, nb(8, 2)), 1.0, 1.0, [0.0, NaN, 0.0], 1.0, 1.0, 3.0, 20.0)
        @test_throws ArgumentError MeshAeroSurrogate(8, 2, zeros(6, nb(8, 2)), zeros(6, nb(8, 2)), 1.0, 1.0, zeros(3), 1.5, 1.0, 3.0, 20.0)
        @test_throws ArgumentError MeshAeroSurrogate(8, 2, zeros(6, nb(8, 2)), zeros(6, nb(8, 2)), 1.0, 1.0, zeros(3), 1.0, -0.1, 3.0, 20.0)
        @test_throws ArgumentError MeshAeroSurrogate(8, 2, zeros(6, nb(8, 2)), zeros(6, nb(8, 2)), 1.0, 1.0, zeros(3), 1.0, 1.0, 20.0, 3.0)
        @test_throws ArgumentError MeshAeroSurrogate(8, 2, zeros(6, nb(8, 2)), zeros(6, nb(8, 2)), 1.0, 1.0, zeros(3), 1.0, 1.0, 0.0, 3.0)
        # fitting checks the limits before tabulating anything
        panels = mesh_aero_panels(_box_triangles(1.0, 1.0, 1.0))
        elapsed = @elapsed @test_throws ArgumentError fit_mesh_aero_surrogate(panels; degree=21, n_directions=100000)
        @test elapsed < 5.0
        @test_throws ArgumentError fit_mesh_aero_surrogate(panels; poly_degree=-1)
        # Enough total samples cannot compensate for an undetermined basis factor.
        @test_throws ArgumentError fit_mesh_aero_surrogate(panels; degree=1, poly_degree=2, n_directions=100, speed_ratios=(4.0, 8.0))
        @test_throws ArgumentError fit_mesh_aero_surrogate(panels; degree=1, poly_degree=2, n_directions=100, speed_ratios=(4.0, 4.0, 8.0))
        @test_throws ArgumentError fit_mesh_aero_surrogate(panels; degree=2, poly_degree=0, n_directions=8, speed_ratios=(4.0, 8.0))
        @test_throws ArgumentError fit_mesh_aero_surrogate(panels; sigma_n=1.2, n_directions=50, degree=2)
        @test_throws ArgumentError fit_mesh_aero_surrogate(panels; speed_ratios=(0.0, 4.0), n_directions=50, degree=2)
        @test_throws ArgumentError fit_mesh_aero_surrogate(panels; speed_ratios=(), n_directions=50, degree=2)
        # effector construction
        sur = MeshAeroSurrogate(good()...)
        @test_throws ArgumentError AerodynamicCoefficientMeshSurrogate(Dict(0 => sur))
        @test_throws ArgumentError AerodynamicCoefficientMeshSurrogate(Dict("1" => sur))
        @test_throws ArgumentError AerodynamicCoefficientMeshSurrogate(Dict(1.0 => sur))
        @test_throws ArgumentError AerodynamicCoefficientMeshSurrogate(Dict(1 => panels))
        @test_throws ArgumentError AerodynamicCoefficientMeshSurrogate(sur; wall_temperature_k=NaN)
        @test_throws ArgumentError AerodynamicCoefficientMeshSurrogate(sur; wall_temperature_k=Inf)
        @test_throws ArgumentError AerodynamicCoefficientMeshSurrogate(sur; wall_temperature_k=-1.0)
        @test AerodynamicCoefficientMeshSurrogate(Dict(3 => sur, 1 => sur)).surrogates[3] === sur
        # the file route validates every field before the evaluator can see the surrogate
        dir = mktempdir()
        fitted = fit_mesh_aero_surrogate(panels; degree=4, poly_degree=1, n_directions=120, speed_ratios=(4.0, 12.0), grid=64, holdout_directions=0)
        path = write_mesh_aero_surrogate(joinpath(dir, "cube.json"), fitted)
        doc = JSON.parsefile(path)
        function written(mutate!)
            d = deepcopy(doc)
            mutate!(d)
            p = joinpath(dir, "variant.json")
            open(p, "w") do io
                JSON.print(io, d)
            end
            return p
        end
        @test read_mesh_aero_surrogate(written(d -> nothing)).degree == 4
        for huge_poly in (typemax(Int), typemax(Int) ÷ 2)
            overflowing = written(d -> begin
                d["degree"] = 1
                d["poly_degree"] = huge_poly
                d["coeff_a"] = [Float64[] for _ in 1:6]
                d["coeff_b"] = [Float64[] for _ in 1:6]
            end)
            err = try; read_mesh_aero_surrogate(overflowing); nothing; catch e; e; end
            @test err isa ArgumentError && occursin("variant.json", err.msg) && occursin("basis size", err.msg)
        end
        oversized = written(d -> (d["degree"] = 25; d["coeff_a"] = [zeros(nb(25, 1)) for _ in 1:6]; d["coeff_b"] = [zeros(nb(25, 1)) for _ in 1:6]))
        err = try; read_mesh_aero_surrogate(oversized); nothing; catch e; e; end
        @test err isa ArgumentError && occursin("variant.json", err.msg) && occursin("degree", err.msg)
        @test_throws ArgumentError read_mesh_aero_surrogate(written(d -> d["degree"] = 4.5))
        @test_throws ArgumentError read_mesh_aero_surrogate(written(d -> d["degree"] = "4"))
        @test_throws ArgumentError read_mesh_aero_surrogate(written(d -> d["poly_degree"] = -1))
        @test_throws ArgumentError read_mesh_aero_surrogate(written(d -> delete!(d, "coeff_b")))
        @test_throws ArgumentError read_mesh_aero_surrogate(written(d -> d["coeff_a"] = d["coeff_a"][1:5]))
        @test_throws ArgumentError read_mesh_aero_surrogate(written(d -> d["coeff_a"][2] = d["coeff_a"][2][1:end-1]))
        @test_throws ArgumentError read_mesh_aero_surrogate(written(d -> d["coeff_a"][1][1] = "x"))
        @test_throws ArgumentError read_mesh_aero_surrogate(written(d -> d["coeff_a"][1][1] = nothing))
        @test_throws ArgumentError read_mesh_aero_surrogate(written(d -> d["moment_reference_m"] = [0.0, 0.0]))
        @test_throws ArgumentError read_mesh_aero_surrogate(written(d -> d["reference_area_m2"] = -1.0))
        @test_throws ArgumentError read_mesh_aero_surrogate(written(d -> d["reference_length_m"] = "1"))
        @test_throws ArgumentError read_mesh_aero_surrogate(written(d -> d["sigma_t"] = 2.0))
        @test_throws ArgumentError read_mesh_aero_surrogate(written(d -> d["speed_ratio_min"] = 40.0))
        @test_throws ArgumentError read_mesh_aero_surrogate(written(d -> d["metadata"] = [1, 2]))
        @test_throws ArgumentError read_mesh_aero_surrogate(written(d -> d["schema"] = "something_else"))
        @test_throws ArgumentError read_mesh_aero_surrogate(joinpath(dir, "absent.json"))
    end

    @testset "engine: propagated attitude integrates the mesh torque, whatever the fitted reference point" begin
        planet = make_no_gram_planet(:mars)
        # a one-sided plate one metre toward body +z of the centre of mass, facing body +x
        tris = _plate_triangles(1.0, 0.0; z0=1.0)
        fitkw = (degree=10, poly_degree=2, n_directions=900, speed_ratios=(4.0, 8.0, 12.0, 24.0), grid=96, holdout_directions=0)
        s0 = fit_mesh_aero_surrogate(mesh_aero_panels(tris; reference_area_m2=1.0, reference_length_m=1.0); fitkw...)
        r0 = SVector(0.3, -0.2, 0.5)
        sr = fit_mesh_aero_surrogate(mesh_aero_panels(tris; reference_area_m2=1.0, reference_length_m=1.0, moment_reference_m=r0); fitkw...)
        rhat = SVector(1.0, 0.0, 0.0)
        r_ii, ω_ii = _mars_state(planet, 120e3, rhat)
        v_ii = sqrt(planet.μ / norm(r_ii)) * normalize(SVector(0.0, cos(0.5), sin(0.5)))
        v_rel = v_ii - cross(ω_ii, r_ii)                       # planet-relative velocity in inertial axes
        x_b = normalize(v_rel)
        z_b = normalize(-r_ii - dot(-r_ii, x_b) * x_b)         # nadir, orthogonalised
        q0 = _quat_from_rows(SMatrix{3, 3, Float64, 9}(vcat(x_b', cross(z_b, x_b)', z_b')))
        @test SM.rot(q0) * x_b ≈ E1 atol=1e-12
        inertia = SMatrix{3, 3, Float64}(400.0 / 12 * 2.0 * I)
        function build(sur)
            bus = SM.Link(root=true, m=400.0, dims=MVector{3, Float64}(1.0, 1.0, 1.0), ref_area=1.0, inertia=inertia)
            ic = SM.CartesianInitialCondition(r_ii, v_ii; q=q0, ang_vel=ZERO3)
            sc = SM.SpacecraftModel(root=bus, initial_condition=ic, inertia_tensor=inertia, id=1)
            return _config(planet, sc, (SM.InverseSquaredGravityModel(), AerodynamicCoefficientMeshSurrogate(sur; wall_temperature_k=200.0)), 6.0, true)
        end
        fields = vcat(STATE_FIELDS, SaveField[
            SaveField(:q, (u, t, integrator) -> [SVector{4, Float64}(sc.q) for sc in u.sc]; per_satellite=true),
            SaveField(:omega, (u, t, integrator) -> [SVector{3, Float64}(sc.ω) for sc in u.sc]; per_satellite=true)])
        t0s, snap0 = _record(build(s0), fields)
        trs, snapr = _record(build(sr), fields)
        @test length(t0s) > 50 && t0s == trs
        # (i) the two fits differ only in the reference point: identical motion
        for (a, b) in zip(snap0, snapr)
            @test a[:q][1] ≈ b[:q][1] atol=1e-7
            @test a[:omega][1] ≈ b[:omega][1] atol=1e-7
            @test a[:pos][1] ≈ b[:pos][1] atol=1e-3
        end
        # (ii) the rate builds up about body -y at the rate the plate torque predicts. The torque
        # per unit drag follows from the surrogate at the initial incidence (body x along the
        # airspeed), so the expectation needs no independent density model.
        s_mid = norm(v_rel) / sqrt(2 * planet.R * 200.0)
        cf, cm = mesh_aero_coefficients(sr, E1, s_mid; tw_ratio=1.0)
        torque_per_drag = (cm + cross(r0, cf)) / (-dot(cf, E1))
        drag_mag = [norm(sn[:drag][1]) for sn in snapr]
        @test drag_mag[2] > 0          # the first save precedes the first RHS evaluation
        impulse = _trapz(trs, drag_mag)
        expected_domega = inertia \ (torque_per_drag * impulse)
        omega_end = snapr[end][:omega][1]
        @test omega_end[2] < 0
        # pitch rate: the plate torque integrated over the recorded drag impulse. The planet-frame
        # airspeed is a few degrees off body x at the start (the engine's planet frame at the initial
        # time is not exactly the pole-aligned frame assumed here), which feeds a few percent of the
        # pitch rate into the other axes; those stay bounded.
        @test isapprox(omega_end[2], expected_domega[2]; rtol=3e-2)
        @test max(abs(omega_end[1]), abs(omega_end[3])) < 0.05 * abs(omega_end[2])
        @test norm(omega_end) > 1e-4
        # the diagnostics are projections of one force: drag, lift and cross add up along the run
        # (the first save precedes the first RHS evaluation and holds empty caches)
        for sn in snapr[2:end]
            @test norm(sn[:lift][1]) + norm(sn[:cross][1]) < 0.2 * norm(sn[:drag][1])
        end
    end

    @testset "engine: radial descent keeps drag and reports an undefined lift basis as zero" begin
        planet = make_no_gram_planet(:mars)
        sphere = _uv_sphere(2.0; nlat=16, nlon=32)
        ss = fit_mesh_aero_surrogate(mesh_aero_panels(sphere; reference_area_m2=4pi); degree=6, poly_degree=2, n_directions=300, speed_ratios=(4.0, 8.0, 16.0), grid=64, holdout_directions=0)
        rhat = normalize(SVector(0.2, 0.9, 0.4))
        r_ii, ω_ii = _mars_state(planet, 112e3, rhat)
        v_ii = cross(ω_ii, r_ii) - 2.0e3 * rhat                   # co-rotating plus straight down
        function build(with_drag)
            bus = SM.Link(root=true, m=400.0, dims=MVector{3, Float64}(4.0, 4.0, 4.0), ref_area=4pi)
            sc = SM.SpacecraftModel(root=bus, initial_condition=SM.CartesianInitialCondition(r_ii, v_ii), id=1)
            effectors = with_drag ? (SM.InverseSquaredGravityModel(), AerodynamicCoefficientMeshSurrogate(ss; wall_temperature_k=200.0)) : (SM.InverseSquaredGravityModel(),)
            return _config(planet, sc, effectors, 6.0, false)
        end
        t, snaps = _record(build(true), STATE_FIELDS)
        tb, base = _record(build(false), STATE_FIELDS)
        energy(sn) = 0.5 * dot(sn[:vel][1], sn[:vel][1]) - planet.μ / norm(sn[:pos][1])
        # the first save precedes the first RHS evaluation (empty caches); from the second on the
        # descent is radial to within the Coriolis drift, and drag acts throughout
        @test all(sn -> norm(sn[:drag][1]) > 0, snaps[2:end])
        # the recorded force does work on the inertial motion; the energy loss matches it
        power = [dot(sn[:drag][1] + sn[:lift][1] + sn[:cross][1], sn[:vel][1]) / 400.0 for sn in snaps[2:end]]
        work = _trapz(t[2:end], power)
        @test work < 0
        @test isapprox(energy(snaps[end]) - energy(snaps[2]), work; rtol=5e-2)
        @test abs(energy(base[end]) - energy(base[1])) < 1e-3 * abs(work)
        # a sphere's force stays along the airspeed, so the lift/cross split stays small once the basis exists
        for sn in snaps[2:end]
            @test norm(sn[:lift][1]) + norm(sn[:cross][1]) < 5e-2 * norm(sn[:drag][1])
        end
    end

    @testset "effector in a simulation tracks the fM box" begin
        planet = make_no_gram_planet(:mars)
        side = 2.0
        ic = SM.InitialCondition(ra=planet.Rp_e + 120e3, rp=planet.Rp_e + 120e3, i=30.0, ω=0.0, Ω=45.0, ν=0.0)
        function config(effector, dir)
            bus = SM.Link(root=true, m=400.0, dims=MVector{3, Float64}(side, side, side), ref_area=side^2, reflection_coefficient=1.0)
            sc = SM.SpacecraftModel(root=bus, initial_condition=ic, id=1)
            return make_example_config(planet=planet, spacecraft=sc, mission_time=300.0,
                initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
                dynamic_effectors=(SM.InverseSquaredGravityModel(), effector), density_model=SM.ExponentialAtmosphereModel(planet),
                ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=false, keplerian=true, EI_km=250.0,
                verbose=false, results=true, results_directory=dir)
        end
        function energy_drop(dir)
            df = DataFrame(Arrow.Table(joinpath(dir, "simulation_results.feather")))
            e(i) = 0.5 * (df.sc1_vel_1[i]^2 + df.sc1_vel_2[i]^2 + df.sc1_vel_3[i]^2) - planet.μ / hypot(df.sc1_pos_1[i], df.sc1_pos_2[i], df.sc1_pos_3[i])
            return e(1) - e(nrow(df))
        end
        panels = mesh_aero_panels(_box_triangles(side, side, side); reference_area_m2=side^2)
        sur = fit_mesh_aero_surrogate(panels; degree=10, poly_degree=2, n_directions=900, speed_ratios=(3.0, 6.0, 12.0, 24.0), grid=96, holdout_directions=0)
        # the fit must be good where the run sits: axial flow on the cube
        cf_axial, _ = panel_aero_coefficients(panels, SVector(1.0, 0.0, 0.0), 14.0; grid=96)
        cf_fit, _ = mesh_aero_coefficients(sur, SVector(1.0, 0.0, 0.0), 14.0)
        @test isapprox(cf_fit[1], cf_axial[1]; rtol=0.06)   # degree-10 harmonics ring a few percent at a cube's edges
        dir_fm = mktempdir(); dir_mesh = mktempdir()
        run_simulation(config(SM.AerodynamicCoefficientfM(), dir_fm); visualization=false)
        # Tw = T for the Hart form; the exponential model reports one fixed temperature
        T_atm = SM.ExponentialAtmosphereModel(planet).temperature_k
        run_simulation(config(AerodynamicCoefficientMeshSurrogate(sur; wall_temperature_k=T_atm), dir_mesh); visualization=false)
        drop_fm = energy_drop(dir_fm)
        drop_mesh = energy_drop(dir_mesh)
        @test drop_fm > 0
        # the cube in axial flow gets shear on all four sides here and on two in the Hart forms:
        # 4 sigma/(s sqrt(pi)) versus 2, a few percent of CD at s ~ 15, so the mesh drags a little more
        @test drop_mesh > 0.98 * drop_fm
        @test isapprox(drop_mesh, drop_fm; rtol=0.1)
        df = DataFrame(Arrow.Table(joinpath(dir_mesh, "simulation_results.feather")))
        @test "sc1_drag_1" in names(df)
        @test maximum(hypot.(df.sc1_drag_1, df.sc1_drag_2, df.sc1_drag_3)) > 0
    end
end
end # module
