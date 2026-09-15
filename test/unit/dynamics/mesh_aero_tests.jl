using Test
using SpaceAGORA
using StaticArrays
using LinearAlgebra
using SpecialFunctions
using DataFrames
using Arrow
using Random

import SpaceAGORA.TelemetryVerification: make_example_config

const SM = SpaceAGORA.SimulationModel
const AE = SM.DynamicEffectors.AerodynamicEffectors

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
        # file path form applies the viewer transform
        magellan = joinpath(@__DIR__, "..", "..", "..", "data", "models", "magellan_nasa_3d_resources.glb")
        if isfile(magellan)
            pm = mesh_aero_panels(magellan; scale=1.0, rotation_deg=(0, 0, 90))
            @test length(pm) > 1000
            @test pm.source == "magellan_nasa_3d_resources.glb"
        end
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
        magellan = joinpath(@__DIR__, "..", "..", "..", "data", "models", "magellan_nasa_3d_resources.glb")
        if isfile(magellan)
            wings = ((region=(x_min=1.9, y_max=1.0), axis=(1.0, 0.0, 0.0), angle_deg=-43.5), (region=(x_max=-1.9, y_max=1.0), axis=(1.0, 0.0, 0.0), angle_deg=-43.5))
            raw = load_model_triangles(magellan)
            posed = load_model_triangles(magellan; articulations=wings)
            sel = raw[1, :] .> 1.9 .&& raw[2, :] .< 1.0
            @test maximum(posed[2, sel]) - minimum(posed[2, sel]) < 0.2        # wing now thin along the flow axis
            @test maximum(posed[3, sel]) - minimum(posed[3, sel]) > 2.5        # and spans model z
        end
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
        run_simulation(config(SM.AerodynamicCoefficientfM(), dir_fm))
        # Tw = T for the Hart form; the exponential model reports one fixed temperature
        T_atm = SM.ExponentialAtmosphereModel(planet).temperature_k
        run_simulation(config(AerodynamicCoefficientMeshSurrogate(sur; wall_temperature_k=T_atm), dir_mesh))
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
