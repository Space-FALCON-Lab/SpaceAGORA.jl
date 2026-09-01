using Test
using LinearAlgebra
using StaticArrays

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
using .SimulationModel

const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
const EARTH = Earth("", SPICE_PATH)

const CB = SimulationModel.SimulationCallbacks
const EM = SimulationModel.EnvironmentModels

# ── Probe density models (no native GRAM required) ──────────────────────────

# Fake GRAM "core" so EnvironmentModels.GRAMAtmosphereModel instances can be
# constructed and evaluated without the native GRAM library being built.
struct FakeGramCore
    rho::Float64
    T::Float64
end

function EM._gram_core_density_state(
    core::FakeGramCore,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    lock_obj,
    vacuum_temperature::Float64
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return core.rho, core.T, SVector{3, Float64}(1.0, 2.0, 3.0)
end

make_fake_gram(; rho::Float64=1.5e-9, T::Float64=210.0) =
    EM.GRAMAtmosphereModel(FakeGramCore(rho, T))

# Smooth deterministic density model for spline-cache tests.
struct SplineProbeDensityModel <: SimulationModel.AbstractDensityModel end

function EM.getDensity(
    ::SplineProbeDensityModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p
)::Tuple{Float64, Float64, SVector{3, Float64}}
    rho = 1.0e-9 * exp(-h / 2.0e6) * (1.0 + 0.05 * sin(el_time / 200.0))
    T = 180.0 + 1.0e-5 * h
    return rho, T, SVector{3, Float64}(1.0, -2.0, 0.5)
end

# Model whose scalar query works but whose batch query throws, so the
# refresh! failure path is reachable while the direct fallback still succeeds.
struct BatchThrowDensityModel <: SimulationModel.AbstractDensityModel end

function EM.getDensity(
    ::BatchThrowDensityModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return 1.23e-7, 199.0, SVector{3, Float64}(0.25, 0.5, 0.75)
end

function EM.getDensityBatch!(
    rhos::AbstractVector{Float64},
    Ts::AbstractVector{Float64},
    winds::AbstractVector{SVector{3, Float64}},
    ::BatchThrowDensityModel,
    hs::AbstractVector{<:Real},
    lats::AbstractVector{<:Real},
    lons::AbstractVector{<:Real},
    el_time::Union{Float64, AbstractVector{<:Real}},
    wind::Bool,
    p
)::Nothing
    error("batch_throw_probe")
end

# Core shaped like GRAMSuite.GRAMAtmosphereModel: exposes `gram` (a trajectory
# driver) and `gram_atmosphere`, so the generate_trajectory branch of the
# track-cache refresh is reached through the wrappers' real property
# forwarding, without the native GRAM library.
struct FakeTrajCore{G, GA}
    gram::G
    gram_atmosphere::GA
end

const FAKE_TRAJ_KWARGS = Ref{Any}(nothing)

function fake_generate_trajectory(gram_atmosphere; kwargs...)
    kw = Dict(kwargs)
    FAKE_TRAJ_KWARGS[] = kw
    n = kw[:n_points]::Int
    pts = Vector{Any}(undef, n)
    for k in 1:n
        # Longitudes deliberately straddle ±180° so both wrap branches of
        # _gram_track_cache_fill_from_trajectory! execute.
        lon_deg = k == 1 ? 190.0 : (k == 2 ? -190.0 : 10.0 * k)
        pts[k] = (
            position = (
                elapsedTime = kw[:initial_elapsed_time] + (k - 1) * kw[:delta_elapsed_time],
                height = kw[:initial_height] + (k - 1) * kw[:delta_height],
                latitude = 5.0 + 0.1 * k,
                longitude = lon_deg,
            ),
            dynamics = (density = 1.0e-8 * k, temperature = 200.0 + k),
            winds = (
                perturbedEWWind = 1.0 * k,
                perturbedNSWind = -1.0 * k,
                perturbedVerticalWind = 0.5 * k,
            ),
        )
    end
    return pts
end

# Module driver, mirroring the native GRAM Julia wrapper: GRAMSuite stores the
# included GRAM module in the core's `gram` field. `generate_trajectory` is
# deliberately not exported — the support check must probe modules with
# `isdefined`, not `hasproperty` (which only sees exported module names).
module FakeTrajDriver end
Base.eval(FakeTrajDriver, :(const generate_trajectory = $fake_generate_trajectory))

make_traj_gram() = EM.GRAMAtmosphereModel(FakeTrajCore(FakeTrajDriver, nothing))

# ── Config / state builders ──────────────────────────────────────────────────

function make_probe_spacecraft(;
    root_area::Float64=12.0,
    panel_area::Float64=6.0,
    links_override=nothing
)
    root = Link(root=true, m=500.0, ref_area=root_area)
    panel = Link(root=false, m=30.0, ref_area=panel_area, r=MVector{3, Float64}(0.0, 1.2, 0.0))
    ic = InitialCondition(
        ra=EARTH.Rp_e + 500e3,
        rp=EARTH.Rp_e + 450e3,
        i=35.0,
        ω=40.0,
        Ω=10.0,
        ν=175.0
    )
    links = links_override === nothing ? [root, panel] : links_override
    dry_mass = root.m + panel.m
    return SpacecraftModel(
        Joint[],
        links,
        root,
        true,
        dry_mass,
        0.0,
        root.inertia,
        0,
        0,
        ic,
        1
    )
end

function probe_config(;
    density_model,
    spacecraft::Vector{SpacecraftModel}=[make_probe_spacecraft()],
    keplerian::Bool=true,
    EI_km::Float64=120.0,
    mission_type=MissionTime
)
    environment_model = EnvironmentModel(
        planet=EARTH,
        EI=EI_km,
        density_model=density_model,
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
        topography=false,
        wind=false
    )
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=MissionConfiguration(
            mission_type=mission_type,
            keplerian=keplerian,
            number_of_orbits=1,
            mission_time=600.0,
            orientation_sim=false,
            num_steps_to_save=1000
        ),
        environment_model=environment_model,
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances()
    )
end

# Inertial states on well-defined conic orbits.
function circular_state(alt_m::Float64)
    r = EARTH.Rp_e + alt_m
    pos = SVector{3, Float64}(r, 0.0, 0.0)
    v = sqrt(EARTH.μ / r)
    vel = SVector{3, Float64}(0.0, v * cosd(30.0), v * sind(30.0))
    return pos, vel
end

function elliptical_apoapsis_state(rp_alt_m::Float64, ra_alt_m::Float64)
    rp = EARTH.Rp_e + rp_alt_m
    ra = EARTH.Rp_e + ra_alt_m
    a = 0.5 * (rp + ra)
    pos = SVector{3, Float64}(ra, 0.0, 0.0)
    v = sqrt(EARTH.μ * (2.0 / ra - 1.0 / a))
    vel = SVector{3, Float64}(0.0, v, 0.0)
    return pos, vel
end

# State with planet-relative speed v and heading angle χ (from north, toward
# east) at a given geocentric latitude/longitude/altitude. Adds back the
# planet-rotation velocity so the relative NED decomposition is exact.
function ned_probe_state(lat_deg::Float64, lon_deg::Float64, alt_m::Float64, v::Float64, chi::Float64)
    lat = deg2rad(lat_deg)
    lon = deg2rad(lon_deg)
    u_up = SVector{3, Float64}(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
    u_e = SVector{3, Float64}(-sin(lon), cos(lon), 0.0)
    u_n = SVector{3, Float64}(-sin(lat) * cos(lon), -sin(lat) * sin(lon), cos(lat))
    pos = (EARTH.Rp_e + alt_m) * u_up
    vel = v * cos(chi) * u_n + v * sin(chi) * u_e + SVector{3, Float64}(cross(EARTH.ω, pos))
    return SVector{3, Float64}(pos), SVector{3, Float64}(vel)
end

function hyperbolic_state(alt_m::Float64)
    r = EARTH.Rp_e + alt_m
    pos = SVector{3, Float64}(r, 0.0, 0.0)
    v_esc = sqrt(2.0 * EARTH.μ / r)
    vel = SVector{3, Float64}(0.0, 1.2 * v_esc, 0.0)
    return pos, vel
end

@testset "Density selection / GRAM cache coverage probes" begin
    # Legacy 3-arg r_intor_p! (used by targeting) reads planet.L_PI, which is
    # all-zeros on a fresh Earth. Pin it to identity for deterministic geometry.
    EARTH.L_PI .= Matrix{Float64}(I, 3, 3)

    fallback_probe = SplineProbeDensityModel()
    gram_a = make_fake_gram(rho=1.0e-9, T=205.0)
    gram_b = make_fake_gram(rho=2.0e-9, T=215.0)

    cfg_gram = probe_config(density_model=make_fake_gram())
    cfg_nongram = probe_config(density_model=fallback_probe)
    cfg_nonkeplerian = probe_config(density_model=make_fake_gram(), keplerian=false)

    @testset "model_selection: per-sat and batch model resolution" begin
        models = EM.GRAMAtmosphereModel[gram_a, gram_b]
        @test CB._density_model_for_sat(models, fallback_probe, 1) === gram_a
        @test CB._density_model_for_sat(models, fallback_probe, 2) === gram_b
        @test CB._density_model_for_sat(models, fallback_probe, 3) === fallback_probe

        p = ODEParams(n_sats=2, args=cfg_gram)
        @test CB._density_model_for_sat(p, 1) === cfg_gram.environment_model.density_model
        push!(p.shared_buffers.density_models, gram_a, gram_b)
        @test CB._density_model_for_sat(p, 2) === gram_b

        # Batch model: empty vector -> shared fallback model.
        empty_models = SimulationModel.AbstractDensityModel[]
        @test CB._density_batch_model_for_callback(empty_models, fallback_probe, 3) === fallback_probe
        # Fewer per-sat models than satellites -> no batch model.
        @test CB._density_batch_model_for_callback(models, fallback_probe, 3) === nothing
        # Mixed per-sat models -> no batch model.
        @test CB._density_batch_model_for_callback(models, fallback_probe, 2) === nothing
        # Uniform per-sat models -> that model.
        uniform = EM.GRAMAtmosphereModel[gram_a, gram_a]
        @test CB._density_batch_model_for_callback(uniform, fallback_probe, 2) === gram_a

        p_uniform = ODEParams(n_sats=2, args=cfg_gram)
        push!(p_uniform.shared_buffers.density_models, gram_a, gram_a)
        @test CB._density_batch_model_for_callback(p_uniform, 2) === gram_a

        @test CB._gram_batch_elapsed_time(5.5, 3) === 5.5
        @test CB._gram_batch_elapsed_time([1.0, 2.0, 3.0], 2) === 2.0
        @test CB._gram_batch_elapsed_time(Float32[4.0, 8.0], 2) === 8.0
    end

    @testset "model_selection: isolated pool gating" begin
        empty_models = SimulationModel.AbstractDensityModel[]
        gram_fallback = make_fake_gram()

        withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => "off") do
            @test CB._gram_isolated_pool_enabled(8) == false
            @test CB._gram_isolated_pool_batch_model_for_callback(empty_models, gram_fallback, 8) === nothing
        end

        withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => "on") do
            @test CB._gram_isolated_pool_enabled(1) == true
            @test CB._gram_isolated_pool_enabled(0) == false
            # Enabled + empty per-sat models + GRAM fallback -> pool model.
            @test CB._gram_isolated_pool_batch_model_for_callback(empty_models, gram_fallback, 4) === gram_fallback
            # Non-GRAM fallback is never pooled.
            @test CB._gram_isolated_pool_batch_model_for_callback(empty_models, fallback_probe, 4) === nothing
            # Per-sat models present -> pool path declined.
            occupied = SimulationModel.AbstractDensityModel[gram_a]
            @test CB._gram_isolated_pool_batch_model_for_callback(occupied, gram_fallback, 4) === nothing

            env = CB._snapshot_callback_env_config()
            @test env.gram_isolated_pool_mode == :on
            @test CB._gram_isolated_pool_enabled(env, 2) == true
            @test CB._gram_isolated_pool_batch_model_for_callback(env, empty_models, gram_fallback, 2) === gram_fallback
            @test CB._gram_isolated_pool_batch_model_for_callback(env, empty_models, fallback_probe, 2) === nothing
            @test CB._gram_isolated_pool_batch_model_for_callback(env, SimulationModel.AbstractDensityModel[gram_a], gram_fallback, 2) === nothing

            p = ODEParams(n_sats=2, args=cfg_gram)
            @test CB._gram_isolated_pool_batch_model_for_callback(p, 2) === cfg_gram.environment_model.density_model
        end

        withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => "auto",
                "SPACEAGORA_GRAM_ISOLATED_POOL_THRESHOLD" => "4") do
            expected = Threads.nthreads() > 1
            @test CB._gram_isolated_pool_enabled(4) == expected
            @test CB._gram_isolated_pool_enabled(3) == false
            env = CB._snapshot_callback_env_config()
            @test CB._gram_isolated_pool_enabled(env, 4) == expected
            @test CB._gram_isolated_pool_enabled(env, 3) == false
        end
    end

    @testset "model_selection: per-sat cache slots" begin
        caches = Vector{Union{Nothing, CB.GramTrackCache}}(nothing, 2)
        c1 = CB._gram_density_cache_for_sat!(caches, 1)
        @test c1 isa CB.GramTrackCache
        @test caches[1] === c1
        @test CB._gram_density_cache_for_sat!(caches, 1) === c1
        # Out-of-range index returns a detached cache and stores nothing.
        c9 = CB._gram_density_cache_for_sat!(caches, 9)
        @test c9 isa CB.GramTrackCache
        @test length(caches) == 2

        vcaches = Vector{Union{Nothing, CB.VacuumPredictedGRAMCache}}(nothing, 2)
        v1 = CB._vacuum_gram_cache_for_sat!(vcaches, 1)
        @test v1 isa CB.VacuumPredictedGRAMCache
        @test v1.valid == false
        @test vcaches[1] === v1
        @test CB._vacuum_gram_cache_for_sat!(vcaches, 1) === v1
        v9 = CB._vacuum_gram_cache_for_sat!(vcaches, 9)
        @test v9 isa CB.VacuumPredictedGRAMCache
        @test length(vcaches) == 2
    end

    @testset "model_selection: isolated pool lifecycle" begin
        p = ODEParams(n_sats=2, args=cfg_gram)
        template = make_fake_gram(rho=3.0e-9, T=222.0)

        models0, locks0 = CB._ensure_gram_isolated_pool!(p, template, 0)
        @test isempty(models0) && isempty(locks0)

        models, locks = CB._ensure_gram_isolated_pool!(p, template, 2)
        @test length(models) == 2 && length(locks) == 2
        @test models === p.shared_buffers.gram_isolated_pool_models
        @test all(m -> m.core isa FakeGramCore && m.core.rho == 3.0e-9, models)
        # Members are deep copies, not aliases of the template.
        @test all(m -> m !== template, models)

        first_model = models[1]
        models_again, _ = CB._ensure_gram_isolated_pool!(p, template, 2)
        @test models_again[1] === first_model  # no rebuild when sizes match

        models_small, locks_small = CB._ensure_gram_isolated_pool!(p, template, 1)
        @test length(models_small) == 1 && length(locks_small) == 1
    end

    @testset "model_selection: isolated pool density state" begin
        p_kep = ODEParams(n_sats=1, args=cfg_gram)
        model = make_fake_gram(rho=4.0e-9, T=233.0)
        lk = ReentrantLock()

        # Above 2000 km: hard vacuum.
        rho, T, wind = CB._gram_isolated_pool_density_state(model, 2.5e6, 0.1, 0.2, 10.0, true, p_kep, lk)
        @test rho == 0.0
        @test T == EARTH.T_ref
        @test wind == SVector{3, Float64}(0.0, 0.0, 0.0)

        # Below EI (drag state): defer to the (fake) GRAM core.
        rho, T, wind = CB._gram_isolated_pool_density_state(model, 100.0e3, 0.1, 0.2, 10.0, true, p_kep, lk)
        @test rho == 4.0e-9
        @test T == 233.0
        @test wind == SVector{3, Float64}(1.0, 2.0, 3.0)

        # Keplerian run above EI still routes to the core.
        rho, _, _ = CB._gram_isolated_pool_density_state(model, 500.0e3, 0.0, 0.0, 0.0, true, p_kep, lk)
        @test rho == 4.0e-9

        # Non-keplerian run above EI uses the analytic polyfit shortcut.
        p_nk = ODEParams(n_sats=1, args=cfg_nonkeplerian)
        rho_nk, T_nk, _ = CB._gram_isolated_pool_density_state(model, 500.0e3, 0.0, 0.0, 0.0, true, p_nk, lk)
        rho_ref, T_ref_poly, _ = EM.density_polyfit(500.0e3, p_nk)
        @test rho_nk == rho_ref
        @test T_nk == T_ref_poly
        @test rho_nk > 0.0
    end

    @testset "model_selection: isolated pool batch eval" begin
        p = ODEParams(n_sats=1, args=cfg_gram)
        model = make_fake_gram(rho=6.0e-9, T=240.0)
        n = 4
        hs = [2.5e6, 2.6e6, 80.0e3, 90.0e3]
        lats = zeros(n)
        lons = zeros(n)
        rhos = zeros(n)
        Ts = zeros(n)
        winds = [SVector{3, Float64}(9.0, 9.0, 9.0) for _ in 1:n]

        # Non-GRAM model dispatches to the generic method: always false.
        withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => "on") do
            @test CB._gram_isolated_pool_batch_eval!(rhos, Ts, winds, fallback_probe, hs, lats, lons, 0.0, true, p) == false
        end

        # Pool disabled: false before any work.
        withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => "off") do
            @test CB._gram_isolated_pool_batch_eval!(rhos, Ts, winds, model, hs, lats, lons, 0.0, true, p) == false
        end

        withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => "on") do
            # Length-mismatch guards, one per buffer.
            @test CB._gram_isolated_pool_batch_eval!(zeros(n - 1), Ts, winds, model, hs, lats, lons, 0.0, true, p) == false
            @test CB._gram_isolated_pool_batch_eval!(rhos, zeros(n - 1), winds, model, hs, lats, lons, 0.0, true, p) == false
            @test CB._gram_isolated_pool_batch_eval!(rhos, Ts, winds[1:n-1], model, hs, lats, lons, 0.0, true, p) == false
            @test CB._gram_isolated_pool_batch_eval!(rhos, Ts, winds, model, hs, zeros(n - 1), lons, 0.0, true, p) == false
            @test CB._gram_isolated_pool_batch_eval!(rhos, Ts, winds, model, hs, lats, zeros(n - 1), 0.0, true, p) == false
            @test CB._gram_isolated_pool_batch_eval!(rhos, Ts, winds, model, hs, lats, lons, zeros(n - 1), true, p) == false
            # Single-worker allotment declines the pool.
            @test CB._gram_isolated_pool_batch_eval!(rhos, Ts, winds, model, hs, lats, lons, 0.0, true, p; allotment_hint=1) == false

            if Threads.nthreads() > 1
                ok = CB._gram_isolated_pool_batch_eval!(rhos, Ts, winds, model, hs, lats, lons, 0.0, true, p; allotment_hint=2)
                @test ok == true
                @test rhos[1] == 0.0 && rhos[2] == 0.0            # vacuum knots
                @test rhos[3] == 6.0e-9 && rhos[4] == 6.0e-9      # fake-core knots
                @test Ts[1] == EARTH.T_ref && Ts[3] == 240.0
                @test winds[1] == SVector{3, Float64}(0.0, 0.0, 0.0)
                @test winds[3] == SVector{3, Float64}(1.0, 2.0, 3.0)
                @test length(p.shared_buffers.gram_isolated_pool_models) >= 2

                # Vector elapsed-time variant.
                fill!(rhos, -1.0)
                ok_vec = CB._gram_isolated_pool_batch_eval!(rhos, Ts, winds, model, hs, lats, lons, collect(1.0:Float64(n)), true, p; allotment_hint=2)
                @test ok_vec == true
                @test rhos[3] == 6.0e-9
            end
        end
    end

    @testset "vacuum_predicted_gram: env config parsing" begin
        withenv("SPACEAGORA_VACUUM_GRAM_CACHE" => nothing,
                "SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS" => nothing,
                "SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S" => nothing,
                "SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M" => nothing) do
            @test CB._vacuum_gram_cache_enabled() == false
            @test CB._vacuum_gram_cache_npoints() == 20
            @test CB._vacuum_gram_cache_horizon_s() == 600.0
            @test CB._vacuum_gram_cache_deviation_m() == 5000.0
        end
        withenv("SPACEAGORA_VACUUM_GRAM_CACHE" => "1") do
            @test CB._vacuum_gram_cache_enabled() == true
        end
        withenv("SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS" => "7") do
            @test CB._vacuum_gram_cache_npoints() == 7
        end
        withenv("SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS" => "2") do
            @test CB._vacuum_gram_cache_npoints() == 4  # floor at 4 knots
        end
        withenv("SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS" => "abc") do
            @test_throws ArgumentError CB._vacuum_gram_cache_npoints()
        end
        withenv("SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S" => "-5.0") do
            @test CB._vacuum_gram_cache_horizon_s() == 10.0  # floor at 10 s
        end
        withenv("SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S" => "1200") do
            @test CB._vacuum_gram_cache_horizon_s() == 1200.0
        end
        withenv("SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M" => "1.0") do
            @test CB._vacuum_gram_cache_deviation_m() == 100.0  # floor at 100 m
        end
        withenv("SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M" => "9000") do
            @test CB._vacuum_gram_cache_deviation_m() == 9000.0
        end
    end

    @testset "vacuum_predicted_gram: spline and interp primitives" begin
        # Natural cubic spline: n <= 2 leaves all second derivatives at zero.
        Ms = Float64[]
        CB._natural_cubic_spline_build!(Ms, [1.0, 2.0], 1.0)
        @test Ms == [0.0, 0.0]
        CB._natural_cubic_spline_build!(Ms, [1.0], 1.0)
        @test Ms == [0.0]

        # Linear data has zero curvature; the spline reproduces it exactly.
        h = 2.0
        ys_lin = [3.0 + 0.5 * (i - 1) * h for i in 1:6]
        Ms_lin = Float64[]
        CB._natural_cubic_spline_build!(Ms_lin, ys_lin, h)
        @test maximum(abs, Ms_lin) < 1e-12
        for t in (0.0, 0.7, 4.4, 9.9)
            @test CB._eval_natural_cubic_spline(t, 0.0, h, ys_lin, Ms_lin) ≈ 3.0 + 0.5 * t atol=1e-10
        end

        # Generic data: the spline interpolates the knot values exactly.
        ys = [0.0, 1.0, -0.5, 2.0, 0.25]
        Ms_g = Float64[]
        CB._natural_cubic_spline_build!(Ms_g, ys, 1.0)
        @test Ms_g[1] == 0.0 && Ms_g[end] == 0.0  # natural boundary
        for (i, y) in enumerate(ys)
            @test CB._eval_natural_cubic_spline(Float64(i - 1), 0.0, 1.0, ys, Ms_g) ≈ y atol=1e-10
        end
        # Clamped extrapolation outside the knot range stays finite.
        @test isfinite(CB._eval_natural_cubic_spline(-0.5, 0.0, 1.0, ys, Ms_g))
        @test isfinite(CB._eval_natural_cubic_spline(99.0, 0.0, 1.0, ys, Ms_g))

        # Piecewise-linear interpolators over a hand-built cache.
        cache = CB.VacuumPredictedGRAMCache()
        cache.t0 = 10.0
        cache.t1 = 12.0
        cache.h = 1.0
        append!(cache.vac_alts, [100.0, 200.0, 400.0])
        append!(cache.vac_positions, [SVector{3, Float64}(1.0, 0.0, 0.0),
                                      SVector{3, Float64}(3.0, 2.0, 0.0),
                                      SVector{3, Float64}(5.0, 4.0, 2.0)])
        append!(cache.winds, [SVector{3, Float64}(0.0, 0.0, 0.0),
                              SVector{3, Float64}(2.0, -2.0, 4.0),
                              SVector{3, Float64}(4.0, -4.0, 8.0)])
        @test CB._interp_vacuum_alt(cache, 10.5) ≈ 150.0
        @test CB._interp_vacuum_alt(cache, 12.0) ≈ 400.0
        @test CB._interp_vacuum_position(cache, 10.5) ≈ SVector{3, Float64}(2.0, 1.0, 0.0)
        @test CB._interp_vacuum_wind(cache, 11.5) ≈ SVector{3, Float64}(3.0, -3.0, 6.0)
        # Clamped beyond the ends: linear continuation from the last segment.
        @test isfinite(CB._interp_vacuum_alt(cache, 20.0))
    end

    @testset "vacuum_predicted_gram: J2 accel and RK4 step" begin
        r = EARTH.Rp_e + 400e3
        pos = SVector{3, Float64}(r, 0.0, 0.0)
        acc = CB._vacuum_j2_accel(pos, EARTH)
        # Manual evaluation of the same expression (equatorial: z = 0).
        # J2 is normalized to the equatorial radius, matching the source.
        j2_scale = 1.5 * EARTH.J2 * EARTH.μ * EARTH.Rp_e^2 / r^4
        expected_x = -EARTH.μ / r^2 + j2_scale * (-1.0)
        @test acc[1] ≈ expected_x rtol=1e-12
        @test abs(acc[2]) < 1e-15 && abs(acc[3]) < 1e-15

        # Polar position: J2 pushes outward along z relative to two-body.
        pos_pole = SVector{3, Float64}(0.0, 0.0, r)
        acc_pole = CB._vacuum_j2_accel(pos_pole, EARTH)
        @test acc_pole[3] ≈ -EARTH.μ / r^2 + j2_scale * 2.0 rtol=1e-12

        vel = SVector{3, Float64}(0.0, sqrt(EARTH.μ / r), 0.0)
        new_pos, new_vel = CB._vacuum_rk4_step(pos, vel, EARTH, 1.0)
        @test norm(new_pos - (pos + vel * 1.0)) < 10.0  # curvature over 1 s is < 10 m
        # Specific orbital energy is conserved to high accuracy over one step.
        e0 = 0.5 * dot(vel, vel) - EARTH.μ / norm(pos)
        e1 = 0.5 * dot(new_vel, new_vel) - EARTH.μ / norm(new_pos)
        @test isapprox(e0, e1; rtol=1e-9)
    end

    @testset "vacuum_predicted_gram: cache build and query" begin
        p = ODEParams(n_sats=1, args=cfg_nongram)
        model = fallback_probe
        pos, vel = circular_state(300e3)
        t = 0.0
        n_pts = 8
        horizon = 300.0
        dev_m = 5000.0

        cache = CB.VacuumPredictedGRAMCache()
        # Degenerate knot count leaves the cache invalid.
        CB._build_vacuum_gram_cache!(cache, model, p, pos, vel, t, 1, horizon)
        @test cache.valid == false

        CB._build_vacuum_gram_cache!(cache, model, p, pos, vel, t, n_pts, horizon)
        @test cache.valid == true
        @test cache.t0 == t
        @test cache.t1 ≈ t + horizon
        @test cache.h ≈ horizon / (n_pts - 1)
        @test length(cache.log_rhos) == n_pts
        @test cache.vac_positions[1] == pos
        # First knot reproduces the direct density query at the initial point.
        l_pi = CB._planet_lpi_at(p, t)
        alt0, lat0, lon0 = CB.rtolatlong(SVector{3, Float64}(l_pi * pos), EARTH)
        rho_direct, T_direct, wind_direct = EM.getDensity(model, alt0, lat0, lon0, t, true, p)
        @test exp(cache.log_rhos[1]) ≈ rho_direct rtol=1e-12
        @test cache.Ts[1] == T_direct
        @test cache.winds[1] == wind_direct
        @test cache.vac_alts[1] ≈ alt0

        # Query on the predicted trajectory: served from the spline.
        t_q = t + 0.5 * cache.h
        pos_q = CB._interp_vacuum_position(cache, t_q)
        rho_q, T_q, wind_q = CB._query_vacuum_gram_cache!(cache, model, p, pos_q, vel, 0.0, t_q, n_pts, horizon, dev_m)
        log_expect = CB._eval_natural_cubic_spline(t_q, cache.t0, cache.h, cache.log_rhos, cache.Ms_rho)
        @test rho_q ≈ exp(log_expect) rtol=1e-12
        @test T_q ≈ CB._eval_natural_cubic_spline(t_q, cache.t0, cache.h, cache.Ts, cache.Ms_T) rtol=1e-12
        @test wind_q ≈ CB._interp_vacuum_wind(cache, t_q)
        @test cache.t0 == t  # hit did not rebuild

        # Deviation beyond threshold forces a rebuild anchored at the query.
        pos_dev = pos_q + SVector{3, Float64}(2.0 * dev_m, 0.0, 0.0)
        rho_rb, T_rb, _ = CB._query_vacuum_gram_cache!(cache, model, p, pos_dev, vel, 0.0, t_q, n_pts, horizon, dev_m)
        @test cache.t0 == t_q  # rebuilt from the query time
        @test rho_rb ≈ exp(cache.log_rhos[1]) rtol=1e-12
        @test T_rb == cache.Ts[1]

        # Query outside the time window also rebuilds.
        t_late = cache.t1 + 100.0
        CB._query_vacuum_gram_cache!(cache, model, p, pos, vel, 0.0, t_late, n_pts, horizon, dev_m)
        @test cache.t0 == t_late

        # Invalid cache + degenerate knot count falls through to a direct query.
        bad = CB.VacuumPredictedGRAMCache()
        rho_fb, T_fb, wind_fb = CB._query_vacuum_gram_cache!(bad, model, p, pos, vel, 0.0, t, 1, horizon, dev_m)
        @test bad.valid == false
        @test rho_fb ≈ rho_direct rtol=1e-12
        @test T_fb == T_direct
        @test wind_fb == wind_direct
    end

    @testset "targeting: env knobs" begin
        withenv("SPACEAGORA_GRAM_TRACK_CACHE_PERIAPSIS_SPLIT" => nothing) do
            @test CB._gram_track_cache_periapsis_split_enabled() == true
        end
        withenv("SPACEAGORA_GRAM_TRACK_CACHE_PERIAPSIS_SPLIT" => "0") do
            @test CB._gram_track_cache_periapsis_split_enabled() == false
        end
        withenv("SPACEAGORA_GRAM_TRACK_CACHE_MAX_NPOS" => nothing) do
            @test CB._gram_track_cache_max_npos() == 512
        end
        withenv("SPACEAGORA_GRAM_TRACK_CACHE_MAX_NPOS" => "1") do
            @test CB._gram_track_cache_max_npos() == 2  # floor at 2
        end
        withenv("SPACEAGORA_GRAM_TRACK_CACHE_MAX_NPOS" => "abc") do
            @test_throws ArgumentError CB._gram_track_cache_max_npos()
        end
    end

    @testset "targeting: geometry helpers" begin
        @test CB._angle_delta_rad(0.1, 0.1 + 2pi + 0.05) ≈ 0.05 atol=1e-12
        @test CB._angle_delta_rad(0.2, 0.2 - 0.3) ≈ -0.3 atol=1e-12
        @test CB._angle_delta_rad(pi - 0.1, -pi + 0.1) ≈ 0.2 atol=1e-12

        # Coincident endpoints degrade to the 1 m floor.
        len, radius = CB._gram_expected_track_length_m(EARTH, 100e3, 0.3, 0.4, 100e3, 0.3, 0.4)
        @test len == 1.0
        @test radius ≈ EARTH.Rp_m + 100e3

        # Pure vertical displacement.
        len_v, _ = CB._gram_expected_track_length_m(EARTH, 100e3, 0.3, 0.4, 101e3, 0.3, 0.4)
        @test len_v ≈ 1000.0 rtol=1e-6

        # Pure horizontal displacement follows the great-circle arc.
        dlat = 0.01
        len_h, radius_h = CB._gram_expected_track_length_m(EARTH, 100e3, 0.0, 0.4, 100e3, dlat, 0.4)
        @test len_h ≈ radius_h * dlat rtol=1e-6

        # Spacing is half the tightest tolerance, floored at 1 m.
        @test CB._gram_track_cache_target_spacing_m(500.0, deg2rad(1.0), 6.5e6) ≈ 250.0
        ang_tight = CB._gram_track_cache_target_spacing_m(1e6, 1e-12, 6.5e6)
        @test ang_tight == 1.0  # angular floor drives spacing to the minimum
    end

    @testset "targeting: Kepler solvers" begin
        # e = 0: E == M exactly.
        @test CB._solve_kepler_elliptic(1.234, 0.0) ≈ 1.234 atol=1e-12
        for (M, e) in ((1.0, 0.5), (0.3, 0.9), (5.8, 0.95))
            E = CB._solve_kepler_elliptic(M, e)
            @test E - e * sin(E) ≈ mod(M, 2pi) atol=1e-9
        end
        for ν in (0.0, 0.7, Float64(pi), 4.5), e in (0.0, 0.3, 0.8)
            E = CB._true_to_eccentric_anomaly(ν, e)
            @test CB._eccentric_to_true_anomaly(E, e) ≈ mod(ν, 2pi) atol=1e-9
        end
    end

    @testset "targeting: conic endpoint targeting" begin
        pos_c, vel_c = circular_state(500e3)
        n_c = sqrt(EARTH.μ / norm(pos_c)^3)
        period = 2pi / n_c

        target = CB._gram_kepler_target(pos_c, vel_c, EARTH, period / 4; include_j2=false)
        @test target !== nothing
        alt_t, lat_t, lon_t = target
        @test isapprox(alt_t, 500e3; atol=40e3)  # circular: altitude preserved (geodetic wobble)
        @test -pi / 2 <= lat_t <= pi / 2
        @test -pi <= lon_t <= pi

        # Degenerate / invalid inputs.
        @test CB._gram_kepler_target(pos_c, vel_c, EARTH, 0.0) === nothing
        @test CB._gram_kepler_target(pos_c, vel_c, EARTH, NaN) === nothing
        pos_h, vel_h = hyperbolic_state(500e3)
        @test CB._gram_kepler_target(pos_h, vel_h, EARTH, 100.0) === nothing
        # Unusable planet object trips the catch fallback.
        @test CB._gram_kepler_target(pos_c, vel_c, nothing, 100.0) === nothing

        # J2 secular rates change the endpoint.
        pos_e, vel_e = elliptical_apoapsis_state(120e3, 800e3)
        t_no_j2 = CB._gram_kepler_target(pos_e, vel_e, EARTH, 1500.0; include_j2=false)
        t_j2 = CB._gram_kepler_target(pos_e, vel_e, EARTH, 1500.0; include_j2=true)
        @test t_no_j2 !== nothing && t_j2 !== nothing
        @test collect(t_no_j2) != collect(t_j2)

        # Periapsis targeting from apoapsis: half the orbital period.
        peri = CB._gram_periapsis_target(pos_e, vel_e, EARTH; include_j2=false)
        @test peri !== nothing
        dt_peri, alt_peri, _, _ = peri
        a_e = 0.5 * (2 * EARTH.Rp_e + 120e3 + 800e3)
        period_e = 2pi / sqrt(EARTH.μ / a_e^3)
        @test dt_peri ≈ period_e / 2 rtol=1e-6
        @test alt_peri < 800e3  # descended from apoapsis
        @test isapprox(alt_peri, 120e3; atol=40e3)
        @test CB._gram_periapsis_target(pos_h, vel_h, EARTH) === nothing
        @test CB._gram_periapsis_target(pos_e, vel_e, nothing) === nothing
        @test CB._gram_descending_periapsis_target(pos_e, vel_e, EARTH; include_j2=false) == peri

        # Orbit-period targeting.
        orbit = CB._gram_orbit_period_target(pos_c, vel_c, EARTH; include_j2=false)
        @test orbit !== nothing
        dt_orbit, alt_end, _, _ = orbit
        @test dt_orbit ≈ period rtol=1e-6
        @test isapprox(alt_end, 500e3; atol=40e3)
        @test CB._gram_orbit_period_target(pos_h, vel_h, EARTH) === nothing
        @test CB._gram_orbit_period_target(pos_c, vel_c, nothing) === nothing

        # Degenerate micro-orbit (r = 0.5 m): the orbital period and the time
        # to periapsis drop below the 1 µs floor.
        pos_tiny = SVector{3, Float64}(0.5, 0.0, 0.0)
        vel_tiny = SVector{3, Float64}(0.0, sqrt(EARTH.μ / 0.5), 0.0)
        @test CB._gram_periapsis_target(pos_tiny, vel_tiny, EARTH) === nothing
        @test CB._gram_orbit_period_target(pos_tiny, vel_tiny, EARTH) === nothing
        # Sub-denormal semi-major axis: a³ underflows and the mean motion is
        # non-finite; every targeting routine must decline.
        pos_eps = SVector{3, Float64}(1.0e-150, 0.0, 0.0)
        vel_eps = SVector{3, Float64}(0.0, sqrt(EARTH.μ / 1.0e-150), 0.0)
        @test CB._gram_kepler_target(pos_eps, vel_eps, EARTH, 100.0) === nothing
        @test CB._gram_periapsis_target(pos_eps, vel_eps, EARTH) === nothing
        @test CB._gram_orbit_period_target(pos_eps, vel_eps, EARTH) === nothing

        # Linear extrapolation target.
        alt_lin, lat_lin, lon_lin = CB._gram_linear_target(pos_c, vel_c, EARTH, 10.0)
        @test isfinite(alt_lin) && isfinite(lat_lin) && isfinite(lon_lin)
        @test isapprox(alt_lin, 500e3; atol=50e3)
    end

    @testset "targeting: entry ballistic model" begin
        p = ODEParams(n_sats=1, args=cfg_gram)

        # Reference area: sum of positive link areas.
        @test CB._gram_entry_reference_area_m2(p, 1) == 18.0
        # All-zero link areas: final floor of 1 m².
        sc_zero = make_probe_spacecraft(root_area=0.0, panel_area=0.0)
        p_zero = ODEParams(n_sats=1, args=probe_config(density_model=make_fake_gram(), spacecraft=[sc_zero]))
        @test CB._gram_entry_reference_area_m2(p_zero, 1) == 1.0
        # Links empty of area but positive root area: root fallback.
        root_pos = Link(root=true, m=500.0, ref_area=7.5)
        panel_zero = Link(root=false, m=30.0, ref_area=0.0, r=MVector{3, Float64}(0.0, 1.2, 0.0))
        sc_root = SpacecraftModel(Joint[], [panel_zero], root_pos, true, 530.0, 0.0, root_pos.inertia, 0, 0,
                                  InitialCondition(ra=EARTH.Rp_e + 500e3, rp=EARTH.Rp_e + 450e3, i=35.0, ω=40.0, Ω=10.0, ν=175.0), 1)
        p_root = ODEParams(n_sats=1, args=probe_config(density_model=make_fake_gram(), spacecraft=[sc_root]))
        @test CB._gram_entry_reference_area_m2(p_root, 1) == 7.5
        # Out-of-range satellite index: guarded fallback.
        @test CB._gram_entry_reference_area_m2(p, 99) == 1.0

        # Entry mass resolution.
        @test CB._gram_entry_mass_kg(p, 1, 250.0) == 250.0
        @test CB._gram_entry_mass_kg(p, 1, NaN) == 530.0   # dry + prop mass
        @test CB._gram_entry_mass_kg(p, 1, -5.0) == 530.0
        @test CB._gram_entry_mass_kg(p, 99, NaN) == 100.0  # guarded fallback

        # Exponential reference density.
        @test CB._gram_entry_reference_density(EARTH, 0.0) ≈ EARTH.ρ_ref
        @test CB._gram_entry_reference_density(EARTH, EARTH.H) ≈ EARTH.ρ_ref * exp(-1.0) rtol=1e-12
        @test CB._gram_entry_reference_density(EARTH, -1.0e9) == 0.0  # overflow guard

        # Allen-Eggers guards.
        pos_en = SVector{3, Float64}(EARTH.Rp_e + 80e3, 0.0, 0.0)
        vel_en = SVector{3, Float64}(0.0, 7400.0, -800.0)
        @test CB._gram_entry_target_allen_eggers(pos_en, vel_en, EARTH, -1.0, 500.0, 10.0) === nothing
        @test CB._gram_entry_target_allen_eggers(pos_en, vel_en, EARTH, 300.0, NaN, 10.0) === nothing
        @test CB._gram_entry_target_allen_eggers(pos_en, vel_en, EARTH, 300.0, 500.0, 0.0) === nothing
        earth_no_atm = Earth(ρ_ref=0.0)
        @test CB._gram_entry_target_allen_eggers(pos_en, vel_en, earth_no_atm, 300.0, 500.0, 10.0) === nothing
        # Co-rotating state has ~zero planet-relative speed.
        vel_corot = SVector{3, Float64}(cross(EARTH.ω, pos_en))
        @test CB._gram_entry_target_allen_eggers(pos_en, vel_corot, EARTH, 300.0, 500.0, 10.0) === nothing

        # Nominal steep entry descends.
        result = CB._gram_entry_target_allen_eggers(pos_en, vel_en, EARTH, 120.0, 500.0, 10.0)
        @test result !== nothing
        h_end, lat_end, lon_end = result
        @test isfinite(h_end) && isfinite(lat_end) && isfinite(lon_end)
        @test h_end < 80e3
        @test -pi < lon_end <= pi

        # Shallow high-altitude arcs exercise the heading/longitude wrap branches.
        # Eastbound across the antimeridian: lon wraps +π -> -π.
        pos_e1, vel_e1 = ned_probe_state(0.0, 179.0, 200e3, 7800.0, pi / 2)
        res_e1 = CB._gram_entry_target_allen_eggers(pos_e1, vel_e1, EARTH, 600.0, 4000.0, 1.0)
        @test res_e1 !== nothing && all(isfinite, collect(res_e1))
        @test -pi < res_e1[3] <= pi
        # Westbound across the antimeridian: lon wraps -π -> +π.
        pos_w1, vel_w1 = ned_probe_state(0.0, -179.0, 200e3, 7800.0, -pi / 2)
        res_w1 = CB._gram_entry_target_allen_eggers(pos_w1, vel_w1, EARTH, 600.0, 4000.0, 1.0)
        @test res_w1 !== nothing && all(isfinite, collect(res_w1))
        @test -pi < res_w1[3] <= pi
        # Near-polar arcs with heading close to ±π: χ wraps in both directions.
        pos_p1, vel_p1 = ned_probe_state(89.85, 0.0, 200e3, 7800.0, 3.1)
        res_p1 = CB._gram_entry_target_allen_eggers(pos_p1, vel_p1, EARTH, 6000.0, 4000.0, 1.0)
        @test res_p1 !== nothing && all(isfinite, collect(res_p1))
        pos_p2, vel_p2 = ned_probe_state(89.85, 0.0, 200e3, 7800.0, -3.1)
        res_p2 = CB._gram_entry_target_allen_eggers(pos_p2, vel_p2, EARTH, 6000.0, 4000.0, 1.0)
        @test res_p2 !== nothing && all(isfinite, collect(res_p2))
    end

    @testset "refresh: fill from trajectory" begin
        cache = CB.GramTrackCache()
        @test_throws ArgumentError CB._gram_track_cache_fill_from_trajectory!(cache, Any[nothing])

        traj = fake_generate_trajectory(nothing;
            initial_height=100.0, initial_latitude=5.0, initial_longitude=190.0,
            initial_elapsed_time=50.0, delta_height=-2.0, delta_latitude=0.1,
            delta_longitude=0.1, delta_elapsed_time=3.0, n_points=3,
            update_initial_perturbations=true)
        CB._gram_track_cache_fill_from_trajectory!(cache, traj)
        @test cache.times == [50.0, 53.0, 56.0]
        @test cache.alts ≈ [100.0e3, 98.0e3, 96.0e3]
        @test cache.lats ≈ deg2rad.([5.1, 5.2, 5.3])
        # Longitudes normalized into (-π, π] with no trig.
        @test cache.lons[1] ≈ deg2rad(190.0) - 2pi
        @test cache.lons[2] ≈ deg2rad(-190.0) + 2pi
        @test cache.lons[3] ≈ deg2rad(30.0)
        @test cache.rhos ≈ [1.0e-8, 2.0e-8, 3.0e-8]
        @test cache.Ts ≈ [201.0, 202.0, 203.0]
        @test cache.winds[2] == SVector{3, Float64}(2.0, -2.0, 1.0)
    end

    @testset "refresh: kepler-or-linear target" begin
        pos_c, vel_c = circular_state(500e3)
        kt = CB._gram_kepler_or_linear_target(pos_c, vel_c, EARTH, 600.0, false)
        @test kt == CB._gram_kepler_target(pos_c, vel_c, EARTH, 600.0; include_j2=false)
        pos_h, vel_h = hyperbolic_state(500e3)
        lt = CB._gram_kepler_or_linear_target(pos_h, vel_h, EARTH, 600.0, false)
        @test lt == CB._gram_linear_target(pos_h, vel_h, EARTH, 600.0)
    end

    @testset "refresh: track cache segments" begin
        model = fallback_probe
        p_time = ODEParams(n_sats=1, args=cfg_nongram)
        horizon = 60.0
        n_points = 8

        run_refresh(p, pos, vel, alt, t; dm=model, seg_end=NaN, kwargs...) = withenv(
            "SPACEAGORA_GRAM_TRACK_CACHE_MAX_NPOS" => "64",
            "SPACEAGORA_GRAM_ISOLATED_POOL" => "off"
        ) do
            cache = CB.GramTrackCache()
            out = CB._gram_track_cache_refresh!(
                cache, dm, p, pos, vel, alt, 0.0, 0.0, t, horizon, n_points,
                500.0, deg2rad(1.0), 0.0, seg_end, get(kwargs, :use_j2, true), 1, NaN
            )
            return cache, out
        end

        # 1. Time mission, out of the atmosphere band, no solver endpoint:
        #    falls back to one full orbital period.
        pos_c, vel_c = circular_state(500e3)
        period_c = 2pi / sqrt(EARTH.μ / norm(pos_c)^3)
        cache, out = run_refresh(p_time, pos_c, vel_c, 500e3, 0.0)
        @test cache.valid == true
        @test cache.t0 == 0.0
        @test isapprox(cache.t1, period_c; rtol=0.05)  # J2 perturbs the rates slightly
        @test length(cache.times) <= 64  # MAX_NPOS cap enforced
        @test all(>(0.0), cache.rhos)
        @test out[1] == cache.rhos[1] && out[2] == cache.Ts[1] && out[3] == cache.winds[1]
        @test cache.index_hint == 1

        # 2. Time mission with a finite solver endpoint: segment ends there.
        cache2, _ = run_refresh(p_time, pos_c, vel_c, 500e3, 10.0; seg_end=60.0)
        @test cache2.valid == true
        @test cache2.t0 == 10.0
        @test cache2.t1 ≈ 60.0

        # 3. Atmosphere band + elliptic orbit: segment runs to periapsis.
        pos_e, vel_e = elliptical_apoapsis_state(60e3, 300e3)
        peri = CB._gram_periapsis_target(pos_e, vel_e, EARTH; include_j2=true)
        @test peri !== nothing
        cache3, _ = run_refresh(p_time, pos_e, vel_e, 100e3, 0.0)
        @test cache3.valid == true
        @test cache3.t1 ≈ peri[1] rtol=1e-6

        # 4. Atmosphere band + hyperbolic entry + Allen-Eggers targeting.
        pos_h, vel_h = hyperbolic_state(100e3)
        cache4, _ = withenv("SPACEAGORA_GRAM_ENTRY_TARGET_MODE" => "allen_eggers") do
            run_refresh(p_time, pos_h, vel_h, 100e3, 0.0; seg_end=45.0)
        end
        @test cache4.valid == true
        @test cache4.t1 ≈ 45.0  # dt_entry = segment_end_t - t

        # 5. Entry model off, solver endpoint available.
        cache5, _ = withenv("SPACEAGORA_GRAM_ENTRY_TARGET_MODE" => "off") do
            run_refresh(p_time, pos_h, vel_h, 100e3, 5.0; seg_end=45.0)
        end
        @test cache5.valid == true
        @test cache5.t1 ≈ 45.0

        # 6. Entry model off, no endpoint: base horizon retained.
        cache6, _ = withenv("SPACEAGORA_GRAM_ENTRY_TARGET_MODE" => "off") do
            run_refresh(p_time, pos_h, vel_h, 100e3, 0.0)
        end
        @test cache6.valid == true
        @test cache6.t1 ≈ horizon

        # 7. Orbit mission out of the band: one full orbital period.
        cfg_orbit = probe_config(density_model=fallback_probe, mission_type=MissionOrbits)
        p_orbit = ODEParams(n_sats=1, args=cfg_orbit)
        cache7, _ = run_refresh(p_orbit, pos_c, vel_c, 500e3, 0.0)
        @test cache7.valid == true
        @test isapprox(cache7.t1, period_c; rtol=0.05)

        # 8. Runtime stats accounting.
        withenv("SPACEAGORA_GRAM_PROFILE" => "1",
                "SPACEAGORA_GRAM_TRACK_CACHE_MAX_NPOS" => "64",
                "SPACEAGORA_GRAM_ISOLATED_POOL" => "off") do
            CB._gram_runtime_stats_reset!()
            cache_s = CB.GramTrackCache()
            CB._gram_track_cache_refresh!(cache_s, model, p_time, pos_c, vel_c, 500e3, 0.0, 0.0, 0.0, horizon, n_points)
            stats = CB._gram_runtime_stats_snapshot()
            @test stats.refresh_calls == 1
            @test stats.refresh_points_total == length(cache_s.times)
            @test stats.refresh_elapsed_s > 0.0
            CB._gram_runtime_stats_reset!()
        end
    end

    @testset "refresh: isolated-pool batch path" begin
        if Threads.nthreads() > 1
            model = make_fake_gram(rho=7.0e-9, T=250.0)
            p = ODEParams(n_sats=1, args=cfg_gram)
            pos_e, vel_e = elliptical_apoapsis_state(60e3, 300e3)
            withenv("SPACEAGORA_GRAM_ISOLATED_POOL" => "on",
                    "SPACEAGORA_GRAM_TRACK_CACHE_MAX_NPOS" => "16") do
                cache = CB.GramTrackCache()
                CB._gram_track_cache_refresh!(cache, model, p, pos_e, vel_e, 100e3, 0.0, 0.0, 0.0, 60.0, 8)
                @test cache.valid == true
                # Every knot at/below 2000 km on a keplerian run routes to the fake core.
                @test all(x -> x == 7.0e-9 || x == 0.0, cache.rhos)
                @test any(==(7.0e-9), cache.rhos)
                @test length(p.shared_buffers.gram_isolated_pool_models) >= 2
            end
        end
    end

    @testset "refresh: GRAM trajectory driver branch" begin
        # Production shape: EM wrapper -> core with :gram/:gram_atmosphere
        # fields -> module driver, reached via the real property forwarding.
        base = make_traj_gram()
        @test CB._gram_track_trajectory_supported(base) == true
        # Surrogate forwards through two wrapper levels to the same core.
        surrogate = EM.GRAMAtmosphereModelSurrogate(base, "", nothing)
        @test CB._gram_track_trajectory_supported(surrogate) == true
        # Non-module drivers are probed with hasproperty.
        nt_driver = EM.GRAMAtmosphereModel(
            FakeTrajCore((generate_trajectory = fake_generate_trajectory,), nothing)
        )
        @test CB._gram_track_trajectory_supported(nt_driver) == true
        # Drivers without generate_trajectory are unsupported.
        no_gen = EM.GRAMAtmosphereModel(FakeTrajCore((initialize! = identity,), nothing))
        @test CB._gram_track_trajectory_supported(no_gen) == false
        @test CB._gram_track_trajectory_supported(fallback_probe) == false
        @test CB._gram_track_trajectory_supported(make_fake_gram()) == false  # no :gram property

        p = ODEParams(n_sats=1, args=probe_config(density_model=surrogate))
        pos_c, vel_c = circular_state(500e3)
        withenv("SPACEAGORA_GRAM_TRACK_CACHE_MAX_NPOS" => "12",
                "SPACEAGORA_GRAM_ISOLATED_POOL" => "off") do
            cache = CB.GramTrackCache()
            FAKE_TRAJ_KWARGS[] = nothing
            out = CB._gram_track_cache_refresh!(cache, surrogate, p, pos_c, vel_c, 500e3, 0.1, 0.2, 7.0, 60.0, 8)
            @test cache.valid == true
            kw = FAKE_TRAJ_KWARGS[]
            @test kw !== nothing
            @test kw[:initial_height] ≈ 500.0      # km
            @test kw[:initial_latitude] ≈ rad2deg(0.1)
            @test kw[:initial_longitude] ≈ rad2deg(0.2)
            @test kw[:initial_elapsed_time] == 7.0
            @test kw[:n_points] == length(cache.times)
            @test cache.times[1] == 7.0
            @test out[1] == cache.rhos[1]

            # Global-lock variant of the trajectory generation.
            old_lock_fn = EM._GRAM_USE_GLOBAL_LOCK_FN[]
            EM._GRAM_USE_GLOBAL_LOCK_FN[] = () -> true
            try
                cache_l = CB.GramTrackCache()
                CB._gram_track_cache_refresh!(cache_l, surrogate, p, pos_c, vel_c, 500e3, 0.1, 0.2, 7.0, 60.0, 8)
                @test cache_l.valid == true
            finally
                EM._GRAM_USE_GLOBAL_LOCK_FN[] = old_lock_fn
            end
        end
    end

    @testset "refresh: failure fallback" begin
        model = BatchThrowDensityModel()
        p = ODEParams(n_sats=1, args=probe_config(density_model=model))
        pos_c, vel_c = circular_state(500e3)

        withenv("SPACEAGORA_GRAM_TRACK_CACHE_MAX_NPOS" => "16",
                "SPACEAGORA_GRAM_ISOLATED_POOL" => "off",
                "SPACEAGORA_GRAM_PROFILE" => "1") do
            CB._gram_runtime_stats_reset!()
            CB._gram_track_cache_warning_emitted[] = false
            cache = CB.GramTrackCache()
            out = @test_logs (:warn, r"GRAM track cache refresh failed") match_mode=:any begin
                CB._gram_track_cache_refresh!(cache, model, p, pos_c, vel_c, 500e3, 0.0, 0.0, 0.0, 60.0, 8)
            end
            @test cache.valid == false
            @test CB._gram_track_cache_warning_emitted[] == true
            # Fallback answer is the direct scalar query.
            @test out[1] == 1.23e-7
            @test out[2] == 199.0
            @test out[3] == SVector{3, Float64}(0.25, 0.5, 0.75)
            stats = CB._gram_runtime_stats_snapshot()
            @test stats.refresh_failures == 1

            # A second failure stays quiet (warning already emitted).
            cache2 = CB.GramTrackCache()
            out2 = @test_logs min_level=Base.CoreLogging.Warn begin
                CB._gram_track_cache_refresh!(cache2, model, p, pos_c, vel_c, 500e3, 0.0, 0.0, 0.0, 60.0, 8)
            end
            @test cache2.valid == false
            @test out2[1] == 1.23e-7
            CB._gram_runtime_stats_reset!()
            CB._gram_track_cache_warning_emitted[] = false
        end
    end
end
