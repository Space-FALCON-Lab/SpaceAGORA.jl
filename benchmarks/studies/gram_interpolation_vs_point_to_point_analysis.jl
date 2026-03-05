const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

using Statistics
using Printf
using Dates
using LinearAlgebra
using StaticArrays
using CSV
using DataFrames

include(joinpath(REPO_ROOT, "src", "simulation_model", "SimulationModel.jl"))
using .SimulationModel

const EM = SimulationModel.EnvironmentModels
const CB = SimulationModel.SimulationCallbacks

const STRICT_RHO_REL_MAX = 0.10
const STRICT_TEMP_ABS_MAX = 0.6
const STRICT_WIND_ABS_MAX = 120.0

@inline function _pctl(x::Vector{Float64}, p::Float64)
    isempty(x) && return 0.0
    xs = sort(x)
    idx = clamp(Int(ceil(p * length(xs))), 1, length(xs))
    return xs[idx]
end

function _make_spacecraft(planet)::SpacecraftModel
    root = Link{0}(root=true, m=120.0, ref_area=1.2)
    ic = InitialCondition(
        ra=planet.Rp_e + 520e3,
        rp=planet.Rp_e + 500e3,
        i=25.0,
        ω=10.0,
        Ω=20.0,
        ν=170.0
    )
    return SpacecraftModel(
        joints=Joint[],
        links=Link[root],
        root=root,
        instant_actuation=true,
        prop_mass=0.0,
        inertia_tensor=root.inertia,
        n_reaction_wheels=0,
        n_thrusters=0,
        initial_condition=ic,
        id=1
    )
end

function _make_args(planet, density_model, initial_time)::SimulationConfiguration
    sc = _make_spacecraft(planet)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false,
            verbose=false,
            generate_plots=false,
            normalize=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=60.0 * 60.0,
            orientation_sim=false,
            num_steps_to_save=1000
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=120.0,
            density_model=density_model,
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=true
        ),
        dynamics_model=DynamicsModel([sc], (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=initial_time,
        integration_tolerances=IntegrationTolerances()
    )
end

function _generate_track(
    planet;
    rp_alt_m::Float64,
    ra_alt_m::Float64,
    i_deg::Float64,
    Ω_deg::Float64,
    ω_deg::Float64,
    ν_start_deg::Float64,
    ν_end_deg::Float64,
    n_samples::Int
)
    rp = planet.Rp_e + rp_alt_m
    ra = planet.Rp_e + ra_alt_m
    a = 0.5 * (ra + rp)
    e = (ra - rp) / (ra + rp)
    i = deg2rad(i_deg)
    Ω = deg2rad(Ω_deg)
    ω = deg2rad(ω_deg)
    νs = collect(range(deg2rad(ν_start_deg), deg2rad(ν_end_deg), length=n_samples))

    pos = Vector{SVector{3, Float64}}(undef, n_samples)
    vel = Vector{SVector{3, Float64}}(undef, n_samples)
    alt = Vector{Float64}(undef, n_samples)
    lat = Vector{Float64}(undef, n_samples)
    lon = Vector{Float64}(undef, n_samples)
    t = Vector{Float64}(undef, n_samples)
    t[1] = 0.0

    for k in 1:n_samples
        oe = SVector{7, Float64}(a, e, i, Ω, ω, νs[k], 0.0)
        r_k, v_k = CB.orbitalelemtorv(oe, planet)
        pos[k] = SVector{3, Float64}(r_k)
        vel[k] = SVector{3, Float64}(v_k)
        rp_k, _ = CB.r_intor_p!(pos[k], vel[k], planet)
        alt[k], lat[k], lon[k] = CB.rtolatlong(rp_k, planet)

        if k > 1
            ds = norm(pos[k] - pos[k - 1])
            vbar = max(1.0, 0.5 * (norm(vel[k]) + norm(vel[k - 1])))
            t[k] = t[k - 1] + ds / vbar
        end
    end

    return (t=t, pos=pos, vel=vel, alt=alt, lat=lat, lon=lon)
end

function _evaluate_direct(track, density_model, p)
    n = length(track.t)
    rho = Vector{Float64}(undef, n)
    temp = Vector{Float64}(undef, n)
    wind = Vector{SVector{3, Float64}}(undef, n)
    elapsed_s = @elapsed begin
        @inbounds for k in 1:n
            rho[k], temp[k], wind[k] = Base.invokelatest(
                getDensity,
                density_model,
                track.alt[k],
                track.lat[k],
                track.lon[k],
                track.t[k],
                true,
                p
            )
        end
    end
    return (rho=rho, temp=temp, wind=wind, elapsed_s=elapsed_s)
end

function _evaluate_interpolated(track, density_model, p, cache_cfg)
    n = length(track.t)
    rho = Vector{Float64}(undef, n)
    temp = Vector{Float64}(undef, n)
    wind = Vector{SVector{3, Float64}}(undef, n)
    cache = CB.GramTrackCache()
    refresh_count = 0
    npos_used = Int[]
    elapsed_s = @elapsed begin
        @inbounds for k in 1:n
            horizon_s, alt_tol_m, ang_tol_rad, n_points = CB._gram_track_cache_profile(cache_cfg, p, track.alt[k])
            seg = CB._gram_track_cache_ready(
                cache,
                track.t[k],
                track.alt[k],
                track.lat[k],
                track.lon[k],
                alt_tol_m,
                ang_tol_rad
            )
            if seg === nothing
                refresh_count += 1
                rho[k], temp[k], wind[k] = Base.invokelatest(
                    CB._gram_track_cache_refresh!,
                    cache,
                    density_model,
                    p,
                    track.pos[k],
                    track.vel[k],
                    track.alt[k],
                    track.lat[k],
                    track.lon[k],
                    track.t[k],
                    horizon_s,
                    n_points,
                    alt_tol_m,
                    ang_tol_rad,
                    cache_cfg.transition_band_m
                )
                push!(npos_used, length(cache.times))
            else
                idx, x = seg
                rho[k], temp[k], wind[k] = CB._gram_track_cache_eval(cache, idx, x)
            end
        end
    end
    return (
        rho=rho,
        temp=temp,
        wind=wind,
        elapsed_s=elapsed_s,
        refresh_count=refresh_count,
        npos_used=npos_used
    )
end

function _summarize_case(case_name::String, track, direct, interp)
    n = length(track.t)
    rho_abs = [abs(interp.rho[k] - direct.rho[k]) for k in 1:n]
    rho_rel = [rho_abs[k] / max(abs(direct.rho[k]), 1e-12) for k in 1:n]
    temp_abs = [abs(interp.temp[k] - direct.temp[k]) for k in 1:n]
    wind_abs = [norm(interp.wind[k] - direct.wind[k]) for k in 1:n]

    npos = interp.npos_used
    npos_mean = isempty(npos) ? 0.0 : mean(Float64.(npos))
    npos_min = isempty(npos) ? 0 : minimum(npos)
    npos_max = isempty(npos) ? 0 : maximum(npos)

    return (
        scenario=case_name,
        samples=n,
        t_span_s=track.t[end] - track.t[1],
        alt_min_km=minimum(track.alt) * 1e-3,
        alt_max_km=maximum(track.alt) * 1e-3,
        rho_abs_max=maximum(rho_abs),
        rho_abs_p95=_pctl(rho_abs, 0.95),
        rho_rel_max=maximum(rho_rel),
        rho_rel_p95=_pctl(rho_rel, 0.95),
        temp_abs_max=maximum(temp_abs),
        temp_abs_p95=_pctl(temp_abs, 0.95),
        wind_abs_max=maximum(wind_abs),
        wind_abs_p95=_pctl(wind_abs, 0.95),
        direct_eval_s=direct.elapsed_s,
        interp_eval_s=interp.elapsed_s,
        speedup=direct.elapsed_s / max(interp.elapsed_s, 1e-9),
        refresh_count=interp.refresh_count,
        refresh_ratio=interp.refresh_count / max(n, 1),
        npos_mean=npos_mean,
        npos_min=npos_min,
        npos_max=npos_max
    )
end

@inline function _regime_for_scenario(name::String)::Symbol
    if name == "orbit"
        return :orbit
    elseif name == "drag_passage"
        return :drag_passage
    elseif name == "entry"
        return :entry
    end
    throw(ArgumentError("Unknown scenario '$name'."))
end

function _candidate_grid(regime::Symbol)
    rows = NamedTuple[]
    if regime == :orbit
        horizons = (6.0, 8.0, 12.0)
        points = (32, 56)
        alt_tols = (1500.0, 3000.0)
        ang_tols = (2.0, 4.0)
    else
        horizons = (0.75, 1.0, 1.5)
        points = (12, 20)
        alt_tols = (300.0, 700.0)
        ang_tols = (0.4, 1.0)
    end

    for h in horizons, n in points, a in alt_tols, g in ang_tols
        push!(rows, (horizon_s=h, points=n, alt_tol_m=a, ang_tol_deg=g))
    end
    return rows
end

function _cfg_with_candidate(base_cfg::CB.GramTrackCacheConfig, regime::Symbol, c)
    if regime == :orbit
        return CB.GramTrackCacheConfig(
            base_cfg.mode,
            base_cfg.entry_horizon_s,
            base_cfg.entry_alt_tol_m,
            base_cfg.entry_ang_tol_rad,
            base_cfg.entry_points,
            Float64(c.horizon_s),
            Float64(c.alt_tol_m),
            deg2rad(Float64(c.ang_tol_deg)),
            Int(c.points),
            base_cfg.transition_band_m
        )
    end

    return CB.GramTrackCacheConfig(
        base_cfg.mode,
        Float64(c.horizon_s),
        Float64(c.alt_tol_m),
        deg2rad(Float64(c.ang_tol_deg)),
        Int(c.points),
        base_cfg.orbit_horizon_s,
        base_cfg.orbit_alt_tol_m,
        base_cfg.orbit_ang_tol_rad,
        base_cfg.orbit_points,
        base_cfg.transition_band_m
    )
end

@inline function _candidate_score(row)::Float64
    # Lower is better: primarily interpolation runtime, with penalties for
    # large interpolation-vs-point deviations and low cache reuse.
    score = row.interp_eval_s
    score += max(0.0, row.rho_rel_max - STRICT_RHO_REL_MAX) * 20.0
    score += max(0.0, row.temp_abs_max - STRICT_TEMP_ABS_MAX) * 2.0
    score += max(0.0, row.wind_abs_max - STRICT_WIND_ABS_MAX) / 60.0
    score += max(0.0, row.refresh_ratio - 0.90) * 2.0
    return score
end

@inline function _is_strict_candidate(row)::Bool
    return row.rho_rel_max <= STRICT_RHO_REL_MAX &&
           row.temp_abs_max <= STRICT_TEMP_ABS_MAX &&
           row.wind_abs_max <= STRICT_WIND_ABS_MAX
end

function _best_strict_candidate(df::DataFrame)
    strict_df = filter(_is_strict_candidate, df)
    nrow(strict_df) == 0 && return nothing
    sort!(strict_df, [:interp_eval_s, :refresh_ratio, :rho_rel_max, :temp_abs_max, :wind_abs_max])
    return strict_df[1, :]
end

@inline function _recommendation_row(scenario::String, row, mode::String, strict_feasible::Bool)
    return (
        scenario=scenario,
        mode=mode,
        strict_feasible=strict_feasible,
        regime=row.regime,
        horizon_s=row.horizon_s,
        points=row.points,
        alt_tol_m=row.alt_tol_m,
        ang_tol_deg=row.ang_tol_deg,
        interp_eval_s=row.interp_eval_s,
        speedup=row.speedup,
        rho_rel_max=row.rho_rel_max,
        temp_abs_max=row.temp_abs_max,
        wind_abs_max=row.wind_abs_max,
        refresh_ratio=row.refresh_ratio,
        npos_mean=row.npos_mean,
        npos_max=row.npos_max,
        score=row.score
    )
end

function _sweep_scenario(
    sc,
    density_model,
    p,
    base_cfg::CB.GramTrackCacheConfig,
    direct
)
    regime = _regime_for_scenario(sc.name)
    candidates = _candidate_grid(regime)
    rows = NamedTuple[]

    for (idx, c) in enumerate(candidates)
        cfg = _cfg_with_candidate(base_cfg, regime, c)
        interp = _evaluate_interpolated(sc.track, density_model, p, cfg)
        summary = _summarize_case(sc.name, sc.track, direct, interp)
        push!(
            rows,
            merge(
                summary,
                (
                    regime=String(regime),
                    candidate_idx=idx,
                    horizon_s=Float64(c.horizon_s),
                    points=Int(c.points),
                    alt_tol_m=Float64(c.alt_tol_m),
                    ang_tol_deg=Float64(c.ang_tol_deg)
                )
            )
        )
    end

    scored = map(row -> merge(row, (score=_candidate_score(row),)), rows)
    df = DataFrame(scored)
    sort!(df, [:score, :interp_eval_s, :rho_rel_max, :temp_abs_max, :wind_abs_max])
    return df, df[1, :]
end

function run_analysis()
    initial_time = InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0)
    spice_path = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
    planet = Mars("", spice_path)
    planet.L_PI .= [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]

    density_model = EM.GRAMAtmosphereModel(
        planet_name="mars",
        initial_time=initial_time
    )
    args = _make_args(planet, density_model, initial_time)
    p = ODEParams{1}(args=args)
    cache_cfg = CB._gram_track_cache_config()

    println("GRAM interpolation vs point-to-point analysis")
    println(@sprintf(
        "Cache config: entry(h=%.3fs,n=%d,alt_tol=%.1fm,ang_tol=%.3fdeg), orbit(h=%.3fs,n=%d,alt_tol=%.1fm,ang_tol=%.3fdeg), transition=%.1fm",
        cache_cfg.entry_horizon_s,
        cache_cfg.entry_points,
        cache_cfg.entry_alt_tol_m,
        rad2deg(cache_cfg.entry_ang_tol_rad),
        cache_cfg.orbit_horizon_s,
        cache_cfg.orbit_points,
        cache_cfg.orbit_alt_tol_m,
        rad2deg(cache_cfg.orbit_ang_tol_rad),
        cache_cfg.transition_band_m
    ))

    scenarios = [
        (
            name="drag_passage",
            track=_generate_track(
                planet;
                rp_alt_m=110e3,
                ra_alt_m=4500e3,
                i_deg=25.0,
                Ω_deg=20.0,
                ω_deg=35.0,
                ν_start_deg=-20.0,
                ν_end_deg=30.0,
                n_samples=220
            )
        ),
        (
            name="entry",
            track=_generate_track(
                planet;
                rp_alt_m=80e3,
                ra_alt_m=1500e3,
                i_deg=22.0,
                Ω_deg=35.0,
                ω_deg=10.0,
                ν_start_deg=-45.0,
                ν_end_deg=-2.0,
                n_samples=220
            )
        ),
        (
            name="orbit",
            track=_generate_track(
                planet;
                rp_alt_m=450e3,
                ra_alt_m=450e3,
                i_deg=30.0,
                Ω_deg=15.0,
                ω_deg=0.0,
                ν_start_deg=0.0,
                ν_end_deg=360.0,
                n_samples=220
            )
        )
    ]

    baseline_rows = NamedTuple[]
    sweep_frames = DataFrame[]
    recommendation_rows = NamedTuple[]

    for sc in scenarios
        println("Running scenario: $(sc.name)")
        direct = _evaluate_direct(sc.track, density_model, p)
        interp_default = _evaluate_interpolated(sc.track, density_model, p, cache_cfg)
        baseline = _summarize_case(sc.name, sc.track, direct, interp_default)
        push!(baseline_rows, merge(baseline, (config_name="default",)))

        println("Sweeping scenario: $(sc.name)")
        sweep_df, best = _sweep_scenario(sc, density_model, p, cache_cfg, direct)
        push!(sweep_frames, sweep_df)
        strict_best = _best_strict_candidate(sweep_df)

        push!(
            recommendation_rows,
            _recommendation_row(sc.name, best, "score_optimal", true)
        )

        if strict_best === nothing
            push!(
                recommendation_rows,
                _recommendation_row(sc.name, best, "strict_fallback", false)
            )
        else
            push!(
                recommendation_rows,
                _recommendation_row(sc.name, strict_best, "strict_accuracy", true)
            )
        end
    end

    baseline_df = DataFrame(baseline_rows)
    sweep_df = vcat(sweep_frames...)
    rec_df = DataFrame(recommendation_rows)

    mkpath(joinpath(REPO_ROOT, "output"))
    out_baseline_csv = joinpath(REPO_ROOT, "output", "gram_interpolation_vs_point_to_point.csv")
    out_sweep_csv = joinpath(REPO_ROOT, "output", "gram_interpolation_parameter_sweep.csv")
    out_rec_csv = joinpath(REPO_ROOT, "output", "gram_interpolation_parameter_recommendations.csv")
    CSV.write(out_baseline_csv, baseline_df)
    CSV.write(out_sweep_csv, sweep_df)
    CSV.write(out_rec_csv, rec_df)

    println("\nBaseline Summary:")
    show(baseline_df, allrows=true, allcols=true)
    println("\n\nRecommended Parameters:")
    show(rec_df, allrows=true, allcols=true)
    println("\n\nSaved CSVs:")
    println("  baseline: $out_baseline_csv")
    println("  sweep:    $out_sweep_csv")
    println("  best:     $out_rec_csv")
    return (baseline=baseline_df, sweep=sweep_df, recommendations=rec_df)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_analysis()
end
