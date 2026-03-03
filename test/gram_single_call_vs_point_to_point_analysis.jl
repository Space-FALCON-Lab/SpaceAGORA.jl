const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

using Statistics
using Printf
using CSV
using DataFrames
using StaticArrays

include(joinpath(@__DIR__, "gram_interpolation_vs_point_to_point_analysis.jl"))

const LARGE_ALT_TOL_M = 1.0e9
const LARGE_ANG_TOL_RAD = Float64(pi)
const SINGLE_CALL_POINTS = 128

function _evaluate_interpolated_single_call(track, density_model, p)
    n = length(track.t)
    rho = Vector{Float64}(undef, n)
    temp = Vector{Float64}(undef, n)
    wind = Vector{SVector{3, Float64}}(undef, n)
    cache = CB.GramTrackCache()
    refresh_count = 0
    npos_used = Int[]
    horizon_s = max(1e-3, track.t[end] - track.t[1])
    segment_end_t = track.t[end]

    elapsed_s = @elapsed begin
        @inbounds for k in 1:n
            seg = CB._gram_track_cache_ready(
                cache,
                track.t[k],
                track.alt[k],
                track.lat[k],
                track.lon[k],
                LARGE_ALT_TOL_M,
                LARGE_ANG_TOL_RAD
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
                    SINGLE_CALL_POINTS,
                    LARGE_ALT_TOL_M,
                    LARGE_ANG_TOL_RAD,
                    0.0,
                    segment_end_t
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

function _trajectory_error_table(scenario::String, track, direct, interp)
    n = length(track.t)
    rho_abs = [abs(interp.rho[k] - direct.rho[k]) for k in 1:n]
    rho_rel = [rho_abs[k] / max(abs(direct.rho[k]), 1e-12) for k in 1:n]
    temp_abs = [abs(interp.temp[k] - direct.temp[k]) for k in 1:n]
    wind_abs = [norm(interp.wind[k] - direct.wind[k]) for k in 1:n]
    lat_deg = rad2deg.(track.lat)
    lon_deg = rad2deg.(track.lon)

    return DataFrame(
        scenario=fill(scenario, n),
        t_s=track.t,
        alt_m=track.alt,
        lat_deg=lat_deg,
        lon_deg=lon_deg,
        rho_direct=direct.rho,
        rho_interp=interp.rho,
        rho_abs_err=rho_abs,
        rho_rel_err=rho_rel,
        temp_direct=direct.temp,
        temp_interp=interp.temp,
        temp_abs_err=temp_abs,
        wind_direct=[norm(w) for w in direct.wind],
        wind_interp=[norm(w) for w in interp.wind],
        wind_abs_err=wind_abs
    )
end

function run_single_call_analysis()
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

    summary_rows = NamedTuple[]
    error_tables = DataFrame[]

    println("Single-call GRAM vs point-to-point analysis")
    println(@sprintf("Forced tolerances: alt=%.3e m, ang=%.3f deg", LARGE_ALT_TOL_M, rad2deg(LARGE_ANG_TOL_RAD)))

    for sc in scenarios
        println("Running scenario: $(sc.name)")
        direct = _evaluate_direct(sc.track, density_model, p)
        single = _evaluate_interpolated_single_call(sc.track, density_model, p)
        summary = _summarize_case(sc.name, sc.track, direct, single)
        push!(
            summary_rows,
            merge(
                summary,
                (
                    mode="single_call_large_tol",
                    forced_alt_tol_m=LARGE_ALT_TOL_M,
                    forced_ang_tol_deg=rad2deg(LARGE_ANG_TOL_RAD),
                    forced_points=SINGLE_CALL_POINTS
                )
            )
        )
        push!(error_tables, _trajectory_error_table(sc.name, sc.track, direct, single))
    end

    summary_df = DataFrame(summary_rows)
    errors_df = vcat(error_tables...)

    mkpath(joinpath(REPO_ROOT, "output"))
    out_summary = joinpath(REPO_ROOT, "output", "gram_single_call_vs_point_to_point_summary.csv")
    out_errors = joinpath(REPO_ROOT, "output", "gram_single_call_vs_point_to_point_trajectory_errors.csv")
    CSV.write(out_summary, summary_df)
    CSV.write(out_errors, errors_df)

    println("\nSummary:")
    show(summary_df, allrows=true, allcols=true)
    println("\n\nSaved CSVs:")
    println("  summary:   $out_summary")
    println("  trajectory $out_errors")
    return (summary=summary_df, errors=errors_df)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_single_call_analysis()
end
