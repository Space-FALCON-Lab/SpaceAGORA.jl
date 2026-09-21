using LinearAlgebra, Plots, Printf, StaticArrays, DelimitedFiles

if !@isdefined(MU);       const MU       = 3.986004418e14;          end
if !@isdefined(C);        const C        = 3.0e8;                   end
if !@isdefined(R_EARTH);  const R_EARTH  = 6_378_137.0;             end
if !@isdefined(R_ATMDEF); const R_ATMDEF = R_EARTH + 100_000.0;     end
if !@isdefined(Ẑ);        const Ẑ        = SVector(0.0, 0.0, 1.0); end
@inline idx(i, off) = 6*(i-1) + off

include("../../functions/1_LOS_Metrics.jl")
include("../../functions/2_Laser_Forces_ver2.jl")
include("../../functions/4_Diagnostics.jl")
include("../../functions/5_OE_Converters.jl")
include("../../functions/6_OE_and_dv_in_RTN.jl")
include("../../functions/12_CSV_Write_Read.jl")

"""
    plot_dv(h, i_deg; nu=nothing) -> NamedTuple

Compute and plot the dΔv/dt slope in RTN components vs number of helpers for
all CSVs matching the given orbital parameters. Slopes and full Δv time series
are cached for fast re-runs.

# Arguments
- `h`:     target altitude in km (e.g. 1000)
- `i_deg`: inclination in degrees (e.g. 0.0)
- `nu`:    true anomaly offset in degrees (optional, e.g. -0.75)

# Returns
A named tuple with fields:
- `plt`           — the dΔv/dt vs N_helpers plot
- `N_helpers_vec` — `Vector{Int}` number of helpers per CSV
- `slope_R`       — `Vector{Float64}` dΔv_R/dt slopes (m/s²)
- `slope_T`       — `Vector{Float64}` dΔv_T/dt slopes (m/s²)
- `slope_N`       — `Vector{Float64}` dΔv_N/dt slopes (m/s²)
- `slope_mat`     — `Matrix{Float64}` (n_files × 3) columns = [R, T, N]
- `dv_timeseries` — `Vector` of NamedTuples `(t, R, T, N)`, one per N_helpers
                    each field is a `Vector{Float64}` of length n_timesteps
- `dv_mag`        — `Vector{Float64}` final total |Δv| magnitude per N_helpers (m/s)
- `dv_change_R`   — `Vector{Float64}` Δv_R(end) − Δv_R(start) per N_helpers (m/s)
- `dv_change_T`   — `Vector{Float64}` Δv_T(end) − Δv_T(start) per N_helpers (m/s)
- `dv_change_N`   — `Vector{Float64}` Δv_N(end) − Δv_N(start) per N_helpers (m/s)
"""
function plot_dv(h::Real, i_deg::Real;
                 nu::Union{Real,Nothing} = nothing)
    # --- Build folder name from parameters ---
    folder = if nu === nothing
        @sprintf("target_h%dkm_i%.1fdeg", h, i_deg)
    else
        @sprintf("target_h%dkm_i%.1fdeg_nu%.2fdeg", h, i_deg, nu)
    end

    data_dir   = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", folder))
    dir_tag    = basename(data_dir)
    output_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "images_reconstructed_from_csv", "delta_v"))
    mkpath(output_dir)

    # Cache layout:
    #   slope      →  delta_v/plot_cache/slope/<dir_tag>/<dir_tag>_slopes.csv
    #   timeseries →  delta_v/plot_cache/timeseries/<dir_tag>/N<N_total>_<dir_tag>_timeseries.csv
    slopes_cache_dir  = joinpath(output_dir, "plot_cache", "slope", dir_tag)
    slopes_cache_path = joinpath(slopes_cache_dir, "$(dir_tag)_slopes.csv")

    isdir(data_dir) || error("Directory not found: $data_dir")
    csv_files = sort(filter(f -> endswith(f, ".csv"), readdir(data_dir; join=true)),
                     by = f -> parse(Int, match(r"_N(\d+)_", basename(f))[1]))

    # All timeseries for the same case share one folder; each file is named by N_total
    ts_cache_dir   = joinpath(output_dir, "plot_cache", "timeseries", dir_tag)
    ts_cache_paths = [joinpath(ts_cache_dir,
                               "N$(parse(Int, match(r"_N(\d+)_", basename(f))[1]))_$(dir_tag)_timeseries.csv")
                      for f in csv_files]

    N_helpers_vec = Int[]
    slope_rows    = Vector{Vector{Float64}}()
    dv_timeseries = NamedTuple[]

    linfit(x, y) = (hcat(x, ones(length(x))) \ y)[1]   # returns slope only

    all_cached = isfile(slopes_cache_path) && all(isfile, ts_cache_paths)

    if all_cached
        raw           = readdlm(slopes_cache_path, ',', Float64; skipstart=1)
        raw           = ndims(raw) == 1 ? reshape(raw, 1, :) : raw  # single-row fix
        N_helpers_vec = Int.(raw[:, 1])
        slope_rows    = [raw[i, 2:4] for i in axes(raw, 1)]
        println("Loaded dv slopes cache: $slopes_cache_path")

        for N_h in N_helpers_vec
            ts_path = joinpath(ts_cache_dir, "N$(N_h + 1)_$(dir_tag)_timeseries.csv")
            ts_raw  = readdlm(ts_path, ',', Float64; skipstart=1)
            push!(dv_timeseries, (t=ts_raw[:, 1], R=ts_raw[:, 2], T=ts_raw[:, 3], N=ts_raw[:, 4]))
        end
        println("Loaded dv time-series caches ($(length(dv_timeseries)) runs)")
    else
        for (csv_path, ts_path) in zip(csv_files, ts_cache_paths)
            local sol, _, p  = load_timeseries_csv(csv_path)
            local t, Δv_hist = delta_v_RTN_time_series(sol, p)
            local dv         = Δv_hist[p[:N]]          # target = last satellite
            local N_h        = p[:N] - 1               # number of helpers
            local slopes     = [linfit(t, dv[row, :]) for row in 1:3]
            push!(N_helpers_vec, N_h)
            push!(slope_rows, slopes)
            push!(dv_timeseries, (t=t, R=dv[1, :], T=dv[2, :], N=dv[3, :]))
            # Save time-series cache for this N_helpers
            mkpath(ts_cache_dir)
            open(ts_path, "w") do io
                writedlm(io, ["t" "dv_R" "dv_T" "dv_N"], ',')
                writedlm(io, hcat(t, dv[1, :], dv[2, :], dv[3, :]), ',')
            end
            @printf("  N_helpers=%d  dΔv/dt: R=%+.4e  T=%+.4e  N=%+.4e  m/s²\n",
                    N_h, slopes[1], slopes[2], slopes[3])
        end
        mkpath(slopes_cache_dir)
        open(slopes_cache_path, "w") do io
            writedlm(io, ["N_helpers" "slope_R" "slope_T" "slope_N"], ',')
            writedlm(io, hcat(N_helpers_vec, stack(slope_rows)'), ',')
        end
        println("Saved dv slopes cache: $slopes_cache_path")
        println("Saved dv time-series caches ($(length(N_helpers_vec)) runs)")
    end

    slope_mat = stack(slope_rows)'   # (n_files × 3)
    slope_R   = slope_mat[:, 1]
    slope_T   = slope_mat[:, 2]
    slope_N   = slope_mat[:, 3]

    # Derived quantities from dv_timeseries
    dv_mag      = [sqrt(ts.R[end]^2 + ts.T[end]^2 + ts.N[end]^2) for ts in dv_timeseries]
    dv_change_R = [ts.R[end] - ts.R[1] for ts in dv_timeseries]
    dv_change_T = [ts.T[end] - ts.T[1] for ts in dv_timeseries]
    dv_change_N = [ts.N[end] - ts.N[1] for ts in dv_timeseries]

    # --- Plot ---
    labels = ["Δv_R", "Δv_T", "Δv_N"]
    styles = [(:blue, :solid), (:orange, :dash), (:green, :dashdot)]

    plt = plot(title="dΔv/dt (RTN) vs Number of Helpers ($dir_tag)",
               xlabel="Number of helpers", ylabel="dΔv/dt (m/s²)",
               legend=:outertopright, yformatter=:scientific,
               xticks=N_helpers_vec)

    for (slopes_vec, lbl, (clr, ls)) in zip([slope_R, slope_T, slope_N], labels, styles)
        plot!(plt, N_helpers_vec, slopes_vec; label=lbl, color=clr, linestyle=ls,
              marker=:circle, ms=0, lw=5)
    end

    return (plt=plt,
            N_helpers_vec=N_helpers_vec,
            slope_R=slope_R,
            slope_T=slope_T,
            slope_N=slope_N,
            slope_mat=slope_mat,
            dv_timeseries=dv_timeseries,
            dv_mag=dv_mag,
            dv_change_R=dv_change_R,
            dv_change_T=dv_change_T,
            dv_change_N=dv_change_N)
end

