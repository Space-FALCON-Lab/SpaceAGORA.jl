using LinearAlgebra, Statistics, Plots, LaTeXStrings, Printf, StaticArrays, DelimitedFiles

if !@isdefined(MU);      const MU      = 3.986004418e14; end
if !@isdefined(C);       const C       = 3.0e8;          end
if !@isdefined(R_EARTH); const R_EARTH = 6_378_137.0;    end
if !@isdefined(R_ATMDEF)
    const R_ATMDEF = 6_478_137.0
end
@inline idx(i, off) = 6*(i-1) + off

include("../functions/1_LOS_Metrics.jl")
include("../functions/2_Laser_Forces_ver3.jl")
include("../functions/4_Diagnostics.jl")
include("../functions/5_OE_Converters.jl")
include("../functions/6_OE_and_dv_in_RTN.jl")
include("../functions/12_CSV_Write_Read.jl")

# ── Target CSV ──────────────────────────────────────────────────────────────
csv_path = normpath(joinpath(@__DIR__, "..", "output", "CSV",
    "target_h1050km_i0.5deg",
    "timeseries_N201_T540000s_h1000km_t1050km_ih0.0deg_it0.5deg_B100_Pin1e+04_rmin0m_rmax2e+05_J2F_toyT.csv"))

output_dir = normpath(joinpath(@__DIR__, "..", "output", "images_reconstructed_from_csv", "delta_v_RTN"))
mkpath(output_dir)
csv_tag    = splitext(basename(csv_path))[1]
cache_path = joinpath(output_dir, "dv_RTN_cache_$(csv_tag).csv")

if isfile(cache_path)
    # ── Load from cache ──────────────────────────────────────────────────────
    raw  = readdlm(cache_path, ',', Float64; skipstart=1)
    t    = raw[:, 1]
    oc   = raw[:, 2]
    dv_R = raw[:, 3]
    dv_T = raw[:, 4]
    dv_N = raw[:, 5]
    println("Loaded Δv RTN cache: $cache_path")
else
    # ── Compute from CSV ─────────────────────────────────────────────────────
    sol, metadata, p = load_timeseries_csv(csv_path)

    # Compute cumulative Δv in RTN for the target satellite (index N)
    t, Δv_hist = delta_v_RTN_time_series(sol, p)
    target_dv  = Δv_hist[p[:N]]   # 3 × length(t)

    dv_R = target_dv[1, :]
    dv_T = target_dv[2, :]
    dv_N = target_dv[3, :]

    # Compute orbit-count axis using instantaneous semi-major axis
    oc = Vector{Float64}(undef, length(t))
    oc[1] = 0.0
    for k in 2:length(t)
        u = sol.u[k]
        r = @SVector [u[idx(p[:N],1)], u[idx(p[:N],2)], u[idx(p[:N],3)]]
        v = @SVector [u[idx(p[:N],4)], u[idx(p[:N],5)], u[idx(p[:N],6)]]
        a  = rv2coe(r, v, MU).a
        T  = 2π * sqrt(a^3 / MU)
        oc[k] = oc[k-1] + (t[k] - t[k-1]) / T
    end

    # Save cache
    open(cache_path, "w") do io
        writedlm(io, ["t" "orbits" "dv_R" "dv_T" "dv_N"], ',')
        writedlm(io, hcat(t, oc, dv_R, dv_T, dv_N), ',')
    end
    println("Saved Δv RTN cache: $cache_path")
end

# ── Two separate figures: (dv_R, dv_T) and dv_N ─────────────────────────────
label_fontsize = 14   # axis label font size
tick_fontsize  = 14   # tick font size
fig_size       = (600, 500)   # figure size in pixels (width, height)
left_margin    = 12Plots.mm   # increase if y-label is clipped
bottom_margin  = 8Plots.mm    # increase if x-label is clipped

common = (lw=2,
          title="", titlefontsize=1,
          formatter=:plain,
          xlabel=L"\mathrm{Orbits}",
          ylabel=L"\Delta v, \; \mathrm{m/s}",
          xguidefontsize=label_fontsize, yguidefontsize=label_fontsize,
          xtickfontsize=tick_fontsize,   ytickfontsize=tick_fontsize,
          xguidefontfamily="Computer Modern",
          yguidefontfamily="Computer Modern",
          tickfontfamily="Computer Modern",
          legendfontsize=label_fontsize,
          legendfontfamily="Computer Modern",
          left_margin=left_margin, bottom_margin=bottom_margin,
          size=fig_size)

p_RT = plot(oc, dv_R; label=L"\Delta v_R", color=:red,   legend=:topleft, common...)
plot!(p_RT, oc, dv_T; label=L"\Delta v_T", color=:blue,  lw=2)

p_N  = plot(oc, dv_N; label=L"\Delta v_N", color=:green, legend=:topleft, common...)

for (sp, tag) in zip((p_RT, p_N), ("RT", "N"))
    display(sp)
    savefig(sp, joinpath(output_dir, "dv_$(tag)_toy_problem_$(csv_tag).png"))
    savefig(sp, joinpath(output_dir, "dv_$(tag)_toy_problem_$(csv_tag).pdf"))
    println("Saved: dv_$(tag)_toy_problem_$(csv_tag).png/.pdf")
end
