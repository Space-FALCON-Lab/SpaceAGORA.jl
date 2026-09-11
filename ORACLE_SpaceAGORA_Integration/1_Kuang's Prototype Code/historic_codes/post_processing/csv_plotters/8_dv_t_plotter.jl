using LinearAlgebra, Plots, Printf, StaticArrays, DelimitedFiles, Serialization

const MU       = 3.986004418e14
const C        = 3.0e8
const R_EARTH  = 6_378_137.0
const R_ATMDEF = R_EARTH + 100_000.0    # default "atmosphere radius" (Kármán line) [m]
@inline idx(i, off) = 6*(i-1) + off

include("../../functions/4_Diagnostics.jl")
include("../../functions/5_OE_Converters.jl")
include("../../functions/1_LOS_Metrics.jl")
include("../../functions/2_Laser_Forces_ver2.jl")
include("../../functions/6_OE_and_dv_in_RTN.jl")
include("../../functions/12_CSV_Write_Read.jl")

# =============================================================================
# --- Config: choose the test case ---
# =============================================================================
h_t    = 1000.0    # target altitude (km)
h_h    = 1000.0    # helper altitude (km)
it_deg = 0.0       # target inclination (deg)
ih_deg = 0.0       # helper inclination (deg)
nu_deg = -0.75     # true-anomaly offset (deg); set 0.0 for no-nu directories
N_h    = 50        # number of helper satellites

# =============================================================================
# --- Resolve directory and CSV file ---
# =============================================================================
csv_base = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV"))

# Build directory name: "target_h<h_t>km_i<it>deg[_nu<nu>deg]"
dir_name = if nu_deg == 0.0
    @sprintf("target_h%.0fkm_i%.1fdeg", h_t, it_deg)
else
    @sprintf("target_h%.0fkm_i%.1fdeg_nu%.2fdeg", h_t, it_deg, nu_deg)
end
data_dir = joinpath(csv_base, dir_name)

@assert isdir(data_dir) "Directory not found: $data_dir"

# Find matching CSV: N = N_h + 1, h_h, ih
N_total  = N_h + 1
pattern  = Regex(string(
    "_N$(N_total)_",
    ".*",
    @sprintf("_h%.0fkm_t%.0fkm_ih%.1fdeg_it%.1fdeg", h_h, h_t, ih_deg, it_deg)
))
matches  = filter(f -> occursin(pattern, f), readdir(data_dir))
@assert length(matches) == 1 "Expected 1 matching CSV, found $(length(matches)) in $data_dir\nPattern: $pattern"

csv_path   = joinpath(data_dir, matches[1])
csv_tag    = replace(matches[1], ".csv" => "")
output_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "images_reconstructed_from_csv", "delta_v"))
dv_cache_dir = joinpath(output_dir, "dv_cache")
mkpath(dv_cache_dir)

println("CSV: ", csv_path)

# =============================================================================
# --- Load or compute Δv time series ---
# =============================================================================
dv_cache = joinpath(dv_cache_dir, "$(csv_tag)_dv_cache.jls")
local t, Δv_hist, p

if isfile(dv_cache)
    cached = deserialize(dv_cache)
    if cached isa Tuple{Vector{Float64}, Vector, Dict}
        t, Δv_hist, p = cached
        println("Loaded from cache: $dv_cache")
    else
        println("Stale cache — recomputing...")
        local sol
        sol, _, p = load_timeseries_csv(csv_path)
        t, Δv_hist = delta_v_RTN_time_series(sol, p)
        serialize(dv_cache, (t, Δv_hist, p))
    end
else
    local sol
    sol, _, p = load_timeseries_csv(csv_path)
    t, Δv_hist = delta_v_RTN_time_series(sol, p)
    serialize(dv_cache, (t, Δv_hist, p))
    println("Computed & cached: $dv_cache")
end

dv = Δv_hist[p[:N]]   # target = last satellite

# =============================================================================
# --- Linear fit: dv ≈ slope·t + intercept ---
# =============================================================================
labels   = ["Δv_R", "Δv_T", "Δv_N"]
suffixes = ["_R_only", "_T_only", "_N_only"]

A    = hcat(t, ones(length(t)))
fits = [A \ dv[row, :] for row in 1:3]

println("\nLinear fit rates of change (Δv = slope·t + intercept):")
for (lbl, f) in zip(labels, fits)
    @printf("  %-6s  slope = %+.6e m/s²   intercept = %+.6e m/s\n", lbl, f[1], f[2])
end

# =============================================================================
# --- Plot ---
# =============================================================================

# individual RTN components with linear fit overlay
for (row, lbl, suf, f) in zip(1:3, labels, suffixes, fits)
    fit_line = f[1] .* t .+ f[2]
    plt = plot(t, dv[row, :]; label=lbl, xlabel="t (s)", ylabel="m/s",
               title="$lbl component (RTN)  [$csv_tag]",
               color=:grey, alpha=0.7, lw=1)
    plot!(plt, t, fit_line;
          label="fit ($(Printf.@sprintf("%+.3e", f[1])) m/s²)",
          color=:red, lw=2, linestyle=:dash)
    savefig(plt, joinpath(output_dir, "dv_RTN_$(csv_tag)$(suf).png"))
end

# all components together
plt_all = plot(t, dv[1, :]; label="Δv_R", xlabel="t (s)", ylabel="m/s",
               title="Δv vs time (RTN)  [$csv_tag]")
plot!(plt_all, t, dv[2, :]; label="Δv_T")
plot!(plt_all, t, dv[3, :]; label="Δv_N")
styles = [(:red, :dash), (:green, :dash), (:blue, :dash)]
for (row, lbl, f, (clr, ls)) in zip(1:3, labels, fits, styles)
    plot!(plt_all, t, f[1] .* t .+ f[2];
          label="$(lbl) fit ($(Printf.@sprintf("%+.3e", f[1])) m/s²)",
          lw=2, linestyle=ls, color=clr)
end
savefig(plt_all, joinpath(output_dir, "dv_RTN_$(csv_tag).png"))

display(plt_all)
println("Saved to: ", output_dir)
