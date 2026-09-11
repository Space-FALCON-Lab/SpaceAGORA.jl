using LinearAlgebra, Statistics, Plots, Printf, StaticArrays, DelimitedFiles

if !@isdefined(MU);      const MU      = 3.986004418e14; end
if !@isdefined(C);       const C       = 3.0e8;          end
if !@isdefined(R_EARTH); const R_EARTH = 6_378_137.0;    end
if !@isdefined(R_ATMDEF)
    const R_ATMDEF = 6_478_137.0
end
@inline idx(i, off) = 6*(i-1) + off

include("../functions/4_Diagnostics.jl")
include("../functions/5_OE_Converters.jl")
include("../functions/6_OE_and_dv_in_RTN.jl")
include("../functions/12_CSV_Write_Read.jl")

const N_PTS_OE_TOY = 2000

# ── Target CSV ──────────────────────────────────────────────────────────────
csv_path = normpath(joinpath(@__DIR__, "..", "output", "CSV",
    "target_h1050km_i0.5deg",
    "timeseries_N201_T540000s_h1000km_t1050km_ih0.0deg_it0.5deg_B100_Pin1e+04_rmin0m_rmax2e+05_J2F_toyT.csv"))

output_dir = normpath(joinpath(@__DIR__, "..", "output", "images_reconstructed_from_csv", "orbital_elements"))
mkpath(output_dir)
csv_tag    = splitext(basename(csv_path))[1]
cache_path = joinpath(output_dir, "delta_OE_cache_$(csv_tag).csv")

# orbit-count helper (same as plot_OE)
function compute_orbits_curve(a_vec, t_vec)
    n = length(t_vec)
    orbits = zeros(n)
    for j in 2:n
        orbits[j] = orbits[j-1] + (t_vec[j] - t_vec[j-1]) / (2π * sqrt(a_vec[j-1]^3 / MU))
    end
    return orbits
end

if isfile(cache_path)
    # ── Load from cache ──────────────────────────────────────────────────────
    raw   = readdlm(cache_path, ',', Float64; skipstart=1)
    t     = raw[:, 1]
    a_vec = raw[:, 2]
    e_vec = raw[:, 3]
    i_vec = raw[:, 4]   # radians
    Ω_vec = raw[:, 5]   # radians
    println("Loaded OE cache: $cache_path")
else
    # ── Compute from CSV using elements_time_series (mirrors plot_OE) ────────
    sol, metadata, p = load_timeseries_csv(csv_path)
    N     = p[:N]

    # Use rv2coe_2pi to avoid ±180° branch-cut jumps in difference plots
    elems     = elements_time_series(sol, MU, rv2coe_2pi)
    oe_series = elems[N]   # target satellite
    t_raw     = sol.t

    a_raw = getfield.(oe_series, :a)
    e_raw = getfield.(oe_series, :e)
    i_raw = getfield.(oe_series, :i)   # radians
    Ω_raw = getfield.(oe_series, :Ω)   # radians

    # Interpolate onto a uniform time grid (same as plot_OE)
    t = collect(LinRange(t_raw[1], t_raw[end], N_PTS_OE_TOY))
    function interp1(vals, tq)
        ii = clamp(searchsortedlast(t_raw, tq), 1, length(t_raw)-1)
        α  = (tq - t_raw[ii]) / (t_raw[ii+1] - t_raw[ii])
        vals[ii] * (1-α) + vals[ii+1] * α
    end
    a_vec = [interp1(a_raw, tq) for tq in t]
    e_vec = [interp1(e_raw, tq) for tq in t]
    i_vec = [interp1(i_raw, tq) for tq in t]
    Ω_vec = [interp1(Ω_raw, tq) for tq in t]

    # Save cache (angles in radians, consistent with plot_OE convention)
    open(cache_path, "w") do io
        writedlm(io, ["t" "a_m" "e" "i_rad" "Omega_rad"], ',')
        writedlm(io, hcat(t, a_vec, e_vec, i_vec, Ω_vec), ',')
    end
    println("Saved OE cache: $cache_path")
end

# ── Compute orbit-count axis (same method as plot_OE) ────────────────────────
oc = compute_orbits_curve(a_vec, t)

# ── Compute Δ from initial value ─────────────────────────────────────────────
Δa = a_vec .- a_vec[1]
Δe = e_vec .- e_vec[1]
Δi = rad2deg.(i_vec .- i_vec[1])   # degrees
ΔΩ = rad2deg.(Ω_vec .- Ω_vec[1])   # degrees

# ── 2×2 subplot: Δa, Δe, Δi, ΔΩ vs orbits ───────────────────────────────────
p_a = plot(oc, Δa, ylabel="Δa (m)",   legend=false, color=:blue,   lw=2, formatter=:scientific)
p_e = plot(oc, Δe, ylabel="Δe",       legend=false, color=:red,    lw=2, formatter=:scientific)
p_i = plot(oc, Δi, ylabel="Δi (deg)", legend=false, color=:green,  lw=2, formatter=:scientific)
p_Ω = plot(oc, ΔΩ, ylabel="ΔΩ (deg)", legend=false, color=:purple, lw=2, formatter=:scientific)

for sp in (p_a, p_e, p_i, p_Ω)
    plot!(sp, xlabel="Orbits")
end

plt = plot(p_a, p_e, p_i, p_Ω,
           layout=(2, 2),
           plot_title="Change in orbital elements — target satellite\n($csv_tag)",
           size=(1000, 700),
           margin=5Plots.mm)

display(plt)

# ── Save ─────────────────────────────────────────────────────────────────────
save_path = joinpath(output_dir, "delta_OE_toy_problem_$(csv_tag).png")
savefig(plt, save_path)
println("Saved figure to: $save_path")
