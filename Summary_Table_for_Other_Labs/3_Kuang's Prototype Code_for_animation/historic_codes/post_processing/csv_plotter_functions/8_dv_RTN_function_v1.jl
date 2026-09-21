# this uses toy problem scheduling technique
using LinearAlgebra, Statistics, Plots, Printf, StaticArrays, DelimitedFiles

if !@isdefined(MU);      const MU      = 3.986004418e14; end
if !@isdefined(C);       const C       = 3.0e8;          end
if !@isdefined(R_EARTH); const R_EARTH = 6_378_137.0;    end
# Unconditional const — re-include just emits a harmless warning (same value).
const R_ATMDEF = 6_478_137.0  # R_EARTH + Kármán line (100 km)
@inline idx(i, off) = 6*(i-1) + off

include("../../functions/1_LOS_Metrics.jl")
include("../../functions/2_Laser_Forces_ver3.jl")
include("../../functions/4_Diagnostics.jl")
include("../../functions/5_OE_Converters.jl")
include("../../functions/6_OE_and_dv_in_RTN.jl")
include("../../functions/12_CSV_Write_Read.jl")

const N_PTS_DV = 2000 # number of points to interpolate Δv time series to for uniformity across CSVs

"""
    plot_dv_RTN(h, i_deg; component=:T, nu=nothing, show_variance=true) -> NamedTuple

Compute and plot the cumulative Δv time series in a chosen RTN component for
the target satellite, across all CSVs matching the given orbital parameters.

# Arguments
- `h`:             target altitude in km (e.g. 1000)
- `i_deg`:         inclination in degrees (e.g. 0.0)
- `component`:     RTN component to plot — `:R`, `:T`, or `:N` (default `:T`)
- `nu`:            true anomaly offset in degrees (optional, e.g. -0.75)
- `show_variance`: overlay ±1σ shaded band (default: true)

# Returns
A named tuple with fields:
- `plt`           — Δv time series plot
- `t`             — shared uniform time vector
- `all_dv`        — `Vector{Vector{Float64}}` interpolated Δv(t) per CSV
- `mean_dv`, `std_dv` — mean and std of Δv(t)
"""
function plot_dv_RTN(h::Real, i_deg::Real;
                     component::Symbol = :T,
                     nu::Union{Real,Nothing} = nothing,
                     show_variance::Bool = true)
    component in (:R, :T, :N) || error("component must be :R, :T, or :N")
    comp_idx = component == :R ? 1 : component == :T ? 2 : 3
    comp_str = string(component)

    # --- Build folder name from parameters ---
    folder = if nu === nothing
        @sprintf("target_h%dkm_i%.1fdeg", h, i_deg)
    else
        @sprintf("target_h%dkm_i%.1fdeg_nu%.2fdeg", h, i_deg, nu)
    end

    data_dir   = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", folder))
    dir_tag    = basename(data_dir)
    output_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "images_reconstructed_from_csv", "delta_v_RTN"))
    mkpath(output_dir)

    cache_path = joinpath(output_dir, "dv_RTN_cache_$(dir_tag).csv")

    all_t    = Vector{Vector{Float64}}()
    all_dv_R = Vector{Vector{Float64}}()
    all_dv_T = Vector{Vector{Float64}}()
    all_dv_N_comp = Vector{Vector{Float64}}()
    N_vals   = Int[]

    if isfile(cache_path)
        # --- Load all three components from cache ---
        header_line = readline(cache_path)
        col_names   = split(header_line, ',')
        N_vals      = [parse(Int, match(r"dv_R_N(\d+)", string(c))[1])
                       for c in col_names if occursin("dv_R_N", string(c))]
        raw    = readdlm(cache_path, ',', Float64; skipstart=1)
        t_vec  = raw[:, 1]
        for k in eachindex(N_vals)
            push!(all_t,       t_vec)
            push!(all_dv_R,    raw[:, 2 + 3*(k-1)])
            push!(all_dv_T,    raw[:, 3 + 3*(k-1)])
            push!(all_dv_N_comp, raw[:, 4 + 3*(k-1)])
        end
        println("Loaded Δv RTN cache: $cache_path")
    else
        # --- Compute all three components from CSVs and cache them ---
        isdir(data_dir) || error("Directory not found: $data_dir")
        csv_files = filter(f -> endswith(f, ".csv"), readdir(data_dir; join=true))

        for csv_path in csv_files
            local sol, _, p  = load_timeseries_csv(csv_path)
            local t, Δv_hist = delta_v_RTN_time_series(sol, p)
            local target_dv  = Δv_hist[p[:N]]   # 3×length(t), target = last satellite

            local t_uniform = collect(LinRange(t[1], t[end], N_PTS_DV))
            local interp = (row, tq) -> begin
                ii = clamp(searchsortedlast(t, tq), 1, length(t)-1)
                α  = (tq - t[ii]) / (t[ii+1] - t[ii])
                target_dv[row, ii] * (1-α) + target_dv[row, ii+1] * α
            end
            push!(all_t,         t_uniform)
            push!(all_dv_R,      [interp(1, tq) for tq in t_uniform])
            push!(all_dv_T,      [interp(2, tq) for tq in t_uniform])
            push!(all_dv_N_comp, [interp(3, tq) for tq in t_uniform])
            push!(N_vals,        p[:N])
        end

        isempty(N_vals) && error("No CSVs loaded from: $data_dir")

        # Save cache: columns t, dv_R_N{N}, dv_T_N{N}, dv_N_N{N} per CSV
        n_csvs = length(N_vals)
        header = Matrix{String}(undef, 1, 1 + 3*n_csvs)
        header[1, 1] = "t"
        for k in 1:n_csvs
            header[1, 2 + 3*(k-1)] = "dv_R_N$(N_vals[k])"
            header[1, 3 + 3*(k-1)] = "dv_T_N$(N_vals[k])"
            header[1, 4 + 3*(k-1)] = "dv_N_N$(N_vals[k])"
        end
        data_mat = hcat(all_t[1],
                        [col for (r, tv, nv) in zip(all_dv_R, all_dv_T, all_dv_N_comp)
                             for col in [r, tv, nv]]...)
        open(cache_path, "w") do io
            writedlm(io, header,   ',')
            writedlm(io, data_mat, ',')
        end
        println("Saved Δv RTN cache: $cache_path")
    end

    # Pick the requested component
    all_dv = comp_idx == 1 ? all_dv_R : comp_idx == 2 ? all_dv_T : all_dv_N_comp

    t_common = all_t[1]
    mean_dv  = mean(stack(all_dv), dims=2) |> vec
    std_dv   = std(stack(all_dv),  dims=2) |> vec

    # --- Plot ---
    plt = plot(title="Cumulative Δv_$comp_str — target satellite ($dir_tag)",
               xlabel="t (s)", ylabel="Δv_$comp_str (m/s)",
               formatter=:scientific, legend=:outertopright)

    for (dv_u, N) in zip(all_dv, N_vals)
        plot!(plt, t_common, dv_u, label="N=$N", color=:grey, alpha=0.5, lw=1)
    end

    if show_variance
        plot!(plt, t_common, mean_dv .+ std_dv,
              fillrange=mean_dv .- std_dv,
              fillalpha=0.2, fillcolor=:blue, linealpha=0, label="±1σ")
    end
    plot!(plt, t_common, mean_dv, label="Mean", color=:blue, lw=2)

    return (plt=plt,
            t=t_common,
            all_dv=all_dv,
            mean_dv=mean_dv,
            std_dv=std_dv)
end


