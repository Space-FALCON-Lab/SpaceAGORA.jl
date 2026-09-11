using LinearAlgebra, Statistics, Plots, Printf, StaticArrays, DelimitedFiles

const MU      = 3.986004418e14
const C       = 3.0e8
const R_EARTH = 6_378_137.0
@inline idx(i, off) = 6*(i-1) + off

include("../../functions/4_Diagnostics.jl")
include("../../functions/12_CSV_Write_Read.jl")

const N_PTS_LM = 2000

"""
    plot_linear_momentum(h, i_deg; nu=nothing, show_variance=true) -> (plt1, plt2)

Plot net satellite momentum and fractional drift for all CSVs matching the
given orbital parameters.

# Arguments
- `h`:             target altitude in km (e.g. 850)
- `i_deg`:         inclination in degrees (e.g. 0.5)
- `nu`:            true anomaly offset in degrees (optional, e.g. -0.75)
- `show_variance`: overlay ±1σ shaded band (default: true)

# Returns
A named tuple with fields:
- `plt1`, `plt2`    — momentum magnitude and fractional drift plots
- `t`               — shared uniform time vector
- `all_Pmag`        — `Vector{Vector{Float64}}` interpolated |P(t)| per CSV
- `all_ΔP`          — `Vector{Vector{Float64}}` interpolated ΔP/P(0) per CSV
- `mean_Pmag`, `std_Pmag` — mean and std of |P(t)|
- `mean_ΔP`,   `std_ΔP`   — mean and std of ΔP/P(0)
"""
function plot_linear_momentum(h::Real, i_deg::Real;
                               nu::Union{Real,Nothing} = nothing,
                               show_variance::Bool = true)
    # --- Build folder name from parameters ---
    folder = if nu === nothing
        @sprintf("target_h%dkm_i%.1fdeg", h, i_deg)
    else
        @sprintf("target_h%dkm_i%.1fdeg_nu%.2fdeg", h, i_deg, nu)
    end

    data_dir   = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", folder))
    dir_tag    = basename(data_dir)
    output_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "images_reconstructed_from_csv", "linear_momentum"))
    mkpath(output_dir)
    cache_path = joinpath(output_dir, "lm_cache_$(dir_tag).csv")

    all_t    = Vector{Vector{Float64}}()
    all_Pmag = Vector{Vector{Float64}}()
    all_ΔP   = Vector{Vector{Float64}}()
    N_vals   = Int[]

    if isfile(cache_path)
        # --- Load from cache ---
        header_line = readline(cache_path)
        col_names   = split(header_line, ',')
        N_vals      = [parse(Int, match(r"Pmag_N(\d+)", string(c))[1])
                       for c in col_names if startswith(string(c), "Pmag")]
        raw         = readdlm(cache_path, ',', Float64; skipstart=1)
        t_vec       = raw[:, 1]
        for k in eachindex(N_vals)
            push!(all_t,    t_vec)
            push!(all_Pmag, raw[:, 2 + 2*(k-1)])
            push!(all_ΔP,   raw[:, 3 + 2*(k-1)])
        end
        println("Loaded linear momentum cache: $cache_path")
    else
        # --- Compute from CSVs ---
        isdir(data_dir) || error("Directory not found: $data_dir")
        csv_files = filter(f -> endswith(f, ".csv"), readdir(data_dir; join=true))

        for csv_path in csv_files
            local sol, _, p = load_timeseries_csv(csv_path)
            local masses    = p[:masses]
            local Pmag      = [total_momentum(u, masses)[2] for u in sol.u]
            local ΔP        = (Pmag .- Pmag[1]) ./ Pmag[1]
            local N         = p[:N]

            local t_uniform = collect(LinRange(sol.t[1], sol.t[end], N_PTS_LM))
            local interp    = (vals, t_query) -> begin
                ii = clamp(searchsortedlast(sol.t, t_query), 1, length(sol.t)-1)
                α  = (t_query - sol.t[ii]) / (sol.t[ii+1] - sol.t[ii])
                vals[ii] * (1-α) + vals[ii+1] * α
            end
            local Pmag_u = [interp(Pmag, t) for t in t_uniform]
            local ΔP_u   = [interp(ΔP,   t) for t in t_uniform]

            push!(all_t,    t_uniform)
            push!(all_Pmag, Pmag_u)
            push!(all_ΔP,   ΔP_u)
            push!(N_vals,   N)
        end

        # Build header and data matrix, then save
        n_csvs = length(N_vals)
        header = Matrix{String}(undef, 1, 1 + 2*n_csvs)
        header[1, 1] = "t"
        for k in 1:n_csvs
            header[1, 2 + 2*(k-1)] = "Pmag_N$(N_vals[k])"
            header[1, 3 + 2*(k-1)] = "dP_N$(N_vals[k])"
        end
        data_mat = hcat(all_t[1],
                        [col for (pm, dp) in zip(all_Pmag, all_ΔP) for col in [pm, dp]]...)
        open(cache_path, "w") do io
            writedlm(io, header,   ',')
            writedlm(io, data_mat, ',')
        end
        println("Saved linear momentum cache: $cache_path")
    end

    # --- Plot ---
    plt1 = plot(title="Net Satellite Momentum ($dir_tag)", xlabel="t (s)",
                ylabel="Total Momentum (kg·m/s)", formatter=:scientific, legend=:outertopright)
    plt2 = plot(title="Linear Momentum Fractional Drift ($dir_tag)", xlabel="t (s)",
                ylabel="ΔP/P(0) = (P(t) - P(0)) / P(0)", formatter=:scientific, legend=:outertopright)

    for (Pmag_u, ΔP_u, N) in zip(all_Pmag, all_ΔP, N_vals)
        plot!(plt1, all_t[1], Pmag_u, label="N=$N", color=:grey, alpha=0.5, lw=1)
        plot!(plt2, all_t[1], ΔP_u,   label="N=$N", color=:grey, alpha=0.5, lw=1)
    end

    mean_Pmag = mean(stack(all_Pmag), dims=2) |> vec
    mean_ΔP   = mean(stack(all_ΔP),   dims=2) |> vec
    std_Pmag  = std(stack(all_Pmag),  dims=2) |> vec
    std_ΔP    = std(stack(all_ΔP),    dims=2) |> vec

    if show_variance
        plot!(plt1, all_t[1], mean_Pmag .+ std_Pmag, fillrange=mean_Pmag .- std_Pmag,
              fillalpha=0.2, fillcolor=:blue, linealpha=0, label="±1σ")
        plot!(plt2, all_t[1], mean_ΔP   .+ std_ΔP,   fillrange=mean_ΔP   .- std_ΔP,
              fillalpha=0.2, fillcolor=:blue, linealpha=0, label="±1σ")
    end
    plot!(plt1, all_t[1], mean_Pmag, label="Mean", color=:blue, lw=2)
    plot!(plt2, all_t[1], mean_ΔP,   label="Mean", color=:blue, lw=2)

    return (plt1=plt1, plt2=plt2,
            t=all_t[1],
            all_Pmag=all_Pmag, all_ΔP=all_ΔP,
            mean_Pmag=mean_Pmag, std_Pmag=std_Pmag,
            mean_ΔP=mean_ΔP,   std_ΔP=std_ΔP)
end


