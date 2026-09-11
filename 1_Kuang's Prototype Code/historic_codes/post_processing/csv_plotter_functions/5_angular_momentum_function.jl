using LinearAlgebra, Statistics, Plots, Printf, StaticArrays, DelimitedFiles

if !@isdefined(MU);      const MU      = 3.986004418e14; end
if !@isdefined(C);       const C       = 3.0e8;          end
if !@isdefined(R_EARTH); const R_EARTH = 6_378_137.0;    end
@inline idx(i, off) = 6*(i-1) + off

include("../../functions/4_Diagnostics.jl")
include("../../functions/12_CSV_Write_Read.jl")

const N_PTS_AM = 2000

"""
    plot_angular_momentum(h, i_deg; nu=nothing, show_variance=true) -> NamedTuple

Plot net satellite angular momentum and fractional drift for all CSVs matching
the given orbital parameters. Results are cached to a CSV for fast re-runs.

# Arguments
- `h`:             target altitude in km (e.g. 850)
- `i_deg`:         inclination in degrees (e.g. 0.5)
- `nu`:            true anomaly offset in degrees (optional, e.g. -0.75)
- `show_variance`: overlay ±1σ shaded band (default: true)

# Returns
A named tuple with fields:
- `plt1`, `plt2`    — momentum magnitude and fractional drift plots
- `t`               — shared uniform time vector
- `all_Hmag`        — `Vector{Vector{Float64}}` interpolated |H(t)| per CSV
- `all_ΔH`          — `Vector{Vector{Float64}}` interpolated ΔH/H(0) per CSV
- `mean_Hmag`, `std_Hmag` — mean and std of |H(t)|
- `mean_ΔH`,   `std_ΔH`   — mean and std of ΔH/H(0)
"""
function plot_angular_momentum(h::Real, i_deg::Real;
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
    output_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "images_reconstructed_from_csv", "angular_momentum"))
    mkpath(output_dir)
    cache_path = joinpath(output_dir, "am_cache_$(dir_tag).csv")

    all_t    = Vector{Vector{Float64}}()
    all_Hmag = Vector{Vector{Float64}}()
    all_ΔH   = Vector{Vector{Float64}}()
    N_vals   = Int[]

    if isfile(cache_path)
        # --- Load from cache ---
        header_line = readline(cache_path)
        col_names   = split(header_line, ',')
        N_vals      = [parse(Int, match(r"Hmag_N(\d+)", string(c))[1])
                       for c in col_names if startswith(string(c), "Hmag")]
        raw         = readdlm(cache_path, ',', Float64; skipstart=1)
        t_vec       = raw[:, 1]
        for k in eachindex(N_vals)
            push!(all_t,    t_vec)
            push!(all_Hmag, raw[:, 2 + 2*(k-1)])
            push!(all_ΔH,   raw[:, 3 + 2*(k-1)])
        end
        println("Loaded angular momentum cache: $cache_path")
    else
        # --- Compute from CSVs ---
        isdir(data_dir) || error("Directory not found: $data_dir")
        csv_files = filter(f -> endswith(f, ".csv"), readdir(data_dir; join=true))

        for csv_path in csv_files
            local sol, _, p = load_timeseries_csv(csv_path)
            local masses    = p[:masses]
            local Hmag      = [angular_momentum(u, masses)[2] for u in sol.u]
            local ΔH        = (Hmag .- Hmag[1]) ./ Hmag[1]
            local N         = p[:N]

            local t_uniform = collect(LinRange(sol.t[1], sol.t[end], N_PTS_AM))
            local interp    = (vals, t_query) -> begin
                ii = clamp(searchsortedlast(sol.t, t_query), 1, length(sol.t)-1)
                α  = (t_query - sol.t[ii]) / (sol.t[ii+1] - sol.t[ii])
                vals[ii] * (1-α) + vals[ii+1] * α
            end
            local Hmag_u = [interp(Hmag, t) for t in t_uniform]
            local ΔH_u   = [interp(ΔH,   t) for t in t_uniform]

            push!(all_t,    t_uniform)
            push!(all_Hmag, Hmag_u)
            push!(all_ΔH,   ΔH_u)
            push!(N_vals,   N)
        end

        # Build header and data matrix, then save
        n_csvs = length(N_vals)
        header = Matrix{String}(undef, 1, 1 + 2*n_csvs)
        header[1, 1] = "t"
        for k in 1:n_csvs
            header[1, 2 + 2*(k-1)] = "Hmag_N$(N_vals[k])"
            header[1, 3 + 2*(k-1)] = "dH_N$(N_vals[k])"
        end
        data_mat = hcat(all_t[1],
                        [col for (hm, dh) in zip(all_Hmag, all_ΔH) for col in [hm, dh]]...)
        open(cache_path, "w") do io
            writedlm(io, header,   ',')
            writedlm(io, data_mat, ',')
        end
        println("Saved angular momentum cache: $cache_path")
    end

    # --- Plot ---
    plt1 = plot(title="Net Satellite Angular Momentum ($dir_tag)", xlabel="t (s)",
                ylabel="Total Angular Momentum (kg·m²/s)", formatter=:scientific, legend=:outertopright)
    plt2 = plot(title="Angular Momentum Fractional Drift ($dir_tag)", xlabel="t (s)",
                ylabel="ΔH/H₀ = (H(t) - H(0)) / H(0)", formatter=:scientific, legend=:outertopright)

    for (Hmag_u, ΔH_u, N) in zip(all_Hmag, all_ΔH, N_vals)
        plot!(plt1, all_t[1], Hmag_u, label="N=$N", color=:grey, alpha=0.5, lw=1)
        plot!(plt2, all_t[1], ΔH_u,   label="N=$N", color=:grey, alpha=0.5, lw=1)
    end

    mean_Hmag = mean(stack(all_Hmag), dims=2) |> vec
    mean_ΔH   = mean(stack(all_ΔH),   dims=2) |> vec
    std_Hmag  = std(stack(all_Hmag),  dims=2) |> vec
    std_ΔH    = std(stack(all_ΔH),    dims=2) |> vec

    if show_variance
        plot!(plt1, all_t[1], mean_Hmag .+ std_Hmag, fillrange=mean_Hmag .- std_Hmag,
              fillalpha=0.2, fillcolor=:blue, linealpha=0, label="±1σ")
        plot!(plt2, all_t[1], mean_ΔH   .+ std_ΔH,   fillrange=mean_ΔH   .- std_ΔH,
              fillalpha=0.2, fillcolor=:blue, linealpha=0, label="±1σ")
    end
    plot!(plt1, all_t[1], mean_Hmag, label="Mean", color=:blue, lw=2)
    plot!(plt2, all_t[1], mean_ΔH,   label="Mean", color=:blue, lw=2)

    return (plt1=plt1, plt2=plt2,
            t=all_t[1],
            all_Hmag=all_Hmag, all_ΔH=all_ΔH,
            mean_Hmag=mean_Hmag, std_Hmag=std_Hmag,
            mean_ΔH=mean_ΔH,   std_ΔH=std_ΔH)
end


