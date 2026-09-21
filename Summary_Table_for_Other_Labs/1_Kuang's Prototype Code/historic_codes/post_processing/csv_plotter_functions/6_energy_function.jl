using LinearAlgebra, Statistics, Plots, Printf, StaticArrays, DelimitedFiles

if !@isdefined(MU);      const MU      = 3.986004418e14; end
if !@isdefined(C);       const C       = 3.0e8;          end
if !@isdefined(R_EARTH); const R_EARTH = 6_378_137.0;    end
@inline idx(i, off) = 6*(i-1) + off

include("../../functions/4_Diagnostics.jl")
include("../../functions/12_CSV_Write_Read.jl")

const N_PTS_E = 2000

"""
    plot_energy(h, i_deg; nu=nothing, show_variance=true) -> NamedTuple

Plot satellite orbital energy variation and fractional drift for all CSVs
matching the given orbital parameters. Results are cached to a CSV for fast
re-runs.

# Arguments
- `h`:             target altitude in km (e.g. 1000)
- `i_deg`:         inclination in degrees (e.g. 0.0)
- `nu`:            true anomaly offset in degrees (optional, e.g. -0.75)
- `show_variance`: overlay ±1σ shaded band (default: true)

# Returns
A named tuple with fields:
- `plt1`, `plt2`        — energy variation and fractional drift plots
- `t`                   — shared uniform time vector
- `all_ΔE`              — `Vector{Vector{Float64}}` interpolated ΔE(t) per CSV
- `all_ΔE_fr`           — `Vector{Vector{Float64}}` interpolated ΔE/|E₀| per CSV
- `mean_ΔE`, `std_ΔE`   — mean and std of ΔE(t)
- `mean_ΔE_fr`, `std_ΔE_fr` — mean and std of ΔE/|E₀|
"""
function plot_energy(h::Real, i_deg::Real;
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
    output_dir = normpath(joinpath(@__DIR__, "..", "..", "output", "images_reconstructed_from_csv", "orbital_energy"))
    mkpath(output_dir)
    cache_path = joinpath(output_dir, "e_cache_$(dir_tag).csv")

    all_t     = Vector{Vector{Float64}}()
    all_ΔE    = Vector{Vector{Float64}}()
    all_ΔE_fr = Vector{Vector{Float64}}()
    N_vals    = Int[]

    if isfile(cache_path)
        # --- Load from cache ---
        header_line = readline(cache_path)
        col_names   = split(header_line, ',')
        N_vals      = [parse(Int, match(r"dE_N(\d+)", string(c))[1])
                       for c in col_names if startswith(string(c), "dE_N")]
        raw         = readdlm(cache_path, ',', Float64; skipstart=1)
        t_vec       = raw[:, 1]
        for k in eachindex(N_vals)
            push!(all_t,     t_vec)
            push!(all_ΔE,    raw[:, 2 + 2*(k-1)])
            push!(all_ΔE_fr, raw[:, 3 + 2*(k-1)])
        end
        println("Loaded energy cache: $cache_path")
    else
        # --- Compute from CSVs ---
        isdir(data_dir) || error("Directory not found: $data_dir")
        csv_files = filter(f -> endswith(f, ".csv"), readdir(data_dir; join=true))

        for csv_path in csv_files
            local sol, _, p = load_timeseries_csv(csv_path)
            local masses    = p[:masses]
            local Etot      = [sum(orbital_energy(u, masses, MU)) for u in sol.u]
            local ΔE        = Etot .- Etot[1]
            local ΔE_frac   = ΔE ./ abs(Etot[1])
            local N         = p[:N]

            local t_uniform = collect(LinRange(sol.t[1], sol.t[end], N_PTS_E))
            local interp    = (vals, t_query) -> begin
                ii = clamp(searchsortedlast(sol.t, t_query), 1, length(sol.t)-1)
                α  = (t_query - sol.t[ii]) / (sol.t[ii+1] - sol.t[ii])
                vals[ii] * (1-α) + vals[ii+1] * α
            end
            local ΔE_u      = [interp(ΔE,      t) for t in t_uniform]
            local ΔE_frac_u = [interp(ΔE_frac, t) for t in t_uniform]

            push!(all_t,     t_uniform)
            push!(all_ΔE,    ΔE_u)
            push!(all_ΔE_fr, ΔE_frac_u)
            push!(N_vals,    N)
        end

        # Build header and data matrix, then save
        n_csvs = length(N_vals)
        header = Matrix{String}(undef, 1, 1 + 2*n_csvs)
        header[1, 1] = "t"
        for k in 1:n_csvs
            header[1, 2 + 2*(k-1)] = "dE_N$(N_vals[k])"
            header[1, 3 + 2*(k-1)] = "dE_fr_N$(N_vals[k])"
        end
        data_mat = hcat(all_t[1],
                        [col for (de, de_fr) in zip(all_ΔE, all_ΔE_fr) for col in [de, de_fr]]...)
        open(cache_path, "w") do io
            writedlm(io, header,   ',')
            writedlm(io, data_mat, ',')
        end
        println("Saved energy cache: $cache_path")
    end

    # --- Plot ---
    plt1 = plot(title="Satellite Orbital Energy Variation ($dir_tag)", xlabel="t (s)",
                ylabel="ΔE = E(t) - E(0) (J)", formatter=:scientific, legend=:outertopright)
    plt2 = plot(title="Orbital Energy Fractional Drift ($dir_tag)", xlabel="t (s)",
                ylabel="ΔE/|E(0)| = (E(t) - E(0)) / |E(0)|", formatter=:scientific, legend=:outertopright)

    for (ΔE_u, ΔE_frac_u, N) in zip(all_ΔE, all_ΔE_fr, N_vals)
        plot!(plt1, all_t[1], ΔE_u,      label="N=$N", color=:grey, alpha=0.5, lw=1)
        plot!(plt2, all_t[1], ΔE_frac_u, label="N=$N", color=:grey, alpha=0.5, lw=1)
    end

    mean_ΔE    = mean(stack(all_ΔE),    dims=2) |> vec
    mean_ΔE_fr = mean(stack(all_ΔE_fr), dims=2) |> vec
    std_ΔE     = std(stack(all_ΔE),     dims=2) |> vec
    std_ΔE_fr  = std(stack(all_ΔE_fr),  dims=2) |> vec

    if show_variance
        plot!(plt1, all_t[1], mean_ΔE    .+ std_ΔE,    fillrange=mean_ΔE    .- std_ΔE,
              fillalpha=0.2, fillcolor=:blue, linealpha=0, label="±1σ")
        plot!(plt2, all_t[1], mean_ΔE_fr .+ std_ΔE_fr, fillrange=mean_ΔE_fr .- std_ΔE_fr,
              fillalpha=0.2, fillcolor=:blue, linealpha=0, label="±1σ")
    end
    plot!(plt1, all_t[1], mean_ΔE,    label="Mean", color=:blue, lw=2)
    plot!(plt2, all_t[1], mean_ΔE_fr, label="Mean", color=:blue, lw=2)

    return (plt1=plt1, plt2=plt2,
            t=all_t[1],
            all_ΔE=all_ΔE, all_ΔE_fr=all_ΔE_fr,
            mean_ΔE=mean_ΔE,    std_ΔE=std_ΔE,
            mean_ΔE_fr=mean_ΔE_fr, std_ΔE_fr=std_ΔE_fr)
end


