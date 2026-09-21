using LinearAlgebra, Statistics, Plots, Printf, StaticArrays, DelimitedFiles

if !@isdefined(MU);      const MU      = 3.986004418e14; end
if !@isdefined(C);       const C       = 3.0e8;          end
if !@isdefined(R_EARTH); const R_EARTH = 6_378_137.0;    end
@inline idx(i, off) = 6*(i-1) + off

include("../../functions/4_Diagnostics.jl")
include("../../functions/5_OE_Converters.jl")
include("../../functions/6_OE_and_dv_in_RTN.jl")
include("../../functions/12_CSV_Write_Read.jl")

const N_PTS_OE = 2000

"""
    plot_OE(h, i_deg; element=:a, sat=-1, nu=nothing, show_variance=true) -> NamedTuple

Compute and plot the time series of a chosen orbital element for the target
satellite, across all CSVs matching the given orbital parameters. All six
elements are cached together so switching element requires no recomputation.

# Arguments
- `h`:             target altitude in km (e.g. 1000)
- `i_deg`:         inclination in degrees (e.g. 0.0)
- `element`:       which OE to plot — `:a`, `:e`, `:i`, `:Ω`, `:ω`, `:ν`
                   (default `:a`). For circular orbits `:ω` falls back to
                   argument of latitude `u`.
- `sat`:           satellite index (default -1 → last satellite = target)
- `nu`:            true anomaly offset in degrees (optional, e.g. -0.75)
- `show_variance`: overlay ±1σ shaded band (default: true)

# Returns
A named tuple with fields:
- `plt`            — orbital element time series plot
- `t`              — shared uniform time vector (seconds)
- `all_vals`       — `Vector{Vector{Float64}}` interpolated OE(t) per CSV
- `mean_vals`, `std_vals` — mean and std of OE(t)
"""
function plot_OE(h::Real, i_deg::Real;
                 element::Symbol = :a,
                 sat::Int = -1,
                 nu::Union{Real,Nothing} = nothing,
                 show_variance::Bool = true,
                 filter_helpers::Union{Vector{Int},Nothing} = nothing)
    element in (:a, :e, :i, :Ω, :ω, :ν) ||
        error("element must be one of :a, :e, :i, :Ω, :ω, :ν")

    # ── folder / path setup ───────────────────────────────────────────────────
    folder = if nu === nothing
        @sprintf("target_h%dkm_i%.1fdeg", h, i_deg)
    else
        @sprintf("target_h%dkm_i%.1fdeg_nu%.2fdeg", h, i_deg, nu)
    end

    data_dir   = normpath(joinpath(@__DIR__, "..", "..", "output", "CSV", folder))
    dir_tag    = basename(data_dir)
    output_dir = normpath(joinpath(@__DIR__, "..", "..", "output",
                                   "images_reconstructed_from_csv", "orbital_elements"))
    mkpath(output_dir)
    cache_path = joinpath(output_dir, "OE_cache_$(dir_tag).csv")

    # storage — all 6 OE cached together
    all_t  = Vector{Vector{Float64}}()
    all_a  = Vector{Vector{Float64}}()
    all_e  = Vector{Vector{Float64}}()
    all_i  = Vector{Vector{Float64}}()
    all_Ω  = Vector{Vector{Float64}}()
    all_ω  = Vector{Vector{Float64}}()   # ω (non-circular) or u (circular)
    all_ν  = Vector{Vector{Float64}}()
    N_vals = Int[]

    if isfile(cache_path)
        # ── load from cache ───────────────────────────────────────────────────
        header_line = readline(cache_path)
        col_names   = split(header_line, ',')
        N_vals = [parse(Int, match(r"a_N(\d+)", string(c))[1])
                  for c in col_names if occursin(r"^a_N\d+$", string(c))]
        raw   = readdlm(cache_path, ',', Float64; skipstart=1)
        t_vec = raw[:, 1]
        for k in eachindex(N_vals)
            base = 2 + 6*(k-1)
            push!(all_t, t_vec)
            push!(all_a, raw[:, base+0])
            push!(all_e, raw[:, base+1])
            push!(all_i, raw[:, base+2])
            push!(all_Ω, raw[:, base+3])
            push!(all_ω, raw[:, base+4])
            push!(all_ν, raw[:, base+5])
        end
        println("Loaded OE cache: $cache_path")

        # ── filter to requested helper counts (N_sat = helpers + 1) ─────────
        if filter_helpers !== nothing
            keep   = [k for k in eachindex(N_vals) if (N_vals[k] - 1) in filter_helpers]
            all_t  = all_t[keep];  all_a = all_a[keep];  all_e = all_e[keep]
            all_i  = all_i[keep];  all_Ω = all_Ω[keep];  all_ω = all_ω[keep]
            all_ν  = all_ν[keep];  N_vals = N_vals[keep]
            isempty(N_vals) && error("No cached runs match filter_helpers=$filter_helpers")
        end
    else
        # ── compute from CSVs ─────────────────────────────────────────────────
        isdir(data_dir) || error("Directory not found: $data_dir")
        csv_files = filter(f -> endswith(f, ".csv"), readdir(data_dir; join=true))

        for csv_path in csv_files
            local sol, _, p = load_timeseries_csv(csv_path)
            local N         = p[:N]
            local sat_idx   = sat == -1 ? N : sat
            local elems     = elements_time_series(sol, MU)
            local oe_series = elems[sat_idx]
            local t         = sol.t

            # extract raw element vectors
            local a_raw = getfield.(oe_series, :a)
            local e_raw = getfield.(oe_series, :e)
            local i_raw = getfield.(oe_series, :i)
            local Ω_raw = getfield.(oe_series, :Ω)

            # ω: use argument of latitude u for circular orbits
            local first_oe = first(oe_series)
            local ω_raw = if hasproperty(first_oe, :ω)
                raw_ω   = getfield.(oe_series, :ω)
                raw_u   = hasproperty(first_oe, :u) ? getfield.(oe_series, :u) :
                          fill(NaN, length(oe_series))
                raw_ν_m = hasproperty(first_oe, :ν) ? getfield.(oe_series, :ν) :
                          fill(NaN, length(oe_series))
                mask_circ = .!(e_raw .> 1e-8) .| isnan.(raw_ν_m)
                [mask_circ[k] ? raw_u[k] : raw_ω[k] for k in eachindex(t)]
            elseif hasproperty(first_oe, :u)
                getfield.(oe_series, :u)
            else
                fill(NaN, length(t))
            end

            local ν_raw = hasproperty(first_oe, :ν) ?
                          getfield.(oe_series, :ν) : fill(NaN, length(oe_series))

            # interpolate all elements to a uniform time grid
            local t_uniform = collect(LinRange(t[1], t[end], N_PTS_OE))
            function interp1(vals, tq)
                ii = clamp(searchsortedlast(t, tq), 1, length(t)-1)
                α  = (tq - t[ii]) / (t[ii+1] - t[ii])
                vals[ii] * (1-α) + vals[ii+1] * α
            end

            push!(all_t, t_uniform)
            push!(all_a, [interp1(a_raw, tq) for tq in t_uniform])
            push!(all_e, [interp1(e_raw, tq) for tq in t_uniform])
            push!(all_i, [interp1(i_raw, tq) for tq in t_uniform])
            push!(all_Ω, [interp1(Ω_raw, tq) for tq in t_uniform])
            push!(all_ω, [interp1(ω_raw, tq) for tq in t_uniform])
            push!(all_ν, [interp1(ν_raw, tq) for tq in t_uniform])
            push!(N_vals, N)
        end

        isempty(N_vals) && error("No CSVs loaded from: $data_dir")

        # ── filter to requested helper counts (N_sat = helpers + 1) ─────────
        if filter_helpers !== nothing
            keep   = [k for k in eachindex(N_vals) if (N_vals[k] - 1) in filter_helpers]
            all_t  = all_t[keep];  all_a = all_a[keep];  all_e = all_e[keep]
            all_i  = all_i[keep];  all_Ω = all_Ω[keep];  all_ω = all_ω[keep]
            all_ν  = all_ν[keep];  N_vals = N_vals[keep]
            isempty(N_vals) && error("No CSVs match filter_helpers=$filter_helpers")
        end

        # ── save cache: t, a/e/i/Ω/ω/ν per CSV ──────────────────────────────
        n_csvs = length(N_vals)
        header = Matrix{String}(undef, 1, 1 + 6*n_csvs)
        header[1, 1] = "t"
        for k in 1:n_csvs
            base = 2 + 6*(k-1)
            header[1, base+0] = "a_N$(N_vals[k])"
            header[1, base+1] = "e_N$(N_vals[k])"
            header[1, base+2] = "i_N$(N_vals[k])"
            header[1, base+3] = "Omega_N$(N_vals[k])"
            header[1, base+4] = "omega_N$(N_vals[k])"
            header[1, base+5] = "nu_N$(N_vals[k])"
        end
        data_mat = hcat(all_t[1],
                        [col for (va, ve, vi, vΩ, vω, vν)
                             in zip(all_a, all_e, all_i, all_Ω, all_ω, all_ν)
                             for col in [va, ve, vi, vΩ, vω, vν]]...)
        open(cache_path, "w") do io
            writedlm(io, header,   ',')
            writedlm(io, data_mat, ',')
        end
        println("Saved OE cache: $cache_path")
    end

    # ── select requested element (convert angles rad → deg) ──────────────────
    rad2deg_vecs(vecs) = [v .* (180/π) for v in vecs]
    all_vals, ylabel_str, elem_label = if element == :a
        all_a, "a, m",       "Semi-major axis"
    elseif element == :e
        all_e, "e",            "Eccentricity"
    elseif element == :i
        rad2deg_vecs(all_i), "i, deg",     "Inclination"
    elseif element == :Ω
        rad2deg_vecs(all_Ω), "Ω, deg",     "RAAN"
    elseif element == :ω
        rad2deg_vecs(all_ω), "ω or u, deg", "ω (non-circ) / u (circ)"
    else  # :ν
        rad2deg_vecs(all_ν), "ν, deg",     "True anomaly"
    end

    t_common  = all_t[1]
    t_hr      = t_common ./ 3600.0
    mean_vals = mean(stack(all_vals), dims=2) |> vec
    std_vals  = std(stack(all_vals),  dims=2) |> vec

    # Compute integrated orbit count x-axis
    # Use instantaneous a(t) for each curve, then average for mean/std
    a_curves = element == :a ? all_vals : all_a
    # Ensure all orbits_curves and orbits_mean have length N_PTS_OE (same as t_common)
    function compute_orbits_curve(a_vec, t_vec)
        n = length(t_vec)
        orbits = zeros(n)
        for j in 2:n
            dt = t_vec[j] - t_vec[j-1]
            a = a_vec[j-1]
            orbits[j] = orbits[j-1] + dt / (2π * sqrt(a^3 / MU))
        end
        return orbits
    end
    orbits_curves = [compute_orbits_curve(ac, t_common) for ac in a_curves]
    # For mean/std, use mean a(t)
    mean_a = mean(stack(all_a), dims=2) |> vec
    orbits_mean = compute_orbits_curve(mean_a, t_common)

    # ── plot ──────────────────────────────────────────────────────────────────
    plt = plot(title="$elem_label — target satellite ($dir_tag)",
               xlabel="Integrated orbits", ylabel=ylabel_str,
               formatter=:scientific, legend=:outertopright)

    for (vals, N, orbits) in zip(all_vals, N_vals, orbits_curves)
        plot!(plt, orbits, vals, label="N=$N", color=:grey, alpha=0.5, lw=1)
    end

    if show_variance
        plot!(plt, orbits_mean, mean_vals .+ std_vals,
              fillrange=mean_vals .- std_vals,
              fillalpha=0.2, fillcolor=:blue, linealpha=0, label="±1σ")
    end
    plot!(plt, orbits_mean, mean_vals, label="Mean", color=:blue, lw=2)

    return (plt=plt,
            t=t_common,
            orbits=orbits_mean,
            N_vals=N_vals,
            all_vals=all_vals,
            mean_vals=mean_vals,
            std_vals=std_vals)
end
