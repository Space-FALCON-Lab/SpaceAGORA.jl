using LinearAlgebra, Statistics, Plots, Printf, StaticArrays, DelimitedFiles

if !@isdefined(MU);      const MU      = 3.986004418e14; end
if !@isdefined(C);       const C       = 3.0e8;          end
if !@isdefined(R_EARTH); const R_EARTH = 6_378_137.0;    end
if !@isdefined(R_ATMDEF)
    const R_ATMDEF = 6_478_137.0  # R_EARTH + Kármán line (100 km)
end
@inline idx(i, off) = 6*(i-1) + off

include("../../functions/1_LOS_Metrics.jl")
include("../../functions/2_Laser_Forces_ver2.jl")
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
    all_orbits = Vector{Vector{Float64}}()  # orbit count arrays per CSV
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
        # Recompute orbit count arrays from sol if needed (not cached)
        # For cache loads, just use nominal period for now (for compatibility)
        for tvec in all_t
            a_m = R_EARTH + h * 1e3
            T_orbit = 2π * sqrt(a_m^3 / MU)
            push!(all_orbits, tvec ./ T_orbit)
        end
        println("Loaded Δv RTN cache: $cache_path")
    else
        # --- Compute all three components from CSVs and cache them ---
        isdir(data_dir) || error("Directory not found: $data_dir")
        csv_files = filter(f -> endswith(f, ".csv"), readdir(data_dir; join=true))

        max_orbits_all = Float64[]
        _loaded = []  # (t, dv_R, dv_T, dv_N, oc) per CSV
        for csv_path in csv_files
            local sol, _, p  = load_timeseries_csv(csv_path)
            local t, Δv_hist = delta_v_RTN_time_series(sol, p)
            local target_dv  = Δv_hist[p[:N]]   # 3×length(t), target = last satellite

            # Compute orbit count array for this CSV using instantaneous a(t)
            oc = Vector{Float64}(undef, length(t))
            oc[1] = 0.0
            for k in 2:length(t)
                u = sol.u[k]
                r = @SVector [u[idx(p[:N],1)], u[idx(p[:N],2)], u[idx(p[:N],3)]]
                v = @SVector [u[idx(p[:N],4)], u[idx(p[:N],5)], u[idx(p[:N],6)]]
                a = rv2coe(r, v, MU).a
                T = 2π * sqrt(a^3 / MU)
                dt = t[k] - t[k-1]
                oc[k] = oc[k-1] + dt / T
            end

            # Interpolate dv and orbit count onto a uniform orbit grid
            oc_max = oc[end]
            push!(max_orbits_all, oc_max)
            push!(_loaded, (t=t, dv_R=target_dv[1,:], dv_T=target_dv[2,:], dv_N=target_dv[3,:], oc=oc))
        end

        isempty(_loaded) && error("No CSVs loaded from: $data_dir")

        # Common orbit-count grid: 0 → min of all max orbit counts (no extrapolation)
        orbit_common = collect(LinRange(0.0, minimum(max_orbits_all), N_PTS_DV))

        # Interpolation helper
        function interp1(xs, ys, xq)
            ii = clamp(searchsortedlast(xs, xq), 1, length(xs)-1)
            α  = (xq - xs[ii]) / (xs[ii+1] - xs[ii])
            ys[ii] * (1-α) + ys[ii+1] * α
        end

        for d in _loaded
            push!(all_orbits, orbit_common)
            push!(all_dv_R,   [interp1(d.oc, d.dv_R, o) for o in orbit_common])
            push!(all_dv_T,   [interp1(d.oc, d.dv_T, o) for o in orbit_common])
            push!(all_dv_N_comp, [interp1(d.oc, d.dv_N, o) for o in orbit_common])
            push!(all_t,      [interp1(d.oc, d.t, o) for o in orbit_common])
            push!(N_vals,     length(d.dv_R))  # N is not used for label here
        end

        # Save cache: columns t, dv_R_N{N}, dv_T_N{N}, dv_N_N{N} per CSV
        n_csvs = length(_loaded)
        header = Matrix{String}(undef, 1, 1 + 3*n_csvs)
        header[1, 1] = "t"
        for k in 1:n_csvs
            header[1, 2 + 3*(k-1)] = "dv_R_N$(k)"
            header[1, 3 + 3*(k-1)] = "dv_T_N$(k)"
            header[1, 4 + 3*(k-1)] = "dv_N_N$(k)"
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
    all_orbit = all_orbits[1]  # all curves are on the same grid
    mean_dv  = mean(stack(all_dv), dims=2) |> vec
    std_dv   = std(stack(all_dv),  dims=2) |> vec

    # --- Plot ---
    plt = plot(title="Cumulative Δv_$comp_str — target satellite ($dir_tag)",
               xlabel="Orbits", ylabel="Δv_$comp_str (m/s)",
               formatter=:scientific, legend=:outertopright)

    for (dv_u, N) in zip(all_dv, N_vals)
        plot!(plt, all_orbit, dv_u, label="N=$N", color=:grey, alpha=0.5, lw=1)
    end

    if show_variance
        plot!(plt, all_orbit, mean_dv .+ std_dv,
              fillrange=mean_dv .- std_dv,
              fillalpha=0.2, fillcolor=:blue, linealpha=0, label="±1σ")
    end
    plot!(plt, all_orbit, mean_dv, label="Mean", color=:blue, lw=2)

    return (plt=plt,
            t=all_t[1],
            orbits=all_orbit,
            all_dv=all_dv,
            mean_dv=mean_dv,
            std_dv=std_dv)
end


