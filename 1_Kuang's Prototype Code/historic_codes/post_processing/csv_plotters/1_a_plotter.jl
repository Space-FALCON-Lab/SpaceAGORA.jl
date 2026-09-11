# Test reading CSV
import Makie
import WGLMakie
WGLMakie.activate!()
using LinearAlgebra
using Plots, Printf
using Statistics 
using DelimitedFiles  # For saving CSVs without extra deps

print("\033c")  # Clears the terminal on Windows

# --------- Constants ---------
const MU       = 3.986004418e14         # Earth μ [m^3/s^2]

include("../functions/5_OE_Converters.jl")
include("../functions/6_OE_and_dv_in_RTN.jl")
include("../functions/7_Plots.jl")
include("../functions/12_CSV_Write_Read.jl")

function _rbf_kernel_2d(X1::AbstractMatrix{<:Real}, X2::AbstractMatrix{<:Real}, ℓ1::Real, ℓ2::Real, σf2::Real)
    # Computes the RBF kernel matrix between two sets of 2D points, with separate length scales for each dimension.
    n1 = size(X1, 1)
    n2 = size(X2, 1)
    K = Matrix{Float64}(undef, n1, n2)
    @inbounds for i in 1:n1
        x1i = Float64(X1[i, 1])
        x2i = Float64(X1[i, 2])
        for j in 1:n2
            d1 = (x1i - Float64(X2[j, 1])) / Float64(ℓ1)
            d2 = (x2i - Float64(X2[j, 2])) / Float64(ℓ2)
            K[i, j] = Float64(σf2) * exp(-0.5 * (d1*d1 + d2*d2))
        end
    end
    return K
end

function fit_gpr_surface_and_plot(N_vals::AbstractVector, t_vals::AbstractVector, a_vals::AbstractVector;
                                  output_dir::AbstractString,
                                  fn_prefix::AbstractString="orbital_elements",
                                  fig_name::AbstractString="a_vs_N_t_gpr_surface.png",
                                  surface_csv_name::AbstractString="a_vs_N_t_gpr_surface.csv",
                                  max_train::Int=2200,
                                  grid_n_N::Int=45,
                                  grid_n_t::Int=60,
                                  ℓN::Real=1.0,
                                  ℓt::Real=1.0,
                                  noise_std::Real=0.03,
                                  scatter_cap::Int=5000,
                                  show_plot::Bool=true,
                                  save_plot::Bool=false)
    ## 1. Fits a GPR surface to the provided (N, t, a) data and generates plots. Also saves the fitted surface values to a CSV.
    n = length(N_vals)

    X = hcat(Float64.(N_vals), Float64.(t_vals))
    y = Float64.(a_vals)

    μx = vec(mean(X; dims=1))
    σx = vec(std(X; dims=1))
    σx .= map(s -> s > 0 ? s : 1.0, σx) # avoid division by zero if any feature has zero variance
    Xs = similar(X)
    Xs[:, 1] .= (X[:, 1] .- μx[1]) ./ σx[1]
    Xs[:, 2] .= (X[:, 2] .- μx[2]) ./ σx[2]

    μy = mean(y)
    σy = std(y)
    σy = σy > 0 ? σy : 1.0
    ys = (y .- μy) ./ σy

    n_train = min(max_train, n)
    train_idx = if n_train == n
        collect(1:n)
    else
        unique(round.(Int, range(1, n; length=n_train)))
    end
    Xtr = Xs[train_idx, :]
    ytr = ys[train_idx]

    K = _rbf_kernel_2d(Xtr, Xtr, ℓN, ℓt, 1.0)
    @inbounds for i in 1:size(K, 1)
        K[i, i] += noise_std^2 + 1e-8
    end

    F = cholesky(Symmetric(K))
    α = F \ ytr

    N_grid = collect(range(minimum(X[:, 1]), maximum(X[:, 1]); length=grid_n_N))
    t_grid = collect(range(minimum(X[:, 2]), maximum(X[:, 2]); length=grid_n_t))

    Xg = Matrix{Float64}(undef, length(N_grid) * length(t_grid), 2)
    p = 1
    for Ni in N_grid
        for tj in t_grid
            Xg[p, 1] = (Ni - μx[1]) / σx[1]
            Xg[p, 2] = (tj - μx[2]) / σx[2]
            p += 1
        end
    end

    Kgs = _rbf_kernel_2d(Xg, Xtr, ℓN, ℓt, 1.0)
    yg = Kgs * α
    ag = yg .* σy .+ μy
    # A_grid[i,j] = a at N_grid[i], t_grid[j]  (shape: n_N × n_t)
    # outer loop was over N, inner over t, so reshape then transpose
    A_grid = Matrix(reshape(ag, length(t_grid), length(N_grid))')

    ytr_hat = (_rbf_kernel_2d(Xtr, Xtr, ℓN, ℓt, 1.0) * α) .* σy .+ μy
    rmse_train = sqrt(mean((y[train_idx] .- ytr_hat).^2))
    println(@sprintf("GPR fit complete: train points=%d, RMSE(train)=%.6g m", length(train_idx), rmse_train))

    ## 2. Output
    # 2.1. Initialize output directory
    outdir = joinpath(output_dir, fn_prefix)
    mkpath(outdir)

    # 2.2. Save the fitted surface values to a CSV for potential future use
    open(joinpath(outdir, surface_csv_name), "w") do io
        writedlm(io, ["N" "t_s" "a_pred_m"], ',')
        surface_rows = Matrix{Float64}(undef, size(Xg, 1), 3)
        q = 1
        for Ni in N_grid
            for tj in t_grid
                surface_rows[q, 1] = Ni
                surface_rows[q, 2] = tj
                surface_rows[q, 3] = ag[q]
                q += 1
            end
        end
        writedlm(io, surface_rows, ',')
    end

    # 2.3. Initialize plotting
    keep_idx = if n <= scatter_cap 
        collect(1:n)
    else
        unique(round.(Int, range(1, n; length=scatter_cap))) # 
    end # prevent plotting to many points, which is very slow

    raw_color_min = minimum(a_vals[keep_idx])
    raw_color_max = maximum(a_vals[keep_idx])
    surface_min = minimum(A_grid)
    surface_max = maximum(A_grid)
    combined_color_range = (min(raw_color_min, surface_min), max(raw_color_max, surface_max))

    # 2.4. plot fig.1. 3D surface with scatter overlay 
    fig = Makie.Figure(size=(1100, 760))
    ax = Makie.Axis3(
        fig[1, 1],
        xlabel="N",
        ylabel="t (s)",
        zlabel="a (m)",
        title="GPR surface fit: a(N,t)"
    )
    Makie.surface!(ax, N_grid, t_grid, A_grid; alpha=0.75)
    Makie.scatter!(ax, Float64.(N_vals[keep_idx]), Float64.(t_vals[keep_idx]), Float64.(a_vals[keep_idx]); markersize=4.0)

    if show_plot
        Makie.display(fig)
    end
    if save_plot
        _save_makie_png(joinpath(outdir, fig_name), fig) # a_vs_N_t_gpr_surface.png
    end

    # 2.5. plot plot fig.2. raw data scatter only
    fig2 = Makie.Figure(size=(1100, 760))
    ax2 = Makie.Axis3(
        fig2[1, 1],
        xlabel="N",
        ylabel="t (s)",
        zlabel="a (m)",
        title="Target satellite semi-major axis: raw data"
    )
    Makie.scatter!(ax2, Float64.(N_vals[keep_idx]), Float64.(t_vals[keep_idx]), Float64.(a_vals[keep_idx]);
                   markersize=4.0)
    if show_plot
        Makie.display(fig2)
    end
    if save_plot
        _save_makie_png(joinpath(outdir, replace(fig_name, ".png" => "_scatter_only.png")), fig2) # a_vs_N_t_gpr_surface_scatter_only.png
    end

    # 2.6. plot fig.3. 2D scatter with color encoding a(N, t)
    fig3 = Makie.Figure(size=(1100, 600))
    ax3 = Makie.Axis(
        fig3[1, 1],
        xlabel="N",
        ylabel="t (s)",
        title="Target satellite a(N,t): raw scatter"
    )
    scatter2d = Makie.scatter!(
        ax3,
        Float64.(N_vals[keep_idx]),
        Float64.(t_vals[keep_idx]);
        color=Float64.(a_vals[keep_idx]),
        markersize=18,
        colormap=:plasma,
        colorrange=combined_color_range
    )
    Makie.Colorbar(fig3[1, 2], scatter2d, label="a (m)")
    if show_plot
        Makie.display(fig3)
    end
    if save_plot
        _save_makie_png(joinpath(outdir, replace(fig_name, ".png" => "_scatter2d.png")), fig3) # a_vs_N_t_gpr_surface_scatter2d.png
    end

    # 2.7. plot fig.4. 2D projection of fitted GPR surface
    fig4 = Makie.Figure(size=(1100, 600))
    ax4 = Makie.Axis(
        fig4[1, 1],
        xlabel="N",
        ylabel="t (s)",
        title="GPR surface projection"
    )
    heatmap_plot = Makie.heatmap!(
        ax4,
        N_grid,
        t_grid,
        A_grid;
        colormap=:plasma,
        colorrange=combined_color_range
    )
    Makie.Colorbar(fig4[1, 2], heatmap_plot, label="a (m)")
    if show_plot
        Makie.display(fig4)
    end
    if save_plot
        _save_makie_png(joinpath(outdir, replace(fig_name, ".png" => "_surface2d.png")), fig4) # a_vs_N_t_gpr_surface_surface2d.png
    end



    return fig
end

## 1. Initialzation
# 1.1. initialize plotting directory
# Plotting
output_dir = joinpath(@__DIR__, "..", "output", "images_reconstructed_from_csv") * "/"
println("Plotting results from CSV reconstruction...")
mkpath(output_dir)
for entry in readdir(output_dir; join=true)
    rm(entry; recursive=true, force=true)
end # clear old images

# 1.2. initialize CSV reading directory
data_dir = joinpath(@__DIR__, "..", "data")
Ns_requested = vcat(2, collect(6:5:161))
scatter_cache_csv = joinpath(data_dir, "a_vs_N_t_target_cache.csv")
rebuild_scatter_cache = false

# 1.3. initialize value arrays
N_vals = Float64[]
t_vals = Float64[]
a_vals = Float64[]

## 2. Load data from CSVs
# 2.1. first try to load from cache if it exists and we're not forcing a rebuild
if isfile(scatter_cache_csv) && !rebuild_scatter_cache # Check if cache file exists and we're not forcing a rebuild
    raw = readdlm(scatter_cache_csv, ',', Float64; skipstart=1)
    append!(N_vals, raw[:, 1])
    append!(t_vals, raw[:, 2])
    append!(a_vals, raw[:, 3])
    println("Loaded cached scatter data from $(scatter_cache_csv)")
end

# 2.2. If cache is empty (either didn't exist or was empty), build it from source CSV files
if isempty(N_vals)
    println("Building scatter cache from source CSV files...")
    for N_req in Ns_requested # Loop over the requested N values to find corresponding CSV files
        pat = Regex("^timeseries_N$(N_req)_.*\\.csv\$") # Pattern to match CSV files for the current N value (e.g., timeseries_N2_..., timeseries_N6_..., etc.)
        matches = filter(path -> occursin(pat, basename(path)), readdir(data_dir; join=true)) # Find all CSV files in the data directory that match the pattern for the current N value

        csv_path_N = first(matches) # Assuming there's only one CSV file per N value, take the first match.
        sol_N, _, _ = load_timeseries_csv(csv_path_N) # load x,y,z,vx,vy,vz time series from the CSV into a solution object structure

        elems = elements_time_series(sol_N, MU) #convert x, y, z, vx, vy, vz (r, v) to orbital elements time series (a,e,i,Ω,ω,ν) for each satellite
        oe_series = elems[end] # target = last satellite

        a = getfield.(oe_series, :a)
        t = sol_N.t
        npts = min(length(t), length(a))

        append!(N_vals, fill(Float64(N_req), npts))
        append!(t_vals, Float64.(t[1:npts]))
        append!(a_vals, Float64.(a[1:npts]))
    end

    open(scatter_cache_csv, "w") do io # Save the combined scatter data to a cache CSV for faster loading next time
        writedlm(io, ["N" "t_s" "a_m"], ',')
        writedlm(io, hcat(N_vals, t_vals, a_vals), ',')
    end
    println("Saved scatter cache to $(scatter_cache_csv)")
end

## 3. Plotting
# 3.1. Plot 3D scatter of a vs N and t with line
plot_a_vs_N_t_3d_from_data(
    N_vals,
    t_vals,
    a_vals;
    IMG_DIR=output_dir,
    fn_prefix="orbital_elements",
    fn="a_vs_N_t_3d.png",
    show_plot=false,
    save_plot=true
) 

# 3.2. Fit GPR surface, generate multiple plots, and save the fitted surface to a CSV
Base.invokelatest(
        fit_gpr_surface_and_plot,
        N_vals,
        t_vals,
        a_vals;
        output_dir=output_dir,
        fn_prefix="orbital_elements",
        fig_name="a_vs_N_t_gpr_surface.png",
        surface_csv_name="a_vs_N_t_gpr_surface.csv",
        show_plot=false,
        save_plot=true
    ) # use Base.invokelatest to avoid world age issues when calling from a script
println("Displayed 3D scatter window.")
