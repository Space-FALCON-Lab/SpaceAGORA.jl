"""
This module provides plotting functions for visualizing satellite orbits and diagnostics.
"""

# Explicitly import plotting functions from Plots to avoid name conflicts
import Plots: plot, plot!, scatter!, savefig, twinx

"""
    Plot satellite orbits in the XY plane.
"""
function plot_orbits(sol; IMG_DIR = nothing, fn="satellite_orbits.png")

    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end

    N = Int(length(sol.u[1]) ÷ 6)
    plt = plot(title="Satellite Orbits (XY)", xlabel="x (m)", ylabel="y (m)", aspect_ratio=1)
    for i in 1:N
        xs = getindex.(sol.u, Ref(idx(i,1))); ys = getindex.(sol.u, Ref(idx(i,2)))
        plot!(plt, xs, ys, label="Sat $i")
    end
    scatter!(plt, [0.0],[0.0], label="Earth")
    outpath = joinpath(IMG_DIR, fn)
    savefig(plt, outpath); return plt
end

"""
    Plot orbital energy per satellite (no total overlay).

    Saves a single plot with one line per satellite.
"""
function plot_orbit_energy_individual_satellites(sol, p; IMG_DIR = nothing, fn="orbit_energy_individual.png")

    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end

    masses = p[:masses]
    mu = p[:mu]

    Evals = [orbital_energy(u, masses, mu) for u in sol.u]

    plt = plot(
        title="Orbital Energy Per Satellite",
        xlabel="t (s)",
        ylabel="Energy (J)",
        lw=1.5,
        legend=:outerright
    )
    N = length(Evals[1])
    for i in 1:N
        plot!(plt, sol.t, [Evals[k][i] for k in eachindex(sol.t)], label="Sat $i")
    end

    outpath = joinpath(IMG_DIR, fn)
    savefig(plt, outpath)
    return plt
end

"""
    Plot total orbital energy (sum over satellites) as a single series.
"""
function plot_orbit_energy_total(sol, p; IMG_DIR = nothing, fn="orbit_energy_total.png")

    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end

    masses = p[:masses]
    mu = p[:mu]

    Evals = [orbital_energy(u, masses, mu) for u in sol.u]
    Etot = [sum(Evals[k]) for k in eachindex(sol.t)]

    plt = plot(
        sol.t,
        Etot,
        title="Total Orbital Energy",
        xlabel="t (s)",
        ylabel="Energy (J)",
        label="Total",
        lw=3,
        ls=:dash,
        color=:red,
        formatter=:scientific
    )

    outpath = joinpath(IMG_DIR, fn)
    savefig(plt, outpath)
    return plt
end

"""
    Plot change in orbital energy over time (relative to t=0), one line per
    satellite plus the total system change.

    ΔE_i(t) = E_i(t) - E_i(0)
    ΔE_tot(t) = ΣᵢΔE_i(t)
"""
function plot_orbital_energy_change(sol, p; IMG_DIR = nothing, fn="orbital_energy_change.png", total_only::Bool=false)

    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end

    masses = p[:masses]
    mu = p[:mu]

    Evals = [orbital_energy(u, masses, mu) for u in sol.u]
    E0 = Evals[1]
    N = length(E0)

    plt = plot(
        title="Change in Orbital Energy",
        xlabel="t (s)",
        ylabel="ΔE (J)",
        lw=1.5,
        legend=:outerright,
        formatter=:scientific
    )
    if !total_only
        for i in 1:N
            plot!(plt, sol.t, [Evals[k][i] - E0[i] for k in eachindex(sol.t)], label="Sat $i")
        end
    end
    # total system ΔE
    plot!(plt, sol.t, [sum(Evals[k]) - sum(E0) for k in eachindex(sol.t)],
          label="Total", lw=3, ls=:dash, color=:red)

    outpath = joinpath(IMG_DIR, fn)
    savefig(plt, outpath)
    return plt
end

"""
    Plot angular momentum over time for each satellite.
"""
# plot individual and total angular momentum on the same y-axis
function plot_angmom(sol, masses; IMG_DIR = nothing, fn="angular_momentum.png")

    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end

    Hs = [angular_momentum(u, masses) for u in sol.u]
    plt = plot(title="Angular Momentum Over Time", xlabel="t (s)", ylabel="|h_i| (kg·m²/s)")
    N = length(masses)
    for i in 1:N
        plot!(plt, sol.t, [Hs[k][1][i] for k in eachindex(sol.t)], label="Sat $i")
    end
    Htot = [Hs[k][2] for k in eachindex(sol.t)]
    plot!(plt, sol.t, Htot, label="Total", lw=3, ls=:dash)
    outpath = joinpath(IMG_DIR, fn)
    savefig(plt, outpath); return plt
end

# plot individual and total angular momentum on separate y-axes, individual on the left and total on the right
function plot_angmom_two_axes(sol, masses; IMG_DIR = nothing, fn="angular_momentum.png")

    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end

    # Compute angular momentum for each satellite and total angular momentum
    Hs = [angular_momentum(u, masses) for u in sol.u]
    N = length(masses)

    # Plot individual satellite angular momentum on the left axis
    plt = plot(
        sol.t,
        [Hs[k][1][1] for k in eachindex(sol.t)],
        title="Angular Momentum Over Time",
        xlabel="t (s)",
        ylabel="|h_i| (kg·m²/s)",
        label="Sat 1",
        legend=:topleft,
        lw=1.5
    )
    for i in 2:N
        plot!(plt, sol.t, [Hs[k][1][i] for k in eachindex(sol.t)], label="Sat $i", lw=1.5)
    end

    # Compute total angular momentum
    Htot = [Hs[k][2] for k in eachindex(sol.t)]

    # Add a secondary axis for total angular momentum with scientific notation
    ax2 = twinx(plt)
    plot!(
        ax2,
        sol.t,
        Htot,
        label="Total",
        lw=3,
        ls=:dash,
        color=:red,
        ylabel="|H_total| (kg·m²/s)",
        ylims=(0, maximum(Htot) + 1e3),
        legend=:topright,
        formatter=:scientific
    )

    # Save the plot
    outpath = joinpath(IMG_DIR, fn)
    savefig(plt, outpath)
    return plt
end

"""
    Plot linear momentum over time for each satellite.
"""
# plot individual and total linear momentum on the same y-axis
function plot_momentum(sol, masses; IMG_DIR = nothing, fn="linear_momentum.png")
        
    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end

    Pvals = [total_momentum(u, masses) for u in sol.u]
    plt = plot(title="Net Satellite Momentum", xlabel="t (s)", ylabel="kg·m/s")
    N = length(masses)
    for i in 1:N
        plot!(plt, sol.t, [Pvals[k][3][i] for k in eachindex(sol.t)], label="Sat $i")
    end
    Pmag = [Pvals[k][2] for k in eachindex(sol.t)] # Magnitude of total momentum
    plot!(plt, sol.t, Pmag, label="Total", lw=3, ls=:dash)
    outpath = joinpath(IMG_DIR, fn)
    savefig(plt, outpath); return plt
end

# plot individual and total linear momentum on separate y-axes, individual on the left and total on the right
function plot_momentum_two_axes(sol, masses; IMG_DIR = nothing, fn="linear_momentum.png")
        
    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end

    # Compute momentum values for each satellite and total momentum
    Pvals = [total_momentum(u, masses) for u in sol.u]
    N = length(masses)

    # Plot individual satellite momentum on the left axis
    plt = plot(
        sol.t,
        [Pvals[k][3][1] for k in eachindex(sol.t)],
        title="Net Satellite Momentum",
        xlabel="t (s)",
        ylabel="Momentum (kg·m/s)",
        label="Sat 1",
        legend=:topleft,
        lw=1.5
    )
    for i in 2:N
        plot!(plt, sol.t, [Pvals[k][3][i] for k in eachindex(sol.t)], label="Sat $i", lw=1.5)
    end

    # Compute total momentum magnitude
    Pmag = [Pvals[k][2] for k in eachindex(sol.t)]

    # Add a secondary axis for total momentum with scientific notation
    ax2 = twinx(plt)
    plot!(
        ax2,
        sol.t,
        Pmag,
        label="Total",
        lw=3,
        ls=:dash,
        color=:red,
        ylabel="Total Momentum (kg·m/s)",
        formatter=:scientific,
        legend=:topright
    )

    # Save the plot
    outpath = joinpath(IMG_DIR, fn)
    savefig(plt, outpath)
    return plt
end

"""
    Plot range between satellites i and j over time, color-coded by in-range/out-of-range status.
    Optionally overlays atmospheric clearance.
"""
# Plot range vs time with color split by in-range/out-of-range; optional clearance overlay
function plot_link_range(sol, p; i::Int=1, j::Int=2, show_clearance::Bool=true,
                        units::Symbol=:km, IMG_DIR = nothing, fn::String="link_range_$(i)-$(j).png")

    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end

    t, rng_m, clr_m, inrng = link_status_time_series(sol, p; i=i, j=j)
    fac = units === :km ? 1e-3 : 1.0
    ulabel = units === :km ? "km" : "m"
    rng_u = rng_m .* fac
    clr_u = clr_m .* fac
    minR  = get(p, :min_range, 0.0) * fac
    maxR  = get(p, :max_range, Inf) * fac
    req_cl = get(p, :atm_clearance, 0.0) * fac

    # Split the curve into two series via NaNs
    rng_in  = copy(rng_u)
    rng_out = copy(rng_u)
    for k in eachindex(t)
        if  inrng[k]; rng_out[k] = NaN end
        if !inrng[k]; rng_in[k]  = NaN end
    end

    plt = plot(title="Link range for pair $(i)-$(j)",
            xlabel="t (s)", ylabel="range ($ulabel)")
    plot!(plt, t, rng_out, label="Out of range")
    plot!(plt, t, rng_in,  label="In range", lw=3)

    # Range gates, if finite
    if isfinite(minR) && minR > 0
        plot!(plt, t, fill(minR, length(t)), label="min_range", ls=:dash)
    end
    if isfinite(maxR)
        plot!(plt, t, fill(maxR, length(t)), label="max_range", ls=:dash)
    end

    if show_clearance
        # Secondary axis for atmosphere clearance
        ax2 = twinx(plt)
        plot!(ax2, t, clr_u, label="Atmos. clearance", ls=:dot,
            ylabel="clearance ($ulabel)")
        plot!(ax2, t, fill(req_cl, length(t)), label="required clearance", ls=:dash)
    end
    
    outpath = joinpath(IMG_DIR, fn)
    savefig(plt, outpath)
    return plt
end

"""
    Print final delta-v in RTN frame and plot time series for a specified satellite.
    
    Inputs:
        sol: ODE solution object
        p: parameter dictionary with keys:
            :N - number of satellites
            :masses - vector of masses for all satellites
            plus all keys required by laser_forces()
        sat: satellite index to report and plot (default 1)
        fn_prefix: filename prefix for saving plot (default "dv_RTN_sat")

    Returns:
        plt: plot object
"""
function report_and_plot_dv_RTN(sol, p; sat::Int=1, IMG_DIR=nothing, fn_prefix="dv_RTN_sat", 
                                R_only::Bool=false, T_only::Bool=false, N_only::Bool=false,
                                RTN_separate::Bool=false)
    
    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end

    t, Δv_hist = delta_v_RTN_time_series(sol, p)
    dv = Δv_hist[sat]

    if RTN_separate
        plt_R = plot(t, dv[1,:], label="Δv_R", xlabel="t (s)", ylabel="m/s",
                title="Δv_R component (RTN) Sat $sat")
        savefig(plt_R, joinpath(IMG_DIR, fn_prefix * "_R_only.png"))
        plt_T = plot(t, dv[2,:], label="Δv_T", xlabel="t (s)", ylabel="m/s",
                title="Δv_T component (RTN) Sat $sat")
        savefig(plt_T, joinpath(IMG_DIR, fn_prefix * "_T_only.png"))
        plt_N = plot(t, dv[3,:], label="Δv_N", xlabel="t (s)", ylabel="m/s",
                title="Δv_N component (RTN) Sat $sat")
        savefig(plt_N, joinpath(IMG_DIR, fn_prefix * "_N_only.png"))
        return plt_R, plt_T, plt_N
    elseif R_only
        plt = plot(t, dv[1,:], label="Δv_R", xlabel="t (s)", ylabel="m/s",
                title="Δv_R component (RTN) Sat $sat")
        outpath = joinpath(IMG_DIR, fn_prefix * "_R_only.png")
    elseif T_only
        plt = plot(t, dv[2,:], label="Δv_T", xlabel="t (s)", ylabel="m/s",
                title="Δv_T component (RTN) Sat $sat")
        outpath = joinpath(IMG_DIR, fn_prefix * "_T_only.png")
    elseif N_only
        plt = plot(t, dv[3,:], label="Δv_N", xlabel="t (s)", ylabel="m/s",
                title="Δv_N component (RTN) Sat $sat")
        outpath = joinpath(IMG_DIR, fn_prefix * "_N_only.png")
    else
        plt = plot(t, dv[1,:], label="Δv_R", xlabel="t (s)", ylabel="m/s",
                title="Δv components (RTN) Sat $sat")
        plot!(plt, t, dv[2,:], label="Δv_T"); plot!(plt, t, dv[3,:], label="Δv_N")
        outpath = joinpath(IMG_DIR, fn_prefix * ".png")
    end
    savefig(plt, outpath)
    return plt
end


"""
    Print final force in RTN frame and plot time series for a specified satellite.

    Inputs:
        sol: ODE solution object
        p: parameter dictionary (must support laser_forces())
        sat: satellite index to report and plot (default 1)
        fn_prefix: filename prefix for saving plot (default "F_RTN_sat")

    Returns:
        plt: plot object
"""
function report_and_plot_F_RTN(sol, p; sat::Int=1, IMG_DIR=nothing, fn_prefix="F_RTN_sat",
                               R_only::Bool=false, T_only::Bool=false, N_only::Bool=false, RTN_separate::Bool=false)

    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "images"))
    end

    t, F_RTN_hist, _ = laser_force_RTN_time_series(sol, p)
    F = F_RTN_hist[sat]

    if R_only
        plt = plot(t, F[1,:], label="F_R", xlabel="t (s)", ylabel="N",
                   title="F_R component (RTN) Sat $sat")
        outpath = joinpath(IMG_DIR, fn_prefix * "_R_only.png")
        savefig(plt, outpath)
    elseif T_only
        plt = plot(t, F[2,:], label="F_T", xlabel="t (s)", ylabel="N",
                   title="F_T component (RTN) Sat $sat")
        outpath = joinpath(IMG_DIR, fn_prefix * "_T_only.png")
        savefig(plt, outpath)
    elseif N_only
        plt = plot(t, F[3,:], label="F_N", xlabel="t (s)", ylabel="N",
                   title="F_N component (RTN) Sat $sat")
        outpath = joinpath(IMG_DIR, fn_prefix * "_N_only.png")
        savefig(plt, outpath)
    elseif RTN_separate
        plt = plot(t, F[1,:], label="F_R", xlabel="t (s)", ylabel="N",
                   title="F_R component (RTN) Sat $sat")
        outpath = joinpath(IMG_DIR, fn_prefix * "_R_only.png")
        savefig(plt, outpath)
        plt = plot(t, F[2,:], label="F_T", xlabel="t (s)", ylabel="N",
                    title="F_T component (RTN) Sat $sat")
        outpath = joinpath(IMG_DIR, fn_prefix * "_T_only.png")
        savefig(plt, outpath)
        plt = plot(t, F[3,:], label="F_N", xlabel="t (s)", ylabel="N",
                   title="F_N component (RTN) Sat $sat")
        outpath = joinpath(IMG_DIR, fn_prefix * "_N_only.png")
        savefig(plt, outpath)
    else
        plt = plot(t, F[1,:], label="F_R", xlabel="t (s)", ylabel="N",
                   title="F components (RTN) Sat $sat")
        plot!(plt, t, F[2,:], label="F_T")
        plot!(plt, t, F[3,:], label="F_N")
        outpath = joinpath(IMG_DIR, fn_prefix * ".png")
        savefig(plt, outpath)
    end


    return plt
end

"""
    Print final delta-P in RTN frame and plot time series for a specified satellite.

    Inputs:
        sol: ODE solution object
        p: parameter dictionary with keys:
            :N - number of satellites
            :masses - vector of masses for all satellites
            plus all keys required by laser_forces()
        sat: satellite index to report and plot (default 1)
        fn_prefix: filename prefix for saving plot (default "dP_RTN_sat")

    Returns:
        plt: plot object
"""
function report_and_plot_delta_P_RTN(sol, p; sat::Int=1, IMG_DIR=nothing, fn_prefix="dP_RTN_sat",
                               R_only::Bool=false, T_only::Bool=false, N_only::Bool=false)
    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end

    t, _, ΔP_RTN_hist = laser_force_RTN_time_series(sol, p)
    dP = ΔP_RTN_hist[sat]

    if R_only
        plt = plot(t, dP[1,:], label="ΔP_R", xlabel="t (s)", ylabel="kg·m/s",
                   title="ΔP_R component (RTN) Sat $sat")
        outpath = joinpath(IMG_DIR, fn_prefix * "_R_only.png")
    elseif T_only
        plt = plot(t, dP[2,:], label="ΔP_T", xlabel="t (s)", ylabel="kg·m/s",
                   title="ΔP_T component (RTN) Sat $sat")
        outpath = joinpath(IMG_DIR, fn_prefix * "_T_only.png")
    elseif N_only
        plt = plot(t, dP[3,:], label="ΔP_N", xlabel="t (s)", ylabel="kg·m/s",
                   title="ΔP_N component (RTN) Sat $sat")
        outpath = joinpath(IMG_DIR, fn_prefix * "_N_only.png")
    else
        plt = plot(t, dP[1,:], label="ΔP_R", xlabel="t (s)", ylabel="kg·m/s",
                   title="ΔP components (RTN) Sat $sat")
        plot!(plt, t, dP[2,:], label="ΔP_T")
        plot!(plt, t, dP[3,:], label="ΔP_N")
        outpath = joinpath(IMG_DIR, fn_prefix * ".png")
    end

    savefig(plt, outpath)
    return plt
end

"""
    Print and plot RTN position components over time for each satellite.

    Inputs:
        sol: ODE solution object
        p: parameter dictionary with keys:
            :N - number of satellites
            :masses - vector of masses for all satellites
            plus all keys required by laser_forces()
        fn_prefix: filename prefix for saving plots (default "r_RTN_sat")

    Returns:
        plt: plot object
"""
function report_and_plot_r_RTN(sol, p; sat::Int=1, IMG_DIR=nothing, fn_prefix="r_RTN_sat",
                                 R_only::Bool=false, T_only::Bool=false, N_only::Bool=false)
    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end

    # Compute RTN position time series
    t, r_RTN_hist, _, _ = r_v_a_RTN_time_series(sol, p)

    r = r_RTN_hist[sat]
    # @printf("\nInitial r (RTN) Sat %d:  R=% .6f  T=% .6f  N=% .6f  |r|=% .6f m/s\n",
    #         sat, r[1,1], r[2,1], r[3,1], norm(r[:,1]))
    # @printf("\nFinal r (RTN) Sat %d:  R=% .6f  T=% .6f  N=% .6f  |r|=% .6f m/s\n",
    #         sat, r[1,end], r[2,end], r[3,end], norm(r[:,end]))
    if R_only
        plt = plot(t, r[1,:], label="r_R", xlabel="t (s)", ylabel="m",
                title="r_R component (RTN) Sat $sat")
        outpath = joinpath(IMG_DIR, fn_prefix * "_R_only.png")
    elseif T_only
        plt = plot(t, r[2,:], label="r_T", xlabel="t (s)", ylabel="m",
                title="r_T component (RTN) Sat $sat")
        outpath = joinpath(IMG_DIR, fn_prefix * "_T_only.png")
    elseif N_only
        plt = plot(t, r[3,:], label="r_N", xlabel="t (s)", ylabel="m",
                title="r_N component (RTN) Sat $sat")
        outpath = joinpath(IMG_DIR, fn_prefix * "_N_only.png")
    else
        plt = plot(t, r[1,:], label="r_R", xlabel="t (s)", ylabel="m",
                title="r components (RTN) Sat $sat")
        plot!(plt, t, r[2,:], label="r_T"); plot!(plt, t, r[3,:], label="r_N")
        outpath = joinpath(IMG_DIR, fn_prefix * ".png")
    end
    savefig(plt, outpath)
    return plt
end

"""
    Print and plot RTN velocity components over time for each satellite.

    Inputs:
        sol: ODE solution object
        p: parameter dictionary with keys:
            :N - number of satellites
            :masses - vector of masses for all satellites
            plus all keys required by laser_forces()
        fn_prefix: filename prefix for saving plots (default "v_RTN_sat")

    Returns:
        plt: plot object
"""
function report_and_plot_v_RTN(sol, p; sat::Int=1, IMG_DIR=nothing, fn_prefix="v_RTN_sat",
                                 R_only::Bool=false, T_only::Bool=false, N_only::Bool=false)
    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end

    # Compute RTN velocity time series
    t, _, v_RTN_hist, _ = r_v_a_RTN_time_series(sol, p)

    v = v_RTN_hist[sat]
    # @printf("\nInitial v (RTN) Sat %d:  R=% .6f  T=% .6f  N=% .6f  |v|=% .6f m/s\n",
    #         sat, v[1,1], v[2,1], v[3,1], norm(v[:,1]))
    # @printf("\nFinal v (RTN) Sat %d:  R=% .6f  T=% .6f  N=% .6f  |v|=% .6f m/s\n",
    #         sat, v[1,end], v[2,end], v[3,end], norm(v[:,end]))
    if R_only
        plt = plot(t, v[1,:], label="v_R", xlabel="t (s)", ylabel="m/s",
                title="v_R component (RTN) Sat $sat")
        outpath = joinpath(IMG_DIR, fn_prefix * "_R_only.png")
    elseif T_only
        plt = plot(t, v[2,:], label="v_T", xlabel="t (s)", ylabel="m/s",
                    title="v_T component (RTN) Sat $sat")
        outpath = joinpath(IMG_DIR, fn_prefix * "_T_only.png")
    elseif N_only
        plt = plot(t, v[3,:], label="v_N", xlabel="t (s)", ylabel="m/s",
                    title="v_N component (RTN) Sat $sat")
        outpath = joinpath(IMG_DIR, fn_prefix * "_N_only.png")
    else
        plt = plot(t, v[1,:], label="v_R", xlabel="t (s)", ylabel="m/s",
                title="v components (RTN) Sat $sat")
        plot!(plt, t, v[2,:], label="v_T"); plot!(plt, t, v[3,:], label="v_N")
        outpath = joinpath(IMG_DIR, fn_prefix * ".png")
    end
    savefig(plt, outpath)
    return plt
end

"""
    Print and plot RTN acceleration components over time for each satellite.
    Optionally overlays gravitational and laser acceleration components.

    Inputs:
        sol: ODE solution object
        p: parameter dictionary with keys:
            :N - number of satellites
            :masses - vector of masses for all satellites
            plus all keys required by laser_forces()
        sat: satellite index to report and plot (default 1)
        fn_prefix: filename prefix for saving plots (default "a_RTN_sat")
        show_a: whether to show total acceleration (default true)
        show_a_gravity: whether to overlay gravitational acceleration (default false)
        show_a_laser: whether to overlay laser acceleration (default false)
        R_only, T_only, N_only: whether to plot only the specified RTN component (default false)

    Returns:
        plt: plot object
"""
function report_and_plot_a_RTN(sol, p; sat::Int=1, IMG_DIR=nothing, fn_prefix="a_RTN_sat",
                                    show_a::Bool=true, show_a_gravity::Bool=false, show_a_laser::Bool=false,
                                 R_only::Bool=false, T_only::Bool=false, N_only::Bool=false)
    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end

    # Compute RTN acceleration time series
    t, _, _, a_RTN_hist = r_v_a_RTN_time_series(sol, p)

    # Compute RTN acceleration time series
    _, a_gravity_RTN_hist, a_laser_RTN_hist = a_gravity_laser_RTN_time_series(sol, p)

    a = a_RTN_hist[sat]
    ag=a_gravity_RTN_hist[sat]
    al=a_laser_RTN_hist[sat]

    plt = plot()  # Initialize empty plot

    # @printf("\nInitial a (RTN) Sat %d:  R=% .6f  T=% .6f  N=% .6f  |a|=% .6f m/s²\n",
    #         sat, a[1,1], a[2,1], a[3,1], norm(a[:,1]))
    # @printf("\nFinal a (RTN) Sat %d:  R=% .6f  T=% .6f  N=% .6f  |a|=% .6f m/s²\n",
    #         sat, a[1,end], a[2,end], a[3,end], norm(a[:,end]))
    if R_only
        if show_a
            plot!(t, a[1,:], label="a_R", xlabel="t (s)", ylabel="m/s²",
                        title="a_R component (RTN) Sat $sat")
        end
        if show_a_gravity
            plot!(t, ag[1,:], label="a_gravity_R", lw=2)
        end
        if show_a_laser
            plot!(t, al[1,:], label="a_laser_R", lw=2)
        end
        outpath = joinpath(IMG_DIR, fn_prefix * "_R_only.png")
    elseif T_only
        if show_a
            plot!(t, a[2,:], label="a_T", xlabel="t (s)", ylabel="m/s²",
                        title="a_T component (RTN) Sat $sat")
        end
        if show_a_gravity
            plot!(t, ag[2,:], label="a_gravity_T", lw=2, ls=:dash)
        end
        if show_a_laser
            plot!(t, al[2,:], label="a_laser_T", lw=2, ls=:dot)
        end
        outpath = joinpath(IMG_DIR, fn_prefix * "_T_only.png")
    elseif N_only
        if show_a
            plot!(t, a[3,:], label="a_N", xlabel="t (s)", ylabel="m/s²",
                        title="a_N component (RTN) Sat $sat")
        end
        if show_a_gravity
            plot!(t, ag[3,:], label="a_gravity_N", lw=2, ls=:dash)
        end
        if show_a_laser
            plot!(t, al[3,:], label="a_laser_N", lw=2, ls=:dot)
        end
        outpath = joinpath(IMG_DIR, fn_prefix * "_N_only.png")
    elseif !R_only && !T_only && !N_only
        if show_a
            plot!(t, a[1,:], label="a_R", xlabel="t (s)", ylabel="m/s²",
                    title="a components (RTN) Sat $sat")
            plot!(plt, t, a[2,:], label="a_T")
            plot!(plt, t, a[3,:], label="a_N")
        end
        if show_a_gravity
            plot!(plt, t, ag[1,:], label="a_gravity_R", lw=2, ls=:dash)
            plot!(plt, t, ag[2,:], label="a_gravity_T", lw=2, ls=:dash)
            plot!(plt, t, ag[3,:], label="a_gravity_N", lw=2, ls=:dash)
        end
        if show_a_laser
            plot!(plt, t, al[1,:], label="a_laser_R", lw=2, ls=:dot)
            plot!(plt, t, al[2,:], label="a_laser_T", lw=2, ls=:dot)
            plot!(plt, t, al[3,:], label="a_laser_N", lw=2, ls=:dot)
        end
        outpath = joinpath(IMG_DIR, fn_prefix * ".png")
    end
    savefig(plt, outpath)
    return plt
end

"""
    Plot orbit energy over time for each satellite and total orbit energy.

    Inputs:
        sol: ODE solution object
        p: parameter dictionary with keys:
            :N - number of satellites
            :masses - vector of masses for all satellites
            plus all keys required by laser_forces()
        IMG_DIR: directory to save the plot (default is "../output/images")
        fn: filename for saving plot (default "orbit_energy.png")

    Returns:
        plt: plot object
"""
function plot_orbit_energy(sol, p; IMG_DIR = nothing, fn="orbit_energy.png", subplot::Bool=false)

    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end

    # Extract parameters from p
    masses = p[:masses]
    mu = p[:mu]

    # Compute orbital energy for each satellite and total orbital energy
    Evals = [orbital_energy(u, masses, mu) for u in sol.u]
    Etot = [sum(Evals[k]) for k in eachindex(sol.t)]

    if subplot
        plt = plot(layout=(2,1), size=(800,600))

        # Individual satellite orbital energy subplot
        N = length(Evals[1])
        plot!(plt[1], title="Orbital Energy Over Time", xlabel="t (s)", ylabel="Energy (J)",
            sol.t, [Evals[k][1] for k in eachindex(sol.t)], label="Sat 1", lw=1.5, legend=:outerright)
        for i in 2:N
            plot!(plt[1], sol.t, [Evals[k][i] for k in eachindex(sol.t)], label="Sat $i", lw=1.5, legend=:outerright)
        end

        # Total orbital energy subplot
        plot!(plt[2], sol.t, Etot,
            title="Total Orbital Energy Over Time", xlabel="t (s)", ylabel="Total Energy (J)",
            label="Total", lw=3, ls=:dash, color=:red, formatter=:scientific)
        
        # Save the plot
        outpath = joinpath(IMG_DIR, fn)
        savefig(plt, outpath)
        return plt
    else
    # Plot individual satellite orbital energy on the left axis
        plt = plot(
            title="Orbital Energy Over Time",
            xlabel="t (s)",
            ylabel="Energy (J)",
            label="Sat 1",
            sol.t,
            [Evals[k][1] for k in eachindex(sol.t)],
            lw=1.5
        )
        N = length(Evals[1])
        for i in 2:N
            plot!(plt, sol.t, [Evals[k][i] for k in eachindex(sol.t)], label="Sat $i", lw=1.5, legend=:topleft)
        end

        # Add a secondary axis for total orbital energy with scientific notation
        ax2 = twinx(plt)
        plot!(
            ax2,
            sol.t,
            Etot,
            label="Total",
            lw=3,
            ls=:dash,
            color=:red,
            ylabel="Total Energy (J)",
            formatter=:scientific,
            legend=:topright
        )

        # Save the plot
        outpath = joinpath(IMG_DIR, fn)
        savefig(plt, outpath)
        return plt
    end
end

function report_and_plot_OE(sol, μ; sat::Int=1, IMG_DIR=nothing, fn_prefix="orbital_elements")
    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end
    outdir = joinpath(IMG_DIR, fn_prefix)
    mkpath(outdir)

    elems = elements_time_series(sol, μ)
    oe_series = elems[sat]

    # extract time series
    a = getfield.(oe_series, :a)
    e = getfield.(oe_series, :e)
    i = getfield.(oe_series, :i)
    Ω = getfield.(oe_series, :Ω)
    first_oe = first(oe_series)
    ω = hasproperty(first_oe, :ω) ? getfield.(oe_series, :ω) :
        hasproperty(first_oe, :u) ? getfield.(oe_series, :u) : fill(NaN, length(oe_series))
    ν = hasproperty(first_oe, :ν) ? getfield.(oe_series, :ν) :
        hasproperty(first_oe, :v) ? getfield.(oe_series, :v) : fill(NaN, length(oe_series))
    u = hasproperty(first_oe, :u) ? getfield.(oe_series, :u) : fill(NaN, length(oe_series))
    t = sol.t

    # Prefer ω for non-circular; fall back to u for circular/equatorial
    mask_circ = .!(e .> 1e-8) .| isnan.(ν)
    ω_or_u = fill(NaN, length(t))
    @inbounds for k in eachindex(t)
        ω_or_u[k] = mask_circ[k] ? u[k] : ω[k]
    end

    plt_a = plot(t, a, label="a",
                 title="Semi-major axis Sat $sat", xlabel="t (s)", ylabel="m")
    savefig(plt_a, joinpath(outdir, "a_sat_$sat.png"))

    plt_e = plot(t, e, label="e",
                 title="Eccentricity Sat $sat", xlabel="t (s)", ylabel="")
    savefig(plt_e, joinpath(outdir, "e_sat_$sat.png"))

    plt_i = plot(t, i, label="i",
                 title="Inclination Sat $sat", xlabel="t (s)", ylabel="rad")
    savefig(plt_i, joinpath(outdir, "i_sat_$sat.png"))

    plt_Ω = plot(t, Ω, label="Ω",
                 title="RAAN Sat $sat", xlabel="t (s)", ylabel="rad")
    savefig(plt_Ω, joinpath(outdir, "Ω_sat_$sat.png"))

    # Plot ω (non-circular) or u (circular)
    plt_ω = plot(t, ω_or_u, label="ω or u",
                 title="ω (non-circ) or u (circ) Sat $sat", xlabel="t (s)", ylabel="rad")
    savefig(plt_ω, joinpath(outdir, "ω_or_u_sat_$sat.png"))

    # Plot ν when defined (will be NaN for circular cases)
    plt_v = plot(t, ν, label="ν",
                 title="True anomaly Sat $sat", xlabel="t (s)", ylabel="rad")
    savefig(plt_v, joinpath(outdir, "ν_sat_$sat.png"))

    return plt_a, plt_e, plt_i, plt_Ω, plt_ω, plt_v
end

function report_and_plot_OE_diff(sol, μ; sat1::Int=1, sat2::Int=2, IMG_DIR=nothing, fn_prefix="orbital_elements_diff")
    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end
    outdir = joinpath(IMG_DIR, fn_prefix)
    mkpath(outdir)

    elems = elements_time_series(sol, μ, rv2coe_2pi)
    oe1 = elems[sat1]
    oe2 = elems[sat2]
    t = sol.t

    # Helper to extract a field or NaN-fill if absent
    function _get(oe_series, sym)
        first_oe = first(oe_series)
        hasproperty(first_oe, sym) ? getfield.(oe_series, sym) : fill(NaN, length(oe_series))
    end

    # Extract per-satellite series
    a1, a2 = _get(oe1, :a), _get(oe2, :a)
    e1, e2 = _get(oe1, :e), _get(oe2, :e)
    i1, i2 = _get(oe1, :i), _get(oe2, :i)
    Ω1, Ω2 = _get(oe1, :Ω), _get(oe2, :Ω)

    # ω or u (circular fallback)
    function _get_ω_or_u(oe_series)
        first_oe = first(oe_series)
        ω  = hasproperty(first_oe, :ω) ? getfield.(oe_series, :ω) : fill(NaN, length(oe_series))
        ν  = hasproperty(first_oe, :ν) ? getfield.(oe_series, :ν) :
             hasproperty(first_oe, :v) ? getfield.(oe_series, :v) : fill(NaN, length(oe_series))
        u  = hasproperty(first_oe, :u) ? getfield.(oe_series, :u) : fill(NaN, length(oe_series))
        e  = _get(oe_series, :e)
        mask_circ = .!(e .> 1e-8) .| isnan.(ν)
        out = fill(NaN, length(oe_series))
        @inbounds for k in eachindex(out)
            out[k] = mask_circ[k] ? u[k] : ω[k]
        end
        return out
    end

    ω_or_u1 = _get_ω_or_u(oe1)
    ω_or_u2 = _get_ω_or_u(oe2)

    # Argument of latitude (ω+ν) — consistent across e thresholds.
    # For circular orbits rv2coe sets ν=u and ω=0, so ω+ν = u.
    # For slightly-eccentric orbits the ν reference shifts with the perigee,
    # so using ω+ν keeps the position angle continuous.
    function _arg_lat(oe_series)
        first_oe = first(oe_series)
        ω_v = hasproperty(first_oe, :ω) ? getfield.(oe_series, :ω) : zeros(length(oe_series))
        ν_v = hasproperty(first_oe, :ν) ? getfield.(oe_series, :ν) :
              hasproperty(first_oe, :v) ? getfield.(oe_series, :v) : fill(NaN, length(oe_series))
        u_v = hasproperty(first_oe, :u) ? getfield.(oe_series, :u) : fill(NaN, length(oe_series))
        e_v = _get(oe_series, :e)
        out = Vector{Float64}(undef, length(oe_series))
        @inbounds for k in eachindex(out)
            out[k] = e_v[k] > 1e-8 ? mod(ω_v[k] + ν_v[k], 2π) : mod(u_v[k], 2π)
        end
        return out
    end

    u1_eff = _arg_lat(oe1)
    u2_eff = _arg_lat(oe2)

    # Differences
    Δa      = a2 .- a1
    Δe      = e2 .- e1
    Δi      = i2 .- i1
    ΔΩ      = Ω2 .- Ω1
    Δω_or_u = ω_or_u2 .- ω_or_u1
    # Wrap angle difference to (-π, π] to avoid ±360° jumps at the atan branch cut
    Δν      = mod.(u2_eff .- u1_eff .+ π, 2π) .- π

    lbl = "Sat $sat2 − Sat $sat1"

    # Convert time axis to number of orbits using sat1's initial semi-major axis
    T_orbit  = 2π * sqrt(a1[1]^3 / μ)
    n_orbits = t ./ T_orbit

    plt_a = plot(n_orbits, Δa, legend=false,
                 xlabel="Number of orbits", ylabel="Δa, m")
    savefig(plt_a, joinpath(outdir, "da_sat$(sat1)_$(sat2).png"))
    savefig(plt_a, joinpath(outdir, "da_sat$(sat1)_$(sat2).pdf"))

    plt_e = plot(n_orbits, Δe, legend=false,
                 xlabel="Number of orbits", ylabel="Δe")
    savefig(plt_e, joinpath(outdir, "de_sat$(sat1)_$(sat2).png"))
    savefig(plt_e, joinpath(outdir, "de_sat$(sat1)_$(sat2).pdf"))

    plt_i = plot(n_orbits, rad2deg.(Δi), legend=false,
                 xlabel="Number of orbits", ylabel="Δi, deg")
    savefig(plt_i, joinpath(outdir, "di_sat$(sat1)_$(sat2).png"))
    savefig(plt_i, joinpath(outdir, "di_sat$(sat1)_$(sat2).pdf"))

    plt_Ω = plot(n_orbits, rad2deg.(ΔΩ), legend=false,
                 xlabel="Number of orbits", ylabel="ΔΩ, deg")
    savefig(plt_Ω, joinpath(outdir, "dΩ_sat$(sat1)_$(sat2).png"))
    savefig(plt_Ω, joinpath(outdir, "dΩ_sat$(sat1)_$(sat2).pdf"))

    plt_ω = plot(n_orbits, rad2deg.(Δω_or_u), legend=false,
                 xlabel="Number of orbits", ylabel="Δω/u, deg")
    savefig(plt_ω, joinpath(outdir, "dω_or_u_sat$(sat1)_$(sat2).png"))
    savefig(plt_ω, joinpath(outdir, "dω_or_u_sat$(sat1)_$(sat2).pdf"))

    plt_v = plot(n_orbits, rad2deg.(Δν), legend=false,
                 xlabel="Number of orbits", ylabel="Δν, deg")
    savefig(plt_v, joinpath(outdir, "dν_sat$(sat1)_$(sat2).png"))
    savefig(plt_v, joinpath(outdir, "dν_sat$(sat1)_$(sat2).pdf"))

    return plt_a, plt_e, plt_i, plt_Ω, plt_ω, plt_v
end

function report_and_plot_OE_diff_v1(sol, μ; sat1::Int=1, sat2::Int=2, IMG_DIR=nothing, fn_prefix="orbital_elements_diff")
    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end
    outdir = joinpath(IMG_DIR, fn_prefix)
    mkpath(outdir)

    elems = elements_time_series(sol, μ, rv2coe_2pi)
    oe1 = elems[sat1]
    oe2 = elems[sat2]
    t = sol.t

    # Helper to extract a field or NaN-fill if absent
    function _get(oe_series, sym)
        first_oe = first(oe_series)
        hasproperty(first_oe, sym) ? getfield.(oe_series, sym) : fill(NaN, length(oe_series))
    end

    # Extract per-satellite series
    a1, a2 = _get(oe1, :a), _get(oe2, :a)
    e1, e2 = _get(oe1, :e), _get(oe2, :e)
    i1, i2 = _get(oe1, :i), _get(oe2, :i)
    Ω1, Ω2 = _get(oe1, :Ω), _get(oe2, :Ω)

    # ω or u (circular fallback)
    function _get_ω_or_u(oe_series)
        first_oe = first(oe_series)
        ω  = hasproperty(first_oe, :ω) ? getfield.(oe_series, :ω) : fill(NaN, length(oe_series))
        ν  = hasproperty(first_oe, :ν) ? getfield.(oe_series, :ν) :
             hasproperty(first_oe, :v) ? getfield.(oe_series, :v) : fill(NaN, length(oe_series))
        u  = hasproperty(first_oe, :u) ? getfield.(oe_series, :u) : fill(NaN, length(oe_series))
        e  = _get(oe_series, :e)
        mask_circ = .!(e .> 1e-8) .| isnan.(ν)
        out = fill(NaN, length(oe_series))
        @inbounds for k in eachindex(out)
            out[k] = mask_circ[k] ? u[k] : ω[k]
        end
        return out
    end

    ω_or_u1 = _get_ω_or_u(oe1)
    ω_or_u2 = _get_ω_or_u(oe2)

    # Argument of latitude (ω+ν) — consistent across e thresholds.
    # For circular orbits rv2coe sets ν=u and ω=0, so ω+ν = u.
    # For slightly-eccentric orbits the ν reference shifts with the perigee,
    # so using ω+ν keeps the position angle continuous.
    function _arg_lat(oe_series)
        first_oe = first(oe_series)
        ω_v = hasproperty(first_oe, :ω) ? getfield.(oe_series, :ω) : zeros(length(oe_series))
        ν_v = hasproperty(first_oe, :ν) ? getfield.(oe_series, :ν) :
              hasproperty(first_oe, :v) ? getfield.(oe_series, :v) : fill(NaN, length(oe_series))
        u_v = hasproperty(first_oe, :u) ? getfield.(oe_series, :u) : fill(NaN, length(oe_series))
        e_v = _get(oe_series, :e)
        out = Vector{Float64}(undef, length(oe_series))
        @inbounds for k in eachindex(out)
            out[k] = e_v[k] > 1e-8 ? mod(ω_v[k] + ν_v[k], 2π) : mod(u_v[k], 2π)
        end
        return out
    end

    u1_eff = _arg_lat(oe1)
    u2_eff = _arg_lat(oe2)

    # Differences
    Δa      = a2 .- a1
    Δe      = e2 .- e1
    Δi      = i2 .- i1
    ΔΩ      = Ω2 .- Ω1
    Δω_or_u = ω_or_u2 .- ω_or_u1
    # Wrap angle difference to (-π, π] to avoid ±360° jumps at the atan branch cut
    Δν      = mod.(u2_eff .- u1_eff .+ π, 2π) .- π

    lbl = "Sat $sat2 − Sat $sat1"

    plt_a = plot(t, Δa, label=lbl,
                 title="Δa ($lbl)", xlabel="t, s", ylabel="Δa, m")
    savefig(plt_a, joinpath(outdir, "da_sat$(sat1)_$(sat2).png"))
    savefig(plt_a, joinpath(outdir, "da_sat$(sat1)_$(sat2).pdf"))

    plt_e = plot(t, Δe, label=lbl,
                 title="Δe ($lbl)", xlabel="t, s", ylabel="Δe")
    savefig(plt_e, joinpath(outdir, "de_sat$(sat1)_$(sat2).png"))
    savefig(plt_e, joinpath(outdir, "de_sat$(sat1)_$(sat2).pdf"))

    plt_i = plot(t, rad2deg.(Δi), label=lbl,
                 title="Δi ($lbl)", xlabel="t, s", ylabel="Δi, deg")
    savefig(plt_i, joinpath(outdir, "di_sat$(sat1)_$(sat2).png"))
    savefig(plt_i, joinpath(outdir, "di_sat$(sat1)_$(sat2).pdf"))

    plt_Ω = plot(t, rad2deg.(ΔΩ), label=lbl,
                 title="ΔΩ ($lbl)", xlabel="t, s", ylabel="ΔΩ, deg")
    savefig(plt_Ω, joinpath(outdir, "dΩ_sat$(sat1)_$(sat2).png"))
    savefig(plt_Ω, joinpath(outdir, "dΩ_sat$(sat1)_$(sat2).pdf"))

    plt_ω = plot(t, rad2deg.(Δω_or_u), label=lbl,
                 title="Δω or Δu ($lbl)", xlabel="t, s", ylabel="Δω/u, deg")
    savefig(plt_ω, joinpath(outdir, "dω_or_u_sat$(sat1)_$(sat2).png"))
    savefig(plt_ω, joinpath(outdir, "dω_or_u_sat$(sat1)_$(sat2).pdf"))

    plt_v = plot(t, rad2deg.(Δν), label=lbl,
                 title="Δν ($lbl)", xlabel="t, s", ylabel="Δν, deg")
    savefig(plt_v, joinpath(outdir, "dν_sat$(sat1)_$(sat2).png"))
    savefig(plt_v, joinpath(outdir, "dν_sat$(sat1)_$(sat2).pdf"))

    return plt_a, plt_e, plt_i, plt_Ω, plt_ω, plt_v
end

"""
    plot_laser_dE_time_series(sol, p; IMG_DIR, fn, total_only)

Plot cumulative laser work ΔE(t) [J] over time for each satellite pair (i→j)
and the system total, using laser_exchanges_time_series.
If total_only=true, only the system total line is plotted.
"""
function plot_laser_dE_time_series(sol, p; IMG_DIR=nothing, fn="laser_dE_timeseries.png", total_only::Bool=false)
    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end

    ΔP_series, ΔE_series, t = laser_exchanges_time_series(sol, p)

    plt = plot(title="Cumulative Laser Work ΔE(t)",
               xlabel="t (s)", ylabel="ΔE (J)",
               lw=1.5, legend=:outerright, formatter=:scientific)

    if !total_only
        sorted_keys = sort(collect(keys(ΔE_series)))
        for key in sorted_keys
            i, j = key
            plot!(plt, t, ΔE_series[key], label="($i→$j)")
        end
    end

    # System total
    ΔE_total = sum(values(ΔE_series))
    plot!(plt, t, ΔE_total, label="Total", lw=3, ls=:dash, color=:red)

    savefig(plt, joinpath(IMG_DIR, fn))
    return plt
end

"""
    plot_laser_dP_time_series(sol, p; IMG_DIR, fn, total_only)

Plot cumulative laser impulse |ΔP|(t) [kg·m/s] over time for each satellite pair (i→j)
and the magnitude of the vector sum, using laser_exchanges_time_series.
If total_only=true, only the magnitude of the total vector sum is plotted.
"""
function plot_laser_dP_time_series(sol, p; IMG_DIR=nothing, fn="laser_dP_timeseries.png", total_only::Bool=false)
    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end

    ΔP_series, ΔE_series, t = laser_exchanges_time_series(sol, p)
    sorted_keys = sort(collect(keys(ΔP_series)))

    plt = plot(title="Cumulative Laser Impulse |ΔP|(t)",
               xlabel="t (s)", ylabel="|ΔP| (kg·m/s)",
               lw=1.5, legend=:outerright, formatter=:scientific)

    if !total_only
        for key in sorted_keys
            i, j = key
            plot!(plt, t, norm.(ΔP_series[key]), label="($i→$j)")
        end
    end

    # Magnitude of vector sum across all pairs
    ΔP_total = [norm(sum(ΔP_series[key][k] for key in sorted_keys)) for k in eachindex(t)]
    plot!(plt, t, ΔP_total, label="|ΔP_total|", lw=3, ls=:dash, color=:red)

    savefig(plt, joinpath(IMG_DIR, fn))
    return plt
end


#######################################!!!!!!!!!!!!!!!!!!!!!!!!
function _save_makie_png(path::AbstractString, fig)
    try
        @eval import CairoMakie
        Makie.save(path, fig; backend=CairoMakie)
    catch
        Makie.save(path, fig)
    end
end

function _dv_component_series(dv_hist::AbstractMatrix{<:Real}, component::Symbol)
    if component === :mag
        return vec(sqrt.(sum(abs2, dv_hist; dims=1)))
    elseif component === :R
        return vec(dv_hist[1, :])
    elseif component === :T
        return vec(dv_hist[2, :])
    elseif component === :N
        return vec(dv_hist[3, :])
    else
        error("component must be one of :mag, :R, :T, :N")
    end
end

function _dv_component_label(component::Symbol)
    component === :mag && return "|Δv|"
    component === :R && return "Δv_R"
    component === :T && return "Δv_T"
    component === :N && return "Δv_N"
    error("component must be one of :mag, :R, :T, :N")
end

function _coerce_dv_component_values(dv_vals, component::Symbol)
    if dv_vals isa AbstractVector
        return Float64.(dv_vals)
    elseif dv_vals isa AbstractMatrix
        nrows, ncols = size(dv_vals)
        if nrows == 3
            return Float64.(_dv_component_series(dv_vals, component))
        elseif ncols == 3
            return Float64.(_dv_component_series(permutedims(dv_vals), component))
        else
            error("dv_vals matrix must have 3 rows or 3 columns to represent RTN components")
        end
    else
        error("dv_vals must be a vector or a matrix with RTN components")
    end
end

function _dv_data_length(dv_vals)
    if dv_vals isa AbstractVector
        return length(dv_vals)
    elseif dv_vals isa AbstractMatrix
        nrows, ncols = size(dv_vals)
        if nrows == 3
            return ncols
        elseif ncols == 3
            return nrows
        else
            error("dv_vals matrix must have 3 rows or 3 columns to represent RTN components")
        end
    else
        error("dv_vals must be a vector or a matrix with RTN components")
    end
end

function plot_a_vs_N_t_3d_from_data(N_vals, t_vals, a_vals; IMG_DIR=nothing,
                                    fn_prefix="orbital_elements", fn="a_vs_N_t_3d.png",
                                    markersize::Real=0.8,
                                    show_plot::Bool=false,
                                    save_plot::Bool=true)
    if !(length(N_vals) == length(t_vals) == length(a_vals))
        error("N_vals, t_vals, and a_vals must have the same length")
    end

    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end
    outdir = joinpath(IMG_DIR, fn_prefix)
    mkpath(outdir)

    try
        @eval import Makie
    catch
        error("Makie is not available in this session.")
    end

    fig = Makie.Figure(size=(1000, 700))
    ax = Makie.Axis3(
        fig[1, 1],
        xlabel="N",
        ylabel="t (s)",
        zlabel="a (m)",
        title="Semi-major axis vs N and time"
    )
    Makie.scatter!(ax, Float64.(N_vals), Float64.(t_vals), Float64.(a_vals); markersize=markersize)
    if show_plot
        Makie.display(fig)
    end
    if save_plot
        _save_makie_png(joinpath(outdir, fn), fig)
    end
    return fig
end

function plot_a_vs_N_t_3d(sols, Ns, μ; sat::Union{Int, Symbol}=:last, IMG_DIR=nothing,
                          fn_prefix="orbital_elements", fn="a_vs_N_t_3d.png",
                          markersize::Real=4.0,
                          show_plot::Bool=false,
                          save_plot::Bool=true)
    if length(sols) != length(Ns)
        error("length(sols) must match length(Ns)")
    end

    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end
    outdir = joinpath(IMG_DIR, fn_prefix)
    mkpath(outdir)

    N_vals = Float64[]
    t_vals = Float64[]
    a_vals = Float64[]

    for (sol, N) in zip(sols, Ns)
        elems = elements_time_series(sol, μ)
        sat_idx = sat === :last ? length(elems) : sat isa Int ? sat : error("sat must be an Int or :last")
        oe_series = elems[sat_idx]

        a = getfield.(oe_series, :a)
        t = sol.t
        npts = min(length(t), length(a))

        append!(N_vals, fill(Float64(N), npts))
        append!(t_vals, Float64.(t[1:npts]))
        append!(a_vals, Float64.(a[1:npts]))
    end
    return plot_a_vs_N_t_3d_from_data(N_vals, t_vals, a_vals;
                                      IMG_DIR=IMG_DIR, fn_prefix=fn_prefix, fn=fn,
                                      markersize=markersize, show_plot=show_plot, save_plot=save_plot)
end

function plot_a_vs_N_t_3d(sol_by_N::AbstractDict, μ; sat::Union{Int, Symbol}=:last, IMG_DIR=nothing,
                          fn_prefix="orbital_elements", fn="a_vs_N_t_3d.png",
                                                    markersize::Real=4.0,
                                                    show_plot::Bool=false,
                                                    save_plot::Bool=true)
    Ns = sort!(collect(keys(sol_by_N)))
    sols = [sol_by_N[N] for N in Ns]
    return plot_a_vs_N_t_3d(sols, Ns, μ; sat=sat, IMG_DIR=IMG_DIR, fn_prefix=fn_prefix, fn=fn,
                                                        markersize=markersize, show_plot=show_plot, save_plot=save_plot)
end

function plot_dv_vs_N_t_3d_from_data(N_vals, t_vals, dv_vals; IMG_DIR=nothing,
                                     fn_prefix="delta_v", fn="dv_vs_N_t_3d.png",
                                     component::Symbol=:mag,
                                     markersize::Real=0.8,
                                     show_plot::Bool=false,
                                     save_plot::Bool=true)
    dv_len = _dv_data_length(dv_vals)
    if !(length(N_vals) == length(t_vals) == dv_len)
        error("N_vals, t_vals, and dv_vals must have the same number of samples")
    end

    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end
    outdir = joinpath(IMG_DIR, fn_prefix)
    mkpath(outdir)

    try
        @eval import Makie
    catch
        error("Makie is not available in this session.")
    end

    dv_component_vals = _coerce_dv_component_values(dv_vals, component)
    label = _dv_component_label(component)

    fig = Makie.Figure(size=(1000, 700))
    ax = Makie.Axis3(
        fig[1, 1],
        xlabel="N",
        ylabel="t (s)",
        zlabel="$(label) (m/s)",
        title="$(label) vs N and time"
    )
    Makie.scatter!(ax, Float64.(N_vals), Float64.(t_vals), dv_component_vals; markersize=markersize)
    if show_plot
        Makie.display(fig)
    end
    if save_plot
        _save_makie_png(joinpath(outdir, fn), fig)
    end
    return fig
end

function plot_dv_vs_N_t_3d(sols, Ns, p; sat::Union{Int, Symbol}=:last, component::Symbol=:mag,
                           IMG_DIR=nothing, fn_prefix="delta_v", fn="dv_vs_N_t_3d.png",
                           markersize::Real=4.0,
                           show_plot::Bool=false,
                           save_plot::Bool=true)
    if length(sols) != length(Ns)
        error("length(sols) must match length(Ns)")
    end

    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end
    outdir = joinpath(IMG_DIR, fn_prefix)
    mkpath(outdir)

    N_vals = Float64[]
    t_vals = Float64[]
    dv_vals = Float64[]

    for (sol, N) in zip(sols, Ns)
        t, dv_hist_all = delta_v_RTN_time_series(sol, p)
        sat_count = length(dv_hist_all)
        sat_idx = sat === :last ? sat_count : sat isa Int ? sat : error("sat must be an Int or :last")
        if sat_idx < 1 || sat_idx > sat_count
            error("sat index $(sat_idx) out of bounds (1:$(sat_count))")
        end

        dv_series = _dv_component_series(dv_hist_all[sat_idx], component)
        npts = min(length(t), length(dv_series))

        append!(N_vals, fill(Float64(N), npts))
        append!(t_vals, Float64.(t[1:npts]))
        append!(dv_vals, Float64.(dv_series[1:npts]))
    end
    return plot_dv_vs_N_t_3d_from_data(N_vals, t_vals, dv_vals;
                                       IMG_DIR=IMG_DIR, fn_prefix=fn_prefix, fn=fn,
                                       component=component,
                                       markersize=markersize, show_plot=show_plot, save_plot=save_plot)
end

function plot_dv_vs_N_t_3d(sol_by_N::AbstractDict, p; sat::Union{Int, Symbol}=:last,
                           component::Symbol=:mag, IMG_DIR=nothing,
                           fn_prefix="delta_v", fn="dv_vs_N_t_3d.png",
                           markersize::Real=4.0,
                           show_plot::Bool=false,
                           save_plot::Bool=true)
    Ns = sort!(collect(keys(sol_by_N)))
    sols = [sol_by_N[N] for N in Ns]
    return plot_dv_vs_N_t_3d(sols, Ns, p; sat=sat, component=component, IMG_DIR=IMG_DIR,
                             fn_prefix=fn_prefix, fn=fn, markersize=markersize,
                             show_plot=show_plot, save_plot=save_plot)
end
#######################################!!!!!!!!!!!!!!!!!!!!!!!!

function report_and_plot_rp_ra(sol, μ, Rbody; sat::Int=1, IMG_DIR=nothing, fn_prefix="apogee_perigee")
    if IMG_DIR === nothing
        IMG_DIR = normpath(joinpath(@__DIR__, "..", "output", "images"))
    end
    outdir = joinpath(IMG_DIR, fn_prefix)
    mkpath(outdir)

    # --- Extract orbital elements ---
    elems = elements_time_series(sol, μ)
    oe_series = elems[sat]

    a = getfield.(oe_series, :a)
    e = getfield.(oe_series, :e)
    t = sol.t

    # --- Compute radii (distance from center) ---
    rp = a .* (1 .- e)         # perigee radius
    ra = a .* (1 .+ e)         # apogee radius

    # --- Convert to altitudes above surface ---
    hp = rp .- Rbody           # perigee altitude
    ha = ra .- Rbody           # apogee altitude

    # ============================
    # Plots
    # ============================

    plt_rp = plot(t, rp, label="r_p",
                  title="Perigee Radius Sat $sat", xlabel="t (s)", ylabel="m")
    savefig(plt_rp, joinpath(outdir, "rp_sat_$sat.png"))

    plt_ra = plot(t, ra, label="r_a",
                  title="Apogee Radius Sat $sat", xlabel="t (s)", ylabel="m")
    savefig(plt_ra, joinpath(outdir, "ra_sat_$sat.png"))

    plt_hp = plot(t, hp, label="h_p",
                  title="Perigee Altitude Sat $sat", xlabel="t (s)", ylabel="m")
    savefig(plt_hp, joinpath(outdir, "hp_sat_$sat.png"))

    plt_ha = plot(t, ha, label="h_a",
                  title="Apogee Altitude Sat $sat", xlabel="t (s)", ylabel="m")
    savefig(plt_ha, joinpath(outdir, "ha_sat_$sat.png"))

    return plt_rp, plt_ra, plt_hp, plt_ha
end
