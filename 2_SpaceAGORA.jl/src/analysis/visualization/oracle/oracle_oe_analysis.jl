# Orbital-element and RTN time-series analysis functions for ORACLE post-processing.
# _rtn_basis is defined in oracle_diagnostics.jl (included before this file in oracle.jl).

using LinearAlgebra
using StaticArrays
using DataFrames
using Printf

# Orbital elements for every satellite at every saved time step.
# coe_fn defaults to rv2coe (angles in (-π,π]); pass rv2coe_2pi for [0,2π).
function elements_time_series(sol::_FlatSol, μ::Real, coe_fn=rv2coe)
    N = Int(length(sol.u[1]) ÷ 6)
    elems = [Vector{NamedTuple}(undef, length(sol.t)) for _ in 1:N]
    for k in eachindex(sol.t)
        uk = sol.u[k]
        for i in 1:N
            r = @SVector [uk[idx(i,1)], uk[idx(i,2)], uk[idx(i,3)]]
            v = @SVector [uk[idx(i,4)], uk[idx(i,5)], uk[idx(i,6)]]
            elems[i][k] = coe_fn(r, v, Float64(μ))
        end
    end
    return elems
end

# Pretty-print initial and final orbital elements for each satellite.
function print_initial_final_elements(sol::_FlatSol, μ::Real; degrees::Bool=true)
    elems = elements_time_series(sol, μ)
    todeg(x) = degrees ? rad2deg(x) : x
    for i in eachindex(elems)
        e0 = elems[i][1]; ef = elems[i][end]
        println("\nSatellite $i:")
        @printf("  a0 = %12.3f m   af = %12.3f m   Δa = %+10.3f m\n", e0.a, ef.a, ef.a-e0.a)
        @printf("  e0 = %12.6g     ef = %12.6g     Δe = %+12.6g\n", e0.e, ef.e, ef.e-e0.e)
        @printf("  i0 = %9.6f°     if = %9.6f°     Δi = %+9.6f°\n", todeg(e0.i), todeg(ef.i), todeg(ef.i-e0.i))
        @printf("  Ω0 = %9.6f°     Ωf = %9.6f°     ΔΩ = %+9.6f°\n", todeg(e0.Ω), todeg(ef.Ω), todeg(ef.Ω-e0.Ω))
        if !isnan(e0.ν)
            @printf("  ω0 = %9.6f°     ωf = %9.6f°\n", todeg(e0.ω), todeg(ef.ω))
            @printf("  ν0 = %9.6f°     νf = %9.6f°\n", todeg(e0.ν), todeg(ef.ν))
        else
            @printf("  u0 = %9.6f°     uf = %9.6f°     Δu = %+9.6f°\n",
                    todeg(e0.u), todeg(ef.u), todeg(ef.u-e0.u))
        end
    end
end

# Position, velocity, and finite-difference acceleration in the RTN frame.
function r_v_a_RTN_time_series(sol::_FlatSol, p::Dict)
    N  = p[:N]
    t  = sol.t; u = sol.u; nT = length(t)

    r_RTN_hist = [zeros(3, nT) for _ in 1:N]
    v_RTN_hist = [zeros(3, nT) for _ in 1:N]
    a_RTN_hist = [zeros(3, nT) for _ in 1:N]

    a_ECI = [zeros(3, nT) for _ in 1:N]
    for i in 1:N
        for k in 1:nT
            if k == 1
                dt = t[2] - t[1]
                for d in 1:3
                    a_ECI[i][d,k] = (u[2][idx(i,3+d)] - u[1][idx(i,3+d)]) / dt
                end
            elseif k == nT
                dt = t[nT] - t[nT-1]
                for d in 1:3
                    a_ECI[i][d,k] = (u[nT][idx(i,3+d)] - u[nT-1][idx(i,3+d)]) / dt
                end
            else
                dt = t[k+1] - t[k-1]
                for d in 1:3
                    a_ECI[i][d,k] = (u[k+1][idx(i,3+d)] - u[k-1][idx(i,3+d)]) / dt
                end
            end
        end
    end

    for k in 1:nT
        uk = u[k]
        for j in 1:N
            r  = @SVector [uk[idx(j,1)], uk[idx(j,2)], uk[idx(j,3)]]
            v  = @SVector [uk[idx(j,4)], uk[idx(j,5)], uk[idx(j,6)]]
            a  = @SVector [a_ECI[j][1,k], a_ECI[j][2,k], a_ECI[j][3,k]]
            er, et, en = _rtn_basis(r, v)
            r_RTN_hist[j][:,k] = [dot(r,er), dot(r,et), dot(r,en)]
            v_RTN_hist[j][:,k] = [dot(v,er), dot(v,et), dot(v,en)]
            a_RTN_hist[j][:,k] = [dot(a,er), dot(a,et), dot(a,en)]
        end
    end
    return t, r_RTN_hist, v_RTN_hist, a_RTN_hist
end

# Analytic gravity + tracker-based laser acceleration in the RTN frame.
function a_gravity_laser_RTN_time_series(sol::_FlatSol, p::Dict)
    N  = p[:N]; mu = p[:mu]
    t  = sol.t; u = sol.u; nT = length(t)

    a_gravity_RTN_hist = [zeros(3, nT) for _ in 1:N]
    a_laser_RTN_hist   = [zeros(3, nT) for _ in 1:N]

    tracker = get(p, :tracker, nothing)

    for k in 1:nT
        uk = u[k]
        for j in 1:N
            r = @SVector [uk[idx(j,1)], uk[idx(j,2)], uk[idx(j,3)]]
            v = @SVector [uk[idx(j,4)], uk[idx(j,5)], uk[idx(j,6)]]
            rmag = norm(r)
            ag_ECI = (-mu / (rmag^3 + 1e-30)) .* r
            er, et, en = _rtn_basis(r, v)
            a_gravity_RTN_hist[j][:,k] = [dot(ag_ECI,er), dot(ag_ECI,et), dot(ag_ECI,en)]
        end
    end

    if !isnothing(tracker) && length(tracker.t_hist) >= 2
        th = tracker.t_hist
        for k_tr in 2:length(th)
            dt_tr = th[k_tr] - th[k_tr-1]
            dt_tr <= 0.0 && continue
            a_R = (tracker.dv_R_hist[k_tr] - tracker.dv_R_hist[k_tr-1]) / dt_tr
            a_T = (tracker.dv_T_hist[k_tr] - tracker.dv_T_hist[k_tr-1]) / dt_tr
            a_N = (tracker.dv_N_hist[k_tr] - tracker.dv_N_hist[k_tr-1]) / dt_tr
            k_sol = clamp(searchsortedlast(t, th[k_tr]), 1, nT)
            a_laser_RTN_hist[1][:,k_sol] .+= [a_R, a_T, a_N]
        end
    end

    return t, a_gravity_RTN_hist, a_laser_RTN_hist
end

# Cumulative ΔV (RTN) time series from the impulse tracker.
function delta_v_RTN_time_series(sol::_FlatSol, p::Dict)
    N  = p[:N]
    tracker = get(p, :tracker, nothing)

    if isnothing(tracker) || isempty(tracker.t_hist)
        t = sol.t
        Δv_RTN_hist = [zeros(3, length(t)) for _ in 1:N]
        return t, Δv_RTN_hist
    end

    t      = tracker.t_hist
    nT     = length(t)
    Δv_RTN_hist = [zeros(3, nT) for _ in 1:N]

    Δv_RTN_hist[1][1,:] = tracker.dv_R_hist
    Δv_RTN_hist[1][2,:] = tracker.dv_T_hist
    Δv_RTN_hist[1][3,:] = tracker.dv_N_hist

    return t, Δv_RTN_hist
end

# Laser force (N) and cumulative ΔP (kg·m/s) in the RTN frame for the target.
function laser_force_RTN_time_series(sol::_FlatSol, p::Dict)
    N      = p[:N]
    masses = p[:masses]
    tracker = get(p, :tracker, nothing)

    t_flat = sol.t; nT_flat = length(t_flat)
    F_RTN_hist  = [zeros(3, nT_flat) for _ in 1:N]
    ΔP_RTN_hist = [zeros(3, nT_flat) for _ in 1:N]

    if isnothing(tracker) || isempty(tracker.t_hist)
        return t_flat, F_RTN_hist, ΔP_RTN_hist
    end

    th    = tracker.t_hist
    nT_tr = length(th)
    m1    = masses[1]

    for k in 2:nT_tr
        dt = th[k] - th[k-1]
        dt <= 0.0 && continue
        F_R = m1 * (tracker.dv_R_hist[k] - tracker.dv_R_hist[k-1]) / dt
        F_T = m1 * (tracker.dv_T_hist[k] - tracker.dv_T_hist[k-1]) / dt
        F_N = m1 * (tracker.dv_N_hist[k] - tracker.dv_N_hist[k-1]) / dt
        k_sol = clamp(searchsortedlast(t_flat, th[k]), 1, nT_flat)
        F_RTN_hist[1][:,k_sol] .+= [F_R, F_T, F_N]
    end
    for k in 1:nT_flat
        k_tr = clamp(searchsortedlast(th, t_flat[k]), 1, nT_tr)
        ΔP_RTN_hist[1][:,k] = [m1*tracker.dv_R_hist[k_tr],
                                m1*tracker.dv_T_hist[k_tr],
                                m1*tracker.dv_N_hist[k_tr]]
    end

    return t_flat, F_RTN_hist, ΔP_RTN_hist
end

# Range, atmospheric clearance, and in-range flag time series for satellite pair (i, j).
function link_status_time_series(sol::_FlatSol, p::Dict; i::Int=1, j::Int=2)
    t  = sol.t; u = sol.u
    R_atm  = get(p, :R_atm,         6_471_000.0)
    atm_cl = get(p, :atm_clearance, 0.0)
    minR   = get(p, :min_range,     0.0)
    maxR   = get(p, :max_range,     Inf)

    nT = length(t)
    rng = zeros(nT); clr = zeros(nT); inrng = falses(nT)

    for k in 1:nT
        uk = u[k]
        ri = @SVector [uk[idx(i,1)], uk[idx(i,2)], uk[idx(i,3)]]
        rj = @SVector [uk[idx(j,1)], uk[idx(j,2)], uk[idx(j,3)]]
        rel = rj - ri
        d   = norm(rel)
        rng[k] = d
        t_min  = -dot(ri, rel) / (dot(rel, rel) + 1e-30)
        t_min  = clamp(t_min, 0.0, 1.0)
        r_mid  = ri + t_min * rel
        clr[k] = norm(r_mid) - R_atm
        inrng[k] = (d >= minR) && (d <= maxR) && (clr[k] >= atm_cl)
    end
    return t, rng, clr, inrng
end

# Build a timeseries DataFrame of cumulative ΔV and OE changes for the target.
function _build_timeseries_dataframe(
    opts::OracleOptions,
    sol,
    oe0_a::Float64,
    oe0_e::Float64,
    oe0_i::Float64,
    oe0_raan::Float64,
    orbit_counts::Vector{Float64},
    mu::Float64,
    tracker,
)::DataFrame
    tf = Float64(sol.t[end])
    times = collect(range(0.0, tf; length=opts.timeseries_points))
    sol_times = Float64.(sol.t)
    rows = NamedTuple[]
    case_id = _case_id(opts)
    for t in times
        r, v = _target_rv_at(sol, Float64(t))
        dv_r_ts, dv_t_ts, dv_n_ts = tracked_dv_at(tracker, Float64(t))
        oe = rv2coe(r, v, mu)
        oc_idx = clamp(searchsortedlast(sol_times, Float64(t)), 1, length(sol_times))
        push!(rows, (
            case_id=case_id,
            time_s=Float64(t),
            orbit=orbit_counts[oc_idx],
            helpers=opts.helpers,
            helper_altitude_km=opts.helper_altitude_km,
            target_altitude_km=opts.target_altitude_km,
            target_inclination_deg=opts.target_inclination_deg,
            helper_inclination_deg=opts.helper_inclination_deg,
            schedule=opts.schedule,
            laser_range_km=opts.laser_range_km,
            laser_power_w=opts.laser_power_w,
            magnification=opts.magnification,
            beta=opts.beta,
            eta=opts.eta,
            mass_kg=opts.mass_kg,
            dv_r_mps=dv_r_ts,
            dv_t_mps=dv_t_ts,
            dv_n_mps=dv_n_ts,
            da_m=oe.a - oe0_a,
            de=oe.e - oe0_e,
            di_deg=rad2deg(oe.i - oe0_i),
            draan_deg=rad2deg(oe.Ω - oe0_raan),
        ))
    end
    return DataFrame(rows)
end
