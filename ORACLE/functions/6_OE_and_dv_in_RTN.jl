# =============================================================================
# RTN basis (public alias for _rtn_basis; used by ported diagnostic functions)
# =============================================================================
@inline function rtn_basis(r::SVector{3,Float64}, v::SVector{3,Float64})
    return _rtn_basis(r, v)  # delegates to 5_OE_Converters._rtn_basis
end

# =============================================================================
# Orbital elements time series (works on _FlatSol or any sol with flat .u)
# =============================================================================

# Elements for every satellite and every time sample.
# coe_fn defaults to rv2coe (angles in (-π, π]); pass rv2coe_2pi for [0, 2π).
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

# =============================================================================
# RTN position / velocity / acceleration time series (works on _FlatSol)
# =============================================================================
# Acceleration is computed via finite differences of velocity (no need for
# ODE right-hand side re-evaluation — avoids coupling to the solver internals).
# =============================================================================
function r_v_a_RTN_time_series(sol::_FlatSol, p::Dict)
    N  = p[:N]
    t  = sol.t; u = sol.u; nT = length(t)

    r_RTN_hist = [zeros(3, nT) for _ in 1:N]
    v_RTN_hist = [zeros(3, nT) for _ in 1:N]
    a_RTN_hist = [zeros(3, nT) for _ in 1:N]

    # Pre-compute finite-difference acceleration for all satellites
    # a[k] ≈ (v[k+1] - v[k-1]) / (t[k+1] - t[k-1])  (central diff, interior)
    #        (v[2] - v[1]) / (t[2] - t[1])              (forward, k=1)
    #        (v[end] - v[end-1]) / (t[end] - t[end-1])  (backward, k=end)
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
            er, et, en = rtn_basis(r, v)
            r_RTN_hist[j][:,k] = [dot(r,er), dot(r,et), dot(r,en)]
            v_RTN_hist[j][:,k] = [dot(v,er), dot(v,et), dot(v,en)]
            a_RTN_hist[j][:,k] = [dot(a,er), dot(a,et), dot(a,en)]
        end
    end
    return t, r_RTN_hist, v_RTN_hist, a_RTN_hist
end

# RTN gravity and laser acceleration time series.
# Gravity is computed analytically; laser acceleration is derived from the
# impulse tracker (stored in p[:tracker]) for the target satellite (index 1).
# All helper satellites have zero laser acceleration.
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
            er, et, en = rtn_basis(r, v)
            a_gravity_RTN_hist[j][:,k] = [dot(ag_ECI,er), dot(ag_ECI,et), dot(ag_ECI,en)]
        end
    end

    # Fill laser acceleration for the target using tracker dV increments / dt
    if !isnothing(tracker) && length(tracker.t_hist) >= 2
        th = tracker.t_hist
        for k_tr in 2:length(th)
            dt_tr = th[k_tr] - th[k_tr-1]
            dt_tr <= 0.0 && continue
            a_R = (tracker.dv_R_hist[k_tr] - tracker.dv_R_hist[k_tr-1]) / dt_tr
            a_T = (tracker.dv_T_hist[k_tr] - tracker.dv_T_hist[k_tr-1]) / dt_tr
            a_N = (tracker.dv_N_hist[k_tr] - tracker.dv_N_hist[k_tr-1]) / dt_tr
            # Map tracker time to flat sol index
            k_sol = clamp(searchsortedlast(t, th[k_tr]), 1, nT)
            a_laser_RTN_hist[1][:,k_sol] .+= [a_R, a_T, a_N]
        end
    end

    return t, a_gravity_RTN_hist, a_laser_RTN_hist
end

# =============================================================================
# Delta-V RTN time series (uses impulse tracker for target; zero for helpers)
# =============================================================================
function delta_v_RTN_time_series(sol::_FlatSol, p::Dict)
    N  = p[:N]
    tracker = get(p, :tracker, nothing)

    # Fall back to zero histories if no tracker is available
    if isnothing(tracker) || isempty(tracker.t_hist)
        t = sol.t
        Δv_RTN_hist = [zeros(3, length(t)) for _ in 1:N]
        return t, Δv_RTN_hist
    end

    t      = tracker.t_hist
    nT     = length(t)
    Δv_RTN_hist = [zeros(3, nT) for _ in 1:N]

    # Target satellite (index 1) dV from tracker (already cumulative RTN)
    Δv_RTN_hist[1][1,:] = tracker.dv_R_hist
    Δv_RTN_hist[1][2,:] = tracker.dv_T_hist
    Δv_RTN_hist[1][3,:] = tracker.dv_N_hist
    # Helper satellites have zero laser dV (Δv_RTN_hist[i] already zeros for i>1)

    return t, Δv_RTN_hist
end

# =============================================================================
# Laser force RTN time series (differentiate tracker dV; target only)
# =============================================================================
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
    # Cumulative ΔP at each flat sol time by interpolation from tracker
    for k in 1:nT_flat
        k_tr = clamp(searchsortedlast(th, t_flat[k]), 1, nT_tr)
        ΔP_RTN_hist[1][:,k] = [m1*tracker.dv_R_hist[k_tr],
                                m1*tracker.dv_T_hist[k_tr],
                                m1*tracker.dv_N_hist[k_tr]]
    end

    return t_flat, F_RTN_hist, ΔP_RTN_hist
end

# =============================================================================
# Link range / clearance / status time series for a satellite pair (i, j)
# =============================================================================
function link_status_time_series(sol::_FlatSol, p::Dict; i::Int=1, j::Int=2)
    t      = sol.t; u = sol.u
    # Default atmosphere top = R_Earth(6371 km) + 100 km = 6 471 km
    _R_atm_default = 6_471_000.0
    R_atm  = get(p, :R_atm,         _R_atm_default)
    atm_cl = get(p, :atm_clearance, 0.0)
    minR   = get(p, :min_range,     0.0)
    maxR   = get(p, :max_range,     Inf)

    nT  = length(t)
    rng = zeros(nT)
    clr = zeros(nT)
    inrng = falses(nT)

    for k in 1:nT
        uk = u[k]
        ri = @SVector [uk[idx(i,1)], uk[idx(i,2)], uk[idx(i,3)]]
        rj = @SVector [uk[idx(j,1)], uk[idx(j,2)], uk[idx(j,3)]]
        rel = rj - ri
        d   = norm(rel)
        rng[k] = d
        # Minimum approach distance to Earth center along the chord ri→rj
        # (for the atmospheric clearance check)
        t_min = -dot(ri, rel) / (dot(rel, rel) + 1e-30)
        t_min = clamp(t_min, 0.0, 1.0)
        r_mid  = ri + t_min * rel
        min_r  = norm(r_mid)
        clr[k] = min_r - R_atm
        in_range = (d >= minR) && (d <= maxR) && (clr[k] >= atm_cl)
        inrng[k] = in_range
    end
    return t, rng, clr, inrng
end

# =============================================================================
# Target r/v accessor kept for backward-compat (SpaceAGORA structured sol)
# =============================================================================
# get rv at any time t from the solution object (sol) for the target satellite (sc[1])
@inline function _target_rv_at(sol, t::Float64)
    u = if isapprox(t, Float64(sol.t[1]); atol=0.0, rtol=0.0)
        sol.u[1]
    elseif isapprox(t, Float64(sol.t[end]); atol=0.0, rtol=0.0)
        sol.u[end]
    else
        sol(t)
    end
    sc = u.sc[1]
    return SVector{3, Float64}(sc.pos), SVector{3, Float64}(sc.vel)
end
