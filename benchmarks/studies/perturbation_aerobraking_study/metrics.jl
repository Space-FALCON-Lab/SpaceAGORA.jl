using LinearAlgebra
using StaticArrays
using Statistics

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

function _radial_velocity(pos::SVector{3,Float64}, vel::SVector{3,Float64})::Float64
    r = norm(pos)
    r <= 0.0 && return NaN
    return dot(pos, vel) / r
end

# Returns (altitude_m, geocentric_lat_rad) for a position vector under a given planet.
function _altitude_and_lat(pos::SVector{3,Float64}, planet)
    r = norm(pos)
    lat = asin(clamp(pos[3] / r, -1.0, 1.0))
    R_surface = _ellipsoid_radius_at_geocentric_lat(planet.Rp_e, planet.Rp_p, lat)
    return r - R_surface, lat
end

# Sample the solution at a set of times and return (times, positions, velocities).
function _sample_trajectory(sol, n_samples::Int=2000)
    t0, tf = sol.t[1], sol.t[end]
    ts = range(t0, tf; length=n_samples)
    positions = SVector{3,Float64}[]
    velocities = SVector{3,Float64}[]
    times = Float64[]
    for t in ts
        u = sol(t)
        pos = SVector{3,Float64}(u.sc[1].pos)
        vel = SVector{3,Float64}(u.sc[1].vel)
        if isfinite(norm(pos)) && isfinite(norm(vel))
            push!(times, Float64(t))
            push!(positions, pos)
            push!(velocities, vel)
        end
    end
    return times, positions, velocities
end

# Locate periapsis crossings (radial velocity ≤0 → >0) using linear interpolation
# over the dense solver time grid. Returns altitude_m at each periapsis.
function _periapsis_altitudes_m(sol, planet)
    ts = sol.t
    length(ts) < 2 && return Float64[]
    altitudes = Float64[]
    for k in 1:(length(ts) - 1)
        t0, t1 = Float64(ts[k]), Float64(ts[k+1])
        t1 > t0 || continue
        u0 = sol(t0); u1 = sol(t1)
        pos0 = SVector{3,Float64}(u0.sc[1].pos)
        vel0 = SVector{3,Float64}(u0.sc[1].vel)
        pos1 = SVector{3,Float64}(u1.sc[1].pos)
        vel1 = SVector{3,Float64}(u1.sc[1].vel)
        rv0 = _radial_velocity(pos0, vel0)
        rv1 = _radial_velocity(pos1, vel1)
        !(isfinite(rv0) && isfinite(rv1)) && continue
        rv0 <= 0.0 && rv1 > 0.0 || continue
        # Linear interpolation to periapsis
        frac = abs(rv0) / (abs(rv0) + abs(rv1))
        t_peri = t0 + frac * (t1 - t0)
        u_peri = sol(t_peri)
        pos_peri = SVector{3,Float64}(u_peri.sc[1].pos)
        alt, _ = _altitude_and_lat(pos_peri, planet)
        isfinite(alt) && push!(altitudes, alt)
    end
    return altitudes
end

# Peak dynamic pressure (Pa) and approximate Chapman heat rate (W/m²) along the
# atmospheric portion of the trajectory. Uses the planet's built-in exponential
# atmosphere for a lightweight density estimate; GRAM density is more accurate but
# requires a lock and is not called here to keep metrics extraction cheap.
function _peak_aero_metrics(sol, planet, EI_m::Float64, n_samples::Int=5000)
    t0, tf = sol.t[1], sol.t[end]
    ts = range(t0, tf; length=n_samples)
    q_dyn_peak = 0.0
    heat_rate_peak = 0.0
    ρ_ref = planet.ρ_ref
    H     = planet.H
    h_ref = planet.Rp_e   # sea-level reference ~ equatorial radius
    k_chapman = hasproperty(planet, :k) ? Float64(planet.k) : 1.9027e-4

    for t in ts
        u = sol(t)
        pos = SVector{3,Float64}(u.sc[1].pos)
        vel = SVector{3,Float64}(u.sc[1].vel)
        alt, _ = _altitude_and_lat(pos, planet)
        alt >= EI_m && continue   # outside atmosphere

        r = norm(pos)
        v = norm(vel)
        ρ = ρ_ref * exp(-(r - h_ref) / H)

        q_dyn = 0.5 * ρ * v^2
        heat  = k_chapman * sqrt(ρ) * v^3

        q_dyn_peak    = max(q_dyn_peak, q_dyn)
        heat_rate_peak = max(heat_rate_peak, heat)
    end
    return q_dyn_peak, heat_rate_peak
end

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

function extract_run_metrics(
    sol,
    args::SimulationConfiguration,
    cfg::BodyStudyConfig,
    ic_params::NamedTuple,
    pert_level::Symbol,
    body_name::Symbol,
    run_id::Int;
    solver_success::Bool=true,
    solve_retcode::String="Success",
    runtime_s::Float64=NaN
)::NamedTuple
    planet = args.environment_model.planet
    EI_m = args.environment_model.EI * 1e3

    spacecraft = args.dynamics_model.spacecraft[1]
    m_init = spacecraft.dry_mass + spacecraft.prop_mass

    # Final state
    u_final = sol(sol.t[end])
    m_final = Float64(u_final.sc[1].mass)

    # Periapsis evolution
    peri_alts = try
        _periapsis_altitudes_m(sol, planet)
    catch
        Float64[]
    end

    n_passes = length(peri_alts)
    h_peri_min = n_passes > 0 ? minimum(peri_alts) : NaN
    h_peri_max = n_passes > 0 ? maximum(peri_alts) : NaN
    h_peri_mean = n_passes > 0 ? mean(peri_alts) : NaN
    h_peri_std  = n_passes > 1 ? std(peri_alts) : 0.0

    n_corridor_violations = count(h -> h < cfg.h_corridor_min_m || h > cfg.h_corridor_max_m, peri_alts)

    # Final apoapsis altitude: highest altitude sampled from the last quarter of trajectory
    tf = sol.t[end]
    t_late = Float64(tf) * 0.75
    h_apo_final = NaN
    try
        ts_late = range(t_late, Float64(tf); length=500)
        alts_late = [begin
            u = sol(t)
            pos = SVector{3,Float64}(u.sc[1].pos)
            alt, _ = _altitude_and_lat(pos, planet)
            alt
        end for t in ts_late]
        h_apo_final = maximum(filter(isfinite, alts_late))
    catch
    end

    # ΔV from Tsiolkovsky: ΔV = Isp * g0 * ln(m_init / m_final)
    g0 = 9.80665
    isp = isempty(args.control_model.control_effectors) ? 300.0 :
          Float64(args.control_model.control_effectors[1].Isp[1])
    dv_total = (isfinite(m_final) && m_final > 0.0 && m_init > m_final) ?
               isp * g0 * log(m_init / m_final) : 0.0

    # Peak aerodynamic loads
    q_dyn_peak, heat_rate_peak = try
        _peak_aero_metrics(sol, planet, EI_m)
    catch
        NaN, NaN
    end

    return (
        run_id=run_id,
        body=String(body_name),
        pert_level=String(pert_level),
        inc_deg=ic_params.inc_deg,
        RAAN_deg=ic_params.RAAN_deg,
        omega_deg=ic_params.omega_deg,
        h_peri_initial_m=ic_params.h_peri_m,
        h_apo_initial_m=ic_params.h_apo_m,
        n_passes=n_passes,
        h_peri_min_m=h_peri_min,
        h_peri_max_m=h_peri_max,
        h_peri_mean_m=h_peri_mean,
        h_peri_std_m=h_peri_std,
        n_corridor_violations=n_corridor_violations,
        h_apo_final_m=h_apo_final,
        dv_total_ms=dv_total,
        q_dyn_peak_Pa=q_dyn_peak,
        heat_rate_peak_Wm2=heat_rate_peak,
        solver_success=solver_success,
        solve_retcode=solve_retcode,
        runtime_s=runtime_s,
    )
end

function failed_run_metrics(
    ::BodyStudyConfig,
    ic_params::NamedTuple,
    pert_level::Symbol,
    body_name::Symbol,
    run_id::Int;
    solve_retcode::String="Failed",
    runtime_s::Float64=NaN
)::NamedTuple
    return (
        run_id=run_id,
        body=String(body_name),
        pert_level=String(pert_level),
        inc_deg=ic_params.inc_deg,
        RAAN_deg=ic_params.RAAN_deg,
        omega_deg=ic_params.omega_deg,
        h_peri_initial_m=ic_params.h_peri_m,
        h_apo_initial_m=ic_params.h_apo_m,
        n_passes=0,
        h_peri_min_m=NaN,
        h_peri_max_m=NaN,
        h_peri_mean_m=NaN,
        h_peri_std_m=NaN,
        n_corridor_violations=0,
        h_apo_final_m=NaN,
        dv_total_ms=NaN,
        q_dyn_peak_Pa=NaN,
        heat_rate_peak_Wm2=NaN,
        solver_success=false,
        solve_retcode=solve_retcode,
        runtime_s=runtime_s,
    )
end
