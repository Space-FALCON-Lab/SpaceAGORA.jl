# Build time-series DataFrame of cumulative ΔV and orbital element changes for the target
# satellite over the simulation.
function _build_timeseries_dataframe(
    opts::OracleCase2Options,
    sol,
    oe0,
    orbit_counts::Vector{Float64},   # cumulative orbit count at each sol.t point (from _orbit_count_from_sol)
    mu::Float64,
    tracker::LaserImpulseTracker,
)::DataFrame
    tf = Float64(sol.t[end])
    times = collect(range(0.0, tf; length=opts.timeseries_points))
    sol_times = Float64.(sol.t)      # saved time points corresponding to orbit_counts
    rows = NamedTuple[]
    case_id = _case_id(opts)
    for t in times
        r, v = _target_rv_at(sol, Float64(t))
        dv_r_ts, dv_t_ts, dv_n_ts = tracked_dv_at(tracker, Float64(t))
        oe = _rv_to_elements(r, v, mu)
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
            da_m=oe.a - oe0.a,
            de=oe.e - oe0.e,
            di_deg=rad2deg(oe.i - oe0.i),
            draan_deg=rad2deg(oe.raan - oe0.raan),
        ))
    end
    return DataFrame(rows)
end

# push!() is a Julia standard library function that appends one or more items to the end of a collection (array, vector, etc.) in-place (the ! convention means it mutates its first argument).
# v = [1, 2, 3]
# push!(v, 4)   
# v is now [1, 2, 3, 4]
