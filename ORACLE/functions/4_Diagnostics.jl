# =============================================================================
# Flat-vector state adapter for SpaceAGORA's structured sol.u format
# =============================================================================
#
# SpaceAGORA sol.u[k] is a structured state with u.sc[i].pos/vel fields.
# The prototype's 7_Plots.jl (and diagnostic functions) expect flat vectors of
# length 6N, indexed via idx(i, off) = 6*(i-1) + off.
#
# _FlatSol wraps a SpaceAGORA solution object, exposing flat-vector .u so that
# all the ported physics helpers work unchanged.
# =============================================================================

# Flat-vector index helper: satellite i (1-based), state offset off (1=x…3=z, 4=vx…6=vz)
@inline idx(i, off) = 6*(i-1) + off

# Convert one SpaceAGORA structured state to a plain Float64 flat vector of length 6N
function _sa_to_flat_u(u, N::Int)::Vector{Float64}
    flat = Vector{Float64}(undef, 6N)
    for i in 1:N
        sc = u.sc[i]
        flat[idx(i,1)] = sc.pos[1]; flat[idx(i,2)] = sc.pos[2]; flat[idx(i,3)] = sc.pos[3]
        flat[idx(i,4)] = sc.vel[1]; flat[idx(i,5)] = sc.vel[2]; flat[idx(i,6)] = sc.vel[3]
    end
    return flat
end

# Lightweight sol adapter that presents flat-vector .u to downstream functions
struct _FlatSol
    t::Vector{Float64}
    u::Vector{Vector{Float64}}
end

# Build a _FlatSol from a SpaceAGORA solution and total spacecraft count N
function _make_flat_sol(sa_sol, N::Int)::_FlatSol
    t = Float64.(sa_sol.t)
    u = [_sa_to_flat_u(sa_sol.u[k], N) for k in eachindex(sa_sol.t)]
    return _FlatSol(t, u)
end

# Build a _FlatSol directly from a feather/Arrow DataFrame written by SpaceAGORA's
# native output pipeline (System 2 — SavedValues callback).
# Column naming follows SpaceAGORA's convention:
#   sc{i}_pos_1/2/3  ← ECI position x/y/z [m] for satellite i
#   sc{i}_vel_1/2/3  ← ECI velocity x/y/z [m/s] for satellite i
# Use this instead of _make_flat_sol when the ODE solution object is not
# needed (e.g. when SPACEAGORA_SOLVER_SAVE_EVERYSTEP=false).
function _make_flat_sol_from_feather(df::DataFrame, N::Int)::_FlatSol
    t = Float64.(df.time)
    nrows = length(t)
    # Pre-extract column arrays to avoid repeated column-name lookups in the inner loop
    pos_cols = [(df[!, "sc$(i)_pos_1"], df[!, "sc$(i)_pos_2"], df[!, "sc$(i)_pos_3"]) for i in 1:N]
    vel_cols = [(df[!, "sc$(i)_vel_1"], df[!, "sc$(i)_vel_2"], df[!, "sc$(i)_vel_3"]) for i in 1:N]
    u = Vector{Vector{Float64}}(undef, nrows)
    for k in 1:nrows
        flat = Vector{Float64}(undef, 6N)
        for i in 1:N
            flat[idx(i,1)] = pos_cols[i][1][k]
            flat[idx(i,2)] = pos_cols[i][2][k]
            flat[idx(i,3)] = pos_cols[i][3][k]
            flat[idx(i,4)] = vel_cols[i][1][k]
            flat[idx(i,5)] = vel_cols[i][2][k]
            flat[idx(i,6)] = vel_cols[i][3][k]
        end
        u[k] = flat
    end
    return _FlatSol(t, u)
end

# Orbit counter that works on _FlatSol (System 2 / feather data).
# Mirrors _orbit_count_from_sol but reads the flat-vector u layout instead of
# SpaceAGORA's structured state.  With 1001 feather points this is as accurate
# as the full-solution version at a fraction of the memory cost.
function _orbit_count_from_flat_sol(flat_sol::_FlatSol, mu::Float64)
    T_series = [begin
        r = SVector{3,Float64}(u[idx(1,1)], u[idx(1,2)], u[idx(1,3)])
        v = SVector{3,Float64}(u[idx(1,4)], u[idx(1,5)], u[idx(1,6)])
        a = _rv_to_elements(r, v, mu).a
        2pi * sqrt(a^3 / mu)
    end for u in flat_sol.u]
    dt    = diff(flat_sol.t)
    T_mid = (T_series[1:end-1] .+ T_series[2:end]) ./ 2
    return cumsum([0.0; dt ./ T_mid])
end

# =============================================================================
# Orbital mechanics diagnostics — work on flat-vector u (compatible with _FlatSol)
# =============================================================================

# Total linear momentum: returns (P_vec, |P|, per-sat_magnitude_vector)
function total_momentum(u::AbstractVector{<:Real}, masses::AbstractVector{<:Real})
    N = length(masses)
    P = zeros(3); pmag = zeros(N)
    @inbounds for i in 1:N
        v = @SVector [u[idx(i,4)], u[idx(i,5)], u[idx(i,6)]]
        pi = masses[i] * v
        P .+= pi; pmag[i] = norm(pi)
    end
    return P, norm(P), pmag
end

# Angular momentum: returns (per-sat |h_i|, total |H|)
function angular_momentum(u::AbstractVector{<:Real}, masses::AbstractVector{<:Real})
    N = length(masses)
    h_each = zeros(N); H = zeros(3)
    @inbounds for i in 1:N
        r = @SVector [u[idx(i,1)], u[idx(i,2)], u[idx(i,3)]]
        p = @SVector [masses[i]*u[idx(i,4)], masses[i]*u[idx(i,5)], masses[i]*u[idx(i,6)]]
        h = cross(r, p)
        h_each[i] = norm(h); H .+= h
    end
    return h_each, norm(H)
end

# Orbital energy: returns per-satellite mechanical energy vector [J]
function orbital_energy(u::AbstractVector{<:Real}, masses::AbstractVector{<:Real}, mu::Real)
    N = length(masses); E = zeros(N)
    @inbounds for i in 1:N
        x,y,z   = u[idx(i,1)], u[idx(i,2)], u[idx(i,3)]
        r       = sqrt(x^2 + y^2 + z^2) + 1e-12
        vx,vy,vz = u[idx(i,4)], u[idx(i,5)], u[idx(i,6)]
        K = 0.5 * masses[i] * (vx^2 + vy^2 + vz^2)
        U = -mu * masses[i] / r
        E[i] = K + U
    end
    return E
end

# =============================================================================
# SpaceAGORA-specific laser exchange diagnostics (use impulse tracker)
# =============================================================================
# These replace the prototype's laser_forces-dependent evaluate_laser_exchanges
# and laser_exchanges_time_series for the SpaceAGORA case-2 runner, where the
# impulse tracker accumulates dV for the target satellite only.
#
# The aggregate laser-pair key is (0, 1): "laser ensemble → target (index 1)".
# =============================================================================

# Evaluate total ΔP and ΔE delivered to the target by the laser ensemble.
# ΔE is approximated as Σ ΔP·v_target_RTN at each tracker step.
# sa_sol: the SpaceAGORA structured sol (for target velocity lookup); or pass
# via p[:sa_sol] when called from 7_Plots.jl with a flat sol as first arg.
function evaluate_laser_exchanges(sol_or_sa, p::Dict)
    tracker = get(p, :tracker, nothing)
    if tracker === nothing || isempty(tracker.t_hist)
        return Dict{Tuple{Int,Int},SVector{3,Float64}}(), Dict{Tuple{Int,Int},Float64}()
    end
    # Resolve SpaceAGORA structured sol for velocity lookup
    sa_sol = sol_or_sa isa _FlatSol ? get(p, :sa_sol, nothing) : sol_or_sa
    mass    = p[:masses][1]               # target mass [kg]
    t_hist  = tracker.t_hist
    dv_R    = tracker.dv_R_hist
    dv_T    = tracker.dv_T_hist
    dv_N    = tracker.dv_N_hist
    nT = length(t_hist)

    ΔP_total = SVector(mass*dv_R[end], mass*dv_T[end], mass*dv_N[end])  # RTN impulse [kg·m/s]

    ΔE_total = 0.0
    if !isnothing(sa_sol)
        # ΔE ≈ Σ m·Δdv · v_target_RTN
        for k in 2:nT
            δdv_R = dv_R[k] - dv_R[k-1]
            δdv_T = dv_T[k] - dv_T[k-1]
            δdv_N = dv_N[k] - dv_N[k-1]
            δP_R = mass * δdv_R; δP_T = mass * δdv_T; δP_N = mass * δdv_N
            t_k  = t_hist[k]
            k_sa = clamp(searchsortedlast(Float64.(sa_sol.t), t_k), 1, length(sa_sol.t))
            sc   = sa_sol.u[k_sa].sc[1]
            r_vec = SVector{3,Float64}(sc.pos); v_vec = SVector{3,Float64}(sc.vel)
            rhat, that, nhat = _rtn_basis(r_vec, v_vec)
            ΔE_total += δP_R*dot(v_vec,rhat) + δP_T*dot(v_vec,that) + δP_N*dot(v_vec,nhat)
        end
    end

    key = (0, 1)
    return Dict(key => ΔP_total), Dict(key => ΔE_total)
end

# Time-series of cumulative ΔP (RTN impulse vector) and ΔE delivered by the laser ensemble.
# Returns (ΔP_series, ΔE_series, t) matching the signature expected by 7_Plots.jl.
# sol_or_sa: may be _FlatSol (called from 7_Plots.jl) or the raw SpaceAGORA sol;
#            SpaceAGORA structured sol is also read from p[:sa_sol] when needed.
function laser_exchanges_time_series(sol_or_sa, p::Dict)
    tracker = get(p, :tracker, nothing)
    # Resolve structured sol for velocity energy computation
    sa_sol = sol_or_sa isa _FlatSol ? get(p, :sa_sol, nothing) : sol_or_sa

    if tracker === nothing || isempty(tracker.t_hist)
        t_fb = sol_or_sa isa _FlatSol ? sol_or_sa.t : Float64.(sol_or_sa.t)
        key  = (0,1)
        return Dict(key => [SVector(0.0,0.0,0.0) for _ in t_fb]),
               Dict(key => zeros(length(t_fb))), t_fb
    end
    mass   = p[:masses][1]
    t_hist = tracker.t_hist
    dv_R   = tracker.dv_R_hist
    dv_T   = tracker.dv_T_hist
    dv_N   = tracker.dv_N_hist
    nT     = length(t_hist)

    ΔP_series = [SVector(mass*dv_R[k], mass*dv_T[k], mass*dv_N[k]) for k in 1:nT]
    ΔE_series = zeros(nT)
    if !isnothing(sa_sol)
        Ecum = 0.0
        for k in 2:nT
            δdv_R = dv_R[k] - dv_R[k-1]
            δdv_T = dv_T[k] - dv_T[k-1]
            δdv_N = dv_N[k] - dv_N[k-1]
            δP_R = mass*δdv_R; δP_T = mass*δdv_T; δP_N = mass*δdv_N
            t_k  = t_hist[k]
            k_sa = clamp(searchsortedlast(Float64.(sa_sol.t), t_k), 1, length(sa_sol.t))
            sc   = sa_sol.u[k_sa].sc[1]
            r_vec = SVector{3,Float64}(sc.pos); v_vec = SVector{3,Float64}(sc.vel)
            rhat, that, nhat = _rtn_basis(r_vec, v_vec)
            Ecum += δP_R*dot(v_vec,rhat) + δP_T*dot(v_vec,that) + δP_N*dot(v_vec,nhat)
            ΔE_series[k] = Ecum
        end
    end

    key = (0,1)
    return Dict(key => ΔP_series), Dict(key => ΔE_series), t_hist
end

# =============================================================================
# Case-ID / CSV / print helpers (runner-level, unchanged)
# =============================================================================

# Generate a unique case ID string based on the options (used for labeling output)
function _case_id(opts::OracleCase2Options)::String
    return @sprintf(
        "Nh%d_h%.0fkm_ti%.1fdeg_hi%.1fdeg",
        opts.helpers,
        opts.helper_altitude_km,
        opts.target_inclination_deg,
        opts.helper_inclination_deg,
    )
end

# Write a DataFrame to CSV, creating the directory if needed, and optionally appending to an existing file.
function _write_csv!(path::String, data; append::Bool=false)
    isempty(path) && return ""
    mkpath(dirname(path))
    do_append = append && isfile(path)
    df = data isa DataFrame ? data : DataFrame(data)
    CSV.write(path, df; append=do_append, header=!do_append)
    return path
end

# Write a DataFrame to Feather (Arrow) format, creating the directory if needed.
function _write_feather!(path::String, data)
    isempty(path) && return ""
    mkpath(dirname(path))
    df = data isa DataFrame ? data : DataFrame(data)
    Arrow.write(path, df)
    return path
end

# Build a descriptive filename stem for one simulation scenario.
# Example: "N4_T504560s_h1050km_t1000km_ih0.0deg_it0.0deg"
function _scenario_stem(s, opts::OracleCase2Options)::String
    T_s = round(Int, s.target_period_s * opts.orbits)
    return @sprintf(
        "N%d_T%ds_h%.0fkm_t%.0fkm_ih%.1fdeg_it%.1fdeg",
        opts.helpers,
        T_s,
        opts.helper_altitude_km,
        opts.target_altitude_km,
        opts.helper_inclination_deg,
        opts.target_inclination_deg,
    )
end

# Print a summary of the simulation results to the console.
function _print_summary(s)
    println("ORACLE Open Cavity Case laser-link run")
    println("  helpers              : $(s.helpers)")
    println("  helper altitude       : $(s.helper_altitude_km) km")
    println("  schedule              : $(s.schedule)")
    println("  solver                : $(s.solver)")
    println("  retcode               : $(s.retcode)")
    @printf("  mission time          : %.2f s (%.4f orbits elapsed)\n", s.orbits * s.target_period_s, s.orbits_elapsed)
    println("  link activations      : $(s.activations)")
    println("  accepted active steps : $(s.active_steps)")
    @printf("  dV_RTN [m/s]          : R=% .6e  T=% .6e  N=% .6e\n", s.dv_r_mps, s.dv_t_mps, s.dv_n_mps)
    @printf("  da [m]                : % .6e\n", s.da_m)
    @printf("  de                    : % .6e\n", s.de)
    @printf("  di [deg]              : % .6e\n", s.di_deg)
    @printf("  dRAAN [deg]           : % .6e\n", s.draan_deg)
end
