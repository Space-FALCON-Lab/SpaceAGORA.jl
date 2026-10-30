# Planner arms (WP-P2; plan: tasks/2026-10-15_periapsis-region-border-targeting/analytic_framework_and_transfer_plan.md R4.5 to R4.10, R7.3, R7.4, R8.1, R8.2; config planner_config.toml).
# ONE planner code path with switches (slack, drift feed-forward, early burn timing, geometry) gives the ladder L0 to L3; the same path flown without feedback (OL) or re-solved at every apoapsis over the whole
# remaining campaign (OL-R) gives the two autonomy-ablation arms; B-ODY (plan section 7, U10) is the Odyssey-style reactive keeper. NOT RUN: written without executing Julia; the selftest (planner_run.jl --phase=selftest)
# is the first thing the Julia runner must execute (Freeze B items F10 to F13).
# Reuses stageg3_common.jl (included, not edited; through it stageo_common.jl, stagen2_common.jl and the Stage F/H files, the estimator v5 and the SafetyStop / U15 helpers). New code only:
#   (i)   pure functions, mirrored line for line in planner_lib.py (cross-checked by the selftest): the exact lever, the lever weight w = H / (ratio x lever) in m/s per ln-q unit, the exact dynamic program over cumulative ln-q
#         offsets (planner_dp / planner_plan), the plan inputs per arm (planner_inputs, planner_ol_schedule), the B-ODY rule (planner_body_dv);
#   (ii)  the ground table (planner_ground_tables.json, sha256 pinned in the config): per geometry the model's burn-free ln-q path, its 15-pass slow part, the lever weights and the apsides, written by planner_predictions.py
#         from geometry_rank.py (the validated Sun-tide closed form, zonal, tesseral and drag channels of analytic_results.md). It is the planner's "predicted drift";
#   (iii) PlannerPolicy: a guidance wrapper that decides the tangential apoapsis burn of each orbit. Coordinates: ln q of the mean-field heat rate at the periapsis event on V_rel (U15). The next pass is predicted with the
#         frozen estimator v5 (stagef_predict_alt via the keeper's make_predictor), the density offset by the P2_10 EWMA of the mean-field offset; passes beyond it by the table's slow increments (feed-forward);
#   (iv)  planner_run_task: one flight (the structure of stageg3_run_task with the planner policy; the apsis callback, the peak callback, the SafetyStop and the output columns are stageg3's);
#   (v)   the Odyssey reference geometry (ODY) injected into the gate geometry list so stageg3_load / stageg3_cal_task serve it (same elements as stageo_config.toml, so the start is the Stage O start).
# Honest limits (also in planner_build.md): superposition of burns is checked for single burns only (R4.1); the table is a model (4 to 5% long in time, field term only where sampled); the heavy v5 predictor is called once
# per orbit (not for OL); no src/ edit; no Julia has been run on this file.
using Dates, TOML, JSON, SHA
include(joinpath(@__DIR__, "stageg3_common.jl"))

# ---- pure core (selftested; mirror of planner_lib.py) -----------------------------------------------------------------------------------------------

"""dr_p / dv for a tangential apoapsis burn (m per m/s), exact vis-viva: 4 a sqrt((1-e)/(1+e)) / sqrt(mu/a)."""
planner_lever_m_per_mps(mu::Real, a::Real, e::Real) = 4a * sqrt((1 - e) / (1 + e)) / sqrt(mu / a)
"""Burn cost per unit of ln q (m/s per ln unit): H / (ratio x lever); `ratio` is burn_response_ratio (1.0, DERIVED from the twin test)."""
planner_weight(H::Real, ratio::Real, mu::Real, a::Real, e::Real) = H / (ratio * planner_lever_m_per_mps(mu, a, e))

"""Exact dynamic program: minimize sum_k L[k] |nu_k - nu_{k-1}| (nu_0 = 0) subject to A[k] <= nu_k <= B[k]. State set {0} U {A} U {B} (an optimal vertex of the piecewise-linear problem). Two relaxation sweeps per step; strict `<` so
ties keep the nearer source (identical to planner_lib.dp_exact). `pen` (m/s) is a fixed cost per burn (a step with nu_k != nu_{k-1}); 0 = none. `lin` (Revision 19 A5, default nothing) is an optional linear cost
lin[k] * nu_k on the level of step k (added after the relaxation, as planner_rev18_offline.dp_lin); a vector with no nonzero entry takes exactly the code path of `lin = nothing` (bit-identical). Returns (nu, obj) or nothing when infeasible."""
function planner_dp(A::AbstractVector{<:Real}, B::AbstractVector{<:Real}, L::AbstractVector{<:Real}; tol::Real=1.0e-9, pen::Real=0.0, lin::Union{Nothing, AbstractVector{<:Real}}=nothing)
    n = length(A); (length(B) == n && length(L) == n) || throw(ArgumentError("A, B, L must have the same length"))
    haslin = lin !== nothing && any(!iszero, lin); haslin && length(lin) != n && throw(ArgumentError("lin must have the same length as A"))
    Af = Float64.(A) .+ 0.0; Bf = Float64.(B) .+ 0.0                      # + 0.0 turns -0.0 into 0.0 (unique() would keep both)
    val = sort(unique(vcat(0.0, Af, Bf))); ns = length(val); i0 = findfirst(==(0.0), val)
    prev = fill(Inf, ns); prev[i0] = 0.0; cur = similar(prev); arg = Vector{Int32}(undef, ns); parent = Matrix{Int32}(undef, n, ns)
    for k in 1:n
        Lk = Float64(L[k]); copyto!(cur, prev); for j in 1:ns; arg[j] = j; end
        for j in 2:ns
            c = cur[j - 1] + Lk * (val[j] - val[j - 1])
            if c < cur[j]; cur[j] = c; arg[j] = arg[j - 1]; end
        end
        for j in (ns - 1):-1:1
            c = cur[j + 1] + Lk * (val[j + 1] - val[j])
            if c < cur[j]; cur[j] = c; arg[j] = arg[j + 1]; end
        end
        if pen > 0.0                                             # staying is free, any move pays pen once: new = min(prev, pen + best move)
            for j in 1:ns
                if prev[j] <= pen + cur[j]; cur[j] = prev[j]; arg[j] = j; else; cur[j] = pen + cur[j]; end
            end
        end
        if haslin; lk = Float64(lin[k]); for j in 1:ns; cur[j] += lk * val[j]; end; end
        for j in 1:ns
            (val[j] < Af[k] - tol || val[j] > Bf[k] + tol) && (cur[j] = Inf)
        end
        parent[k, :] = arg; copyto!(prev, cur)
    end
    j = argmin(prev); isfinite(prev[j]) || return nothing
    nu = zeros(n); obj = prev[j]
    for k in n:-1:1; nu[k] = val[j]; j = Int(parent[k, j]); end
    return (; nu, obj)
end

"""Lazy tie-break weights: step m (0-based) of n costs L_m (1 + eps (n - m) / n)."""
planner_tie_weights(L::AbstractVector{<:Real}, eps::Real) = (n = length(L); Float64.(L) .* (1 .+ eps .* (n .- (0:(n - 1))) ./ n))

"""Offsets nu (cumulative ln-q) with level + nu inside [lo, hi]; (nu, du, dv_cost) or nothing. dv_cost uses the untouched weights."""
function planner_plan(level, lo, hi, w; eps::Real=1.0e-6, tol::Real=1.0e-9, pen::Real=0.0, lin=nothing)
    A = Float64.(lo) .- Float64.(level); B = Float64.(hi) .- Float64.(level)
    r = lin === nothing ? planner_dp(A, B, planner_tie_weights(w, eps); tol=tol, pen=pen) : planner_dp(A, B, planner_tie_weights(w, eps); tol=tol, pen=pen, lin=Float64.(lin)); r === nothing && return nothing
    du = diff(vcat(0.0, r.nu))
    return (; nu=r.nu, du, dv_cost=sum(Float64.(w) .* abs.(du)))
end

"""Table value at 1-based pass p, clamped to the table (the last values are held beyond it)."""
planner_tget(arr::AbstractVector, p::Integer) = Float64(arr[clamp(p, 1, length(arr))])
"""Slow ln-q increment from pass p to p + 1 (0 outside the table)."""
planner_dslow(tab, p::Integer) = (1 <= p <= length(tab.d_slow)) ? Float64(tab.d_slow[p]) : 0.0

"""DP inputs for the decision before pass j (modes online and olr). sw = (slack, ff, early, horizon); tab = ground table; P = planner constants; z0 = predicted ln q of pass j without a further burn (full);
w0 = live weight of this burn; bias = OL-R drift-update slope (ln q per pass). FF off: persistence of z0, one step (the DP is a clamp to the nearest tube edge). FF on: z0, then the slow anchor z0 - fast_j plus the table's slow
increments (plus gain x bias for the next cap_passes steps) inside the tube shrunk by half the fast reserve at each edge."""
function planner_inputs(sw, tab, P, j::Integer, z0::Real, w0::Real; bias::Real=0.0, lower_margin::Real=0.0, fb=nothing, fb_N::Integer=0, mu::Real=NaN)
    hi0 = log(P.q_cap_mean_w_cm2); W = sw.slack ? tab.W : P.l0_width_ln; lo0 = hi0 - W + lower_margin
    sw.ff || return (; level=[Float64(z0)], lo=[lo0], hi=[hi0], w=[Float64(w0)])
    n = sw.horizon > 0 ? Int(sw.horizon) : max(1, tab.n - j + 1 + Int(P.tail_passes))
    res = sw.slack ? P.fast_reserve_ln : 0.0
    s_j = z0 - planner_tget(tab.fast, j); level = Vector{Float64}(undef, n); level[1] = z0; acc = 0.0
    cap = Int(P.update_cap_passes); gain = Float64(P.update_gain)
    for m in 1:(n - 1)
        acc += planner_dslow(tab, j + m - 1) + (m <= cap ? gain * bias : 0.0)
        level[m + 1] = s_j + acc
    end
    lo = vcat(lo0, fill(lo0 + res / 2, n - 1)); hi = vcat(hi0, fill(hi0 - res / 2, n - 1))
    w = sw.early ? vcat(Float64(w0), [planner_tget(tab.w, j + m) for m in 1:(n - 1)]) : fill(Float64(w0), n)
    fb === nothing && return (; level, lo, hi, w)
    return (; level, lo, hi, w, lin=planner_term_lin(tab, fb, P, mu, P.ra_target_m, j, fb_N, n))          # Revision 19 A5: a linear terminal-periapsis cost at the arrival pass fb_N (nothing when it lies outside the plan)
end

"""Open-loop schedule (arm OL): passes p = 2 .. n + tail, level = ln q_first + the table's slow increments, tube shrunk by half the fast reserve at both edges, the table's true weights. Returns Dict pass => signed dv
(m/s, + raises the periapsis); bursts below min_burn are dropped. Computed once, before the first decision, from the table and the calibrated first-pass q only (no measurement). Terminal-periapsis term: planner_ol_plan."""
planner_ol_schedule(tab, P, q_first_ln::Real) = planner_ol_plan(tab, P, q_first_ln).sched

# ---- terminal-periapsis term (Revision 19 A5; default off; mirror of planner_lib.py) -----------------------------------------------------------------

"""Exact vis-viva tangential periapsis-raise burn at apoapsis ra from periapsis rp to rp_des (m/s): v_a(a') - v_a(a), v_a = sqrt(mu (2/ra - 1/a)); negative when it lowers the periapsis."""
planner_first_burn_dv(mu::Real, ra::Real, rp::Real, rp_des::Real) = sqrt(mu * (2 / ra - 2 / (ra + rp_des))) - sqrt(mu * (2 / ra - 2 / (ra + rp)))
"""c_T = H dFB/drp (m/s per ln-q unit): price of a unit of final level at the closing apoapsis (ra_f, rp_f) with the desired periapsis radius fixed. dFB/drp = 1 / lever(a_f, e_f) exactly (lever = 4 a^2 v_a / mu = d rp / d v_a),
so the term needs neither the desired altitude nor a finite difference."""
planner_first_burn_slope(mu::Real, ra_f::Real, rp_f::Real, H::Real) = H / planner_lever_m_per_mps(mu, (ra_f + rp_f) / 2, (ra_f - rp_f) / (ra_f + rp_f))

"""Per-step linear cost vector of length n for a DP whose step 1 is the burn before pass j: fb.weight x c_T at the step of the arrival pass N (step N - j + 1), zero elsewhere; c_T from the table periapsis and scale height at
pass min(N, tab.n) and the closing apoapsis ra_target. nothing when that step lies outside 1..n (no self-consistent arrival inside the plan: the term is dropped)."""
function planner_fb_lin(tab, fb, mu::Real, ra_target::Real, j::Integer, N::Integer, n::Integer)
    m = N - j + 1; (1 <= m <= n) || return nothing
    kk = clamp(N, 1, tab.n); rp_t = tab.a[kk] * (1 - tab.e[kk]); lin = zeros(n); lin[m] = fb.weight * planner_first_burn_slope(mu, ra_target, rp_t, tab.H[kk])
    return lin
end

"""Linear cost (m/s per ln-q unit) of the shortfall of the closing apoapsis below the target circle (Revision 22 item 3, switch apo_one_sided): shortfall = target - closing apoapsis grows with the last drop, d drop / d ln q = drop,
so c = slope (m/s per km of shortfall) x drop / 1000 (m -> km)."""
planner_apo_slope(slope_mps_per_km::Real, drop_m::Real) = slope_mps_per_km * drop_m / 1000.0
"""Linear cost vector of the terminal terms at the arrival pass N of a DP whose step 1 is before pass j: the FB term (planner_fb_lin) plus, when P.apo_one_sided, the shortfall term with the table's drop at pass min(N, n) (drop at q_req,
ASSUMED to stand for the last drop). nothing when the step lies outside the plan."""
function planner_term_lin(tab, fb, P, mu::Real, ra_target::Real, j::Integer, N::Integer, n::Integer)
    lin = planner_fb_lin(tab, fb, mu, ra_target, j, N, n); lin === nothing && return nothing
    P.apo_one_sided || return lin
    kk = clamp(N, 1, tab.n); ra_tab = tab.a .* (1 .+ tab.e); drop = kk < tab.n ? ra_tab[kk] - ra_tab[kk + 1] : ra_tab[tab.n - 1] - ra_tab[tab.n]
    lin[N - j + 1] += planner_apo_slope(P.apo_shortfall_slope_mps_per_km, drop); return lin
end

"""1-D linear interpolation clamped at the ends (numpy.interp for ascending xp)."""
function planner_interp1(x::Real, xp::AbstractVector{<:Real}, fp::AbstractVector{<:Real})
    x <= xp[1] && return Float64(fp[1]); x >= xp[end] && return Float64(fp[end])
    i = searchsortedlast(xp, x); t = (x - xp[i]) / (xp[i + 1] - xp[i]); return Float64(fp[i]) + t * (Float64(fp[i + 1]) - Float64(fp[i]))
end

"""In-model closing-apoapsis pass and time of a planned ln-q path x (x[k] = ln q of pass k, held at its last value beyond its length): the drop per unit q is read at the CURRENT apoapsis radius (table drop interpolated on ra),
scaled by kappa (MEASURED per geometry), the period is Kepler on the osculating a; the burn-induced periapsis shift is -H nu (Revision 18b / 19 B; planner_rev18_offline.arrival). Returns (N, days, ra_end); N = 3 n when the target is not reached."""
function planner_arrival(tab, mu::Real, ra_target::Real, x::AbstractVector{<:Real}, q_first_ln::Real, kappa::Real)
    n = tab.n; ra_tab = tab.a .* (1 .+ tab.e); rp_tab = tab.a .* (1 .- tab.e); ra = ra_tab[1]; t = tab.t[1]
    d0 = vcat(ra_tab[1:(end - 1)] .- ra_tab[2:end], ra_tab[end - 1] - ra_tab[end]); ra_r = reverse(ra_tab); d0_r = reverse(d0); kmax = 3 * n
    for k in 1:kmax
        kk = min(k, n); xk = Float64(x[min(k, length(x))]); nu = xk - (q_first_ln + (tab.path[kk] - tab.path[1])); rp = rp_tab[kk] - tab.H[kk] * nu; a = 0.5 * (ra + rp)
        drop = kappa * planner_interp1(ra, ra_r, d0_r) * exp(xk) / tab.q_req; ra_next = ra - drop; amid = 0.5 * (a + 0.5 * (ra_next + rp)); Tk = 2π * sqrt(amid^3 / mu)
        ra_next <= ra_target && return (; N=k, days=(t + 0.5 * Tk) / 86400.0, ra_end=ra_next)
        ra = ra_next; t += Tk
    end
    return (; N=kmax, days=NaN, ra_end=ra)
end

"""Open-loop plan on the ground table. fb = nothing: the existing schedule (same DP call as before the term existed; N = 0). fb = (; mu, kappa, weight): the terminal-periapsis term of Revision 19 as in planner_rev19_offline.model_plan: solve
without the term, find the arrival pass N of that plan (planner_arrival), re-solve with the linear cost fb.weight x c_T at the step of pass N, and repeat until N is self-consistent (at most `iters` re-solves). When no arrival pass lies
inside the plan the term is dropped (dropped = true) and the plan without it is returned. Returns (sched, N, dropped, converged, nu)."""
function planner_ol_plan(tab, P, q_first_ln::Real; fb=nothing, iters::Integer=3)
    hi0 = log(P.q_cap_mean_w_cm2); lo0 = hi0 - tab.W; res = P.fast_reserve_ln
    N = tab.n + Int(P.tail_passes); passes = collect(2:N); level = Float64[]; acc = 0.0
    for p in passes; acc += planner_dslow(tab, p - 1); push!(level, q_first_ln + acc); end
    lo = fill(lo0 + res / 2, length(passes)); hi = fill(hi0 - res / 2, length(passes)); w = [planner_tget(tab.w, p) for p in passes]
    plain() = planner_plan(level, lo, hi, w; eps=P.dp_tie_eps, pen=P.dp_burn_penalty_mps)
    pl = nothing; Narr = 0; dropped = false; converged = true
    if fb === nothing
        pl = plain()
    else
        lin = nothing
        for it in 0:iters
            pl = planner_plan(level, lo, hi, w; eps=P.dp_tie_eps, pen=P.dp_burn_penalty_mps, lin=lin); pl === nothing && break
            Nn = planner_arrival(tab, fb.mu, P.ra_stop_m, vcat(q_first_ln, level .+ pl.nu), q_first_ln, fb.kappa).N
            Nn == Narr && break
            Narr = Nn; lin = planner_term_lin(tab, fb, P, fb.mu, P.ra_target_m, 2, Narr, length(level))
            if lin === nothing; dropped = true; pl = plain(); break; end
            it == iters && (converged = false)
        end
    end
    sched = Dict{Int, Float64}(); pl === nothing && return (; sched, N=Narr, dropped, converged, nu=Float64[], level)
    for (i, p) in enumerate(passes)
        abs(pl.du[i]) * w[i] >= P.dp_min_burn_mps && (sched[p] = -pl.du[i] * w[i])
    end
    return (; sched, N=Narr, dropped, converged, nu=pl.nu, level)
end

"""B-ODY: the next pass's Q_s is predicted as exp(z0) x qs_over_qmean; burn (re-center on the phase target) only when it is outside [q_low, q_high] and the last burn is at least cadence_s ago."""
function planner_body_dv(B, P, j::Integer, t::Real, last_burn_t::Real, z0::Real, w0::Real)
    qs = exp(z0) * B.qs_over_qmean
    (B.q_low_qs_w_cm2 <= qs <= B.q_high_qs_w_cm2 || (t - last_burn_t) < B.cadence_s) && return 0.0
    tgt = j <= B.mp1_last_pass ? B.q_target_mp1_qs_w_cm2 : B.q_target_mp2_qs_w_cm2
    dv = -(log(tgt / B.qs_over_qmean) - z0) * w0
    return abs(dv) >= P.dp_min_burn_mps ? dv : 0.0
end

# ---- terminal aim (Revision 14; only arms with the terminal switch) ----------------------------------------------------------------------------------

planner_median(v::AbstractVector{<:Real}) = (s = sort(Float64.(v)); n = length(s); isodd(n) ? s[(n + 1) ÷ 2] : 0.5 * (s[n ÷ 2] + s[n ÷ 2 + 1]))
"""Terminal aim: R = osculating apoapsis radius at the decision minus the target (m); c_hist = observed (apoapsis drop of a pass) / (its measured q_mean) in m per W/cm2. Active when at least 3 samples exist and
0 < R + aim <= activate_steps x step. Splits R + aim into N = ceil((R + aim) / step) equal drops d*, returns the ln q whose predicted drop (c x q) equals d*, clipped to [ln(min_factor x q_cap), ln q_cap]. nothing when inactive."""
function planner_terminal_aim(R::Real, c_hist::AbstractVector{<:Real}, P)
    length(c_hist) < 3 && return nothing
    Reff = R + P.term_aim_below_m; Reff <= 0 && return nothing; Reff > P.term_activate_steps * P.term_step_m && return nothing
    K = min(length(c_hist), P.term_c_window); c = planner_median(c_hist[(end - K + 1):end])
    N = max(1, ceil(Int, Reff / P.term_step_m)); dstar = Reff / N; hi = log(P.q_cap_mean_w_cm2); lo = hi + log(P.term_min_q_factor)
    return (; lnq=clamp(log(dstar / c), lo, hi), N, dstar, c)
end

# ---- ground table ----------------------------------------------------------------------------------------------------------------------------------

"""Ground table of one geometry id from planner_ground_tables.json (path and sha256 from the config; a mismatch is an ERROR; an empty sha256 is an ERROR unless allow_unhashed). Fields: n, q_req, W, path, slow, d_slow,
fast, w, H, a, e, t (arrays by 1-based pass)."""
function planner_table(meta, gid::AbstractString; allow_unhashed::Bool=false, file::Union{Nothing, AbstractString}=nothing, section::AbstractString="ground_table")
    pt = meta["planner"][section]; fp = file === nothing ? joinpath(stageg3_task_dir(), String(pt["path"])) : String(file); raw = read(fp); h = bytes2hex(sha256(raw)); want = String(pt["sha256"])
    if file === nothing && !allow_unhashed
        isempty(want) && error("planner_config.toml [planner.ground_table] sha256 is empty: run planner_predictions.py --tables, record the hash in the config, then fly")
        h == want || error("ground table sha256 mismatch: $fp is $h, config has $want")
    end
    T = JSON.parse(String(raw))
    if !haskey(T["geometries"], gid) && file === nothing && section == "ground_table" && haskey(meta["planner"], "ground_table_generic")      # Revision 16: the generic orbits live in their own pinned file
        return planner_table(meta, gid; allow_unhashed=allow_unhashed, section="ground_table_generic")
    end
    if !haskey(T["geometries"], gid) && file === nothing && section == "ground_table_generic" && haskey(meta["planner"], "ground_table_rev23")      # Revision 23b: the mean-field tables (keys G0528, G0890, G0975M, G0063M) live in their own pinned file
        return planner_table(meta, gid; allow_unhashed=allow_unhashed, section="ground_table_rev23")
    end
    if !haskey(T["geometries"], gid) && file === nothing && section == "ground_table_rev23" && haskey(meta["planner"], "ground_table_rev25")      # Revision 25: the tables of the low-altitude targets (G0486, G1230, G1038)
        return planner_table(meta, gid; allow_unhashed=allow_unhashed, section="ground_table_rev25")
    end
    haskey(T["geometries"], gid) || error("ground table has no geometry $gid")
    g = T["geometries"][gid]; f = k -> Float64.(g[k])
    return (; n=Int(g["n"]), q_req=Float64(g["q_req"]), W=Float64(g["W"]), path=f("path"), slow=f("slow"), d_slow=f("d_slow"), fast=f("fast"), w=f("w"), H=f("H"), a=f("a"), e=f("e"), t=f("t"))
end

# ---- Revision 20 guard v2 (time-to-go; switch ttg_guard, default off; mirror of planner_rev20_replay.py guard_solve with GUARD_V2; spec results_rev20_guard_v2.md section 5) -------------------------------

"""Guard constants ([planner.guard], all ASSUMED and frozen in Revision 20c): (; pen_g, rel_slack, s_frac, aim_frac, anchor_k, anchor_from_pass, bias_steps, k_window, ema_beta, eps_max, max_solves, dev_thr_g, start_pass,
w_max_periods, w_rp_m, t_offset_d (Dict gid => days)). The actuator (tube), the infeasible mode (pin), the kappa statistic (median) and the edge margin (0) are fixed in the code."""
function planner_guard_constants(meta; set::AbstractString="v2", target::Union{Nothing, AbstractString}=nothing)
    g = set == "v2" ? meta["planner"]["guard"] : meta["planner"]["guard_sets"][set]; f = k -> Float64(g[k])
    node_check = Bool(get(g, "node_check", false)); node_tol = get(g, "node_tol_deg", nothing)
    node_check && node_tol === nothing && throw(ErrorException("guard node_check is enabled but node_tol_deg has no value: no tolerance has been given (Revision 25); set node_tol_deg explicitly"))
    w_rp = Float64(get(planner_target_table(meta, target), "guard_w_rp_m", g["w_rp_m"]))         # "v2" = [planner.guard] (the default); other sets: [planner.guard_sets.<name>] (same keys)
    return (; pen_g=f("pen_g"), rel_slack=f("rel_slack"), s_frac=f("s_frac"), aim_frac=f("aim_frac"), anchor_k=Int(g["anchor_k"]), anchor_from_pass=Int(g["anchor_from_pass"]), bias_steps=Int(g["bias_steps"]),
        k_window=Int(g["k_window"]), ema_beta=f("ema_beta"), eps_max=f("eps_max"), max_solves=Int(g["max_solves"]), dev_thr_g=f("dev_thr_g"), start_pass=Int(g["start_pass"]), w_max_periods=f("w_max_periods"),
        w_rp_m=w_rp, node_check, node_tol_deg=node_tol === nothing ? NaN : Float64(node_tol), term_model=Bool(get(g, "term_model", false)), t_offset_d=Dict{String, Float64}(String(k) => Float64(v) for (k, v) in g["t_offset_d"]))
end

"""Kepler period (s) of the orbit with apoapsis ra and periapsis rp."""
planner_kepler_period(mu::Real, ra::Real, rp::Real) = (a = 0.5 * (ra + rp); 2π * sqrt(a^3 / mu))
"""Measured drop / model drop (state dependent, unit kappa) of the last k passes: median over passes p = max(2, n - k + 1) .. n of (ras[p-1] - ras[p]) / (d0(ras[p-1]) q_p / q_req); nothing for n < 4 or fewer than 3 usable values.
xs[p] = realized ln q of pass p, ras[p] = apoapsis radius after pass p."""
function planner_ttg_kappa_hat(xs::AbstractVector{<:Real}, ras::AbstractVector{<:Real}, tab, k::Integer)
    n = length(xs); n < 4 && return nothing
    ra_tab = tab.a .* (1 .+ tab.e); d0 = vcat(ra_tab[1:(end - 1)] .- ra_tab[2:end], ra_tab[end - 1] - ra_tab[end]); ra_r = reverse(ra_tab); d0_r = reverse(d0); vals = Float64[]
    for p in max(2, n - k + 1):n
        r_ = ras[p - 1]; dr = r_ - ras[p]; m = planner_interp1(r_, ra_r, d0_r) * exp(xs[p]) / tab.q_req
        (dr > 0 && m > 0) && push!(vals, dr / m)
    end
    return length(vals) >= 3 ? planner_median(vals) : nothing
end

"""Anchored level forecast of the guard at decision j (nothing before pass anchor_from_pass): lvl_alt[1] = lvl[1], lvl_alt[m+1] = u + cbar + slow(j + m), cbar = mean over the last min(anchor_k, j - 2) passes of
xs_i - (u0 + nuh_i) - path_i (path and slow of the table by pass index, clamped), u = executed offset now including the initial offset u0."""
function planner_ttg_lvl_alt(G, tab, j::Integer, lvl::AbstractVector{<:Real}, xs, nuh, u0::Real, u::Real)
    (G.anchor_k > 0 && j >= G.anchor_from_pass) || return nothing
    Kk = min(G.anchor_k, j - 2); ii = (j - Kk):(j - 1)
    cbar = sum(xs[q] - (u0 + nuh[q]) - planner_tget(tab.path, q) for q in ii) / length(ii)
    alt = Float64.(lvl); for m in 1:(length(lvl) - 1); alt[m + 1] = u + cbar + planner_tget(tab.slow, j + m); end
    return alt
end

"""Wrap an angle in degrees to (-180, 180]."""
planner_wrap_deg(x::Real) = (y = mod(Float64(x) + 180.0, 360.0) - 180.0; y == -180.0 ? 180.0 : y)
"""Secular J2 rate of the ascending node (rad/s): dOmega/dt = -1.5 n J2 (Re / p)^2 cos i, n = sqrt(mu / a^3), p = a (1 - e^2) (Kaula; as planner_rev21_position.j2_orbit per orbit)."""
planner_j2_node_rate(mu::Real, J2::Real, Re::Real, a::Real, e::Real, inc_rad::Real) = -1.5 * sqrt(mu / a^3) * J2 * (Re / (a * (1 - e^2)))^2 * cos(inc_rad)
"""Node error (deg) at the arrival time t_arr: wrap(Omega_now + rate (t_arr - t_now) - Omega_target); rate in deg/s."""
planner_node_error_deg(om_now_deg::Real, rate_deg_s::Real, t_now::Real, t_arr::Real, om_target_deg::Real) = planner_wrap_deg(om_now_deg + rate_deg_s * (t_arr - t_now) - om_target_deg)
"""Arrival time (s) at which the node error is zero, the instance of the node cycle (period 360 / |rate|) nearest to tt."""
function planner_node_aim(om_now_deg::Real, rate_deg_s::Real, t_now::Real, om_target_deg::Real, tt::Real)
    rate_deg_s == 0.0 && return Float64(tt)
    t0 = t_now + planner_wrap_deg(om_target_deg - om_now_deg) / rate_deg_s; per = 360.0 / abs(rate_deg_s)
    return t0 + per * round((tt - t0) / per)
end

"""Closing-apoapsis time (s) predicted from the state after pass j - 1 (apoapsis ra0, periapsis rp0, periapsis time t0) for planned ln q x[m+1] and periapsis rp[m+1] of pass j + m (held at the last entry), the drop per unit q
kappa d0(ra) exp(x) / q_req, Kepler periods. Returns (t, steps) with t = NaN after 400 steps (planner_rev20_replay.predict_arrival)."""
function planner_ttg_predict(tab, mu::Real, ra_target::Real, kap::Real, ra0::Real, rp0::Real, t0::Real, x::AbstractVector{<:Real}, rpp::AbstractVector{<:Real}; term=nothing, ra_aim=nothing)
    ra_tab = tab.a .* (1 .+ tab.e); d0 = vcat(ra_tab[1:(end - 1)] .- ra_tab[2:end], ra_tab[end - 1] - ra_tab[end]); ra_r = reverse(ra_tab); d0_r = reverse(d0)
    ra = Float64(ra0); rp_prev = Float64(rp0); t = Float64(t0); nx = length(x)
    for m in 0:399
        xm = x[min(m + 1, nx)]; rpm = rpp[min(m + 1, nx)]
        t += 0.5 * planner_kepler_period(mu, ra, rp_prev) + 0.5 * planner_kepler_period(mu, ra, rpm)
        dr = kap * planner_interp1(ra, ra_r, d0_r) * exp(xm) / tab.q_req
        if term !== nothing                                                                   # guard v3 (term_model): the terminal aim takes over inside activate_steps x step of the target: N equal drops of (R + aim) / N
            Reff = ra - (ra_aim === nothing ? ra_target : ra_aim) + term.below
            (0 < Reff <= term.steps * term.step) && (dr = Reff / max(1, ceil(Int, Reff / term.step)))
        end
        ra -= dr
        ra <= ra_target && return (t + 0.5 * planner_kepler_period(mu, ra, rpm), m + 1)
        rp_prev = rpm
    end
    return (NaN, 400)
end

"""Guard core (planner_rev20_replay.guard_solve with GUARD_V2). Inputs of decision j: the planner's plan inputs (lvl, A_ = lo - lvl, B_ = hi - lvl, L_ = tie weights) and its DP offsets nu0, the state after pass j - 1 (ra_now, rp_now,
t_now = time of that periapsis), the target time tstar_s (aim offset included), kap_meas (kappa_hat or nothing), the drift slope `bias`, the local scale height, the anchored forecast lvl_alt (or nothing), and the smoothed arrival
of the previous decision (NaN = none). Returns (; nu, engaged, infeasible, ema, ta_raw, ta, s, steps, solves, late): nu = nu0 while the predicted arrival is inside [t* - s, t* + s]; otherwise the tube is tightened on the side
the arrival needs (the upper edge is never raised) by the bisected fraction u of the tube width, at most max_solves DP solves in all."""
function planner_ttg_decide(G, tab, P, mu::Real, kap::Real, j::Integer, lvl, A_, B_, L_, nu0, ra_now::Real, rp_now::Real, t_now::Real, tstar_s::Real, kap_meas, bias::Real, Hloc::Real, lvl_alt, ema_prev::Real; node=nothing)
    kh = kap_meas === nothing ? kap : kap_meas; n = length(lvl); cap = Int(P.update_cap_passes); gain = Float64(P.update_gain); K2 = max(G.bias_steps, cap)
    corr = [gain * bias * (min(m, K2) - min(m, cap)) for m in 0:(n - 1)]
    lvl0 = lvl_alt === nothing ? lvl : lvl_alt; lvl_g = lvl0 .+ corr
    A2 = A_ .+ lvl .- lvl_g; B2 = B_ .+ lvl .- lvl_g; corr = lvl_g .- lvl
    pred(nu) = (x = lvl_g .+ nu; k0 = min(j - 1, tab.n); rp_tab = tab.a .* (1 .- tab.e);
        rpp = [rp_now + (rp_tab[min(j + m, tab.n)] - rp_tab[k0]) - Hloc * nu[m + 1] for m in 0:(n - 1)]; planner_ttg_predict(tab, mu, P.ra_stop_m, kh, ra_now, rp_now, t_now, x, rpp; term=term, ra_aim=P.ra_target_m))
    term = G.term_model ? (; below=P.term_aim_below_m, steps=P.term_activate_steps, step=P.term_step_m) : nothing
    W3 = G.w_max_periods * planner_kepler_period(mu, P.ra_target_m, G.w_rp_m)
    dpg(A, B) = (r = planner_dp(A, B, L_; pen=G.pen_g); r === nothing ? nothing : r.nu)
    nu0g = any(!iszero, corr) ? dpg(A2, B2) : nu0
    nu0g === nothing && return (; nu=nu0, engaged=false, infeasible=false, ema=Float64(ema_prev), ta_raw=NaN, ta=NaN, s=NaN, steps=0, solves=0, late=false)
    ta0, steps0 = pred(nu0g); solves = 1
    isfinite(ta0) || (ta0 = t_now + 1.0e9)
    ta_raw = ta0
    if G.ema_beta > 0
        ta0 = isnan(ema_prev) ? ta0 : G.ema_beta * ema_prev + (1.0 - G.ema_beta) * ta0
    end
    ema = ta0; dlt = ta_raw - ta0
    s = G.s_frac * W3 + G.rel_slack * max(ta0 - t_now, 0.0); a_s = G.aim_frac * s
    node_ok = true
    if G.node_check                                                                                  # Revision 22 item 2 (disabled by default, NO tolerance given): a closing passage counts only if its node error is within node_tol_deg; otherwise aim at the arrival that zeroes it
        node === nothing && throw(ErrorException("guard node_check is enabled but no node state was supplied"))
        isfinite(G.node_tol_deg) || throw(ErrorException("guard node_check is enabled but node_tol_deg has no value"))
        node_ok = abs(planner_node_error_deg(node.omega_now_deg, node.rate_deg_s, t_now, ta0, node.omega_target_deg)) <= G.node_tol_deg
        node_ok || (tstar_s = planner_node_aim(node.omega_now_deg, node.rate_deg_s, t_now, node.omega_target_deg, tstar_s))
    end
    if node_ok && tstar_s - s <= ta0 <= tstar_s + s
        return (; nu=nu0, engaged=false, infeasible=false, ema, ta_raw, ta=ta0, s, steps=steps0, solves, late=false)
    end
    late = ta0 > tstar_s + s; width = maximum(B2 .- A2)
    function solve(u)
        e_ = u * G.eps_max * width
        late ? dpg(min.(A2 .+ e_, B2), B2) : dpg(A2, max.(B2 .- e_, A2))
    end
    nu_hi = solve(1.0); solves += 1
    if nu_hi === nothing
        return (; nu=nu0, engaged=true, infeasible=true, ema, ta_raw, ta=ta0, s, steps=steps0, solves, late)
    end
    th, _ = pred(nu_hi); th = isfinite(th) ? th - dlt : t_now + 1.0e9
    if (late && th > tstar_s + s) || (!late && th < tstar_s - s)
        return (; nu=nu_hi, engaged=true, infeasible=true, ema, ta_raw, ta=ta0, s, steps=steps0, solves, late)
    end
    lo_u = 0.0; hi_u = 1.0; best = nu_hi
    while solves < G.max_solves
        mid = 0.5 * (lo_u + hi_u); nu_m = solve(mid); solves += 1
        nu_m === nothing && break
        tm, _ = pred(nu_m); tm = isfinite(tm) ? tm - dlt : t_now + 1.0e9
        if abs(tm - tstar_s) <= a_s; best = nu_m; break; end
        too_late = tm > tstar_s
        if (late && too_late) || (!late && !too_late); lo_u = mid; else; hi_u = mid; best = nu_m; end
    end
    return (; nu=best, engaged=true, infeasible=false, ema, ta_raw, ta=ta0, s, steps=steps0, solves, late)
end

"""Guard setup of a flight: constants, kappa of the geometry (planner_fb_params), and the target time. t* = the model-world nominal arrival of the a priori plan from the first-pass level x0 (planner_ol_plan, then planner_arrival with
kappa; planner_rev20_replay.tstar_for), plus the geometry's aim offset t_offset_d, in seconds from the start."""
function planner_ttg_setup(meta, gid::AbstractString, tab, P, mu::Real, x0::Real; set::AbstractString="v2", tstar_d=nothing)
    G = planner_guard_constants(meta; set, target=planner_geometry_target(meta, gid)); kap = planner_fb_params(meta, gid, mu).kappa
    days = if tstar_d === nothing                                                                    # Revision 25 (default off, [planner.ic] use_plan_t_star): t* from the IC plan instead of the model-world nominal
        op = planner_ol_plan(tab, P, x0); planner_arrival(tab, mu, P.ra_stop_m, vcat(x0, op.level .+ op.nu), x0, kap).days
    else Float64(tstar_d) end
    return (; G, mu=Float64(mu), kappa=kap, tstar_nominal_d=days, tstar_s=(days + get(G.t_offset_d, String(gid), 0.0)) * 86400.0)
end

# ---- the policy ---------------------------------------------------------------------------------------------------------------------------------------

mutable struct PlannerPolicy <: SimulationModel.AbstractGuidanceModel
    f::StageFPolicy               # carries the inner CorridorGuidancePolicy, make_predictor, fit, measurements and the edge log (decisions.csv)
    mu::Float64
    omega_fn::Any                 # t_s -> rotation vector in the inertial frame
    vrel0::Float64                # first-pass relative speed from the calibration
    vr0::Float64                  # osculating V_rel at the first call; NaN until then
    mode::Symbol                  # :body | :online | :ol | :olr
    sw::NamedTuple                # (slack, ff, early, update, horizon)
    P::NamedTuple
    B::NamedTuple
    tab::Any                      # ground table or nothing (B-ODY)
    vlog::Vector{NamedTuple}      # (pass, vrel_p_mps, vp_inertial_mps) filled by the apsis callback
    h_local_m::Float64            # H_local of the calibration (live lever weight)
    q_first::Float64              # calibrated first-pass q (W/cm2)
    nu_cum::Float64               # executed cumulative ln-q offset
    nu_hist::Vector{Float64}      # nu_hist[p] = cumulative offset in effect at pass p
    apriori::Union{Nothing, Dict{Int, Float64}}   # a priori schedule from the ground table (pass => signed dv); executed as is by OL, logged against the replanned burn by the others
    last_burn_t::Float64
    last_ra::Float64              # osculating apoapsis radius at the previous decision (Revision 14 terminal aim)
    c_hist::Vector{Float64}       # observed apoapsis drop / measured q_mean per pass (m per W/cm2)
    apriori_build_s::Float64      # wall time to build the a priori schedule (one DP)
    rlog::Vector{NamedTuple}      # per-apoapsis replan log (replan.csv, Revision 11)
    u_s::Float64                  # Revision 15 remedy 4: smoothed burn-free level (NaN until the first call)
    last_dv_sign::Float64         # Revision 15 remedy 1: sign of the previous executed burn (0 = none)
    last_burn_j::Int
    fb::Any                       # Revision 19 A5 terminal-periapsis term: (; mu, kappa, weight) or nothing (switch off)
    apriori_N::Int                # arrival pass of the a priori plan (0 unless the term is on); the replanner prices the term at that pass until the Revision 20 time-to-go replaces it
    ic::Any                       # Revision 20 initial-condition plan (planner_ic_plan) or nothing ("ref"); its x0 anchors the a priori schedule
    ttg::Any                      # Revision 20 guard v2 setup (planner_ttg_setup) or nothing (switch ttg_guard off)
    ra_hist::Vector{Float64}      # osculating apoapsis radius at every decision j >= 2 (= after pass j - 1), index j - 1; filled only with the guard
    ta_ema::Float64               # guard: smoothed predicted arrival time (s), NaN until the first prediction
    ttg_N::Int                    # guard: absolute arrival pass of its last prediction (0 = none); prices the FB term when both switches are on
    node_ctx::Any                 # guard node check (default off): (; rot, omega_target_deg, J2, Re) or nothing
    node_now::Any                 # guard node check: (; omega_now_deg, rate_deg_s, omega_target_deg) of the current decision or nothing
end

"""OL-R drift update: slope (ln q per pass) of the residual r_i = realized ln q_i - (ln q_1 + path_i - path_1 + executed offset_i) over the last `update_window` flown passes (least squares); 0 until there are window + 1."""
function planner_bias(w::PlannerPolicy, j::Int)
    K = Int(w.P.update_window); meas = w.f.meas; vl = w.vlog; n = min(length(meas), length(vl), j - 1)
    (n >= K + 1 && w.tab !== nothing) || return 0.0
    lq(i) = log(0.5 * meas[i].rho_mean * vl[i].vrel_p_mps^3 / w.P.w2m)
    q1 = lq(1); p1 = w.tab.path[1]; idx = (n - K + 1):n
    r = Float64[lq(i) - (q1 + planner_tget(w.tab.path, i) - p1 + (i <= length(w.nu_hist) ? w.nu_hist[i] : w.nu_cum)) for i in idx]
    x = collect(Float64, idx); xm = sum(x) / K; rm = sum(r) / K; den = sum((x .- xm) .^ 2)
    return den > 0 ? sum((x .- xm) .* (r .- rm)) / den : 0.0
end

"""Terminal-aim bookkeeping of one decision (Revision 14): appends the observed (apoapsis drop of the last pass) / (its measured q_mean) to c_hist, records the apoapsis radius, and returns the aim (or nothing). Reads the
measurements only; the c_hist / last_ra it updates feed nothing but the aim."""
function planner_terminal_track!(w::PlannerPolicy, a_now::Float64, e_now::Float64)
    meas = w.f.meas; P = w.P; ra_now = a_now * (1 + e_now); nm = min(length(meas), length(w.vlog))
    if isfinite(w.last_ra) && nm >= 1
        drop = w.last_ra - ra_now; qprev = 0.5 * meas[nm].rho_mean * w.vlog[nm].vrel_p_mps^3 / P.w2m
        (drop > 0 && qprev > 0) && push!(w.c_hist, drop / qprev)
    end
    w.last_ra = ra_now
    return planner_terminal_aim(ra_now - P.ra_target_m, w.c_hist, P)
end

"""Predicted burn-free ln q of the next pass (z0) and the density offset yhat (the P2_10 EWMA of the measured mean-field offset): the quantity every measurement-driven decision starts from."""
function planner_predict_z0(w::PlannerPolicy, predictor, vrel_now::Float64)
    f = w.f; P = w.P; meas = f.meas
    ym = Float64[log(m.rho_mean / corridor_density(f.fit, m.alt_m)) for m in meas]
    yhat = isempty(ym) ? 0.0 : stagef_predict_next("P2_10", ym, Float64[m.lon_rad for m in meas])
    z0 = log(0.5 * corridor_density(f.fit, predictor(0.0)) * exp(yhat) * vrel_now^3 / P.w2m)
    return z0, yhat
end

"""Decision before pass j. Returns (dv signed m/s, z0, w0, yhat, bias, nu0, n_plan, reason). Mode :ol reads the ground schedule (no predictor call, no measurement); with the terminal switch (arm OLT, Revision 20b) the terminal
aim is checked first, from the measurements alone, and only once it is active are the predictor and the aim used: before activation the decision, and every logged column, equal OL's."""
function planner_decide(w::PlannerPolicy, j::Int, t::Float64, predictor, vrel_now::Float64, a_now::Float64, e_now::Float64)
    f = w.f; P = w.P
    w0 = planner_weight(w.h_local_m, P.burn_response_ratio, w.mu, a_now, e_now)
    if w.sw.ttg_guard                                                  # guard bookkeeping: the apoapsis radius after pass j - 1 (index j - 1)
        while length(w.ra_hist) < j - 2; push!(w.ra_hist, NaN); end
        length(w.ra_hist) == j - 2 ? push!(w.ra_hist, a_now * (1 + e_now)) : (w.ra_hist[j - 1] = a_now * (1 + e_now))
    end
    if w.mode === :ol
        if w.sw.terminal
            ta = planner_terminal_track!(w, a_now, e_now)
            if ta !== nothing
                z0, yhat = planner_predict_z0(w, predictor, vrel_now)
                nu0 = ta.lnq - z0; dv = -nu0 * w0; abs(dv) < P.dp_min_burn_mps && (dv = 0.0)
                return (; dv, z0, w0, yhat, bias=0.0, nu0, n_plan=ta.N, reason=dv == 0.0 ? "terminal_hold" : "terminal_aim", dp_s=0.0)
            end
        end
        dv = get(w.apriori, j, 0.0)
        return (; dv, z0=NaN, w0, yhat=NaN, bias=0.0, nu0=NaN, n_plan=length(w.apriori), reason=dv == 0.0 ? "ol_hold" : "ol_scheduled", dp_s=0.0)
    end
    meas = f.meas
    z0, yhat = planner_predict_z0(w, predictor, vrel_now)
    if w.sw.terminal                                                  # Revision 14: terminal aim replaces the tube DP once the target is within reach
        ta = planner_terminal_track!(w, a_now, e_now)
        if ta !== nothing
            nu0 = ta.lnq - z0; dv = -nu0 * w0; abs(dv) < P.dp_min_burn_mps && (dv = 0.0)
            return (; dv, z0, w0, yhat, bias=0.0, nu0, n_plan=ta.N, reason=dv == 0.0 ? "terminal_hold" : "terminal_aim", dp_s=0.0)
        end
    end
    if w.mode === :body
        dv = planner_body_dv(w.B, P, j, t, w.last_burn_t, z0, w0)
        return (; dv, z0, w0, yhat, bias=0.0, nu0=NaN, n_plan=0, reason=dv == 0.0 ? "body_hold" : "body_recenter", dp_s=0.0)
    end
    hi0 = log(P.q_cap_mean_w_cm2); lo0f = hi0 - (w.sw.slack ? w.tab.W : P.l0_width_ln); outside = z0 > hi0 || z0 < lo0f      # predicted level outside the FULL tube (forced)
    if w.sw.smooth_beta > 0                                                                                                      # Revision 15 remedy 4: smooth the burn-free level
        u = z0 - w.nu_cum; w.u_s = isnan(w.u_s) ? u : w.sw.smooth_beta * w.u_s + (1 - w.sw.smooth_beta) * u
        outside || (z0 = w.u_s + w.nu_cum)
    end
    bias = (w.mode === :olr && w.sw.update) ? planner_bias(w, j) : 0.0
    fbN = (w.sw.ttg_guard && w.ttg_N > 0) ? w.ttg_N : w.apriori_N                                                                  # the guard's arrival pass (previous decision) prices the term when both switches are on
    fbon = (w.sw.fb_term || w.P.apo_one_sided) && w.fb !== nothing && w.sw.ff && fbN > 0                                                        # Revision 19 A5 (default off): terminal-periapsis cost at the a priori plan's arrival pass
    inp = fbon ? planner_inputs(w.sw, w.tab, P, j, z0, w0; bias=bias, lower_margin=w.sw.lower_margin, fb=w.fb, fb_N=fbN, mu=w.mu) : planner_inputs(w.sw, w.tab, P, j, z0, w0; bias=bias, lower_margin=w.sw.lower_margin)
    t_dp = time_ns(); pl = fbon ? planner_plan(inp.level, inp.lo, inp.hi, inp.w; eps=P.dp_tie_eps, pen=P.dp_burn_penalty_mps, lin=inp.lin) : planner_plan(inp.level, inp.lo, inp.hi, inp.w; eps=P.dp_tie_eps, pen=P.dp_burn_penalty_mps); dp_s = (time_ns() - t_dp) * 1.0e-9
    pl === nothing && return (; dv=0.0, z0, w0, yhat, bias, nu0=NaN, n_plan=length(inp.level), reason="dp_infeasible", dp_s)
    dv = -pl.nu[1] * w0
    if w.sw.dev_thr > 0 && w.mode === :olr && !outside                                                                           # Revision 15 remedy 3: follow the a priori burn unless the replan differs by more than the threshold
        ap = get(w.apriori, j, 0.0); abs(dv - ap) <= w.sw.dev_thr && (dv = ap)
    end
    if w.sw.rev_deadband > 0 && dv != 0.0 && !outside && w.last_dv_sign != 0.0 && sign(dv) != w.last_dv_sign && j - w.last_burn_j <= w.sw.rev_window && abs(pl.nu[1]) < w.sw.rev_deadband      # remedy 1
        dv = 0.0
    end
    abs(dv) < max(P.dp_min_burn_mps, w.sw.min_burn_arm) && (dv = 0.0)                                                            # remedy 2 (minimum burn size; the penalty is the arm's burn_penalty_mps)
    if w.sw.ttg_guard && w.ttg !== nothing && j >= w.ttg.G.start_pass                                                           # Revision 20 guard v2 (default off)
        g = planner_ttg_policy_step!(w, j, inp, pl, w0, bias, dv, a_now, e_now)
        return (; dv=g.dv, z0, w0, yhat, bias, nu0=pl.nu[1], n_plan=length(inp.level), reason=g.dv == 0.0 ? "plan_hold" : (g.engaged ? "guard_burn" : "plan_burn"), dp_s,
            ttg_engaged=g.engaged, ttg_infeasible=g.infeasible, ttg_ta_s=g.ta, ttg_s_s=g.s, ttg_tt_s=w.ttg.tstar_s, ttg_solves=g.solves, ttg_dv_ref=dv)
    end
    return (; dv, z0, w0, yhat, bias, nu0=pl.nu[1], n_plan=length(inp.level), reason=dv == 0.0 ? "plan_hold" : "plan_burn", dp_s)
end

"""Guard step of one planner decision (after the planner's own burn dv_ref is known): gathers the flown history (realized ln q, apoapsis radii, executed offsets), runs planner_ttg_decide, applies the guard's burn when the
tube actuator engaged (minimum burn as the planner's, no a priori snap) and keeps the planner's own burn unless the guard differs by more than dev_thr_g. Updates ta_ema and ttg_N. Returns (; dv, engaged, infeasible, ta, s, solves)."""
function planner_ttg_policy_step!(w::PlannerPolicy, j::Int, inp, pl, w0::Float64, bias::Float64, dv_ref::Float64, a_now::Float64, e_now::Float64)
    T = w.ttg; G = T.G; P = w.P; tab = w.tab; meas = w.f.meas; vl = w.vlog; nf = min(length(meas), length(vl), j - 1)
    xs = Float64[log(0.5 * meas[i].rho_mean * vl[i].vrel_p_mps^3 / P.w2m) for i in 1:nf]; ras = w.ra_hist[1:min(nf, length(w.ra_hist))]
    nuh = Float64[i <= length(w.nu_hist) ? w.nu_hist[i] : w.nu_cum for i in 1:nf]; u0 = planner_ic_u0(w.ic); u = u0 + w.nu_cum
    lvl = inp.level; A_ = inp.lo .- lvl; B_ = inp.hi .- lvl; L_ = planner_tie_weights(inp.w, P.dp_tie_eps)
    ra_now = a_now * (1 + e_now); rp_now = a_now * (1 - e_now)
    t_now = nf >= 1 ? Float64(meas[nf].t_s) : 0.0
    kap_meas = length(ras) == nf ? planner_ttg_kappa_hat(xs, ras, tab, G.k_window) : nothing
    r = planner_ttg_decide(G, tab, P, T.mu, T.kappa, j, lvl, A_, B_, L_, pl.nu, ra_now, rp_now, t_now, T.tstar_s, kap_meas, bias, w.h_local_m,
        length(ras) == nf ? planner_ttg_lvl_alt(G, tab, j, lvl, xs, nuh, u0, u) : nothing, w.ta_ema; node=G.node_check ? w.node_now : nothing)
    w.ta_ema = r.ema; r.steps > 0 && isfinite(r.ta_raw) && (w.ttg_N = j + r.steps - 1)
    dv = dv_ref
    if r.engaged
        dg = -r.nu[1] * w0; abs(dg) < max(P.dp_min_burn_mps, w.sw.min_burn_arm) && (dg = 0.0)
        dv = abs(dg - dv_ref) <= G.dev_thr_g ? dv_ref : dg
    end
    return (; dv, engaged=r.engaged, infeasible=r.infeasible, ta=r.ta, s=r.s, solves=r.solves)
end

function SimulationModel.GuidanceHooks.calcGuidanceEffect!(w::PlannerPolicy, u, p, t::Float64, sat_idx::Int64)
    f = w.f; policy = f.inner; commands = p.shared_buffers.maneuver_commands
    (sat_idx < 1 || sat_idx > length(commands)) && return nothing
    if !policy.armed
        commands[sat_idx] = SimulationModel.CommandTypes.PropulsiveManeuverCommand(valid=true, source_orbit=policy.periapsis_count + 1)
        return nothing
    end
    r = SVector{3, Float64}(u.sc[sat_idx].pos); v = SVector{3, Float64}(u.sc[sat_idx].vel); mass = Float64(u.sc[sat_idx].mass)
    vehicle = policy.config.values["vehicle"]; exhaust_velocity = Float64(vehicle["specific_impulse_s"]) * 9.80665
    consumed_at_limit = Float64(vehicle["thrust_n"]) * policy.max_burn_duration_s / exhaust_velocity
    duration_limited_dv = consumed_at_limit < mass ? exhaust_velocity * log(mass / (mass - consumed_at_limit)) : Inf
    wet_mass = Float64(vehicle["dry_mass_kg"]) + Float64(vehicle["propellant_mass_kg"])
    delivered_so_far = mass < wet_mass ? exhaust_velocity * log(wet_mass / mass) : 0.0
    policy.remaining_budget_mps = max(0.0, Float64(policy.config.values["comparison"]["total_delta_v_budget_mps"]) - delivered_so_far)
    max_allowed_dv = min(policy.max_dv_mps, policy.remaining_budget_mps, duration_limited_dv)
    j = policy.periapsis_count + 1
    vr = periapsis_vrel(r, v, w.mu, w.omega_fn(t)); isnan(w.vr0) && (w.vr0 = vr); vrel_now = w.vrel0 * vr / w.vr0      # U15, as StageGPolicy
    pv = periapsis_vectors(r, v, w.mu); a_now = pv.rp / (1 - pv.e)
    predictor = f.make_predictor(p, r, v, mass, Float64(t))
    if w.mode !== :body && w.apriori === nothing                     # the a priori schedule: one DP on the ground table and the calibrated first-pass q, no measurement
        tb = time_ns(); op = planner_ol_plan(w.tab, w.P, planner_anchor_ln(w.q_first, w.ic); fb=(w.sw.fb_term || w.P.apo_one_sided) ? w.fb : nothing); w.apriori = op.sched; w.apriori_N = op.N; w.apriori_build_s = (time_ns() - tb) * 1.0e-9      # Revision 20: the planned x0 anchors the schedule when an IC scenario is on
    end
    if w.ttg !== nothing && w.ttg.G.node_check                                  # guard node check (default off): node of the orbit in the body-equator inertial frame and its secular J2 rate
        nc = w.node_ctx; nc === nothing && error("guard node_check is on but the target node is not set ([planner.targets.*] node_deg)")
        hb = cross(transpose(nc.rot) * r, transpose(nc.rot) * v); inc_rad = acos(clamp(hb[3] / norm(hb), -1.0, 1.0))
        w.node_now = (; omega_now_deg=rad2deg(atan(hb[1], -hb[2])), rate_deg_s=rad2deg(planner_j2_node_rate(w.mu, nc.J2, nc.Re, a_now, pv.e, inc_rad)), omega_target_deg=nc.omega_target_deg)
    end
    d = planner_decide(w, j, Float64(t), predictor, vrel_now, a_now, pv.e)
    dv = clamp(d.dv, -max_allowed_dv, max_allowed_dv)
    while length(w.nu_hist) < j - 1; push!(w.nu_hist, w.nu_hist[end]); end
    dv != 0.0 && (w.nu_cum += -dv / d.w0; w.last_burn_t = Float64(t); w.last_dv_sign = sign(dv); w.last_burn_j = j)
    length(w.nu_hist) == j - 1 ? push!(w.nu_hist, w.nu_cum) : (w.nu_hist[j] = w.nu_cum)
    push!(f.edge_log, (; request_time_s=Float64(t), target_pass=j, n_meas=length(f.meas), yhat=d.yhat, lower_m=NaN, upper_m=NaN, collapsed=false, delta_v_mps=dv, reason=d.reason,
        z0=d.z0, w0=d.w0, bias=d.bias, nu0=d.nu0, n_plan=d.n_plan, nu_cum=w.nu_cum))
    nm = min(length(f.meas), length(w.vlog)); lq_meas = nm >= 1 ? log(0.5 * f.meas[nm].rho_mean * w.vlog[nm].vrel_p_mps^3 / w.P.w2m) : NaN
    push!(w.rlog, (; pass=j, t_s=Float64(t), apriori_dv_mps=(w.apriori === nothing ? NaN : get(w.apriori, j, 0.0)), planned_dv_mps=d.dv, executed_dv_mps=dv, reason=d.reason,
        a_m=a_now, e=pv.e, r_m=norm(r), vrel_mps=vrel_now, lnq_last_measured=lq_meas, z0_pred_lnq=d.z0, density_offset_yhat=d.yhat, drift_bias_per_pass=d.bias, nu_cum=w.nu_cum, w0=d.w0,
        n_plan=d.n_plan, dp_solve_s=d.dp_s, apriori_build_s=(j == 2 ? w.apriori_build_s : 0.0)))
    if w.sw.ttg_guard                                                                                                             # guard columns, every row of a guard arm (Revision 20)
        w.rlog[end] = merge(w.rlog[end], (; ttg_engaged=get(d, :ttg_engaged, false), ttg_infeasible=get(d, :ttg_infeasible, false), ttg_ta_s=get(d, :ttg_ta_s, NaN), ttg_s_s=get(d, :ttg_s_s, NaN),
            ttg_tt_s=(w.ttg === nothing ? NaN : w.ttg.tstar_s), ttg_solves=get(d, :ttg_solves, 0), ttg_dv_ref=get(d, :ttg_dv_ref, NaN), ttg_arrival_pass=w.ttg_N))
    end
    push!(policy.decisions, (request_time_s=Float64(t), source_orbit=j, target_pass=j, policy="planner_" * String(w.mode), delta_v_mps=dv,
        predicted_periapsis_altitude_m=NaN, target_periapsis_altitude_m=NaN, trigger_reason=d.reason, feasible=true, predictor_evaluations=1))
    # the pre-apoapsis callback (examples/AGORA_Earth_Periapsis_Corridor.jl) reads `last_selection` to validate the finite-thrust plan, record it in burn_plans and stop on a rejected plan
    policy.last_selection = (; delta_v_mps=dv, feasible=true, reason=Symbol(d.reason), error="", target_periapsis_altitude_m=NaN, predicted_periapsis_altitude_m=NaN, evaluations=1)
    commands[sat_idx] = SimulationModel.CommandTypes.PropulsiveManeuverCommand(valid=true, delta_v_mps=abs(dv), direction_rad=dv >= 0.0 ? 0.0 : π, source_orbit=j)
    policy.armed = false
    return nothing
end

# ---- terminal wrap for the Stage O keeper (Revision 16): N2-I-T -------------------------------------------------------------------------------------------

"""Wrapper around the UNCHANGED Stage O keeper (StageOPolicy over StageN2Policy). With `terminal` false it only delegates (reproduction gate: the flight must equal the Stage O flight). With `terminal` true, once the target
is within reach (planner_terminal_aim, the SAME function and constants as the planner arms) it takes over the burn decision of that orbit with the terminal aim; before that and whenever the aim is inactive it delegates."""
mutable struct TerminalWrap <: SimulationModel.AbstractGuidanceModel
    inner::Any                    # StageOPolicy
    f::StageFPolicy
    mu::Float64
    omega_fn::Any
    vrel0::Float64
    vr0::Float64
    P::NamedTuple
    vlog::Vector{NamedTuple}
    h_local_m::Float64
    terminal::Bool
    last_ra::Float64
    c_hist::Vector{Float64}
    nu_cum::Float64
    apriori_build_s::Float64
    rlog::Vector{NamedTuple}
    n2::NamedTuple                # (t_mult, w_m, aim) of the keeper, only to fill the same decision-log columns
end

function SimulationModel.GuidanceHooks.calcGuidanceEffect!(w::TerminalWrap, u, p, t::Float64, sat_idx::Int64)
    f = w.f; policy = f.inner; commands = p.shared_buffers.maneuver_commands
    (w.terminal && sat_idx >= 1 && sat_idx <= length(commands) && policy.armed) || return SimulationModel.GuidanceHooks.calcGuidanceEffect!(w.inner, u, p, t, sat_idx)
    r = SVector{3, Float64}(u.sc[sat_idx].pos); v = SVector{3, Float64}(u.sc[sat_idx].vel); mass = Float64(u.sc[sat_idx].mass)
    j = policy.periapsis_count + 1; pv = periapsis_vectors(r, v, w.mu); a_now = pv.rp / (1 - pv.e); ra_now = a_now * (1 + pv.e)
    meas = f.meas; nm = min(length(meas), length(w.vlog))
    if isfinite(w.last_ra) && nm >= 1
        drop = w.last_ra - ra_now; qprev = 0.5 * meas[nm].rho_mean * w.vlog[nm].vrel_p_mps^3 / w.P.w2m
        (drop > 0 && qprev > 0) && push!(w.c_hist, drop / qprev)
    end
    w.last_ra = ra_now
    ta = planner_terminal_aim(ra_now - w.P.ra_target_m, w.c_hist, w.P)
    ta === nothing && return SimulationModel.GuidanceHooks.calcGuidanceEffect!(w.inner, u, p, t, sat_idx)
    vehicle = policy.config.values["vehicle"]; exhaust_velocity = Float64(vehicle["specific_impulse_s"]) * 9.80665
    consumed_at_limit = Float64(vehicle["thrust_n"]) * policy.max_burn_duration_s / exhaust_velocity
    duration_limited_dv = consumed_at_limit < mass ? exhaust_velocity * log(mass / (mass - consumed_at_limit)) : Inf
    wet_mass = Float64(vehicle["dry_mass_kg"]) + Float64(vehicle["propellant_mass_kg"]); delivered_so_far = mass < wet_mass ? exhaust_velocity * log(wet_mass / mass) : 0.0
    policy.remaining_budget_mps = max(0.0, Float64(policy.config.values["comparison"]["total_delta_v_budget_mps"]) - delivered_so_far)
    max_allowed_dv = min(policy.max_dv_mps, policy.remaining_budget_mps, duration_limited_dv)
    vr = periapsis_vrel(r, v, w.mu, w.omega_fn(t)); isnan(w.vr0) && (w.vr0 = vr); vrel_now = w.vrel0 * vr / w.vr0
    predictor = f.make_predictor(p, r, v, mass, Float64(t))
    ym = Float64[log(m.rho_mean / corridor_density(f.fit, m.alt_m)) for m in meas]; yhat = isempty(ym) ? 0.0 : stagef_predict_next("P2_10", ym, Float64[m.lon_rad for m in meas])
    z0 = log(0.5 * corridor_density(f.fit, predictor(0.0)) * exp(yhat) * vrel_now^3 / w.P.w2m); w0 = planner_weight(w.h_local_m, w.P.burn_response_ratio, w.mu, a_now, pv.e)
    dv = -(ta.lnq - z0) * w0; abs(dv) < w.P.dp_min_burn_mps && (dv = 0.0); dv = clamp(dv, -max_allowed_dv, max_allowed_dv)
    dv != 0.0 && (w.nu_cum += -dv / w0)
    push!(f.edge_log, (; request_time_s=Float64(t), target_pass=j, n_meas=length(meas), yhat, lower_m=NaN, upper_m=NaN, collapsed=false, delta_v_mps=dv, reason="terminal_aim", m=NaN, s=NaN, c=NaN,
        t_mult=w.n2.t_mult, w_m=w.n2.w_m, aim=w.n2.aim, target_m=NaN))
    push!(w.rlog, (; pass=j, t_s=Float64(t), reason="terminal_aim", executed_dv_mps=dv, z0_pred_lnq=z0, dp_solve_s=0.0, apriori_build_s=0.0))
    push!(policy.decisions, (request_time_s=Float64(t), source_orbit=j, target_pass=j, policy="terminal_aim", delta_v_mps=dv, predicted_periapsis_altitude_m=NaN, target_periapsis_altitude_m=NaN,
        trigger_reason="terminal_aim", feasible=true, predictor_evaluations=1))
    policy.last_selection = (; delta_v_mps=dv, feasible=true, reason=:terminal_aim, error="", target_periapsis_altitude_m=NaN, predicted_periapsis_altitude_m=NaN, evaluations=1)
    commands[sat_idx] = SimulationModel.CommandTypes.PropulsiveManeuverCommand(valid=true, delta_v_mps=abs(dv), direction_rad=dv >= 0.0 ? 0.0 : π, source_orbit=j)
    policy.armed = false
    return nothing
end

# ---- geometry list (ODY injected), constants ------------------------------------------------------------------------------------------------------

"""stageg3_meta plus the Odyssey reference geometry (id ODY, from [planner.geometry.ODY]) appended to meta["stageg3"]["geometry"], so stageg3_load / stageg3_cal_task / stageg3_run_task machinery accepts gid "ODY".
Idempotent (the dict is cached by stageg3_meta per config path). Call before any stageg3_load in a process."""
function planner_meta(path::AbstractString)
    meta = stageg3_meta(path); G = meta["stageg3"]["geometry"]
    for (key, o) in meta["planner"]["geometry"]            # ODY (Odyssey reference) and, Revision 16, the generic orbits GA to GD
        gid = String(o["id"])
        any(g -> String(g["id"]) == gid, G) || push!(G, Dict{String, Any}("id" => gid, "role" => String(o["role"]), "inc_deg" => Float64(o["inc_deg"]), "omega_deg" => Float64(o["omega_deg"]), "raan_deg" => Float64(o["raan_deg"]),
            "epoch_off_d" => Float64(o["epoch_off_d"]), "epoch_utc" => String(o["epoch_utc"]), "field_id" => gid))
    end
    return meta
end

"""Geometry id of a cell token: ODY -> ODY; G1 / G2 / Gworst -> the id of that role in the gate geometry table."""
function planner_gid(meta, token::AbstractString)
    token == "ODY" && return "ODY"
    i = findfirst(g -> String(g["role"]) == token, meta["stageg3"]["geometry"]); i === nothing && throw(ArgumentError("unknown geometry token $token"))
    return String(meta["stageg3"]["geometry"][i]["id"])
end

"""Target table of [planner.target] (the default target) merged with [planner.targets.<name>] (Revision 25: the apoapsis target is a parameter). name = nothing / "" / "default" = the default table. Keys (all optional; the default is
the current value): ra_target_m, tail_passes, horizon_days, term_aim_below_m, term_step_m, term_activate_steps, term_min_q_factor, term_c_window, guard_w_rp_m, event, apo_one_sided, apo_tol_m, apo_aim_above_frac,
apo_shortfall_slope_mps_per_km, node_deg."""
function planner_target_table(meta, name::Union{Nothing, AbstractString}=nothing)
    PL = meta["planner"]; base = Dict{String, Any}(get(PL, "target", Dict{String, Any}()))
    (name === nothing || isempty(name) || name == "default") && return base
    haskey(get(PL, "targets", Dict{String, Any}()), String(name)) || error("unknown target $name ([planner.targets.$name] missing)")
    return merge(base, Dict{String, Any}(PL["targets"][String(name)]))
end
"""Target name of a geometry id ([planner.geometry.<key>] field `target`), nothing = the default target."""
function planner_geometry_target(meta, gid::AbstractString)
    for (_, o) in get(meta["planner"], "geometry", Dict{String, Any}()); String(o["id"]) == String(gid) && return haskey(o, "target") ? String(o["target"]) : nothing; end
    return nothing
end
"""Horizon (days) of a geometry: its target's `horizon_days` when set, else [planner] horizon_days."""
planner_horizon_days(meta, gid::AbstractString) = Float64(get(planner_target_table(meta, planner_geometry_target(meta, gid)), "horizon_days", meta["planner"]["horizon_days"]))

"""Planner constants as NamedTuples (P: arithmetic constants; B: B-ODY; sw: switches of arm `key`). Revision 25: `target` selects [planner.targets.<name>] (nothing = the default target)."""
function planner_constants(meta, w2m::Real, h_local_m::Real, arm::AbstractString=""; target::Union{Nothing, AbstractString}=nothing)
    PL = meta["planner"]; D = PL["dp"]; U = PL["update"]; f = k -> Float64(PL[k])
    P = (; q_cap_mean_w_cm2=f("q_cap_mean_w_cm2"), l0_width_ln=f("l0_width_ln"), fast_reserve_ln=f("fast_reserve_ln"), tail_passes=Int(PL["tail_passes"]), burn_response_ratio=f("burn_response_ratio"),
        q_never_exceed_qs_w_cm2=f("q_never_exceed_qs_w_cm2"), dp_tie_eps=Float64(D["tie_eps"]), dp_min_burn_mps=Float64(D["min_burn_mps"]), dp_burn_penalty_mps=Float64(get(get(PL["arms"], arm, Dict{String, Any}()), "burn_penalty_mps", D["burn_penalty_mps"])), update_gain=Float64(U["gain"]), update_cap_passes=Int(U["cap_passes"]),
        update_window=Int(U["window"]), w2m=Float64(w2m), h_local_m=Float64(h_local_m),
        ra_target_m=Float64(meta["stageo"]["ra_target_m"]), term_aim_below_m=Float64(get(get(PL, "terminal", Dict()), "term_aim_below_m", 12500.0)), term_step_m=Float64(get(get(PL, "terminal", Dict()), "term_step_m", 40000.0)),
        term_activate_steps=Int(get(get(PL, "terminal", Dict()), "term_activate_steps", 3)), term_min_q_factor=Float64(get(get(PL, "terminal", Dict()), "term_min_q_factor", 0.25)), term_c_window=Int(get(get(PL, "terminal", Dict()), "term_c_window", 6)))
    T = planner_target_table(meta, target); String(get(T, "event", "apoapsis")) == "apoapsis" || throw(ErrorException("Not implemented: decision event $(T["event"]) (Revision 22 item 5: apoapsis is the only event)"))
    one_sided = Bool(get(T, "apo_one_sided", false)); tol = Float64(get(T, "apo_tol_m", 25000.0)); frac = Float64(get(T, "apo_aim_above_frac", 0.5))
    P = merge(P, (; ra_target_m=Float64(get(T, "ra_target_m", P.ra_target_m)), tail_passes=Int(get(T, "tail_passes", P.tail_passes)), term_step_m=Float64(get(T, "term_step_m", P.term_step_m)),
        term_activate_steps=Int(get(T, "term_activate_steps", P.term_activate_steps)), term_min_q_factor=Float64(get(T, "term_min_q_factor", P.term_min_q_factor)), term_c_window=Int(get(T, "term_c_window", P.term_c_window)),
        term_aim_below_m=one_sided ? -tol * frac : Float64(get(T, "term_aim_below_m", P.term_aim_below_m)),                     # one-sided tolerance (Revision 22 item 3): the aim sits apo_aim_above_frac of the band ABOVE the target (a negative "below")
        ra_stop_m=Float64(get(T, "ra_target_m", P.ra_target_m)) + (one_sided ? tol : 0.0),                                          # closing pass = first apoapsis at or below target + tol (default: the target itself)
        apo_one_sided=one_sided, apo_tol_m=tol, apo_aim_above_frac=frac, apo_shortfall_slope_mps_per_km=Float64(get(T, "apo_shortfall_slope_mps_per_km", 0.056))))
    b = PL["body"]
    B = (; q_low_qs_w_cm2=Float64(b["q_low_qs_w_cm2"]), q_high_qs_w_cm2=Float64(b["q_high_qs_w_cm2"]), q_target_mp1_qs_w_cm2=Float64(b["q_target_mp1_qs_w_cm2"]), q_target_mp2_qs_w_cm2=Float64(b["q_target_mp2_qs_w_cm2"]),
        qs_over_qmean=Float64(b["qs_over_qmean"]), cadence_s=Float64(b["cadence_s"]), mp1_last_pass=Int(b["mp1_last_pass"]))
    return P, B
end
planner_switches(meta, key::AbstractString) = (a = meta["planner"]["arms"][key]; (; slack=Bool(get(a, "slack", false)), ff=Bool(get(a, "ff", false)), early=Bool(get(a, "early", false)),
    update=Bool(get(a, "update", false)), horizon=Int(get(a, "horizon", 0)), terminal=Bool(get(a, "terminal", false)),
    rev_window=Int(get(a, "rev_window_orbits", 0)), rev_deadband=Float64(get(a, "rev_deadband_ln", 0.0)), min_burn_arm=Float64(get(a, "min_burn_arm_mps", 0.0)), dev_thr=Float64(get(a, "dev_thr_mps", 0.0)),
    smooth_beta=Float64(get(a, "smooth_beta", 0.0)), lower_margin=Float64(get(a, "lower_margin_ln", 0.0)),
    fb_term=Bool(get(a, "fb_term", false)), ic_scenario=String(get(a, "ic_scenario", "ref")), ttg_guard=Bool(get(a, "ttg_guard", false))))                                  # Revision 19 A5 / Revision 20 switches, default off / "ref"
planner_arm_code(key::AbstractString) = "PL_" * key

"""Stage O N2 keeper (kept unchanged: Kalman constants, edge multiplier, aim, band cap from the arm's config fields, SOURCED stagen_flight kalman_params.csv, stageo_data calibrated_arms.csv and stageo_config [stageo.n2]) inside the
terminal wrapper. Arm fields: n2_kr, n2_kq, n2_km0, n2_kv0, n2_t_mult, n2_aim, n2_w_m, n2_track_vrel; the wrapper's `terminal` follows the arm's terminal switch."""
function planner_n2it_policy(meta, A, fpol, ms, omega_fn, G, P, vlog, sw)
    n2 = (; t_mult=Float64(A["n2_t_mult"]), w_m=Float64(A["n2_w_m"]), aim=Float64(A["n2_aim"]))
    pol = StageN2Policy(fpol, Float64(A["n2_kr"]), Float64(A["n2_kq"]), Float64(A["n2_km0"]), Float64(A["n2_kv0"]), n2.t_mult, n2.w_m, n2.aim)
    so = StageOPolicy(pol, fpol, Float64(ms.harm.gm_m3s2), Bool(A["n2_track_vrel"]), NaN, Float64(G.vrel_first_mps))
    return TerminalWrap(so, fpol, Float64(ms.harm.gm_m3s2), omega_fn, Float64(G.vrel_first_mps), NaN, P, vlog, Float64(G.hloc_at_hg_m), sw.terminal, NaN, Float64[], 0.0, 0.0, NamedTuple[], n2)
end

"""kappa of a geometry id from the calibration file [planner.terminal_periapsis] kappa_file (path relative to the task directory; entry kappa_med, used only when `accepted`), else `default`. Revision 21: G0975 (MEASURED by planner_rev18_calibrate from the gen1 flights)."""
function planner_kappa_file(t, gid::AbstractString, default::Real)
    haskey(t, "kappa_file") || return Float64(default)
    fp = joinpath(stageg3_task_dir(), String(t["kappa_file"])); isfile(fp) || return Float64(default)
    J = JSON.parsefile(fp); (haskey(J, gid) && Bool(get(J[gid], "accepted", false))) || return Float64(default)
    return Float64(J[gid]["kappa_med"])
end

"""Terminal-periapsis parameters of geometry id `gid` ([planner.terminal_periapsis]): (; mu, kappa, weight). kappa is MEASURED per geometry (Revision 19 B), kappa_default ASSUMED for the others."""
function planner_fb_params(meta, gid::AbstractString, mu::Real)
    t = meta["planner"]["terminal_periapsis"]; kt = get(t, "kappa", Dict{String, Any}())
    kappa = haskey(kt, String(gid)) ? Float64(kt[String(gid)]) : planner_kappa_file(t, String(gid), Float64(t["kappa_default"]))
    return (; mu=Float64(mu), kappa, weight=Float64(t["fb_weight"]))
end

# ---- Revision 20 initial conditions (WP1): ic_scenario ref | S-fixed | S-free -------------------------------------------------------------------------------

"""Initial-condition plan of a pinned planner_rev20_ic.select_ic dump (path and sha256 in [planner.ic]; a mismatch or an empty hash is an ERROR unless `file` is given). The file is a JSON object with a list `plans` (as
rev20_data/ic/headline.json) or the list itself; the entry is selected by geom == token, scenario, and the [planner.ic] selectors fb_term, horizon, delta_ln, h_des_m (t_star_d must be null). "ref" returns nothing (the unchanged flight).
Returns (; scenario, x0, nu_pre, w_pre, dv_pre_signed (m/s, + raises the periapsis: -nu_pre w_pre), dv_pre_charged (m/s; S-fixed only), h_p0_m)."""
function planner_ic_plan(meta, token::AbstractString, scenario::AbstractString; file::Union{Nothing, AbstractString}=nothing)
    scenario == "ref" && return nothing
    scenario in ("S-fixed", "S-free") || error("ic_scenario must be ref, S-fixed or S-free, got $scenario")
    pt = meta["planner"]["ic"]; fp = file === nothing ? joinpath(stageg3_task_dir(), String(pt["path"])) : String(file); raw = read(fp)
    if file === nothing
        want = String(pt["sha256"]); isempty(want) && error("planner_config.toml [planner.ic] sha256 is empty: pin the select_ic dump before flying an IC scenario")
        h = bytes2hex(sha256(raw)); h == want || error("IC plan file sha256 mismatch: $fp is $h, config has $want")
    end
    J = JSON.parse(String(raw)); plans = collect(J isa AbstractDict ? J["plans"] : J)
    if file === nothing                                                                              # Revision 25: further pinned IC files ([planner.ic_extra.<name>] path, sha256), plans appended
        for (nm, ex) in get(meta["planner"], "ic_extra", Dict{String, Any}())
            fe = joinpath(stageg3_task_dir(), String(ex["path"])); rawe = read(fe); wante = String(get(ex, "sha256", "")); isempty(wante) && error("[planner.ic_extra.$nm] sha256 is empty")
            bytes2hex(sha256(rawe)) == wante || error("IC plan file sha256 mismatch: $fe (config has $wante)")
            Je = JSON.parse(String(rawe)); append!(plans, Je isa AbstractDict ? Je["plans"] : Je)
        end
    end
    tsel = get(pt, "t_star_d", nothing)                                                              # Revision 25 (default null): select the plan computed for the arrival target t* (days)
    fbt = Bool(get(pt, "fb_term", false)); hz = String(get(pt, "horizon", "arrival")); dl = Float64(get(pt, "delta_ln", 0.0)); hd = Float64(get(pt, "h_des_m", 250000.0))
    cand = [p for p in plans if String(p["geom"]) == token && String(p["scenario"]) == scenario && Bool(p["fb_term"]) == fbt && String(p["horizon"]) == hz && isapprox(Float64(p["delta_ln"]), dl; atol=1.0e-12) &&
        isapprox(Float64(p["h_des_m"]), hd; atol=1.0e-6) && (tsel === nothing ? p["t_star_d"] === nothing : (p["t_star_d"] !== nothing && isapprox(Float64(p["t_star_d"]), Float64(tsel); atol=1.0e-9)))]
    length(cand) == 1 || error("IC plan file has $(length(cand)) entries for geom=$token scenario=$scenario fb_term=$fbt horizon=$hz delta_ln=$dl h_des_m=$hd (need exactly one)")
    p = cand[1]; Bool(p["feasible"]) || error("IC plan for $token $scenario is infeasible")
    nu_pre = Float64(p["nu_pre"]); w_pre = Float64(p["w_pre"])
    return (; scenario, x0=Float64(p["x0"]), nu_pre, w_pre, dv_pre_signed=-nu_pre * w_pre, dv_pre_charged=scenario == "S-fixed" ? w_pre * abs(nu_pre) : 0.0, h_p0_m=Float64(p["h_p0_m"]),
        insertion_node_deg=haskey(p, "insertion_node_deg") && p["insertion_node_deg"] !== nothing ? Float64(p["insertion_node_deg"]) : nothing, t_star_d=p["t_star_d"] === nothing ? nothing : Float64(p["t_star_d"]))
end
"""Calibration directory name of an S-free start (the first pass calibrated to exp(x0))."""
planner_ic_cal_name(gid::AbstractString, x0::Real) = string(gid, "__x0_", round(Float64(x0); digits=5))
"""JSON token of a geometry in the IC file: ODY and G1 by role token, the generic orbits by id."""
planner_ic_token(gid::AbstractString, grole::AbstractString) = grole in ("ODY", "G1") ? String(grole) : String(gid)
"""First-pass ln q the a priori schedule (and the guard's target) is anchored at: the calibrated q_first for "ref"; S-fixed shifts it by the pre-burn offset nu_pre (planner_rev20_replay: x_ref + u0); S-free was calibrated to exp(x0)."""
planner_anchor_ln(q_first::Real, ic) = ic === nothing ? log(q_first) : ic.scenario == "S-fixed" ? log(q_first) + ic.nu_pre : ic.x0
""" Initial offset u0 of the guard's bookkeeping (executed offset excluding the planner's burns): nu_pre for S-fixed, else 0."""
planner_ic_u0(ic) = (ic !== nothing && ic.scenario == "S-fixed") ? ic.nu_pre : 0.0

"""Node-check context of a flight: nothing unless the arm's guard set has node_check; then (; rot, omega_target_deg, J2, Re) with the target's node_deg ([planner.target(s)] node_deg, no default) and the Kaula Mars J2 / R of [planner.node_check]."""
function planner_node_ctx(meta, target, rot, sw, A)
    sw.ttg_guard || return nothing; G = planner_guard_constants(meta; set=String(get(A, "guard_set", "v2")), target=target); G.node_check || return nothing
    nd = get(planner_target_table(meta, target), "node_deg", nothing); nd === nothing && error("guard node_check is on but the target has no node_deg")
    nc = meta["planner"]["node_check"]; return (; rot, omega_target_deg=Float64(nd), J2=Float64(nc["J2"]), Re=Float64(nc["R_m"]))
end

# ---- phase cal --------------------------------------------------------------------------------------------------------------------------------------

"""Calibration case (stageg3_cal_task) for a geometry id including ODY. seed fields as stageg3_cal_task: config_path, gid, gram_tag, out."""
function planner_cal_task(seed)
    planner_meta(seed.config_path)
    return stageg3_cal_task(seed)
end

"""Calibration case of an S-free start: stageg3_cal_task (copied, not edited; Revision 20) with the calibration target heat rate exp(x0) instead of the band's Q_FIRST. seed fields as stageg3_cal_task plus `cal_name`
(planner_ic_cal_name) and `x0`; writes <seed.out>/<cal_name>/{guess,nodes,calibration}.csv. The flight reads that directory instead of <gid>."""
function planner_cal_x0_task(seed)
    t0 = time()
    try
        planner_meta(seed.config_path)
        meta, sc, g, ep = stageg3_load(seed.config_path, seed.gid); scc = stageo_cal_view(sc); ms = stagef_setup(sc)
        b = scc["band"]; b["q_low_w_cm2"] = exp(Float64(seed.x0))^2 / Float64(b["q_high_w_cm2"])                    # gram_calibrate targets sqrt(q_low q_high) = exp(x0)
        epoch = parse_epoch(sc["grid"]["epochs_utc"][1]); rot = mars_rot(ms.planet, epoch)
        et0 = SimulationModel.ephemerides_time_seconds(epoch, SpiceEphemeridesModel()); stagef_env!(false, meta)
        model = stagef_gram_model(sc, epoch, stagef_seed("mars" * seed.gram_tag)); pargs = gram_pargs(ms.planet, Float64(sc["atmosphere_gram"]["entry_interface_km"]))
        inc, raan, omega = Float64(g["inc_deg"]), Float64(g["raan_deg"]), Float64(g["omega_deg"])
        gg = pi_rule_altitude(scc, ms, rot, et0, model, pargs, inc, raan, omega)
        f = gram_exp_fit(sc, model, pargs, gg.h_g, gg.site)
        temp = Float64(sc["atmosphere"]["temperature_k"]); top = Float64(sc["atmosphere"]["fit_top_altitude_m"])
        pwf = pw_fit(f.pw_h, f.pw_rho, temp, "piecewise"); msc = with_fit(ms, exp_fit(f.rho_g, gg.h_g, f.H, top, temp, "cell fit"))
        dir = joinpath(seed.out, seed.cal_name); mkpath(dir)
        seed_alt, recs = gram_calibrate(msc, scc, model, rot, epoch, dir, gg.h_g; inc_deg=inc, raan_deg=raan, omega_deg=omega, pwfit=pwf)
        write_rows_m(joinpath(dir, "calibration.csv"), recs)
        write_rows_m(joinpath(dir, "nodes.csv"), [(h_m=a, rho_kg_m3=b_) for (a, b_) in zip(f.pw_h, f.pw_rho)])
        core_t = model.core.initial_time; stageg3_gram_guard(core_t, ep)
        row = (; planet="mars", tag=seed.cal_name, status="ok", h_guess_m=gg.h_g, seed_alt_m=seed_alt, vrel_first_mps=recs[end].vrel_mps, q_first_w_cm2=recs[end].q_w_cm2, q_target_w_cm2=recs[end].q_target,
            site_lat_rad=gg.site.lat, site_lon_rad=gg.site.lon, fit_H_m=f.H, pw_max_abs_lnrho=f.pw_max_abs_lnrho, hloc_at_hg_m=f.hloc_at_hg, n_cal_iter=length(recs),
            gram_core_epoch=string(core_t.year, "-", core_t.month, "-", core_t.day, " ", core_t.hour, ":", core_t.minute, ":", core_t.second), wall_s=time() - t0, error="", epoch_utc=ep)
        write_rows_m(joinpath(dir, "guess.csv"), [row]); return row
    catch err
        return (; planet="mars", tag=seed.cal_name, status="cal_error", h_guess_m=NaN, seed_alt_m=NaN, vrel_first_mps=NaN, q_first_w_cm2=NaN, q_target_w_cm2=NaN, site_lat_rad=NaN, site_lon_rad=NaN, fit_H_m=NaN,
            pw_max_abs_lnrho=NaN, hloc_at_hg_m=NaN, n_cal_iter=0, gram_core_epoch="", wall_s=time() - t0, error=replace(first(sprint(showerror, err), 400), "\n" => " "), epoch_utc="")
    end
end

# ---- phase run: one flight ------------------------------------------------------------------------------------------------------------------------------

"""Revision 17a: the GRAM constructor keywords of a dispersion level ([planner.dispersion.<name>]; keys gram_perturbation_scales (4 numbers) and mars_mgcm_dust_levels (3 numbers)). "" or "L0" = no keywords (the unchanged flight)."""
function planner_disp_kwargs(meta, name::AbstractString; gid::AbstractString="")
    (isempty(name) || name == "L0") && return (;)
    d = meta["planner"]["dispersion"][name]; kw = Pair{Symbol, Any}[]
    if !isempty(gid) && haskey(d, gid); d = d[gid]; end                                            # Revision 21: levels of the generic geometries are sub-tables [planner.dispersion.<level>.<gid>]
    (haskey(d, "gram_perturbation_scales") || haskey(d, "mars_mgcm_dust_levels")) || error("dispersion level $name has no entry for geometry $gid")
    haskey(d, "gram_perturbation_scales") && push!(kw, :gram_perturbation_scales => Tuple(Float64.(d["gram_perturbation_scales"])))
    haskey(d, "mars_mgcm_dust_levels") && push!(kw, :mars_mgcm_dust_levels => Tuple(Float64.(d["mars_mgcm_dust_levels"])))
    return (; kw...)
end

"""One planner flight. seed fields: config_path, gid, arm (B-ODY key BODY, L0..L3, OL, OLR), tag (GRAM seed string), run_id, out, cal_root, horizon_days, dt_atm_s (NaN = config), smoke_force, decision_speed ("relative"),
grole (role token of the cell). Plain-data row (stageg3 columns plus mode and switches); never throws."""
function planner_run_task(seed)
    t0 = time()
    try
        meta = planner_meta(seed.config_path); PL = meta["planner"]; A = PL["arms"][seed.arm]; mode = Symbol(String(A["mode"]))
        mode in (:body, :online, :ol, :olr, :n2it) || error("arm $(seed.arm) is not a planner arm (mode $mode)")
        meta, sc, geom, ep = stageg3_load(seed.config_path, seed.gid; dt_atm_s=seed.dt_atm_s); G3 = meta["stageg3"]; SO = meta["stageo"]; ms0 = stagef_setup(sc); sw0 = planner_switches(meta, seed.arm)
        inc, raan, omega = Float64(geom["inc_deg"]), Float64(geom["raan_deg"]), Float64(geom["omega_deg"])
        epoch = parse_epoch(ep); rot = mars_rot(ms0.planet, epoch)
        ic_name = String(sw0.ic_scenario); ic = planner_ic_plan(meta, planner_ic_token(seed.gid, seed.grole), ic_name)      # Revision 20: nothing for "ref" (the unchanged flight)
        (ic !== nothing && ic.insertion_node_deg !== nothing && Bool(get(PL["ic"], "use_insertion_node", false))) && (raan = ic.insertion_node_deg)      # Revision 25 (default off): the plan's insertion node replaces the geometry's Omega0 (the first-pass calibration must have been made at it)
        cdir = joinpath(seed.cal_root, ic_name == "S-free" ? planner_ic_cal_name(seed.gid, ic.x0) : seed.gid); G = CSV.read(joinpath(cdir, "guess.csv"), DataFrame)[1, :]; N = CSV.read(joinpath(cdir, "nodes.csv"), DataFrame)
        pwf = pw_fit(N.h_m, N.rho_kg_m3, Float64(sc["atmosphere"]["temperature_k"]), "piecewise"); ms = with_fit(ms0, pwf)
        stagef_env!(true, meta)
        disp = hasproperty(seed, :disp) ? String(seed.disp) : ""; model = stagef_gram_model(sc, epoch, stagef_seed("mars" * seed.tag); planner_disp_kwargs(meta, disp; gid=seed.gid)...); core_t = model.core.initial_time; stageg3_gram_guard(core_t, ep)
        pargs = gram_pargs(ms.planet, Float64(sc["atmosphere_gram"]["entry_interface_km"]))
        dir = joinpath(seed.out, seed.run_id); mkpath(dir); horizon = Float64(seed.horizon_days) * SECONDS_PER_DAY
        state0 = mars_initial_state(ms, rot, Float64(G.seed_alt_m), inc, raan, omega)
        if ic_name == "S-fixed" && ic.dv_pre_signed != 0.0                                                           # Revision 20: impulsive tangential burn at the starting apoapsis before pass 1 (propellant mass not debited: ~2e-2 kg)
            state0 = (state0[1], state0[2] * (1 + ic.dv_pre_signed / norm(state0[2])))
        end
        base = gram_args(ms, sc, model, epoch, horizon, dir; state=state0)
        meas = NamedTuple[]; vlog = NamedTuple[]
        margin = Float64(G3["keeper_margin_m"]); cfg, hq, hl = stagef_corridor_config(ms, Float64(G.vrel_first_mps), margin); c = cfg.values["comparison"]
        inner = CorridorGuidancePolicy(cfg, :favorable, Symbol(meta["stagef"]["policy"]), Float64(c["max_burn_delta_v_mps"]), Float64(c["max_burn_duration_s"]),
            Float64(c["total_delta_v_budget_mps"]), Float64(c["near_passive_deadband_m"]), Float64(c["target_periapsis_altitude_m"]), deg2rad(Float64(c["planning_lead_deg"])),
            0, false, nothing, NamedTuple[], NamedTuple[], String[])
        st = stagef_est_state(ms, epoch); et0e = SimulationModel.ephemerides_time_seconds(epoch, SpiceEphemeridesModel())
        maker = (p, r, v, mass, t) -> (dv -> stagef_predict_alt(st, r, v, Float64(dv), et0e + t))
        w2m = Float64(sc["band"]["w_cm2_to_w_m2"])
        fpol = StageFPolicy(inner, maker, mode === :n2it ? :adaptive : :fixed, pwf, Float64(G.vrel_first_mps), Float64(sc["band"]["q_low_w_cm2"]) * w2m, Float64(sc["band"]["q_high_w_cm2"]) * w2m,
            margin, meas, "", 0.0, 0.0, Float64(meta["analysis"]["adapt_min_band_m"]), NamedTuple[])
        ds = Symbol(seed.decision_speed); ds === :relative || error("planner arms use decision_speed = relative (U15)")
        et0p = SimulationModel.ephemerides_time_seconds(epoch, SpiceEphemeridesModel()); planet = ms.planet
        omega_fn = t -> (SimulationModel.planet_frame_lpi(planet, et0p + t, SpiceEphemeridesModel()))' * SVector{3, Float64}(planet.ω)
        gtarget = planner_geometry_target(meta, seed.gid); P, B = planner_constants(meta, w2m, Float64(G.hloc_at_hg_m), seed.arm; target=gtarget); sw = planner_switches(meta, seed.arm)
        tab = (mode === :body || mode === :n2it) ? nothing : planner_table(meta, seed.gid)
        wrap = mode === :n2it ? planner_n2it_policy(meta, A, fpol, ms, omega_fn, G, P, vlog, sw) : PlannerPolicy(fpol, Float64(ms.harm.gm_m3s2), omega_fn, Float64(G.vrel_first_mps), NaN, mode, sw, P, B, tab, vlog, Float64(G.hloc_at_hg_m), Float64(G.q_first_w_cm2),
            0.0, [0.0], nothing, -Inf, NaN, Float64[], 0.0, NamedTuple[], NaN, 0.0, 0, (sw.fb_term || P.apo_one_sided) ? (fbp = planner_fb_params(meta, seed.gid, ms.harm.gm_m3s2); sw.fb_term ? fbp : (; fbp..., weight=0.0)) : nothing, 0, ic,
            sw.ttg_guard ? planner_ttg_setup(meta, seed.gid, tab, P, ms.harm.gm_m3s2, planner_anchor_ln(Float64(G.q_first_w_cm2), ic); set=String(get(A, "guard_set", "v2")), tstar_d=(ic !== nothing && Bool(get(PL["ic"], "use_plan_t_star", false))) ? ic.t_star_d : nothing) : nothing, Float64[], NaN, 0,
            planner_node_ctx(meta, gtarget, rot, sw, A), nothing)
        v = cfg.values["vehicle"]; thruster = BaseThrusterModel([Float64(v["thrust_n"])], [0.0], [0.0], [-1.0], [-1.0], [Float64(v["specific_impulse_s"])])
        args = SimulationConfiguration(file_paths=base.file_paths, simulation_settings=base.simulation_settings, mission_configuration=base.mission_configuration,
            environment_model=base.environment_model, dynamics_model=base.dynamics_model, guidance_model=GuidanceModel((wrap,), [1.0e12]), navigation_model=base.navigation_model,
            control_model=ControlModel((thruster,), [1.0e12]), initial_time=base.initial_time, integration_tolerances=base.integration_tolerances,
            solver_config=SolverConfig(solver_mode=base.solver_config.solver_mode, maxiters=Int(SO["solver_maxiters"])))
        qstop = seed.smoke_force ? Float64(G3["smoke_forced_safety_q_w_cm2"]) : P.q_never_exceed_qs_w_cm2
        nh = String(A["safety"]) == "q_alt" ? Float64(PL["safety"]["n_hloc"]) : Inf       # n_hloc = Inf: the altitude test (alt0 - n_hloc H) never fires
        ss = SafetyStop(true, qstop, nh, Float64(G.hloc_at_hg_m))
        obs = NamedTuple[]; pk = PassBook(P.ra_stop_m)
        extra = (stageg3_apsis_callback(obs, meas, vlog, args, inner, model, pargs, pwf, pk, ss), stageo_peak_callback(args, model, pargs, pk, Float64(SO["peak_gate_altitude_m"])), _preapoapsis_guidance_callback(inner))
        tr = time(); sol = run_simulation(args; return_solution=true, isolate_state=false, extra_callbacks=extra); run_wall = time() - tr
        book_flush!(pk); by = Dict(r.pass => r for r in pk.recs); vl = Dict(r.pass => r for r in vlog)
        logpath = joinpath(dir, "gram_density_perturbation_log.csv")
        rows0 = stagef_table(sol, args, ms, rot, sc, pwf, logpath); isfile(logpath) && rm(logpath; force=true); isfile(logpath * ".summary.toml") && rm(logpath * ".summary.toml"; force=true)
        rows = [merge(r, (; rho_peak_kg_m3=haskey(by, r.pass) ? by[r.pass].rho_peak_kg_m3 : NaN, rp_m=Float64(r.radius_m), ra_m=haskey(by, r.pass) ? by[r.pass].ra_m : NaN,
                           vp_inertial_mps=haskey(vl, r.pass) ? vl[r.pass].vp_inertial_mps : NaN)) for r in rows0]
        write_rows_m(joinpath(dir, "passes.csv"), rows)
        ve = Float64(sc["vehicle"]["specific_impulse_s"]) * 9.80665; dv = ve * log(ms.wet_mass_kg / Float64(sol.u[end].sc[1].mass)); ic === nothing || (dv += ic.dv_pre_charged); burns = length(inner.burn_plans)
        stopped = !isempty(inner.failures); failures = join(inner.failures, " | ")
        write_rows_m(joinpath(dir, "decisions.csv"), fpol.edge_log); write_rows_m(joinpath(dir, "replan.csv"), wrap.rlog); write_rows_m(joinpath(dir, "measured.csv"), meas)
        ret = string(sol.retcode)
        stop_reason = pk.stopped ? "end_radius" : !isempty(ss.reason) ? ss.reason : stopped ? "controller_stopped" : ret == "Success" ? "horizon" : "retcode_" * ret
        return (; planet="mars", tag=seed.tag, run_id=seed.run_id, arm=seed.arm, ctrl=String(mode), code=planner_arm_code(seed.arm), status="ok", retcode=ret,
            passes=length(rows), t_end_days=sol.t[end] / SECONDS_PER_DAY, horizon_days=seed.horizon_days, completed=pk.stopped && !stopped, controller_stopped=stopped, controller_failures=failures, total_dv_mps=dv, burns=burns,
            collapsed_decisions=0, gram_core_epoch=string(core_t.year, "-", core_t.month, "-", core_t.day, " ", core_t.hour, ":", core_t.minute, ":", core_t.second),
            run_wall_s=run_wall, wall_s=time() - t0, error="", stop_reason=stop_reason,
            geometry=seed.gid, role=seed.grole, keeper=seed.arm, epoch_utc=ep, inc_deg=inc, omega_deg=omega, raan_deg=raan, epoch_off_d=Float64(geom["epoch_off_d"]), decision_speed=seed.decision_speed,
            safety_q_w_cm2=qstop, safety_alt0_m=ss.alt0_m, safety_h_local_m=ss.h_local_m, safety_reason=ss.reason,
            planner_dp_total_s=sum((r.dp_solve_s for r in wrap.rlog); init=0.0) + wrap.apriori_build_s, planner_dp_fraction=(sum((r.dp_solve_s for r in wrap.rlog); init=0.0) + wrap.apriori_build_s) / max(run_wall, 1.0e-9), n_replan_rows=length(wrap.rlog),
            mode=String(mode), slack=sw.slack, ff=sw.ff, early=sw.early, update=sw.update, horizon_passes=sw.horizon, terminal=sw.terminal, nu_cum=wrap.nu_cum, n_decisions=length(fpol.edge_log), disp,
            fb_term=sw.fb_term, ttg_guard=sw.ttg_guard, ic_scenario=ic_name, ic_x0=ic === nothing ? NaN : ic.x0, ic_dv_pre_mps=ic === nothing ? 0.0 : ic.dv_pre_charged)
    catch err
        return (; planet="mars", tag=seed.tag, run_id=seed.run_id, arm=seed.arm, ctrl="", code=planner_arm_code(seed.arm), status="run_error", retcode="", passes=0, t_end_days=NaN, horizon_days=seed.horizon_days,
            completed=false, controller_stopped=false, controller_failures="", total_dv_mps=NaN, burns=0, collapsed_decisions=0, gram_core_epoch="", run_wall_s=NaN, wall_s=time() - t0,
            error=replace(first(sprint(showerror, err), 500), "\n" => " "), stop_reason="error",
            geometry=seed.gid, role=seed.grole, keeper=seed.arm, epoch_utc="", inc_deg=NaN, omega_deg=NaN, raan_deg=NaN, epoch_off_d=NaN, decision_speed=seed.decision_speed,
            safety_q_w_cm2=NaN, safety_alt0_m=NaN, safety_h_local_m=NaN, safety_reason="", planner_dp_total_s=NaN, planner_dp_fraction=NaN, n_replan_rows=0, mode="", slack=false, ff=false, early=false, update=false, horizon_passes=0, terminal=false, nu_cum=NaN, n_decisions=0, disp=hasproperty(seed, :disp) ? String(seed.disp) : "",
            fb_term=false, ttg_guard=false, ic_scenario="", ic_x0=NaN, ic_dv_pre_mps=NaN)
    end
end
