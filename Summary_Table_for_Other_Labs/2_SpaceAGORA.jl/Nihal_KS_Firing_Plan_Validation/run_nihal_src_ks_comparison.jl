#!/usr/bin/env julia
#
# Reproduces the Nihal_new_firing_plan schedulers using Nihal's own top-level
# src/ KS-regularized (Kustaanheimo-Stiefel) propagator instead of SpaceAGORA's
# OpenCavityLaserLinkModel, so the two can be compared directly on identical
# initial conditions and schedules. Two-body motion is captured exactly by the
# KS regularization; only J2 + laser are numerically integrated perturbations.
#
# Laser force uses Kuang's canonical open-cavity formula F = B*Pin/c (no extra
# factor), calibrated here via cr=magnification=100 so it matches ORACLE's
# eta=beta=1, magnification=100, power_w=10_000, mass_kg=227 exactly.

using LinearAlgebra, StaticArrays, OrdinaryDiffEq, CSV, DataFrames, Printf

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SRC_ROOT  = joinpath(REPO_ROOT, "src")
const PLAN_DIR  = joinpath(@__DIR__, "Input", "Schedule_v2")

include(joinpath(SRC_ROOT, "constants.jl"))
include(joinpath(SRC_ROOT, "ks_common.jl"))
include(joinpath(SRC_ROOT, "ks_forces.jl"))
include(joinpath(SRC_ROOT, "ks_propagation.jl"))
include(joinpath(SRC_ROOT, "orbital_setup.jl"))

# Canonical Kuang production laser parameters (matches ORACLE's magnification=100,
# power_w=10_000, mass_kg=227; Kuang's real Fint=B*Pin/c has no extra factor).
const CR = 100.0
const LASER_POWER_W = 10_000.0
const MASS_KG = 227.0
const LINK_MAX_RANGE_KM = 200.0
const ATMOSPHERE_TOP_KM = 100.0

# Solves Kepler's equation M = E - e sin(E) by Newton iteration, then converts E to true anomaly.
function mean_to_true_anomaly_deg(M_deg::Float64, e::Float64; tol::Float64=1e-13, max_iter::Int=30)
    M_norm = mod(deg2rad(M_deg) + pi, 2pi) - pi
    E = e < 0.8 ? M_norm : pi
    for _ in 1:max_iter
        f  = E - e * sin(E) - M_norm
        fp = 1.0 - e * cos(E)
        Δ  = f / fp
        E -= Δ
        abs(Δ) < tol && break
    end
    ν = 2.0 * atan(sqrt((1.0 + e) / (1.0 - e)) * tan(E / 2.0))
    return rad2deg(mod(ν, 2pi))
end

# Standalone rv->classical-elements (avoids pulling in orbital_element_tracking_scheduling.jl's
# ForwardDiff/JuMP/SCIP dependencies for a single utility function).
function oe_from_rv(x::AbstractVector{<:Real}; mu::Real)
    position = x[1:3]
    velocity = x[4:6]
    r = norm(position)
    v2 = dot(velocity, velocity)
    energy = 0.5 * v2 - mu / r
    a = -mu / (2.0 * energy)
    h_vec = cross(position, velocity)
    e_vec = (cross(velocity, h_vec) ./ mu) .- (position ./ r)
    e = norm(e_vec)
    return a, e
end

function load_plan(sched_name::String)
    ic_df = CSV.read(joinpath(PLAN_DIR, "initial_conditions.csv"), DataFrame)
    ic_df.satellite = Int.(round.(ic_df.satellite))
    sched_df = CSV.read(joinpath(PLAN_DIR, "schedule_$(sched_name).csv"), DataFrame)

    helper_cols = names(sched_df)[2:end]
    col_matches = match.(r"^(\d+)_to_(\d+)$", helper_cols)
    col_helper_ids = [parse(Int, m.captures[1]) for m in col_matches]
    target_csv_id = parse(Int, col_matches[1].captures[2])
    helper_csv_ids = sort(unique(col_helper_ids))
    ordered_csv_ids = vcat([target_csv_id], helper_csv_ids)
    pos_of = Dict(cid => k for (k, cid) in enumerate(ordered_csv_ids))

    n_intervals = nrow(sched_df)
    schedule_starts = Float64.(sched_df.start_time_s)
    active_helper_pos = zeros(Int, n_intervals)   # 0 = no helper active; else satellite position
    for (ci, cid) in enumerate(col_helper_ids)
        col = sched_df[!, helper_cols[ci]]
        slot_pos = pos_of[cid]
        for k in 1:n_intervals
            col[k] == 1.0 && (active_helper_pos[k] = slot_pos)
        end
    end
    dt = n_intervals >= 2 ? (schedule_starts[end] - schedule_starts[1]) / (n_intervals - 1) : 0.0
    mission_time_s = schedule_starts[end] + dt

    return (ic_df=ic_df, ordered_csv_ids=ordered_csv_ids,
            schedule_starts=schedule_starts, active_helper_pos=active_helper_pos,
            mission_time_s=mission_time_s, n_helpers=length(helper_csv_ids))
end

function build_u0(plan)
    N = plan.n_helpers + 1
    u0 = zeros(Float64, 10 * N)
    for (pos, cid) in enumerate(plan.ordered_csv_ids)
        row = only(findall(==(cid), plan.ic_df.satellite))
        r = plan.ic_df[row, :]
        nu_deg = mean_to_true_anomaly_deg(r.M_deg, r.e)
        r_eci, v_eci = coe_to_rv(
            r.a_km, r.e, deg2rad(r.i_deg), deg2rad(r.raan_deg), deg2rad(r.omega_deg), deg2rad(nu_deg);
            mu=MU_EARTH_KM3_S2,
        )
        u0[ks_state_slice(pos)] .= cartesian_to_ks_state(r_eci, v_eci; mu=MU_EARTH_KM3_S2, t0=0.0)
    end
    return u0, N
end

function make_rhs(plan, N::Int; enable_laser::Bool)
    schedule_starts = plan.schedule_starts
    active_helper_pos = plan.active_helper_pos
    return function rhs!(du, u, p, t)
        active_pos = 0
        if enable_laser && !isempty(schedule_starts)
            _, _, _, t_target = ks_state_components(u, 1)
            k = clamp(searchsortedlast(schedule_starts, t_target), 1, length(schedule_starts))
            active_pos = active_helper_pos[k]
        end
        pair_list = active_pos > 0 ? [(active_pos, 1)] : Tuple{Int,Int}[]
        pair_data = ks_pair_geometry_data(
            u; pair_list=pair_list, link_max_range_km=LINK_MAX_RANGE_KM, atmosphere_top_km=ATMOSPHERE_TOP_KM,
        )
        laser_accels = ks_total_laser_accelerations(
            pair_data, N; control_amplitudes=ones(N), cr=CR,
            laser_power_w=LASER_POWER_W, satellite_mass_kg=MASS_KG,
        )
        for sat in 1:N
            p_sat, q_sat, h_sat, _ = ks_state_components(u, sat)
            r_sat = ks_position(p_sat)
            total_accel = ks_j2_acceleration(r_sat; mu=MU_EARTH_KM3_S2, j2=J2_EARTH, earth_radius_km=R_EARTH_KM) .+
                          laser_accels[:, sat]
            ks_satellite_rhs_sundman_time!(view(du, ks_state_slice(sat)), p_sat, q_sat, h_sat, total_accel)
        end
    end
end

function propagate(plan, u0, N::Int; enable_laser::Bool)
    tf = plan.mission_time_s
    initial_min_radius = minimum(ks_radius_scalar(ks_state_components(u0, sat)[1]) for sat in 1:N)
    s_upper = max(tf / initial_min_radius * 2.5, 1e-6)
    rhs! = make_rhs(plan, N; enable_laser=enable_laser)
    condition(u, s, integrator) = ks_min_physical_time(u) - tf
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!; save_positions=(false, false))
    prob = ODEProblem(rhs!, copy(u0), (0.0, s_upper))
    return solve(prob, Tsit5(); reltol=1e-12, abstol=1e-12, callback=cb, save_everystep=false)
end

function run_scheduler(name::String)
    plan = load_plan(name)
    u0, N = build_u0(plan)
    sol_ctrl = propagate(plan, u0, N; enable_laser=true)
    sol_ref  = propagate(plan, u0, N; enable_laser=false)

    r_c, v_c, _, _ = ks_cartesian_components(sol_ctrl.u[end], 1)
    r_r, v_r, _, _ = ks_cartesian_components(sol_ref.u[end], 1)
    a_c, e_c = oe_from_rv(vcat(r_c, v_c); mu=MU_EARTH_KM3_S2)
    a_r, e_r = oe_from_rv(vcat(r_r, v_r); mu=MU_EARTH_KM3_S2)

    Δa_km = a_c - a_r
    Δe    = e_c - e_r
    @printf("[%s] Nihal-src KS propagation: Δa = %+.4f km   Δe = %+.6e   (ctrl=%s, ref=%s)\n",
        name, Δa_km, Δe, sol_ctrl.retcode, sol_ref.retcode)
    return Δa_km, Δe
end

for name in ["Max_a", "Max_combined", "Max_i", "Max_omega"]
    run_scheduler(name)
end
