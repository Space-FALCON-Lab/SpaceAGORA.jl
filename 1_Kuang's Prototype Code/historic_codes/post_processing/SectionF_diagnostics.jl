"""
SectionF_diagnostics.jl

Conservation-law validation diagnostics for Section F of the paper:
"Assessing Small Satellite Maneuverability in Proliferated LEO Constellations
 with Open Cavity Laser-Interlinks"

Computes four validation residuals from a saved ODE solution and parameter dict:
  1. Net OCL force        F_net^(OCL)(t)   – should be ≈ 0 when link is active (Eq. 21)
  2. Net OCL torque       τ_net^(OCL)(t)   – should be ≈ 0 when link is active  (Eq. 21)
  3. Energy residual      ε_E(t)           – should stay at integration error   (Eq. 22)
  4. Angular-momentum residual ε_H^(OCL)(t) – should stay ≈ 0                  (Eq. 23)

The helper satellite is assumed to be satellite 1 and the target satellite 2
(index convention: helper=1, target=2).
"""

"""
    compute_sectionF_diagnostics(sol, p; target_idx=2, helper_idx=1) -> NamedTuple

Post-process an ODE solution to produce all four Section-F validation time series.

# Arguments
- `sol`        : ODE solution object with fields `.t` (time) and `.u` (state vector at each step)
- `p`          : parameter dictionary (same `p` passed to `nbody_photon!`)
- `target_idx` : satellite index of the target   (default 2)
- `helper_idx` : satellite index of the helper   (default 1)

# Returns a NamedTuple with fields:
- `t`            : time vector (s)
- `orbit_count`  : cumulative orbit count of the target (dimensionless)
- `F_net_mag`    : magnitude of net OCL force vector at each time step (N)
- `tau_net_mag`  : magnitude of net OCL torque vector at each time step (N·m)
- `delta_Eorb`   : ΔE_orb(t) = E_orb(t) – E_orb(0)  (J)
- `W_OCL`        : accumulated OCL work W_OCL(t)       (J)
- `eps_E`        : energy residual ε_E(t) = ΔE_orb – W_OCL  (J)
- `eps_H_mag`    : magnitude of angular-momentum residual ε_H^(OCL)(t) (kg·m²/s)
- `Hz_total`     : z-component of total satellite angular momentum H_z(t) (kg·m²/s)
"""
function compute_sectionF_diagnostics(sol, p; target_idx::Int=2, helper_idx::Int=1)

    N       = p[:N]
    masses  = p[:masses]
    mu      = p[:mu]
    nsteps  = length(sol.t)

    # ── Helpers ────────────────────────────────────────────────────────────────
    _r(u, i) = @SVector [u[idx(i,1)], u[idx(i,2)], u[idx(i,3)]]
    _v(u, i) = @SVector [u[idx(i,4)], u[idx(i,5)], u[idx(i,6)]]

    # ── Orbit count for the target satellite ──────────────────────────────────
    function _target_period(u)
        r = _r(u, target_idx)
        v = _v(u, target_idx)
        rv_norm = norm(r)
        v2      = dot(v, v)
        ε       = 0.5*v2 - mu/rv_norm          # specific orbital energy
        a       = -mu / (2*ε)
        return 2π * sqrt(max(a, 1.0)^3 / mu)
    end

    T_period = [_target_period(sol.u[k]) for k in 1:nsteps]
    dt_vec   = diff(sol.t)
    T_mid    = 0.5 .* (T_period[1:end-1] .+ T_period[2:end])
    orbit_increments = dt_vec ./ T_mid
    orbit_count = cumsum([0.0; orbit_increments])

    # ── Initial total orbital energy ──────────────────────────────────────────
    _Eorb(u) = sum(orbital_energy(u, masses, mu))
    Eorb_0   = _Eorb(sol.u[1])

    # ── Pre-allocate output arrays ─────────────────────────────────────────────
    F_net_mag   = zeros(nsteps)
    tau_net_mag = zeros(nsteps)
    delta_Eorb  = zeros(nsteps)
    W_OCL       = zeros(nsteps)          # cumulative OCL work
    eps_H_x     = zeros(nsteps)          # angular-momentum residual components
    eps_H_y     = zeros(nsteps)
    eps_H_z     = zeros(nsteps)
    Hz_total    = zeros(nsteps)

    # ── Evaluate at every saved time step ─────────────────────────────────────
    # W_OCL is integrated with the trapezoidal rule over saved time steps.
    # The angular-momentum residual uses the same trapezoidal rule.

    # Power and torque integrands at each step
    power_integrand      = zeros(nsteps)   # sum_i Q_i · v_i
    dH_target_integrand  = zeros(3, nsteps) # r_target × F_{target←helper}
    dH_helper_integrand  = zeros(3, nsteps) # r_helper × F_{helper←target}

    for k in 1:nsteps
        u  = sol.u[k]
        Fdict, _ = laser_forces(u, p)      # Dict (i,j) => force on j due to i

        # --- Net OCL force and torque (Eq. 21) --------------------------------
        F_net   = @MVector zeros(3)
        tau_net = @MVector zeros(3)
        for i in 1:N
            # Sum all forces acting ON satellite i
            F_on_i = @MVector zeros(3)
            for ((src, tgt), Fvec) in Fdict
                if tgt == i
                    F_on_i .+= Fvec
                end
            end
            F_net   .+= F_on_i
            r_i      = _r(u, i)
            tau_net .+= cross(r_i, SVector(F_on_i[1], F_on_i[2], F_on_i[3]))
        end
        F_net_mag[k]   = norm(F_net)
        tau_net_mag[k] = norm(tau_net)

        # --- Power integrand (for W_OCL) --------------------------------------
        pow = 0.0
        for i in 1:N
            v_i = _v(u, i)
            for ((src, tgt), Fvec) in Fdict
                if tgt == i
                    pow += dot(Fvec, v_i)
                end
            end
        end
        power_integrand[k] = pow

        # --- Angular-momentum-residual integrands (Eq. 23) --------------------
        # Force on target due to helper
        key_th = (helper_idx, target_idx)   # force on target_idx due to helper_idx
        key_ht = (target_idx, helper_idx)   # force on helper_idx due to target_idx

        r_target = _r(u, target_idx)
        r_helper = _r(u, helper_idx)

        F_on_target = haskey(Fdict, key_th) ? Fdict[key_th] : SVector(0.0, 0.0, 0.0)
        F_on_helper = haskey(Fdict, key_ht) ? Fdict[key_ht] : SVector(0.0, 0.0, 0.0)

        torque_target = cross(r_target, F_on_target)
        torque_helper = cross(r_helper, F_on_helper)

        dH_target_integrand[:, k] = torque_target
        dH_helper_integrand[:, k] = torque_helper

        # --- ΔE_orb -----------------------------------------------------------
        delta_Eorb[k] = _Eorb(u) - Eorb_0

        # --- Total Hz (z-component of angular momentum) -----------------------
        H_z = 0.0
        for i in 1:N
            r_i = _r(u, i)
            v_i = _v(u, i)
            H_z += masses[i] * (r_i[1]*v_i[2] - r_i[2]*v_i[1])
        end
        Hz_total[k] = H_z
    end

    # ── Integrate W_OCL and ε_H cumulatively (trapezoidal rule) ───────────────
    cumW   = 0.0
    cumH_x = 0.0; cumH_y = 0.0; cumH_z = 0.0

    for k in 1:nsteps-1
        dt = sol.t[k+1] - sol.t[k]

        # Work
        cumW      += 0.5 * (power_integrand[k] + power_integrand[k+1]) * dt
        W_OCL[k+1] = cumW

        # Angular-momentum residual
        for dim in 1:3
            avg = 0.5 * (dH_target_integrand[dim, k] + dH_target_integrand[dim, k+1] +
                         dH_helper_integrand[dim, k] + dH_helper_integrand[dim, k+1])
            if dim == 1; cumH_x += avg * dt
            elseif dim == 2; cumH_y += avg * dt
            else; cumH_z += avg * dt
            end
        end
        eps_H_x[k+1] = cumH_x
        eps_H_y[k+1] = cumH_y
        eps_H_z[k+1] = cumH_z
    end

    # ── Energy residual ε_E(t) = ΔE_orb(t) – W_OCL(t) ────────────────────────
    eps_E = delta_Eorb .- W_OCL

    # ── Angular-momentum residual magnitude ────────────────────────────────────
    eps_H_mag = sqrt.(eps_H_x.^2 .+ eps_H_y.^2 .+ eps_H_z.^2)

    return (
        t            = sol.t,
        orbit_count  = orbit_count,
        F_net_mag    = F_net_mag,
        tau_net_mag  = tau_net_mag,
        delta_Eorb   = delta_Eorb,
        W_OCL        = W_OCL,
        eps_E        = eps_E,
        eps_H_mag    = eps_H_mag,
        Hz_total     = Hz_total,
    )
end
