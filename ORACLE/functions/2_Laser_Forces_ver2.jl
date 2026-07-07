# Laser impulse tracker: integrates F/m·dt at every accepted ODE step in RTN frame
Base.@kwdef mutable struct _LaserImpulseTracker
    # ---------------------------------------------------------------------------
    # Laser impulse tracker: integrates F/m·dt at every accepted ODE step so that
    # the reported dV_RTN is the true accumulated laser impulse, not an ECI
    # velocity difference that is contaminated by orbital-phase drift.
    # ---------------------------------------------------------------------------
    t_prev::Float64            = 0.0
    dv_R::Float64              = 0.0
    dv_T::Float64              = 0.0
    dv_N::Float64              = 0.0
    t_hist::Vector{Float64}    = Float64[]
    dv_R_hist::Vector{Float64} = Float64[]
    dv_T_hist::Vector{Float64} = Float64[]
    dv_N_hist::Vector{Float64} = Float64[]
end

# Creates a DiscreteCallback that accumulates the laser impulse in RTN frame at every accepted ODE step.
function _make_laser_impulse_callback(
    model::OpenCavityLaserLinkModel, # the live laser model (isolate_state=false, so no deep copy)
    tracker::_LaserImpulseTracker,   # mutable struct that accumulates dV; lives outside the integrator
    mass_kg::Float64,                # target satellite mass, needed to convert force → acceleration
)
    function affect!(integrator)            # called by DiffEq at every accepted ODE step
        dt = integrator.t - tracker.t_prev  # time elapsed since the last accepted step
        if dt > 0.0                         # skip the very first call where t_prev == t (dt = 0)
            helper_idx = model.active_helper_idx  # which helper is currently firing (0 = none)
            if helper_idx > 0                     # only accumulate when a link is active
                sc = integrator.u.sc              # array of all spacecraft states at this step
                tgt_pos = SVector{3, Float64}(sc[model.target_idx].pos)  # target ECI position [m]
                tgt_vel = SVector{3, Float64}(sc[model.target_idx].vel)  # target ECI velocity [m/s]
                hlp_pos = SVector{3, Float64}(sc[helper_idx].pos)        # active helper ECI position [m]

                # Recompute the laser force inline (laser_link_pair_force is not exported)
                rel = tgt_pos - hlp_pos   # vector from helper → target [m]
                rho = norm(rel)            # separation distance [m]
                if rho > 0.0 && rho <= model.range_m  # link is physically valid
                    # Scalar force magnitude: F = η·β·M·P / c  [N]
                    F_mag = model.eta * model.beta *
                            model.magnification * model.power_w /
                            299_792_458.0
                    force = F_mag * rel / rho              # force vector along line-of-sight [N]
                    rhat, that, nhat = _rtn_basis(tgt_pos, tgt_vel)  # RTN unit vectors at target
                    accel = force / mass_kg                # acceleration vector [m/s²]
                    # Trapezoidal-style: a·dt ≈ ΔV for this step (left-endpoint rectangle rule)
                    tracker.dv_R += dot(accel, rhat) * dt  # radial component [m/s]
                    tracker.dv_T += dot(accel, that) * dt  # along-track component [m/s]
                    tracker.dv_N += dot(accel, nhat) * dt  # cross-track component [m/s]
                end
            end
        end
        tracker.t_prev = integrator.t  # advance the "last step" timestamp for the next call

        # Append current time and cumulative dV to the history vectors for post-run lookup
        push!(tracker.t_hist,    integrator.t)
        push!(tracker.dv_R_hist, tracker.dv_R)
        push!(tracker.dv_T_hist, tracker.dv_T)
        push!(tracker.dv_N_hist, tracker.dv_N)
    end
    return DiscreteCallback(
        (u, t, integrator) -> true,  # condition: fire at every accepted step
        affect!;
        save_positions=(false, false), # don't save an extra solution snapshot at trigger/affect time
    )
end

# Look up the accumulated (dv_R, dv_T, dv_N) at a given time from the history.
function _tracked_dv_at(tracker::_LaserImpulseTracker, t::Float64)
    isempty(tracker.t_hist) && return (0.0, 0.0, 0.0)
    idx = clamp(searchsortedlast(tracker.t_hist, t), 1, length(tracker.t_hist))
    return (tracker.dv_R_hist[idx], tracker.dv_T_hist[idx], tracker.dv_N_hist[idx])
end
