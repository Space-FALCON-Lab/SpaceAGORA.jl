"""
Prospective adapter planning reserve. `fraction` reduces enabled physical limits
before planning; it never changes the shared validator. At most one quarter of
this reserve may be consumed by the scale-based Float64 resolution estimate.
The default one-percent reserve is a policy, not a total-acceleration theorem.
Every delivered reference must still pass independent validation.
"""
struct RPOPlanningHeadroom
    fraction::Float64
    function RPOPlanningHeadroom(; fraction=0.01)
        f = Float64(fraction)
        isfinite(f) && 0 < f < 1 || throw(ArgumentError("Planning headroom must be strictly between zero and one."))
        new(f)
    end
end

"""
Compute reduced adapter limits and screen finite-difference output resolution.
The estimates `32eps()*position_scale/dt` and `32eps()*velocity_scale/dt`
reserve budget for stored samples and arithmetic. They are not allowances in
validation and do not certify arbitrary planner calculations. Unsupported
resolution rejects instead of increasing a physical limit. Recheck using the
actual returned positions and velocities before delivering a candidate.
"""
function rpo_planning_budget(request::RPOPlanningRequest, policy::RPOPlanningHeadroom;
                             position_scale_m, velocity_scale_mps)
    c = request.constraints
    vlim = c.max_speed_mps
    alim = c.max_acceleration_mps2
    ps, vs = Float64(position_scale_m), Float64(velocity_scale_mps)
    scales_ok = isfinite(ps) && ps >= 0 && isfinite(vs) && vs >= 0
    ve = (32eps(Float64) * ps) / request.reference_dt_s
    ae = (32eps(Float64) * vs) / request.reference_dt_s
    reserve_v = isnothing(vlim) ? nothing : policy.fraction * vlim
    reserve_a = isnothing(alim) ? nothing : policy.fraction * alim
    planned_v = isnothing(vlim) ? nothing : (1 - policy.fraction) * vlim
    planned_a = isnothing(alim) ? nothing : (1 - policy.fraction) * alim
    speed_ok = isnothing(vlim) || (isfinite(ve) && reserve_v > 0 && 0 < planned_v < vlim && ve <= reserve_v / 4)
    accel_ok = isnothing(alim) || (isfinite(ae) && reserve_a > 0 && 0 < planned_a < alim && ae <= reserve_a / 4)
    supported = scales_ok && speed_ok && accel_ok
    return (supported=supported, reason=supported ? :supported : :insufficient_reference_precision,
        fraction=policy.fraction, max_speed_mps=planned_v, max_acceleration_mps2=planned_a,
        position_scale_m=ps, velocity_scale_mps=vs,
        speed_resolution_mps=ve, acceleration_resolution_mps2=ae,
        speed_reserve_mps=reserve_v, acceleration_reserve_mps2=reserve_a)
end
