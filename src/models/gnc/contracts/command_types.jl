module CommandTypes

export PropulsiveManeuverCommand, PropulsiveBurnPlan, AerobrakingControlCommand

Base.@kwdef struct PropulsiveManeuverCommand
    valid::Bool = false
    delta_v_mps::Float64 = 0.0
    direction_rad::Float64 = 0.0
    source_orbit::Int64 = 0
end

"""
    PropulsiveBurnPlan(; valid=false, delta_v_mps=0.0, direction_rad=0.0, source_orbit=0, start_burn_s=-1.0,
                       stop_burn_s=-1.0, thrust_n=0.0, isp_s=0.0, commanded_impulse_n_s=0.0,
                       propellant_required_kg=0.0)

One constant-thrust burn for the burn scheduler (`BaseThrusterModel`): a burn of
`delta_v_mps` (m/s), prograde for `direction_rad = 0` and retrograde for `π`, fired
at `thrust_n` (N) with specific impulse `isp_s` (s) on the window
`[start_burn_s, stop_burn_s)` (s into the run), taking `propellant_required_kg`
(kg) for an impulse of `commanded_impulse_n_s` (N·s). `source_orbit` is the orbit
number of the command it came from, `-1` when none. Only a `valid` plan is executed.
The scheduler plans and executes these; a planner of the burn-plan family
(`BurnPlanFamily`) proposes one, with the propellant, impulse and window duration
of `burn_plan_sizing`.
"""
Base.@kwdef struct PropulsiveBurnPlan
    valid::Bool = false
    delta_v_mps::Float64 = 0.0
    direction_rad::Float64 = 0.0
    source_orbit::Int64 = 0
    start_burn_s::Float64 = -1.0
    stop_burn_s::Float64 = -1.0
    thrust_n::Float64 = 0.0
    isp_s::Float64 = 0.0
    commanded_impulse_n_s::Float64 = 0.0
    propellant_required_kg::Float64 = 0.0
end

Base.@kwdef struct AerobrakingControlCommand
    alpha_command::Float64 = 0.0
end

end # module CommandTypes
