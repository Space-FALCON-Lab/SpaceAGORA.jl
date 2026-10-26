# MIT-licensed regression expressions from SpaceAGORA
# bfa32f6a44c5f242714b8cf182a451e863119a05,
# benchmarks/studies/cygnss_reconstruction/cygnss_slew.jl:258-288.
# The complete producer and its tests live in SpaceAGORAMissionReconstructions.
# This synthetic fixture keeps generic output/actuator assertions independent
# of mission data and drivers. No campaign is run by this file.
function slew_save_fields(args, wheels)
    SC = SM.SimulationCallbacks
    ω_of(u) = SVector{3, Float64}(u.sc[1].ω)
    return [
        SC.SaveField(:body_rate_rad_s, (u, t, integ) -> [ω_of(u)]; per_satellite=true, column_prefix="body_rate"),
        SC.SaveField(:wheel_momentum_nms,
            (u, t, integ) -> [SM.wheel_momentum_body(wheels, t)];
            per_satellite=true, column_prefix="wheel_momentum_nms"),
        SC.SaveField(:wheel_torque_nm,
            (u, t, integ) -> [norm(SM.wheel_reaction_torque(SM.wheel_momentum_body(wheels, t),
                SM.wheel_momentum_rate_body(wheels, t), ω_of(u)))];
            per_satellite=true),
    ]
end
