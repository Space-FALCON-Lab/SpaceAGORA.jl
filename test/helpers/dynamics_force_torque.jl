# A dynamic effector's force and torque on spacecraft `i` with state view `x`, as
# the dynamics evaluate it at the run's current time.
function dynamics_force_torque(model, x, p, i::Integer)
    engine = SpaceAGORA.SimulationEngine
    sample = engine.build_state_sample(x, p.args.dynamics_model.spacecraft[i],
        p.args.mission_configuration.orientation_sim, Int(i))
    return engine._evaluate_dynamic_effector(model, x, sample, p, Int(i),
        Float64(p.shared_buffers.current_time[]))
end
