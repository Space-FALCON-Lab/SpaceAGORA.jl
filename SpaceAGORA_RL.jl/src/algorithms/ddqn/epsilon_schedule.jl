Base.@kwdef struct EpsilonSchedule
    start::Float64 = 1.0
    stop::Float64 = 0.01
    decay_steps::Int = 500_000
    decay_start_step::Int = 10_000
end

function epsilon_value(schedule::EpsilonSchedule, step::Integer)
    step <= schedule.decay_start_step && return schedule.start
    schedule.decay_steps <= 0 && return schedule.stop
    frac = min(1.0, (step - schedule.decay_start_step) / schedule.decay_steps)
    return max(schedule.stop, schedule.start + (schedule.stop - schedule.start) * frac)
end
