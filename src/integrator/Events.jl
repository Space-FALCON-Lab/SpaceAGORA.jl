

function event(solution)
    # Impact definition
    breaker = true
    if length(solution.t_events[end]) != 0
        println("IMPACT!")
        breaker = false
    elseif length(solution.t_events[end-1]) != 0
        breaker = false
        println("PERIAPSIS GREATER THAN APOAPSIS!")
    end
    
    return breaker
end