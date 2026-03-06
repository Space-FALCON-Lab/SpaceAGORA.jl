
function event(count_impact, count_apoapsisgreaterpreiapsis)
    # Impact definition
    breaker = true
    if count_impact != 0
        println("IMPACT!")
        breaker = false
    elseif count_apoapsisgreaterpreiapsis != 0
        breaker = false
        println("PERIAPSIS GREATER THAN APOAPSIS!")
    end
    
    return breaker
end