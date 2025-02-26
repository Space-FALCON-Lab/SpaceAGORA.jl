import .config

function Earth_firing_plan(planet, apoapsis, periapsis)
    max_q = 0.6 # Max dynamic pressure
    min_q = 0.5 # Min dynamic pressure
    estimated_ρ = planet.ρ_ref * exp(-(periapsis-6378e3) / planet.H)
    periapsis_velocity = sqrt(planet.μ *(2/periapsis - 2/(apoapsis+periapsis)))
    q = 0.5 * estimated_ρ * periapsis_velocity^2
    if q > max_q
        rp_new = -planet.H*log(2*min_q/(planet.ρ_ref*periapsis_velocity^2)) + 6378e3
        println("Aerobraking at ", periapsis, " m")
        println("New periapsis: ", rp_new, " m")
        args[:delta_v] = abs(sqrt(2*planet.μ/periapsis - planet.μ*2/(periapsis+apoapsis)) - sqrt(2*planet.μ/rp_new - planet.μ*2/(rp_new+apoapsis)))
        args[:phi] = deg2rad(180.0)
    else
        args[:delta_v] = 0.0
        args[:phi] = deg2rad(0.0)
    end
    return args
end