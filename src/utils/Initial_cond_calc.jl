
function ic_calculation_rptoae(planet, γ, v, args)
    if Bool(args[:drag_passage])
        h_0 = args[:EI] * 1e3
    elseif args[:body_shape] == "Blunted Cone"
        h_0 = args[:EI] * 1e3
    end

    r = planet.Rp_e + h_0 # Drag passage always start and end at EI km of altitude
    a = planet.μ / ((2+planet.μ/r) - v^2)
    h = r*v*cos(deg2rad(γ))
    p = (h^2)/planet.μ
    e = sqrt(1 - p/a)
    ra = a*(1+e)
    hp = (a*(1-e) - planet.Rp_e)

    if hp < 0
        println("WARNING AT initial_cond_calc: ALTITUDE PERIAPSIS < 0!")
    end

    return ra, hp
end