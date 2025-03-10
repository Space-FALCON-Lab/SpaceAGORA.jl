
function ic_calculation_rptoae(planet, γ, v, args)
    if Bool(args[:drag_passage])
        h_0 = args[:EI] * 1e3
    elseif args[:body_shape] == "Blunted Cone"
        h_0 = args[:EI] * 1e3
    end

    # function func(lat)
    #     f = (planet.Rp_e - planet.Rp_p) / planet.Rp_e
    #     e = 1 - (1 - f)^2
    #     N = p_class.Rp_e / sqrt(1 - e*sin(lat)^2)

        

    #     return 
    # end

    # find_zero(lat -> func(lat), (-pi/2, pi/2), Bisection())

    r = planet.Rp_e + h_0                   # Drag passage always start and end at EI km of altitude
    a = planet.μ / ((2*planet.μ/r) - v^2)
    h = r*v*cos(deg2rad(γ))
    p = (h^2)/planet.μ
    e = sqrt(1 - p/a)
    ra = a*(1+e)
    hp = (a*(1-e) - planet.Rp_e)

    # println("hp1: $(hp)")

    # Ri, Vi = orbitalelemtorv([a, e, args[:inclination], args[:Ω], args[:ω], 0.0], planet)

    # date_initial = from_utc(DateTime(args[:year], args[:month], args[:day], args[:hours], args[:minutes], args[:secs]))

    # Rp, Vp = r_intor_p(Ri, Vi, planet, 0, 0, date_initial, 0)

    # hp, _, _ = rtolatlong(Rp, planet)

    # println("hp2: $(hp)")

    if hp < 0
        println("WARNING AT initial_cond_calc: ALTITUDE PERIAPSIS < 0!")
    end

    return ra, hp
end