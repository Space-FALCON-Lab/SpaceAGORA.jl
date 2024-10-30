include("../config.jl")


function interp(a, b, x)
    """

    """
    
    # check delta == diff b and a
    if (abs(b-a) > 20.0)
        if b <= 360.0 && b >= 350.0
            b = 360.0 - b
        elseif a <= 360.0 && a >= 350.0
            a = 360.0 - a
        end
    end

    value = x * (b - a) + a

    return value
end

function temperature_linear(h, p)
    """

    """

    # into atmosphere
    if config.drag_state == true
        T = p.T
    else
        T = p.T
    end

    return T
end

function wind_def(model, indexes, x, table)
    """

    """

    # boundary 1
    index_b1 = indexes[0]
    index_b2 = indexes[1]

    if model[] == 1
        wE, wN = interp(table[][index_b1], table[][index_b2], x), interp(table[][index_b1], table[][index_b2], x)
        wind = [wE, wN, 0]
    elseif model[] == 2
        wE, wN, vN = interp(table[][index_b1], table[][index_b2], x), interp(table[][index_b1], table[][index_b2], x), interp(table[][index_b1], table[][index_b2], x)
        wind = [wE, wN, vN] 
    end

    return wind
end

function define_model_parameters(OE, p, t0, tf_prev, model, final_state_angle)
    """

    """


end


function density_constant(h, p, OE=0, lat=0, lon=0, timereal=0, t0=0, tf_prev=0, montecarlo=0, Wind=0, args=0, version=[])
    """

    """

    if config.drag_state == false
        ρ = 0.0
    else
        ρ = p.ρ_ref
    end
    
    T = temperature_linear(h, p)

    wind = [0, 0, 0]

    return ρ, T, wind
end

function density_no(h, p, OE=0, lat=0, lon=0, timereal=0, t0=0, tf_prev=0, montecarlo=0, Wind=0, args=0, version=[])
    """

    """

    T = temperature_linear(h, p)

    wind = [0, 0, 0]

    return ρ, T, wind
end

function density_exp(h, ρ, OE=0, lat=0, lon=0, timereal=0, t0=0, tf_prev=0, montecarlo=0, Wind=0, args=0, version=[])
    """

    """

    ρ = p.ρ_ref * exp((p.h_ref - h)./p.H)

    T = temperature_linear(h, p)

    wind = [0, 0, 0]

    return ρ, T, wind
end

function mars_gram(h, p, OE, lat, lon, timereal, t0, tf_prev, montecarlo, Wind, args, version=[])
    """

    """

    

end