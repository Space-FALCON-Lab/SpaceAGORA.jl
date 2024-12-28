include("../utils/Ref_system_conf.jl")
include("../utils/Reference_system.jl")

using PythonCall

sys = pyimport("sys")
os = pyimport("os")

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(@__FILE__)), "GRAMpy"))

gram = pyimport("gram")

import .config

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
    if config.cnf.drag_state == true
        T = p.T
    else
        T = p.T
    end

    return T
end

# function wind_def(model, indexes, x, table)
#     """

#     """

#     # NEED TO FINISH

#     # boundary 1
#     index_b1 = indexes[0]
#     index_b2 = indexes[1]

#     if model[] == 1
#         wE, wN = interp(table[][index_b1], table[][index_b2], x), interp(table[][index_b1], table[][index_b2], x)
#         wind = [wE, wN, 0]
#     elseif model[] == 2
#         wE, wN, vN = interp(table[][index_b1], table[][index_b2], x), interp(table[][index_b1], table[][index_b2], x), interp(table[][index_b1], table[][index_b2], x)
#         wind = [wE, wN, vN] 
#     end

#     return wind
# end

# function define_model_parameters(OE, p, t0, tf_prev, model, final_state_angle)
#     """

#     """


# end


function density_constant(h, p, OE=0, lat=0, lon=0, timereal=0, t0=0, tf_prev=0, montecarlo=0, Wind=0, args=0, version=[])
    """

    """

    if config.cnf.drag_state == false
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

function density_exp(h, p, OE=0, lat=0, lon=0, timereal=0, t0=0, tf_prev=0, montecarlo=0, Wind=0, args=0, version=[])
    """

    """

    ρ = p.ρ_ref * exp.((p.h_ref .- h)/p.H)

    T = temperature_linear(h, p)

    wind = [0, 0, 0]

    return ρ, T, wind
end

function density_gram(h, p, OE, lat, lon, timereal, t0, tf_prev, montecarlo, Wind, args, version=[], atmosphere=nothing)
    """

    """

    if not version
        version = args.MarsGram_version
    end

    if config.drag_state == False
        rho , T , wind = density_exp(h, p)
    else
        position = gram.Position()
        position.height = h*1e-3
        lat = rad2deg(lat)
        lon = rad2deg(lon)
        position.latitude = lat
        position.longitude = lon
        position.elapsedTime = 0
        atmosphere.setPosition(position)
        # print('set planet position', position.latitude, position.longitude, position.height)
        atmosphere.update()
        # print('update')
        atmos = atmosphere.getAtmosphereState()
        # print('get atmo state')
        rho = atmos.density
        T = atmos.temperature
        wind = [montecarlo ? atmos.perturbedEWWind : atmos.ewWind,
                montecarlo ? atmos.perturbedNSWind : atmos.nsWind,
                montecarlo ? atmos.verticalWind : 0]
    end

    return rho, T, wind
end