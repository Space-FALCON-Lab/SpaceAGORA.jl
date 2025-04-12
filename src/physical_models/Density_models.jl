include("../utils/Reference_system.jl")

using PythonCall

sys = pyimport("sys")

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

    return 0.0, T, wind
end

function density_exp(h, p, OE=0, lat=0, lon=0, timereal=0, t0=0, tf_prev=0, montecarlo=0, Wind=0, args=0, version=[])
    """

    """

    ρ = p.ρ_ref * exp.((p.h_ref .- h)/p.H)

    T = temperature_linear(h, p)

    wind = [0, 0, 0]

    return ρ, T, wind
end

function density_poly(h, p, OE=0, lat=0, lon=0, timereal=0, t0=0, tf_prev=0, montecarlo=0, Wind=0, args=0, version=[])
    """

    """

    if lowercase(p.name) == "mars"
        polyfit = [4.66697190e-30, -2.43123427e-26,  5.57787983e-23, -7.40000386e-20,
                   6.26120522e-17, -3.50802246e-14,  1.30275628e-11, -3.08201353e-09,
                   4.07513114e-07, -1.47695990e-05, -2.78993562e-03,  1.68510666e-01,
                   -1.10000029e+01];
    elseif lowercase(p.name) == "venus"
        polyfit = [5.70439037e-29, -2.99330044e-25,  6.95564356e-22, -9.42033628e-19,
                   8.23282558e-16, -4.85359445e-13,  1.95869719e-10, -5.37212138e-08,
                   9.72466921e-06, -1.09642572e-03,  7.00391655e-02, -2.38273462e+00,
                   3.47398091e+01]
    elseif lowercase(p.name) == "earth"
        polyfit = [3.10016712e-50, -2.17809613e-46,  6.63624570e-43, -1.11280155e-39,
                   1.02387057e-36, -2.99949753e-34, -4.14125039e-31,  5.08494877e-28,
                   -1.08434957e-25, -2.60445523e-22,  3.36457782e-19, -2.24280648e-16, 
                   9.99503253e-14, -3.18276016e-11,  7.37905396e-09, -1.24091526e-06, 
                   1.48245690e-04, -1.20911350e-02,  6.29822909e-01, -1.87314469e+01,
                   2.33755356e+02]
    elseif lowercase(p.name) == "titan"
        polyfit = [-9.70621744e-31,  4.69037546e-27, -9.80609148e-24,  1.15380789e-20,
                   -8.26238796e-18,  3.55022155e-15, -7.45004592e-13, -5.77924371e-11,
                   8.39394837e-08, -2.43825123e-05,  3.64437354e-03, -3.14678393e-01,
                   8.17796838e+00]
    end

    if length(h) == 1
        power = zeros(length(polyfit))
        for i = 1:length(polyfit)
            power[i] = h^(length(polyfit) - i)
        end

        exponent = dot(polyfit, power)
        ρ = exp(exponent)
    else
        ρ = zeros(length(h))

        for j = 1:length(h)
            power = zeros(length(polyfit))
            for i = 1:length(polyfit)
                power[i] = h[j]^(length(polyfit) - i)
            end

            exponent = dot(polyfit, power)
            ρ[j] = exp(exponent)
        end
    end

    return ρ
end

# @pyexec """
# def density_gram(h, p, OE, lat, lon, timereal, t0, tf_prev, montecarlo, Wind, args, elapsed time, version=[], atmosphere=None):

#     if not version:
#         version = args[:MarsGram_version]

#     if config.drag_state == False:
#         rho , T , wind = density_exp(h, p)
#     else:
#         position = gram.Position()
#         position.height = h*1e-3
#         lat = rad2deg(lat)
#         lon = rad2deg(lon)
#         position.latitude = lat
#         position.longitude = lon
#         position.elapsedTime = 0 
#         atmosphere.setPosition(position)
#         atmosphere.update()
#         atmos = atmosphere.getAtmosphereState()
#         rho = atmos.density
#         T = atmos.temperature
#         wind = [atmos.perturbedEWWind if montecarlo else atmos.ewWind,
#                 atmos.perturbedNSWind if montecarlo else atmos.nsWind,
#                 atmos.verticalWind if montecarlo else 0]

#     return rho, T, wind """ => density_gram

function density_gram(h::Float64, p, lat::Float64, lon::Float64, montecarlo::Bool, Wind::Bool, args::Dict, el_time::Float64, atmosphere=nothing, gram=nothing)
    """

    """

    # sys.path.append(args[:directory_Gram])

    # gram = pyimport("gram")
    # println("Lat: $lat, Lon: $lon, Alt: $h")
    if config.cnf.drag_state == false && args[:keplerian] == false
        rho , T , wind = density_exp(h, p)
        rho = 0.0
    elseif config.cnf.drag_state == true || args[:keplerian] == true
        position = gram.Position()
        position.height = h * 1e-3
        # println("Height: ", position.height)
        # println("Lat: $lat, Lon: $lon, Alt: $h")
        lat = rad2deg(lat)
        lon = rad2deg(lon)
        position.latitude = lat
        # position.longitude = 165 + 2/(24*60*60)*el_time
        position.longitude = lon
        
        position.elapsedTime = el_time # Time since start in s
        atmosphere.setPosition(position)
        # if p.name == "mars"   
        #     position.height += atmosphere.getPosition().surfaceHeight
        #     atmosphere.setPosition(position)
        # end
        # println("set planet position ", position.latitude, position.longitude, position.height)
        # sleep(1.0)
        atmosphere.update()
        # println("update")
        # sleep(1.0)
        atmos = atmosphere.getAtmosphereState()
        # println("get atmo state")
        # sleep(1.0)
        rho = atmos.density
        T = atmos.temperature
        wind = [montecarlo ? atmos.perturbedEWWind : atmos.ewWind,
                montecarlo ? atmos.perturbedNSWind : atmos.nsWind,
                atmos.verticalWind]
    end

    return rho, T, wind
end