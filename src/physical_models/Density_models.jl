include("../utils/Reference_system.jl")

using PythonCall
using SatelliteToolbox
# SpaceIndices.init()
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
    if config.cnf().drag_state == true
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

    if config.cnf().drag_state == false
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

function density_polyfit(h, p, lat::Float64=0.0, lon::Float64=0.0, montecarlo::Bool=false, Wind::Bool=false, args::Dict=Dict(), el_time::Float64=0.0, atmosphere=nothing, gram=nothing)
    """
    Calculate the density using polynomial fit coefficients.

    Parameters
    ----------
    h : Float64
        Height in meters.
    p : Planet struct
        Planet data
    OE : Any
        Orbital elements (not used in this function).
    lat : Float64
        Latitude in radians (not used in this function).
    lon : Float64
        Longitude in radians (not used in this function).
    timereal : Any
        Real time (not used in this function).
    t0 : Any
        Initial time (not used in this function).
    tf_prev : Any
        Previous final time (not used in this function).
    montecarlo : Bool
        Whether to use Monte Carlo simulation (not used in this function).
    Wind : Bool
        Whether to include wind effects (not used in this function).
    args : Any
        Additional arguments (not used in this function).
    version : Any
        GRAM Version information (not used in this function).
    
    Returns
    -------
    ρ : Float64
        Density in kg/m³.
    T : Float64
        Temperature in Kelvin.
    wind : Vector{Float64}
        Wind vector in m/s.
    """

    if typeof(h) != Float64
        polyfit = p.polyfit_coeffs
        power = zeros(length(polyfit),length(h))
        # Convert height from meters to kilometers
        h = h * 1e-3
        # Calculate the polynomial value at height h
        for i=1:length(polyfit)
            power[i,:] = (h).^(length(polyfit)-i)
        end

        # Calculate the exponent term of the density using the polynomial coefficients
        exponent = zeros(length(h))
        for j=1:length(h)
            exponent[j] = sum(polyfit .* power[:,j])
        end

        # Calculate the density
        ρ = exp.(exponent)
        T = temperature_linear(h, p)

        wind = [0, 0, 0]

        return ρ, T, wind
    else
        polyfit = p.polyfit_coeffs
        power = zeros(length(polyfit))
        # Convert height from meters to kilometers
        h = h * 1e-3
        # Calculate the polynomial value at height h
        for i=1:length(polyfit)
            power[i] = (h)^(length(polyfit)-i)
        end
        # Calculate the exponent term of the density using the polynomial coefficients
        exponent = sum(polyfit .* power)
        # Calculate the density
        ρ = exp(exponent)
        T = temperature_linear(h, p)

        wind = [0, 0, 0]

        return ρ, T, wind
    end
end

function density_gram(h::Float64, p, lat::Float64, lon::Float64, montecarlo::Bool, Wind::Bool, args::Dict, el_time::Float64, atmosphere=nothing, gram=nothing)
    """

    """
    if config.cnf().drag_state == false && args[:keplerian] == false
        if h > 2000.0e3

            rho = 0.0
            T = temperature_linear(h, p)
            wind = [0.0, 0.0, 0.0]
        else
            rho, T, wind = density_polyfit(h, p)
        end
    elseif config.cnf().drag_state == true || args[:keplerian] == true
        if h > 2000.0e3

            rho = 0.0
            T = temperature_linear(h, p)
            wind = [0.0, 0.0, 0.0]
        else
            position = gram.Position()
            position.height = h * 1e-3
            position.latitude = rad2deg(lat)
            position.longitude = rad2deg(lon)

            position.elapsedTime = el_time # Time since start in s
            atmosphere.setPosition(position)
            atmosphere.update()
            atmos = atmosphere.getAtmosphereState()
            rho = pyconvert(Float64, atmos.density)
            T = pyconvert(Float64, atmos.temperature)
            wind = SVector{3, Float64}([pyconvert(Float64, montecarlo ? atmos.perturbedEWWind : atmos.ewWind),
                    pyconvert(Float64, montecarlo ? atmos.perturbedNSWind : atmos.nsWind),
                    pyconvert(Float64, atmos.verticalWind)])
        end
    end

    return rho, T, wind
end

function density_nrlmsise(h::Float64, p, lat::Float64, lon::Float64, montecarlo::Bool, Wind::Bool, args::Dict, current_time::DateTime)
    """

    """

    if config.cnf().drag_state == false && args[:keplerian] == false
        rho , T , wind = density_exp(h, p)
        rho = 0.0
    elseif config.cnf().drag_state == true || args[:keplerian] == true
        jd = datetime2julian(current_time)
        atmo = SatelliteToolbox.AtmosphericModels.nrlmsise00(jd, h, lat, lon, 150, 150, 3)
        rho = atmo.total_density
        T = atmo.temperature
        wind = [0.0,0.0,0.0]
    end

    return rho, T, wind
end