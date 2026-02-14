# include("../utils/Reference_system.jl")

# using .SimulationModel
using PythonCall
using SatelliteToolbox
# using .AbstractTypes: AbstractPlanet, AbstractDensityModel

struct NoAtmosphereModel <: AbstractDensityModel
    # No additional fields needed for a simple no-atmosphere model
end

struct ExponentialAtmosphereModel <: AbstractDensityModel
    ρ_ref::Float64 # Reference density at reference height (kg/m³)
    h_ref::Float64 # Reference height (m)
    H::Float64     # Scale height (m)
end

struct GRAMAtmosphereModel{G, GA} <: AbstractDensityModel
    # GRAM Version information (e.g., "GRAM2016", "GRAM2010", etc.)
    gram::G
    gram_atmosphere::GA
end

struct PolynomialFitAtmosphereModel <: AbstractDensityModel
    polyfit_coeffs::Vector{Float64} # Coefficients for the polynomial fit (ordered from highest degree to lowest)
end

struct NRLMSISE00AtmosphereModel <: AbstractDensityModel
    # No additional fields needed for NRLMSISE-00, but we can add any necessary parameters or configurations here if needed in the future
end

function GRAMAtmosphereModel(; gram_directory::String="GRAMpy", gram_data_directory::String="GRAM_Data", planet_name::String="earth", seed::Int=1001, initial_time::InitialTime=InitialTime())
    # Initialize the GRAM model and atmosphere object
    sys = pyimport("sys")
    os = pyimport("os")
    if !(gram_directory in pyconvert(Vector{String}, sys.path))
        sys.path.append(gram_directory)
    end
    gram = pyimport("gram")
    inputParameters = Dict("earth" => gram.EarthInputParameters(),
                            "mars" => gram.MarsInputParameters(),
                            "venus" => gram.VenusInputParameters(),
                            "titan" => gram.TitanInputParameters())
        
    namelistReaders = Dict("earth" => gram.EarthNamelistReader(),
                        "mars" => gram.MarsNamelistReader(),
                        "venus" => gram.VenusNamelistReader(),
                        "titan" => gram.TitanNamelistReader())
        
    atmospheres = Dict("earth" => gram.EarthAtmosphere(),
                    "mars" => gram.MarsAtmosphere(),
                    "venus" => gram.VenusAtmosphere(),
                    "titan" => gram.TitanAtmosphere())

    input_parameters = inputParameters[planet_name]

    # Mars has some weird specific parameters, so this line is just to check to make sure the it doesn't do it for the other planets
    if planet_name == "mars"
        # input_parameters.dataPath = os.path.join(os.path.dirname(os.path.abspath(@__FILE__)),"..", "GRAM_Data", "Mars", "data", "")
        input_parameters.dataPath = gram_data_directory * "/Mars/data/"
        if !Bool(os.path.exists(input_parameters.dataPath))
            throw(ArgumentError("GRAM data path not found: " * input_parameters.dataPath))
        end
    end

    if planet_name == "earth"
        # input_parameters.dataPath = os.path.join(os.path.dirname(os.path.abspath(@__FILE__)),"..", "GRAM_Data", "Mars", "data", "")
        input_parameters.dataPath = gram_data_directory * "/Earth/data/"
        if !Bool(os.path.exists(input_parameters.dataPath))
            throw(ArgumentError("GRAM data path not found: " * input_parameters.dataPath))
        end
    end

    reader = namelistReaders[planet_name]
    reader.tryGetSpicePath(input_parameters)

    gram_atmosphere = atmospheres[planet_name]
    gram_atmosphere.setInputParameters(input_parameters)
    
    if planet_name == "earth"
        gram_atmosphere.setMERRA2Parameters(0, -90.0, 90.0, 0.0, 359.99999)
    end
    gram_atmosphere.setPerturbationScales(1.5)
    gram_atmosphere.setMinRelativeStepSize(0.5)
    gram_atmosphere.setSeed(seed)

    if planet_name == "mars"
        gram_atmosphere.setMOLAHeights(false)
    end

    ttime = gram.GramTime()
    ttime.setStartTime(initial_time.year, initial_time.month, initial_time.day, initial_time.hour, initial_time.minute, initial_time.second, gram.UTC, gram.PET)
    gram_atmosphere.setStartTime(ttime)
    return GRAMAtmosphereModel(gram, gram_atmosphere)
end
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
    # if config.cnf.drag_state == true
    #     T = p.T
    # else
    #     T = p.T
    # end

    return p.T
end

# function density_constant(h, p, OE=0, lat=0, lon=0, timereal=0, t0=0, tf_prev=0, montecarlo=0, Wind=0, args=0, version=[])
#     """

#     """

#     if config.cnf.drag_state == false
#         ρ = 0.0
#     else
#         ρ = p.ρ_ref
#     end
    
#     T = temperature_linear(h, p)

#     wind = [0, 0, 0]

#     return ρ, T, wind
# end

# function density_no(h, p, OE=0, lat=0, lon=0, timereal=0, t0=0, tf_prev=0, montecarlo=0, Wind=0, args=0, version=[])
#     """

#     """

#     T = temperature_linear(h, p)

#     wind = [0, 0, 0]

#     return 0.0, T, wind
# end

"""
    Get the density using the Exponential atmosphere model.

    Parameters
    ----------
    model : ExponentialAtmosphereModel
        Model of the atmosphere.
    h : Float64
        Height in meters.
    lat : Float64
        Latitude in radians.
    lon : Float64
        Longitude in radians.
    el_time : Float64
        Elapsed time since the start of the simulation in seconds.
    wind : Bool
        Whether to include wind effects in the density calculation.

    Returns
    -------
    ρ : Float64
        Density in kg/m³.
    T : Float64
        Temperature in Kelvin.
    wind : SVector{3, Float64}
        Wind vector in m/s.
"""
function getDensity(model::ExponentialAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool)::Tuple{Float64, Float64, SVector{3, Float64}}
    """

    """

    ρ = p.ρ_ref * exp.((p.h_ref .- h)/p.H)

    T = temperature_linear(h, p)

    wind_vec = SVector{3, Float64}(0, 0, 0)

    return ρ, T, wind_vec
end

function getDensity(model::PolynomialFitAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool)::Tuple{Float64, Float64, SVector{3, Float64}}

    if typeof(h) != Float64
        # polyfit = model.polyfit_coeffs
        power = zeros(length(model.polyfit_coeffs),length(h))
        # Convert height from meters to kilometers
        h *= 1e-3
        # Calculate the polynomial value at height h
        for i=1:length(polyfit)
            power[i,:] = (h).^(length(polyfit)-i)
        end

        # Calculate the exponent term of the density using the polynomial coefficients
        exponent = zeros(length(h))
        @inbounds for j=eachindex(h)
            exponent[j] = sum(polyfit .* power[:,j])
        end

        # Calculate the density
        ρ = exp.(exponent)
        T = temperature_linear(h, p)

        wind_vec = SVector{3, Float64}(0, 0, 0)

        return ρ, T, wind_vec
    else
        # polyfit = model.polyfit_coeffs
        power = MVector(zeros(length(model.polyfit_coeffs)))
        # Convert height from meters to kilometers
        h *= 1e-3
        # Calculate the polynomial value at height h
        @inbounds for i=eachindex(model.polyfit_coeffs)
            power[i] = h^(length(model.polyfit_coeffs)-i)
        end
        # Calculate the exponent term of the density using the polynomial coefficients
        exponent = sum(model.polyfit_coeffs .* power)
        # Calculate the density
        ρ = exp(exponent)
        T = temperature_linear(h, p)

        wind_vec = SVector{3, Float64}(0, 0, 0)

        return ρ, T, wind_vec
    end
end

function getDensity(model::GRAMAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool)::Tuple{Float64, Float64, SVector{3, Float64}}
    cnf = params.cnf
    if cnf.drag_state == false && args[:keplerian] == false
        if h > 2000.0e3
            rho = 0.0
            T = p.T
            wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
        else
            rho, T, wind_vec = density_polyfit(h, p)
        end
    else
        if h > 2000.0e3
            rho = 0.0
            T = temperature_linear(h, p)
            wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
        else
            position = model.gram.Position()
            position.height = h * 1e-3
            position.latitude = rad2deg(lat)
            position.longitude = rad2deg(lon)

            position.elapsedTime = el_time # Time since start in s
            model.gram_atmosphere.setPosition(position)
            model.gram_atmosphere.update()
            atmos = model.gram_atmosphere.getAtmosphereState()
            rho = pyconvert(Float64, atmos.density)
            T = pyconvert(Float64, atmos.temperature)
            wind_vec = SVector{3, Float64}([pyconvert(Float64, wind ? atmos.perturbedEWWind : atmos.ewWind),
                    pyconvert(Float64, wind ? atmos.perturbedNSWind : atmos.nsWind),
                    pyconvert(Float64, atmos.verticalWind)])
        end
    end

    return rho, T, wind_vec
end

function getDensity(model::NRLMSISE00AtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool)::Tuple{Float64, Float64, SVector{3, Float64}}

    if config.cnf.drag_state == false && args[:keplerian] == false
        rho , T , wind_vec = density_exp(h, p)
        rho = 0.0
    elseif config.cnf.drag_state == true || args[:keplerian] == true
        jd = datetime2julian(current_time)
        atmo = SatelliteToolbox.AtmosphericModels.nrlmsise00(jd, h, lat, lon, 150, 150, 3)
        rho = atmo.total_density
        T = atmo.temperature
        wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    return rho, T, wind_vec
end