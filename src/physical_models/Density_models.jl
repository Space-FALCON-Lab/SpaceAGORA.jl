# include("../utils/Reference_system.jl")

# using .SimulationModel
using SatelliteToolbox
using StaticArrays
using LinearAlgebra
using ..SimulationModel: GRAM_LOCK
# export NoAtmosphereModel, ExponentialAtmosphereModel, GRAMAtmosphereModel, PolynomialFitAtmosphereModel, NRLMSISE00AtmosphereModel
# export getDensity
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
    gram::G
    gram_atmosphere::GA
    gram_root::String
    gram_data_root::String
    spice_root::String
    planet_name::String
    initial_time::InitialTime
end

struct PolynomialFitAtmosphereModel <: AbstractDensityModel
    polyfit_coeffs::Vector{Float64} # Coefficients for the polynomial fit (ordered from highest degree to lowest)
end

struct NRLMSISE00AtmosphereModel <: AbstractDensityModel
    # No additional fields needed for NRLMSISE-00, but we can add any necessary parameters or configurations here if needed in the future
end

const _GRAM_WRAPPER = Ref{Any}(nothing)
const _GRAM_WRAPPER_FILE = Ref{String}("")
const _GRAM_SEED_WARNING_EMITTED = Ref(false)
const _GRAM_WIND_WARNING_EMITTED = Ref(false)
const _GRAM_LOCK_OFF_WARNING_EMITTED = Ref(false)

@inline _spaceagora_repo_root() = normpath(joinpath(@__DIR__, "..", ".."))
@inline _gram_lib_extension() = Sys.iswindows() ? "dll" : (Sys.isapple() ? "dylib" : "so")
@inline function _gram_use_global_lock()::Bool
    mode = lowercase(strip(get(ENV, "SPACEAGORA_GRAM_GLOBAL_LOCK", "on")))
    if mode in ("on", "true", "1", "yes")
        return true
    elseif mode in ("off", "false", "0", "no")
        if Threads.nthreads() > 1
            if !_GRAM_LOCK_OFF_WARNING_EMITTED[]
                _GRAM_LOCK_OFF_WARNING_EMITTED[] = true
                @warn "SPACEAGORA_GRAM_GLOBAL_LOCK=off is unsafe with Threads.nthreads()>1. Keeping GRAM global lock enabled."
            end
            return true
        end
        return false
    end
    throw(ArgumentError("Unsupported SPACEAGORA_GRAM_GLOBAL_LOCK='$mode'. Use one of: on, off."))
end

function _as_abspath(path::AbstractString)::String
    if isempty(path)
        return ""
    end
    return isabspath(path) ? normpath(path) : normpath(joinpath(_spaceagora_repo_root(), path))
end

@inline function _is_gram_root(path::AbstractString)::Bool
    return isdir(joinpath(path, "Build")) && isdir(joinpath(path, "Julia"))
end

function _resolve_gram_root(gram_root_directory::String, gram_directory::String)::String
    candidates = String[]

    !isempty(gram_root_directory) && push!(candidates, gram_root_directory)
    !isempty(get(ENV, "GRAM_ROOT", "")) && push!(candidates, ENV["GRAM_ROOT"])

    if !isempty(gram_directory)
        push!(candidates, gram_directory)
        push!(candidates, joinpath(gram_directory, ".."))
        push!(candidates, joinpath(gram_directory, "..", ".."))
    end

    push!(candidates, "GRAM Suite 2.0")

    for candidate in candidates
        candidate_abs = _as_abspath(candidate)
        _is_gram_root(candidate_abs) && return candidate_abs
    end

    throw(ArgumentError(
        "Unable to locate GRAM Suite root. Set `gram_root_directory` or `ENV[\"GRAM_ROOT\"]` to a path containing `Build/` and `Julia/`."
    ))
end

function _load_gram_wrapper!(gram_root::String)::Tuple{Any, Bool}
    wrapper_file = normpath(joinpath(gram_root, "Julia", "GRAM.jl"))
    isfile(wrapper_file) || throw(ArgumentError("GRAM Julia wrapper not found at: $wrapper_file"))

    loaded_now = false
    if _GRAM_WRAPPER[] === nothing
        Base.include(@__MODULE__, wrapper_file)
        _GRAM_WRAPPER[] = Base.invokelatest(getfield, @__MODULE__, :GRAM)
        _GRAM_WRAPPER_FILE[] = wrapper_file
        loaded_now = true
    elseif _GRAM_WRAPPER_FILE[] != wrapper_file
        throw(ArgumentError(
            "GRAM wrapper already loaded from $(_GRAM_WRAPPER_FILE[]). Restart Julia to load a different wrapper path."
        ))
    end

    return _GRAM_WRAPPER[], loaded_now
end

function _resolve_spice_directory(gram_root::String, spice_directory::String, gram_data_directory::String)::String
    candidates = String[]
    !isempty(spice_directory) && push!(candidates, spice_directory)
    push!(candidates, joinpath(gram_root, "SPICE"))
    !isempty(gram_data_directory) && push!(candidates, joinpath(gram_data_directory, "SPICE"))
    push!(candidates, "GRAM_Data/SPICE")

    for candidate in candidates
        candidate_abs = _as_abspath(candidate)
        isdir(candidate_abs) && return candidate_abs
    end

    throw(ArgumentError("Unable to find a valid SPICE directory for GRAM initialization."))
end

function _resolve_gram_data_root(gram_root::String, gram_data_directory::String)::String
    candidates = String[]
    push!(candidates, gram_root)
    !isempty(gram_data_directory) && push!(candidates, gram_data_directory)
    push!(candidates, "GRAM_Data")

    for candidate in candidates
        candidate_abs = _as_abspath(candidate)
        if isdir(joinpath(candidate_abs, "Earth", "data")) || isdir(joinpath(candidate_abs, "Mars", "data"))
            return candidate_abs
        end
    end

    throw(ArgumentError("Unable to find GRAM data directories (expected `<root>/Earth/data` or `<root>/Mars/data`)."))
end

function _gram_body(gram, planet_name::String)::Int
    key = lowercase(strip(planet_name))
    key == "earth" && return gram.BODY_EARTH
    key == "mars" && return gram.BODY_MARS
    key == "venus" && return gram.BODY_VENUS
    key == "titan" && return gram.BODY_TITAN
    key == "jupiter" && return gram.BODY_JUPITER
    key == "uranus" && return gram.BODY_URANUS
    key == "neptune" && return gram.BODY_NEPTUNE
    throw(ArgumentError("Unsupported GRAM planet_name '$planet_name'."))
end

function _gram_data_path(gram_data_root::String, planet_name::String)::Union{Nothing, String}
    key = lowercase(strip(planet_name))
    if key == "earth" || key == "mars"
        planet_dir = uppercasefirst(key)
        data_path = joinpath(gram_data_root, planet_dir, "data")
        isdir(data_path) || throw(ArgumentError("GRAM data path not found: $data_path"))
        return data_path
    end
    return nothing
end

function GRAMAtmosphereModel(;
    gram_directory::String="GRAM Suite 2.0",
    gram_data_directory::String="GRAM_Data",
    gram_root_directory::String="",
    gram_library_path::String="",
    spice_directory::String="",
    planet_name::String="earth",
    seed::Int=1001,
    initial_time::InitialTime=InitialTime()
)
    gram_root = _resolve_gram_root(gram_root_directory, gram_directory)
    gram, loaded_now = _load_gram_wrapper!(gram_root)
    if loaded_now
        return Base.invokelatest(
            GRAMAtmosphereModel;
            gram_directory=gram_directory,
            gram_data_directory=gram_data_directory,
            gram_root_directory=gram_root_directory,
            gram_library_path=gram_library_path,
            spice_directory=spice_directory,
            planet_name=planet_name,
            seed=seed,
            initial_time=initial_time
        )
    end

    if !isempty(gram_library_path)
        gram.set_library!(_as_abspath(gram_library_path))
    else
        local_lib = joinpath(gram_root, "Build", "lib", "libGRAM.$(_gram_lib_extension())")
        if isfile(local_lib)
            gram.set_library!(local_lib)
        end
    end

    spice_root = _resolve_spice_directory(gram_root, spice_directory, gram_data_directory)
    gram.initialize!(spice_root)

    gram_data_root = _resolve_gram_data_root(gram_root, gram_data_directory)
    body = _gram_body(gram, planet_name)
    data_path = _gram_data_path(gram_data_root, planet_name)
    gram_atmosphere = data_path === nothing ? gram.create_atmosphere(body) : gram.create_atmosphere(body; data_path=data_path)

    gram.set_start_time!(
        gram_atmosphere;
        year=Int(initial_time.year),
        month=Int(initial_time.month),
        day=Int(initial_time.day),
        hour=Int(initial_time.hour),
        minute=Int(initial_time.minute),
        seconds=Float64(initial_time.second),
        scale=1,
        frame=1
    )

    if seed != 1001 && !_GRAM_SEED_WARNING_EMITTED[]
        _GRAM_SEED_WARNING_EMITTED[] = true
        @warn "GRAMAtmosphereModel seed is ignored by the current Julia GRAM wrapper."
    end

    return GRAMAtmosphereModel(
        gram,
        gram_atmosphere,
        gram_root,
        gram_data_root,
        spice_root,
        lowercase(strip(planet_name)),
        initial_time
    )
end

function Base.deepcopy_internal(model::GRAMAtmosphereModel, stackdict::IdDict)
    if haskey(stackdict, model)
        return stackdict[model]
    end

    copied = GRAMAtmosphereModel(
        gram_root_directory=model.gram_root,
        gram_data_directory=model.gram_data_root,
        spice_directory=model.spice_root,
        planet_name=model.planet_name,
        initial_time=model.initial_time
    )
    stackdict[model] = copied
    return copied
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

    return p.T_ref
end

function getDensity(model::NoAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    # Return zero density and a default temperature for the no-atmosphere model
    ρ = 0.0
    T = p.args.environment_model.planet.T_ref
    wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
    return ρ, T, wind_vec
end

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

function getDensity(model::PolynomialFitAtmosphereModel, h::Union{Float64, Vector{Float64}}, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params

    if typeof(h) != Float64 || length(h) != 1
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
        T = p.args.environment_model.planet.T_ref

        wind_vec = SVector{3, Float64}(0, 0, 0)

        return ρ, T, wind_vec
    else
        # polyfit = model.polyfit_coeffs
        power = MVector{length(model.polyfit_coeffs)}(zeros(length(model.polyfit_coeffs)))
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
        T = p.args.environment_model.planet.T_ref

        wind_vec = SVector{3, Float64}(0, 0, 0)

        return ρ, T, wind_vec
    end
end

@inline function _gram_density_state(
    model::GRAMAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool
)::Tuple{Float64, Float64, SVector{3, Float64}}
    model.gram.set_position!(
        model.gram_atmosphere;
        height=h * 1e-3,
        latitude=rad2deg(lat),
        longitude=rad2deg(lon),
        elapsed_time=el_time
    )

    err = model.gram.update!(model.gram_atmosphere)
    err == 0 || throw(ErrorException("GRAM update failed (code=$err): $(model.gram.get_error_message())"))

    atmos = model.gram.get_dynamics_state(model.gram_atmosphere)
    rho_local = Float64(atmos.density)
    T_local = Float64(atmos.temperature)
    if isdefined(model.gram, :get_winds_state)
        winds = model.gram.get_winds_state(model.gram_atmosphere)
        wind_vec_local = if wind
            SVector{3, Float64}(
                Float64(winds.perturbedEWWind),
                Float64(winds.perturbedNSWind),
                Float64(winds.perturbedVerticalWind)
            )
        else
            SVector{3, Float64}(
                Float64(winds.ewWind),
                Float64(winds.nsWind),
                Float64(winds.verticalWind)
            )
        end
        return rho_local, T_local, wind_vec_local
    end

    if !_GRAM_WIND_WARNING_EMITTED[]
        _GRAM_WIND_WARNING_EMITTED[] = true
        @warn "Loaded GRAM wrapper does not expose wind state. Returning zero wind."
    end
    return rho_local, T_local, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function getDensity(model::GRAMAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    # cnf = params.cnf
    # println("In GRAM density model")
    # Check atmospheric interface conditions
    # planet_radius = p.args.environment_model.planet.Rp_e
    EI = p.args.environment_model.EI * 1e3
    drag_state = h - EI <= 0.0
    if !drag_state && !p.args.mission_configuration.keplerian
        if h > 2000.0e3
            rho = 0.0
            T = p.args.environment_model.planet.T_ref
            wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
        else
            rho, T, wind_vec = density_polyfit(h, p)
        end
    else
        if h > 2000.0e3
            rho = 0.0
            T = p.args.environment_model.planet.T_ref
            wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
        else
            rho, T, wind_vec = if _gram_use_global_lock()
                lock(GRAM_LOCK) do
                    _gram_density_state(model, h, lat, lon, el_time, wind)
                end
            else
                _gram_density_state(model, h, lat, lon, el_time, wind)
            end
        end
    end
    # println("el_time: ", el_time, "h: ", h, " lat: ", rad2deg(lat), " lon: ", rad2deg(lon), " rho: ", rho, " T: ", T, " wind_vec: ", wind_vec)
    # println("rho: ", rho)
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

function density_polyfit(h::Float64, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    # Load the polynomial coefficients for the planet
    coeffs = p.args.environment_model.planet.polyfit_coeffs
    poly_model = PolynomialFitAtmosphereModel(vec(coeffs))
    ρ, T, wind_vec = getDensity(poly_model, h, 0.0, 0.0, 0.0, false, p) # Latitude, longitude, elapsed time, and wind don't affect the density in this model
    return ρ, T, wind_vec
end
# end # module
