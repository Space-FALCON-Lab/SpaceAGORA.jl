mutable struct MarsAtmosphereAdapter
    planet
    actual_mode::Symbol
    gram_model
    fallback_model
    cache::Dict{NTuple{4, Int}, Tuple{Float64, Float64, SVector{3, Float64}}}
    cache_altitude_step_m::Float64
    cache_longitude_step_deg::Float64
    cache_time_step_s::Float64
    lat_deg::Float64
    initial_elapsed_time_s::Float64
end

struct OfflineGRAMSurrogateTag
    planet_name::String
end

@inline function _offline_gram_surrogate_file(config::PrelimConfig)
    candidate = joinpath(
        config.repo_root,
        "data",
        "GRAMSuite.jl",
        "GRAM Suite 2.0",
        "simulation",
        "GRAM",
        "static_grids",
        "mars_surrogate.jls",
    )
    return isfile(candidate) ? candidate : ""
end

function build_atmosphere_adapter(config::PrelimConfig)
    fallback_model = SpaceAGORA.ExponentialAtmosphereModel(config.planet; temperature_k=config.planet.T_ref)
    gram_model = nothing
    actual_mode = :exponential
    if config.atmosphere_mode == :gram
        try
            gram_model = if config.use_gram_surrogate
                try
                    ENV.GRAMAtmosphereModelSurrogate(planet_name="mars", initial_time=config.initial_time)
                catch
                    GRAMSuite.GRAMAtmosphereModel(planet_name="mars", initial_time=config.initial_time)
                end
            else
                try
                    ENV.GRAMAtmosphereModel(planet_name="mars", initial_time=config.initial_time)
                catch
                    GRAMSuite.GRAMAtmosphereModel(planet_name="mars", initial_time=config.initial_time)
                end
            end
            actual_mode = :gram
        catch err
            surrogate_file = _offline_gram_surrogate_file(config)
            if !isempty(surrogate_file)
                @warn "Compiled GRAM initialization failed; using bundled offline Mars surrogate grid for this study." exception=(err, catch_backtrace()) surrogate_file
                gram_model = (
                    base=OfflineGRAMSurrogateTag("mars"),
                    surrogate_file=surrogate_file,
                )
                actual_mode = :gram_offline
            else
                @warn "Falling back to exponential Mars atmosphere because GRAM initialization failed." exception=(err, catch_backtrace())
            end
        end
    end
    return MarsAtmosphereAdapter(
        config.planet,
        actual_mode,
        gram_model,
        fallback_model,
        Dict{NTuple{4, Int}, Tuple{Float64, Float64, SVector{3, Float64}}}(),
        config.gram_cache_altitude_m,
        config.gram_cache_longitude_deg,
        config.gram_cache_time_s,
        config.lat_deg,
        config.initial_elapsed_time_s,
    )
end

@inline function atmosphere_label(adapter::MarsAtmosphereAdapter)
    return adapter.actual_mode in (:gram, :gram_offline) ? "marsgram" : "exponential_fallback"
end

@inline function _wrapped_longitude_deg(theta_rad::Float64)
    return mod(rad2deg(theta_rad), 360.0)
end

@inline function _atmosphere_cache_key(
    adapter::MarsAtmosphereAdapter,
    h_m::Float64,
    lat_deg::Float64,
    lon_deg::Float64,
    elapsed_time_s::Float64,
)
    return (
        round(Int, h_m / adapter.cache_altitude_step_m),
        round(Int, lat_deg / adapter.cache_longitude_step_deg),
        round(Int, lon_deg / adapter.cache_longitude_step_deg),
        round(Int, elapsed_time_s / adapter.cache_time_step_s),
    )
end

function atmosphere_state(
    adapter::MarsAtmosphereAdapter,
    h_m::Float64,
    lat_rad::Float64,
    lon_rad::Float64,
    elapsed_time_s::Float64;
    include_wind::Bool=false,
)::Tuple{Float64, Float64, SVector{3, Float64}}
    h_clamped = max(h_m, 0.0)
    lat_deg = rad2deg(lat_rad)
    lon_deg = mod(rad2deg(lon_rad), 360.0)
    key = _atmosphere_cache_key(adapter, h_clamped, lat_deg, lon_deg, elapsed_time_s)
    cached = get(adapter.cache, key, nothing)
    if cached !== nothing && (!include_wind || cached[3] != SVector{3, Float64}(0.0, 0.0, 0.0))
        return cached
    end

    if adapter.actual_mode == :gram && adapter.gram_model !== nothing
        try
            rho, temperature_k, wind_vec = if adapter.gram_model isa ENV.GRAMAtmosphereModel || adapter.gram_model isa ENV.GRAMAtmosphereModelSurrogate
                ENV._gram_point_density(
                    adapter.gram_model,
                    h_clamped,
                    lat_rad,
                    deg2rad(lon_deg),
                    adapter.initial_elapsed_time_s + elapsed_time_s,
                    include_wind,
                )
            else
                GRAMSuite.point_density_state(
                    adapter.gram_model,
                    h_clamped,
                    lat_rad,
                    deg2rad(lon_deg),
                    adapter.initial_elapsed_time_s + elapsed_time_s,
                    include_wind,
                )
            end
            rho = max(Float64(rho), 0.0)
            temperature_k = max(Float64(temperature_k), 1.0)
            wind_state = SVector{3, Float64}(Float64(wind_vec[1]), Float64(wind_vec[2]), Float64(wind_vec[3]))
            adapter.cache[key] = (rho, temperature_k, wind_state)
            return rho, temperature_k, wind_state
        catch err
            @warn "GRAM point-density query failed; using exponential fallback for this sample." exception=(err, catch_backtrace())
            adapter.actual_mode = :exponential
        end
    elseif adapter.actual_mode == :gram_offline && adapter.gram_model !== nothing
        try
            rho, temperature_k, wind_vec = GRAMSuite.surrogate_density_state(
                adapter.gram_model.base,
                adapter.gram_model.surrogate_file,
                nothing,
                h_clamped,
                lat_rad,
                deg2rad(lon_deg),
                adapter.initial_elapsed_time_s + elapsed_time_s,
                include_wind;
                vacuum_temperature=adapter.planet.T_ref,
            )
            rho = max(Float64(rho), 0.0)
            temperature_k = max(Float64(temperature_k), 1.0)
            wind_state = SVector{3, Float64}(Float64(wind_vec[1]), Float64(wind_vec[2]), Float64(wind_vec[3]))
            adapter.cache[key] = (rho, temperature_k, wind_state)
            return rho, temperature_k, wind_state
        catch err
            @warn "Offline GRAM surrogate query failed; using exponential fallback for this sample." exception=(err, catch_backtrace())
            adapter.actual_mode = :exponential
        end
    end

    rho, temperature_k, wind_vec = SpaceAGORA.getDensity(
        adapter.fallback_model,
        h_clamped,
        lat_deg,
        lon_deg,
        adapter.initial_elapsed_time_s + elapsed_time_s,
        include_wind,
    )
    rho = max(Float64(rho), 0.0)
    temperature_k = max(Float64(temperature_k), 1.0)
    wind_state = SVector{3, Float64}(Float64(wind_vec[1]), Float64(wind_vec[2]), Float64(wind_vec[3]))
    adapter.cache[key] = (rho, temperature_k, wind_state)
    return rho, temperature_k, wind_state
end

function density_temperature(
    adapter::MarsAtmosphereAdapter,
    h_m::Float64,
    theta_rad::Float64,
    elapsed_time_s::Float64,
)::Tuple{Float64, Float64}
    rho, temperature_k, _ = atmosphere_state(
        adapter,
        h_m,
        deg2rad(adapter.lat_deg),
        theta_rad,
        elapsed_time_s;
        include_wind=false,
    )
    return rho, temperature_k
end

@inline function mach_number(adapter::MarsAtmosphereAdapter, velocity_mps::Float64, temperature_k::Float64)::Float64
    speed_of_sound = sqrt(max(adapter.planet.γ * adapter.planet.R * max(temperature_k, 1.0), 1.0))
    return velocity_mps / speed_of_sound
end
