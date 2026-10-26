struct SpaceAGORADensityAdapter{M,C,P} <: Function
    model::M
    context::C
    planet::P
    latitude::Float64
    longitude::Float64
    wind::Bool
    epoch_s::Float64
    rotation_rate_rad_s::Float64
end

# The KS prediction runs in a non-rotating frame aligned with the planet-fixed
# frame at `epoch_s`. After `elapsed_time_s` the planet has turned by Ωτ about
# its pole, so planet-fixed coordinates are the prediction ones turned by -Ωτ.
@inline function _mpc_prediction_to_planet_fixed(
    source::SpaceAGORADensityAdapter, position, elapsed_time_s)
    θ = source.rotation_rate_rad_s * Float64(elapsed_time_s)
    c, s = cos(θ), sin(θ)
    x, y, z = position
    return SVector{3,Float64}(c * x + s * y, -s * x + c * y, z)
end

@inline function _mpc_planet_fixed_to_prediction(
    source::SpaceAGORADensityAdapter, vector, elapsed_time_s)
    θ = source.rotation_rate_rad_s * Float64(elapsed_time_s)
    c, s = cos(θ), sin(θ)
    x, y, z = vector
    return SVector{3,Float64}(c * x - s * y, s * x + c * y, z)
end

function (source::SpaceAGORADensityAdapter)(altitude_m, elapsed_time_s,
    position_m=nothing)
    if position_m === nothing
        sampled_altitude_m = Float64(altitude_m)
        sampled_latitude = source.latitude
        sampled_longitude = source.longitude
    else
        position_pf = _mpc_prediction_to_planet_fixed(
            source, SVector{3,Float64}(position_m), elapsed_time_s)
        sampled_altitude_m, sampled_latitude, sampled_longitude =
            rtolatlong(position_pf, source.planet)
    end
    return Float64(getDensity(
        source.model,
        Float64(sampled_altitude_m),
        Float64(sampled_latitude),
        Float64(sampled_longitude),
        source.epoch_s + Float64(elapsed_time_s),
        source.wind,
        source.context,
    )[1])
end

function (source::SpaceAGORADensityAdapter)(::Val{:gradient}, position_m,
    elapsed_time_s)
    position = SVector{3,Float64}(position_m)
    model = source.model
    if applicable(density_altitude_derivative, model, 0.0)
        position_pf = _mpc_prediction_to_planet_fixed(source, position, elapsed_time_s)
        altitude_m, latitude, longitude = rtolatlong(position_pf, source.planet)
        derivative_per_m = density_altitude_derivative(model, altitude_m)
        # Geodetic up direction in the planet-fixed frame, back to prediction axes.
        up_pf = @SVector [
            cos(latitude) * cos(longitude),
            cos(latitude) * sin(longitude),
            sin(latitude),
        ]
        return derivative_per_m *
            _mpc_planet_fixed_to_prediction(source, up_pf, elapsed_time_s)
    end
    step_m = 1.0
    center = source(0.0, elapsed_time_s, position)
    return SVector{3,Float64}(ntuple(axis -> begin
        offset = SVector{3,Float64}(ntuple(
            index -> index == axis ? step_m : 0.0, 3))
        ks_density_difference(
            source(0.0, elapsed_time_s, position + offset),
            center,
            source(0.0, elapsed_time_s, position - offset),
            step_m)
    end, 3))
end

"""
    density_function_from_spaceagora(args_or_params; latitude=0, longitude=0,
        wind=false, epoch_s=0) -> callable

Return the configured SpaceAGORA atmosphere as an MPC density source
`(altitude_m, elapsed_s[, position_m])`.

With a position, the callable samples `getDensity` at the geodetic altitude,
latitude and longitude of that position. Positions are expressed in the MPC
prediction frame: non-rotating, aligned with the planet-fixed frame at
`epoch_s` (s into the run), so after `elapsed_s` the planet has turned by
`Ω * elapsed_s` about its pole. The atmosphere is sampled at run time
`epoch_s + elapsed_s`. Without a position it uses `altitude_m` and the given
`latitude`/`longitude` (rad).

The callable also supplies the Cartesian density gradient the drag Jacobian
needs: analytic along the geodetic vertical for models with a
[`density_altitude_derivative`](@ref) method (the polynomial fit gives the exact
derivative of its log-density polynomial), otherwise a 1 m finite difference.
"""
function density_function_from_spaceagora(
    obj;
    latitude::Real=0.0,
    longitude::Real=0.0,
    wind::Bool=false,
    epoch_s::Real=0.0,
)
    args = _mpc_scenario_args(obj)
    context = hasproperty(obj, :args) ? obj : args
    planet = args.environment_model.planet
    return SpaceAGORADensityAdapter(
        args.environment_model.density_model,
        context,
        planet,
        Float64(latitude),
        Float64(longitude),
        wind,
        Float64(epoch_s),
        _planet_rotation_rate_rad_s(planet),
    )
end
