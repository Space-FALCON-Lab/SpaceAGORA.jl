struct SpaceAGORADensityAdapter{M,C,P} <: Function
    model::M
    context::C
    planet::P
    latitude::Float64
    longitude::Float64
    wind::Bool
end

function (source::SpaceAGORADensityAdapter)(altitude_m, elapsed_time_s,
    position_ii_m=nothing)
    if position_ii_m === nothing
        sampled_altitude_m = Float64(altitude_m)
        sampled_latitude = source.latitude
        sampled_longitude = source.longitude
    else
        sampled_altitude_m, sampled_latitude, sampled_longitude =
            rtolatlong(SVector{3,Float64}(position_ii_m), source.planet)
    end
    return Float64(getDensity(
        source.model,
        Float64(sampled_altitude_m),
        Float64(sampled_latitude),
        Float64(sampled_longitude),
        Float64(elapsed_time_s),
        source.wind,
        source.context,
    )[1])
end

function (source::SpaceAGORADensityAdapter)(::Val{:gradient}, position_ii_m,
    elapsed_time_s)
    position = SVector{3,Float64}(position_ii_m)
    altitude_m, latitude, longitude = rtolatlong(position, source.planet)
    model = source.model
    if applicable(density_altitude_derivative, model, altitude_m)
        derivative_per_m = density_altitude_derivative(model, altitude_m)
        up = @SVector [
            cos(latitude) * cos(longitude),
            cos(latitude) * sin(longitude),
            sin(latitude),
        ]
        return derivative_per_m * up
    end
    step_m = 1.0
    return SVector{3,Float64}(ntuple(axis -> begin
        offset = SVector{3,Float64}(ntuple(
            index -> index == axis ? step_m : 0.0, 3))
        (source(0.0, elapsed_time_s, position + offset) -
         source(0.0, elapsed_time_s, position - offset)) / (2.0 * step_m)
    end, 3))
end

"""
Return the configured SpaceAGORA atmosphere as an MPC density source.

The returned callable evaluates `getDensity` directly. It also supplies the
Cartesian density gradient required by the analytical drag Jacobian. For a
`PolynomialFitAtmosphereModel`, that gradient is the exact derivative of the
configured log-density polynomial.
"""
function density_function_from_spaceagora(
    obj;
    latitude::Real=0.0,
    longitude::Real=0.0,
    wind::Bool=false,
)
    args = _mpc_scenario_args(obj)
    context = hasproperty(obj, :args) ? obj : args
    return SpaceAGORADensityAdapter(
        args.environment_model.density_model,
        context,
        args.environment_model.planet,
        Float64(latitude),
        Float64(longitude),
        wind,
    )
end
