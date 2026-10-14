struct GramSample
    mean_density::Float64
    std_density::Float64
end

struct GramSeedPlan
    truth_seed::Int
    prior_nominal_seed::Int
    prior_dispersion_seeds::Vector{Int}
end

function gram_seed_plan(truth_seed::Int, prior_seed::Int, n_dispersion::Int)::GramSeedPlan
    n_dispersion >= 2 || throw(ArgumentError("At least two dispersed GRAM members are required."))
    truth_seed != prior_seed ||
        throw(ArgumentError("GRAM truth and prior nominal seeds must differ."))
    members = collect((prior_seed + 1):(prior_seed + n_dispersion))
    truth_seed in members &&
        throw(ArgumentError("GRAM truth seed must differ from every prior ensemble seed."))
    return GramSeedPlan(truth_seed, prior_seed, members)
end

function _prior_member_seeds(prior_seed::Int, n_dispersion::Int)::Vector{Int}
    n_dispersion >= 2 || throw(ArgumentError("At least two dispersed GRAM members are required."))
    return collect((prior_seed + 1):(prior_seed + n_dispersion))
end

function _to_initial_time(dt::DateTime)
    return SM.InitialTime(
        year=Dates.year(dt),
        month=Dates.month(dt),
        day=Dates.day(dt),
        hour=Dates.hour(dt),
        minute=Dates.minute(dt),
        second=Dates.second(dt) + Dates.millisecond(dt) / 1000.0
    )
end

function _build_gram_model(initial_time::SM.InitialTime, planet_name::String, seed::Int)
    return GRAMSuite.GRAMAtmosphereModel(
        ;
        planet_name=planet_name,
        initial_time=initial_time,
        seed=seed,
    )
end

function _gram_density_at_point(
    model,
    point::TrajectoryPoint;
    elapsed_time_s::Float64=point.elapsed_s,
    perturbed::Bool=false,
)::Float64
    gram = model.gram
    atmos = model.gram_atmosphere
    set_position! = Base.invokelatest(getproperty, gram, :set_position!)
    update! = Base.invokelatest(getproperty, gram, :update!)
    get_dynamics_state = Base.invokelatest(getproperty, gram, :get_dynamics_state)
    get_density_state = Base.invokelatest(getproperty, gram, :get_density_state)
    return lock(RuntimeServices.GRAM_LOCK) do
        Base.invokelatest(
            set_position!,
            atmos;
            height=point.alt_m * 1e-3,
            latitude=point.lat_deg,
            longitude=point.lon_deg,
            elapsed_time=elapsed_time_s
        )
        err = Base.invokelatest(update!, atmos)
        if err != 0
            get_error_message = Base.invokelatest(getproperty, gram, :get_error_message)
            throw(ErrorException("GRAM update failed (code=$err): $(Base.invokelatest(get_error_message))"))
        end
        dyn = Base.invokelatest(get_dynamics_state, atmos)
        if perturbed
            density = Base.invokelatest(get_density_state, atmos)
            return Float64(density.perturbedDensity)
        end
        return Float64(dyn.density)
    end
end

@inline function _elapsed_from_initial_s(point::TrajectoryPoint, initial_dt::DateTime)::Float64
    return Dates.value(point.dt - initial_dt) * 1.0e-3
end

function gram_truth_profiles(
    pair::TrajectoryPair,
    initial_dt::DateTime;
    planet_name::String="earth",
    seed::Int,
)::Tuple{Vector{Float64}, Vector{Float64}}
    initial_time = _to_initial_time(initial_dt)
    truth_model = _build_gram_model(initial_time, planet_name, seed)
    aero_truth = Vector{Float64}(undef, length(pair.aerocapture))
    edl_truth = Vector{Float64}(undef, length(pair.edl))
    @inbounds for i in eachindex(pair.aerocapture)
        point = pair.aerocapture[i]
        aero_truth[i] = _gram_density_at_point(
            truth_model,
            point;
            elapsed_time_s=_elapsed_from_initial_s(point, initial_dt),
            perturbed=true,
        )
    end
    @inbounds for i in eachindex(pair.edl)
        point = pair.edl[i]
        edl_truth[i] = _gram_density_at_point(
            truth_model,
            point;
            elapsed_time_s=_elapsed_from_initial_s(point, initial_dt),
            perturbed=true,
        )
    end
    return aero_truth, edl_truth
end

function gram_prior_samples(
    points::Vector{TrajectoryPoint},
    initial_dt::DateTime;
    planet_name::String="earth",
    prior_seed::Int=20_042,
    n_dispersion::Int=24
)::Vector{GramSample}
    dispersion_seeds = _prior_member_seeds(prior_seed, n_dispersion)
    initial_time = _to_initial_time(initial_dt)
    nominal_model = _build_gram_model(initial_time, planet_name, prior_seed)
    dispersed_models = [
        _build_gram_model(initial_time, planet_name, seed)
        for seed in dispersion_seeds
    ]

    out = Vector{GramSample}(undef, length(points))
    vals = Vector{Float64}(undef, n_dispersion)
    @inbounds for i in eachindex(points)
        elapsed_time_s = _elapsed_from_initial_s(points[i], initial_dt)
        mu = _gram_density_at_point(nominal_model, points[i]; elapsed_time_s)
        for j in 1:n_dispersion
            vals[j] = _gram_density_at_point(
                dispersed_models[j],
                points[i];
                elapsed_time_s,
                perturbed=true,
            )
        end
        sigma = max(std(vals), max(mu * 1.0e-6, 1.0e-12))
        out[i] = GramSample(mu, sigma)
    end
    return out
end
