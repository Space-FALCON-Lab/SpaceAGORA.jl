include(joinpath(@__DIR__, "common.jl"))

function _mc_int_env(name::String, default::Int)::Int
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return default
    parsed = tryparse(Int, raw)
    parsed === nothing && throw(ArgumentError("$(name) must be an integer; got $(repr(raw))."))
    parsed > 0 || throw(ArgumentError("$(name) must be > 0; got $(parsed)."))
    return parsed
end

const MC_SMOKE = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"
const MC_SAMPLES = _mc_int_env("SPACEAGORA_MC_SAMPLES", MC_SMOKE ? 1 : 8)
const MC_THREADS = _mc_int_env("SPACEAGORA_MC_THREADS", min(Threads.nthreads(), MC_SAMPLES))

function _mc_mission_time_s(default_time_s::Float64)::Float64
    !MC_SMOKE && return default_time_s
    raw = strip(get(ENV, "SPACEAGORA_EXAMPLE_SMOKE_MISSION_TIME", "120.0"))
    parsed = tryparse(Float64, raw)
    if parsed === nothing || !(parsed > 0.0)
        return min(default_time_s, 120.0)
    end
    return min(default_time_s, parsed)
end

function make_config_for_seed(seed::Int)
    planet = make_no_gram_planet(:earth)
    phase_deg = 175.0 + 0.05 * (seed - 1)
    periapsis_offset_m = 50.0 * (seed - 1)

    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.05, 2.05, 2.8),
        panel_dims=(0.01, 5.7 / 2.0, 1.0),
        bus_mass=620.0,
        panel_mass_each=10.0,
        panel_offset_y=2.05 / 2.0 + 5.7 / 4.0,
        ic=InitialCondition(
            ra=56_378.7978559e3,
            rp=planet.Rp_e + 200_590.0 + periapsis_offset_m,
            i=89.876,
            ω=75.505,
            Ω=104.115,
            ν=phase_deg
        ),
        prop_mass=200.0,
        id=1
    )

    args = make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=_mc_mission_time_s(30.0 * 60.0),
        initial_time=InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0),
        dynamic_effectors=(InverseSquaredJ2GravityModel(),),
        density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(),
        orientation_sim=false,
        keplerian=true,
        EI_km=120.0,
        verbose=false,
        results=false
    )

    return args
end

seeds = collect(1:MC_SAMPLES)
println("Running $(length(seeds)) Monte Carlo samples with $(MC_THREADS) worker task(s)")

result = run_monte_carlo(seeds; threads=MC_THREADS) do seed
    run_simulation(make_config_for_seed(seed))
    return (; seed)
end

println("Successful samples: $(length(result.successful))")
println("Failed samples: $(length(result.failed))")
println("Elapsed wall time: $(round(result.elapsed_s; digits=3)) s")

if !isempty(result.failed)
    for sample in result.failed
        println("  seed $(sample.seed) failed: $(sample.error)")
    end
    error("Monte Carlo example failed for $(length(result.failed)) sample(s).")
end
