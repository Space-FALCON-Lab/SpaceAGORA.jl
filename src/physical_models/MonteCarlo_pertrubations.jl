## Note: HOW THE MONTECARLO PERTURBATIONS ARE SET UP
## the seed is the monte carlo number, the random number generator generates config.counter_random number, which increase of +1 every time the generator is called.
## the generates output is always the last number of the list. When config.counter_random number reaches 100, it is re-initialized to 1 to avoid storage problem and slowing down the function.
## Changing the config.counter_random number results in changing the random number generates, otherwise this would be the same everytime, biasing the results.

include("../config.jl")

using Distributions

function unifrom_distribution(dispersion, number=0)
    """

    """

    Random.seed!(config.index_MonteCarlo)

    return rand(Uniform(-dispersion, dispersion), number)[end]
end

function gaussian_distribution(μ, σ, number = 0)
    """

    """

    config.counter += 1

    Random.seed!(config.counter_random)

    number = rand(1:500)

    Random.seed!(config.counter_random)

    return rand(Normal(μ, σ), number)[end]
end

function monte_carlo_aerodynamics(CL_body, CD_body, args)
    """

    """

    uncertainty_CD, uncertainty_CL = args.CD_dispersion/100, args.CL_ispersion/100
    CD_body += unifrom_distribution(CD_body*uncertainty_CD, 1)
    CL_body += unifrom_distribution(CL_body*uncertainty_CL, 1)

    return CL_body, CD_body
end

function monte_carlo_density(density, args)
    """

    """

    Random.seed!(config.index_MonteCarlo)

    density = rand(Uniform(density*0.5, density*2), 3)

    return density
end

function monte_carlo_initial_condition(state, args)
    """

    """

    state[:Apoapsis] += unifrom_distribution(args[:ra_dispersion], 3)
    state[:Periapsis] += unifrom_distribution(args[:rp_dispersion], 4)
    state[:Inclination] += unifrom_distribution(args[:i_dispersion], 5)
    state[:Ω] += unifrom_distribution(args[:Ω_dispersion], 6)
    state[:ω] += unifrom_distribution(args[:ω_dispersion], 7)
    state[:vi] += unifrom_distribution(args[:vi_dispersion], 8)

    return state
end

function monte_carlo_true_anomaly(state, args)
    """

    """
    
    state[:vi] += unifrom_distribution(args[:vi_dispersion], 9)

    return state
end

function monte_carlo_guidance_closedform(state, args)
    """

    """

    state[:ra] += unifrom_distribution(args[:ra_dispersion_gnc], 10)
    state[:rp] += unifrom_distribution(args[:rp_dispersion_gnc], 11)
    state[:i] += unifrom_distribution(args[:i_dispersion_gnc], 12)
    state[:Ω] += unifrom_distribution(args[:Ω_dispersion_gnc], 13)
    state[:ω] += unifrom_distribution(args[:ω_dispersion_gnc], 14)
    state[:vi] += unifrom_distribution(args[:vi_dispersion_gnc], 15)

    return state
end

function monte_carlo_guidance_environment(ρ, T, S, args)
    """

    """

    mean_ρ, mean_T, mean_S = args[:ρ_mudispersion_gnc]/100, args[:T_mudispersion_gnc]/100, args[:S_mudispersion_gnc]/100
    std_ρ, std_T, std_S = args[:ρ_sigmadispersion_gnc]/100, args[:T_sigmadispersion_gnc]/100, args[:S_sigmadispersion_gnc]/100

    ρ += gaussian_distribution(ρ*mean_ρ, ρ*std_ρ, 16)
    T += gaussian_distribution(T*mean_T, T*std_T, 17)
    S += gaussian_distribution(S*mean_S, S*std_S, 18)

    return ρ, T, S
end