import .config

function planet_data(ip)

    # try
    #     if haskey(ip, :planet)
    #         ip = ip[:planet]
    #     end
    # catch
    #     nothing
    # end

    if (ip == 0 || (typeof(ip) == String && cmp(lowercase(ip), "earth") == 0)) # Earth
        Rp_e = 6.3781e6            # equatorial radius, m
        Rp_p = 6.3568e6            # polar radius, m
        Rp_m = 6.3710e6            # volumetric mean radius, m
        mass = 5.97219e24          # mass, kg
        g_ref = 9.798              # acceleration due to gravity, m/s²
        ρ_ref = 1.225              # density, kg/m³
        μ = 3.9860e14              # gravitational parameter, m³/s²
        h_ref = 0 * 1e3            # reference altitude, m
        H = 8.5 * 1e3              # scale height, m
        R = 287.1                  # specific gas constant, J/(kg·K)
        γ = 1.4005                 # specific heat ratio
        T = 300.0                  # temperature, K
        p = 101400.0               # surface pressure, km/(m⋅s)
        J2 = 1.08263e-3            # Earth's dynamic form factor
        k = 1.83e-4                # Chapman heating coefficient, kg^0.5/m
        # k = 1.7623e-4            # Sutton - Graves heating coefficient, kg^0.5/m
        ω = [0.0, 0.0, 7.2921066e-5]   # Earth's rotation rate, rad/s
        μ_fluid = 1.5*10e-5        # kinematic viscosity, m²/s
        Lz = -9.8/1e3              # vertical temperature gradient, K/m
        name = "earth"
    elseif (ip == 1 || (typeof(ip) == String && cmp(lowercase(ip), "mars") == 0)) # Mars
        Rp_e = 3.3962e6            # equatorial radius, m
        Rp_p = 3.3762e6            # polar radius, m
        Rp_m = 3.3895e6            # volumetric mean radius, m
        mass = 6.4185e23           # mass, kg
        g_ref = 3.71               # acceleration due to gravity, m/s²
        ρ_ref = 8.7489231e-07      # density, kg/m³
        μ = 4.2828e13              # gravitational parameter, m³/s²
        h_ref = 90 * 1e3           # reference altitude, m
        H = 6.308278108 * 1e3      # scale height, m
        R = 188.92                 # specific gas constant, J/(kg·K)
        γ = 1.33                   # specific heat ratio
        T = 150.0                  # temperature, K
        p = 636.0                  # surface pressure, km/(m⋅s)
        J2 = 1.96045e-3            # Mars' dynamic form factor
        k = 1.898e-4               # Chapman heating coefficient, kg^0.5/m
        # k = 1.7623e-4            # Sutton - Graves heating coefficient, kg^0.5/m
        ω = [0.0, 0.0, 7.088236e-5]    # Mars' rotation rate, rad/s
        μ_fluid = 13.06*10e-6      # kinematic viscosity, m²/s
        Lz = -4.5/1e3              # vertical temperature gradient, K/m
        name = "mars"
    elseif (ip == 2 || (typeof(ip) == String && cmp(lowercase(ip), "venus") == 0)) # Venus
        Rp_e = 6.0518e6            # equatorial radius, m
        Rp_p = 6.0518e6            # polar radius, m
        Rp_m = 6.0518e6            # volumetric mean radius, m
        mass = 4.8685e24           # mass, kg
        g_ref = 8.87               # acceleration due to gravity, m/s²
        ρ_ref = 65.0               # density, kg/m³
        μ = 3.24858592e14          # gravitational parameter, m³/s²
        h_ref = 0 * 1e3           # reference altitude, m
        H = 15.9 * 1e3            # scale height, m
        R = 188.92                 # specific gas constant, J/(kg·K)
        γ = 1.2857                 # specific heat ratio
        T = 100.0                  # temperature, K
        p = 9200000.0              # surface pressure, km/(m⋅s)
        J2 = 4.458e-6              # Venus' dynamic form factor
        k = 1.896e-4               # Chapman heating coefficient, kg^0.5/m
        # k = 1.7623e-4            # Sutton - Graves heating coefficient, kg^0.5/m
        ω = [0.0, 0.0, -2.99e-7]  # Venus' rotation rate, rad/s
        μ_fluid = 2.0*10e-6      # kinematic viscosity, m²/s
        Lz = -10.7/1e3             # vertical temperature gradient, K/m
        name = "venus"
    elseif (ip == 3 || (typeof(ip) == String && cmp(lowercase(ip), "sun") == 0)) # Sun
        Rp_e = 6.9634e8            # equatorial radius, m
        Rp_p = 6.9634e8            # polar radius, m
        Rp_m = 6.9634e8            # volumetric mean radius, m
        mass = 1.9891e30           # mass, kg
        g_ref = 274                # m/s^2
        ρ_ref = 0
        μ = 1.32712440018e20       # gravitational parameter, m^3/s^2
        h_ref = 0
        H = 0
        R = 0
        γ = 0
        T = 0
        p = 0
        J2 = 0
        k = 0
        ω = [0, 0, 0]
        μ_fluid = 0                # kinematic viscosity, m²/s
        Lz = 0
        name = "sun"
    end

    # println(model.planet)
    config.model.planet = config.Planet(Rp_e, Rp_p, Rp_m, mass, p, k, ω, g_ref, ρ_ref, h_ref, H, R, γ, T, J2, μ, μ_fluid, Lz, name)
    
    return config.model.planet
end