include("../physical_models/Planet_shapes.jl")

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
        μ = 3.986004418e14              # gravitational parameter, m³/s²
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
        α = 0.0                    # Right ascension of the north pole of rotation, radians
        SF = 1361.0                # Solar irradiance, W/m²
        δ = pi/2                    # Declination of the north pole of rotation, radians
        topography_function = Earth_elevation! # Earth topography function
        name = "earth"
    elseif (ip == 1 || (typeof(ip) == String && cmp(lowercase(ip), "mars") == 0)) # Mars
        Rp_e = 3.3962e6 #3.3962           # equatorial radius, m
        Rp_p = 3.3762e6 #3.3762           # polar radius, m
        Rp_m = 3.3895e6            # volumetric mean radius, m
        mass = 6.4169e23           # mass, kg
        g_ref = 3.73               # acceleration due to gravity, m/s²
        ρ_ref = 8.7489231e-07      # density, kg/m³
        μ = 4.282837362069909e13 # 4.2828314258067e13              # gravitational parameter, m³/s²
        h_ref = 90 * 1e3           # reference altitude, m
        H = 6.308278108 * 1e3      # scale height, m
        R = 188.92                 # specific gas constant, J/(kg·K)
        γ = 1.33                   # specific heat ratio
        T = 150.0                  # temperature, K
        p = 636.0                  # surface pressure, km/(m⋅s)
        J2 = 1.96045e-3            # Mars' dynamic form factor
        k = 1.898e-4               # Chapman heating coefficient, kg^0.5/m
        # k = 1.7623e-4            # Sutton - Graves heating coefficient, kg^0.5/m
        ω = [0.0, 0.0, 7.08823596e-5]    # Mars' rotation rate, rad/s
        μ_fluid = 13.06*10e-6      # kinematic viscosity, m²/s
        Lz = -4.5/1e3              # vertical temperature gradient, K/m
        α = deg2rad(317.68143)     # Right ascension of the north pole of rotation, radians
        δ = deg2rad(52.88650)      # Declination of the north pole of rotation, radians
        SF = 586.2                # Solar irradiance, W/m²
        topography_function = Mars_elevation! # Mars topography function
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
        α = deg2rad(272.76)        # Right ascension of the north pole of rotation, radians
        δ = deg2rad(67.16)         # Declination of the north pole of rotation, radians
        SF = 2601.3                # Solar irradiance, W/m²
        topography_function = Venus_elevation! # Venus topography function
        name = "venus"
    elseif (ip == 3 || (typeof(ip) == String && cmp(lowercase(ip), "sun") == 0)) # Sun
        Rp_e = 6.9634e8            # equatorial radius, m
        Rp_p = 6.9634e8            # polar radius, m
        Rp_m = 6.9634e8            # volumetric mean radius, m
        mass = 1.9891e30           # mass, kg
        g_ref = 274                # m/s^2
        ρ_ref = 0
        μ = 1.3271244001799e20       # gravitational parameter, m^3/s^2
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
        α = 0.0
        δ = 0.0
        SF = 0
        topography_function = (args, Clm, Slm, latitude, longitude) -> 0.0
        name = "sun"
    elseif (ip == 4 || (typeof(ip) == String && cmp(lowercase(ip), "moon") == 0)) # Moon
        Rp_e = 1.7381e6            # equatorial radius, m
        Rp_p = 1.7360e6            # polar radius, m
        Rp_m = 1.7374e6            # volumetric mean radius, m
        mass = 0.07346e24           # mass, kg
        g_ref = 1.62                # acceleration due to gravity, m/s²
        ρ_ref = 0
        μ = 4.9028005821478e12              # gravitational parameter, m³/s²
        h_ref = 0
        H = 0
        R = 0
        γ = 0
        T = 0
        p = 0
        J2 = 202.7e-6
        k = 0
        ω = [0, 0, 0]
        μ_fluid = 0                # kinematic viscosity, m²/s
        Lz = 0
        α = 0.0
        δ = 0.0
        SF = 1361.0                # Solar irradiance, W/m²
        topography_function = (args, Clm, Slm, latitude, longitude) -> 0.0
        name = "moon"
    elseif (ip == 5 || (typeof(ip) == String && cmp(lowercase(ip), "jupiter") == 0))
        Rp_e = 7.1492e7
        Rp_p = 6.6854e7
        Rp_m = 6.9911e7
        mass = 1.89813e27
        g_ref = 25.92 # m/s^2
        ρ_ref = 0
        μ = 1.26686534e17 # gravitational parameter, m^3/s^2
        h_ref = 0 * 10e3
        H = 0
        R = 0
        γ = 0
        T = 0
        p = 0
        J2 = 14736e-6
        k = 0
        ω = [0, 0, 1.758e-4]
        μ_fluid = 0                # kinematic viscosity, m²/s
        Lz = 0
        α = 0.0
        δ = 0.0
        SF = 50.26                # Solar irradiance, W/m²
        name = "jupiter"
    elseif (ip == 6 || (typeof(ip) == String && cmp(lowercase(ip), "saturn") == 0))
        Rp_e = 6.0268e7
        Rp_p = 5.4364e7
        Rp_m = 5.8232e7
        mass = 5.68319e26
        g_ref = 11.19 # m/s^2
        ρ_ref = 0
        μ = 3.7931187e16 # gravitational parameter, m^3/s^2
        h_ref = 0 * 10e3
        H = 0
        R = 0
        γ = 0
        T = 0
        p = 0
        J2 = 16290e-6
        k = 0
        ω = [0, 0, 1.65e-4]
        μ_fluid = 0                # kinematic viscosity, m²/s
        Lz = 0
        α = deg2rad(40.589)
        δ = deg2rad(83.537)
        SF = 14.82                # Solar irradiance, W/m²
        name = "saturn"
    elseif (ip == 7 || (typeof(ip) == String && cmp(lowercase(ip), "titan") == 0))
        Rp_e = 2.575e6
        Rp_p = 2.575e6
        Rp_m = 2.575e6
        mass = 1.3452e23
        g_ref = 1.352 # m/s^2
        ρ_ref = 5.3
        μ = 8.981e12 # gravitational parameter, m^3/s^2
        h_ref = 0 * 10e3
        H = 21.0e3
        R = 290.0
        γ = 1.3846
        T = 94
        p = 146.7
        J2 = 3.15e-5
        k = 1.74e-4
        ω = [0, 0, 7.37e-6]
        μ_fluid = 0                # kinematic viscosity, m²/s
        Lz = -1.352/1e3
        α = deg2rad(39.4827)
        δ = deg2rad(83.4279)
        SF = 14.82                # Solar irradiance, W/m², probably same as Saturn but couldn't find real data
        name = "titan"
    end

    # Derived in References/J2000_to_pci.mlx(.m)
    # Converts from J2000 to Planet Centered Inertial (PCI) 
    # frame based on the planet's North pole of rotation
    # α = Right ascension of the north pole of rotation, radians
    # δ = Declination of the north pole of rotation, radians
    if name == "earth"
        J2000_to_pci = [1 0 0; 0 1 0; 0 0 1]
    else
        σ1 = sqrt(cos(δ)^4 + cos(δ)^2*sin(δ)^2)
        J2000_to_pci = SMatrix{3, 3, Float64}([-sin(α) cos(α) 0;
                        -cos(δ)*cos(α)*sin(δ)/σ1 -cos(δ)*sin(α)*sin(δ)/σ1 cos(δ)^2/σ1;
                        cos(δ)*cos(α) cos(δ)*sin(α) sin(δ)])   
    end 
    planet = config.Planet(Rp_e, 
                                        Rp_p, 
                                        Rp_m, 
                                        mass, 
                                        p, 
                                        k, 
                                        ω, 
                                        g_ref, 
                                        ρ_ref, 
                                        h_ref, 
                                        H, 
                                        R, 
                                        γ, 
                                        T, 
                                        J2, 
                                        μ, 
                                        μ_fluid, 
                                        Lz, 
                                        α, 
                                        δ, 
                                        J2000_to_pci,
                                        MMatrix{3, 3, Float64}(zeros(3,3)), 
                                        zeros(3, 3), 
                                        zeros(3, 3), 
                                        zeros(3, 3),
                                        zeros(3, 3),
                                        name,
                                        zeros(3, 3),
                                        zeros(3, 3),
                                        zeros(3, 3),
                                        zeros(3, 3),
                                        zeros(3, 3),
                                        zeros(3, 3),
                                        zeros(3, 3),
                                        zeros(3),
                                        zeros(3),
                                        topography_function)
    
    return planet
end