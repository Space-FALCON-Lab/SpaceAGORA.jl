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
        # polyfit_coeffs = [-1.7539409645214832e-57, 2.735656076315809e-53, -1.8243490769488347e-49, 6.504765617793163e-46, -1.1637408657034938e-42, 8.044884138893168e-41, 4.264962263039017e-36, -7.651115834387683e-33, -3.188248308052816e-30, 3.8370830656820503e-26, -8.557502178008995e-23, 1.137879849173412e-19, -1.0408232216096158e-16, 6.834085016894604e-14, -3.2506596548183e-11, 1.1089006707870246e-08, -2.639423772958483e-06, 0.0004165844083994442, -0.03967261693733797, 1.8349343859319074, -38.14918904018883]
        
        # 50 km to 500 km
        polyfit_coeffs = [-2.05897273251208688e+01,
                            -8.56737992269037107e+00,
                            3.45019353462444567e+00,
                            -1.76485096058359936e+00,
                            6.01801952591328626e-01,
                            7.24429484383713046e-02,
                            -3.30361927785729759e-01,
                            3.07509524503914389e-01,
                            -1.58408358148278777e-01,
                            7.77052279665210299e-03,
                            7.81913274877799047e-02,
                            -9.03365513831131534e-02,
                            5.60317820649417070e-02,
                            -1.22200602431319916e-02,
                            -1.47378888840894624e-02,
                            1.75437183414442650e-02,
                            -4.28819900201210177e-03,
                            -1.10920129900025574e-02,
                            1.82659646193858577e-02,
                            -1.51388155198109841e-02,
                            6.19700872461139817e-03]
        name = "earth"
    elseif (ip == 1 || (typeof(ip) == String && cmp(lowercase(ip), "mars") == 0)) # Mars
        Rp_e = 3.3962e6 #3.3962    # equatorial radius, m
        Rp_p = 3.3762e6 #3.3762    # polar radius, m
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
        SF = 586.2                 # Solar irradiance, W/m²
        topography_function = Mars_elevation! # Mars topography function

        # -1 km to 2000 km
        # polyfit_coeffs = [-8.278592174668491e-43, 1.2598495030132498e-38, -8.634065871212132e-35, 3.5185552646901455e-31, -9.480197229347404e-28, 1.7753104600795092e-24, -2.3622107295909874e-21, 2.2393603867716714e-18, -1.487031340144351e-15, 6.592111911218399e-13, -1.714014789283248e-10, 1.3556252797088945e-08, 5.196239221937857e-06, -0.0012393556758398866, -0.0500835105059738, -4.213431227716942]
        
        # 100 km to 300 km
        # polyfit_coeffs = [-3.691310097181554e-58, 5.819173546214448e-54, -3.9285937578286423e-50, 1.4222601230188116e-46, -2.606951392190571e-43, 3.2943551967480965e-41, 9.394166176413728e-37, -1.7651753457891617e-33, -5.79069281873952e-31, 8.639557954110502e-27, -1.991207114225621e-23, 2.7207390647640917e-20, -2.5611296697872007e-17, 1.7386922029136165e-14, -8.619727907575625e-12, 3.1040218147963276e-09, -7.949080301839893e-07, 0.00013834108975291533, -0.014729001168514675, 0.6707044510751348, -19.414578139119545]
        
        # 50 km to 300 km
        polyfit_coeffs = [-22.904887634115337,
                                    -11.841971600850178,
                                    2.449273531108809,
                                    0.07122528331622979,
                                    -0.3143478095124491,
                                    0.1734745257301972,
                                    -0.14396272246410283,
                                    0.06104818577407287,
                                    0.0738268254164405,
                                    -0.11983534975043726,
                                    0.052506209635111954,
                                    0.0022260853200571835,
                                    -0.0245022530854842]

        name = "mars"
    elseif (ip == 2 || (typeof(ip) == String && cmp(lowercase(ip), "venus") == 0)) # Venus
        Rp_e = 6.0518e6            # equatorial radius, m
        Rp_p = 6.0518e6            # polar radius, m
        Rp_m = 6.0518e6            # volumetric mean radius, m
        mass = 4.8685e24           # mass, kg
        g_ref = 8.87               # acceleration due to gravity, m/s²
        ρ_ref = 65.0               # density, kg/m³
        μ = 3.24858592e14          # gravitational parameter, m³/s²
        h_ref = 0 * 1e3            # reference altitude, m
        H = 15.9 * 1e3             # scale height, m
        R = 188.92                 # specific gas constant, J/(kg·K)
        γ = 1.2857                 # specific heat ratio
        T = 100.0                  # temperature, K
        p = 9200000.0              # surface pressure, km/(m⋅s)
        J2 = 4.458e-6              # Venus' dynamic form factor
        k = 1.896e-4               # Chapman heating coefficient, kg^0.5/m
        # k = 1.7623e-4            # Sutton - Graves heating coefficient, kg^0.5/m
        ω = [0.0, 0.0, -2.99e-7]   # Venus' rotation rate, rad/s
        μ_fluid = 2.0*10e-6        # kinematic viscosity, m²/s
        Lz = -10.7/1e3             # vertical temperature gradient, K/m
        α = deg2rad(272.76)        # Right ascension of the north pole of rotation, radians
        δ = deg2rad(67.16)         # Declination of the north pole of rotation, radians
        SF = 2601.3                # Solar irradiance, W/m²
        topography_function = Venus_elevation! # Venus topography function
        # polyfit_coeffs = [1.295014716586507e-57, -1.920381283790201e-53, 1.2024671159968765e-49, -3.931503383921753e-46, 5.985870736864543e-43, 2.115956905107091e-40, -2.4659597875857534e-36, 3.0591710987549437e-33, 3.951465781537392e-30, -1.8949093746237393e-26, 3.123829612747949e-23, -2.928033666820754e-20, 1.5168683041510048e-17, -1.5135241597177884e-15, -3.865230229956326e-12, 3.1328117105612896e-09, -1.2501690556294552e-06, 0.00028978339946121796, -0.03741075092352375, 2.149847471180469, -43.08275565785116]
        
        # 50 km to 250 km
        polyfit_coeffs = [-17.17603495187312,
                            -15.971487996237625,
                            3.1916728583530163,
                            0.8296122739608919,
                            -0.891116626156514,
                            -0.018054118127076735,
                            0.2016004439406335,
                            0.0537350005662433,
                            -0.06184791826147179,
                            -0.002221184567720138,
                            -0.03852447322489234,
                            0.11389766896205009,
                            -0.023720859333601994]
        
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
        polyfit_coeffs = zeros(1)
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
        polyfit_coeffs = zeros(1)
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
        topography_function = (args, Clm, Slm, latitude, longitude) -> 0.0 # Jupiter topography function
        polyfit_coeffs = zeros(1)
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
        topography_function = (args, Clm, Slm, latitude, longitude) -> 0.0 # Saturn topography function
        polyfit_coeffs = zeros(1)
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
        topography_function = (args, Clm, Slm, latitude, longitude) -> 0.0 # Titan topography function
        # polyfit_coeffs = [1.7989756686197253e-58, -2.7298975030491325e-54, 1.7620522402686604e-50, -6.025021166267467e-47, 1.0056316643424087e-43, 9.494104496406468e-42, -3.8472088727076255e-37, 6.051435602297366e-34, 4.074478639170247e-31, -3.244699052533356e-27, 6.66877802035039e-24, -8.360025139024445e-21, 7.301165978344981e-18, -4.650857357165472e-15, 2.1978197729328097e-12, -7.705014392936314e-10, 1.9713879988437584e-07, -3.551889476633975e-05, 0.004248542489215875, -0.3277965440319509, 8.128293001726805] 
        
        # 50 km to 2000 km
        polyfit_coeffs = [-1.94278696179288133e+01,
                            -1.44426058071846857e+01,
                            2.69604850426823095e+00,
                            -3.36162464509102088e-01,
                            9.00778126111024535e-02,
                            -1.36262683746282226e-01,
                            2.32875124056659416e-01,
                            -2.26385345471503546e-01,
                            1.18549197857537245e-01,
                            -2.89659248379974324e-02,
                            1.65054791585107323e-02,
                            -4.78606382329083910e-02,
                            6.50600435627339962e-02,
                            -5.03265704035375031e-02,
                            2.26198266668309125e-02,
                            -6.45321729948887882e-03,
                            5.09967575456080676e-03,
                            -9.43672721749253159e-03,
                            9.48561230043162319e-03,
                            -6.20859370479269826e-03,
                            1.20768474131884517e-03]
        
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
                            polyfit_coeffs,
                            topography_function)
    
    return planet
end