module Planets
    include("../physical_models/Planet_shapes.jl")
    using ..AbstractTypes: AbstractPlanet
    using StaticArrays
    using CSV
    using SPICE
    export Earth, Mars, Venus, Titan
    δ(i, j) = ==(i, j)

    @kwdef mutable struct TopographyHarmonicsWorkspace
        Clm::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        Slm::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        Plm::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        fn::Vector{Float64} = [0.0, 0.0, 0.0]
        G::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        H::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    end # struct HarmonicsWorkspace

    @kwdef struct Earth <: AbstractPlanet
        name::String = "Earth" # Name of the planet
        Rp_e::Float64 = 6.3781e6 # Equatorial radius in meters
        Rp_p::Float64 = 6.3568e6 # Polar radius in meters
        Rp_m::Float64 = 6.371e6  # Mean radius in meters
        mass::Float64  = 5.972e24 # Mass in kg
        p::Float64 = 101325.0 # Surface pressure in Pascals
        k::Float64 = 1.83e-4 # Chapman heating coefficient, kg^0.5/m
        ω::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 7.2921066e-5) # Angular velocity vector in rad/s
        μ::Float64 = 3.986004418e14 # Standard gravitational parameter in m^3/s^2
        J2::Float64 = 1.08263e-3 # J2 coefficient
        g_ref::Float64 = 9.80665 # Standard gravity in m/s^2
        ρ_ref::Float64 = 1.225 # Sea level atmospheric density in kg/m^3
        h_ref::Float64 = 0.0 # Reference altitude for exponential atmosphere in meters
        H::Float64 = 8.5e3 # Scale height for exponential atmosphere in meters
        R::Float64 = 287.1 # Specific gas constant for dry air in J/(kg*K)
        T_ref::Float64 = 288.15 # Reference temperature for exponential atmosphere in K
        γ::Float64 = 1.4 # Ratio of specific heats for air
        T::Float64 = 300.0 # Surface temperature in K
        μ_fluid::Float64 = 1.5e-5 # Dynamic viscosity of air in Pa*s
        Lz::Float64 = -9.8e-3 # Vertical temperature gradient in K/m for calculating temperature at altitude
        α::Float64 = 0.0 # Right-ascension of north pole relative to J2000 in degrees
        δ::Float64 = 0.0 # Declination of north pole relative to J2000 in degrees
        J2000_to_pci::SMatrix{3, 3, Float64} = SMatrix{3, 3, Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]) # Rotation matrix from J2000 to planet-centered inertial frame
        L_PI::MMatrix{3, 3, Float64} = MMatrix{3, 3, Float64}([0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]) # Rotation matrix from planet-centered inertial frame to planet-centered, planet-fixed frame, function of time
        topography_workspace::TopographyHarmonicsWorkspace = TopographyHarmonicsWorkspace() # Workspace for calculating topography harmonics
        polyfit_coeffs::Vector{Float64} = [-1.7539409645214832e-57, 2.735656076315809e-53, -1.8243490769488347e-49, 6.504765617793163e-46, -1.1637408657034938e-42, 8.044884138893168e-41, 4.264962263039017e-36, -7.651115834387683e-33, -3.188248308052816e-30, 3.8370830656820503e-26, -8.557502178008995e-23, 1.137879849173412e-19, -1.0408232216096158e-16, 6.834085016894604e-14, -3.2506596548183e-11, 1.1089006707870246e-08, -2.639423772958483e-06, 0.0004165844083994442, -0.03967261693733797, 1.8349343859319074, -38.14918904018883]
        topography_function::Function = Earth_elevation! # Function to calculate elevation based on topography
    end # struct Earth

    @kwdef struct Mars <: AbstractPlanet
        name::String = "Mars" # Name of the planet
        Rp_e::Float64 = 3.3962e6 # Equatorial radius in meters
        Rp_p::Float64 = 3.3762e6 # Polar radius in meters
        Rp_m::Float64 = 3.3895e6 # Mean radius in meters
        mass::Float64  = 0.64171e24 # Mass in kg
        p::Float64 = 636.0 # Surface pressure in Pascals
        k::Float64 = 1.898e-4 # Chapman heating coefficient, kg^0.5/m
        ω::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 7.08823596e-5) # Angular velocity vector in rad/s
        μ::Float64 = 4.2828314e13 # Standard gravitational parameter in m^3/s^2
        J2::Float64 = 1.96045e-3 # J2 coefficient
        g_ref::Float64 = 3.72076 # Standard gravity in m/s^2
        ρ_ref::Float64 = 8.7489231e-7 # Sea level atmospheric density in kg/m^3
        h_ref::Float64 = 90.0e3 # Reference altitude for exponential atmosphere in meters
        H::Float64 = 6.308278108e3 # Scale height for exponential atmosphere in meters
        R::Float64 = 188.92 # Specific gas constant for CO2 in J/(kg*K)
        T_ref::Float64 = 150.0 # Reference temperature for exponential atmosphere in K
        γ::Float64 = 1.33 # Ratio of specific heats for CO2
        T::Float64 = 150.0 # Surface temperature in K
        μ_fluid::Float64 = 13.06e-6 # Dynamic viscosity of CO2 in Pa*s
        Lz::Float64 = -4.5e-3 # Vertical temperature gradient in K/m for calculating temperature at altitude
        α::Float64 = deg2rad(317.68143) # Right-ascension of north pole
        δ::Float64 = deg2rad(52.88650) # Declination of north pole
        J2000_to_pci::SMatrix{3, 3, Float64} = SMatrix{3, 3, Float64}([0.67330 0.7394 0.0; -0.5896 0.5369 0.6034; 0.4462 -0.4062 0.7974]) # Rotation matrix from J2000 to planet-centered inertial frame # Rotation matrix from J2000 to planet-centered inertial frame
        L_PI::MMatrix{3, 3, Float64} = MMatrix{3, 3, Float64}([0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]) # Rotation matrix from planet-centered inertial frame to planet-centered, planet-fixed frame, function of time
        topography_workspace::TopographyHarmonicsWorkspace = TopographyHarmonicsWorkspace() # Workspace for calculating topography harmonics
        polyfit_coeffs::SVector{21, Float64} = SVector{21, Float64}(-3.691310097181554e-58, 5.819173546214448e-54, -3.9285937578286423e-50, 1.4222601230188116e-46, -2.606951392190571e-43, 3.2943551967480965e-41, 9.394166176413728e-37, -1.7651753457891617e-33, -5.79069281873952e-31, 8.639557954110502e-27, -1.991207114225621e-23, 2.7207390647640917e-20, -2.5611296697872007e-17, 1.7386922029136165e-14, -8.619727907575625e-12, 3.1040218147963276e-09, -7.949080301839893e-07, 0.00013834108975291533, -0.014729001168514675, 0.6707044510751348, -19.414578139119545)
    end # struct Mars

    @kwdef struct Venus <: AbstractPlanet
        name::String = "Venus" # Name of the planet
        Rp_e::Float64 = 6.0518e6 # Equatorial radius in meters
        Rp_p::Float64 = 6.0518e6 # Polar radius in meters
        Rp_m::Float64 = 6.0518e6 # Mean radius in meters
        mass::Float64  = 4.8685e24 # Mass in kg
        p::Float64 = 9.2e6 # Surface pressure in Pascals
        k::Float64 = 1.896e-4 # Chapman heating coefficient, kg^0.5/m
        ω::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, -2.99e-7) # Angular velocity vector in rad/s
        μ::Float64 = 3.24858592e14 # Standard gravitational parameter in m^3/s^2
        J2::Float64 = 4.458e-6 # J2 coefficient
        g_ref::Float64 = 8.87 # Standard gravity in m/s^2
        ρ_ref::Float64 = 65.0 # Sea level atmospheric density in kg/m^3
        h_ref::Float64 = 0.0 # Reference altitude for exponential atmosphere in meters
        H::Float64 = 15.9e3 # Scale height for exponential atmosphere in meters
        R::Float64 = 188.92 # Specific gas constant in J/(kg*K)
        T_ref::Float64 = 100.0 # Reference temperature in K
        γ::Float64 = 1.2857 # Ratio of specific heats
        T::Float64 = 100.0 # Surface temperature in K
        μ_fluid::Float64 = 2.0e-6 # Dynamic viscosity in Pa*s
        Lz::Float64 = -10.7e-3 # Vertical temperature gradient in K/m
        α::Float64 = deg2rad(272.76) # Right-ascension of north pole
        δ::Float64 = deg2rad(67.16) # Declination of north pole
        J2000_to_pci::SMatrix{3, 3, Float64} = SMatrix{3, 3, Float64}(
            [0.9988399975 0.0481524597 0.0;
             -0.0443769404 0.9205233406 0.3881590739;
             0.0186908142 -0.3877088084 0.92159239]
        )
        L_PI::MMatrix{3, 3, Float64} = MMatrix{3, 3, Float64}([0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0])
        topography_workspace::TopographyHarmonicsWorkspace = TopographyHarmonicsWorkspace()
        polyfit_coeffs::Vector{Float64} = [1.295014716586507e-57, -1.920381283790201e-53, 1.2024671159968765e-49, -3.931503383921753e-46, 5.985870736864543e-43, 2.115956905107091e-40, -2.4659597875857534e-36, 3.0591710987549437e-33, 3.951465781537392e-30, -1.8949093746237393e-26, 3.123829612747949e-23, -2.928033666820754e-20, 1.5168683041510048e-17, -1.5135241597177884e-15, -3.865230229956326e-12, 3.1328117105612896e-9, -1.2501690556294552e-6, 0.00028978339946121796, -0.03741075092352375, 2.149847471180469, -43.08275565785116]
    end # struct Venus

    @kwdef struct Titan <: AbstractPlanet
        name::String = "Titan" # Name of the planet
        Rp_e::Float64 = 2.575e6 # Equatorial radius in meters
        Rp_p::Float64 = 2.575e6 # Polar radius in meters
        Rp_m::Float64 = 2.575e6 # Mean radius in meters
        mass::Float64  = 1.3452e23 # Mass in kg
        p::Float64 = 146.7 # Surface pressure in Pascals
        k::Float64 = 1.74e-4 # Chapman heating coefficient, kg^0.5/m
        ω::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 7.37e-6) # Angular velocity vector in rad/s
        μ::Float64 = 8.981e12 # Standard gravitational parameter in m^3/s^2
        J2::Float64 = 3.15e-5 # J2 coefficient
        g_ref::Float64 = 1.352 # Standard gravity in m/s^2
        ρ_ref::Float64 = 5.3 # Sea level atmospheric density in kg/m^3
        h_ref::Float64 = 0.0 # Reference altitude for exponential atmosphere in meters
        H::Float64 = 21.0e3 # Scale height for exponential atmosphere in meters
        R::Float64 = 290.0 # Specific gas constant in J/(kg*K)
        T_ref::Float64 = 94.0 # Reference temperature in K
        γ::Float64 = 1.3846 # Ratio of specific heats
        T::Float64 = 94.0 # Surface temperature in K
        μ_fluid::Float64 = 0.0 # Dynamic viscosity in Pa*s
        Lz::Float64 = -1.352e-3 # Vertical temperature gradient in K/m
        α::Float64 = deg2rad(39.4827) # Right-ascension of north pole
        δ::Float64 = deg2rad(83.4279) # Declination of north pole
        J2000_to_pci::SMatrix{3, 3, Float64} = SMatrix{3, 3, Float64}(
            [-0.6358452054 0.7718166069 0.0;
             -0.7667447037 -0.6316668225 0.1144534171;
             0.0883370481 0.0727746565 0.9934286161]
        )
        L_PI::MMatrix{3, 3, Float64} = MMatrix{3, 3, Float64}([0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0])
        topography_workspace::TopographyHarmonicsWorkspace = TopographyHarmonicsWorkspace()
        polyfit_coeffs::Vector{Float64} = [1.7989756686197253e-58, -2.7298975030491325e-54, 1.7620522402686604e-50, -6.025021166267467e-47, 1.0056316643424087e-43, 9.494104496406468e-42, -3.8472088727076255e-37, 6.051435602297366e-34, 4.074478639170247e-31, -3.244699052533356e-27, 6.66877802035039e-24, -8.360025139024445e-21, 7.301165978344981e-18, -4.650857357165472e-15, 2.1978197729328097e-12, -7.705014392936314e-10, 1.9713879988437584e-7, -3.551889476633975e-5, 0.004248542489215875, -0.3277965440319509, 8.128293001726805]
    end # struct Titan

    @inline function _furnsh_required(spice_path::String, relpath::String)
        kernel_path = joinpath(spice_path, relpath)
        isfile(kernel_path) || throw(ArgumentError("Required SPICE kernel not found: $kernel_path"))
        furnsh(kernel_path)
        return kernel_path
    end

    function _furnsh_first_existing(spice_path::String, relpaths::NTuple{N, String}) where {N}
        for relpath in relpaths
            kernel_path = joinpath(spice_path, relpath)
            if isfile(kernel_path)
                furnsh(kernel_path)
                return kernel_path
            end
        end
        throw(ArgumentError("Unable to find required SPICE kernel in $(spice_path). Tried: $(join(relpaths, ", "))"))
    end

    # Constructors
    function Earth(topo_harmonics_file::String, spice_path::String="GRAM Suite 2.0/SPICE")
        earth = Earth()
        # TopographyHarmonicsWorkspace!(topo_harmonics_file, earth)
        _furnsh_required(spice_path, "pck/pck00011.tpc")
        _furnsh_required(spice_path, "lsk/naif0012.tls")
        _furnsh_first_existing(spice_path, ("spk/planets/de440s.bsp", "spk/planets/de440_GRAM.bsp"))
        # calc_J2000_to_pci_rotation_matrix!(earth.α, earth.δ, earth)
        return earth
    end

    function Mars(topo_harmonics_file::String, spice_path::String="GRAM Suite 2.0/SPICE")
        mars = Mars()
        # TopographyHarmonicsWorkspace!(topo_harmonics_file, mars)
        _furnsh_required(spice_path, "pck/pck00011.tpc")
        _furnsh_required(spice_path, "lsk/naif0012.tls")
        _furnsh_first_existing(spice_path, ("spk/planets/de440s.bsp", "spk/planets/de440_GRAM.bsp"))
        # calc_J2000_to_pci_rotation_matrix!(mars.α, mars.δ, mars)
        return mars
    end

    function Venus(topo_harmonics_file::String, spice_path::String="GRAM Suite 2.0/SPICE")
        venus = Venus()
        # TopographyHarmonicsWorkspace!(topo_harmonics_file, venus)
        _furnsh_required(spice_path, "pck/pck00011.tpc")
        _furnsh_required(spice_path, "lsk/naif0012.tls")
        _furnsh_first_existing(spice_path, ("spk/planets/de440s.bsp", "spk/planets/de440_GRAM.bsp"))
        # calc_J2000_to_pci_rotation_matrix!(venus.α, venus.δ, venus)
        return venus
    end

    function Titan(topo_harmonics_file::String, spice_path::String="GRAM Suite 2.0/SPICE")
        titan = Titan()
        # TopographyHarmonicsWorkspace!(topo_harmonics_file, titan)
        _furnsh_required(spice_path, "pck/pck00011.tpc")
        _furnsh_required(spice_path, "lsk/naif0012.tls")
        _furnsh_first_existing(spice_path, ("spk/planets/de440s.bsp", "spk/planets/de440_GRAM.bsp"))
        _furnsh_first_existing(spice_path, ("spk/satellites/sat441.bsp", "spk/satellites/sat441_GRAM.bsp"))
        # calc_J2000_to_pci_rotation_matrix!(titan.α, titan.δ, titan)
        return titan
    end

    # Helper functions
    
    function TopographyHarmonicsWorkspace!(topo_harmonics_file::String, planet::P) where P <: AbstractPlanet
        if topo_harmonics_file != ""
            Clm_topo, Slm_topo = read_topography_harmonics(topo_harmonics_file)
            # planet.topography_workspace = TopographyHarmonicsWorkspace()
            planet.topography_workspace.Clm = Clm_topo
            planet.topography_workspace.Slm = Slm_topo
            # Precompute the fn, gnm, hnm values for the topography function workspace
            n = 1:size(Clm_topo, 1)-1
            N = size(Clm_topo, 1)-1
            M = size(Clm_topo, 2)-1
            planet.topography_workspace.fn = @. √( ((1 + (n == 1))*(2.0*n + 1)) / (2.0*n) )
            # Preallocate matrices with 'undef' for maximum speed
            planet.topography_workspace.G = Matrix{Float64}(undef, N + 1, M + 1)
            planet.topography_workspace.H = Matrix{Float64}(undef, N + 1, M + 1)
            
            # Initialize with zeros to handle unused indices safely
            planet.topography_workspace.G .= 0.0
            planet.topography_workspace.H .= 0.0

            @inbounds for j = 1:M+1
                m = j - 1
                @inbounds for i = j+1:N+1
                    n = i - 1
                    
                    # Common denominator factor
                    den = (n + m) * (n - m)
                    
                    # Precompute Gnm
                    planet.topography_workspace.G[i, j] = sqrt(((2n + 1) * (2n - 1)) / den)

                    # Precompute Hnm (only if i > j + 1, which means n > m + 1)
                    if i > j + 1
                        num_h = (2n + 1) * (n - m - 1) * (n + m - 1)
                        den_h = (2n - 3) * den
                        planet.topography_workspace.H[i, j] = sqrt(num_h / den_h)
                    end
                end
            end

            planet.topography_workspace.Plm = Matrix{Float64}(undef, N + 1, M + 1) # Preallocate the matrix for the Associated Legendre Polynomial evaluations
        end
    end

    read_topography_harmonics(filepath::String)::Tuple{Matrix{Float64}, Matrix{Float64}} = begin
        harmonics_data = CSV.File(filepath)
        
        # Pre-initialize the Clm and Slm arrays
        total_data_size = size(harmonics_data, 1)
        degree = maximum(harmonics_data.degree[end]) + 1

        # p_class.A_topo = zeros(degree+1, degree+1) # Preallocate the matrix for the Associated Legendre Polynomial evaluations
        Clm_topo = zeros(Float64, degree, degree)
        Slm_topo = zeros(Float64, degree, degree)

        # Read in all the data from the DataFrame
        for i=1:total_data_size
            l = harmonics_data.degree[i] + 1 # Get the degree, l, from the data and convert to an index (subtract 1 because the data starts at 2nd degree coefficient)
            m = harmonics_data.order[i] + 1 # Get the order, m, from the data and convert to an index (add 1 because the data starts at 0th order coefficient)
            Clm_topo[l, m] = harmonics_data.C[i]
            Slm_topo[l, m] = harmonics_data.S[i]
        end
        return Clm_topo, Slm_topo
    end

    calc_J2000_to_pci_rotation_matrix!(α::Float64, δ::Float64, planet::P) where P <: AbstractPlanet = begin
        σ1 = sqrt(cos(δ)^4 + cos(δ)^2*sin(δ)^2)
        planet.J2000_to_pci .= SMatrix{3, 3, Float64}([-sin(α) cos(α) 0;
                        -cos(δ)*cos(α)*sin(δ)/σ1 -cos(δ)*sin(α)*sin(δ)/σ1 cos(δ)^2/σ1;
                        cos(δ)*cos(α) cos(δ)*sin(α) sin(δ)])
    end
   
end # module Planets
