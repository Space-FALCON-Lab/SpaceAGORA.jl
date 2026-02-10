module Planets
    include("../physical_models/Planet_shapes.jl")
    using ..AbstractTypes: Planet
    using StaticArrays
    using CSV
    export Earth, Mars
    δ(i, j) = ==(i, j)

    @kwdef mutable struct TopographyHarmonicsWorkspace
        Clm::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        Slm::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        Plm::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        fn::Vector{Float64} = [0.0, 0.0, 0.0]
        G::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
        H::Matrix{Float64} = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    end # struct HarmonicsWorkspace

    @kwdef struct Earth <: Planet
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

    @kwdef struct Mars <: Planet
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

    # Constructors
    function Earth(topo_harmonics_file::String)
        earth = Earth()
        TopographyHarmonicsWorkspace!(topo_harmonics_file, earth)
        return earth
    end

    function Mars(topo_harmonics_file::String)
        mars = Mars()
        TopographyHarmonicsWorkspace!(topo_harmonics_file, mars)
        return mars
    end

    # Helper functions
    function TopographyHarmonicsWorkspace!(topo_harmonics_file::String, planet::Planet)
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
                for i = j+1:N+1
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

    calc_J2000_to_pci_rotation_matrix!(α::Float64, δ::Float64, planet::Planet) = begin
        σ1 = sqrt(cos(δ)^4 + cos(δ)^2*sin(δ)^2)
        planet.J2000_to_pci = SMatrix{3, 3, Float64}([-sin(α) cos(α) 0;
                        -cos(δ)*cos(α)*sin(δ)/σ1 -cos(δ)*sin(α)*sin(δ)/σ1 cos(δ)^2/σ1;
                        cos(δ)*cos(α) cos(δ)*sin(α) sin(δ)])
    end
   
end # module Planets