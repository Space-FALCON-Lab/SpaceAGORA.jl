module Planets
    include(joinpath(@__DIR__, "..", "..", "environment", "ephemerides", "planet_shapes.jl"))
    using ..AbstractTypes: AbstractPlanet
    using StaticArrays
    using CSV
    using SPICE
    export Earth, Mars, Venus, Moon, Titan
    const SPICE_LOCK = parentmodule(parentmodule(@__MODULE__)).RuntimeServices.SPICE_LOCK
    # Every lock site in this file is planet construction, kernel furnishing or a
    # body-constant lookup -- setup, not solve. They are attributed so that a
    # solve's occupancy is not silently missing the setup that preceded it, and
    # so that a run whose "lock cost" is really kernel loading reads as such.
    const _SPICE_BODY_LOCK =
        parentmodule(parentmodule(@__MODULE__)).RuntimeServices.tracked_lock(:spice_body)
    const MARS_MU_M3S2 = 0.4282837285418775e5 * 1e9

    # Planet constructors (Earth(), Mars(), ...) are called once per case setup, but
    # some callers (e.g. Monte Carlo sample loops that build a fresh mission config
    # per trial) call them many times per process. `furnsh` has no cheap "already
    # loaded" check of its own and CSPICE's kernel table is a fixed-size global
    # resource with no automatic eviction, so reloading the same kernel path
    # thousands of times in one process without ever unloading eventually exhausts
    # it ("Insufficient dynamic kernel table space" / "too many kernels loaded").
    # Track furnished paths here so repeat construction is a cheap no-op instead of
    # a repeat furnsh call.
    const _FURNISHED_KERNELS = Set{String}()

    # Two things can desynchronise this cache from CSPICE's actual pool.
    #
    # `kclear()` wipes CSPICE's own kernel pool but has no way to know about this
    # cache (or the planet-instance caches below), so anything that calls `kclear()`
    # must call this too — otherwise a kernel already in `_FURNISHED_KERNELS` gets
    # skipped on the next furnish even though `kclear()` just unloaded it from CSPICE,
    # leaving lookups (e.g. utc2et needing the leapseconds kernel) failing against an
    # empty pool. Also drops the cached Earth()/Mars()/etc. instances themselves
    # (defined further down, near their constructors) since those were built from
    # kernels that no longer exist post-kclear.
    #
    # PRECOMPILATION is the second route, and it is worse because it is silent
    # and permanent. `_FURNISHED_KERNELS` is a module-level `const`, so a
    # precompile workload that constructs a SPICE-backed planet bakes the
    # furnished paths into the pkgimage; every later process then starts with a
    # cache claiming kernels are loaded and a CSPICE pool that is empty, and
    # nothing furnishes them again. `src/precompile_workload.jl` calls this at
    # the end of `@setup_workload` for exactly that reason.
    @inline function _reset_furnished_kernels!()
        lock(_SPICE_BODY_LOCK) do
            empty!(_FURNISHED_KERNELS)
            empty!(_EARTH_CACHE)
            empty!(_MARS_CACHE)
            empty!(_VENUS_CACHE)
            empty!(_TITAN_CACHE)
            empty!(_MOON_CACHE)
        end
        return nothing
    end

    @inline function _furnsh_once(kernel_path::String)
        resolved = abspath(kernel_path)
        lock(_SPICE_BODY_LOCK) do
            resolved in _FURNISHED_KERNELS && return nothing
            furnsh(resolved)
            push!(_FURNISHED_KERNELS, resolved)
            return nothing
        end
        return nothing
    end
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
        Rp_e::Float64 = 6.3781366e6 # Equatorial radius in meters (SPICE pck00011)
        Rp_p::Float64 = 6.3567519e6 # Polar radius in meters (SPICE pck00011)
        Rp_m::Float64 = 6.371008366666667e6  # Mean radius in meters (SPICE pck00011)
        mass::Float64  = 5.972e24 # Mass in kg
        p::Float64 = 101325.0 # Surface pressure in Pascals
        k::Float64 = 1.83e-4 # Chapman heating coefficient, kg^0.5/m
        ω::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 7.2921066e-5) # Angular velocity vector in rad/s
        # μ::Float64 = 3.986004415e14 # Standard gravitational parameter in m^3/s^2, GMAT default value
        μ::Float64 = 3.98600436233e14 # Standard gravitational parameter in m^3/s^2, DE421 value ends with 233 after 436
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
        L_PI::MMatrix{3, 3, Float64} = MMatrix{3, 3, Float64}(zeros(3, 3)) # Rotation matrix from J2000 to planet-fixed frame, function of time
        topography_workspace::TopographyHarmonicsWorkspace = TopographyHarmonicsWorkspace() # Workspace for calculating topography harmonics
        polyfit_coeffs::Vector{Float64} = [-1.7539409645214832e-57, 2.735656076315809e-53, -1.8243490769488347e-49, 6.504765617793163e-46, -1.1637408657034938e-42, 8.044884138893168e-41, 4.264962263039017e-36, -7.651115834387683e-33, -3.188248308052816e-30, 3.8370830656820503e-26, -8.557502178008995e-23, 1.137879849173412e-19, -1.0408232216096158e-16, 6.834085016894604e-14, -3.2506596548183e-11, 1.1089006707870246e-08, -2.639423772958483e-06, 0.0004165844083994442, -0.03967261693733797, 1.8349343859319074, -38.14918904018883]
        topography_function::Function = Earth_elevation! # Function to calculate elevation based on topography
    end # struct Earth

    @kwdef struct Mars <: AbstractPlanet
        name::String = "Mars" # Name of the planet
        Rp_e::Float64 = 3.396190e6 # Equatorial radius in meters
        Rp_p::Float64 = 3.3762e6 # Polar radius in meters
        Rp_m::Float64 = 3.3895e6 # Mean radius in meters
        mass::Float64  = 0.64171e24 # Mass in kg
        p::Float64 = 636.0 # Surface pressure in Pascals
        k::Float64 = 1.898e-4 # Chapman heating coefficient, kg^0.5/m
        ω::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 7.08823596e-5) # Angular velocity vector in rad/s
        μ::Float64 = MARS_MU_M3S2 # Standard gravitational parameter in m^3/s^2
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
        L_PI::MMatrix{3, 3, Float64} = MMatrix{3, 3, Float64}(zeros(3, 3)) # Rotation matrix from J2000 to planet-fixed frame, function of time
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
        μ::Float64 = 3.24858599e14 # Standard gravitational parameter in m^3/s^2
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
        L_PI::MMatrix{3, 3, Float64} = MMatrix{3, 3, Float64}(zeros(3, 3))
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
        L_PI::MMatrix{3, 3, Float64} = MMatrix{3, 3, Float64}(zeros(3, 3))
        topography_workspace::TopographyHarmonicsWorkspace = TopographyHarmonicsWorkspace()
        polyfit_coeffs::Vector{Float64} = [1.7989756686197253e-58, -2.7298975030491325e-54, 1.7620522402686604e-50, -6.025021166267467e-47, 1.0056316643424087e-43, 9.494104496406468e-42, -3.8472088727076255e-37, 6.051435602297366e-34, 4.074478639170247e-31, -3.244699052533356e-27, 6.66877802035039e-24, -8.360025139024445e-21, 7.301165978344981e-18, -4.650857357165472e-15, 2.1978197729328097e-12, -7.705014392936314e-10, 1.9713879988437584e-7, -3.551889476633975e-5, 0.004248542489215875, -0.3277965440319509, 8.128293001726805]
    end # struct Titan

    @kwdef struct Moon <: AbstractPlanet
        name::String = "Moon"
        Rp_e::Float64 = 1.7374e6
        Rp_p::Float64 = 1.7360e6
        Rp_m::Float64 = 1.7374e6
        mass::Float64 = 7.346e22
        p::Float64 = 0.0
        k::Float64 = 0.0
        ω::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 2.6617e-6)
        μ::Float64 = 4.902799e12
        J2::Float64 = 2.027e-4
        g_ref::Float64 = 1.62
        ρ_ref::Float64 = 0.0
        h_ref::Float64 = 0.0
        H::Float64 = 0.0
        R::Float64 = 0.0
        T_ref::Float64 = 0.0
        γ::Float64 = 0.0
        T::Float64 = 0.0
        μ_fluid::Float64 = 0.0
        Lz::Float64 = 0.0
        α::Float64 = 0.0
        δ::Float64 = 0.0
        L_PI::MMatrix{3, 3, Float64} = MMatrix{3, 3, Float64}(zeros(3, 3))
        topography_workspace::TopographyHarmonicsWorkspace = TopographyHarmonicsWorkspace()
        polyfit_coeffs::Vector{Float64} = [0.0]
    end # struct Moon

    @inline function _furnsh_required(spice_path::String, relpath::String)
        kernel_path = joinpath(spice_path, relpath)
        isfile(kernel_path) || throw(ArgumentError("Required SPICE kernel not found: $kernel_path"))
        _furnsh_once(kernel_path)
        return kernel_path
    end

    function _furnsh_first_existing(spice_path::String, relpaths::NTuple{N, String}) where {N}
        for relpath in relpaths
            kernel_path = joinpath(spice_path, relpath)
            if isfile(kernel_path)
                _furnsh_once(kernel_path)
                return kernel_path
            end
        end
        throw(ArgumentError("Unable to find required SPICE kernel in $(spice_path). Tried: $(join(relpaths, ", "))"))
    end

    function _furnsh_first_existing_if_available(spice_path::String, relpaths::NTuple{N, String}) where {N}
        for relpath in relpaths
            kernel_path = joinpath(spice_path, relpath)
            if isfile(kernel_path)
                _furnsh_once(kernel_path)
                return kernel_path
            end
        end
        return nothing
    end

    @inline function _spice_body_pool_name(planet_name::String)::String
        return planet_name == "Moon" ? "MOON" : uppercase(planet_name)
    end

    @inline function _spice_body_radii_m(planet_name::String)::NTuple{3, Float64}
        radii_km = lock(_SPICE_BODY_LOCK) do
            bodvrd(_spice_body_pool_name(planet_name), "RADII")
        end
        return (radii_km[1] * 1e3, radii_km[2] * 1e3, radii_km[3] * 1e3)
    end

    @inline function _spice_body_gm_m3s2(planet_name::String)::Float64
        gm_km3s2 = lock(_SPICE_BODY_LOCK) do
            bodvrd(_spice_body_pool_name(planet_name), "GM")
        end
        return gm_km3s2[1] * 1e9
    end

    function _gravity_constants_kernel_if_available(spice_path::String)
        for relpath in (
            "pck/de_403_masses.tpc",
            # "spk/planets/de_403_masses.tpc",
            "pck/gm_de440.tpc",
            "pck/gm_de441.tpc",
            "pck/gm_de431.tpc",
            "pck/gm_de430.tpc"
        )
            kernel_path = joinpath(spice_path, relpath)
            if isfile(kernel_path)
                _furnsh_once(kernel_path)
                return kernel_path
            end
        end
        return nothing
    end

    function _spice_backed_planet_kwargs(planet_name::String)
        rp_e, rp_p_2, rp_p = _spice_body_radii_m(planet_name)
        kwargs = Dict{Symbol, Any}(
            :Rp_e => rp_e,
            :Rp_p => rp_p,
            :Rp_m => (rp_e + rp_p_2 + rp_p) / 3.0
        )
        if planet_name == "Mars"
            kwargs[:μ] = MARS_MU_M3S2
        else
            try
                kwargs[:μ] = _spice_body_gm_m3s2(planet_name)
            catch
                # GM constants require a dedicated kernel in addition to the standard PCK.
            end
        end
        return kwargs
    end

    @inline function _planetary_kernel_override_relpath()::String
        return strip(get(ENV, "SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH", ""))
    end

    function _furnsh_planetary_kernel(spice_path::String)
        override_relpath = _planetary_kernel_override_relpath()
        if !isempty(override_relpath)
            return _furnsh_required(spice_path, override_relpath)
        end
        return _furnsh_first_existing(
            spice_path,
            (
                "spk/planets/de430.bsp",
                "spk/planets/de421.bsp",
                "spk/planets/de442s.bsp",
                "spk/planets/de442.bsp",
                "spk/planets/de440s.bsp",
                "spk/planets/de440_GRAM.bsp"
            )
        )
    end

    @inline function _furnsh_mars_system_kernel(spice_path::String)
        return _furnsh_first_existing(
            spice_path,
            (
                "spk/satellites/mar097_GRAM.bsp",
                "spk/satellites/mar097.bsp",
            )
        )
    end

    function _furnsh_mars_pck(spice_path::String)
        # Always try to load a modern generic text PCK first so the standard
        # Mars frame constants are available from the starter-pack bundle.
        for relpath in ("pck/pck00011.tpc", "pck/pck00010.tpc")
            kernel_path = joinpath(spice_path, relpath)
            if isfile(kernel_path)
                _furnsh_once(kernel_path)
                return kernel_path
            end
        end

        # Backward-compatible fallback for older SPICE bundles that only ship the
        # pck00008 Mars data plus the local quadratic M2 correction patch.
        _furnsh_required(spice_path, "pck/pck00008.tpc")
        return _furnsh_required(spice_path, "pck/mars_iau2000_m2_quadratic_patch.tpc")
    end

    @inline function _spice_body_fixed_frame(planet_name::String)::String
        return planet_name == "Moon" ? "MOON_PA_DE421" : lock(_SPICE_BODY_LOCK) do
            _, resolved_frame = cnmfrm(planet_name)
            resolved_frame
        end
    end

    # Constructors
    #
    # `Earth("", spice_path)` / `Mars(...)` / etc. are called pervasively across the
    # codebase (case builders, Monte Carlo sample loops, test setup) with the same
    # `(topo_harmonics_file, spice_path)` pair every time within a process. The
    # `topo_harmonics_file` argument is currently inert everywhere: `TopographyHarmonicsWorkspace!`
    # (below) is the only thing that ever mutates a planet after construction, and no
    # caller in this codebase invokes it (every call site passes `topo_harmonics_file=""`
    # and the call is commented out in each constructor). That makes these constructors
    # pure for any given key, so — on top of `_furnsh_once` making repeat kernel loads a
    # no-op — cache the constructed planet itself: repeat calls become a `Dict` lookup
    # instead of a `SPICE_LOCK`-guarded `bodvrd` FFI round trip (Mars/Venus/Titan/Moon) or
    # kernel-existence probing (all five). This matters most for Monte Carlo sample loops,
    # which build a fresh mission config (and therefore a fresh planet) per trial, all
    # serialized on the same global lock when run multi-threaded.
    #
    # NOTE: this assumes the returned planet is never mutated after construction. If
    # `TopographyHarmonicsWorkspace!` is ever wired back up for a real caller, this cache
    # needs revisiting — every consumer of a given key would be sharing one mutable
    # workspace.
    const _EARTH_CACHE = Dict{Tuple{String, String}, Earth}()
    const _MARS_CACHE  = Dict{Tuple{String, String}, Mars}()
    const _VENUS_CACHE = Dict{Tuple{String, String}, Venus}()
    const _TITAN_CACHE = Dict{Tuple{String, String}, Titan}()
    const _MOON_CACHE  = Dict{Tuple{String, String}, Moon}()

    function Earth(topo_harmonics_file::String, spice_path::String="data/GRAMSuite.jl/GRAM Suite 2.0/SPICE")
        key = (topo_harmonics_file, spice_path)
        return lock(_SPICE_BODY_LOCK) do
            haskey(_EARTH_CACHE, key) && return _EARTH_CACHE[key]
            _furnsh_required(spice_path, "pck/pck00011.tpc")
            _furnsh_required(spice_path, "lsk/naif0012.tls")
            _furnsh_planetary_kernel(spice_path)
            _gravity_constants_kernel_if_available(spice_path)
            # The starter-pack SPICE bundle shipped in-repo may omit the high-precision
            # Earth orientation kernels. When they are absent, runtime frame transforms
            # fall back to the generic IAU_EARTH frame from pck00011.tpc.
            _furnsh_first_existing_if_available(
                spice_path,
                ("pck/earth_latest_high_prec.bpc", "pck/earth_200101_990628_predict.bpc")
            )
            _furnsh_first_existing_if_available(
                spice_path,
                (
                    "tf/earth_assoc_itrf93.tf",
                    # "fk/planets/earth_assoc_itrf93.tf",
                    # "fk/planets/earth_fixed.tf"
                )
            )
            earth = Earth()
            # TopographyHarmonicsWorkspace!(topo_harmonics_file, earth)
            _EARTH_CACHE[key] = earth
            return earth
        end
    end

    function Mars(topo_harmonics_file::String, spice_path::String="data/GRAMSuite.jl/GRAM Suite 2.0/SPICE")
        key = (topo_harmonics_file, spice_path)
        return lock(_SPICE_BODY_LOCK) do
            haskey(_MARS_CACHE, key) && return _MARS_CACHE[key]
            _furnsh_mars_pck(spice_path)
            _furnsh_required(spice_path, "lsk/naif0012.tls")
            _furnsh_planetary_kernel(spice_path)
            _furnsh_mars_system_kernel(spice_path)
            _gravity_constants_kernel_if_available(spice_path)
            mars = Mars(; _spice_backed_planet_kwargs("Mars")...)
            # TopographyHarmonicsWorkspace!(topo_harmonics_file, mars)
            _MARS_CACHE[key] = mars
            return mars
        end
    end

    function Venus(topo_harmonics_file::String, spice_path::String="data/GRAMSuite.jl/GRAM Suite 2.0/SPICE")
        key = (topo_harmonics_file, spice_path)
        return lock(_SPICE_BODY_LOCK) do
            haskey(_VENUS_CACHE, key) && return _VENUS_CACHE[key]
            _furnsh_required(spice_path, "pck/pck00011.tpc")
            _furnsh_required(spice_path, "lsk/naif0012.tls")
            _furnsh_planetary_kernel(spice_path)
            _gravity_constants_kernel_if_available(spice_path)
            venus = Venus(; _spice_backed_planet_kwargs("Venus")...)
            # TopographyHarmonicsWorkspace!(topo_harmonics_file, venus)
            _VENUS_CACHE[key] = venus
            return venus
        end
    end

    function Titan(topo_harmonics_file::String, spice_path::String="data/GRAMSuite.jl/GRAM Suite 2.0/SPICE")
        key = (topo_harmonics_file, spice_path)
        return lock(_SPICE_BODY_LOCK) do
            haskey(_TITAN_CACHE, key) && return _TITAN_CACHE[key]
            _furnsh_required(spice_path, "pck/pck00010.tpc")
            _furnsh_required(spice_path, "lsk/naif0012.tls")
            _furnsh_planetary_kernel(spice_path)
            _gravity_constants_kernel_if_available(spice_path)
            _furnsh_first_existing(spice_path, ("spk/satellites/sat441.bsp", "spk/satellites/sat441_GRAM.bsp"))
            titan = Titan(; _spice_backed_planet_kwargs("Titan")...)
            # TopographyHarmonicsWorkspace!(topo_harmonics_file, titan)
            _TITAN_CACHE[key] = titan
            return titan
        end
    end

    function Moon(topo_harmonics_file::String, spice_path::String="data/GRAMSuite.jl/GRAM Suite 2.0/SPICE")
        key = (topo_harmonics_file, spice_path)
        return lock(_SPICE_BODY_LOCK) do
            haskey(_MOON_CACHE, key) && return _MOON_CACHE[key]
            _furnsh_required(spice_path, "pck/pck00011.tpc")
            _furnsh_required(spice_path, "lsk/naif0012.tls")
            _furnsh_planetary_kernel(spice_path)
            _gravity_constants_kernel_if_available(spice_path)
            _furnsh_required(spice_path, "spk/satellites/SPICELunaCurrentKernel.bpc")
            _furnsh_required(spice_path, "tf/SPICELunaFrameKernel.tf")
            moon = Moon(; _spice_backed_planet_kwargs("Moon")...)
            _MOON_CACHE[key] = moon
            return moon
        end
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

end # module Planets
