include("../physical_models/MonteCarlo_pertrubations.jl")
include("../physical_models/Planet_data.jl")
include("../physical_models/Mission.jl")
include("../utils/Save_csv.jl")
include("../utils/Plot_data.jl")
include("Aerobraking.jl")
include("../utils/Reference_system.jl")

using SPICE
using StaticArrays
using AstroTime
using PythonCall
using Arrow

function aerobraking_campaign(args, state)
    save_res = args[:results]
    config.cnf.Gram_directory = args[:directory_Gram]

    # Descent towards Mars
    purpose = "Aerobraking around Mars"

    mission = Dict(:Purpose => purpose,
                   :Planet => args[:planet],
                   :Gravity_Model => args[:gravity_model], 
                   :Density_Model => args[:density_model], 
                   :Wind => args[:wind],
                   :Aerodynamic_Model => args[:aerodynamic_model],
                   :Thermal_Model => args[:thermal_model],
                   :Control => args[:control_mode],
                   :Firings => args[:thrust_control],
                   :Shape => args[:body_shape],
                   :Monte_Carlo => args[:montecarlo])

    if args[:print_res] == true
        println("Mission is: ", mission)
    end

    ip = mission_def(mission)
    p_class = planet_data(ip.M.planet)
    config.model.planet = p_class

    # Load SPICE kernels if required
    furnsh(args[:directory_Spice] * "/pck/pck00011.tpc")
    # furnsh(args[:directory_Spice] * "/spk/planets/de440_GRAM.bsp")
    furnsh(args[:directory_Spice] * "/lsk/naif0012.tls")
    # furnsh(args[:directory_Spice] * "/spk/planets/de440s.bsp")
    furnsh(args[:directory_Spice] * "/spk/planets/de430.bsp")
    furnsh(args[:directory_Spice] * "/spk/satellites/sat441_GRAM.bsp")
    furnsh(args[:directory_Spice] * "/spk/satellites/mar097_GRAM.bsp")
    
    # If using lat/lon initial conditions, correct the initial orbital elements
    if args[:orientation_type] == 1
        # Get the latitude and longitude of the initial conditions
        lat = args[:latitude]
        lon = args[:longitude]

        # Convert latitude and longitude to radians
        lat_rad = deg2rad(lat)
        lon_rad = deg2rad(lon)
        α_rad = deg2rad(args[:azimuth])
        γ_rad = deg2rad(args[:γ_initial_a])
        config.cnf.et = utc2et(to_utc(DateTime(args[:year], args[:month], args[:day], args[:hours], args[:minutes], args[:secs])))
        p_class.L_PI .= SMatrix{3, 3, Float64}(pxform("J2000", "IAU_" * uppercase(p_class.name), config.cnf.et))*p_class.J2000_to_pci'
        # will have to rethink this to use the gamma/v step initial conditions
        OE = latlongtoOE([lat_rad, lon_rad, args[:EI]*1e3], p_class, γ_rad, α_rad, args[:v_initial_a])
 
        OE[3:5] = rad2deg.(OE[3:5])
        args[:inclination] = OE[3]
        args[:Ω] = OE[4]
        args[:ω] = OE[5]
        state[:Inclination] = OE[3]
        state[:Ω] = OE[4]
        state[:ω] = OE[5]
    end
    # Set up n-body gravity
    if length(args[:n_bodies]) != 0
        for i=1:length(args[:n_bodies])
            push!(config.cnf.n_bodies_list, planet_data(args[:n_bodies][i]))
        end
    end

    # Set up spherical harmonics coefficients
    if args[:gravity_harmonics] == true
        # Read in the gravity harmonics data
        harmonics_data = CSV.read(args[:gravity_harmonics_file], DataFrame)
        
        # Pre-initialize the Clm and Slm arrays
        total_data_size = size(harmonics_data, 1)
        degree = maximum(harmonics_data[:, 1]) + 1

        p_class.A_grav = zeros(degree+1, degree+1) # Preallocate the matrix for the Associated Legendre Polynomial evaluations
        p_class.Clm = zeros(degree, degree)
        p_class.Slm = zeros(degree, degree)

        # Read in all the data from the DataFrame
        for i=1:total_data_size
            l = harmonics_data[i, 1] + 1 # Get the degree, l, from the data and convert to an index (subtract 1 because the data starts at 2nd degree coefficient)
            m = harmonics_data[i, 2] + 1 # Get the order, m, from the data and convert to an index (add 1 because the data starts at 0th order coefficient)
            p_class.Clm[l, m] = harmonics_data[i, 3]
            p_class.Slm[l, m] = harmonics_data[i, 4]
        end

        # Precalculate N1, N2
        L = args[:L]
        M = args[:M]
        N1 = zeros(L+4, L+4)
        N2 = zeros(L+4, L+4)
        VR01 = zeros(L+1, L+1)
        VR11 = zeros(L+1, L+1)
        sqrt_2 = sqrt(2)
        for m = 0:M+2
            j = m + 1
            for l = m+2:L+2
                i = l + 1
                N1[i, j] = √((2*l+1)*(2*l-1)/(l+m)/(l-m))
                N2[i, j] = √((l+m-1)*(2*l+1)*(l-m-1)/(2*l-3)/(l+m)/(l-m))
            end
        end

        for l = 0:L
            i = l + 1
            for m = 0:min(M, l)
                j = m + 1
                divisor = m == 0 ? sqrt_2 : 1
                VR01[i, j] = sqrt((l-m)*(l+m+1)) / divisor
                VR11[i, j] = sqrt((2*l+1)*(l+m+2)*(l+m+1)/(2*l+3)) / divisor
            end
        end
        p_class.N1 = N1
        p_class.N2 = N2
        p_class.VR01 = VR01
        p_class.VR11 = VR11

        A = zeros(L+4, L+4)
        R = zeros(L+4)
        I = zeros(L+4)
        A[1, 1] = 1
        # Fill the diagonal elements of A
        for l = 1:L+2
            i = l + 1
            A[i, i] = sqrt((2*l+1)/(2*l))*A[i-1, i-1]
        end
        p_class.Re = R
        p_class.Im = I
        p_class.A = A

    end

    # Set up the planet shape
    if args[:topography_model] == "Spherical Harmonics"
        harmonics_data = CSV.read(args[:topography_harmonics_file], DataFrame)
        
        # Pre-initialize the Clm and Slm arrays
        total_data_size = size(harmonics_data, 1)
        degree = maximum(harmonics_data[:, 1]) + 1

        p_class.A_topo = zeros(degree+1, degree+1) # Preallocate the matrix for the Associated Legendre Polynomial evaluations
        p_class.Clm_topo = zeros(degree, degree)
        p_class.Slm_topo = zeros(degree, degree)

        # Read in all the data from the DataFrame
        for i=1:total_data_size
            l = harmonics_data[i, 1] + 1 # Get the degree, l, from the data and convert to an index (subtract 1 because the data starts at 2nd degree coefficient)
            m = harmonics_data[i, 2] + 1 # Get the order, m, from the data and convert to an index (add 1 because the data starts at 0th order coefficient)
            p_class.Clm_topo[l, m] = harmonics_data[i, 3]
            p_class.Slm_topo[l, m] = harmonics_data[i, 4]
        end

    end

    if args[:gravity_model] == "Inverse Squared"
        p_class.Rp_p = p_class.Rp_e
    end

    # Vehicle - calculation notebook page1
    # Mass
    # dry_mass = args[:spacecraft_model].dry_mass
    # prop_mass = args[:spacecraft_model].prop_mass 
    # mass = dry_mass + prop_mass
    mass = config.get_spacecraft_mass(args[:spacecraft_model])


    # Spacecraft Shape
    # if args[:body_shape] == "Spacecraft"
    #     area_body = args[:length_sp] * args[:height_sp]   # 33.38 # 7.26# This is recalculated for the new sc config. 11 (look notes)# m^2 2001 Mars Odyssey Aerobraking, Smith & Bell paper
    #     length_sp = args[:length_sp]                      # 11.4 # 3.7617#5.7 # m solar array length https://www.jpl.nasa.gov/news/press_kits/odysseyarrival.pdf
    #     height_sp = area_body / length_sp

    #     # Main Body
    #     length_ody = args[:length_sat]   # m 
    #     height_ody = args[:height_sat]   # m
    #     width_ody = args[:width_sat]     # m

    # # Blunted Body Shape
    # elseif args[:body_shape] == "Blunted Cone"
    #     δ = args[:cone_angle] # deg
    #     nose_radius = args[:nose_radius] # 0.6638 # m
    #     base_radius = args[:base_radius] # 2.65/2 # m
    # end

    apoapsis = state[:Apoapsis]

    state[:Periapsis] = p_class.Rp_e + state[:Periapsis]*1e3
    state[:vi] = (lowercase(args[:type_of_mission]) == "time" && args[:ν] >= 0.0) ? deg2rad(args[:ν]) : deg2rad(180.0001)

    if args[:montecarlo] == true
        state = monte_carlo_initial_condition(state, args)
    end

    semimajoraxis_in = (state[:Apoapsis] + state[:Periapsis])/2
    eccentricity_in = (state[:Apoapsis] - state[:Periapsis]) / (state[:Apoapsis] + state[:Periapsis])
    apoapsis = state[:Apoapsis]
    periapsis = state[:Periapsis]

    # Initial Condition
    if args[:drag_passage] == true
        h_0 = args[:EI] * 1e3
    elseif args[:body_shape] == "Blunted Cone"
        h_0 = args[:EI] * 1e3
        args[:AE] = h_0/1e3
        args[:EI] = h_0/1e3
    end

    if Bool(args[:drag_passage]) || args[:body_shape] == "Blunted Cone"
        r = p_class.Rp_e + h_0
        
        state[:vi] = -acos(1 / eccentricity_in * (semimajoraxis_in * (1 - eccentricity_in^2) / r - 1))
        
        if args[:montecarlo] == true
            state = monte_carlo_true_anomaly(state, args)
            apoapsis = state[:Apoapsis]
            periapsis = state[:Periapsis]
        end
    end

    # Initial Model Definition
    # Body
    if args[:body_shape] == "Spacecraft"
        b_class = args[:spacecraft_model]
        println("Area: " * string(config.get_spacecraft_reference_area(b_class)) * " m^2")

    elseif args[:body_shape] == "Blunted Cone" # TODO: Change this to new spacecraft model method

        Mass = mass
        Delta = δ
        NoseRadius = nose_radius
        BaseRadius = base_radius
        Area_tot = pi * BaseRadius^2

        b_class = config.Body(Mass, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Area_tot, Delta, NoseRadius, BaseRadius)

    end

    function initialconditions()
        a = semimajoraxis_in
        e = eccentricity_in
        i = deg2rad(state[:Inclination])
        Ω = deg2rad(state[:Ω])
        ω = deg2rad(state[:ω])
        vi = state[:vi]
        m0 = mass
        year = args[:year]
        month = args[:month]
        day = args[:day]
        hour = args[:hours]
        min = args[:minutes]
        second = args[:secs]
        time_rot = args[:planettime]
        el_time = 0.0 # Elapsed time in seconds
        DateTimeIC = from_utc(DateTime(year, month, day, hour, min, second))
        DateTimeJ2000 = from_utc(DateTime(2000, 1, 1, 12, 0, 0))
        ic = config.Initial_condition(a, e, i, Ω, ω, vi, m0, year, month, day, hour, min, second, time_rot, el_time, DateTimeIC, DateTimeJ2000)

        return ic
    end

    ic_class = initialconditions()

    function aerodynamics()
        δ = deg2rad(0)
        α = deg2rad(args[:α])
        thermal_accomodation_factor = args[:thermal_accomodation_factor]
        reflection_coefficient = args[:reflection_coefficient]
        thermal_contact = 0
        heat_rate_limit = args[:max_heat_rate]
        heat_load_limit = args[:max_heat_load]

        a = config.Aerodynamics(δ, α, thermal_accomodation_factor, reflection_coefficient, thermal_contact, heat_rate_limit, heat_load_limit)

        return a        
    end

    a_class = aerodynamics()

    # Engine
    args[:phi] = deg2rad(args[:phi])
    function engine()
        ϕ = args[:phi]
        g_e = 9.81
        T = args[:thrust]
        Isp = 200

        e = config.Engines(ϕ, g_e, T, Isp)

        return e
    end

    e_class = engine()

    function model()
        body = b_class
        planet = p_class
        initialcondition = ic_class
        aerodynamics = a_class
        engine = e_class

        m = config.Model(body, planet, aerodynamics, engine, initialcondition)

        return m
    end

    m = model()

    # Define gram atmosphere
    gram = nothing
    gram_atmosphere = nothing
    if uppercase(args[:density_model]) == "GRAM"
        sys = pyimport("sys")
        os = pyimport("os")
        if !(args[:directory_Gram] in pyconvert(Vector{String}, sys.path))
            sys.path.append(args[:directory_Gram])
        end
        gram = pyimport("gram")
        inputParameters = Dict("earth" => gram.EarthInputParameters(),
                            "mars" => gram.MarsInputParameters(),
                            "venus" => gram.VenusInputParameters(),
                            "titan" => gram.TitanInputParameters())
        
        namelistReaders = Dict("earth" => gram.EarthNamelistReader(),
                            "mars" => gram.MarsNamelistReader(),
                            "venus" => gram.VenusNamelistReader(),
                            "titan" => gram.TitanNamelistReader())
            
        atmospheres = Dict("earth" => gram.EarthAtmosphere(),
                        "mars" => gram.MarsAtmosphere(),
                        "venus" => gram.VenusAtmosphere(),
                        "titan" => gram.TitanAtmosphere())

        planet_name = m.planet.name
        input_parameters = inputParameters[planet_name]

        # Mars has some weird specific parameters, so this line is just to check to make sure the it doesn't do it for the other planets
        if planet_name == "mars"
            # input_parameters.dataPath = os.path.join(os.path.dirname(os.path.abspath(@__FILE__)),"..", "GRAM_Data", "Mars", "data", "")
            input_parameters.dataPath = args[:directory_Gram_data] * "/Mars/data/"
            if !Bool(os.path.exists(input_parameters.dataPath))
                throw(ArgumentError("GRAM data path not found: " * input_parameters.dataPath))
            end
        end

        if planet_name == "earth"
            # input_parameters.dataPath = os.path.join(os.path.dirname(os.path.abspath(@__FILE__)),"..", "GRAM_Data", "Mars", "data", "")
            input_parameters.dataPath = args[:directory_Gram_data] * "/Earth/data/"
            if !Bool(os.path.exists(input_parameters.dataPath))
                throw(ArgumentError("GRAM data path not found: " * input_parameters.dataPath))
            end
        end

        reader = namelistReaders[planet_name]
        reader.tryGetSpicePath(input_parameters)

        gram_atmosphere = atmospheres[planet_name]
        gram_atmosphere.setInputParameters(input_parameters)
        
        if planet_name == "earth"
            gram_atmosphere.setMERRA2Parameters(0, -90.0, 90.0, 0.0, 359.99999)
        end
        gram_atmosphere.setPerturbationScales(1.5)
        gram_atmosphere.setMinRelativeStepSize(0.5)
        if args[:montecarlo] == 1
            gram_atmosphere.setSeed(Int(round(rand()*10000)))
        else
            gram_atmosphere.setSeed(1001)
        end

        if planet_name == "mars"
            gram_atmosphere.setMOLAHeights(false)
        end

        ttime = gram.GramTime()
        ttime.setStartTime(args[:year], args[:month], args[:day], args[:hours], args[:minutes], args[:secs], gram.UTC, gram.PET)
        gram_atmosphere.setStartTime(ttime)
    end
    # Initialization - Reset all the config index for new simulation
    config.cnf.count_aerobraking = 0
    config.cnf.count_overcome_hr = 0
    config.cnf.save_index_heat = 0
    config.cnf.index_propellant_mass = 1
    config.cnf.counter_random = 0
    config.cnf.DU = Bool(args[:normalize]) ? semimajoraxis_in : 1
    config.cnf.TU = Bool(args[:normalize]) ? sqrt(config.cnf.DU^3 / m.planet.μ) : 1
    config.cnf.MU = Bool(args[:normalize]) ? mass : 1

    # Save results
    if save_res == 1
        if args[:filename] == 1
            if args[:montecarlo] == true
                folder_name = args[:simulation_filename][1:findfirst!(args[:simulation_filename], "_nMC)")]
            else
                folder_name = args[:simulation_filename]
            end
            
            name = args[:directory_results] * "/" * folder_name

            if !isdir(args[:directory_results])
                mkpath(args[:directory_results])
            end
            
            filename = name * ".csv"
        else
            name = args[:directory_results] * "/" * "GRAMver_" * string(args[:Gram_version])

            if !isdir(args[:directory_results])
                mkpath(args[:directory_results])
            end

            filename = name * ".csv"
        end
        # if the file already exists, clear the current data
        if filesize(filename) > 0
            file = open(filename, "w")
            truncate(file, 0)
        end
        # save_csv(filename, args)
    end

    # Initialize the arrow writer for plotting
    arrow_writer = nothing
    temp_name = nothing
    if args[:plot] == true
        temp_name = tempname()
        println("Temporary directory created for plotting: " * temp_name)
        arrow_writer = open(Arrow.Writer, temp_name)
    end

    ##########################################################
    # RUN SIMULATION
    config.cnf.heat_rate_limit = args[:max_heat_rate]
    t_el = @elapsed begin
        aerobraking(ip, m, args, gram, gram_atmosphere, filename, temp_name)
    end
    ##########################################################

    # Finalize the arrow writer if plotting is enabled
    if args[:plot] == true
        close(arrow_writer)
        println("Arrow writer closed. Data saved to: " * temp_name)
    end
    # Print final results
    if Bool(args[:print_res])
        println("ρ: " * string(maximum(config.solution.physical_properties.ρ)) * " kg/m^3")
        println("heat rate: " * string(maximum(maximum.(config.solution.performance.heat_rate))) * " W/cm^2")
    end

    

    if Bool(args[:print_res])
        println("Elapsed time: " * string(t_el) * " s")
    end

    if args[:plot] == true
        plots(state, m, name, args, temp_name)
    end

    # rm(temp_name, recursive=true, force=true) # Remove the temporary directory used for plotting
end