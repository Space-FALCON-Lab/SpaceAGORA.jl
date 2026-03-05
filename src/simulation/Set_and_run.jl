using .SimulationModel

# include("../physical_models/MonteCarlo_pertrubations.jl")
# include("../physical_models/Planet_data.jl")
# include("../physical_models/Mission.jl")
# include("../utils/Save_csv.jl")
# include("../utils/Plot_data.jl")
include(joinpath(@__DIR__, "execution", "simulation_execution.jl"))
# include("../utils/Reference_system.jl")

using SPICE
using StaticArrays
using AstroTime
using Arrow

function execute_campaign(args; isolate_state::Bool=true, state=nothing)
    throw(ArgumentError("execute_campaign expects SimulationConfiguration; got $(typeof(args))."))
end

function execute_campaign(args::SimulationConfiguration; isolate_state::Bool=true, state=nothing)
    save_res = args.simulation_settings.results
    filename = ""
    arrow_filename = ""
    arrow_writer = nothing

    # mission = Dict(:Purpose => purpose,
    #                :Planet => args[:planet],
    #                :Gravity_Model => args[:gravity_model], 
    #                :Density_Model => args[:density_model], 
    #                :Wind => args[:wind],
    #                :Aerodynamic_Model => args[:aerodynamic_model],
    #                :Thermal_Model => args[:thermal_model],
    #                :Control => args[:control_mode],
    #                :Firings => args[:thrust_control],
    #                :Shape => args[:body_shape],
    #                :Monte_Carlo => args[:montecarlo])

    # if args.simulation_settings.verbose
    #     println("Mission is: ", mission)
    # end

    # Create config parameters
    # m = Model()
    # cnf = Cnf()
    # solution = Solution()

    # ip = mission_def(mission)
    # p_class = planet_data(ip.M.planet)

    # Load SPICE kernels if required
    # furnsh(args.file_paths.SPICE * "/pck/pck00011.tpc")
    # # furnsh(args.file_paths.SPICE * "/spk/planets/de440_GRAM.bsp")
    # furnsh(args.file_paths.SPICE * "/lsk/naif0012.tls")
    # furnsh(args.file_paths.SPICE * "/spk/planets/de440s.bsp")
    # # furnsh(args.file_paths.SPICE * "/spk/planets/de430.bsp")
    # furnsh(args.file_paths.SPICE * "/spk/satellites/sat441_GRAM.bsp")
    # furnsh(args.file_paths.SPICE * "/spk/satellites/mar097_GRAM.bsp")
    
    # If using lat/lon initial conditions, correct the initial orbital elements
    # if args[:orientation_type] == 1
    #     # Get the latitude and longitude of the initial conditions
    #     lat = args[:latitude]
    #     lon = args[:longitude]

    #     # Convert latitude and longitude to radians
    #     lat_rad = deg2rad(lat)
    #     lon_rad = deg2rad(lon)
    #     α_rad = deg2rad(args[:azimuth])
    #     γ_rad = deg2rad(args[:γ_initial_a])
    #     cnf.et = utc2et(to_utc(DateTime(args[:year], args[:month], args[:day], args[:hours], args[:minutes], args[:secs])))
    #     p_class.L_PI .= SMatrix{3, 3, Float64}(pxform("J2000", "IAU_" * uppercase(p_class.name), config.cnf.et))*p_class.J2000_to_pci'
    #     # will have to rethink this to use the gamma/v step initial conditions
    #     OE = latlongtoOE([lat_rad, lon_rad, args[:EI]*1e3], p_class, γ_rad, α_rad, args[:v_initial_a])
 
    #     OE[3:5] = rad2deg.(OE[3:5])
    #     args[:inclination] = OE[3]
    #     args[:Ω] = OE[4]
    #     args[:ω] = OE[5]
    #     state[:Inclination] = OE[3]
    #     state[:Ω] = OE[4]
    #     state[:ω] = OE[5]
    # end

    # Set up the planet shape (move all this initialization to planet initialization?)
    # if args.physical_models.topography
    #     harmonics_data = CSV.read(args.file_paths.topography_harmonics*"/"*args.physical_models.harmonics_data, DataFrame)
        
    #     # Pre-initialize the Clm and Slm arrays
    #     total_data_size = size(harmonics_data, 1)
    #     degree = maximum(harmonics_data[:, 1]) + 1

    #     p_class.A_topo = zeros(degree+1, degree+1) # Preallocate the matrix for the Associated Legendre Polynomial evaluations
    #     p_class.Clm_topo = zeros(degree, degree)
    #     p_class.Slm_topo = zeros(degree, degree)

    #     # Read in all the data from the DataFrame
    #     for i=1:total_data_size
    #         l = harmonics_data[i, 1] + 1 # Get the degree, l, from the data and convert to an index (subtract 1 because the data starts at 2nd degree coefficient)
    #         m = harmonics_data[i, 2] + 1 # Get the order, m, from the data and convert to an index (add 1 because the data starts at 0th order coefficient)
    #         p_class.Clm_topo[l, m] = harmonics_data[i, 3]
    #         p_class.Slm_topo[l, m] = harmonics_data[i, 4]
    #     end

    # end

    # if args[:gravity_model] == "Inverse Squared"
    #     p_class.Rp_p = p_class.Rp_e
    # end

    # Vehicle - calculation notebook page1
    # Mass
    # dry_mass = args[:spacecraft_model].dry_mass
    # prop_mass = args[:spacecraft_model].prop_mass 
    # mass = dry_mass + prop_mass
    # mass = get_spacecraft_mass(args[:spacecraft_model])

    # apoapsis = state[:Apoapsis]

    # state[:Periapsis] = p_class.Rp_e + state[:Periapsis]*1e3
    # state[:vi] = (lowercase(args[:type_of_mission]) == "time" && args[:ν] >= 0.0) ? deg2rad(args[:ν]) : deg2rad(180.0001)

    # if args[:montecarlo] == true
    #     state = monte_carlo_initial_condition(state, args)
    # end

    # semimajoraxis_in = (state[:Apoapsis] + state[:Periapsis])/2
    # eccentricity_in = (state[:Apoapsis] - state[:Periapsis]) / (state[:Apoapsis] + state[:Periapsis])
    # apoapsis = state[:Apoapsis]
    # periapsis = state[:Periapsis]

    # Initial Condition
    # if args[:drag_passage] == true
    #     h_0 = args[:EI] * 1e3
    # elseif args[:body_shape] == "Blunted Cone"
    #     h_0 = args[:EI] * 1e3
    #     args[:AE] = h_0/1e3
    #     args[:EI] = h_0/1e3
    # end

    # TODO: Make this part of the initial condition definition in config file, e.g., constructor for drag passage/entry specifically to set initial true anomaly and radius
    # if Bool(args[:drag_passage]) || args[:body_shape] == "Blunted Cone"
    #     r = p_class.Rp_e + h_0
        
    #     state[:vi] = -acos(1 / eccentricity_in * (semimajoraxis_in * (1 - eccentricity_in^2) / r - 1))
        
    #     if args[:montecarlo] == true
    #         state = monte_carlo_true_anomaly(state, args)
    #         apoapsis = state[:Apoapsis]
    #         periapsis = state[:Periapsis]
    #     end
    # end

    # Initial Model Definition
    # Body
    # if args[:body_shape] == "Spacecraft"
    #     b_class = args[:spacecraft_model]
    #     if args[:print_res]
    #         println("Area: " * string(get_spacecraft_reference_area(b_class)) * " m^2")
    #     end

    # elseif args[:body_shape] == "Blunted Cone" # TODO: Change this to new spacecraft model method

    #     Mass = mass
    #     Delta = δ
    #     NoseRadius = nose_radius
    #     BaseRadius = base_radius
    #     Area_tot = pi * BaseRadius^2

    #     b_class = Body(Mass, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Area_tot, Delta, NoseRadius, BaseRadius)

    # end

    # This is all defined in args already, just need to make sure it is all passed into 
    # function initialconditions()
    #     a = semimajoraxis_in
    #     e = eccentricity_in
    #     i = deg2rad(state[:Inclination])
    #     Ω = deg2rad(state[:Ω])
    #     ω = deg2rad(state[:ω])
    #     vi = state[:vi]
    #     println("vi: ", vi)
    #     m0 = mass
    #     year = args[:year]
    #     month = args[:month]
    #     day = args[:day]
    #     hour = args[:hours]
    #     min = args[:minutes]
    #     second = args[:secs]
    #     time_rot = args[:planettime]
    #     el_time = 0.0 # Elapsed time in seconds
    #     DateTimeIC = from_utc(DateTime(year, month, day, hour, min, second))
    #     DateTimeJ2000 = from_utc(DateTime(2000, 1, 1, 12, 0, 0))
    #     ic = Initial_condition(a, e, i, Ω, ω, vi, m0, year, month, day, hour, min, second, time_rot, el_time, DateTimeIC, DateTimeJ2000)

    #     return ic
    # end

    # ic_class = initialconditions()

    # function aerodynamics()
    #     δ = deg2rad(0)
    #     α = deg2rad(args[:α])
    #     thermal_accomodation_factor = args[:thermal_accomodation_factor]
    #     reflection_coefficient = args[:reflection_coefficient]
    #     thermal_contact = 0
    #     heat_rate_limit = args[:max_heat_rate]
    #     heat_load_limit = args[:max_heat_load]

    #     a = Aerodynamics(δ, α, thermal_accomodation_factor, reflection_coefficient, thermal_contact, heat_rate_limit, heat_load_limit)

    #     return a        
    # end

    # a_class = aerodynamics()

    # Engine
    # args[:phi] = deg2rad(args[:phi])
    # function engine()
    #     ϕ = args[:phi]
    #     g_e = 9.81
    #     T = args[:thrust]
    #     Isp = 200

    #     e = Engines(ϕ, g_e, T, Isp)

    #     return e
    # end

    # e_class = engine()

    # function model()
    #     body = b_class
    #     planet = p_class
    #     initialcondition = ic_class
    #     aerodynamics = a_class
    #     engine = e_class

    #     # m = Model(body, planet, aerodynamics, engine, initialcondition)

    #     return m
    # end

    # m = model()

    # Define gram atmosphere
    # gram = nothing
    # gram_atmosphere = nothing
    # if uppercase(args[:density_model]) == "GRAM"
    #     sys = pyimport("sys")
    #     os = pyimport("os")
    #     if !(args[:directory_Gram] in pyconvert(Vector{String}, sys.path))
    #         sys.path.append(args[:directory_Gram])
    #     end
    #     gram = pyimport("gram")
    #     inputParameters = Dict("earth" => gram.EarthInputParameters(),
    #                         "mars" => gram.MarsInputParameters(),
    #                         "venus" => gram.VenusInputParameters(),
    #                         "titan" => gram.TitanInputParameters())
        
    #     namelistReaders = Dict("earth" => gram.EarthNamelistReader(),
    #                         "mars" => gram.MarsNamelistReader(),
    #                         "venus" => gram.VenusNamelistReader(),
    #                         "titan" => gram.TitanNamelistReader())
            
    #     atmospheres = Dict("earth" => gram.EarthAtmosphere(),
    #                     "mars" => gram.MarsAtmosphere(),
    #                     "venus" => gram.VenusAtmosphere(),
    #                     "titan" => gram.TitanAtmosphere())

    #     planet_name = m.planet.name
    #     input_parameters = inputParameters[planet_name]

    #     # Mars has some weird specific parameters, so this line is just to check to make sure the it doesn't do it for the other planets
    #     if planet_name == "mars"
    #         # input_parameters.dataPath = os.path.join(os.path.dirname(os.path.abspath(@__FILE__)),"..", "GRAM_Data", "Mars", "data", "")
    #         input_parameters.dataPath = args[:directory_Gram_data] * "/Mars/data/"
    #         if !Bool(os.path.exists(input_parameters.dataPath))
    #             throw(ArgumentError("GRAM data path not found: " * input_parameters.dataPath))
    #         end
    #     end

    #     if planet_name == "earth"
    #         # input_parameters.dataPath = os.path.join(os.path.dirname(os.path.abspath(@__FILE__)),"..", "GRAM_Data", "Mars", "data", "")
    #         input_parameters.dataPath = args[:directory_Gram_data] * "/Earth/data/"
    #         if !Bool(os.path.exists(input_parameters.dataPath))
    #             throw(ArgumentError("GRAM data path not found: " * input_parameters.dataPath))
    #         end
    #     end

    #     reader = namelistReaders[planet_name]
    #     reader.tryGetSpicePath(input_parameters)

    #     gram_atmosphere = atmospheres[planet_name]
    #     gram_atmosphere.setInputParameters(input_parameters)
        
    #     if planet_name == "earth"
    #         gram_atmosphere.setMERRA2Parameters(0, -90.0, 90.0, 0.0, 359.99999)
    #     end
    #     gram_atmosphere.setPerturbationScales(1.5)
    #     gram_atmosphere.setMinRelativeStepSize(0.5)
    #     if args[:montecarlo] == 1
    #         gram_atmosphere.setSeed(Int(round(rand()*10000)))
    #     else
    #         gram_atmosphere.setSeed(1001)
    #         # gram_atmosphere.setSeed(Int(round(rand()*10000)))
    #     end

    #     if planet_name == "mars"
    #         gram_atmosphere.setMOLAHeights(false)
    #     end

    #     ttime = gram.GramTime()
    #     ttime.setStartTime(args[:year], args[:month], args[:day], args[:hours], args[:minutes], args[:secs], gram.UTC, gram.PET)
    #     gram_atmosphere.setStartTime(ttime)
    # end
    # Initialization - Reset all the config index for new simulation
    # cnf.count_aerobraking = 0
    # cnf.count_overcome_hr = 0
    # cnf.save_index_heat = 0
    # cnf.index_propellant_mass = 1
    # cnf.counter_random = 0
    # cnf.DU = args.simulation_settings.normalize ? args.initial_conditions.a : 1
    # cnf.TU = args.simulation_settings.normalize ? sqrt(cnf.DU^3 / args.environment_model.planet.μ) : 1
    # cnf.MU = args.simulation_settings.normalize ? args.dynamics_model.roots[1].dry_mass + args.dynamics_model.roots[1].prop_mass : 1

    # Save results
    if args.simulation_settings.results
        name = args.simulation_settings.results_directory * "/results"

        if !isdir(args.simulation_settings.results_directory)
            mkpath(args.simulation_settings.results_directory)
        end

        if args.simulation_settings.save_csv
            filename = name * ".csv"
            # if the file already exists, clear the current data
            if isfile(filename) && filesize(filename) > 0
                open(filename, "w") do file
                    truncate(file, 0)
                end
            end
        end

        # Initialize the arrow writer for plotting
        # arrow_writer = nothing
        # temp_name = nothing
    # if args[:plot] == true
        arrow_filename = name * ".feather"
        # if args.simulation_settings.verbose
        #     println("Temporary directory created for plotting: " * arrow_filename)
        # end
        arrow_writer = open(Arrow.Writer, arrow_filename)
    # end
        # save_csv(filename, args)
    end

    

    ##########################################################
    # RUN SIMULATION
    # cnf.heat_rate_limit = args[:max_heat_rate]
    # params = (cnf, m, solution)
    t_el = @elapsed begin
        execute_case(args, filename, arrow_filename; isolate_state=isolate_state)
    end
    # cnf = params[1]
    # m = params[2]
    # solution = params[3]
    ##########################################################

    # Finalize the arrow writer if plotting is enabled
    # if args[:plot] == true
    if arrow_writer !== nothing
        close(arrow_writer)
        if args.simulation_settings.verbose
            println("Arrow writer closed. Data saved to: " * arrow_filename)
        end
    end
    # end
    # Print final results
    # if Bool(args[:print_res])
    #     println(solution.physical_properties.ρ)
    #     println("ρ: " * string(maximum(solution.physical_properties.ρ)) * " kg/m^3")
    #     println("heat rate: " * string(maximum(maximum.(solution.performance.heat_rate))) * " W/cm^2")
    # end

    

    # if Bool(args[:print_res])
    #     println("Elapsed time: " * string(t_el) * " s")
    # end

    # if args[:plot] == true
    #     plots(state, m, name, args, params, arrow_filename)
    # end

    # rm(temp_name, recursive=true, force=true) # Remove the temporary directory used for plotting
end

function execute_campaign(args::SimulationConfiguration, state; isolate_state::Bool=true)
    return execute_campaign(args; isolate_state=isolate_state, state=state)
end
