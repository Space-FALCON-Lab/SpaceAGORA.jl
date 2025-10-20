include("../simulation/Run.jl")
# include("../config.jl") #TODO:Figure out how to run multiple times without having to comment this line out
include("../utils/maneuver_plans.jl")
include("../utils/attitude_control_plans.jl")
# include("SpacecraftModel.jl")

import .config
import .ref_sys
using Profile
using Random
using WGLMakie
using PlotlyJS
using Base.Threads
using Dash
# using Dash
using DashCoreComponents
using DashHtmlComponents
using OrdinaryDiffEq
# import .SpacecraftModel
# Define spacecraft model
spacecraft = config.SpacecraftModel()
# Add bodies to the spacecraft model
main_bus = config.Link(root=true, 
                        r=SVector{3, Float64}(0.0, 0.0, 0.0), 
                        q=SVector{4, Float64}([0, -0.6321683, -0.07370895, 0.7713171]), # true odyssey
                        # q=SVector{4, Float64}([0, 0, 0, 1]),
                        # q=SVector{4, Float64}([0, 0, 0.6371, 0.7707]), # 0 inc
                        # q=SVector{4, Float64}([0, -0.6358, -0.0338, 0.7711]), # 90 inc
                        ṙ=SVector{3, Float64}([0,0,0]), 
                        dims=SVector{3, Float64}([2.2,2.6,1.7]), 
                        ref_area=2.6*1.7,
                        m=391.0, 
                        gyro=4,
                        max_torque=2.0,
                        max_h=200.0,
                        attitude_control_rate=1.0/10.0, # 3 Hz
                        J_rw=MMatrix{3, 4, Float64}([1.0 0.0 0.0 0.57735; 0.0 1.0 0.0 0.57735; 0.0 0.0 1.0 0.57735]),#0.57735
                        attitude_control_function=lqr_constant_α_β)

L_panel = config.Link(r=SVector{3, Float64}(0.0, -2.6/2 - 3.89/4, 0.0), 
                        # q=SVector{4, Float64}([0, 0.4617, 0, 0.8870]),
                        q=SVector{4, Float64}([0, 0, 0, 1]),
                        ṙ=SVector{3, Float64}([0,0,0]), 
                        dims=SVector{3, Float64}([0.01, 3.89/2, 1.7]), 
                        ref_area=3.89*1.7/2,
                        m=10.0, 
                        gyro=0)
R_panel = config.Link(r=SVector{3, Float64}(0.0, 2.6/2 + 3.89/4, 0.0),
                        # q=SVector{4, Float64}([0, 0.4617, 0, 0.8870]),
                        q=SVector{4, Float64}([0, 0, 0, 1]),
                        ṙ=SVector{3, Float64}([0,0,0]), 
                        dims=SVector{3, Float64}([0.01, 3.89/2, 1.7]), 
                        ref_area=3.89*1.7/2,
                        m=10.0, 
                        gyro=0)

config.add_body!(spacecraft, main_bus, prop_mass=50.0)
config.add_body!(spacecraft, L_panel)
config.add_body!(spacecraft, R_panel)



L_panel_joint = config.Joint(main_bus, L_panel)
R_panel_joint = config.Joint(R_panel, main_bus)
config.add_joint!(spacecraft, L_panel_joint)
config.add_joint!(spacecraft, R_panel_joint)



println("Spacecraft model initialized with $(length(spacecraft.links)) bodies.")
# println("Spacecraft roots: $spacecraft.roots")
println("Spacecraft COM: $(config.get_COM(spacecraft, main_bus))")
println("Spacecraft MOI: $(config.get_inertia_tensor(spacecraft, main_bus))")
# config.model.body = spacecraft
# println("Number of reaction wheels: ", config.model.body.n_reaction_wheels)
# println("Number of reaction wheels true: $(spacecraft.n_reaction_wheels)")
args = Dict(# Misc Simulation
            :results => 1,                                                                                      # Generate csv file for results True=1, False=0
            :passresults => 1,                                                                                  # Pass results as output True=1, False=0
            :print_res => 1,                                                                                    # Print some lines True=1, False=0
            :directory_results => "/workspaces/SpaceAGORA.jl/output/odyssey_quat_test",                # Directory where to save the results
            :directory_Gram => "/workspaces/SpaceAGORA.jl/GRAMpy",                                                    # Directory where Gram is
            :directory_Gram_data => "/workspaces/SpaceAGORA.jl/GRAM_Data",                                            # Directory where Gram data is
            :directory_Spice => "/workspaces/SpaceAGORA.jl/GRAM_Data/SPICE",                                          # Directory where SPICE files are located
            :Gram_version => 0,                                                                                 # MarsGram x file to use
            :montecarlo_analysis => 0,                                                                          # Generate csv file for Montecarlo results True=1, False=0
            :plot => 1,
            :plot_type => "multi", #multi or single                                                                                         # Generate pdf plots of results True=1, False=0
            :dashboard => true,                                                                                     # Start dashboard True=1, False=0
            :filename => 1,                                         # Filename with specifics of simulation, True =1, False=0
            :machine => "",                                         # choices=['Laptop' , 'Cluster' , 'Aero' , 'Desktop_Home','Karnap_Laptop']
            :integrator => "Julia",                                 # choices=['Costumed', 'Julia'] Costumed customed integrator, Julia DifferentialEquations.jl library integrator, only for drag passage, others phases use RK4
            :normalize => 0,                                       # Normalize the integration True=1, False=0
            :closed_form => 0,                                     # Closed form solution True=1, False=0
            # Type of Mission
            :type_of_mission => "Time",                           # choices=['Drag Passage' , 'Orbits' , 'Aerobraking Campaign']
            :keplerian => 1,                                        # Do not include drag passage: True=1, False=0
            :number_of_orbits => 1,                                 # Number of aerobraking passage
            :mission_time => 6000.0, #60000.0,                                  # Mission time in seconds, used only for Time mission type
            :orientation_sim => false,                               # Orientation simulation True=1, False=0, if false, will only propagate position
            :rng => MersenneTwister(12345),                               # Random number generator for reproducibility
            :synchronized_threads => true,                           # Synchronized threads for multi spacecraft simulation
            :orb_picture_viz => true,                               # Visualize orbit picture during simulation

        
            
            #swarm simulation configuration
            :swarm_config => 1,                                      # Swarm configuration: 0 = no swarm, 1 = yes swarm
            :n_spacecraft => 1,                                       # Number of spacecraft to simulate
            :n_target_obj => 1,                                       # Number of target objects to simulate

            # Physical Model
            :planet => 0,                                           # Earth = 0, Mars = 1, Venus = 2
            :planettime => 0.0,                                     # Initial time of the mission, sec. Importalmnt for J2 effect and rotation of the planet
            :gravity_model => "Inverse Squared and J2 effect",      # choices=['Constant' , 'Inverse Squared' , 'Inverse Squared and J2 effect', 'GRAM']
            :density_model => "Exponential",                               # choices=['Constant' , 'Exponential' , 'Gram']
            :topography_model => "Spherical Harmonics",                             # choices=['None' , 'Spherical Harmonics']
            :topography_harmonics_file => "/workspaces/SpaceAGORA.jl/Topography_harmonics_data/MOLA.csv", # File with the topography harmonics coefficients
            :topo_degree => 90,                                     # Maximum degree of the topography harmonics (Defined in the file)
            :topo_order => 90,                                      # Maximum order of the topography harmonics (Defined in the file)
            :wind => 1,                                             # Wind calculation only if density model is Gram True=1, False=0
            :aerodynamic_model => "Mach-dependent",                 # choices=['Cd and Cl Constant' , 'Mach-dependent' , 'No-Ballistic flight with axial coefficient']: "Mach-dependent" specific for spacecraft shape, "No-Ballistic flight" specific for blunted-cone shape
            :thermal_model => "Maxwellian Heat Transfer",    # choices=['Maxwellian Heat Transfer' , 'Convective and Radiative']: "Maxwellian Heat Transfer" specific for spacecraft shape, "Convective and Radiative" specific for blunted-cone shape
            :interactive_forces => true,                     # Interactive forces between spacecraft toggle. True= include interactive forces, False= do not include interactive forces
            :spice_call => false,                                 # Call SPICE for every time step
            # Perturbations
            :n_bodies => ["Sun"],                                        # Add names of bodies you want to simulate the gravity of to a list. Keep list empty if not required to simulate extra body gravity.
            :srp => 0,                                             # Solar Radiation Pressure True=1, False=0
            :gravity_harmonics => 1,                                            # Gravity Spherical harmonics True=1, False=0
            :gravity_harmonics_file => "/workspaces/SpaceAGORA.jl/Gravity_harmonics_data/Mars50c.csv", # File with the gravity harmonics coefficients
            :L => 50,                                              # Maximum degree of the gravity harmonics (Defined in the file)
            :M => 50,                                              # Maximum order of the gravity harmonics (Defined in the file)

            # Rates
            :trajectory_rate => 100.0,                              # Rate at which the trajectory in drag passage integrate using RK4
            :flash1_rate => 3.0,                                    # Rate at which Control Mode-1 is called
            :save_rate => 5.0,                                      # Rate at which the data trajectory are saved
            :simulation_filename => "filename",
            
            # Body
            :body_shape => "Spacecraft",                            # choices=['Spacecraft' , 'Blunted Cone']
            :max_heat_rate => 0.15,                                 # Max heat rate the heat rate control will start to react to
            :max_heat_load => 30.0,                                 # Max heat load the heat load control will not be overcomed
            # :dry_mass => 411.0,                                     # Initial dry mass of body in kg
            # :prop_mass => 50.0,                                     # Initial propellant mass of body in kg
            :reflection_coefficient => 0.9,                         # Diffuse reflection sigma =0, for specular reflection sigma = 1
            :thermal_accomodation_factor => 1.0,                    # Thermal accomodation factor, Shaaf and Chambre
            :α => 90.0,                                             # Max angle of attack of solar panels

            # # Fill for Spacecraft body shape only
            # :length_sat => 2.2,                                     # Length of the satellite in m
            # :height_sat => 1.7,                                     # Height of the satellite in m
            # :width_sat => 2.6,                                      # Width of the satellite in m
            # :length_sp => 3.76,                                     # Length of the solar panels in m
            # :height_sp => 1.93,                                     # Height of the solar panels in m

            # # Fill for Blunted Cone body shape only
            # :cone_angle => 70.0,                                    # Cone angle of the blunted cone in deg
            # :base_radius => 2.65/2,                                 # Base radius of the blunted cone in m
            # :nose_radius => 0.6638,                                 # Nose radius of the blunted cone in m
            :spacecraft_model => spacecraft,                            # Spacecraft model with bodies and joints
            
            # Engine
            :thrust => 4.0,                                         # Maximum magnitude thrust in N
            
            # Control Mode
            :control_mode => 0,                                     # Use Rotative Solar Panels Control:  False=0, Only heat rate=1, Only heat load=2, Heat rate and Heat load = 3
            :security_mode => 1,                                    # Security mode that set the angle of attack to 0 deg if predicted heat load exceed heat load limit
            :second_switch_reevaluation => 1,                       # Reevaluation of the second switch time when the time is closer to it
            :control_in_loop => 1,                                  # Control in loop, control called during integration of trajectory, full state knowledge
            :flash2_through_integration => 0,                       # Integration of the equations of motion and lambda to define time switches and revaluation second time switch
            :solar_panel_control_rate => 1.0/3.0,                        # Rate at which the solar panel controller is called

            # Initial Conditions
            :initial_condition_type => 3,                          # Initial Condition ra,hp = 0, Initial Condition v, gamma = 1, Initial Condition a, e = 2, type=3 is constellation simulation
        :ra_initial_a => 6571e3,                # Initial Apoapsis Radius for ISS orbit in m (Earth radius + ~400 km)
        :ra_initial_b => 6571e3,                # Final Apoapsis Radius for ISS orbit in m
        :ra_step => 5e10,                       # Step Apoapsis Radius for for-loop in m (single value)
        :hp_initial_a => 6771e3,                # Initial Periapsis Altitude for ISS orbit in m (circular orbit)
        :hp_initial_b => 6771e3,                # Final Periapsis Altitude for ISS orbit in m
        :hp_step => 10000000.0,                        # Step Periapsis Radius for for-loop in m (single value)
        :v_initial_a => 7660.0,                 # Initial Velocity (m/s) for ISS orbit
        :v_initial_b => 7660.0,                 # Final Velocity (m/s) for ISS orbit
        :v_step => 10000.0,                         # Step Velocity (m/s) for for-loop (single value)
        :a_initial_a => 6771e3,                 # Initial Semi-major axis (m) for ISS orbit
        :a_initial_b => 6771e3,                 # Final Semi-major axis (m) for ISS orbit
        :a_step => 1.0,                         # Step Semi-major axis (m) for for-loop (single value)
        :e_initial_a => 0.01,                    # Initial Eccentricity for ISS orbit (circular)
        :e_initial_b => 0.01,                    # Final Eccentricity for ISS orbit
        :e_step => 0.1,                        # Step Eccentricity for for-loop (single value)
            
            :orientation_type => 0,                                   # Initial Condition orientation = 0, Initial Condition orientation and velocity = 1
            :γ_initial_a => -2.5,                                    # Initial Gamma (deg) for for-loop if initial conditions are in v and gamma
            :γ_initial_b => 7.0,                                    # Final Gamma (deg) for for-loop if initial conditions are in v and gamma
            :γ_step => 100,                                         # Step Gamma (deg) for for-loop if initial conditions are in v and gamma
            :inclination => 93.522,                                   # Inclination Orbit, deg
            :ω => 109.7454,                                              # AOP, deg
            :Ω => 28.1517,                                              # RAAN, deg
            :ν => 320.0,                                              # True Anomaly, deg, set to negative value to start at apoapsis if type_of_mission is "Time"
            :EI => 250.0,                                           # Entry Interface, km
            :AE => 250.0,                                           # Atmospheric Exit, km
            :year => 2001,                                          # Mission year
            :month => 11,                                           # Mission month
            :day => 6,                                             # Mission day
            :hours => 19,                                           # Mission hour
            :minutes => 0,                                         # Mission minute
            :secs => 32.0,                                          # Mission second
            
            # Final Conditions
            :final_apoapsis => 3390.0e3+503e3, # 5088116.837416616, # 4905.974818462152e3                  # Final apoapsis radius if aerobraking campaign

            # Do not change
            :heat_load_sol => 0,                                    # Heat load solution leave it to 0 and change it only for control mode = 2:  Max energy depletaion=0, Min energy depletion=1, One switch max-min=2, One switch min-max = 3
            :thrust_control => "None",                              # choices=['None' , 'Aerobraking Maneuver' , 'Drag Passage Firing']
            :phi => 180.0,                                          # Thrust Angle, deg
            :delta_v => 0,                                          # Delta-v of Aerobraking Manuver,m/s
            :apoapsis_targeting => 0,                               # Apoapsis Targeting Enabled
            :ra_fin_orbit => 25000e3,                               # Target final apoapsis for the orbit, m
            :maneuver_plan => Odyssey_firing_plan,                # Maneuver plan function
            
            # Monte Carlo Simulations
            :montecarlo => 0,                                       # Run Monte Carlo simulation True=1, False=0
            :initial_montecarlo_number => 1,                        # Initial Monte Carlo sample number
            :montecarlo_size => 1000,                               # number of Monte Carlo samples
            
            # Monte Carlo Perturbations
            :CD_dispersion => 10.0,                                 # Max dispersion of CD for Uniform Distribution, %
            :CL_dispersion => 10.0,                                 # Max dispersion of CL for Uniform Distribution, %
            :rp_dispersion => 87.0*0.05/3,                                  # Max dispersion for initial vacuum periapsis radius following uniform distribution, km
            :ra_dispersion => 28559.0*0.05/3,                                  # Max dispersion for initial apoapsis radius following uniform distribution, km
            :i_dispersion => 0.25,                                  # Max dispersion for initial inclination following uniform distribution, deg
            :Ω_dispersion => 0.25,                                  # Max dispersion for initial right ascension of the ascending node following uniform distribution, deg
            :ω_dispersion => 0.25,                                  # Max dispersion for initial argument of periapsis following uniform distribution, deg
            :vi_dispersion => 0.025,                                # Max dispersion for initial true anomaly following uniform distribution, deg
            
            # MonteCarlo Perturbation Guidance - Closed Form Solution (only for online)
            :ra_dispersion_gnc => 0.25,                             # Max dispersion for initial apoapsis radius used by gnc following uniform distribution, km
            :rp_dispersion_gnc => 0.25,                             # Max dispersion for initial periapsis radius used by gnc following uniform distribution, km
            :i_dispersion_gnc => 0.025,                             # Max dispersion for initial inclination used by gnc following uniform distribution, deg
            :Ω_dispersion_gnc => 0.025,                             # Max dispersion for initial right ascension of the ascending node used by gnc following uniform distribution, deg
            :ω_dispersion_gnc => 0.0,                               # Max dispersion for initial argument of periapsis used by gnc following uniform distribution, deg
            :vi_dispersion_gnc => 0.0,                              # Max dispersion for initial true anomaly used by gnc following uniform distribution, deg
            
            # Online trajectory control (heat rate)
            :ρ_mudispersion_gnc => 0.0,                             # Mean dispersion of rho for Gaussian Distribution, %
            :ρ_sigmadispersion_gnc => 1.0,                          # Std dispersion of rho for Gaussian Distribution, %
            :T_mudispersion_gnc => 0.0,                             # Mean dispersion of T for Gaussian Distribution, %
            :T_sigmadispersion_gnc => 1.0,                          # Std dispersion of T for Gaussian Distribution, %
            :S_mudispersion_gnc => 0.0,                             # Mean dispersion of S for Gaussian Distribution, %
            :S_sigmadispersion_gnc => 1.0,                          # Std dispersion of S for Gaussian Distribution, %
            :multiplicative_factor_heatload => 1.0,                 # Multiplicative factor for heat rate prediction when calculated heat load

            :a_tol => 1e-5,                                         # Absolute tolerance for integration
            :r_tol => 1e-3,                                         # Relative tolerance for integration
            :a_tol_orbit => 1e-8,                                    # Absolute tolerance for orbit integration (outside atmosphere, i.e., step 1 and step 3)
            :r_tol_orbit => 1e-6,                                    # Relative tolerance for orbit integration (outside atmosphere, i.e., step 1 and step 3)
            :a_tol_drag => 1e-8,                                       # Absolute tolerance for drag passage integration (inside atmosphere, i.e., step 2)
            :r_tol_drag => 1e-6,                                       # Relative tolerance for drag passage integration (inside atmosphere, i.e., step 2)
            :a_tol_quaternion => 1e-9,                                  # Absolute tolerance for quaternion integration (inside atmosphere, i.e., step 2)
            :r_tol_quaternion => 1e-7,                                  # Relative tolerance for quaternion integration (inside atmosphere, i.e., step 2)
            :dt_max => 1.0,                                         # Maximum time step for integration, s
            :dt_max_orbit => 30.0,                                   # Maximum time step for orbit integration (outside atmosphere, i.e., step 1 and step 3), s
            :dt_max_drag => 1.0,                                    # Maximum time step for drag passage
            :Odyssey_sim => 0                                      # Simulate Odyssey Mission
            )

        # Create a dictionary of spacecraft buses using the main_bus as a template
        spacecraft_buses = Dict{Int, typeof(config.SpacecraftModel())}()

        for i in 1:args[:n_spacecraft]
                # bus = config.Link(root=true, 
                #                                   r=main_bus.r, 
                #                                   q=main_bus.q, 
                #                                   ṙ=main_bus.ṙ, 
                #                                   dims=main_bus.dims, 
                #                                   ref_area=main_bus.ref_area,
                #                                   m=main_bus.m, 
                #                                   gyro=main_bus.gyro,
                #                                   max_torque=main_bus.max_torque,
                #                                   max_h=main_bus.max_h,
                #                                   attitude_control_rate=main_bus.attitude_control_rate,
                #                                   J_rw=main_bus.J_rw,
                spacecraft_buses[i] = deepcopy(args[:spacecraft_model])
                spacecraft_buses[i].laser_effector = true
                spacecraft_buses[i].laser_power = 100e3 # 100 kW PLACEHOLDER VALUE
                spacecraft_buses[i].laser_range = 1000e3 # 1000 km PLACEHOLDER VALUE
                spacecraft_buses[i].uid = i #setting spacecraft id

                # Initialize access logging and storage
                spacecraft_buses[i].access_log = Vector{Tuple{Float64, Float64, Int64}}() # (start_time, end_time, target_id)
                spacecraft_buses[i].access_storage = Dict{Int64, config.Access}() # Dictionary storing access windows by target id
                
        end

        # Create a dictionary of target objects using the main_bus as a template
        target_objects = Dict{Int, typeof(config.SpacecraftModel())}()

        for i in 1:args[:n_target_obj]
                # target = config.Link(root=true, 
                #                                          r=main_bus.r, 
                #                                          q=main_bus.q, 
                #                                          ṙ=main_bus.ṙ, 
                #                                          dims=main_bus.dims, 
                #                                          ref_area=main_bus.ref_area,
                #                                          m=main_bus.m, 
                #                                          gyro=main_bus.gyro,
                #                                          max_torque=main_bus.max_torque,
                #                                          max_h=main_bus.max_h,
                #                                          attitude_control_rate=main_bus.attitude_control_rate,
                #                                          J_rw=main_bus.J_rw,
                #                                          attitude_control_function=main_bus.attitude_control_function)
                target_objects[i] = deepcopy(args[:spacecraft_model])
                target_objects[i].uid = i #setting debris object id
                target_objects[i].laser_effector = false

                # Initialize access logging and storage
                target_objects[i].access_log = Vector{Tuple{Float64, Float64, Int64}}() #(start_time, end_time, spacecraft_id)
                target_objects[i].access_storage = Dict{Int64, config.Access}() # Dictionary storing access windows by spacecraft id
        end

        # Create a dictionary of initial conditions for each spacecraft
        spacecraft_initial_conditions = Dict{Int, Dict{Symbol, Any}}()
        rng = args[:rng]

        for i in 1:args[:n_spacecraft]
                spacecraft_initial_conditions[i] = Dict(
                        :ra_initial_a => args[:ra_initial_a]+1000.0, # Adding 1 km to initial apoapsis radius of debris
                        :ra_initial_b => args[:ra_initial_a]+1000.0,
                        :ra_step => args[:ra_step],
                        :hp_initial_a => args[:hp_initial_a],
                        :hp_initial_b => args[:hp_initial_a],
                        :hp_step => args[:hp_step],
                        :v_initial_a => args[:v_initial_a],
                        :v_initial_b => args[:v_initial_a],
                        :v_step => args[:v_step],
                        :a_initial_a => args[:a_initial_a]+1000.0,
                        :a_initial_b => args[:a_initial_a]+1000.0,
                        :a_step => args[:a_step],
                        :e_initial_a => args[:e_initial_a],
                        :e_initial_b => args[:e_initial_a],
                        :e_step => args[:e_step],
                        :γ_initial_a => args[:γ_initial_a],
                        :γ_initial_b => args[:γ_initial_a],
                        :γ_step => args[:γ_step],
                        :inclination => args[:inclination],
                        :ω => args[:ω],
                        :Ω => args[:Ω],
                        :ν => args[:ν]
                )

        end
        # Create a dictionary of initial conditions for each debris object
        target_initial_conditions = Dict{Int, Dict{Symbol, Any}}()

        for i in 1:args[:n_target_obj]
                target_initial_conditions[i] = Dict(
                        :ra_initial_a => args[:ra_initial_a], # Adding 1 km to initial apoapsis radius of debris
                        :ra_initial_b => args[:ra_initial_a],
                        :ra_step => args[:ra_step],
                        :hp_initial_a => args[:hp_initial_a],
                        :hp_initial_b => args[:hp_initial_a],
                        :hp_step => args[:hp_step],
                        :v_initial_a => args[:v_initial_a],
                        :v_initial_b => args[:v_initial_a],
                        :v_step => args[:v_step],
                        :a_initial_a => args[:a_initial_a],
                        :a_initial_b => args[:a_initial_a],
                        :a_step => args[:a_step],
                        :e_initial_a => args[:e_initial_a],
                        :e_initial_b => args[:e_initial_a],
                        :e_step => args[:e_step],
                        :γ_initial_a => args[:γ_initial_a],
                        :γ_initial_b => args[:γ_initial_a],
                        :γ_step => args[:γ_step],
                        :inclination => args[:inclination],
                        :ω => args[:ω],
                        :Ω => args[:Ω],
                        :ν => args[:ν]
                )
        end
        args[:target_objects] = target_objects
        args[:spacecraft_buses] = spacecraft_buses
        args[:target_initial_conditions] = target_initial_conditions
        args[:spacecraft_initial_conditions] = spacecraft_initial_conditions

        println("Initialized $(length(target_objects)) target objects.")

        println("Initialized $(length(spacecraft_buses)) spacecraft buses.")

# @profview run_analysis(args)
t = @elapsed begin          
        # Run the simulation
        sol = run_analysis(args)
        # if Bool(args[:passresults])
        #     println("Ra initial = " * string((sol.orientation.oe[1][1] * (1 + sol.orientation.oe[2][1]))* 1e-3) * " km, Ra new = " * string((sol.orientation.oe[1][end] * (1 + sol.orientation.oe[2][end]))* 1e-3) * " km - Actual periapsis altitude = " * string(minimum(sol.orientation.alt) * 1e-3) * " km - Target Ra = " * string(args[:final_apoapsis] * 1e-3) * " km")
        # end
end

        println("COMPUTATIONAL TIME = " * string(t) * " s") 

# end
        # end
    # end
# end

# mc_runs = 50
# nominal_ra = args[:ra_initial_a]
# nominal_rp = args[:hp_initial_a]
# nominal_i = args[:inclination]
# nominal_Ω = args[:Ω]
# nominal_ω = args[:ω]
# for i in 24:mc_runs
#     a = @allocated begin
#         t = @elapsed begin
#             args[:directory_results] = "/workspaces/ABTS.jl/output/odyssey_MC_polyfit_atmo_orig_disp/" * string(i)
#             println("Monte Carlo Run: ", i)
#             if i == 1
#                 args[:print_res] = 1
#             else
#                 args[:print_res] = 0
#             end
#             # args[:print_res] = 1
#             # println(randn()*sqrt(args[:ra_dispersion]) * 1e3)
#             args[:ra_initial_a] = nominal_ra + randn()*sqrt(args[:ra_dispersion]) * 1e3
#             args[:hp_initial_a] = nominal_rp + randn()*sqrt(args[:rp_dispersion]) * 1e3
#             args[:inclination] = nominal_i + randn()*sqrt(args[:i_dispersion])
#             args[:Ω] = nominal_Ω + randn()*sqrt(args[:Ω_dispersion])
#             args[:ω] = nominal_ω + randn()*sqrt(args[:ω_dispersion])

#             # Run the simulation
#             sol = run_analysis(args)

#             if Bool(args[:passresults])
#                 println("Ra initial = " * string((sol.orientation.oe[1][1] * (1 + sol.orientation.oe[2][1]))* 1e-3) * " km, Ra new = " * string((sol.orientation.oe[1][end] * (1 + sol.orientation.oe[2][end]))* 1e-3) * " km - Actual periapsis altitude = " * string(minimum(sol.orientation.alt) * 1e-3) * " km - Target Ra = " * string(args[:final_apoapsis] * 1e-3) * " km")
#             end
#         end

#         println("COMPUTATIONAL TIME = " * string(t) * " s")
#         config.reset_config()
#     end
#     println("Memory allocated = " * string(a / 1e6) * " MB")
# end