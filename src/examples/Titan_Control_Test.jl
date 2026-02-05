include("../simulation/Run.jl")
include("../config.jl") #TODO:Figure out how to run multiple times without having to comment this line out
include("../utils/maneuver_plans.jl")
include("../utils/attitude_control_plans.jl")
# include("SpacecraftModel.jl")

import .config
import .ref_sys
using Profile

spacecraft = config.SpacecraftModel()
# Add bodies to the spacecraft model
main_bus = config.Link(root=true, 
                        r=SVector{3, Float64}(0.0, 0.0, 0.0), 
                        q=SVector{4, Float64}([0, -0.6321683, -0.07370895, 0.7713171]),
                        ṙ=SVector{3, Float64}([0,0,0]), 
                        dims=SVector{3, Float64}([2.2,2.6,1.7]), 
                        ref_area=2.6*1.7,
                        m=2400.0, 
                        gyro=4,
                        max_torque=2.0,
                        max_h=200.0,
                        attitude_control_rate=1.0/10.0, # 3 Hz
                        J_rw=MMatrix{3, 4, Float64}([1.0 0.0 0.0 0.57735; 0.0 1.0 0.0 0.57735; 0.0 0.0 1.0 0.57735]),#0.57735
                        attitude_control_function=lqr_constant_α_β)

L_panel = config.Link(r=SVector{3, Float64}(0.0, -2.6/2 - 3.89/4, 0.0), 
                        q=SVector{4, Float64}([0, 0, 0, 1]),
                        ṙ=SVector{3, Float64}([0,0,0]), 
                        dims=SVector{3, Float64}([0.01, 3.89/2, 1.7]), 
                        ref_area=3.89*1.7/2,
                        m=10.0, 
                        gyro=0)
R_panel = config.Link(r=SVector{3, Float64}(0.0, 2.6/2 + 3.89/4, 0.0),
                        q=SVector{4, Float64}([0, 0, 0, 1]),
                        ṙ=SVector{3, Float64}([0,0,0]), 
                        dims=SVector{3, Float64}([0.01, 3.89/2, 1.7]), 
                        ref_area=3.89*1.7/2,
                        m=10.0, 
                        gyro=0)

config.add_body!(spacecraft, main_bus, prop_mass=1000.0)
config.add_body!(spacecraft, L_panel)
config.add_body!(spacecraft, R_panel)

L_panel_joint = config.Joint(main_bus, L_panel)
R_panel_joint = config.Joint(R_panel, main_bus)
config.add_joint!(spacecraft, L_panel_joint)
config.add_joint!(spacecraft, R_panel_joint)

args = Dict(# Misc Simulation
            :results => true,                                                                                      # Generate csv file for results True=1, False=0
            :passresults => true,                                                                                  # Pass results as output True=1, False=0
            :print_res => true,                                                                                    # Print some lines True=1, False=0
            :directory_results => "output/titan_control_test",            # Directory where to save the results
            :directory_Gram => "GRAMpy",                   # Directory where Gram is
            :directory_Gram_data => "GRAM_Data",           # Directory where Gram data is
            :directory_Spice => "GRAM_Data/SPICE",         # Directory where SPICE files are located
            :Gram_version => 0,                                                                                 # MarsGram x file to use
            :montecarlo_analysis => 0,                                                                          # Generate csv file for Montecarlo results True=1, False=0
            :plot => true,                                                                                         # Generate pdf plots of results True=1, False=0
            :filename => true,                                         # Filename with specifics of simulation, True =1, False=0
            :machine => "",                                         # choices=['Laptop' , 'Cluster' , 'Aero' , 'Desktop_Home','Karnap_Laptop']
            :integrator => "Julia",                                 # choices=['Costumed', 'Julia'] Costumed customed integrator, Julia DifferentialEquations.jl library integrator, only for drag passage, others phases use RK4
            :normalize => true,                                       # Normalize the results True=1, False=0
            :closed_form => true,                                   # Closed form solution for the drag passage True=1, False=0
            :save_csv => false,                                      # Save csv files of the trajectory True=1, False=0

            # Type of Mission
            :type_of_mission => "Orbits",                           # choices=['Drag Passage' , 'Orbits' , 'Aerobraking Campaign']
            :keplerian => false,                                        # Do not include drag passage: True=1, False=0
            :number_of_orbits => 5,                                 # Number of aerobraking passage
            :orientation_sim => false,                              # Orientation simulation True=1, False=0
            :save_steps => 10000,                                   # Number of steps to store in memory before saving to file

            # Physical Model
            :planet => 7,                                           # Earth = 0, Mars = 1, Venus = 2, Titan = 7
            :planettime => 0.0,                                  # Initial time of the mission, sec. Important for J2 effect and rotation of the planet
            :gravity_model => "Inverse Squared and J2 effect",      # choices=['Constant' , 'Inverse Squared' , 'Inverse Squared and J2 effect']
            
            :density_model => "Gram",                               # choices=['Constant' , 'Exponential' , 'Gram']
            :topography_model => "None",                             # choices=['None' , 'Spherical Harmonics']
            :topography_harmonics_file => "", # File with the topography harmonics coefficients
            :topo_degree => 90,                                     # Maximum degree of the topography harmonics (Defined in the file)
            :topo_order => 90,  
            :wind => 1,                                             # Wind calculation only if density model is Gram True=1, False=0
            :aerodynamic_model => "Mach-dependent",                 # choices=['Cd and Cl Constant' , 'Mach-dependent' , 'No-Ballistic flight with axial coefficient']: "Mach-dependent" specific for spacecraft shape, "No-Ballistic flight" specific for blunted-cone shape
            :thermal_model => "Maxwellian Heat Transfer",           # choices=['Maxwellian Heat Transfer' , 'Convective and Radiative']: "Maxwellian Heat Transfer" specific for spacecraft shape, "Convective and Radiative" specific for blunted-cone shape
            
            # Perturbations
            :n_bodies => ["Saturn"],                                        # Add names of bodies you want to simulate the gravity of to a list. Keep list empty if not required to simulate extra body gravity.
            :srp => 1,                                             # Solar Radiation Pressure True=1, False=0
            :eclipse => false,                                   # Eclipse model for SRP True=1, False=0
            :gravity_harmonics => 0,                                            # Gravity Spherical harmonics True=1, False=0
            :gravity_harmonics_file => "Gravity_harmonics_data/titan5.csv", # File with the gravity harmonics coefficients
            :L => 5,                                              # Maximum degree of the gravity harmonics (Defined in the file)
            :M => 5,                                              # Maximum order of the gravity harmonics (Defined in the file)
            :magnetic_field => false,                             # Magnetic Field True=1, False=0

            # Rates
            :trajectory_rate => 100.0,                              # Rate at which the trajectory in drag passage integrate using RK4
            :flash1_rate => 3.0,                                    # Rate at which Control Mode-1 is called
            :save_rate => 5.0,                                      # Rate at which the data trajectory are saved
            :control_in_loop => 1,                                  # Control in loop, control called during integration of trajectory, full state knowledge
            :flash2_through_integration => 0,                       # Integration of the equations of motion and lambda to define time switches and revaluation second time switch
            :struct_ctrl => 0,                                    # Structural thermal control, True=1, False=0
            :targeting_ctrl => 0,                                   # Targeting control True=1, False=0
            
            # Body
            :body_shape => "Spacecraft",                            # choices=['Spacecraft' , 'Blunted Cone']
            :max_heat_rate => 0.15,                                 # Max heat rate the heat rate control will start to react to
            :max_heat_load => 30.0,                                 # Max heat load the heat load control will not be overcomed
            :reflection_coefficient => 0.9,                         # Diffuse reflection sigma =0, for specular reflection sigma = 1
            :thermal_accomodation_factor => 1.0,                    # Thermal accomodation factor, Shaaf and Chambre
            :α => 90.0,                                             # Max angle of attack of solar panels
            :spacecraft_model => spacecraft,                         # Spacecraft model defined in SpacecraftModel.jl

            # Engine
            :thrust => 4.0,                                         # Maximum magnitude thrust in N
            
            # Control Mode
            :control_mode => 3,                                     # Use Rotative Solar Panels Control:  False=0, Only heat rate=1, Only heat load=2, Heat rate and Heat load = 3
            :security_mode => 0,                                    # Security mode that set the angle of attack to 0 deg if predicted heat load exceed heat load limit
            :second_switch_reevaluation => 1,                       # Reevaluation of the second switch time when the time is closer to it
            
            # Initial Conditions
            :initial_condition_type => 0,                           # Initial Condition ra,hp = 0, Initial Condition v, gamma = 1
            :orientation_type => 0,                               # Orientation type: 0 = Keplerian elements, 1 = Cartesian coordinates
            :ra_initial_a => 15000e3 + 2.5755e6, # 28523.95e3,                # Initial Apoapsis Radius for for-loop in m
            :ra_initial_b => 1e21,                               # Final Apoapsis Radius for for-loop in m
            :ra_step => 5e21,                                       # Step Apoapsis Radius for for-loop in m
            :hp_initial_a => 720_000.0,#176590.0,#188140.0                                 # Initial Periapsis Altitude for for-loop in m
            :hp_initial_b => 1590000.0,                              # Final Periapsis Altitude for for-loop in m
            :hp_step => 10000000.0,                                 # Step Periapsis Radius for for-loop in m
            :v_initial_a => 2100.0,                                 # Initial Velocity (m/s) for for-loop if initial conditions are in v and gamma
            :v_initial_b => 5000.0,                                 # Final Velocity (m/s) for for-loop if initial conditions are in v and gamma
            :v_step => 100.0e3,                                       # Step Velocity (m/s) for for-loop if initial conditions are in v and gamma
            :γ_initial_a => 7.5,                                    # Initial Gamma (deg) for for-loop if initial conditions are in v and gamma
            :γ_initial_b => 8.0,                                    # Final Gamma (deg) for for-loop if initial conditions are in v and gamma
            :γ_step => 1.0,                                         # Step Gamma (deg) for for-loop if initial conditions are in v and gamma
            :inclination => 85.37,                                   # Inclination Orbit, deg
            :ω => 90.0,                                              # AOP, deg
            :Ω => 64.495,                                              # RAAN, deg
            :EI => 1500.0,                                           # Entry Interface, km
            :AE => 1500.0,                                           # Atmospheric Exit, km
            :year => 2031,                                          # Mission year
            :month => 10,                                           # Mission month
            :day => 15,                                             # Mission day
            :hours => 14,                                           # Mission hour
            :minutes => 21,                                         # Mission minute
            :secs => 0.0,                                          # Mission second
            
            # Final Conditions
            :final_apoapsis => 1500e3 + 2.5755e6, # 4905.974818462152e3                  # Final apoapsis radius if aerobraking campaign

            # Do not change
            :heat_load_sol => 0,                                    # Heat load solution #leave it to 0 and change it only for control mode = 2:  Max energy depletaion=0, Min energy depletion=1, One switch max-min=2, One switch min-max = 3
            :thrust_control => "Aerobraking Maneuver",                              # choices=['None' , 'Aerobraking Maneuver' , 'Drag Passage Firing']
            :phi => 180.0,                                          # Thrust Angle, deg
            :delta_v => 0,                                          # Delta-v of Aerobraking Manuver,m/s
            :apoapsis_targeting => 0,                               # Apoapsis Targeting Enabled
            :ra_fin_orbit => 5000e3,                               # Target final apoapsis for the orbit, m
            :maneuver_plan => titan_firing_plan,                    # Maneuver Plan for the mission
            
            # Monte Carlo Simulations
            :montecarlo => 0,                                       # Run Monte Carlo simulation True=1, False=0
            :initial_montecarlo_number => 1,                        # Initial Monte Carlo sample number
            :montecarlo_size => 1000,                               # number of Monte Carlo samples
            
            # Monte Carlo Perturbations
            :CD_dispersion => 10.0,                                 # Max dispersion of CD for Uniform Distribution, %
            :CL_dispersion => 10.0,                                 # Max dispersion of CL for Uniform Distribution, %
            :rp_dispersion => 720*0.05/3,                                  # Max dispersion for initial vacuum periapsis radius following uniform distribution, km
            :ra_dispersion => 17000*0.05/3,                                  # Max dispersion for initial apoapsis radius following uniform distribution, km
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
            :a_tol_drag => 1e-10,                                       # Absolute tolerance for drag passage integration (inside atmosphere, i.e., step 2)
            :r_tol_drag => 1e-8,                                       # Relative tolerance for drag passage integration (inside atmosphere, i.e., step 2)
            :a_tol_quaternion => 1e-8,                                  # Absolute tolerance for quaternion integration (inside atmosphere, i.e., step 2)
            :r_tol_quaternion => 1e-5,                                  # Relative tolerance for quaternion integration (inside atmosphere, i.e., step 2)
            :dt_max => 1.0,                                         # Maximum time step for integration, s
            :dt_max_orbit => 10.0,                                   # Maximum time step for orbit integration (outside atmosphere, i.e., step 1 and step 3), s
            :dt_max_drag => 0.1,                                    # Maximum time step for drag passage

            :Odyssey_sim => 0                                       # Simulate Odyssey Mission
            )


    t = @elapsed begin
        # Run the simulation
        sol = run_analysis(args)

        if Bool(args[:passresults])
            println("Ra initial = " * string((sol.orientation.oe[1][1] * (1 + sol.orientation.oe[2][1]))* 1e-3) * " km, Ra new = " * string((sol.orientation.oe[1][end] * (1 + sol.orientation.oe[2][end]))* 1e-3) * " km - Actual periapsis altitude = " * string(minimum(sol.orientation.alt) * 1e-3) * " km - Target Ra = " * string(args[:final_apoapsis] * 1e-3) * " km")
        end
    end

    println("COMPUTATIONAL TIME = " * string(t) * " s")