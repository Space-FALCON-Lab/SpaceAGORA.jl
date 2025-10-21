include("../simulation/Run.jl")
include("../config.jl")
include("../utils/maneuver_plans.jl")

using SPICE
args = Dict(# Misc Simulation
            :results => 1,                                                                                      # Generate csv file for results True=1, False=0
            :passresults => 1,                                                                                  # Pass results as output True=1, False=0
            :print_res => 0,                                                                                    # Print some lines True=1, False=0
            :directory_results => "/workspaces/ABTS.jl/output/earth_GHSM_low_hp",            # Directory where to save the results
            :directory_Gram => "/workspaces/ABTS.jl/GRAMpy",                   # Directory where Gram is
            :directory_Gram_data => "/workspaces/ABTS.jl/GRAM_Data",           # Directory where Gram data is
            :directory_Spice => "/workspaces/ABTS.jl/GRAM_Data/SPICE",         # Directory where SPICE files are located
            :Gram_version => 0,                                                                                 # MarsGram x file to use
            :montecarlo_analysis => 0,                                                                          # Generate csv file for Montecarlo results True=1, False=0
            :plot => 0,                                                                                         # Generate pdf plots of results True=1, False=0
            :filename => 1,                                         # Filename with specifics of simulation, True =1, False=0
            :machine => "",                                         # choices=['Laptop' , 'Cluster' , 'Aero' , 'Desktop_Home','Karnap_Laptop']
            :integrator => "Julia",                                 # choices=['Costumed', 'Julia'] Costumed customed integrator, Julia DifferentialEquations.jl library integrator, only for drag passage, others phases use RK4
            :normalize => 0,                                       # Normalize the state during integration True=1, False=0
<<<<<<< HEAD
=======
            :closed_form => 0,                                     # Closed form solution for the drag passage True=1, False=0

>>>>>>> targeting_attiude_merge
            # Type of Mission
            :type_of_mission => "Orbits",                           # choices=['Drag Passage' , 'Orbits' , 'Aerobraking Campaign']
            :keplerian => 0,                                        # Do not include drag passage: True=1, False=0, NOTE: Can't be used with aerobraking campaign, must be 'Orbits' in type_of_mission
            :number_of_orbits => 100,                                 # Number of aerobraking passage

            # Physical Model
            :planet => 0,                                           # Earth = 0, Mars = 1, Venus = 2
            :planettime => 0.0,                                     # Initial time of the mission, sec. Important for J2 effect and rotation of the planet
            :gravity_model => "Inverse Squared",      # choices=['Constant' , 'Inverse Squared' , 'Inverse Squared and J2 effect']
            :density_model => "Gram",                               # choices=['No-Density', 'Constant' , 'Exponential' , 'nrlmsise', 'Gram']
            :topography_model => "Spherical Harmonics",                             # choices=['None' , 'Spherical Harmonics']
            :topography_harmonics_file => "/workspaces/ABTS.jl/Topography_harmonics_data/Earth2012.csv", # File with the topography harmonics coefficients
            :topo_degree => 50,                                     # Maximum degree of the topography harmonics (Defined in the file)
            :topo_order => 50,                                      # Maximum order of the topography harmonics (Defined in the file)

            :wind => 1,                                             # Wind calculation only if density model is Gram True=1, False=0
            :aerodynamic_model => "Mach-dependent",                 # choices=['Cd and Cl Constant' , 'Mach-dependent' , 'No-Ballistic flight with axial coefficient']: "Mach-dependent" specific for spacecraft shape, "No-Ballistic flight" specific for blunted-cone shape
            :thermal_model => "Maxwellian Heat Transfer",           # choices=['Maxwellian Heat Transfer' , 'Convective and Radiative']: "Maxwellian Heat Transfer" specific for spacecraft shape, "Convective and Radiative" specific for blunted-cone shape
            
            # Perturbations
            :n_bodies => ["Sun", "Moon"],                                        # Add names of bodies you want to simulate the gravity of to a list. Keep list empty if not required to simulate extra body gravity.
            :srp => 0,                                             # Solar Radiation Pressure True=1, False=0
            :gravity_harmonics => 1,                               # Gravity Harmonics True=1, False=0
            :gravity_harmonics_file => "/workspaces/ABTS.jl/Gravity_harmonics_data/egm96.csv",                  # Gravity Harmonics file to use
            :L => 50,                                              # Maximum degree of gravity harmonics
            :M => 50,                                              # Maximum order of gravity harmonics

            # Rates
            :trajectory_rate => 100.0,                              # Rate at which the trajectory in drag passage integrate using RK4
            :flash1_rate => 3.0,                                    # Rate at which Control Mode-1 is called
            :save_rate => 5.0,                                      # Rate at which the data trajectory are saved
            :control_in_loop => 0,                                  # Control in loop, control called during integration of trajectory, full state knowledge
            :flash2_through_integration => 0,                       # Integration of the equations of motion and lambda to define time switches and revaluation second time switch
            :struct_ctrl => 0,                                    # Structural thermal control, True=1, False=0
            :targeting_ctrl => 0,                                   # Targeting control True=1, False=0
            # Body
            :body_shape => "Spacecraft",                            # choices=['Spacecraft' , 'Blunted Cone']
            :max_heat_rate => 0.15,                                 # Max heat rate the heat rate control will start to react to
            :max_heat_load => 30.0,                                 # Max heat load the heat load control will not be overcomed
            :dry_mass => 640.0,                                     # Initial dry mass of body in kg
            :prop_mass => 200.0,                                     # Initial propellant mass of body in kg
            :reflection_coefficient => 0.9,                         # Diffuse reflection sigma =0, for specular reflection sigma = 1
            :thermal_accomodation_factor => 1.0,                    # Thermal accomodation factor, Shaaf and Chambre
            :α => 90.0,                                             # Max angle of attack of solar panels

            # Fill for Spacecraft body shape only
            :length_sat => 2.05,                                     # Length of the satellite in m
            :height_sat => 2.8,                                     # Height of the satellite in m
            :width_sat => 3.7,                                      # Width of the satellite in m
            :length_sp => 5.7,                                     # Length of the solar panels in m
            :height_sp => 1.0,                                     # Height of the solar panels in m

            # Fill for Blunted Cone body shape only
            :cone_angle => 70.0,                                    # Cone angle of the blunted cone in deg
            :base_radius => 2.65/2,                                 # Base radius of the blunted cone in m
            :nose_radius => 0.6638,                                 # Nose radius of the blunted cone in m
            
            # Engine
            :thrust => 4.0,                                         # Maximum magnitude thrust in N
            
            # Control Mode
            :control_mode => 0,                                     # Use Rotative Solar Panels Control:  False=0, Only heat rate=1, Only heat load=2, Heat rate and Heat load = 3
            :security_mode => 0,                                    # Security mode that set the angle of attack to 0 deg if predicted heat load exceed heat load limit
            :second_switch_reevaluation => 1,                       # Reevaluation of the second switch time when the time is closer to it
            
            # Initial Conditions
            :initial_condition_type => 0,                           # Initial Condition ra,hp = 0, Initial Condition v, gamma = 1
            :ra_initial_a => 56378.7978559e3, #56094.7978559e3, # 56378e3,                # Initial Apoapsis Radius for for-loop in m
            :ra_initial_b => 1e21,                               # Final Apoapsis Radius for for-loop in m
            :ra_step => 5e21,                                       # Step Apoapsis Radius for for-loop in m
            :hp_initial_a => 200590.0,#150590.14 - 8e3,#200590.0,#                                 # Initial Periapsis Altitude for for-loop in m
            :hp_initial_b => 9990000.0,                              # Final Periapsis Altitude for for-loop in m
            :hp_step => 1e20,                                 # Step Periapsis Radius for for-loop in m
            :v_initial_a => 3700.0,                                 # Initial Velocity (m/s) for for-loop if initial conditions are in v and gamma
            :v_initial_b => 5000.0,                                 # Final Velocity (m/s) for for-loop if initial conditions are in v and gamma
            :v_step => 100.0,                                       # Step Velocity (m/s) for for-loop if initial conditions are in v and gamma
            :γ_initial_a => 2.5,                                    # Initial Gamma (deg) for for-loop if initial conditions are in v and gamma
            :γ_initial_b => 7.0,                                    # Final Gamma (deg) for for-loop if initial conditions are in v and gamma
            :γ_step => 0.5,                                         # Step Gamma (deg) for for-loop if initial conditions are in v and gamma
            
            :orientation_type => 0,                                   # Initial Condition in Keplerian Elements = 0, Initial Condition in Cartesian Elements = 1
            :inclination => 89.876,                                   # Inclination Orbit, deg
            :ω => 75.505,                                              # AOP, deg
            :Ω => 104.115,                                              # RAAN, deg
            :EI => 300.0,                                           # Entry Interface, km
            :AE => 300.0,                                           # Atmospheric Exit, km
            :year => 2014,                                          # Mission year
            :month => 5,                                           # Mission month
            :day => 27,                                             # Mission day
            :hours => 5,                                           # Mission hour
            :minutes => 0,                                         # Mission minute
            :secs => 0.0,                                          # Mission second
            
            # Final Conditions
            :final_apoapsis => 7378.0e3,#35000e3,#7378.0e3,  # 4905.974818462152e3                  # Final apoapsis radius if aerobraking campaign

            # Do not change
            :heat_load_sol => 0,                                    # Heat load solution #leave it to 0 and change it only for control mode = 2:  Max energy depletaion=0, Min energy depletion=1, One switch max-min=2, One switch min-max = 3
            :thrust_control => "Aerobraking Maneuver",                              # choices=['None' , 'Aerobraking Maneuver' , 'Drag Passage Firing']
            :phi => 180.0,                                          # Thrust Angle, deg
            :delta_v => 0,                                          # Delta-v of Aerobraking Manuver,m/s
            :apoapsis_targeting => 0,                               # Apoapsis Targeting Enabled
            :ra_fin_orbit => 25000e3,                               # Target final apoapsis for the orbit, m
            :maneuver_plan => Earth_firing_plan,                    # Firing plan function
            
            # Monte Carlo Simulations
            :montecarlo => 0,                                       # Run Monte Carlo simulation True=1, False=0
            :initial_montecarlo_number => 1,                        # Initial Monte Carlo sample number
            :montecarlo_size => 1000,                               # number of Monte Carlo samples
            
            # Monte Carlo Perturbations
            :CD_dispersion => 10.0,                                 # Max dispersion of CD for Uniform Distribution, %
            :CL_dispersion => 10.0,                                 # Max dispersion of CL for Uniform Distribution, %
            :rp_dispersion => 200.590*0.05/3,                                  # Max dispersion for initial vacuum periapsis radius following uniform distribution, km
            :ra_dispersion => 56000.0*0.05/3,                                  # Max dispersion for initial apoapsis radius following uniform distribution, km
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
            :Odyssey_sim => 0                                       # Simulate Odyssey Mission
            )

# Calculating time of simulation
# for hour = 0:23
#     args[:hours] = hour
#     args[:directory_results] = "/home/space-falcon-1/Documents/ABTS.jl/output/venus_express_" * string(hour)
#     t = @elapsed begin
#         sol = run_analysis(args)
#     end
#     println("COMPUTATIONAL TIME = " * string(t) * " s")
# end

#####################
# For comparison with GMAT with different perturbations #
#####################
# gravity_models = ["Inverse Squared", "Inverse Squared and J2 effect", "Inverse Squared and J2 effect"]
# gravity_models_str = ["I", "IJ2", "IJ2GH"]
# # density_models = ["No-Density", "nrlmsise", "GRAM"]
# n_bodies = [[], ["Moon"], ["Sun"], ["Moon", "Sun"]]
# n_bodies_str = ["", "M", "S", "SM"]

# furnsh(args[:directory_Spice] * "/pck/pck00011.tpc")
# furnsh(args[:directory_Spice] * "/spk/planets/de440_GRAM.bsp")
# furnsh(args[:directory_Spice] * "/lsk/naif0012.tls")
# furnsh(args[:directory_Spice] * "/spk/planets/de440s.bsp")
# furnsh(args[:directory_Spice] * "/spk/satellites/sat441_GRAM.bsp")

# # for i in 1:length(density_models)
#     for j in 1:length(gravity_models)
#         for k in 1:length(n_bodies)
#             args[:gravity_model] = gravity_models[j]
#             if j == 3
#                 args[:gravity_harmonics] = 1
#             else
#                 args[:gravity_harmonics] = 0
#             end
#             # args[:density_model] = density_models[i]
#             args[:n_bodies] = n_bodies[k]
#             args[:directory_results] = "/workspaces/ABTS.jl/output/earth_G"*gravity_models_str[j]*n_bodies_str[k]
#             # args[:filename_results] = "Earth_Aerobraking_" * string(i) * "_" * string(j) * "_" * string(k)
#             println("Running simulation with gravity model: " * gravity_models[j] * ", n_bodies: " * string(n_bodies[k]))
#             t = @elapsed begin
#                 sol = run_analysis(args)
#             end
#             println("COMPUTATIONAL TIME = " * string(t) * " s")
#         end
#     end
# # end

mc_runs = 0
nominal_ra = args[:ra_initial_a]
nominal_rp = args[:hp_initial_a]
nominal_i = args[:inclination]
nominal_Ω = args[:Ω]
nominal_ω = args[:ω]
for i in 1:mc_runs+1
    t = @elapsed begin
        if i == 1
            args[:ra_initial_a] = nominal_ra
            args[:hp_initial_a] = nominal_rp
            args[:inclination] = nominal_i
            args[:Ω] = nominal_Ω
            args[:ω] = nominal_ω
            args[:directory_results] = "/workspaces/ABTS.jl/output/earth_MC_ab_campaign/Nominal"
            args[:print_res] = 1
        else
            # Generate random values for the parameters
            args[:ra_initial_a] = nominal_ra + randn()*sqrt(args[:ra_dispersion]) * 1e3
            args[:hp_initial_a] = nominal_rp + randn()*sqrt(args[:rp_dispersion]) * 1e3
            args[:inclination] = nominal_i + randn()*sqrt(args[:i_dispersion])
            args[:Ω] = nominal_Ω + randn()*sqrt(args[:Ω_dispersion])
            args[:ω] = nominal_ω + randn()*sqrt(args[:ω_dispersion])
            # args[:directory_results] = "/workspaces/ABTS.jl/output/earth_MC_ab_campaign/" * string(i-1)
            args[:print_res] = 0
        end
        # args[:ra_initial_a] = nominal_ra + randn()*sqrt(args[:ra_dispersion]) * 1e3
        # args[:hp_initial_a] = nominal_rp + randn()*sqrt(args[:rp_dispersion]) * 1e3
        # args[:inclination] = nominal_i + randn()*sqrt(args[:i_dispersion])
        # args[:Ω] = nominal_Ω + randn()*sqrt(args[:Ω_dispersion])
        # args[:ω] = nominal_ω + randn()*sqrt(args[:ω_dispersion])

        # Run the simulation
        sol = run_analysis(args)

        if Bool(args[:passresults])
            println("Ra initial = " * string((sol.orientation.oe[1][1] * (1 + sol.orientation.oe[2][1]))* 1e-3) * " km, Ra new = " * string((sol.orientation.oe[1][end] * (1 + sol.orientation.oe[2][end]))* 1e-3) * " km - Actual periapsis altitude = " * string(minimum(sol.orientation.alt) * 1e-3) * " km - Target Ra = " * string(args[:final_apoapsis] * 1e-3) * " km")
        end
    end

    println("COMPUTATIONAL TIME = " * string(t) * " s")
end