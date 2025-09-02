include("../simulation/Run.jl")
# include("../config.jl") #TODO:Figure out how to run multiple times without having to comment this line out
include("../utils/maneuver_plans.jl")
include("../utils/attitude_control_plans.jl")
include("../utils/quaternion_utils.jl")
# include("SpacecraftModel.jl")

import .config
import .ref_sys
using Profile
# import .SpacecraftModel
# Define spacecraft model
spacecraft = config.SpacecraftModel()
# Add bodies to the spacecraft model
p = SVector{3, Float64}([0.1, 0.2, -0.3])
q = 1/(1+norm(p)^2)*SVector{4, Float64}([2*p; 1-norm(p)^2])
skew = (ω) -> SMatrix{3, 3, Float64}([0 -ω[3] ω[2];
                                   ω[3] 0 -ω[1];
                                   -ω[2] ω[1] 0])
dcm = (q[4]^2 - norm(q[1:3])^2)*I(3) - 2*q[4]*skew(q[1:3]) + 2*q[1:3]*q[1:3]' # DCM from quaternion
ω_body = SVector{3, Float64}([0.001, -0.01, 0.03]) # Reference angular velocity
ω_ref = dcm'*ω_body
h = sqrt(30.0/7.0)
w = sqrt(6.0)
d = sqrt(66.0/7.0)

# q = SVector{4, Float64}([0.0, 0.0, sin(pi/4), cos(pi/4)]) # Quaternion for the main bus
main_bus = config.Link(root=true, 
                        r=SVector{3, Float64}(0.0, 0.0, 0.0), 
                        q=SVector{4, Float64}(q),
                        # q=SVector{4, Float64}([0.0, 0.0, 0.0, 1.0]),
                        ṙ=SVector{3, Float64}([0.0, 0.0, 0.0]), 
                        ω=SVector{3, Float64}(ω_body),
                        dims=SVector{3, Float64}([1.5, 1.8, 2.86]), 
                        ref_area=1.5*2.86,
                        m=200.0, 
                        gyro=0)

config.add_body!(spacecraft, main_bus, prop_mass=1.0)
# config.add_body!(spacecraft, L_panel)
# config.add_body!(spacecraft, R_panel)

# L_panel_joint = config.Joint(main_bus, L_panel)
# R_panel_joint = config.Joint(R_panel, main_bus)
# config.add_joint!(spacecraft, L_panel_joint)
# config.add_joint!(spacecraft, R_panel_joint)
config.set_inertia_tensor!(spacecraft, main_bus, 
                        SMatrix{3, 3, Float64}(Diagonal([900.0, 800.0, 600.0])))
lenXHub = 1.50 # m
lenYHub = 1.80 # m
lenZHub = 2.86 # m
array_width = 7.262 # m
area2 = pi*array_width/2*array_width/2

println("Spacecraft model initialized with $(length(spacecraft.links)) bodies.")
# println("Spacecraft roots: $spacecraft.roots")
println("Spacecraft COM: $(config.get_COM(spacecraft, main_bus))")
println("Spacecraft MOI: $(config.get_inertia_tensor(spacecraft, main_bus))")

args = Dict(# Misc Simulation
            :results => 1,                                                                                      # Generate csv file for results True=1, False=0
            :passresults => 1,                                                                                  # Pass results as output True=1, False=0
            :print_res => 1,                                                                                    # Print some lines True=1, False=0
            :directory_results => "/workspaces/ABTS.jl/output/basilisk_torque_free_comparison",                # Directory where to save the results
            :directory_Gram => "/workspaces/ABTS.jl/GRAMpy",                                                    # Directory where Gram is
            :directory_Gram_data => "/workspaces/ABTS.jl/GRAM_Data",                                            # Directory where Gram data is
            :directory_Spice => "/workspaces/ABTS.jl/GRAM_Data/SPICE",                                          # Directory where SPICE files are located
            :Gram_version => 0,                                                                                 # MarsGram x file to use
            :montecarlo_analysis => 0,                                                                          # Generate csv file for Montecarlo results True=1, False=0
            :plot => 1,                                                                                         # Generate pdf plots of results True=1, False=0
            :filename => 1,                                         # Filename with specifics of simulation, True =1, False=0
            :machine => "",                                         # choices=['Laptop' , 'Cluster' , 'Aero' , 'Desktop_Home','Karnap_Laptop']
            :integrator => "Julia",                                 # choices=['Costumed', 'Julia'] Costumed customed integrator, Julia DifferentialEquations.jl library integrator, only for drag passage, others phases use RK4
            :normalize => 1,                                       # Normalize the integration True=1, False=0
            :closed_form => 0,                                     # Closed form solution True=1, False=0
            # Type of Mission
            :type_of_mission => "Time",                           # choices=['Drag Passage' , 'Orbits' , 'Aerobraking Campaign']
            :keplerian => 1,                                        # Do not include drag passage: True=1, False=0
            :number_of_orbits => 10,                                 # Number of aerobraking passage
            :mission_time => 1000000.0,                                  # Mission time in seconds, used only for Time mission type
            :orientation_sim => true,                                  # Orientation simulation True=1, False=0, if false, will only propagate position

            # Physical Model
            :planet => 0,                                           # Earth = 0, Mars = 1, Venus = 2
            :planettime => 0.0,                                     # Initial time of the mission, sec. Important for J2 effect and rotation of the planet
            :gravity_model => "Inverse Squared",      # choices=['Constant' , 'Inverse Squared' , 'Inverse Squared and J2 effect', 'GRAM']
            :density_model => "Gram",                               # choices=['Constant' , 'Exponential' , 'Gram']
            :topography_model => "None",                             # choices=['None' , 'Spherical Harmonics']
            :topography_harmonics_file => "/workspaces/ABTS.jl/Topography_harmonics_data/Earth2012.csv", # File with the topography harmonics coefficients
            :topo_degree => 90,                                     # Maximum degree of the topography harmonics (Defined in the file)
            :topo_order => 90,                                      # Maximum order of the topography harmonics (Defined in the file)
            :wind => 1,                                             # Wind calculation only if density model is Gram True=1, False=0
            :aerodynamic_model => "Mach-dependent",                 # choices=['Cd and Cl Constant' , 'Mach-dependent' , 'No-Ballistic flight with axial coefficient']: "Mach-dependent" specific for spacecraft shape, "No-Ballistic flight" specific for blunted-cone shape
            :thermal_model => "Maxwellian Heat Transfer",           # choices=['Maxwellian Heat Transfer' , 'Convective and Radiative']: "Maxwellian Heat Transfer" specific for spacecraft shape, "Convective and Radiative" specific for blunted-cone shape
            
            # Perturbations
            :n_bodies => [],                                        # Add names of bodies you want to simulate the gravity of to a list. Keep list empty if not required to simulate extra body gravity.
            :srp => false,                                             # Solar Radiation Pressure true/false
            :eclipse => false,                                         # Whether to include eclipse conditions in SRP calculation
            :gravity_gradient => false,                                   # Gravity Gradient true/false
            :gravity_harmonics => 0,                                            # Gravity Spherical harmonics True=1, False=0
            :gravity_harmonics_file => "/workspaces/ABTS.jl/Gravity_harmonics_data/EarthGGM05C.csv", # File with the gravity harmonics coefficients
            :L => 50,                                              # Maximum degree of the gravity harmonics (Defined in the file)
            :M => 50,                                              # Maximum order of the gravity harmonics (Defined in the file)

            # Rates
            :trajectory_rate => 100.0,                              # Rate at which the trajectory in drag passage integrate using RK4
            :flash1_rate => 3.0,                                    # Rate at which Control Mode-1 is called
            :save_rate => 5.0,                                      # Rate at which the data trajectory are saved
            
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
            
            # Initial Conditions
            :initial_condition_type => 2,                           # Initial Condition ra,hp = 0, Initial Condition v, gamma = 1
            :ra_initial_a => 15000.0e3,                # Initial Apoapsis Radius for for-loop in m
            :ra_initial_b => 50000e3,                               # Final Apoapsis Radius for for-loop in m
            :ra_step => 5e10,                                       # Step Apoapsis Radius for for-loop in m
            :hp_initial_a => 145.0e3,                                 # Initial Periapsis Altitude for for-loop in m
            :hp_initial_b => 1590000.0e3,                              # Final Periapsis Altitude for for-loop in m
            :hp_step => 1e12,                                 # Step Periapsis Radius for for-loop in m
            :v_initial_a => 4500.0,                                 # Initial Velocity (m/s) for for-loop if initial conditions are in v and gamma
            :v_initial_a => 4500.0,                                 # Initial Velocity (m/s) for for-loop if initial conditions are in v and gamma
            :v_initial_b => 5000.0,                                 # Final Velocity (m/s) for for-loop if initial conditions are in v and gamma
            :v_step => 1000.0,                                       # Step Velocity (m/s) for for-loop if initial conditions are in v and gamma
            :a_initial_a => 10000.0e3,                # Initial Semi-major axis for for-loop in m
            :a_initial_b => 10001.0e3,                               # Final Semi-major axis for for-loop in m
            :a_step => 5e10,                                       # Step Semi-major axis for for-loop in m
            :e_initial_a => 0.1,                                   # Initial Eccentricity for for-loop in m
            :e_initial_b => 0.11,                                   # Final Eccentricity for for-loop in m
            :e_step => 0.1,                                       # Step Eccentricity for for-loop in m
            
            :orientation_type => 0,                                   # Initial Condition orientation = 0, Initial Condition orientation and velocity = 1
            :γ_initial_a => -2.5,                                    # Initial Gamma (deg) for for-loop if initial conditions are in v and gamma
            :γ_initial_b => 7.0,                                    # Final Gamma (deg) for for-loop if initial conditions are in v and gamma
            :γ_step => 100,                                         # Step Gamma (deg) for for-loop if initial conditions are in v and gamma
            :inclination => 33.3,                                   # Inclination Orbit, deg
            :ω => 347.8,                                              # AOP, deg
            :Ω => 48.2,                                              # RAAN, deg
            :ν => 85.3,                                               # True Anomaly, deg
            :EI => 160.0,                                           # Entry Interface, km
            :AE => 160.0,                                           # Atmospheric Exit, km
            :year => 2001,                                          # Mission year
            :month => 11,                                           # Mission month
            :day => 6,                                             # Mission day
            :hours => 19,                                           # Mission hour
            :minutes => 0,                                         # Mission minute
            :secs => 0.0,                                          # Mission second
            
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
            :a_tol_orbit => 1e-11,                                    # Absolute tolerance for orbit integration (outside atmosphere, i.e., step 1 and step 3)
            :r_tol_orbit => 1e-9,                                    # Relative tolerance for orbit integration (outside atmosphere, i.e., step 1 and step 3)
            :a_tol_drag => 1e-8,                                       # Absolute tolerance for drag passage integration (inside atmosphere, i.e., step 2)
            :r_tol_drag => 1e-6,                                       # Relative tolerance for drag passage integration (inside atmosphere, i.e., step 2)
            :a_tol_quaternion => 1e-11,                                  # Absolute tolerance for quaternion integration (inside atmosphere, i.e., step 2)
            :r_tol_quaternion => 1e-9,                                  # Relative tolerance for quaternion integration (inside atmosphere, i.e., step 2)
            :dt_max => 1.0,                                         # Maximum time step for integration, s
            :dt_max_orbit => 10.0,                                   # Maximum time step for orbit integration (outside atmosphere, i.e., step 1 and step 3), s
            :dt_max_drag => 1.0,                                    # Maximum time step for drag passage

            :Odyssey_sim => 0                                      # Simulate Odyssey Mission
            )

# # Calculating time of simulation
# @profview run_analysis(args)
t = @elapsed begin          
    # Run the simulation
    sol = run_analysis(args)
    if Bool(args[:passresults])
        println("Ra initial = " * string((sol.orientation.oe[1][1] * (1 + sol.orientation.oe[2][1]))* 1e-3) * " km, Ra new = " * string((sol.orientation.oe[1][end] * (1 + sol.orientation.oe[2][end]))* 1e-3) * " km - Actual periapsis altitude = " * string(minimum(sol.orientation.alt) * 1e-3) * " km - Target Ra = " * string(args[:final_apoapsis] * 1e-3) * " km")
    end
end

println("COMPUTATIONAL TIME = " * string(t) * " s")