# The complete script of "Adding a Force or Torque of Your Own"
# (docs/src/user/custom_effector.md), run by CI.
# The region between `# snippet: script` and its `# end snippet` is the page's
# script, and each nested region one of its shorter excerpts; the page shows each
# without the marker lines, and the HPC and extensibility docs gate compares them.
# snippet: script
using SpaceAGORA
using StaticArrays
using LinearAlgebra
using CSV, DataFrames

# --- 1. The effector -------------------------------------------------------------
# snippet: effector
struct AlongTrackThrust <: SpaceAGORA.AbstractForceTorqueModel
    force_n::Float64          # thrust magnitude, newtons, always along the velocity
end

function SpaceAGORA.wrench(model::AlongTrackThrust, x::SpaceAGORA.StateSample,
                           env::SpaceAGORA.EnvironmentSample, t::Float64)
    direction = x.vel_ii / norm(x.vel_ii)          # unit vector along the inertial velocity
    force_ii = model.force_n * direction            # newtons, inertial frame
    torque_body = SVector{3, Float64}(0.0, 0.0, 0.0)  # newton-metres, body frame
    return force_ii, torque_body
end
# end snippet

# --- 2. A simulation that uses it ---------------------------------------------------
planet = make_no_gram_planet(:earth)
altitude_m = 400e3
a0 = planet.Rp_e + altitude_m
period_s = 2π * sqrt(a0^3 / planet.μ)

function build_config(effectors; results_directory::String)
    ic = InitialCondition(a=a0, e=0.001, i=51.6, ω=0.0, Ω=0.0, ν=0.0)
    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.0, 2.0, 2.0), panel_dims=(0.02, 2.0, 1.0),
        bus_mass=500.0, panel_mass_each=10.0, panel_offset_y=2.0,
        ic=ic, prop_mass=0.0, id=1)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=true, verbose=false, results_directory=results_directory,
            generate_plots=false, generate_filenames=false, normalize=false, save_csv=true),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime, mission_time=period_s, number_of_orbits=1,
            keplerian=true, orientation_sim=false, num_steps_to_save=1000, data_rate=10.0),
        environment_model=EnvironmentModel(
            planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
            ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false),
        dynamics_model=DynamicsModel([spacecraft], effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances())
end

thrust = AlongTrackThrust(0.052)     # 0.052 N on 520 kg: 1e-4 m/s^2
output_root = joinpath(@__DIR__, "output", "effector_walkthrough")
config_without = build_config((InverseSquaredGravityModel(),); results_directory=joinpath(output_root, "without"))
config_with    = build_config((InverseSquaredGravityModel(), thrust); results_directory=joinpath(output_root, "with"))

results_without = run_simulation(config_without; return_results=true)
results_with = run_simulation(config_with; return_results=true)
results_with.table          # named state columns, also written to the CSV
results_with.files          # files this call wrote
results_with.configuration  # isolated configuration that actually ran

# --- 3. Verify the effect --------------------------------------------------------
function semi_major_axis(csv_path)
    df = CSV.read(csv_path, DataFrame)
    last = df[end, :]
    r = norm((last.sc1_pos_1, last.sc1_pos_2, last.sc1_pos_3))
    v = norm((last.sc1_vel_1, last.sc1_vel_2, last.sc1_vel_3))
    energy = v^2 / 2 - planet.μ / r          # specific orbital energy, J/kg
    return -planet.μ / (2 * energy), last.time, last.sc1_mass
end

a_without, t_end, mass_kg = semi_major_axis(joinpath(output_root, "without", "simulation_results.csv"))
a_with, _, _ = semi_major_axis(joinpath(output_root, "with", "simulation_results.csv"))
accel = thrust.force_n / mass_kg
v0 = sqrt(planet.μ / a0)
expected_da = 2 * a0^2 / planet.μ * accel * v0 * t_end
observed_da = a_with - a_without
println("mass used by the engine: $(mass_kg) kg, thrust acceleration: $(accel) m/s^2")
println("semi-major axis after one orbit without thrust: $(round(a_without / 1e3, digits=3)) km")
println("semi-major axis after one orbit with thrust:    $(round(a_with / 1e3, digits=3)) km")
println("expected rise: $(round(expected_da, digits=1)) m, observed rise: $(round(observed_da, digits=1)) m, ratio $(round(observed_da / expected_da, digits=3))")
# end snippet

using Test
# CI's checks of the page's claims, outside the region the page shows. The ratio to
# the first-order estimate stays a diagnostic, as the page says; CI checks that both
# runs reach the end of the orbit, that the engine uses the 520 kg vehicle, and that
# the thrust raises the orbit.
@test t_end ≈ period_s rtol=1e-6
@test mass_kg == 520.0
@test observed_da > 0
