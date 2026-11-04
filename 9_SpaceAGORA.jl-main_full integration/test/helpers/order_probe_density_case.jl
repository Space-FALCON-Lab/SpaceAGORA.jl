# A pure-Julia stand-in for GRAM's perturbed winds, used by
# gram_wind_determinism_tests.jl and the child processes it starts. Defines
# types and functions only; the including module supplies `SM`, `EM` and `CB`.
#
# Every query advances a shared counter and the returned density depends on
# it, so, like GRAM's random walk, a value depends on how many queries came
# before it and in what order. The dependence is tiny (1e-12 relative), so the
# step-size controller does not react to it, but it is visible bit for bit.
mutable struct OrderProbeModel <: SM.AbstractDensityModel
    calls::Int
    lock::ReentrantLock
end
OrderProbeModel() = OrderProbeModel(0, ReentrantLock())

# Thread ids that issued queries, across every copy of the model (runs deep
# copy their configuration).
const ORDER_PROBE_THREADS = Set{Int}()
const ORDER_PROBE_THREADS_LOCK = ReentrantLock()

EM.density_model_history_dependent(::OrderProbeModel) = true
EM.reset_density_model_history!(m::OrderProbeModel) = (lock(() -> m.calls = 0, m.lock); true)
CB.density_model_threadsafe(::OrderProbeModel) = true
CB.density_model_work_is_heavy(::OrderProbeModel) = true

function EM.getDensity(m::OrderProbeModel, h::Float64, lat::Float64, lon::Float64,
        t::Float64, wind::Bool, p)
    k = lock(m.lock) do
        m.calls += 1
    end
    lock(() -> push!(ORDER_PROBE_THREADS, Threads.threadid()), ORDER_PROBE_THREADS_LOCK)
    rho = 1.0e-9 * exp(-(h - 180.0e3) / 30.0e3) * (1.0 + 1.0e-12 * (k % 997))
    return rho, 700.0, SVector{3, Float64}(5.0, -3.0, 0.0)
end

function order_probe_args(; n::Int=8, mission_time::Float64=60.0)
    planet = SM.Earth()
    spacecraft = SM.SpacecraftModel[]
    for i in 1:n
        root = SM.Link(root=true, m=500.0, ref_area=12.0)
        ic = SM.InitialCondition(ra=planet.Rp_e + 400e3, rp=planet.Rp_e + 150e3,
            i=35.0, ω=40.0, Ω=360.0 * (i - 1) / n, ν=-5.0)
        push!(spacecraft, SM.SpacecraftModel(SM.Joint[], SM.Link[root], root, true,
            root.m, 0.0, root.inertia, 0, 0, ic, i))
    end
    return SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false,
            generate_plots=false, normalize=false, save_csv=false),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime,
            keplerian=true, number_of_orbits=1, mission_time=mission_time,
            orientation_sim=false, num_steps_to_save=10, data_rate=10.0),
        environment_model=SM.EnvironmentModel(planet=planet, EI=600.0,
            density_model=OrderProbeModel(),
            thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=true, ephemerides_model=SM.SimpleEphemeridesModel()),
        dynamics_model=SM.DynamicsModel(spacecraft,
            (SM.InverseSquaredGravityModel(), SM.AerodynamicCoefficientfM())),
        guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=SM.ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=SM.IntegrationTolerances(reltol_orbit=1e-9,
            abstol_orbit=1e-9, dt_max_orbit=5.0))
end

function order_probe_final(args)
    empty!(ORDER_PROBE_THREADS)
    sol = SpaceAGORA.run_simulation(args; return_solution=true)
    return (t=copy(sol.t), u=collect(sol.u[end]), threads=sort!(collect(ORDER_PROBE_THREADS)))
end
