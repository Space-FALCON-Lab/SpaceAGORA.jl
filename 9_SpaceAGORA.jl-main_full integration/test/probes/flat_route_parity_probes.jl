# The flat constellation route must reproduce the serial route bit for bit:
# per-effector slots summed in effector order, and one compiled kernel per
# effector shared by every route. This probe runs the determinism smoke's
# two-satellite plunge through both routes with a mixed effector set and
# compares every saved column exactly.
using Test
using CSV
using DataFrames
using SpaceAGORA
using SpaceAGORA.SimulationModel

const REPO_ROOT_PARITY = normpath(joinpath(@__DIR__, "..", ".."))
const SPICE_PATH_PARITY = joinpath(REPO_ROOT_PARITY, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")

@testset "Flat constellation route matches the serial route bit for bit" begin
    planet = Earth("", SPICE_PATH_PARITY)
    harmonics_file = joinpath(REPO_ROOT_PARITY, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
    function make_sc(id::Int64, ν_deg::Float64)
        root = Link(root=true, m=140.0, ref_area=1.2)
        ic = InitialCondition(ra=planet.Rp_e + 520e3, rp=planet.Rp_e + 500e3, i=28.0, ω=15.0, Ω=20.0, ν=ν_deg)
        return SpacecraftModel(
            joints=Joint[], links=Link[root], root=root, instant_actuation=true, prop_mass=15.0,
            inertia_tensor=root.inertia, n_reaction_wheels=0, n_thrusters=0, initial_condition=ic, id=id,
        )
    end
    effectors = (
        InverseSquaredJ2GravityModel(),
        AerodynamicCoefficientfM(),
        NBodyGravityModel(["Sun", "Moon"], "Earth", SPICE_PATH_PARITY),
        SolarRadiationPressureModel(1.2, 12.0; direct=true, albedo=true, ir=true),
        GravitationalHarmonicsModel(6, 6, harmonics_file, planet),
        InverseSquaredGravityModel(),
    )
    cfg = SimulationConfiguration(
        simulation_settings=SimulationSettings(results=true, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime, keplerian=true, number_of_orbits=1, mission_time=300.0, orientation_sim=false, num_steps_to_save=50),
        environment_model=EnvironmentModel(planet=planet, EI=120.0, density_model=ExponentialAtmosphereModel(planet), thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet), topography=false, wind=false),
        dynamics_model=DynamicsModel([make_sc(1, 170.0), make_sc(2, 182.0)], effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0),
    )
    function run_route(env_pairs)
        withenv(env_pairs...) do
            mktempdir() do tmp
                cd(tmp) do
                    run_simulation(cfg)
                    df = CSV.read(joinpath(cfg.simulation_settings.results_directory, "simulation_results.csv"), DataFrame)
                    cols = [c for c in names(df) if eltype(df[!, c]) <: Real]
                    return Matrix{Float64}(df[:, cols]), cols
                end
            end
        end
    end
    pinned = ["SPACEAGORA_RHS_CALIBRATE" => "off", "SPACEAGORA_EFFECTOR_PARALLEL_HEAVY_ONLY" => "0",
              "SPACEAGORA_EFFECTOR_THREAD_THRESHOLD" => "1", "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "1"]
    ref, cols = run_route(vcat(pinned, ["SPACEAGORA_RHS_EXECUTION_MODE" => "serial", "SPACEAGORA_EFFECTOR_PARALLEL" => "off", "SPACEAGORA_MULTIBODY_PARALLEL" => "off"]))
    @test size(ref, 1) >= 10
    for (label, env_pairs) in (
        "flat, 2 workers" => vcat(pinned, ["SPACEAGORA_RHS_EXECUTION_MODE" => "flat", "SPACEAGORA_EFFECTOR_PARALLEL" => "on", "SPACEAGORA_EFFECTOR_MAX_THREADS" => "2", "SPACEAGORA_MULTIBODY_PARALLEL" => "off"]),
        "flat, 4 workers, n-body threaded" => vcat(pinned, ["SPACEAGORA_RHS_EXECUTION_MODE" => "flat", "SPACEAGORA_EFFECTOR_PARALLEL" => "on", "SPACEAGORA_EFFECTOR_MAX_THREADS" => "4", "SPACEAGORA_MULTIBODY_PARALLEL" => "on", "SPACEAGORA_MULTIBODY_MAX_THREADS" => "2"]),
        "satellite batch" => vcat(pinned, ["SPACEAGORA_RHS_EXECUTION_MODE" => "satellite", "SPACEAGORA_RHS_BATCH_PARALLEL" => "on", "SPACEAGORA_EFFECTOR_PARALLEL" => "off", "SPACEAGORA_MULTIBODY_PARALLEL" => "off"]),
    )
        out, _ = run_route(env_pairs)
        @test size(out) == size(ref)
        if size(out) == size(ref)
            d = abs.(out .- ref)
            bad = findall(>(0.0), vec(maximum(d; dims=1)))
            @test isempty(bad)
            isempty(bad) || @info "route differs" label columns=cols[bad][1:min(end, 8)] max=maximum(d)
        end
    end
end
