const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

using CSV
using DataFrames
using SPICE
using SpaceAGORA
using SpaceAGORA.SimulationModel

const SimulationModel = SpaceAGORA.SimulationModel
const quat_mult = SimulationModel.quat_mult
const run_simulation = SpaceAGORA.run_simulation

if Threads.nthreads() < 2
    error("Threaded smoke requires at least 2 Julia threads; got $(Threads.nthreads())")
end

spice_path = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
planet = Earth("", spice_path)

function make_sc(id::Int64, ν_deg::Float64)
    root = Link{0}(root=true, m=140.0, ref_area=1.2)
    ic = InitialCondition(
        ra=planet.Rp_e + 520e3,
        rp=planet.Rp_e + 500e3,
        i=28.0,
        ω=15.0,
        Ω=20.0,
        ν=ν_deg
    )
    return SpacecraftModel(
        joints=Joint[],
        links=Link[root],
        root=root,
        instant_actuation=true,
        prop_mass=15.0,
        inertia_tensor=root.inertia,
        n_reaction_wheels=0,
        n_thrusters=0,
        initial_condition=ic,
        id=id
    )
end

sc1 = make_sc(1, 170.0)
sc2 = make_sc(2, 160.0)

thruster = BaseThrusterModel(
    thrust=[0.9, 0.8],
    direction=[0.0, π],
    Δv=[0.0, 0.0],
    start_burn_time=[0.0, 0.0],
    stop_burn_time=[15.0, 15.0],
    Isp=[300.0, 280.0]
)

args = SimulationConfiguration(
    simulation_settings=SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        normalize=false
    ),
    mission_configuration=MissionConfiguration(
        mission_type=MissionTime,
        keplerian=true,
        number_of_orbits=1,
        mission_time=60.0,
        orientation_sim=false,
        num_steps_to_save=200
    ),
    environment_model=EnvironmentModel(
        planet=planet,
        EI=120.0,
        density_model=NoAtmosphereModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false,
        wind=false
    ),
    dynamics_model=DynamicsModel([sc1, sc2], (InverseSquaredGravityModel(),)),
    guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
    navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
    control_model=ControlModel(control_effectors=(thruster,), control_rates=[1.0]),
    initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
    integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
)

mktempdir() do tmp
    cd(tmp) do
        run_simulation(args)
        csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
        if !isfile(csv_path)
            error("Expected simulation_results.csv to be written by threaded smoke run")
        end
        df = CSV.read(csv_path, DataFrame)
        if nrow(df) < 6
            error("Threaded smoke produced too few rows: $(nrow(df))")
        end

        required_cols = [
            "sc1_pos_1", "sc1_vel_1", "sc1_mass",
            "sc2_pos_1", "sc2_vel_1", "sc2_mass"
        ]
        for col in required_cols
            if !(col in names(df))
                error("Missing required column $col in threaded smoke output")
            end
            if !all(isfinite, Float64.(df[!, col]))
                error("Non-finite values found in $col during threaded smoke")
            end
        end

        if !(minimum(Float64.(df.sc1_mass)) < Float64(df.sc1_mass[1]))
            error("Expected sc1 mass to decrease during threaded smoke run")
        end
        if !(minimum(Float64.(df.sc2_mass)) < Float64(df.sc2_mass[1]))
            error("Expected sc2 mass to decrease during threaded smoke run")
        end
    end
end

println("threaded_smoke_ok")

# ---------------------------------------------------------------------------
# Thread-count independence. The effector loop, the n-body loop and the
# multibody aero loop collect per-item results and sum them in a fixed order on
# one thread, so the same configuration must give bit-identical trajectories
# whether the effectors are evaluated serially or on 2 or 4 workers, and
# whether the n-body bodies are evaluated serially or on 2 workers. Three
# effectors with different costs make the partition non-trivial.
#
# The RHS execution route is pinned for every run below. By default the engine
# times the candidate routes at start-up (SPACEAGORA_RHS_CALIBRATE=auto) and
# keeps the fastest, and the batched routes (satellite_batch, flat) evaluate the
# harmonics through @fastmath/@turbo kernels that are not bit-identical to the
# per-satellite scalar path on every CPU. That choice is a wall-clock decision,
# so leaving it unpinned turned this check into a lottery on the CI runners.
# The route drift itself is measured and printed further down, not asserted.
# ---------------------------------------------------------------------------
harmonics_file = joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
det_effectors = (
    InverseSquaredGravityModel(),
    NBodyGravityModel(["Sun", "Moon"], "Earth", spice_path),
    GravitationalHarmonicsModel(4, 4, harmonics_file, planet),
    AerodynamicCoefficientfM(),   # exercises the per-link aero loop (legacy path) as well
)
det_environment = EnvironmentModel(
    planet=planet,
    EI=120.0,
    density_model=ExponentialAtmosphereModel(planet),
    thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
    topography=false,
    wind=false
)
det_args = SimulationConfiguration(
    simulation_settings=args.simulation_settings,
    mission_configuration=MissionConfiguration(
        mission_type=MissionTime,
        keplerian=true,
        number_of_orbits=1,
        mission_time=300.0,
        orientation_sim=false,
        num_steps_to_save=50
    ),
    environment_model=det_environment,
    dynamics_model=DynamicsModel([sc1, sc2], det_effectors),
    guidance_model=args.guidance_model,
    navigation_model=args.navigation_model,
    control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
    initial_time=args.initial_time,
    integration_tolerances=args.integration_tolerances
)

function det_run(env_pairs::Vector{Pair{String, String}})::Tuple{Matrix{Float64}, Vector{String}}
    return withenv(env_pairs...) do
        mktempdir() do tmp
            cd(tmp) do
                run_simulation(det_args)
                df = CSV.read(joinpath(det_args.simulation_settings.results_directory, "simulation_results.csv"), DataFrame)
                cols = [c for c in names(df) if eltype(df[!, c]) <: Real]
                return Matrix{Float64}(df[:, cols]), cols
            end
        end
    end
end

# Describe a mismatch precisely enough to act on from a CI log alone.
function det_describe_difference(out::Matrix{Float64}, ref::Matrix{Float64}, cols::Vector{String})::String
    d = abs.(out .- ref)
    d[isnan.(d)] .= Inf   # NaN in one run only counts as a difference
    bad = findall(>(0.0), vec(maximum(d; dims=1)))
    rows = unique(getindex.(findall(>(0.0), d), 1))
    io = IOBuffer()
    print(io, "max abs difference ", maximum(d), " in ", length(bad), " column(s) over ", length(rows), "/", size(d, 1), " row(s)")
    for j in bad[1:min(end, 8)]
        i = argmax(view(d, :, j))
        print(io, "; ", cols[j], "[", i, "] ref=", repr(ref[i, j]), " out=", repr(out[i, j]))
    end
    return String(take!(io))
end

det_machine = let info = Sys.cpu_info()
    "cpu=$(isempty(info) ? "unknown" : String(info[1].model)) threads=$(Threads.nthreads())"
end

det_pinned = [
    "SPACEAGORA_RHS_CALIBRATE" => "off",
]
det_force_on = [
    "SPACEAGORA_EFFECTOR_PARALLEL_HEAVY_ONLY" => "0",
    "SPACEAGORA_EFFECTOR_THREAD_THRESHOLD" => "1",
    "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "1",
]
det_serial_route = vcat(det_pinned, ["SPACEAGORA_RHS_EXECUTION_MODE" => "serial"])
det_per_sat_route = vcat(det_pinned, ["SPACEAGORA_RHS_EXECUTION_MODE" => "per_satellite"])
det_ref, det_cols = det_run(vcat(det_serial_route, ["SPACEAGORA_EFFECTOR_PARALLEL" => "off", "SPACEAGORA_MULTIBODY_PARALLEL" => "off"]))
size(det_ref, 1) >= 10 || error("Determinism check produced too few rows: $(size(det_ref, 1))")
det_variants = [
    "effectors on 2 workers" => vcat(det_per_sat_route, det_force_on, ["SPACEAGORA_EFFECTOR_PARALLEL" => "on", "SPACEAGORA_EFFECTOR_MAX_THREADS" => "2", "SPACEAGORA_MULTIBODY_PARALLEL" => "off"]),
    "effectors on 4 workers" => vcat(det_per_sat_route, det_force_on, ["SPACEAGORA_EFFECTOR_PARALLEL" => "on", "SPACEAGORA_EFFECTOR_MAX_THREADS" => "4", "SPACEAGORA_MULTIBODY_PARALLEL" => "off"]),
    "n-body on 2 workers"    => vcat(det_serial_route, det_force_on, ["SPACEAGORA_EFFECTOR_PARALLEL" => "off", "SPACEAGORA_MULTIBODY_PARALLEL" => "on", "SPACEAGORA_MULTIBODY_MAX_THREADS" => "2"]),
    "everything on"          => vcat(det_per_sat_route, det_force_on, ["SPACEAGORA_EFFECTOR_PARALLEL" => "on", "SPACEAGORA_EFFECTOR_MAX_THREADS" => "4", "SPACEAGORA_MULTIBODY_PARALLEL" => "on", "SPACEAGORA_MULTIBODY_MAX_THREADS" => "2"]),
]
for (label, env_pairs) in det_variants
    out, _ = det_run(env_pairs)
    if size(out) != size(det_ref)
        error("Threaded determinism: $label produced $(size(out)) rows/cols, serial produced $(size(det_ref)) [$det_machine]")
    end
    if !isequal(out, det_ref)
        error("Threaded determinism: $label differs from the serial run: $(det_describe_difference(out, det_ref, det_cols)) [$det_machine]")
    end
end
println("threaded_determinism_ok variants=$(length(det_variants)) rows=$(size(det_ref, 1)) route=per_satellite $det_machine")

# Route drift report. The batched routes may legitimately differ from the
# scalar route at the last bit on a given CPU; print how much so the number is
# on record for every runner, without failing the job.
det_routes = [
    # Two spacecraft sit below SPACEAGORA_RHS_BATCH_THREAD_THRESHOLD (default 16),
    # so the batch loop only runs when the switch is forced on.
    "satellite_batch" => vcat(det_pinned, ["SPACEAGORA_RHS_EXECUTION_MODE" => "satellite", "SPACEAGORA_RHS_BATCH_PARALLEL" => "on", "SPACEAGORA_EFFECTOR_PARALLEL" => "off", "SPACEAGORA_MULTIBODY_PARALLEL" => "off"]),
    "flat"            => vcat(det_pinned, det_force_on, ["SPACEAGORA_RHS_EXECUTION_MODE" => "flat", "SPACEAGORA_EFFECTOR_PARALLEL" => "on", "SPACEAGORA_EFFECTOR_MAX_THREADS" => "4", "SPACEAGORA_MULTIBODY_PARALLEL" => "off"]),
]
for (label, env_pairs) in det_routes
    out, _ = det_run(env_pairs)
    if size(out) != size(det_ref)
        println("threaded_route_report route=$label result=size_mismatch $(size(out)) vs $(size(det_ref)) [$det_machine]")
    elseif isequal(out, det_ref)
        println("threaded_route_report route=$label result=identical [$det_machine]")
    else
        println("threaded_route_report route=$label result=differs $(det_describe_difference(out, det_ref, det_cols)) [$det_machine]")
    end
end
