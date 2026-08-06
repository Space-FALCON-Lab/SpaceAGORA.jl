using Test
using LinearAlgebra
using StaticArrays

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
using .SimulationModel

include(joinpath(REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
const run_simulation = SimulationEngine.run_simulation

const EARTH = make_no_gram_planet(:earth)

# Custom drag-style effector shaped like an extension's: a wrench model that
# declares its atmosphere dependence through the public requirements hook
# rather than being one of the built-in AerodynamicCoefficient* types.
struct ProbeCustomDragModel <: SimulationModel.AbstractForceTorqueModel end

SimulationModel.environment_requirements(::ProbeCustomDragModel) =
    SimulationModel.EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)

function SimulationModel.wrench(
    ::ProbeCustomDragModel,
    x::SimulationModel.StateSample,
    env::SimulationModel.EnvironmentSample,
    t::Float64,
)
    return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end

# ── Config builders ──────────────────────────────────────────────────────────

function make_probe_spacecraft(alt_km::Float64)
    root = Link{0}(root=true, m=100.0, ref_area=10.0)
    ic = InitialCondition(
        ra=EARTH.Rp_e + alt_km * 1e3,
        rp=EARTH.Rp_e + alt_km * 1e3,
        i=30.0,
        ω=0.0,
        Ω=0.0,
        ν=0.0
    )
    return SpacecraftModel(Joint[], [root], root, true, 100.0, 0.0, root.inertia, 0, 0, ic, 1)
end

function make_drag_probe_config(;
    density_model,
    alt_km::Float64,
    EI_km::Float64,
    solver_config::Union{Nothing, SolverConfig}=nothing,
    dynamic_effectors::Tuple=(InverseSquaredJ2GravityModel(), AerodynamicCoefficientfM()),
    mission_time::Float64=900.0
)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=false,
            number_of_orbits=1,
            mission_time=mission_time,
            orientation_sim=false,
            num_steps_to_save=1000
        ),
        environment_model=EnvironmentModel(
            planet=EARTH,
            EI=EI_km,
            density_model=density_model,
            ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=false,
            wind=false
        ),
        dynamics_model=DynamicsModel([make_probe_spacecraft(alt_km)], dynamic_effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(year=2019, month=6, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-8, abstol_orbit=1e-8, dt_max_orbit=20.0,
            reltol_atmosphere=1e-8, abstol_atmosphere=1e-8, dt_max_atmosphere=0.2
        ),
        solver_config=solver_config
    )
end

function final_state(cfg)
    sol = run_simulation(cfg; return_solution=true)
    @test string(sol.retcode) == "Success"
    sc = sol.u[end].sc[1]
    pos = SVector{3, Float64}(sc.pos[1], sc.pos[2], sc.pos[3])
    vel = SVector{3, Float64}(sc.vel[1], sc.vel[2], sc.vel[3])
    return pos, vel, Float64(sol.t[end])
end

specific_energy(pos, vel) = 0.5 * dot(vel, vel) - EARTH.μ / norm(pos)

const PROBE_RHO = 1.0e-9

drag_model() = ConstantDensityModel(density_kg_m3=PROBE_RHO, temperature_k=900.0)

# Drag must both separate the trajectory from the vacuum twin and dissipate
# specific orbital energy relative to it.
function assert_drag_engaged(; alt_km::Float64, EI_km::Float64, solver_config=nothing)
    pos_vac, vel_vac, t_vac = final_state(make_drag_probe_config(
        density_model=NoAtmosphereModel(), alt_km=alt_km, EI_km=EI_km, solver_config=solver_config
    ))
    pos_atm, vel_atm, t_atm = final_state(make_drag_probe_config(
        density_model=drag_model(), alt_km=alt_km, EI_km=EI_km, solver_config=solver_config
    ))
    @test t_vac == t_atm
    @test norm(pos_vac - pos_atm) > 100.0
    @test specific_energy(pos_atm, vel_atm) < specific_energy(pos_vac, vel_vac)
    return nothing
end

@testset "EI-partition drag engagement probes" begin
    @testset "circular LEO entirely below EI decays" begin
        # The regression this pins: a non-vacuum density model with the whole
        # orbit inside the entry interface must produce continuous drag — no
        # downward EI crossing is required to engage aero.
        assert_drag_engaged(alt_km=250.0, EI_km=300.0)
        assert_drag_engaged(alt_km=250.0, EI_km=300.0, solver_config=SolverConfig(solver_mode=:split_imex))
    end

    @testset "orbit entirely above EI decays on the split/implicit path" begin
        # The split/implicit partition used to zero aero strictly above EI for
        # every density model, so an orbit that never crossed EI downward flew
        # drag-free. Non-vanishing density models must now stay engaged.
        assert_drag_engaged(alt_km=400.0, EI_km=120.0, solver_config=SolverConfig(solver_mode=:split_imex))
        # The default path always engaged drag above EI; keep that pinned too.
        assert_drag_engaged(alt_km=400.0, EI_km=120.0)
    end

    @testset "vacuum model keeps the implicit fast path" begin
        # NoAtmosphereModel is the only model that may take the above-EI
        # shortcut; the trait is the gate.
        @test SimulationModel.EnvironmentModels.density_vanishes_above_entry_interface(NoAtmosphereModel())
        @test !SimulationModel.EnvironmentModels.density_vanishes_above_entry_interface(drag_model())
        @test !SimulationModel.EnvironmentModels.density_vanishes_above_entry_interface(ExponentialAtmosphereModel(EARTH))
        @test !SimulationModel.EnvironmentModels.density_vanishes_above_entry_interface(NRLMSISE00AtmosphereModel())
    end

    @testset "custom atmosphere-consuming effector suppresses the warning" begin
        # Extensions implement drag as custom wrench effectors that declare
        # environment_requirements(model).atmosphere = true; the diagnostic
        # must recognize the public hook, not just the built-in aero types.
        # Runs BEFORE the positive warning test so the @warn call site's
        # maxlog=1 budget is provably still available there.
        cfg_custom = make_drag_probe_config(
            density_model=drag_model(),
            alt_km=250.0,
            EI_km=300.0,
            dynamic_effectors=(InverseSquaredJ2GravityModel(), ProbeCustomDragModel()),
            mission_time=60.0
        )
        @test_logs min_level=Base.CoreLogging.Warn run_simulation(cfg_custom)
    end

    @testset "density model without aero effector warns loudly" begin
        cfg_no_aero = make_drag_probe_config(
            density_model=drag_model(),
            alt_km=250.0,
            EI_km=300.0,
            dynamic_effectors=(InverseSquaredJ2GravityModel(),),
            mission_time=60.0
        )
        @test_logs (:warn, r"NO aerodynamic force will ever be applied") match_mode=:any run_simulation(cfg_no_aero)
    end
end
