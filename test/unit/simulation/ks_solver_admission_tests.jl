module KSSolverAdmissionTests

# Admission boundary and optional-package handling of the KS solver route,
# exercised without the KSPropagator extension. The extension is a weak
# dependency that ordinary environments and CI do not install, so these tests
# cover the core-side validation, the runtime checks and the trial RHS service.

using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA
using SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const T0 = SM.InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0)
const PLANET = SM.make_no_gram_planet(:earth)

mklink(; root=false, m, dims, r=(0.0, 0.0, 0.0), q=(0.0, 0.0, 0.0, 1.0)) =
    SM.Link(root=root, m=m, dims=MVector{3, Float64}(dims...), r=MVector{3, Float64}(r...), q=MVector{4, Float64}(q...))
ic_at(r, v) = SM.CartesianInitialCondition(SVector{3, Float64}(r), SVector{3, Float64}(v);
    q=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0), ang_vel=SVector{3, Float64}(0.0, 0.0, 0.0))

# One rigid link, no mounted hardware: the shape the initial KS route admits.
function rigid_bus()
    bus = mklink(root=true, m=20.0, dims=(1.0, 0.8, 0.6))
    return SM.SpacecraftModel(; joints=SM.Joint[], links=[bus], root=bus,
        initial_condition=ic_at((7.0e6, 0.0, 0.0), (0.0, 7.5e3, 0.0)),
        inertia_tensor=bus.inertia, attachments=SM.CompliantAttachment[])
end

function ks_args(; spacecraft=rigid_bus(), effectors=(SM.InverseSquaredJ2GravityModel(),),
        density=SM.NoAtmosphereModel(), orientation=false, ephem=SM.SimpleEphemeridesModel(), mode=:ks)
    return make_example_config(planet=PLANET, spacecraft=spacecraft, mission_time=10.0, initial_time=T0,
        dynamic_effectors=effectors, density_model=density, ephemerides_model=ephem,
        orientation_sim=orientation, keplerian=true, verbose=false, results=false,
        results_directory=mktempdir(), solver_config=SM.SolverConfig(solver_mode=mode))
end

function refusal_message(f)
    err = try
        f()
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    return err isa ArgumentError ? err.msg : ""
end

@testset "KS solver mode is a named, optional solver" begin
    @test SM.SolverConfig(solver_mode=:ks).solver_mode === :ks
    @test_throws ArgumentError SM.SolverConfig(solver_mode=:kss)
    @test SM.SolverConfig(solver_mode=:tsit5).solver_mode === :tsit5
end

@testset "trial RHS safety declarations" begin
    @test SE.trial_rhs_safe(SM.InverseSquaredGravityModel())
    @test SE.trial_rhs_safe(SM.InverseSquaredJ2GravityModel())
    @test SE.trial_rhs_safe(SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=PLANET))
    # Unknown models are conservative by default.
    @test !SE.trial_rhs_safe(nothing)
    @test !SE.trial_rhs_safe("not a model")
    @test !SE.trial_rhs_safe(SM.AerodynamicCommandedAreaDragModel())
end

const EXTENSION_ABSENT = Base.get_extension(SpaceAGORA, :SpaceAGORAKSPropagatorExt) === nothing

@testset "admitted configuration stops at the optional package" begin
    if !EXTENSION_ABSENT
        @test_skip "KS extension is loaded in this environment"
    else
        message = refusal_message(() -> SE._validate_ks_config!(ks_args()))
        @test occursin("KSPropagator", message)
        # The same refusal reaches run_simulation before any integration.
        run_message = refusal_message(() -> run_simulation(ks_args()))
        @test occursin("KSPropagator", run_message)
        # A deterministic exponential atmosphere is inside the boundary.
        atmosphere_message = refusal_message(() -> SE._validate_ks_config!(ks_args(density=SM.ExponentialAtmosphereModel(PLANET))))
        @test occursin("KSPropagator", atmosphere_message)
        # Ordinary solvers are untouched by the KS validation.
        @test SM.SolverConfig(solver_mode=:tsit5) === ks_args(mode=:tsit5).solver_config || ks_args(mode=:tsit5).solver_config.solver_mode === :tsit5
    end
end

@testset "admission boundary refusals name their rule" begin
    @test occursin("orientation_sim=false", refusal_message(() -> SE._validate_ks_config!(ks_args(orientation=true))))
    @test occursin("NoAtmosphereModel or deterministic ExponentialAtmosphereModel",
        refusal_message(() -> SE._validate_ks_config!(ks_args(density=SM.PolynomialFitAtmosphereModel(PLANET)))))
    @test occursin("exactly one configured inverse-square gravity model",
        refusal_message(() -> SE._validate_ks_config!(ks_args(effectors=()))))
    @test occursin("exactly one configured inverse-square gravity model",
        refusal_message(() -> SE._validate_ks_config!(ks_args(effectors=(SM.InverseSquaredJ2GravityModel(), SM.InverseSquaredGravityModel())))))
    @test occursin("has not declared trial_rhs_safe",
        refusal_message(() -> SE._validate_ks_config!(ks_args(effectors=(SM.InverseSquaredJ2GravityModel(), SM.AerodynamicCommandedAreaDragModel())))))
    three_links = make_three_body_spacecraft(bus_dims=(1.2, 1.1, 0.9), panel_dims=(0.01, 0.8, 0.4),
        bus_mass=50.0, panel_mass_each=2.0, panel_offset_y=1.0, prop_mass=0.0,
        ic=SM.InitialCondition(ra=PLANET.Rp_e + 500.0e3, rp=PLANET.Rp_e + 400.0e3, i=10.0, Ω=0.0, ω=0.0, ν=0.0))
    @test occursin("single rigid spacecraft link", refusal_message(() -> SE._validate_ks_config!(ks_args(spacecraft=three_links))))
end

@testset "runtime checks and the trial RHS service" begin
    args = ks_args()
    p = SM.ODEParams(n_sats=1, args=args)
    SE._initialize_runtime_env_config!(p)
    @test SE._validate_ks_runtime!(p) === nothing
    u = SE.build_initial_conditions(args)
    du = zero(u)
    SE._trial_cartesian_rhs!(du, u, p, 0.0)
    @test all(isfinite, du)
    @test any(!iszero, du.sc[1].vel)                       # gravity acts on the trial state
    @test du.sc[1].pos ≈ u.sc[1].vel                        # kinematics of the configured RHS
    # A second trial at the same time with another state must re-evaluate, not reuse a stamp.
    u2 = copy(u); u2.sc[1].pos .*= 1.001
    du2 = zero(u)
    SE._trial_cartesian_rhs!(du2, u2, p, 0.0)
    @test !(du2.sc[1].vel ≈ du.sc[1].vel)
    # Frozen per-step atmosphere samples are refused for a state-dependent atmosphere.
    withenv("SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "on") do
        frozen = ks_args(density=SM.ExponentialAtmosphereModel(PLANET))
        pf = SM.ODEParams(n_sats=1, args=frozen)
        SE._initialize_runtime_env_config!(pf)
        @test occursin("SPACEAGORA_DENSITY_FREEZE_PER_STEP=false", refusal_message(() -> SE._validate_ks_runtime!(pf)))
        uf = SE.build_initial_conditions(frozen)
        @test occursin("density_freeze_per_step=false", refusal_message(() -> SE._trial_cartesian_rhs!(zero(uf), uf, pf, 0.0)))
    end
end

end
