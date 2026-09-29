module SplitIMEXAccuracyTests
using Test
using SpaceAGORA
using OrdinaryDiffEq
using SparseArrays
const SE = SpaceAGORA.SimulationEngine
const SM = SpaceAGORA.SimulationModel

# Nonlinear stiff relaxation with the manufactured exact solution
# u(t) = [cos(t), 2cos(t)]. The explicit forcing is nonzero, so this
# exercises the split path and Newton solve rather than an unsplit surrogate.
function implicit!(du, u, stiffness, t)
    c = cos(t)
    du[1] = -stiffness * (u[1]^3 - c^3)
    du[2] = -stiffness * (u[2]^3 - 8c^3)
    return nothing
end
function explicit!(du, u, p, t)
    du[1] = -sin(t)
    du[2] = -2sin(t)
    return nothing
end
exact(t) = [cos(t), 2cos(t)]

@testset "Nonlinear split accuracy with component tolerances and reuse" begin
    rtol, atol = [1e-9, 1e-10], [1e-10, 1e-11]
    @test SM.SolverConfig().split_imex_solver === :kencarp4
    for sparse_jac in (false, true), mode in (:kencarp4,)
        f = SE._split_component_function(implicit!, sparse_jac ? spdiagm(ones(2)) : nothing)
        prob = SplitODEProblem(f, explicit!, exact(0.0), (0.0, 0.8), 100.0)
        alg = SE._split_imex_solver_spec(SM.SolverConfig(split_imex_solver=mode), sparse_jac).alg
        integrator = init(prob, alg; reltol=rtol, abstol=atol, dtmax=0.01)
        sol = solve!(integrator)
        @test string(sol.retcode) == "Success"
        @test sol.u[end] ≈ exact(0.8) rtol=1e-7 atol=1e-10
        @test integrator.opts.reltol == rtol
        @test integrator.opts.abstol == atol

        integrator.p = 75.0
        reinit!(integrator, exact(0.2); t0=0.2, tf=1.0)
        reused = solve!(integrator)
        @test string(reused.retcode) == "Success"
        @test reused.u[end] ≈ exact(1.0) rtol=1e-7 atol=1e-10
        @test integrator.opts.reltol == rtol
        @test integrator.opts.abstol == atol
    end
end

@testset "Alternative split algorithms retain their nonlinear policy" begin
    # Applying fresh Jacobians to every KenCarp method regressed KenCarp47
    # on this shifted-time exact-solution case. Keep its established behavior.
    for sparse_jac in (false, true), mode in (:kencarp47, :kencarp58)
        f = SE._split_component_function(implicit!, sparse_jac ? spdiagm(ones(2)) : nothing)
        prob = SplitODEProblem(f, explicit!, exact(0.2), (0.2, 1.0), 75.0)
        alg = SE._split_imex_solver_spec(SM.SolverConfig(split_imex_solver=mode), sparse_jac).alg
        sol = solve(prob, alg; reltol=[1e-9, 1e-10], abstol=[1e-10, 1e-11], dtmax=0.01)
        @test string(sol.retcode) == "Success"
        @test sol.u[end] ≈ exact(1.0) rtol=1e-7 atol=1e-10
    end
end

@testset "Default split IMEX preserves bounded atmospheric accuracy" begin
    # This nonlinear public-API case detects the dependency migration's
    # KenCarp4 regression that a linear cancellation probe cannot expose.
    # The endpoint limits below apply to this 600-second fixture only.
    planet = SM.make_no_gram_planet(:earth)
    function configuration(mode, directory)
        body = SM.Link(root=true, m=400.0, ref_area=4.0)
        ic = SM.InitialCondition(ra=planet.Rp_e + 400_000.0, rp=planet.Rp_e + 130_000.0,
            i=40.0, ω=0.0, Ω=0.0, ν=-35.0)
        sc = SM.SpacecraftModel(SM.Joint[], [body], body, true, 400.0, 0.0,
            body.inertia, 0, 0, ic, 1)
        args = SpaceAGORA.TelemetryVerification.make_example_config(
            planet=planet, spacecraft=sc, mission_time=600.0,
            initial_time=SM.InitialTime(year=2021, month=3, day=4, hour=5),
            dynamic_effectors=(SM.InverseSquaredJ2GravityModel(), SM.AerodynamicCoefficientfM()),
            density_model=SM.ExponentialAtmosphereModel(planet),
            ephemerides_model=SM.SimpleEphemeridesModel(), EI_km=300.0,
            results=false, verbose=false, results_directory=directory,
            solver_config=SM.SolverConfig(solver_mode=mode))
        mission = SM.SimConfig.MissionConfiguration(mission_type=SM.MissionTime,
            keplerian=false, number_of_orbits=1, mission_time=600.0,
            orientation_sim=false, num_steps_to_save=1000, data_rate=5.0)
        return SM.SimConfig._with_configuration(args; mission_configuration=mission)
    end
    mktempdir() do directory
        ref = SpaceAGORA.run_simulation(configuration(:dp8, directory); return_solution=true)
        candidate = SpaceAGORA.run_simulation(configuration(:split_imex, directory); return_solution=true)
        @test string(ref.retcode) == string(candidate.retcode) == "Success"
        @test ref.t[end] == candidate.t[end] == 600.0
        delta = collect(candidate.u[end]) .- collect(ref.u[end])
        @test sqrt(sum(abs2, delta[1:3])) < 0.01
        @test sqrt(sum(abs2, delta[4:6])) < 5e-5
    end
end
end
