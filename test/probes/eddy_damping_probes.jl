using Test
using LinearAlgebra
using Logging
using StaticArrays
using Random

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

using SpaceAGORA
const SimulationModel = SpaceAGORA.SimulationModel
using .SimulationModel

const SimulationEngine = SpaceAGORA.SimulationEngine
const run_simulation = SimulationEngine.run_simulation
const build_initial_conditions = SimulationEngine.build_initial_conditions
const ODEParams = SimulationModel.ODEParams

const EARTH = make_no_gram_planet(:earth)
const PE = SimulationModel.DynamicEffectors.PerturbationEffectors

# ── builders ─────────────────────────────────────────────────────────────────
function make_magnet_spacecraft(; q0, w0, m_dipole=1.0)
    root = Link{0}(root=true, m=4.0,
        dims=MVector{3, Float64}(0.1, 0.1, 0.3), ref_area=0.03)
    push!(root.magnets, Magnet(m=MVector{3, Float64}(0.0, 0.0, m_dipole)))
    ra = EARTH.Rp_e + 520e3
    ic = InitialCondition(ra, 0.0, 60.0, 0.0, 0.0, 0.0, q0, w0)
    return SpacecraftModel(Joint[], [root], root, true, 4.0, 0.0, root.inertia, 0, 0, ic, 1)
end

function make_attitude_config(; effectors, mission_time=600.0, orientation_sim=true,
                              q0=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
                              w0=SVector{3, Float64}(0.02, 0.01, 0.0628))
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime, keplerian=false, number_of_orbits=1,
            mission_time=mission_time, orientation_sim=orientation_sim,
            num_steps_to_save=500),
        environment_model=EnvironmentModel(
            planet=EARTH, EI=300.0,
            density_model=NoAtmosphereModel(),
            ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=false, wind=false),
        dynamics_model=DynamicsModel([make_magnet_spacecraft(q0=q0, w0=w0)], effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(year=2015, month=6, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-8, abstol_orbit=1e-8, dt_max_orbit=5.0,
            reltol_atmosphere=1e-8, abstol_atmosphere=1e-8, dt_max_atmosphere=1.0),
    )
end

rot_ke(sc) = begin
    ω = SVector{3, Float64}(sc.ω[1], sc.ω[2], sc.ω[3])
    0.5 * dot(ω, ω)   # inertia-weighted constant factor irrelevant for trend
end

@testset "eddy-current damping probes" begin
    @testset "torque law and validation" begin
        rng = MersenneTwister(7)
        for _ in 1:50
            B = SVector{3, Float64}(randn(rng, 3)...) .* 3e-5
            ω = SVector{3, Float64}(randn(rng, 3)...) .* 0.05
            τ = eddy_damping_torque(1.0e3, B, ω)
            @test τ ≈ 1.0e3 .* cross(B, cross(B, ω))
            # independent BAC-CAB formulation — not a transcription of the code
            @test τ ≈ 1.0e3 .* (B .* dot(B, ω) .- ω .* dot(B, B))
            @test dot(τ, ω) <= 1e-18   # dissipative
        end
        Bz = SVector{3, Float64}(2e-5, -1e-5, 3e-5)
        @test eddy_damping_torque(1.0e3, Bz, SVector{3, Float64}(0.0, 0.0, 0.0)) ==
              SVector{3, Float64}(0.0, 0.0, 0.0)
        B = SVector{3, Float64}(2e-5, -1e-5, 3e-5)
        ω_par = 0.05 .* B ./ norm(B)
        @test norm(eddy_damping_torque(1.0e3, B, ω_par)) ≈ 0.0 atol = 1e-18
        @test_throws ArgumentError EddyCurrentDampingModel(k_e=-1.0)
        @test_throws ArgumentError EddyCurrentDampingModel(k_e=0.0)
        @test_throws ArgumentError EddyCurrentDampingModel(k_e=Inf)
        @test_throws ArgumentError EddyCurrentDampingModel(k_e=NaN)
        @test_throws ArgumentError EddyCurrentDampingModel(k_e=1.0, field_model=:igrf)
        @test_throws ArgumentError EddyCurrentDampingModel(k_e=1.0, field_model=:sideways)
        m = EddyCurrentDampingModel(k_e=5.0, field_model=:igrf, igrf_year=2015.4)
        @test m.k_e == 5.0 && m.igrf_year == 2015.4
    end

    @testset "post-2030 IGRF epoch: one construction warning, silent hot loop" begin
        # The library's reduced-accuracy warning has no maxlog and would fire
        # on every RHS call; the model warns exactly once at construction and
        # _magnetic_field_inertial evaluates warning-free.
        m32 = @test_logs (:warn, r"reduced for epochs past 2030") EddyCurrentDampingModel(
            k_e=1.0, field_model=:igrf, igrf_year=2032.0)
        @test m32.igrf_year == 2032.0
        earth = PE.Earth()
        l_pi = SMatrix{3, 3, Float64, 9}(I)
        r = 6898e3 * SVector(cosd(25.0) * cosd(40.0), cosd(25.0) * sind(40.0), sind(25.0))
        alt, lat, lon = PE.rtolatlong(r, earth)
        B_ii = @test_logs min_level = Logging.Warn PE._magnetic_field_inertial(
            m32, l_pi, r, lat, lon, alt)
        @test 1.5e-5 < norm(B_ii) < 7e-5   # the epoch actually evaluated
        # Pre-2031 epochs stay warning-free at construction too.
        @test_logs min_level = Logging.Warn EddyCurrentDampingModel(
            k_e=1.0, field_model=:igrf, igrf_year=2025.4)
    end

    @testset "ODE-path torque: dissipative when attitude simulated, zero otherwise" begin
        model = EddyCurrentDampingModel(k_e=1.0e4)
        args_on = make_attitude_config(
            effectors=(InverseSquaredGravityModel(), model))
        p = ODEParams(n_sats=1, args=args_on)
        u = build_initial_conditions(args_on)
        p.shared_buffers.current_time[] = 0.0
        f, τ = SimulationModel.calcForceTorque(model, u.sc[1], p, 1)
        @test f == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test all(isfinite, τ) && norm(τ) > 0.0
        ω0 = SVector{3, Float64}(Float64.(getproperty(u.sc[1], :ω))...)
        @test dot(τ, ω0) < 0.0

        args_off = make_attitude_config(
            effectors=(InverseSquaredGravityModel(), model), orientation_sim=false)
        p2 = ODEParams(n_sats=1, args=args_off)
        u2 = build_initial_conditions(args_off)
        p2.shared_buffers.current_time[] = 0.0
        f2, τ2 = SimulationModel.calcForceTorque(model, u2.sc[1], p2, 1)
        @test f2 == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test τ2 == SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    @testset "wrench path: the live engine route, both field sources, guards" begin
        # Real simulations reach this effector through wrench(model, StateSample,
        # EnvironmentSample, t) — _wrench_method_available short-circuits
        # calcForceTorque — so the translational no-op guard that matters is the
        # one in wrench, pinned here directly and end-to-end below.
        ES = parentmodule(PE.StateSample)
        earth = PE.Earth()
        l_pi = SMatrix{3, 3, Float64, 9}(I)
        r = 6898e3 * SVector(cosd(25.0) * cosd(40.0), cosd(25.0) * sind(40.0), sind(25.0))
        alt, lat, lon = PE.rtolatlong(r, earth)
        pf = ES.PlanetFrameSample(l_pi, r, SVector(0.0, 0.0, 0.0), alt, lat, lon)
        env = PE.EnvironmentSample(earth; planet_frame=pf)
        q = SVector{4, Float64}(sind(15.0), 0.0, 0.0, cosd(15.0))
        ω = SVector{3, Float64}(0.02, -0.01, 0.03)
        x = PE.StateSample(r, SVector(7.6e3, 0.0, 0.0), 4.0; q_ib=q, ω_body=ω)
        k_e = 2.0e3
        for model in (EddyCurrentDampingModel(k_e=k_e),
                      EddyCurrentDampingModel(k_e=k_e, field_model=:igrf, igrf_year=2015.12))
            B_ii = PE._magnetic_field_inertial(model, l_pi, r, lat, lon, alt)
            @test 1.5e-5 < norm(B_ii) < 7e-5          # a real LEO-band field was evaluated
            B_b = PE.rot(q) * B_ii
            f, τ = PE.wrench(model, x, env, 0.0)
            @test f == SVector(0.0, 0.0, 0.0)
            @test τ ≈ k_e .* (B_b .* dot(B_b, ω) .- ω .* dot(B_b, B_b))
            @test dot(τ, ω) < 0.0
            # translational guard: any missing attitude channel -> exact zero
            for (qg, ωg) in ((nothing, ω), (q, nothing), (nothing, nothing))
                xg = PE.StateSample(r, SVector(7.6e3, 0.0, 0.0), 4.0; q_ib=qg, ω_body=ωg)
                @test PE.wrench(model, xg, env, 0.0) ==
                      (SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0))
            end
            @test_throws ArgumentError PE.wrench(model, x, PE.EnvironmentSample(earth), 0.0)
        end

        # ODE-path :igrf branch (the rtolatlong + _magnetic_field_inertial copy)
        model_igrf = EddyCurrentDampingModel(k_e=1.0e4, field_model=:igrf, igrf_year=2015.12)
        args = make_attitude_config(effectors=(InverseSquaredGravityModel(), model_igrf))
        p = ODEParams(n_sats=1, args=args)
        u = build_initial_conditions(args)
        p.shared_buffers.current_time[] = 0.0
        f, τ = SimulationModel.calcForceTorque(model_igrf, u.sc[1], p, 1)
        @test f == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test all(isfinite, τ) && norm(τ) > 0.0
        ω0 = SVector{3, Float64}(Float64.(getproperty(u.sc[1], :ω))...)
        @test dot(τ, ω0) < 0.0

        # end-to-end translational no-op: with orientation_sim=false the engine's
        # wrench route must leave the trajectory bit-identical.
        base = make_attitude_config(
            effectors=(InverseSquaredGravityModel(),), orientation_sim=false)
        plus = make_attitude_config(
            effectors=(InverseSquaredGravityModel(), EddyCurrentDampingModel(k_e=1.0e4)),
            orientation_sim=false)
        sol_b = run_simulation(base; return_solution=true)
        sol_p = run_simulation(plus; return_solution=true)
        @test string(sol_b.retcode) == "Success" && string(sol_p.retcode) == "Success"
        @test SVector{3, Float64}(sol_p.u[end].sc[1].pos...) ==
              SVector{3, Float64}(sol_b.u[end].sc[1].pos...)
        @test SVector{3, Float64}(sol_p.u[end].sc[1].vel...) ==
              SVector{3, Float64}(sol_b.u[end].sc[1].vel...)
    end

    @testset "end-to-end spin-down: damping dissipates, magnet alone conserves" begin
        damped = make_attitude_config(effectors=(
            InverseSquaredGravityModel(),
            MagneticTorqueRodModel(),
            EddyCurrentDampingModel(k_e=1.0e4)))
        undamped = make_attitude_config(effectors=(
            InverseSquaredGravityModel(),
            MagneticTorqueRodModel()))
        sol_d = run_simulation(damped; return_solution=true)
        sol_u = run_simulation(undamped; return_solution=true)
        @test string(sol_d.retcode) == "Success"
        @test string(sol_u.retcode) == "Success"
        ke0_d = rot_ke(sol_d.u[1].sc[1]);  ke1_d = rot_ke(sol_d.u[end].sc[1])
        ke0_u = rot_ke(sol_u.u[1].sc[1]);  ke1_u = rot_ke(sol_u.u[end].sc[1])
        @test ke1_d < 0.9 * ke0_d                       # damping visibly dissipates
        @test abs(ke1_u - ke0_u) < 0.05 * ke0_u         # magnet torque alone ~conserves |ω|²
        @test ke1_d < ke1_u                             # damped ends below undamped
    end
end
