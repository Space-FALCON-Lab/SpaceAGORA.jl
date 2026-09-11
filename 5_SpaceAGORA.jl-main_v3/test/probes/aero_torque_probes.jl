using Test
using LinearAlgebra
using StaticArrays

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
using .SimulationModel

include(joinpath(REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
const run_simulation = SimulationEngine.run_simulation
const ODEParams = SimulationModel.ODEParams

const EARTH = make_no_gram_planet(:earth)
const PE = SimulationModel.DynamicEffectors.PerturbationEffectors
const AE = SimulationModel.DynamicEffectors.AerodynamicEffectors
const ES = parentmodule(PE.StateSample)

# ── shared sample fixtures ───────────────────────────────────────────────────
const EARTH_T = PE.Earth()
const L_PI_I = SMatrix{3, 3, Float64, 9}(I)
const R_SAMPLE = 6698e3 * SVector(1.0, 0.0, 0.0)                 # ~320 km over the equator
const V_PP = SVector(0.0, 7.7e3, 0.0)                            # circular-ish, flow along -y in the airframe
const ATM = ES.AtmosphereSample(1.0e-11, 800.0, SVector(0.0, 0.0, 0.0))

function sample_env()
    alt, lat, lon = AE.rtolatlong(R_SAMPLE, EARTH_T)
    pf = ES.PlanetFrameSample(L_PI_I, R_SAMPLE, V_PP, alt, lat, lon)
    return PE.EnvironmentSample(EARTH_T; planet_frame=pf, atmosphere=ATM)
end

function box_link(; cop=(0.0, 0.0, 0.0), r=(0.0, 0.0, 0.0), q=(0.0, 0.0, 0.0, 1.0), root=true)
    return Link{0}(root=root, m=4.0,
        dims=MVector{3, Float64}(0.1, 0.1, 0.3), ref_area=0.03,
        r=MVector{3, Float64}(r...),
        q=MVector{4, Float64}(q...),
        cop_offset_b=MVector{3, Float64}(cop...))
end

function box_spacecraft(links; q0=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
                        w0=SVector{3, Float64}(0.0, 0.0, 0.0))
    ic = InitialCondition(6698e3, 0.0, 45.0, 0.0, 0.0, 0.0, q0, w0)
    root = links[1]
    return SpacecraftModel(Joint[], collect(links), root, true, 4.0, 0.0, root.inertia, 0, 0, ic, 1)
end

state_for(sc; q_ib=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
          ω=SVector{3, Float64}(0.0, 0.0, 0.0)) =
    PE.StateSample(R_SAMPLE, V_PP, 4.0; q_ib=q_ib, ω_body=ω, spacecraft=sc)

@testset "aero attitude-tracking and CoP torque probes" begin
    env = sample_env()
    fm = AerodynamicCoefficientfM()
    cc = AerodynamicCoefficientConstant()

    @testset "aero geometry tracks the PROPAGATED attitude" begin
        sc = box_spacecraft((box_link(),))
        q_id = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
        q_rot = SVector{4, Float64}(sind(45.0), 0.0, 0.0, cosd(45.0))   # 90 deg about x
        f_id, _ = AE.wrench(fm, state_for(sc; q_ib=q_id), env, 0.0)
        f_rot, _ = AE.wrench(fm, state_for(sc; q_ib=q_rot), env, 0.0)
        @test norm(f_id) > 0.0 && norm(f_rot) > 0.0
        # Rotating the propagated attitude MUST change the fM incidence and
        # therefore the force. The pre-fix code read the static Link.q and
        # returned identical wrenches for any q_ib.
        @test !isapprox(norm(f_id), norm(f_rot); rtol=1e-3)
        # Static Link.q untouched throughout:
        @test SVector{4, Float64}(sc.links[1].q...) == q_id
    end

    @testset "root-link CoP lever arm" begin
        d = 0.05
        cop = SVector(0.0, 0.0, d)
        sc_cop = box_spacecraft((box_link(cop=(cop...,)),))
        sc_zero = box_spacecraft((box_link(),))
        for model in (cc, fm)
            f0, τ0 = AE.wrench(model, state_for(sc_zero), env, 0.0)
            f1, τ1 = AE.wrench(model, state_for(sc_cop), env, 0.0)
            @test τ0 == SVector(0.0, 0.0, 0.0)          # no offset -> exactly no torque
            @test f1 == f0                              # CoP never changes the force
            # identity attitude: root frame == inertial frame here (l_pi = I)
            @test τ1 ≈ cross(cop, f1)
        end
        # ROTATED attitude: the torque frame pin that kills a rot/rot'
        # transpose mutant — force must be rotated INTO the root body frame
        # (rot(q_ib), inertial -> body) before crossing with the body-frame
        # lever. Rotate about z (perpendicular to BOTH the flow and the z
        # lever): the drag force lands along body +/-x depending on the
        # rotation direction, so the two frame conventions give opposite-sign
        # torques. (A rotation about x is degenerate here — it maps the drag
        # axis onto the lever axis and both conventions give zero.)
        q_rot = SVector{4, Float64}(0.0, 0.0, sind(45.0), cosd(45.0))
        for model in (cc, fm)
            f, τ = AE.wrench(model, state_for(sc_cop; q_ib=q_rot), env, 0.0)
            @test norm(f) > 0.0
            @test τ ≈ cross(cop, PE.rot(q_rot) * f) atol = 1e-16
            @test !(τ ≈ cross(cop, PE.rot(q_rot)' * f))  # the transpose is wrong physics
        end
    end

    @testset "panel lever arms: frame-consistent, mirrored pair cancels" begin
        panel = (; y) -> box_link(root=false, r=(0.0, y, 0.0))
        sc_pair = box_spacecraft((box_link(), panel(y=0.5), panel(y=-0.5)))
        _, τ_pair = AE.wrench(cc, state_for(sc_pair), env, 0.0)
        @test norm(τ_pair) < 1e-18                      # symmetric pair -> exact cancel

        # single canted panel with a link-frame CoP: lever arm must compose
        # through the child attitude (rot(q_child)') into the root frame.
        q_child = SVector{4, Float64}(sind(45.0), 0.0, 0.0, cosd(45.0))  # 90 deg about x
        cop_link = SVector(0.0, 0.0, 0.04)
        one_panel = box_link(root=false, r=(0.0, 0.5, 0.0),
            q=(q_child...,), cop=(cop_link...,))
        sc_one = box_spacecraft((box_link(), one_panel))
        f_one, τ_one = AE.wrench(cc, state_for(sc_one), env, 0.0)
        sc_root_only = box_spacecraft((box_link(),))
        f_root, τ_root = AE.wrench(cc, state_for(sc_root_only), env, 0.0)
        f_panel = f_one - f_root                        # the panel's own force share
        lever = SVector(0.0, 0.5, 0.0) + PE.rot(q_child)' * cop_link
        @test τ_one ≈ cross(lever, f_panel) atol = 1e-12
    end

    @testset "constant-model drag never thrusts (folded incidence)" begin
        # The constant CD(alpha) line is physical on [0, pi/2] only; the signed
        # atan2 incidence a propagated attitude produces must be folded, or
        # alpha = -pi/2 yields CD = -0.6 and drag pumps orbital energy
        # (Codex, PR #86). Sweep attitudes: drag must always oppose the flow,
        # and reversed flow must see the same face (fold symmetry).
        sc = box_spacecraft((box_link(),))
        v_hat = normalize(V_PP)
        for ax in (SVector(1.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0),
                   normalize(SVector(1.0, 1.0, 1.0)))
            for deg in (-150.0, -90.0, -45.0, 30.0, 90.0, 135.0, 180.0)
                q = SVector{4, Float64}((sind(deg / 2) .* ax)..., cosd(deg / 2))
                f, _ = AE.wrench(cc, state_for(sc; q_ib=q), env, 0.0)
                @test dot(f, v_hat) < 0.0                # always drag, never thrust
            end
        end
        # fold symmetry: flow along body -x sees the max-drag face exactly
        # like flow along body +x (90 deg vs -90 deg about z).
        q_p = SVector{4, Float64}(0.0, 0.0, sind(45.0), cosd(45.0))
        q_m = SVector{4, Float64}(0.0, 0.0, -sind(45.0), cosd(45.0))
        f_p, _ = AE.wrench(cc, state_for(sc; q_ib=q_p), env, 0.0)
        f_m, _ = AE.wrench(cc, state_for(sc; q_ib=q_m), env, 0.0)
        @test norm(f_p) ≈ norm(f_m)
        @test AE._fold_constant_incidence(-pi / 2) == pi / 2
        @test AE._fold_constant_incidence(Float64(pi)) == 0.0
        @test AE._fold_constant_incidence(pi / 2) == pi / 2   # historical range untouched
    end

    @testset "child composition order under a rotated root (fM)" begin
        # A child at attitude q_child on a root at q_ib must see exactly the
        # incidence of a root-only body at the COMPOSED attitude. Different
        # rotation axes make the composition order observable — the
        # wrong-order mutant rot(q_child)' * rot(q_ib)' fails this.
        qmul(a, b) = SVector{4, Float64}(
            a[4] * b[1] + b[4] * a[1] + (a[2] * b[3] - a[3] * b[2]),
            a[4] * b[2] + b[4] * a[2] + (a[3] * b[1] - a[1] * b[3]),
            a[4] * b[3] + b[4] * a[3] + (a[1] * b[2] - a[2] * b[1]),
            a[4] * b[4] - (a[1] * b[1] + a[2] * b[2] + a[3] * b[3]))
        q_ib = SVector{4, Float64}(sind(45.0), 0.0, 0.0, cosd(45.0))       # 90 deg about x
        q_child = SVector{4, Float64}(0.0, sind(45.0), 0.0, cosd(45.0))    # 90 deg about y
        R_expect = PE.rot(q_child) * PE.rot(q_ib)
        cand1, cand2 = qmul(q_child, q_ib), qmul(q_ib, q_child)
        q_comb = isapprox(PE.rot(cand1), R_expect; atol=1e-12) ? cand1 : cand2
        @test PE.rot(q_comb) ≈ R_expect atol = 1e-12                        # self-validating helper

        child = box_link(root=false, r=(0.0, 0.7, 0.0), q=(q_child...,))
        f_both, _ = AE.wrench(fm, state_for(box_spacecraft((box_link(), child)); q_ib=q_ib), env, 0.0)
        f_root, _ = AE.wrench(fm, state_for(box_spacecraft((box_link(),)); q_ib=q_ib), env, 0.0)
        f_child = f_both - f_root
        f_equiv, _ = AE.wrench(fm, state_for(box_spacecraft((box_link(),)); q_ib=q_comb), env, 0.0)
        @test norm(f_child) > 0.0
        @test f_child ≈ f_equiv atol = 1e-12
    end

    @testset "per-link atmosphere samples the PROPAGATED link position" begin
        q_ib = SVector{4, Float64}(sind(45.0), 0.0, 0.0, cosd(45.0))       # 90 deg about x
        r_child = SVector(0.0, 0.7, 0.0)
        child = box_link(root=false, r=(r_child...,))
        sc = box_spacecraft((box_link(), child))
        sampled = SVector{3, Float64}[]
        fn = pos -> (push!(sampled, pos); (ATM.rho_kg_m3, ATM.temperature_k, nothing))
        AE._aero_pure_wrench(:fm, state_for(sc; q_ib=q_ib), env, fn)
        @test length(sampled) == 1                                          # non-root links only
        # body.r is a root-frame offset: it must rotate with the propagated
        # root attitude (l_pi = I here, so planet frame == inertial frame).
        @test sampled[1] ≈ R_SAMPLE + PE.rot(q_ib)' * r_child atol = 1e-9
        @test !(sampled[1] ≈ R_SAMPLE + r_child)                            # unrotated is wrong
    end

    @testset "missing propagated attitude in the orientation branch throws" begin
        sc = box_spacecraft((box_link(),))
        @test_throws ArgumentError AE._aero_link_angles(
            sc, sc.links[1], 1, true, SVector(0.0, 1.0, 0.0), 0.0, :max_drag, nothing)
    end

    @testset "orientation_sim=false stays force-only and CoP-blind" begin
        sc_cop = box_spacecraft((box_link(cop=(0.0, 0.0, 0.05)),))
        sc_zero = box_spacecraft((box_link(),))
        for model in (cc, fm)
            f1, τ1 = AE.wrench(model, PE.StateSample(R_SAMPLE, V_PP, 4.0; spacecraft=sc_cop), env, 0.0)
            f0, τ0 = AE.wrench(model, PE.StateSample(R_SAMPLE, V_PP, 4.0; spacecraft=sc_zero), env, 0.0)
            @test τ1 == SVector(0.0, 0.0, 0.0) && τ0 == SVector(0.0, 0.0, 0.0)
            @test f1 == f0                              # cop invisible without attitude
        end
    end

    @testset "end-to-end: CoP torque drives the attitude, not the orbit" begin
        atmo = ExponentialAtmosphereModel(1.0e-11, 300e3, 50e3;
            temperature_k=800.0, valid_max_altitude_m=1000e3)
        function config(sc)
            return SimulationConfiguration(
                simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
                mission_configuration=MissionConfiguration(
                    mission_type=MissionTime, keplerian=false, number_of_orbits=1,
                    mission_time=600.0, orientation_sim=true, num_steps_to_save=500),
                environment_model=EnvironmentModel(
                    planet=EARTH, EI=300.0,
                    density_model=atmo,
                    ephemerides_model=SimpleEphemeridesModel(),
                    thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
                    topography=false, wind=false),
                dynamics_model=DynamicsModel([sc], (InverseSquaredGravityModel(), AerodynamicCoefficientfM())),
                guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
                navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
                control_model=ControlModel((), Float64[]),
                initial_time=InitialTime(year=2015, month=6, day=1, hour=0, minute=0, second=0.0),
                integration_tolerances=IntegrationTolerances(
                    reltol_orbit=1e-8, abstol_orbit=1e-8, dt_max_orbit=5.0,
                    reltol_atmosphere=1e-8, abstol_atmosphere=1e-8, dt_max_atmosphere=1.0),
            )
        end
        mk(cop) = box_spacecraft((box_link(cop=cop),);
            w0=SVector{3, Float64}(0.0, 0.0, 0.01))
        sol_cop = run_simulation(config(mk((0.0, 0.0, 0.05))); return_solution=true)
        sol_zero = run_simulation(config(mk((0.0, 0.0, 0.0))); return_solution=true)
        @test string(sol_cop.retcode) == "Success" && string(sol_zero.retcode) == "Success"
        ω_cop = SVector{3, Float64}(sol_cop.u[end].sc[1].ω...)
        ω_zero = SVector{3, Float64}(sol_zero.u[end].sc[1].ω...)
        @test !(ω_cop ≈ ω_zero)                          # the torque acted on the attitude
        p_cop = SVector{3, Float64}(sol_cop.u[end].sc[1].pos...)
        p_zero = SVector{3, Float64}(sol_zero.u[end].sc[1].pos...)
        @test norm(p_cop - p_zero) < 1.0                 # ...but barely on the orbit (< 1 m over 600 s)
    end
end
