using Test
using LinearAlgebra
using StaticArrays

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
using .SimulationModel

include(joinpath(REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
const run_simulation = SimulationEngine.run_simulation

const EARTH = make_no_gram_planet(:earth)

# The fM wrench internals live in the nested effector module, not on
# SimulationModel itself (the top-level shim only re-exports the coefficient fn).
const AeroFM = SimulationModel.DynamicEffectors.AerodynamicEffectors

# ── Config builders ──────────────────────────────────────────────────────────
# Single-prism spacecraft with a configurable root attitude and Hart-relevant
# dims, so the fixed-attitude incidence modes are separable analytically.

function make_prism_spacecraft(alt_km::Float64; dims=(0.1, 0.1, 0.3), q_root=MVector{4, Float64}(0, 0, 0, 1))
    root = Link{0}(
        root=true, m=4.0,
        dims=MVector{3, Float64}(dims[1], dims[2], dims[3]),
        ref_area=dims[2] * dims[3],
        q=q_root,
    )
    ic = InitialCondition(
        ra=EARTH.Rp_e + alt_km * 1e3,
        rp=EARTH.Rp_e + alt_km * 1e3,
        i=30.0, ω=0.0, Ω=0.0, ν=0.0
    )
    return SpacecraftModel(Joint[], [root], root, true, 4.0, 0.0, root.inertia, 0, 0, ic, 1)
end

function make_incidence_probe_config(;
    incidence::Symbol=:max_drag,
    alt_km::Float64=250.0,
    EI_km::Float64=300.0,
    q_root=MVector{4, Float64}(0, 0, 0, 1),
    mission_time::Float64=900.0,
    aero::Bool=true,
)
    effectors = aero ?
        (InverseSquaredJ2GravityModel(), AerodynamicCoefficientfM(fixed_attitude_incidence=incidence)) :
        (InverseSquaredJ2GravityModel(),)
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
            density_model=ConstantDensityModel(density_kg_m3=1.0e-9, temperature_k=900.0),
            ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=false,
            wind=false
        ),
        dynamics_model=DynamicsModel(
            [make_prism_spacecraft(alt_km; q_root=q_root)],
            effectors,
        ),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(year=2019, month=6, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0,
            reltol_atmosphere=1e-9, abstol_atmosphere=1e-9, dt_max_atmosphere=0.2
        ),
    )
end

function final_state(cfg)
    sol = run_simulation(cfg; return_solution=true)
    @test string(sol.retcode) == "Success"
    sc = sol.u[end].sc[1]
    pos = SVector{3, Float64}(sc.pos[1], sc.pos[2], sc.pos[3])
    vel = SVector{3, Float64}(sc.vel[1], sc.vel[2], sc.vel[3])
    return pos, vel
end

specific_energy(pos, vel) = 0.5 * dot(vel, vel) - EARTH.μ / norm(pos)

# Energy dissipated relative to a drag-free J2-only twin at the SAME final time
# is drag-proportional, so dissipation ratios between incidence modes isolate
# CD*A ratios. The twin must share the mission time: the two-body osculating
# energy carries a J2 oscillation much larger than the drag signal, and only a
# phase-locked baseline cancels it (an early-time baseline aliases it instead).
function dissipation(; incidence::Symbol, q_root=MVector{4, Float64}(0, 0, 0, 1))
    pos, vel = final_state(make_incidence_probe_config(incidence=incidence, q_root=q_root))
    return specific_energy(pos, vel)
end

function j2_baseline_energy(; q_root=MVector{4, Float64}(0, 0, 0, 1))
    pos, vel = final_state(make_incidence_probe_config(aero=false, q_root=q_root))
    return specific_energy(pos, vel)
end

# Direct Hart coefficient evaluation at a controlled incidence for expected
# ratios (mirrors the wrench's 3-arg call reading body.α).
function hart_cd(alpha::Float64; dims=(0.1, 0.1, 0.3))
    b = Link{0}(root=true, m=4.0,
        dims=MVector{3, Float64}(dims[1], dims[2], dims[3]),
        ref_area=dims[2] * dims[3])
    b.α = alpha
    b.β = 0.0
    v = 7500.0
    T = 900.0
    S = sqrt(EARTH.γ * 0.5) * v / sqrt(EARTH.γ * EARTH.R * T)
    _, CD, _, _, _, _ = AeroFM.aerodynamic_coefficient_fM(b, T, S)
    return CD
end

@testset "fM fixed-attitude incidence probes" begin
    @testset "helpers and defaults" begin
        @test AerodynamicCoefficientfM().fixed_attitude_incidence === :max_drag
        @test_throws ArgumentError AeroFM._validate_fm_incidence(:junk)
        # Cauchy mean projected area: box and thin panel
        box = Link{0}(root=true, m=1.0, dims=MVector{3, Float64}(0.1, 0.1, 0.3), ref_area=0.03)
        @test AeroFM._aero_link_area(box, :tumbling_average) ≈ 0.5 * (0.01 + 0.03 + 0.03)
        @test AeroFM._aero_link_area(box, :max_drag) == 0.03
        panel = Link{0}(m=0.2, dims=MVector{3, Float64}(0.0, 0.1, 0.3), ref_area=0.03)
        @test AeroFM._aero_link_area(panel, :tumbling_average) ≈ 0.015  # one-sided area / 2, zero-thickness limit
        # :attitude composes child incidence with the root attitude (child
        # quaternions are root-relative): an identity-quaternion child rigidly
        # mounted on a 45-degree-pitched root shares the root's incidence.
        q45 = MVector{4, Float64}(0.0, sin(pi / 8), 0.0, cos(pi / 8))
        root45 = Link{0}(root=true, m=4.0, dims=MVector{3, Float64}(0.1, 0.1, 0.3), ref_area=0.03, q=q45)
        child_id = Link{0}(m=1.0, dims=MVector{3, Float64}(0.1, 0.1, 0.3), ref_area=0.03)
        α_root = AeroFM._attitude_link_alpha(root45, root45)
        @test abs(α_root - pi / 2) ≈ pi / 4 atol = 1e-12   # 45-degree incidence off flow-normal
        @test AeroFM._attitude_link_alpha(child_id, root45) ≈ α_root atol = 1e-12
        # ...whereas the historical uncomposed child formula reads flow-normal.
        @test AeroFM._quaternion_link_alpha(child_id) ≈ pi / 2 atol = 1e-12
    end

    @testset "attitude mode coincides with max_drag for identity quaternions" begin
        e_max = dissipation(incidence=:max_drag)
        e_att = dissipation(incidence=:attitude)
        @test isapprox(e_max, e_att; rtol=1e-10)
    end

    @testset "a configured 45-degree attitude reduces drag per Hart incidence" begin
        q45 = MVector{4, Float64}(0.0, sin(pi / 8), 0.0, cos(pi / 8))  # 45 deg about body y
        e_j2 = j2_baseline_energy(q_root=q45)                  # drag-free, phase-locked baseline
        e_max = dissipation(incidence=:max_drag, q_root=q45)   # root pinned flow-normal regardless
        e_att = dissipation(incidence=:attitude, q_root=q45)   # honors the 45-degree attitude
        @test e_att > e_max  # less dissipation (higher remaining energy) at oblique incidence
        # quantitative: dissipation ratio tracks the Hart CD ratio (same area)
        ratio_measured = (e_j2 - e_att) / (e_j2 - e_max)
        root45 = Link{0}(root=true, m=4.0, dims=MVector{3, Float64}(0.1, 0.1, 0.3), ref_area=0.03, q=q45)
        α45 = AeroFM._attitude_link_alpha(root45, root45)
        ratio_expected = hart_cd(α45) / hart_cd(pi / 2)
        @test isapprox(ratio_measured, ratio_expected; rtol=0.08)
    end

    @testset "tumbling_average charges the Cauchy mean projected area" begin
        e_j2 = j2_baseline_energy()
        e_max = dissipation(incidence=:max_drag)
        e_tum = dissipation(incidence=:tumbling_average)
        ratio_measured = (e_j2 - e_tum) / (e_j2 - e_max)
        # identical normal-incidence CD; ratio is pure area: (S_total/4)/ref_area
        ratio_expected = (0.5 * (0.01 + 0.03 + 0.03)) / 0.03
        @test isapprox(ratio_measured, ratio_expected; rtol=0.02)
    end
end
