using Test
using CSV
using DataFrames
using LinearAlgebra
using StaticArrays
using SPICE

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(REPO_ROOT, "src", "simulation_model", "SimulationModel.jl"))
using .SimulationModel

# run_simulation.jl expects quat_mult in the including scope.
const quat_mult = SimulationModel.quat_mult
include(joinpath(REPO_ROOT, "src", "simulation", "run_simulation.jl"))

struct ConstantDensityModel <: SimulationModel.AbstractDensityModel
    rho::Float64
    temp::Float64
end

function SimulationModel.EnvironmentModels.getDensity(
    model::ConstantDensityModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p
)
    return model.rho, model.temp, SVector{3, Float64}(0.0, 0.0, 0.0)
end

const SPICE_PATH = joinpath(REPO_ROOT, "GRAM_Data", "SPICE")
const EARTH = Earth("", SPICE_PATH)

angle_distance(a::Float64, b::Float64) = abs(mod(a - b + pi, 2pi) - pi)

function specific_energy(df::DataFrame, mu::Float64)
    r = sqrt.(df.sc1_pos_1.^2 .+ df.sc1_pos_2.^2 .+ df.sc1_pos_3.^2)
    v2 = df.sc1_vel_1.^2 .+ df.sc1_vel_2.^2 .+ df.sc1_vel_3.^2
    return 0.5 .* v2 .- mu ./ r
end

function make_spacecraft(;
    ra_alt_m::Float64,
    rp_alt_m::Float64,
    i_deg::Float64=35.0,
    ω_deg::Float64=40.0,
    Ω_deg::Float64=10.0,
    ν_deg::Float64=175.0,
    orientation_state::Union{Nothing, Tuple{SVector{4, Float64}, SVector{3, Float64}}}=nothing
)
    root = Link{0}(root=true, m=500.0, ref_area=12.0)
    panel = Link{0}(root=false, m=30.0, ref_area=6.0, r=MVector{3, Float64}(0.0, 1.2, 0.0))

    if isnothing(orientation_state)
        ic = InitialCondition(
            ra=EARTH.Rp_e + ra_alt_m,
            rp=EARTH.Rp_e + rp_alt_m,
            i=i_deg,
            ω=ω_deg,
            Ω=Ω_deg,
            ν=ν_deg
        )
    else
        q0, w0 = orientation_state
        ra = EARTH.Rp_e + ra_alt_m
        rp = EARTH.Rp_e + rp_alt_m
        a = (ra + rp) / 2.0
        e = (ra - rp) / (ra + rp)
        ic = InitialCondition(a, e, i_deg, ω_deg, Ω_deg, ν_deg, q0, w0)
    end

    dry_mass = root.m + panel.m
    return SpacecraftModel(
        Joint[],
        [root, panel],
        root,
        true,
        dry_mass,
        0.0,
        root.inertia,
        0,
        0,
        ic,
        1
    )
end

function build_config(;
    spacecraft::SpacecraftModel,
    density_model,
    orientation_sim::Bool,
    mission_time::Float64,
    EI_km::Float64,
    dynamic_effectors::Tuple,
    keplerian::Bool=true,
    tolerances::IntegrationTolerances=IntegrationTolerances()
)
    environment_model = EnvironmentModel(
        planet=EARTH,
        EI=EI_km,
        density_model=density_model,
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
        topography=false,
        wind=false
    )

    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=true, verbose=false, generate_plots=false),
        mission_configuration=MissionConfiguration(
            mission_type="Time",
            keplerian=keplerian,
            number_of_orbits=1,
            mission_time=mission_time,
            orientation_sim=orientation_sim,
            num_steps_to_save=1000
        ),
        environment_model=environment_model,
        dynamics_model=DynamicsModel([spacecraft], dynamic_effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=tolerances
    )
end

function run_case(args::SimulationConfiguration)
    return mktempdir() do tmp
        cd(tmp) do
            run_simulation(args)
            @test isfile("simulation_results.csv")
            return CSV.read("simulation_results.csv", DataFrame)
        end
    end
end

@testset "Deterministic Smoke + No-Drag Energy Invariant" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=400e3, ν_deg=175.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=1200.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    df = run_case(args)
    @test nrow(df) > 100
    @test all(isfinite, df.sc1_pos_1)
    @test all(isfinite, df.sc1_vel_1)

    eps = specific_energy(df, EARTH.μ)
    energy_drift = maximum(abs.(eps .- first(eps))) / abs(first(eps))
    @test energy_drift < 1e-8
end

@testset "Orbital Elements Round-Trip Invariant" begin
    ic = InitialCondition(
        ra=EARTH.Rp_e + 550e3,
        rp=EARTH.Rp_e + 300e3,
        i=37.0,
        ω=55.0,
        Ω=25.0,
        ν=130.0
    )
    r, v = orbitalelemtorv(ic, EARTH)
    oe = rvtoorbitalelement(SVector{3, Float64}(r), SVector{3, Float64}(v), EARTH)

    @test isapprox(oe[1], ic.a; rtol=1e-10, atol=1e-6)
    @test isapprox(oe[2], ic.e; rtol=1e-10, atol=1e-10)
    @test angle_distance(oe[3], ic.i) < 1e-8
    @test angle_distance(oe[4], ic.Ω) < 1e-8
    @test angle_distance(oe[5], ic.ω) < 1e-8
    @test angle_distance(oe[6], ic.ν) < 1e-8
end

@testset "Circular Orbit Robustness" begin
    ic_eq = InitialCondition(
        ra=EARTH.Rp_e + 500e3,
        rp=EARTH.Rp_e + 500e3,
        i=0.0,
        ω=0.0,
        Ω=0.0,
        ν=30.0
    )
    r, v = orbitalelemtorv(ic_eq, EARTH)
    oe = rvtoorbitalelement(SVector{3, Float64}(r), SVector{3, Float64}(v), EARTH)
    @test all(isfinite, oe)
    @test abs(oe[2]) < 1e-10

    sc = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=500e3,
        i_deg=0.0,
        ω_deg=0.0,
        Ω_deg=0.0,
        ν_deg=30.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=300.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    df = run_case(args)
    @test nrow(df) > 10
end

@testset "Quaternion Norm Invariant" begin
    q0 = normalize(SVector{4, Float64}(0.1, -0.2, 0.3, 0.9))
    w0 = SVector{3, Float64}(0.01, -0.015, 0.02)
    sc = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=420e3,
        i_deg=40.0,
        ω_deg=10.0,
        Ω_deg=20.0,
        ν_deg=175.0,
        orientation_state=(q0, w0)
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-9,
            abstol_orbit=1e-9,
            reltol_quaternion=1e-9,
            abstol_quaternion=1e-9,
            dt_max_orbit=5.0
        )
    )

    df = run_case(args)
    qnorm = sqrt.(df.sc1_q_1.^2 .+ df.sc1_q_2.^2 .+ df.sc1_q_3.^2 .+ df.sc1_q_4.^2)
    @test maximum(abs.(qnorm .- 1.0)) < 1e-6
end

@testset "Drag Dissipates Specific Orbital Energy" begin
    sc = make_spacecraft(
        ra_alt_m=220e3,
        rp_alt_m=100e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=180.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=ConstantDensityModel(1e-6, 240.0),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=140.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        keplerian=false,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_orbit=5.0,
            dt_max_atmosphere=0.2
        )
    )

    df = run_case(args)
    eps = specific_energy(df, EARTH.μ)
    @test last(eps) < first(eps) - 1e5
end
