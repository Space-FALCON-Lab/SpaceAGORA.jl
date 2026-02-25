using Test
using CSV
using DataFrames
using LinearAlgebra
using StaticArrays
using SPICE
using JET
using Aqua

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

function make_single_link_spacecraft(;
    ra_alt_m::Float64,
    rp_alt_m::Float64,
    i_deg::Float64=35.0,
    ω_deg::Float64=40.0,
    Ω_deg::Float64=10.0,
    ν_deg::Float64=175.0
)
    root = Link{0}(root=true, m=500.0, ref_area=12.0)
    ic = InitialCondition(
        ra=EARTH.Rp_e + ra_alt_m,
        rp=EARTH.Rp_e + rp_alt_m,
        i=i_deg,
        ω=ω_deg,
        Ω=Ω_deg,
        ν=ν_deg
    )

    return SpacecraftModel(
        Joint[],
        [root],
        root,
        true,
        root.m,
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
    tolerances::IntegrationTolerances=IntegrationTolerances(),
    initial_time::InitialTime=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0)
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
        initial_time=initial_time,
        integration_tolerances=tolerances
    )
end

function interp_linear(times::AbstractVector{<:Real}, values::AbstractVector{<:Real}, t::Float64)
    idx = searchsortedlast(times, t)
    if idx <= 0
        return Float64(values[1])
    elseif idx >= length(times)
        return Float64(values[end])
    end
    t1 = Float64(times[idx])
    t2 = Float64(times[idx + 1])
    y1 = Float64(values[idx])
    y2 = Float64(values[idx + 1])
    α = (t - t1) / (t2 - t1)
    return (1.0 - α) * y1 + α * y2
end

function make_agora_earth_spacecraft()
    main_bus = Link{0}(root=true, m=620.0, ref_area=2.05 * 2.8)
    left_panel = Link{0}(root=false, m=10.0, ref_area=5.7 * 1.0 / 2.0, r=MVector{3, Float64}(0.0, -2.05 / 2.0 - 5.7 / 4.0, 0.0))
    right_panel = Link{0}(root=false, m=10.0, ref_area=5.7 * 1.0 / 2.0, r=MVector{3, Float64}(0.0, 2.05 / 2.0 + 5.7 / 4.0, 0.0))
    ic = InitialCondition(
        ra=56_378.7978559e3,
        rp=EARTH.Rp_e + 200_590.0,
        i=89.876,
        ω=75.505,
        Ω=104.115,
        ν=175.0
    )

    return SpacecraftModel(
        Joint[],
        [main_bus, left_panel, right_panel],
        main_bus,
        true,
        main_bus.m + left_panel.m + right_panel.m,
        200.0,
        main_bus.inertia,
        0,
        0,
        ic,
        1
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

@testset "API Convenience Constructors" begin
    ic = InitialCondition()
    @test ic isa InitialCondition
    @test ic.a == 0.0
    @test ic.e == 0.0

    link = Link()
    @test link isa Link{0}

    joint = Joint()
    @test joint isa Joint

    sc = SpacecraftModel()
    @test sc isa SpacecraftModel
    @test sc.root.root

    custom_root = Link{0}(root=true, m=100.0)
    sc_custom = SpacecraftModel(root=custom_root, id=42)
    @test sc_custom.id == 42
    @test sc_custom.root === custom_root
    @test any(link -> link === custom_root, sc_custom.links)

    nbody = NBodyGravityModel(["Sun"], "Earth", SPICE_PATH)
    @test nbody isa NBodyGravityModel
    @test nbody.primary_body_name == "Earth"
    @test nbody.body_names == ("Sun",)
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

@testset "AGORA Earth Regression (Golden)" begin
    golden_path = joinpath(@__DIR__, "golden", "agora_earth_regression.csv")
    @test isfile(golden_path)
    golden = CSV.read(golden_path, DataFrame)

    sc = make_agora_earth_spacecraft()
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=6.0 * 3600.0,
        EI_km=300.0,
        dynamic_effectors=(InverseSquaredJ2GravityModel(),),
        keplerian=true,
        initial_time=InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0),
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            dt_max_orbit=15.0
        )
    )

    df = run_case(args)
    @test nrow(df) > 1000
    times = Vector{Float64}(df.time)

    for row in eachrow(golden)
        t = Float64(row.time)
        pos_atol = Float64(row.pos_atol_m)
        vel_atol = Float64(row.vel_atol_mps)

        pos1 = interp_linear(times, df.sc1_pos_1, t)
        pos2 = interp_linear(times, df.sc1_pos_2, t)
        pos3 = interp_linear(times, df.sc1_pos_3, t)
        vel1 = interp_linear(times, df.sc1_vel_1, t)
        vel2 = interp_linear(times, df.sc1_vel_2, t)
        vel3 = interp_linear(times, df.sc1_vel_3, t)

        @test isapprox(pos1, Float64(row.pos1); atol=pos_atol, rtol=0.0)
        @test isapprox(pos2, Float64(row.pos2); atol=pos_atol, rtol=0.0)
        @test isapprox(pos3, Float64(row.pos3); atol=pos_atol, rtol=0.0)
        @test isapprox(vel1, Float64(row.vel1); atol=vel_atol, rtol=0.0)
        @test isapprox(vel2, Float64(row.vel2); atol=vel_atol, rtol=0.0)
        @test isapprox(vel3, Float64(row.vel3); atol=vel_atol, rtol=0.0)
    end
end

@testset "N-Body Gravity Typed API Smoke" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), NBodyGravityModel(["Sun"], "Earth", SPICE_PATH)),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    df = run_case(args)
    @test nrow(df) > 10
    @test all(isfinite, df.sc1_pos_1)
    @test all(isfinite, df.sc1_vel_1)
end

@testset "Single-Link Drag Dissipates Specific Orbital Energy" begin
    sc = make_single_link_spacecraft(
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

@testset "JET Static Analysis" begin
    JET.@test_opt InitialCondition()
    JET.@test_opt Link()
    JET.@test_opt Joint()
    JET.@test_opt SpacecraftModel()

    sc = SpacecraftModel()
    JET.@test_opt rotate_to_inertial(sc, sc.root, 1)
end

@testset "Aqua Package Quality" begin
    Aqua.test_all(
        SimulationModel;
        ambiguities=false,
        stale_deps=false,
        deps_compat=false,
        project_extras=false,
        piracies=false,
        persistent_tasks=false,
        undocumented_names=false
    )
end
