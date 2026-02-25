using Test
using CSV
using DataFrames
using LinearAlgebra
using StaticArrays
using ComponentArrays
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

struct TimedTangentialThrusterModel <: SimulationModel.AbstractControlEffectorModel
    thrust::Float64
    direction_sign::Float64 # +1.0 prograde, -1.0 retrograde
    start_time::Float64
    stop_time::Float64
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

function SimulationModel.calcControlForceTorque(
    model::TimedTangentialThrusterModel,
    u::AbstractVector{Float64},
    p::ODEParams,
    i::Int64,
    t::Float64
)
    if t < model.start_time || t > model.stop_time
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    v = SVector{3, Float64}(u.vel)
    vm = norm(v)
    if vm == 0.0 || !isfinite(vm)
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    force = model.thrust * model.direction_sign * (v / vm)
    return force, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.calcControlEffect!(
    model::TimedTangentialThrusterModel,
    u::ComponentVector,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    return nothing
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

function make_base_thruster_model(;
    thrust::Float64,
    direction::Float64=0.0,
    Δv::Float64=0.0,
    start_burn_time::Float64=0.0,
    stop_burn_time::Float64=0.0,
    Isp::Float64=300.0
)
    return BaseThrusterModel(
        thrust=[thrust],
        direction=[direction],
        Δv=[Δv],
        start_burn_time=[start_burn_time],
        stop_burn_time=[stop_burn_time],
        Isp=[Isp]
    )
end

function build_config(;
    spacecraft::SpacecraftModel,
    density_model,
    orientation_sim::Bool,
    mission_time::Float64,
    EI_km::Float64,
    dynamic_effectors::Tuple,
    control_effectors::Tuple=(),
    control_rates::Vector{Float64}=Float64[],
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
        control_model=ControlModel(control_effectors=control_effectors, control_rates=control_rates),
        initial_time=initial_time,
        integration_tolerances=tolerances
    )
end

function build_config_multi(;
    spacecraft::Vector{SpacecraftModel},
    density_model,
    orientation_sim::Bool,
    mission_time::Float64,
    EI_km::Float64,
    dynamic_effectors::Tuple,
    control_effectors::Tuple=(),
    control_rates::Vector{Float64}=Float64[],
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
        dynamics_model=DynamicsModel(spacecraft, dynamic_effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=control_effectors, control_rates=control_rates),
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

function run_case_silent(args::SimulationConfiguration)
    return mktempdir() do tmp
        cd(tmp) do
            redirect_stdout(devnull) do
                run_simulation(args)
            end
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

    @testset "Effector Rate Validation" begin
        @test GuidanceModel((), Float64[]) isa GuidanceModel
        @test NavigationModel((), Float64[]) isa NavigationModel
        @test ControlModel((), Float64[]) isa ControlModel

        @test_throws ArgumentError GuidanceModel((:g1,), Float64[])
        @test_throws ArgumentError NavigationModel((:n1,), Float64[])
        @test_throws ArgumentError ControlModel((:c1,), Float64[])

        @test_throws ArgumentError GuidanceModel((:g1,), [0.0])
        @test_throws ArgumentError NavigationModel((:n1,), [-1.0])
        @test_throws ArgumentError ControlModel((:c1,), [Inf])
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

@testset "Deterministic Replay (No-Drag)" begin
    sc = make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    df1 = run_case_silent(args)
    df2 = run_case_silent(args)
    @test nrow(df1) == nrow(df2)
    @test nrow(df1) > 10

    sample_idxs = round.(Int, range(1, nrow(df1), length=8))
    for idx in sample_idxs
        t = Float64(df1.time[idx])
        p1 = SVector{3, Float64}(Float64(df1.sc1_pos_1[idx]), Float64(df1.sc1_pos_2[idx]), Float64(df1.sc1_pos_3[idx]))
        v1 = SVector{3, Float64}(Float64(df1.sc1_vel_1[idx]), Float64(df1.sc1_vel_2[idx]), Float64(df1.sc1_vel_3[idx]))
        p2 = SVector{3, Float64}(
            interp_linear(df2.time, df2.sc1_pos_1, t),
            interp_linear(df2.time, df2.sc1_pos_2, t),
            interp_linear(df2.time, df2.sc1_pos_3, t)
        )
        v2 = SVector{3, Float64}(
            interp_linear(df2.time, df2.sc1_vel_1, t),
            interp_linear(df2.time, df2.sc1_vel_2, t),
            interp_linear(df2.time, df2.sc1_vel_3, t)
        )
        @test norm(p1 - p2) < 0.1
        @test norm(v1 - v2) < 1e-4
    end
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

@testset "Two-Spacecraft Isolation vs Single-Craft Baselines" begin
    sc_a = make_spacecraft(ra_alt_m=520e3, rp_alt_m=430e3, ν_deg=150.0)
    sc_b = make_spacecraft(ra_alt_m=700e3, rp_alt_m=650e3, ν_deg=110.0)

    args_multi = build_config_multi(
        spacecraft=[sc_a, sc_b],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    args_a = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=430e3, ν_deg=150.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    args_b = build_config(
        spacecraft=make_spacecraft(ra_alt_m=700e3, rp_alt_m=650e3, ν_deg=110.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    df_multi = run_case_silent(args_multi)
    df_a = run_case_silent(args_a)
    df_b = run_case_silent(args_b)

    @test nrow(df_multi) > 10
    sample_idxs = round.(Int, range(1, nrow(df_multi), length=8))
    for idx in sample_idxs
        t = Float64(df_multi.time[idx])

        pa_m = SVector{3, Float64}(Float64(df_multi.sc1_pos_1[idx]), Float64(df_multi.sc1_pos_2[idx]), Float64(df_multi.sc1_pos_3[idx]))
        va_m = SVector{3, Float64}(Float64(df_multi.sc1_vel_1[idx]), Float64(df_multi.sc1_vel_2[idx]), Float64(df_multi.sc1_vel_3[idx]))
        pa_s = SVector{3, Float64}(
            interp_linear(df_a.time, df_a.sc1_pos_1, t),
            interp_linear(df_a.time, df_a.sc1_pos_2, t),
            interp_linear(df_a.time, df_a.sc1_pos_3, t)
        )
        va_s = SVector{3, Float64}(
            interp_linear(df_a.time, df_a.sc1_vel_1, t),
            interp_linear(df_a.time, df_a.sc1_vel_2, t),
            interp_linear(df_a.time, df_a.sc1_vel_3, t)
        )
        # Multi-satellite adaptive stepping can introduce modest trajectory differences vs single-body runs.
        @test norm(pa_m - pa_s) < 200.0
        @test norm(va_m - va_s) < 0.2

        pb_m = SVector{3, Float64}(Float64(df_multi.sc2_pos_1[idx]), Float64(df_multi.sc2_pos_2[idx]), Float64(df_multi.sc2_pos_3[idx]))
        vb_m = SVector{3, Float64}(Float64(df_multi.sc2_vel_1[idx]), Float64(df_multi.sc2_vel_2[idx]), Float64(df_multi.sc2_vel_3[idx]))
        pb_s = SVector{3, Float64}(
            interp_linear(df_b.time, df_b.sc1_pos_1, t),
            interp_linear(df_b.time, df_b.sc1_pos_2, t),
            interp_linear(df_b.time, df_b.sc1_pos_3, t)
        )
        vb_s = SVector{3, Float64}(
            interp_linear(df_b.time, df_b.sc1_vel_1, t),
            interp_linear(df_b.time, df_b.sc1_vel_2, t),
            interp_linear(df_b.time, df_b.sc1_vel_3, t)
        )
        @test norm(pb_m - pb_s) < 200.0
        @test norm(vb_m - vb_s) < 0.2
    end
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

@testset "Thruster Edge Cases (AI Second-Pass)" begin
    @testset "BaseThrusterModel Validation" begin
        model_ok = BaseThrusterModel(
            thrust=[1.0, 2.0],
            direction=[0.0, π],
            Δv=[10.0, 20.0],
            start_burn_time=[0.0, 0.0],
            stop_burn_time=[100.0, 120.0],
            Isp=[300.0, 300.0]
        )
        @test length(model_ok.thrust) == 2
        @test length(model_ok.direction) == 2

        @test_throws ArgumentError BaseThrusterModel(
            thrust=[1.0, 2.0],
            direction=[0.0],
            Δv=[10.0, 20.0],
            start_burn_time=[0.0, 0.0],
            stop_burn_time=[100.0, 120.0],
            Isp=[300.0, 300.0]
        )
    end

    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=120.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=600.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    p = ODEParams{1}(args=args)
    u = build_initial_conditions(args).sc[1]

    @testset "calcControlForceTorque" begin
        model = make_base_thruster_model(thrust=2.0, direction=0.0, start_burn_time=50.0, stop_burn_time=150.0)

        force_pre, torque_pre = calcControlForceTorque(model, u, p, 1, 49.9)
        @test force_pre == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test torque_pre == SVector{3, Float64}(0.0, 0.0, 0.0)

        expected_dir = normalize(SVector{3, Float64}(u.vel))
        force_start, torque_start = calcControlForceTorque(model, u, p, 1, 50.0)
        @test isapprox(force_start, 2.0 .* expected_dir; atol=1e-12, rtol=0.0)
        @test torque_start == SVector{3, Float64}(0.0, 0.0, 0.0)

        force_stop, _ = calcControlForceTorque(model, u, p, 1, 150.0)
        @test isapprox(force_stop, 2.0 .* expected_dir; atol=1e-12, rtol=0.0)

        force_post, _ = calcControlForceTorque(model, u, p, 1, 150.1)
        @test force_post == SVector{3, Float64}(0.0, 0.0, 0.0)

        model.direction[1] = π
        force_retro, _ = calcControlForceTorque(model, u, p, 1, 100.0)
        @test dot(force_retro, expected_dir) < 0.0

        u_zero = ComponentVector(pos=[0.0, 0.0, 0.0], vel=[0.0, 0.0, 0.0], mass=1.0, heat_loads=[0.0])
        force_zero, torque_zero = calcControlForceTorque(model, u_zero, p, 1, 100.0)
        @test force_zero == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test torque_zero == SVector{3, Float64}(0.0, 0.0, 0.0)

        force_oob, torque_oob = calcControlForceTorque(model, u, p, 2, 100.0)
        @test force_oob == SVector{3, Float64}(0.0, 0.0, 0.0)
        @test torque_oob == SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    @testset "calcControlEffect!" begin
        model = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
        state = build_initial_conditions(args)
        calcControlEffect!(model, state, p, 100.0, 1)
        @test isfinite(model.start_burn_time[1])
        @test isfinite(model.stop_burn_time[1])
        @test model.start_burn_time[1] < model.stop_burn_time[1]
        expected_burn_duration = state.sc[1].mass * model.Δv[1] / model.thrust[1]
        @test isapprox(
            model.stop_burn_time[1] - model.start_burn_time[1],
            expected_burn_duration;
            atol=1e-9,
            rtol=0.0
        )

        # Pre-ignition tracking: a future scheduled window should be retimed.
        model_track = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=10_000.0, stop_burn_time=11_000.0)
        calcControlEffect!(model_track, state, p, 100.0, 1)
        @test model_track.start_burn_time[1] != 10_000.0
        @test model_track.stop_burn_time[1] != 11_000.0

        # Post-ignition lock: once burn start time has been reached, keep fixed.
        s_sched = model.start_burn_time[1]
        e_sched = model.stop_burn_time[1]
        calcControlEffect!(model, state, p, 120.0, 1)
        @test model.start_burn_time[1] == s_sched
        @test model.stop_burn_time[1] == e_sched

        sc_ineligible = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=210.0)
        args_ineligible = build_config(
            spacecraft=sc_ineligible,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=600.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true
        )
        p_ineligible = ODEParams{1}(args=args_ineligible)
        state_ineligible = build_initial_conditions(args_ineligible)
        model_ineligible = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=11.0, stop_burn_time=22.0)
        calcControlEffect!(model_ineligible, state_ineligible, p_ineligible, 100.0, 1)
        @test model_ineligible.start_burn_time[1] == 11.0
        @test model_ineligible.stop_burn_time[1] == 22.0

        model_zero_thrust = make_base_thruster_model(thrust=0.0, Δv=20.0, start_burn_time=33.0, stop_burn_time=44.0)
        calcControlEffect!(model_zero_thrust, state, p, 100.0, 1)
        @test model_zero_thrust.start_burn_time[1] == 33.0
        @test model_zero_thrust.stop_burn_time[1] == 44.0

        sc_edge_block = make_spacecraft(ra_alt_m=600e3, rp_alt_m=400e3, ν_deg=180.0)
        r_edge_block, _ = orbitalelemtorv(sc_edge_block.initial_condition, EARTH)
        ei_edge_block_km = (norm(SVector{3, Float64}(r_edge_block)) - EARTH.Rp_e) / 1e3
        args_edge_block = build_config(
            spacecraft=sc_edge_block,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=600.0,
            EI_km=ei_edge_block_km,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true
        )
        p_edge_block = ODEParams{1}(args=args_edge_block)
        state_edge_block = build_initial_conditions(args_edge_block)
        model_edge_block = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
        calcControlEffect!(model_edge_block, state_edge_block, p_edge_block, 100.0, 1)
        @test model_edge_block.start_burn_time[1] == -1.0
        @test model_edge_block.stop_burn_time[1] == -1.0

        sc_edge_allow = make_spacecraft(ra_alt_m=600e3, rp_alt_m=400e3, ν_deg=170.0)
        r_edge_allow, _ = orbitalelemtorv(sc_edge_allow.initial_condition, EARTH)
        ei_edge_allow_km = (norm(SVector{3, Float64}(r_edge_allow)) - EARTH.Rp_e) / 1e3
        args_edge_allow = build_config(
            spacecraft=sc_edge_allow,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=600.0,
            EI_km=ei_edge_allow_km,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true
        )
        p_edge_allow = ODEParams{1}(args=args_edge_allow)
        state_edge_allow = build_initial_conditions(args_edge_allow)
        model_edge_allow = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
        calcControlEffect!(model_edge_allow, state_edge_allow, p_edge_allow, 100.0, 1)
        @test model_edge_allow.start_burn_time[1] != -1.0
        @test model_edge_allow.stop_burn_time[1] != -1.0

        state_hyperbolic = build_initial_conditions(args)
        rmag = norm(SVector{3, Float64}(state_hyperbolic.sc[1].pos))
        escape_speed = sqrt(2.0 * EARTH.μ / rmag)
        vhat = normalize(SVector{3, Float64}(state_hyperbolic.sc[1].vel))
        state_hyperbolic.sc[1].vel .= (1.2 * escape_speed) .* vhat
        model_hyperbolic = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=91.0, stop_burn_time=92.0)
        @test_nowarn calcControlEffect!(model_hyperbolic, state_hyperbolic, p, 100.0, 1)
        @test model_hyperbolic.start_burn_time[1] == 91.0
        @test model_hyperbolic.stop_burn_time[1] == 92.0

        state_near_parabolic = build_initial_conditions(args)
        rmag_parabolic = norm(SVector{3, Float64}(state_near_parabolic.sc[1].pos))
        vhat_parabolic = normalize(SVector{3, Float64}(state_near_parabolic.sc[1].vel))
        v_near_parabolic = 0.999999 * sqrt(2.0 * EARTH.μ / rmag_parabolic)
        state_near_parabolic.sc[1].vel .= v_near_parabolic .* vhat_parabolic
        model_near_parabolic = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=95.0, stop_burn_time=96.0)
        @test_nowarn calcControlEffect!(model_near_parabolic, state_near_parabolic, p, 100.0, 1)
        @test isfinite(model_near_parabolic.start_burn_time[1])
        @test isfinite(model_near_parabolic.stop_burn_time[1])

        state_singular = build_initial_conditions(args)
        state_singular.sc[1].pos .= 0.0
        state_singular.sc[1].vel .= 0.0
        model_singular = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=93.0, stop_burn_time=94.0)
        @test_nowarn calcControlEffect!(model_singular, state_singular, p, 100.0, 1)
        @test model_singular.start_burn_time[1] == 93.0
        @test model_singular.stop_burn_time[1] == 94.0

        state_tiny_a = build_initial_conditions(args)
        state_tiny_a.sc[1].pos .= SVector{3, Float64}(1.0, 0.0, 0.0)
        state_tiny_a.sc[1].vel .= SVector{3, Float64}(0.0, 0.1, 0.0)
        model_tiny_a = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=97.0, stop_burn_time=98.0)
        @test_nowarn calcControlEffect!(model_tiny_a, state_tiny_a, p, 100.0, 1)
        @test isfinite(model_tiny_a.start_burn_time[1])
        @test isfinite(model_tiny_a.stop_burn_time[1])

        sc_multi_1 = make_spacecraft(ra_alt_m=650e3, rp_alt_m=450e3, ν_deg=120.0)
        sc_multi_2 = make_spacecraft(ra_alt_m=700e3, rp_alt_m=500e3, ν_deg=150.0)
        args_multi = build_config_multi(
            spacecraft=[sc_multi_1, sc_multi_2],
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=600.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true
        )
        p_multi = ODEParams{2}(args=args_multi)
        state_multi = build_initial_conditions(args_multi)
        model_multi = BaseThrusterModel(
            thrust=[2.0, 3.0],
            direction=[0.0, π],
            Δv=[20.0, 10.0],
            start_burn_time=[-1.0, -2.0],
            stop_burn_time=[-3.0, -4.0],
            Isp=[300.0, 300.0]
        )
        calcControlEffect!(model_multi, state_multi, p_multi, 100.0, 1)
        @test model_multi.start_burn_time[1] != -1.0
        @test model_multi.stop_burn_time[1] != -3.0
        @test model_multi.start_burn_time[2] == -2.0
        @test model_multi.stop_burn_time[2] == -4.0

        s1 = model_multi.start_burn_time[1]
        e1 = model_multi.stop_burn_time[1]
        calcControlEffect!(model_multi, state_multi, p_multi, 100.0, 2)
        @test model_multi.start_burn_time[1] == s1
        @test model_multi.stop_burn_time[1] == e1
        @test model_multi.start_burn_time[2] != -2.0
        @test model_multi.stop_burn_time[2] != -4.0

        model_oob = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
        @test_nowarn calcControlEffect!(model_oob, state, p, 100.0, 2)
        @test model_oob.start_burn_time[1] == -1.0
        @test model_oob.stop_burn_time[1] == -1.0
    end

    @testset "schmitt_trigger" begin
        @test schmitt_trigger(0.80, 0.75, 0.25) == 1.0
        @test schmitt_trigger(0.20, 0.75, 0.25) == 0.0
        @test schmitt_trigger(0.75, 0.75, 0.25) == 0.0
        @test schmitt_trigger(0.50, 0.75, 0.25) == 0.0

        # Sequential behavior check: current implementation is memoryless.
        state_high = schmitt_trigger(0.80, 0.75, 0.25)
        state_mid_after_high = schmitt_trigger(0.50, 0.75, 0.25)
        _ = schmitt_trigger(0.20, 0.75, 0.25)
        state_mid_after_low = schmitt_trigger(0.50, 0.75, 0.25)
        @test state_high == 1.0
        @test state_mid_after_high == 0.0
        @test state_mid_after_low == 0.0
        @test state_mid_after_high == state_mid_after_low
    end

    @testset "integrate_impulse!" begin
        link = Link{0}(root=true, attitude_control_rate=0.1)
        ω = 20.0

        thr_zero = Thruster(max_thrust=10.0, cutoff_frequency=ω, κ=0.0)
        impulse_zero = integrate_impulse!(link, thr_zero, 0.0, 0.0)
        @test impulse_zero == 0.0
        @test thr_zero.κ == 0.0

        thr_full = Thruster(max_thrust=10.0, cutoff_frequency=ω, κ=0.0)
        dt = link.attitude_control_rate
        expected_full = thr_full.max_thrust * (dt + (-1.0) / ω * (1 - exp(-ω * dt)))
        impulse_full = integrate_impulse!(link, thr_full, dt, 0.0)
        @test isapprox(impulse_full, expected_full; atol=1e-12, rtol=0.0)
        @test isapprox(thr_full.κ, 1 - exp(-ω * dt); atol=1e-12, rtol=0.0)

        thr_decay = Thruster(max_thrust=10.0, cutoff_frequency=ω, κ=1.0)
        expected_decay = thr_decay.max_thrust / ω * (1 - exp(-ω * dt))
        impulse_decay = integrate_impulse!(link, thr_decay, 0.0, 0.0)
        @test isapprox(impulse_decay, expected_decay; atol=1e-12, rtol=0.0)
        @test isapprox(thr_decay.κ, exp(-ω * dt); atol=1e-12, rtol=0.0)

        thr_small_ω = Thruster(max_thrust=10.0, cutoff_frequency=1e-12, κ=0.3)
        impulse_small_ω = integrate_impulse!(link, thr_small_ω, 0.05, 0.0)
        @test isfinite(impulse_small_ω)
        @test impulse_small_ω >= 0.0
        @test isfinite(thr_small_ω.κ)

        thr_ref_neg = Thruster(max_thrust=10.0, cutoff_frequency=ω, κ=0.4)
        thr_neg = deepcopy(thr_ref_neg)
        thr_zero_ref = deepcopy(thr_ref_neg)
        impulse_neg = integrate_impulse!(link, thr_neg, -0.01, 0.0)
        impulse_zero_ref = integrate_impulse!(link, thr_zero_ref, 0.0, 0.0)
        @test isapprox(impulse_neg, impulse_zero_ref; atol=1e-12, rtol=0.0)
        @test isapprox(thr_neg.κ, thr_zero_ref.κ; atol=1e-12, rtol=0.0)

        thr_ref_hi = Thruster(max_thrust=10.0, cutoff_frequency=ω, κ=0.4)
        thr_hi = deepcopy(thr_ref_hi)
        thr_dt_ref = deepcopy(thr_ref_hi)
        impulse_hi = integrate_impulse!(link, thr_hi, 10.0 * link.attitude_control_rate, 0.0)
        impulse_dt_ref = integrate_impulse!(link, thr_dt_ref, link.attitude_control_rate, 0.0)
        @test isapprox(impulse_hi, impulse_dt_ref; atol=1e-12, rtol=0.0)
        @test isapprox(thr_hi.κ, thr_dt_ref.κ; atol=1e-12, rtol=0.0)
    end

    @testset "thrust_calculation_schmitt_trigger!" begin
        mktempdir() do tmp
            cd(tmp) do
                link = Link{0}(root=true, attitude_control_rate=0.1)
                thr = Thruster(
                    max_thrust=10.0,
                    min_firing_time=0.05,
                    level_on=0.75,
                    level_off=0.25,
                    cutoff_frequency=100.0,
                    κ=0.0
                )

                thrust_calculation_schmitt_trigger!(link, thr, 1.0, 0.0)
                @test thr.thrust == 0.0
                @test isfile("thruster_debug.csv")

                thrust_calculation_schmitt_trigger!(link, thr, 9.0, 0.1)
                @test thr.thrust > 0.0
                @test thr.thrust <= thr.max_thrust + 1e-9

                thr_zero_max = Thruster(max_thrust=0.0, cutoff_frequency=100.0, κ=0.0)
                @test_nowarn thrust_calculation_schmitt_trigger!(link, thr_zero_max, 1.0, 0.2)
                @test thr_zero_max.thrust == 0.0
            end
        end
    end

    @testset "update_thrusters!" begin
        mktempdir() do tmp
            cd(tmp) do
                link = Link{0}(root=true, attitude_control_rate=0.1)
                update_thrusters!(link, SVector{3, Float64}(0.0, 0.0, 0.0), 0.0)
                @test isempty(link.thrusters)
                @test size(link.J_thruster) == (3, 0)

                thr = Thruster(
                    max_thrust=1.0,
                    cutoff_frequency=100.0,
                    min_firing_time=0.0,
                    location=MVector{3, Float64}(0.0, 1.0, 0.0),
                    direction=MVector{3, Float64}(0.0, 0.0, 1.0)
                )
                push!(link.thrusters, thr)

                update_thrusters!(link, SVector{3, Float64}(-1.0, 0.0, 0.0), 0.1)
                @test link.thrusters[1].thrust >= 0.0

                link_full = Link{0}(root=true, attitude_control_rate=0.1)
                push!(link_full.thrusters, Thruster(max_thrust=100.0, cutoff_frequency=1e6, min_firing_time=0.0, location=MVector{3, Float64}(0.0, 1.0, 0.0), direction=MVector{3, Float64}(0.0, 0.0, 1.0)))
                push!(link_full.thrusters, Thruster(max_thrust=100.0, cutoff_frequency=1e6, min_firing_time=0.0, location=MVector{3, Float64}(0.0, 0.0, 1.0), direction=MVector{3, Float64}(1.0, 0.0, 0.0)))
                push!(link_full.thrusters, Thruster(max_thrust=100.0, cutoff_frequency=1e6, min_firing_time=0.0, location=MVector{3, Float64}(1.0, 0.0, 0.0), direction=MVector{3, Float64}(0.0, 1.0, 0.0)))
                τ_req_full = SVector{3, Float64}(1.0, 2.0, 3.0)
                update_thrusters!(link_full, τ_req_full, 0.2)
                @test rank(link_full.J_thruster) == 3
                @test all(isfinite, link_full.J_thruster)
                @test all(thr_i -> isfinite(thr_i.thrust) && thr_i.thrust >= 0.0, link_full.thrusters)
                @test any(thr_i -> thr_i.thrust > 0.0, link_full.thrusters)
                thrust_full = [thr_i.thrust for thr_i in link_full.thrusters]
                τ_ach_full = link_full.J_thruster * thrust_full
                @test norm(τ_ach_full - τ_req_full) / norm(τ_req_full) < 1e-6

                link_singular = Link{0}(root=true, attitude_control_rate=0.1)
                push!(link_singular.thrusters, Thruster(max_thrust=100.0, cutoff_frequency=1e6, min_firing_time=0.0, location=MVector{3, Float64}(0.0, 1.0, 0.0), direction=MVector{3, Float64}(0.0, 0.0, 1.0)))
                push!(link_singular.thrusters, Thruster(max_thrust=100.0, cutoff_frequency=1e6, min_firing_time=0.0, location=MVector{3, Float64}(0.0, 2.0, 0.0), direction=MVector{3, Float64}(0.0, 0.0, 1.0)))
                τ_req_singular = SVector{3, Float64}(0.0, 1.0, 0.0)
                update_thrusters!(link_singular, τ_req_singular, 0.3)
                @test rank(link_singular.J_thruster) == 1
                @test all(isfinite, link_singular.J_thruster)
                @test all(thr_i -> isfinite(thr_i.thrust) && thr_i.thrust >= 0.0, link_singular.thrusters)
                thrust_singular = [thr_i.thrust for thr_i in link_singular.thrusters]
                τ_ach_singular = link_singular.J_thruster * thrust_singular
                τ_lsq_singular = link_singular.J_thruster * (pinv(link_singular.J_thruster) * τ_req_singular)
                @test norm(τ_ach_singular - τ_lsq_singular) < 1e-6

                link_degenerate = Link{0}(root=true, attitude_control_rate=0.1)
                push!(
                    link_degenerate.thrusters,
                    Thruster(
                        max_thrust=50.0,
                        cutoff_frequency=100.0,
                        min_firing_time=0.0,
                        location=MVector{3, Float64}(0.5, 0.0, 0.0),
                        direction=MVector{3, Float64}(0.0, 0.0, 0.0)
                    )
                )
                @test_nowarn update_thrusters!(link_degenerate, SVector{3, Float64}(0.5, 0.0, 0.0), 0.4)
                @test all(isfinite, link_degenerate.J_thruster)
                @test link_degenerate.J_thruster[:, 1] == zeros(3)
                @test isfinite(link_degenerate.thrusters[1].thrust)
                @test link_degenerate.thrusters[1].thrust >= 0.0
            end
        end
    end
end

@testset "Control Burn Energy Sign (End-to-End)" begin
    sc0 = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
    args0 = build_config(
        spacecraft=sc0,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=250.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(),
        control_rates=Float64[],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    df0 = run_case_silent(args0)
    eps0 = specific_energy(df0, EARTH.μ)
    Δeps0 = last(eps0) - first(eps0)
    @test abs(Δeps0) < 2e3

    scp = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
    thr_pro = TimedTangentialThrusterModel(800.0, +1.0, 100.0, 101.0)
    argsp = build_config(
        spacecraft=scp,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=250.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(thr_pro,),
        control_rates=[1.0],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    dfp = run_case_silent(argsp)
    epsp = specific_energy(dfp, EARTH.μ)
    Δepsp = last(epsp) - first(epsp)
    @test Δepsp > 5e3

    scr = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
    thr_retro = TimedTangentialThrusterModel(800.0, -1.0, 100.0, 101.0)
    argsr = build_config(
        spacecraft=scr,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=250.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(thr_retro,),
        control_rates=[1.0],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    dfr = run_case_silent(argsr)
    epsr = specific_energy(dfr, EARTH.μ)
    Δepsr = last(epsr) - first(epsr)
    @test Δepsr < -5e3
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
