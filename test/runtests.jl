using Test
using CSV
using DataFrames
using LinearAlgebra
using Logging
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
# Legacy wrappers (run_analysis/aerobraking_campaign).
include(joinpath(REPO_ROOT, "src", "simulation", "Run.jl"))

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
    simulation_settings::SimulationSettings=SimulationSettings(results=true, verbose=false, generate_plots=false, normalize=false),
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
        simulation_settings=simulation_settings,
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
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
    simulation_settings::SimulationSettings=SimulationSettings(results=true, verbose=false, generate_plots=false, normalize=false),
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
        simulation_settings=simulation_settings,
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
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

function run_case(args::SimulationConfiguration; isolate_state::Bool=true)
    return mktempdir() do tmp
        cd(tmp) do
            run_simulation(args; isolate_state=isolate_state)
            @test isfile("simulation_results.csv")
            return CSV.read("simulation_results.csv", DataFrame)
        end
    end
end

function run_case_silent(args::SimulationConfiguration; isolate_state::Bool=true)
    return mktempdir() do tmp
        cd(tmp) do
            redirect_stdout(devnull) do
                run_simulation(args; isolate_state=isolate_state)
            end
            @test isfile("simulation_results.csv")
            return CSV.read("simulation_results.csv", DataFrame)
        end
    end
end

function run_case_capture_stdout(args::SimulationConfiguration; expect_results_csv::Bool=true, isolate_state::Bool=true)
    return mktempdir() do tmp
        cd(tmp) do
            output = ""
            mktemp() do path, io
                redirect_stdout(io) do
                    run_simulation(args; isolate_state=isolate_state)
                end
                flush(io)
                seekstart(io)
                output = read(io, String)
            end
            if expect_results_csv
                @test isfile("simulation_results.csv")
                return CSV.read("simulation_results.csv", DataFrame), output
            else
                @test !isfile("simulation_results.csv")
                return DataFrame(), output
            end
        end
    end
end

function run_case_via_run_analysis(args::SimulationConfiguration; expect_results_csv::Bool=true, isolate_state::Bool=true)
    return mktempdir() do tmp
        cd(tmp) do
            redirect_stdout(devnull) do
                run_analysis(args; isolate_state=isolate_state)
            end
            if expect_results_csv
                @test isfile("simulation_results.csv")
                return CSV.read("simulation_results.csv", DataFrame)
            else
                @test !isfile("simulation_results.csv")
                return DataFrame()
            end
        end
    end
end

function run_case_via_campaign(args::SimulationConfiguration; expect_results_csv::Bool=true, isolate_state::Bool=true, state=nothing)
    return mktempdir() do tmp
        cd(tmp) do
            redirect_stdout(devnull) do
                if state === nothing
                    aerobraking_campaign(args; isolate_state=isolate_state)
                else
                    aerobraking_campaign(args, state; isolate_state=isolate_state)
                end
            end
            if expect_results_csv
                @test isfile("simulation_results.csv")
                return CSV.read("simulation_results.csv", DataFrame)
            else
                @test !isfile("simulation_results.csv")
                return DataFrame()
            end
        end
    end
end

function run_case_via_run_orbitalelements(args::SimulationConfiguration; expect_results_csv::Bool=true, isolate_state::Bool=true)
    return mktempdir() do tmp
        cd(tmp) do
            redirect_stdout(devnull) do
                run_orbitalelements(args; isolate_state=isolate_state)
            end
            if expect_results_csv
                @test isfile("simulation_results.csv")
                return CSV.read("simulation_results.csv", DataFrame)
            else
                @test !isfile("simulation_results.csv")
                return DataFrame()
            end
        end
    end
end

function run_case_via_run_vgamma(args::SimulationConfiguration; expect_results_csv::Bool=true, isolate_state::Bool=true)
    return mktempdir() do tmp
        cd(tmp) do
            redirect_stdout(devnull) do
                run_vgamma(args; isolate_state=isolate_state)
            end
            if expect_results_csv
                @test isfile("simulation_results.csv")
                return CSV.read("simulation_results.csv", DataFrame)
            else
                @test !isfile("simulation_results.csv")
                return DataFrame()
            end
        end
    end
end

function seed_solution_for_save_csv!(
    solution::Solution;
    n_bodies::Int,
    n_reaction_wheels::Int,
    n_thrusters::Int,
    base::Float64=1.0,
    closed_form::Bool=false
)
    solution.physical_properties.α = [Float64[] for _ in 1:n_bodies]
    solution.physical_properties.β = [Float64[] for _ in 1:n_bodies]
    solution.performance.heat_rate = [Float64[] for _ in 1:n_bodies]
    solution.performance.heat_load = [Float64[] for _ in 1:n_bodies]
    solution.physical_properties.rw_h = [Float64[] for _ in 1:n_reaction_wheels]
    solution.physical_properties.rw_τ = [Float64[] for _ in 1:n_reaction_wheels]
    solution.physical_properties.thruster_forces = [Float64[] for _ in 1:n_thrusters]

    function _push_sample!(field)
        if field isa Vector{Float64}
            push!(field, base)
        elseif field isa Vector{Int64}
            push!(field, 1)
        elseif field isa Vector{Vector{Float64}}
            for subfield in field
                push!(subfield, base)
            end
        end
        return nothing
    end

    for group in (solution.orientation, solution.physical_properties, solution.performance, solution.forces)
        for fname in fieldnames(typeof(group))
            _push_sample!(getfield(group, fname))
        end
    end

    if closed_form
        solution.closed_form.t_cf = [base + 10.0]
        solution.closed_form.h_cf = [base + 20.0]
        solution.closed_form.γ_cf = [base + 30.0]
        solution.closed_form.v_cf = [base + 40.0]
    end

    return solution
end

const LEGACY_TARGETING_LOADED = Ref(false)
module LegacyTargetingSandbox
using ..SimulationModel
end
const LEGACY_TARGETING_SANDBOX = LegacyTargetingSandbox

const LEGACY_CONFIG_LOADED = Ref(false)
module LegacyConfigSandbox
using ..SimulationModel
end
const LEGACY_CONFIG_SANDBOX = LegacyConfigSandbox

module LegacyRuntimeSandbox
using ..SimulationModel
end
const LEGACY_RUNTIME_SANDBOX = LegacyRuntimeSandbox

const LEGACY_COMPLETE_PASSAGE_LOADED = Ref(false)
module LegacyCompletePassageSandbox
using ..SimulationModel
end
const LEGACY_COMPLETE_PASSAGE_SANDBOX = LegacyCompletePassageSandbox

const GUIDANCE_SANDBOX_LOADED = Ref(false)
module GuidanceSandbox
using ..SimulationModel
using ..SimulationModel.AbstractTypes: AbstractGuidanceModel
using ComponentArrays
end
const GUIDANCE_SANDBOX = GuidanceSandbox

const LEGACY_MONTE_CARLO_SANDBOX_LOADED = Ref(false)
module LegacyMonteCarloSandbox
using ..SimulationModel
using CSV
using DataFrames
using Statistics
end
const LEGACY_MONTE_CARLO_SANDBOX = LegacyMonteCarloSandbox

const LEGACY_CLOSED_FORM_SANDBOX_LOADED = Ref(false)
module LegacyClosedFormSandbox
using ..SimulationModel
end
const LEGACY_CLOSED_FORM_SANDBOX = LegacyClosedFormSandbox

module IncludeOrderSandbox
end
const INCLUDE_ORDER_SANDBOX = IncludeOrderSandbox

module ExportImportSandbox
end
const EXPORT_IMPORT_SANDBOX = ExportImportSandbox

@testset "SimulationModel Export Contract" begin
    sandbox = EXPORT_IMPORT_SANDBOX
    @test_nowarn Base.include(sandbox, joinpath(REPO_ROOT, "src", "simulation_model", "SimulationModel.jl"))
    @test_nowarn Core.eval(sandbox, :(using .SimulationModel))

    required_public_names = [
        :SimulationConfiguration,
        :InitialCondition,
        :InitialTime,
        :MissionConfiguration,
        :EnvironmentModel,
        :IntegrationTolerances,
        :ControlModel,
        :GuidanceModel,
        :NavigationModel,
        :DynamicsModel
    ]

    for sym in required_public_names
        @test Core.eval(sandbox, :(isdefined(@__MODULE__, $(QuoteNode(sym)))))
    end
end

@testset "Include-Order + Name Ambiguity Smoke" begin
    sandbox = INCLUDE_ORDER_SANDBOX

    Base.include_string(sandbox, """
    module ConflictingExports
        export SimulationConfiguration
        struct SimulationConfiguration end
    end
    using .ConflictingExports
    """)

    @test_nowarn Base.include(sandbox, joinpath(REPO_ROOT, "src", "simulation_model", "SimulationModel.jl"))
    Core.eval(sandbox, :(const quat_mult = SimulationModel.quat_mult))
    @test_nowarn Base.include(sandbox, joinpath(REPO_ROOT, "src", "simulation", "run_simulation.jl"))
    @test isdefined(sandbox, :run_simulation)
end

@testset "Complete Passage Typed Entry Contract Smoke" begin
    module_name = gensym(:CompletePassageContractSandbox)
    Core.eval(Main, :(module $module_name end))
    sandbox = getfield(Main, module_name)

    @test_nowarn Base.include(sandbox, joinpath(REPO_ROOT, "src", "simulation_model", "SimulationModel.jl"))
    Core.eval(sandbox, quote
        const _asim_contract_isolate_state = Ref{Union{Nothing, Bool}}(nothing)
        function run_simulation(args::SimulationModel.SimulationConfiguration; isolate_state::Bool=true)
            _asim_contract_isolate_state[] = isolate_state
            return :asim_forwarded
        end
    end)

    complete_src = read(joinpath(REPO_ROOT, "src", "simulation", "Complete_passage.jl"), String)
    start_idx = findfirst("function asim(args::SimulationConfiguration; isolate_state::Bool=true)", complete_src)
    next_idx = findfirst("function asim(initial_state, numberofpassage, args, params)", complete_src)
    @test start_idx !== nothing
    @test next_idx !== nothing

    typed_asim_src = strip(complete_src[first(start_idx):(first(next_idx) - 1)])
    @test occursin("return run_simulation(args; isolate_state=isolate_state)", typed_asim_src)

    Core.eval(sandbox, :(using .SimulationModel))
    Base.include_string(sandbox, typed_asim_src)
    @test isdefined(sandbox, :asim)

    typed_args = Core.eval(sandbox, quote
        planet = SimulationModel.Earth("", joinpath(Main.REPO_ROOT, "GRAM_Data", "SPICE"))
        env = SimulationModel.EnvironmentModel(
            planet=planet,
            EI=120.0,
            density_model=SimulationModel.NoAtmosphereModel(),
            topography=false,
            topo_degree=0,
            topo_order=0,
            wind=false,
            thermal_model=SimulationModel.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet)
        )
        sc = SimulationModel.SpacecraftModel()
        SimulationModel.SimulationConfiguration(
            simulation_settings=SimulationModel.SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
            mission_configuration=SimulationModel.MissionConfiguration(
                mission_type=SimulationModel.MissionTime,
                keplerian=true,
                number_of_orbits=1,
                mission_time=60.0,
                orientation_sim=false,
                num_steps_to_save=10
            ),
            environment_model=env,
            dynamics_model=SimulationModel.DynamicsModel([sc], (SimulationModel.InverseSquaredGravityModel(),)),
            guidance_model=SimulationModel.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
            navigation_model=SimulationModel.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
            control_model=SimulationModel.ControlModel(control_effectors=(), control_rates=Float64[]),
            initial_time=SimulationModel.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
            integration_tolerances=SimulationModel.IntegrationTolerances(reltol_orbit=1e-8, abstol_orbit=1e-8)
        )
    end)

    result = getfield(sandbox, :asim)(typed_args; isolate_state=false)
    @test result == :asim_forwarded
    @test Core.eval(sandbox, :(_asim_contract_isolate_state[])) == false
end

function ensure_legacy_targeting_loaded!()
    if LEGACY_TARGETING_LOADED[]
        return
    end

    Core.eval(LEGACY_TARGETING_SANDBOX, :(include(joinpath(Main.REPO_ROOT, "src", "control", "Control.jl"))))
    Core.eval(LEGACY_TARGETING_SANDBOX, :(include(joinpath(Main.REPO_ROOT, "src", "control", "targeting_control", "targeting.jl"))))
    Core.eval(LEGACY_TARGETING_SANDBOX, :(include(joinpath(Main.REPO_ROOT, "src", "control", "targeting_control", "Eom_targeting.jl"))))

    LEGACY_TARGETING_LOADED[] = true
end

function ensure_legacy_config_loaded!()
    if LEGACY_CONFIG_LOADED[]
        return
    end

    Core.eval(LEGACY_CONFIG_SANDBOX, :(include(joinpath(Main.REPO_ROOT, "src", "physical_models", "Mission.jl"))))
    Core.eval(LEGACY_CONFIG_SANDBOX, :(include(joinpath(Main.REPO_ROOT, "src", "physical_models", "Planet_data.jl"))))
    Core.eval(LEGACY_CONFIG_SANDBOX, :(include(joinpath(Main.REPO_ROOT, "src", "utils", "Define_mission.jl"))))

    LEGACY_CONFIG_LOADED[] = true
end

function ensure_legacy_complete_passage_loaded!()
    if LEGACY_COMPLETE_PASSAGE_LOADED[]
        return
    end

    source = read(joinpath(REPO_ROOT, "src", "simulation", "Complete_passage.jl"), String)
    typed_start = findfirst("function asim(args::SimulationConfiguration; isolate_state::Bool=true)", source)
    legacy_start = findfirst("function asim(initial_state, numberofpassage, args, params)", source)
    typed_start === nothing && throw(ArgumentError("Typed asim entrypoint missing in Complete_passage.jl"))
    legacy_start === nothing && throw(ArgumentError("Legacy asim entrypoint missing in Complete_passage.jl"))

    typed_method = strip(source[first(typed_start):(first(legacy_start) - 1)])
    legacy_method = strip(source[first(legacy_start):end])

    Core.eval(LEGACY_COMPLETE_PASSAGE_SANDBOX, :(using .SimulationModel))
    Base.include_string(LEGACY_COMPLETE_PASSAGE_SANDBOX, typed_method * "\n\n" * legacy_method)
    LEGACY_COMPLETE_PASSAGE_LOADED[] = true
end

function ensure_guidance_sandbox_loaded!()
    if GUIDANCE_SANDBOX_LOADED[]
        return
    end

    Core.eval(GUIDANCE_SANDBOX, :(include(joinpath(Main.REPO_ROOT, "src", "guidance", "Thruster_guidance_models.jl"))))
    Core.eval(GUIDANCE_SANDBOX, :(include(joinpath(Main.REPO_ROOT, "src", "guidance", "Thruster_guidance_functions.jl"))))
    GUIDANCE_SANDBOX_LOADED[] = true
end

function ensure_legacy_monte_carlo_loaded!()
    if LEGACY_MONTE_CARLO_SANDBOX_LOADED[]
        return
    end

    Core.eval(LEGACY_MONTE_CARLO_SANDBOX, :(include(joinpath(Main.REPO_ROOT, "src", "utils", "MonteCarlo_set.jl"))))
    LEGACY_MONTE_CARLO_SANDBOX_LOADED[] = true
end

function ensure_legacy_closed_form_loaded!()
    if LEGACY_CLOSED_FORM_SANDBOX_LOADED[]
        return
    end

    Core.eval(LEGACY_CLOSED_FORM_SANDBOX, :(include(joinpath(Main.REPO_ROOT, "src", "utils", "Closed_form_solution.jl"))))
    LEGACY_CLOSED_FORM_SANDBOX_LOADED[] = true
end

@testset "Complete Passage Legacy Entrypoint Smoke" begin
    ensure_legacy_complete_passage_loaded!()
    sandbox = LEGACY_COMPLETE_PASSAGE_SANDBOX

    args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )

    @test_throws ArgumentError sandbox.asim(args; isolate_state=true)
    @test_throws UndefVarError sandbox.asim(nothing, 1, args, ())
end

@testset "API Convenience Constructors" begin
    @testset "Mission/Environment Config Validation" begin
        mc_default = MissionConfiguration()
        @test mc_default.mission_type == MissionTime

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            SimulationModel.SimConfig._deprecated_mission_type_input_warned[] = false
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            mc_str = MissionConfiguration(mission_type="Time", mission_time=600.0, number_of_orbits=1, num_steps_to_save=10)
            @test mc_str.mission_type == MissionTime
            @test mc_str.mission_type == "Time"
            @test "Time" == mc_str.mission_type
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            mc_sym = MissionConfiguration(mission_type=:orbits, mission_time=600.0, number_of_orbits=2, num_steps_to_save=10)
            @test mc_sym.mission_type == MissionOrbits
            @test mc_sym.mission_type == :orbits
            @test :orbits == mc_sym.mission_type
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            mc_enum = MissionConfiguration(mission_type=MissionOrbits, mission_time=600.0, number_of_orbits=3, num_steps_to_save=10)
            @test mc_enum.mission_type == MissionOrbits
            @test mc_enum.mission_type == "Orbits"
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "1") do
            SimulationModel.SimConfig._deprecated_mission_type_input_warned[] = false
            @test_logs (:warn, r"mission_type=.*deprecated") MissionConfiguration(mission_type="Time")
            @test SimulationModel.SimConfig._deprecated_mission_type_input_warned[] == true
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            SimulationModel.SimConfig._deprecated_mission_type_input_warned[] = false
            @test_logs MissionConfiguration(mission_type="Time")
            @test SimulationModel.SimConfig._deprecated_mission_type_input_warned[] == false
        end

        @test_throws ArgumentError MissionConfiguration(mission_type="invalid")
        @test_throws ArgumentError MissionConfiguration(mission_time=0.0)
        @test_throws ArgumentError MissionConfiguration(number_of_orbits=0)
        @test_throws ArgumentError MissionConfiguration(num_steps_to_save=0)

        env_ok = EnvironmentModel(
            planet=EARTH,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=true,
            topo_degree=8,
            topo_order=8,
            wind=false
        )
        @test env_ok.EI == 120.0

        @test_throws ArgumentError EnvironmentModel(
            planet=EARTH,
            EI=-1.0,
            density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=false,
            topo_degree=8,
            topo_order=8,
            wind=false
        )

        @test_throws ArgumentError EnvironmentModel(
            planet=EARTH,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=true,
            topo_degree=-1,
            topo_order=8,
            wind=false
        )

        @test_throws ArgumentError EnvironmentModel(
            planet=EARTH,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=true,
            topo_degree=8,
            topo_order=-1,
            wind=false
        )
    end

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

        sc1 = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0)
        sc2 = make_spacecraft(ra_alt_m=550e3, rp_alt_m=500e3, ν_deg=160.0)
        args_bad_slots = build_config_multi(
            spacecraft=[sc1, sc2],
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=60.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            control_effectors=(make_base_thruster_model(thrust=1.0, Δv=1.0, start_burn_time=-1.0, stop_burn_time=-1.0),),
            control_rates=[1.0],
            keplerian=true
        )
        @test_throws ArgumentError run_case_silent(args_bad_slots)
    end
end

@testset "RHS Completeness: Mass Derivative" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=120.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )

    u = build_initial_conditions(args)
    du = copy(u)
    du.sc[1].mass = 789.0
    p = ODEParams{1}(args=args)
    spacecraft_dynamics!(du, u, p, 0.0)
    @test du.sc[1].mass == 0.0

    du_inactive = copy(u)
    du_inactive.sc[1].mass = 123.0
    p_inactive = ODEParams{1}(args=args, is_active=[false])
    spacecraft_dynamics!(du_inactive, u, p_inactive, 0.0)
    @test du_inactive.sc[1].mass == 0.0
end

@testset "Legacy Targeting Smoke" begin
    ensure_legacy_targeting_loaded!()
    sandbox = LEGACY_TARGETING_SANDBOX

    @test isdefined(sandbox, :_legacy_control_log_enabled)
    @test isdefined(sandbox, :_legacy_targeting_log_enabled)
    @test isdefined(sandbox, :_legacy_sim_targeting_log_enabled)
    @test isdefined(sandbox, :_legacy_eom_targeting_log_enabled)
    @test isdefined(sandbox, :_legacy_eom_ctrl_log_enabled)
    @test isdefined(sandbox, :_legacy_get_cnf)
    @test isdefined(sandbox, :_legacy_get_solution)
    @test isdefined(sandbox, :_legacy_control_strict_exceptions)
    @test isdefined(sandbox, :_legacy_control_exception_fallback)

    args_quiet = Dict{Symbol, Any}(:print_res => false, :verbose => false)
    args_print = Dict{Symbol, Any}(:print_res => true, :verbose => false)
    args_verbose = Dict{Symbol, Any}(:print_res => false, :verbose => true)

    @test sandbox._legacy_control_log_enabled(args_quiet) == false
    @test sandbox._legacy_control_log_enabled(args_print) == true
    @test sandbox._legacy_control_log_enabled(args_verbose) == true

    @test sandbox._legacy_targeting_log_enabled(args_quiet) == false
    @test sandbox._legacy_targeting_log_enabled(args_print) == true
    @test sandbox._legacy_targeting_log_enabled(args_verbose) == true

    @test sandbox._legacy_sim_targeting_log_enabled(args_quiet) == false
    @test sandbox._legacy_sim_targeting_log_enabled(args_print) == true
    @test sandbox._legacy_sim_targeting_log_enabled(args_verbose) == true

    @test sandbox._legacy_eom_targeting_log_enabled(args_quiet) == false
    @test sandbox._legacy_eom_targeting_log_enabled(args_print) == true
    @test sandbox._legacy_eom_targeting_log_enabled(args_verbose) == true

    @test sandbox._legacy_eom_ctrl_log_enabled(args_quiet) == false
    @test sandbox._legacy_eom_ctrl_log_enabled(args_print) == true
    @test sandbox._legacy_eom_ctrl_log_enabled(args_verbose) == true

    debug_key = "SPACEAGORA_DEBUG_LEGACY_CONTROL"
    old_debug = get(ENV, debug_key, nothing)
    try
        ENV[debug_key] = "1"
        @test sandbox._legacy_control_log_enabled(args_quiet) == true
        @test sandbox._legacy_targeting_log_enabled(args_quiet) == true
        @test sandbox._legacy_sim_targeting_log_enabled(args_quiet) == true
        @test sandbox._legacy_eom_targeting_log_enabled(args_quiet) == true
        @test sandbox._legacy_eom_ctrl_log_enabled(args_quiet) == true
    finally
        if old_debug === nothing
            if haskey(ENV, debug_key)
                delete!(ENV, debug_key)
            end
        else
            ENV[debug_key] = old_debug
        end
    end

    m_dummy = (aerodynamics=(α=1.234,),)
    @test sandbox.no_control(nothing, m_dummy) == 1.234

    heat_rate = sandbox.heat_rate_calc(1.0, 1e-6, 250.0, 250.0, 287.0, 1.4, 3.0, 0.3)
    @test isfinite(heat_rate)
    @test heat_rate >= 0.0

    local_cnf = (α=0.314, α_past=0.25)
    @test sandbox._legacy_get_cnf(Dict{Symbol, Any}(:cnf => local_cnf)).α == local_cnf.α
    local_solution = (orientation=(time=[1.0],),)
    @test sandbox._legacy_get_solution(Dict{Symbol, Any}(:solution => local_solution)).orientation.time[1] == 1.0
    @test sandbox._legacy_control_strict_exceptions(args_quiet) == false
    @test sandbox._legacy_control_strict_exceptions(Dict{Symbol, Any}(:strict_legacy_control_exceptions => true)) == true
    withenv("SPACEAGORA_STRICT_LEGACY_CONTROL_EXCEPTIONS" => "1") do
        @test sandbox._legacy_control_strict_exceptions(args_quiet) == true
    end
    err = ErrorException("legacy-control-fallback")
    @test sandbox._legacy_control_exception_fallback(args_quiet, "unit-test", err, stacktrace(), 7.0) == 7.0
    withenv("SPACEAGORA_DEBUG_LEGACY_CONTROL" => "1") do
        @test_logs (:warn, r"Legacy control fallback in unit-test") sandbox._legacy_control_exception_fallback(args_quiet, "unit-test", err, stacktrace(), 8.0)
    end
    withenv("SPACEAGORA_STRICT_LEGACY_CONTROL_EXCEPTIONS" => "1") do
        @test_throws ErrorException sandbox._legacy_control_exception_fallback(args_quiet, "unit-test", err, stacktrace(), 9.0)
    end
    @test sandbox.control_solarpanels_heatrate(nothing, nothing, Dict{Symbol, Any}(), [0], nothing; cnf=local_cnf) == local_cnf.α

    # Exercise legacy control branch behavior with deterministic local stubs.
    Core.eval(sandbox, :(module config
        get_spacecraft_reference_area(body) = 1.0
    end))
    Core.eval(sandbox, :(aerodynamic_coefficient_fM(ang, body, T_p, S, aero, MonteCarlo=false) = (0.1, 0.5 + 0.2*sin(ang))))

    m_struct = (aerodynamics=(α=1.2,), body=(dummy=1,))
    α_drag_max = sandbox.control_struct_load(nothing, m_struct, Dict{Symbol, Any}(:max_dyn_press => 100.0), 3.0, 250.0, 1.0, false)
    α_drag_min = sandbox.control_struct_load(nothing, m_struct, Dict{Symbol, Any}(:max_dyn_press => 0.0), 3.0, 250.0, 1.0, false)
    α_drag_mid = sandbox.control_struct_load(nothing, m_struct, Dict{Symbol, Any}(:max_dyn_press => 0.8), 3.0, 250.0, 1.0, false)
    @test α_drag_max == 1.2
    @test α_drag_min == 0.0001
    @test 0.0001 < α_drag_mid < 1.2

    m_heat_hi = (aerodynamics=(thermal_accomodation_factor=1.0, α=1.2, heat_rate_limit=1e9), planet=(R=287.0, γ=1.4))
    m_heat_lo = (aerodynamics=(thermal_accomodation_factor=1.0, α=1.2, heat_rate_limit=-1.0), planet=(R=287.0, γ=1.4))
    args_heat = Dict{Symbol, Any}()
    α_heat_hi = sandbox.control_solarpanels_heatrate(nothing, m_heat_hi, args_heat, [1], [250.0, 1e-6, 3.0], 0.0, 0.0, 0.0; cnf=local_cnf)
    α_heat_lo = sandbox.control_solarpanels_heatrate(nothing, m_heat_lo, args_heat, [1], [250.0, 1e-6, 3.0], 0.0, 0.0, 0.0; cnf=local_cnf)
    @test α_heat_hi == 1.2
    @test α_heat_lo == 0.0001

    heat_rate_max = sandbox.heat_rate_calc(1.0, 1e-6, 250.0, 250.0, 287.0, 1.4, 3.0, 1.2)
    heat_rate_min = sandbox.heat_rate_calc(1.0, 1e-6, 250.0, 250.0, 287.0, 1.4, 3.0, 0.0001)
    heat_limit_mid = (heat_rate_max + heat_rate_min) / 2 + 1e-5
    m_heat_mid = (aerodynamics=(thermal_accomodation_factor=1.0, α=1.2, heat_rate_limit=heat_limit_mid), planet=(R=287.0, γ=1.4))
    α_heat_mid = sandbox.control_solarpanels_heatrate(nothing, m_heat_mid, args_heat, [1], [250.0, 1e-6, 3.0], 0.0, 0.0, 0.0; cnf=local_cnf)
    @test 0.0001 < α_heat_mid < 1.2
end

@testset "Legacy CNF Threading Guard" begin
    function strip_comments(src::String)
        # Remove block comments first, then trim line comments (including trailing inline comments).
        no_block = replace(src, r"#=.*?=#"s => "")
        no_line = map(line -> first(split(line, '#'; limit=2)), split(no_block, '\n'; keepempty=true))
        return join(no_line, "\n")
    end

    scoped_files = [
        "src/control/Control.jl",
        "src/control/targeting_control/targeting.jl",
        "src/control/heatload_control/Time_switch_calcs.jl",
        "src/control/heatload_control/Second_tsw_calcs.jl",
        "src/control/heatload_control/Security_mode.jl",
        "src/control/targeting_control/sim_targeting.jl",
        "src/control/utils/Eoms.jl",
        "src/control/utils/Eom_ctrl.jl",
        "src/control/utils/Propulsive_maneuvers.jl",
        "src/control/Eoms.jl",
        "src/control/Eom_ctrl.jl",
        "src/control/targeting_control/Eom_targeting.jl"
    ]
    for relpath in scoped_files
        src = strip_comments(read(joinpath(REPO_ROOT, relpath), String))
        @test !occursin("config.cnf", src)
        @test !occursin("config.solution", src)
    end
end

@testset "Legacy Config Parsing Cleanup" begin
    ensure_legacy_config_loaded!()
    sandbox = LEGACY_CONFIG_SANDBOX

    @test isdefined(sandbox, :_legacy_parse_planet_id)
    @test sandbox._legacy_parse_planet_id("earth") == 0
    @test sandbox._legacy_parse_planet_id(:Mars) == 1
    @test sandbox._legacy_parse_planet_id(7) == 7
    @test_throws ArgumentError sandbox._legacy_parse_planet_id("pluto")
    @test_throws ArgumentError sandbox.planet_data("earth")

    @test isdefined(sandbox, :_legacy_mission_kind)
    @test sandbox._legacy_mission_kind("Time") == :time
    @test sandbox._legacy_mission_kind("Drag Passage") == :drag_passage
    @test sandbox._legacy_mission_kind(MissionTime) == :time
    @test sandbox._legacy_mission_kind(MissionOrbits) == :orbits

    mission = Dict{Symbol, Any}(
        :Planet => "Venus",
        :Gravity_Model => "Inverse Squared and J2 Effect",
        :Density_Model => "No-Density",
        :Wind => true,
        :Aerodynamic_Model => "Mach-dependent",
        :Shape => "Spacecraft",
        :Control => 2,
        :Firings => "Drag Passage Firing",
        :Thermal_Model => "Shaaf and Chambre",
        :Monte_Carlo => 0
    )
    ip = withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
        sandbox.mission_def(mission)
    end
    @test ip.gm isa LegacyGravityModelCode
    @test ip.dm isa LegacyDensityModelCode
    @test ip.am isa LegacyAerodynamicModelCode
    @test ip.tc isa LegacyThrustControlCode
    @test ip.tm isa LegacyThermalModelCode
    @test ip.M.planet == 2
    @test ip.gm == 2
    @test ip.dm == 2
    @test ip.wm == 1
    @test ip.am == 1
    @test ip.cm == 2
    @test ip.tc == 2
    @test ip.tm == 2
    @test ip.mc == 0

    mission_default = deepcopy(mission)
    mission_default[:Planet] = "UnknownPlanet"
    ip_default = withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
        sandbox.mission_def(mission_default)
    end
    @test ip_default.M.planet == 1

    @test isdefined(sandbox, :_deprecated_mission_dict_input_warned)
    empty!(sandbox._deprecated_mission_dict_input_warned[])
    withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "1") do
        @test_logs (:warn, r"mission_def: legacy string/symbol value for :Gravity_Model") sandbox.mission_def(Dict{Symbol, Any}(:Gravity_Model => "Constant"))
        @test_logs sandbox.mission_def(Dict{Symbol, Any}(:Gravity_Model => "Inverse Squared"))
    end
    empty!(sandbox._deprecated_mission_dict_input_warned[])
    withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
        @test_logs sandbox.mission_def(Dict{Symbol, Any}(:Density_Model => "Exponential"))
    end

    ip_from_ints = InitialParameters(Mission(), 1, 3, 0, 2, 2, 0, 1, 0)
    @test ip_from_ints.gm == LegacyGravityInverseSquared
    @test ip_from_ints.dm == LegacyDensityGRAM
    @test ip_from_ints.am == LegacyAeroNoBallisticAxial
    @test ip_from_ints.tm == LegacyThermalMaxwellian
    @test ip_from_ints.tc == LegacyThrustAerobrakingManeuver
    @test ip_from_ints.gm == 1
    @test ip_from_ints.dm == 3
    @test_throws ArgumentError InitialParameters(Mission(), 99, 1, 0, 0, 1, 0, 0, 0)

    @test _legacy_run_planet_id("earth") == 0
    @test _legacy_run_planet_id(:Mars) == 1
    @test _legacy_run_planet_id(7) == 7
    @test_throws ArgumentError _legacy_run_planet_id("pluto")

    planet_earth = _legacy_typed_planet("earth", Dict{Symbol, Any}())
    @test getfield(planet_earth, :name) == "Earth"
    planet_jupiter = _legacy_typed_planet("jupiter", Dict{Symbol, Any}())
    @test isapprox(planet_jupiter.Rp_e, 7.1492e7; atol=0.0, rtol=0.0)
    @test isapprox(planet_jupiter.μ, 1.26686534e17; atol=0.0, rtol=0.0)

    args_time = Dict{Symbol, Any}(
        :type_of_mission => "time",
        :drag_passage => 1,
        :number_of_orbits => 3,
        :body_shape => "Spacecraft",
        :aerodynamic_model => "No-Ballistic flight with axial coefficient",
        :thermal_model => "Convective and Radiative",
        :control_mode => 1,
        :thrust_control => "Aerobraking Maneuver",
        :thrust => 0.0,
        :delta_v => 5.0,
        :Odyssey_sim => false
    )
    out_time = redirect_stdout(devnull) do
        sandbox.def_miss(deepcopy(args_time))
    end
    @test out_time[:drag_passage] == 0
    @test out_time[:aerodynamic_model] == "Mach-dependent"
    @test out_time[:thermal_model] == "Maxwellian Heat Transfer"
    @test out_time[:thrust] == 0.1

    args_entry = deepcopy(args_time)
    args_entry[:type_of_mission] = "Entry"
    out_entry = redirect_stdout(devnull) do
        sandbox.def_miss(args_entry)
    end
    @test out_entry[:drag_passage] == 1
    @test out_entry[:number_of_orbits] == 1

    args_typed = deepcopy(args_time)
    args_typed[:type_of_mission] = MissionTime
    out_typed = redirect_stdout(devnull) do
        sandbox.def_miss(args_typed)
    end
    @test out_typed[:drag_passage] == 0

    args_campaign = deepcopy(args_time)
    args_campaign[:type_of_mission] = "Aerobraking Campaign"
    args_campaign[:number_of_orbits] = 12
    out_campaign = redirect_stdout(devnull) do
        sandbox.def_miss(args_campaign)
    end
    @test out_campaign[:drag_passage] == 0
    @test out_campaign[:number_of_orbits] == 1000

    args_blunted = deepcopy(args_time)
    args_blunted[:body_shape] = "Blunted Cone"
    args_blunted[:aerodynamic_model] = "Mach-dependent"
    args_blunted[:thermal_model] = "Maxwellian Heat Transfer"
    args_blunted[:control_mode] = 2
    args_blunted[:thrust_control] = "Aerobraking Maneuver"
    args_blunted[:thrust] = 0.0
    args_blunted[:delta_v] = 7.5
    out_blunted = redirect_stdout(devnull) do
        sandbox.def_miss(args_blunted)
    end
    @test out_blunted[:aerodynamic_model] == "No-Ballistic flight with axial coefficient"
    @test out_blunted[:thermal_model] == "Convective and Radiative"
    @test out_blunted[:control_mode] == 0
    @test out_blunted[:thrust] == 0.1
    @test out_blunted[:delta_v] == 7.5

    args_drag_passage_maneuver = deepcopy(args_time)
    args_drag_passage_maneuver[:type_of_mission] = "Drag Passage"
    args_drag_passage_maneuver[:thrust_control] = "Aerobraking Maneuver"
    out_drag_passage_maneuver = redirect_stdout(devnull) do
        sandbox.def_miss(args_drag_passage_maneuver)
    end
    @test out_drag_passage_maneuver[:drag_passage] == 1
    @test out_drag_passage_maneuver[:thrust_control] == "None"

    args_odyssey = deepcopy(args_time)
    args_odyssey[:Odyssey_sim] = true
    args_odyssey[:gravity_model] = "Inverse Squared"
    out_odyssey = redirect_stdout(devnull) do
        sandbox.def_miss(args_odyssey)
    end
    @test out_odyssey[:type_of_mission] == "Aerobraking Campaign"
    @test out_odyssey[:number_of_orbits] == 350
    @test out_odyssey[:planet] == 1
    @test out_odyssey[:control_mode] == 0
    @test out_odyssey[:thrust_control] == "Aerobraking Maneuver"
    @test out_odyssey[:drag_passage] == 0
    @test out_odyssey[:hp_initial_a] == 87000
    @test out_odyssey[:final_apoapsis] == 3390.0e3 + 503e3
    @test haskey(out_odyssey, :inital_condition_type)
    @test out_odyssey[:inital_condition_type] == 0

    # When a legacy Planet constructor is provided, cover all legacy planet branches.
    Core.eval(sandbox, :(Planet(args...) = (Rp_e=args[1], μ=args[16], name=args[27], topography_function=args[end])))
    legacy_planet_expected = Dict(
        0 => (6.3781e6, 3.98600436000000e14, "earth"),
        1 => (3.3962e6, 4.28283140000000e13, "mars"),
        2 => (6.0518e6, 3.24858592e14, "venus"),
        3 => (6.9634e8, 1.3271244002331e20, "sun"),
        4 => (1.7381e6, 4.9028005821478e12, "moon"),
        5 => (7.1492e7, 1.26686534e17, "jupiter"),
        6 => (6.0268e7, 3.7931187e16, "saturn"),
        7 => (2.575e6, 8.981e12, "titan")
    )
    for pid in 0:7
        pdat = sandbox.planet_data(pid)
        rp_ref, mu_ref, name_ref = legacy_planet_expected[pid]
        @test isapprox(pdat.Rp_e, rp_ref; atol=0.0, rtol=0.0)
        @test isapprox(pdat.μ, mu_ref; atol=0.0, rtol=0.0)
        @test lowercase(String(pdat.name)) == name_ref
    end
end

@testset "Legacy Integrator/Event/Output Smoke" begin
    sandbox = LEGACY_RUNTIME_SANDBOX

    Core.eval(sandbox, :(include(joinpath(Main.REPO_ROOT, "src", "integrator", "Events.jl"))))
    Core.eval(sandbox, :(include(joinpath(Main.REPO_ROOT, "src", "integrator", "Integrators.jl"))))
    Core.eval(sandbox, :(include(joinpath(Main.REPO_ROOT, "src", "utils", "Save_results.jl"))))

    @test sandbox.event(0, 0) == true
    @test redirect_stdout(devnull) do
        sandbox.event(1, 0)
    end == false
    @test redirect_stdout(devnull) do
        sandbox.event(0, 1)
    end == false

    m = (planet=(Rp_e=6.3781e6, μ=3.986004418e14),)
    solution_state = (t_events=[Any[]],)
    y_impact = [m.planet.Rp_e + 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0]
    breaker_impact, sol_after_impact = redirect_stdout(devnull) do
        sandbox.impact(0.0, y_impact, m, solution_state, Dict{Symbol, Any}())
    end
    @test breaker_impact == true
    @test sol_after_impact.t_events[end] == ["true"]

    y_250_in = [m.planet.Rp_e + 240e3, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0]
    y_250_out = [m.planet.Rp_e + 260e3, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0]
    y_120_in = [m.planet.Rp_e + 110e3, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0]
    y_120_out = [m.planet.Rp_e + 130e3, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0]
    @test sandbox.eventsecondstep(1.0, y_250_out, 0.0, y_250_in, m, Dict{Symbol, Any}()) == true
    @test sandbox.heat_check(1.0, y_120_out, 0.0, y_120_in, m, Dict{Symbol, Any}()) == true
    @test sandbox.out_drag_passage(1.0, [m.planet.Rp_e + 11e3, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0], 0.0, [m.planet.Rp_e + 9e3, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0], m, Dict{Symbol, Any}(:AE => 10e3)) == true

    runtime_ok = sandbox._legacy_get_save_results_runtime_state(Dict{Symbol, Any}(:cnf => 1, :solution => 2, :model => 3))
    @test runtime_ok.cnf == 1
    @test runtime_ok.solution == 2
    @test runtime_ok.model == 3
    @test_throws ArgumentError sandbox._legacy_get_save_results_runtime_state(Dict{Symbol, Any}(:cnf => 1))
end

@testset "Legacy Save_csv Direct Smoke" begin
    sandbox = LEGACY_RUNTIME_SANDBOX
    if !isdefined(sandbox, :save_csv)
        Core.eval(sandbox, :(include(joinpath(Main.REPO_ROOT, "src", "utils", "Save_csv.jl"))))
    end

    body = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, i_deg=35.0, ω_deg=20.0, Ω_deg=15.0, ν_deg=120.0)
    body.n_reaction_wheels = 2
    body.n_thrusters = 2
    model = Model(body=body)
    cnf = with_logger(Logging.NullLogger()) do
        Cnf()
    end

    mktempdir() do tmp
        csv_path = joinpath(tmp, "legacy_save_csv.csv")
        arrow_path = joinpath(tmp, "legacy_save_csv.arrow")

        args_no_cf = Dict{Symbol, Any}(:save_csv => true, :closed_form => 0, :print_res => false)
        solution_no_cf = seed_solution_for_save_csv!(
            Solution();
            n_bodies=length(model.body.links),
            n_reaction_wheels=model.body.n_reaction_wheels,
            n_thrusters=model.body.n_thrusters,
            base=1.0,
            closed_form=false
        )
        redirect_stdout(devnull) do
            sandbox.save_csv(csv_path, args_no_cf, arrow_path, (cnf, model, solution_no_cf))
        end

        @test isfile(csv_path)
        @test isfile(arrow_path)
        @test filesize(csv_path) > 0
        @test filesize(arrow_path) > 0
        df_first = CSV.read(csv_path, DataFrame)
        @test nrow(df_first) == 1
        @test "link_1_aoa" in names(df_first)
        @test "link_2_aoa" in names(df_first)
        @test "rw_h_1" in names(df_first)
        @test "rw_h_2" in names(df_first)
        @test "thruster_force_1" in names(df_first)
        @test "thruster_force_2" in names(df_first)

        arrow_size_before = filesize(arrow_path)
        args_with_cf = Dict{Symbol, Any}(:save_csv => true, :closed_form => 1, :print_res => false)
        solution_with_cf = seed_solution_for_save_csv!(
            Solution();
            n_bodies=length(model.body.links),
            n_reaction_wheels=model.body.n_reaction_wheels,
            n_thrusters=model.body.n_thrusters,
            base=2.0,
            closed_form=true
        )
        redirect_stdout(devnull) do
            sandbox.save_csv(csv_path, args_with_cf, arrow_path, (cnf, model, solution_with_cf))
        end

        @test filesize(arrow_path) > arrow_size_before
        df_second = CSV.read(csv_path, DataFrame)
        @test nrow(df_second) == 2
        @test isapprox(Float64(df_second.t_cf[end]), 12.0; atol=0.0, rtol=0.0)
        @test isapprox(Float64(df_second.h_cf[end]), 22.0; atol=0.0, rtol=0.0)
        @test isapprox(Float64(df_second.gamma_cf[end]), 32.0; atol=0.0, rtol=0.0)
        @test isapprox(Float64(df_second.v_cf[end]), 42.0; atol=0.0, rtol=0.0)
    end
end

@testset "Guidance Thruster Campaign Helpers" begin
    ensure_guidance_sandbox_loaded!()
    sandbox = GUIDANCE_SANDBOX

    guidance_model = sandbox.AerobrakingCampaignPropulsiveManeuverGuidanceModel(
        maneuver_orbit_number=[2, 3],
        maneuver_Δv=[-5.0, 20.0]
    )

    thruster = BaseThrusterModel(
        thrust=[500.0],
        direction=[0.0],
        Δv=[0.0],
        start_burn_time=[-1.0],
        stop_burn_time=[-1.0],
        Isp=[300.0]
    )
    args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(thruster,),
        control_rates=[1.0],
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p = ODEParams{1}(args=args)
    u = ComponentVector(pos=[0.0, 0.0, 0.0], vel=[0.0, 0.0, 0.0], mass=1.0, heat_loads=[0.0])

    p.orbit_counter[1] = 3
    sandbox.calcGuidanceEffect!(guidance_model, u, p, 0.0, 1)
    @test thruster.Δv[1] == 20.0
    @test isapprox(thruster.direction[1], π; atol=1e-12, rtol=0.0)

    p.orbit_counter[1] = 2
    sandbox.calcGuidanceEffect!(guidance_model, u, p, 0.0, 1)
    @test thruster.Δv[1] == 5.0
    @test isapprox(thruster.direction[1], 0.0; atol=1e-12, rtol=0.0)

    p.orbit_counter[1] = 1
    sandbox.calcGuidanceEffect!(guidance_model, u, p, 0.0, 1)
    @test thruster.Δv[1] == 0.0
end

@testset "Legacy Monte Carlo Helpers Smoke" begin
    ensure_legacy_monte_carlo_loaded!()
    sandbox = LEGACY_MONTE_CARLO_SANDBOX

    cnf = with_logger(Logging.NullLogger()) do
        Cnf()
    end
    solution = Solution()

    runtime = sandbox._legacy_get_mc_runtime_state(Dict{Symbol, Any}(:cnf => cnf, :solution => solution))
    @test runtime.cnf === cnf
    @test runtime.solution === solution
    @test_throws ArgumentError sandbox._legacy_get_mc_runtime_state(Dict{Symbol, Any}(:cnf => cnf))

    mc_dict, count, args_mc = sandbox.MonteCarlo_setting(Dict{Symbol, Any}(:montecarlo => false))
    @test count == 0
    @test args_mc[:montecarlo_size] == 1
    @test args_mc[:intial_montecarlo_number] == 1
    @test all(k -> haskey(mc_dict, k), (:N_passages, :Duration, :Median_heat, :Periapsis_min, :Periapsis_max))

    args_passage = Dict{Symbol, Any}(
        :cnf => cnf,
        :solution => solution,
        :print_res => false,
        :simulation_filename => "campaign"
    )
    cnf.altitude_periapsis = [1.0]
    cnf.max_heatrate = [2.0]
    sandbox.MonteCarlo_setting_passage(3, args_passage)
    @test cnf.counter_random == 3
    @test cnf.index_MonteCarlo == 3
    @test isempty(cnf.altitude_periapsis)
    @test isempty(cnf.max_heatrate)
    @test endswith(args_passage[:simulation_filename], "_nMC=3")

    solution.orientation.number_of_passage = [7]
    solution.orientation.time = [123.0]
    solution.performance.heat_rate = [[10.0, 12.0], [1.0]]
    cnf.max_heatrate = [4.0, 16.0]
    cnf.altitude_periapsis = [80.0, 100.0]
    args_append = Dict{Symbol, Any}(
        :cnf => cnf,
        :solution => solution,
        :type_of_mission => "Time",
        :max_heat_rate => 11.0,
        :print_res => false
    )
    count_after = sandbox.MonteCarlo_append(mc_dict, args_append, 0)
    @test count_after == 1
    @test mc_dict[:N_passages][1] == 7
    @test isapprox(mc_dict[:Duration][1], 123.0; atol=0.0, rtol=0.0)
    @test isapprox(mc_dict[:Median_heat][1], 10.0; atol=0.0, rtol=0.0)
    @test mc_dict[:Periapsis_min][1] == 80.0
    @test mc_dict[:Periapsis_max][1] == 100.0

    mktempdir() do tmp
        args_save = Dict{Symbol, Any}(
            :save_results => false,
            :simulation_filename => "campaign_nMC=1",
            :directory_results => tmp * "/",
            :control_mode => 1,
            :max_heat_rate => 100.0,
            :montecarlo_size => 1
        )
        state = Dict{Symbol, Any}(:Apoapsis => 7000e3, :Periapsis => 120.0)
        @test_nowarn sandbox.MonteCarlo_save(args_save, state, mc_dict)
        @test isempty(readdir(tmp))
    end
end

@testset "Legacy Closed-Form Helper Smoke" begin
    ensure_legacy_closed_form_loaded!()
    sandbox = LEGACY_CLOSED_FORM_SANDBOX

    solution = Solution()
    sandbox.results(solution, [1.0], [2.0], [3.0], [4.0])
    @test solution.closed_form.t_cf == [1.0]
    @test solution.closed_form.h_cf == [2.0]
    @test solution.closed_form.γ_cf == [3.0]
    @test solution.closed_form.v_cf == [4.0]

    solution_blunted = Solution()
    solution_blunted.orientation.time = [0.0, 1.0, 2.0]
    cnf_blunted = with_logger(Logging.NullLogger()) do
        Cnf()
    end
    params = (cnf_blunted, nothing, solution_blunted)
    mission_stub = (initial_condition=(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),)
    args_blunted = Dict{Symbol, Any}(:body_shape => "Blunted Cone")

    t_cf, h_cf, γ_cf, v_cf = sandbox.closed_form(args_blunted, mission_stub, params)
    @test t_cf == zeros(3)
    @test h_cf == zeros(3)
    @test γ_cf == zeros(3)
    @test v_cf == zeros(3)
    @test solution_blunted.closed_form.t_cf == zeros(3)
    @test solution_blunted.closed_form.h_cf == zeros(3)
    @test solution_blunted.closed_form.γ_cf == zeros(3)
    @test solution_blunted.closed_form.v_cf == zeros(3)
end

@testset "Typed Planet Constructors + Topography Workspace" begin
    mars = Mars("", SPICE_PATH)
    venus = Venus("", SPICE_PATH)
    titan = Titan("", SPICE_PATH)

    @test mars.name == "Mars"
    @test venus.name == "Venus"
    @test titan.name == "Titan"
    @test mars.μ > 0.0
    @test venus.μ > 0.0
    @test titan.μ > 0.0

    mktempdir() do tmp
        topo_file = joinpath(tmp, "topo_harmonics.csv")
        write(topo_file, "degree,order,C,S\n0,0,1.0,0.0\n1,0,0.1,0.0\n1,1,0.05,0.02\n")
        earth = Earth("", SPICE_PATH)
        SimulationModel.Planets.TopographyHarmonicsWorkspace!(topo_file, earth)
        @test size(earth.topography_workspace.Clm) == (2, 2)
        @test size(earth.topography_workspace.Slm) == (2, 2)
        @test isapprox(earth.topography_workspace.Clm[2, 2], 0.05; atol=0.0, rtol=0.0)
        @test isapprox(earth.topography_workspace.Slm[2, 2], 0.02; atol=0.0, rtol=0.0)
    end
end

@testset "Run Simulation State Isolation" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=45.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    l_pi_before = args.environment_model.planet.L_PI
    @test_nowarn run_simulation(args)
    @test args.environment_model.planet.L_PI == l_pi_before

    args_no_isolation = deepcopy(args)
    @test_nowarn run_simulation(args_no_isolation; isolate_state=false)

    if Threads.nthreads() >= 2
        args_parallel = deepcopy(args)
        t1 = Threads.@spawn run_simulation(args_parallel)
        t2 = Threads.@spawn run_simulation(args_parallel)
        @test_nowarn fetch(t1)
        @test_nowarn fetch(t2)
    end
end

@testset "Legacy Global State Guard (Remaining Modules)" begin
    function strip_comments(src::String)
        no_block = replace(src, r"#=.*?=#"s => "")
        no_line = map(line -> first(split(line, '#'; limit=2)), split(no_block, '\n'; keepempty=true))
        return join(no_line, "\n")
    end

    remaining_files = [
        "src/utils/Save_results.jl",
        "src/utils/MonteCarlo_set.jl",
        "src/physical_models/Propulsive_maneuvers.jl"
    ]

    for relpath in remaining_files
        src = strip_comments(read(joinpath(REPO_ROOT, relpath), String))
        @test !occursin("config.", src)
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

@testset "Legacy Wrapper Parity + Robustness" begin
    sc = make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0)
    settings = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=settings,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    df_direct = run_case_silent(args)
    df_wrapper = run_case_via_run_analysis(args)

    @test nrow(df_direct) == nrow(df_wrapper)
    @test nrow(df_wrapper) > 10
    sample_idxs = round.(Int, range(1, nrow(df_direct), length=8))
    for idx in sample_idxs
        t = Float64(df_direct.time[idx])
        p_direct = SVector{3, Float64}(Float64(df_direct.sc1_pos_1[idx]), Float64(df_direct.sc1_pos_2[idx]), Float64(df_direct.sc1_pos_3[idx]))
        v_direct = SVector{3, Float64}(Float64(df_direct.sc1_vel_1[idx]), Float64(df_direct.sc1_vel_2[idx]), Float64(df_direct.sc1_vel_3[idx]))
        p_wrapper = SVector{3, Float64}(
            interp_linear(df_wrapper.time, df_wrapper.sc1_pos_1, t),
            interp_linear(df_wrapper.time, df_wrapper.sc1_pos_2, t),
            interp_linear(df_wrapper.time, df_wrapper.sc1_pos_3, t)
        )
        v_wrapper = SVector{3, Float64}(
            interp_linear(df_wrapper.time, df_wrapper.sc1_vel_1, t),
            interp_linear(df_wrapper.time, df_wrapper.sc1_vel_2, t),
            interp_linear(df_wrapper.time, df_wrapper.sc1_vel_3, t)
        )
        @test norm(p_direct - p_wrapper) < 0.1
        @test norm(v_direct - v_wrapper) < 1e-4
    end

    settings_no_csv = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=false,
        normalize=false
    )
    args_no_csv = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=settings_no_csv,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    mktempdir() do tmp
        cd(tmp) do
            redirect_stdout(devnull) do
                run_analysis(args_no_csv)
            end
            @test isfile("simulation_results.csv")
            @test isfile(joinpath("output", "results.feather"))
        end
    end

    settings_no_results = SimulationSettings(
        results=false,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=false,
        normalize=false
    )
    args_no_results = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=settings_no_results,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    _ = run_case_via_run_analysis(args_no_results; expect_results_csv=false)
end

@testset "Legacy Entry Dispatch Parity" begin
    sc = make_spacecraft(ra_alt_m=540e3, rp_alt_m=430e3, ν_deg=168.0)
    settings = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=400.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=settings,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    function assert_df_parity(df_ref::DataFrame, df_cmp::DataFrame)
        @test nrow(df_ref) == nrow(df_cmp)
        @test nrow(df_ref) > 10
        sample_idxs = round.(Int, range(1, nrow(df_ref), length=8))
        for idx in sample_idxs
            t = Float64(df_ref.time[idx])
            p_ref = SVector{3, Float64}(Float64(df_ref.sc1_pos_1[idx]), Float64(df_ref.sc1_pos_2[idx]), Float64(df_ref.sc1_pos_3[idx]))
            v_ref = SVector{3, Float64}(Float64(df_ref.sc1_vel_1[idx]), Float64(df_ref.sc1_vel_2[idx]), Float64(df_ref.sc1_vel_3[idx]))
            p_cmp = SVector{3, Float64}(
                interp_linear(df_cmp.time, df_cmp.sc1_pos_1, t),
                interp_linear(df_cmp.time, df_cmp.sc1_pos_2, t),
                interp_linear(df_cmp.time, df_cmp.sc1_pos_3, t)
            )
            v_cmp = SVector{3, Float64}(
                interp_linear(df_cmp.time, df_cmp.sc1_vel_1, t),
                interp_linear(df_cmp.time, df_cmp.sc1_vel_2, t),
                interp_linear(df_cmp.time, df_cmp.sc1_vel_3, t)
            )
            @test norm(p_ref - p_cmp) < 0.1
            @test norm(v_ref - v_cmp) < 1e-4
        end
    end

    df_direct = run_case_silent(args; isolate_state=false)
    df_analysis = run_case_via_run_analysis(args; isolate_state=false)
    df_campaign = run_case_via_campaign(args; isolate_state=false)
    df_campaign_with_state = run_case_via_campaign(args; isolate_state=false, state=Dict(:legacy => true))
    df_run_oe = run_case_via_run_orbitalelements(args; isolate_state=false)
    df_run_vg = run_case_via_run_vgamma(args; isolate_state=false)

    assert_df_parity(df_direct, df_analysis)
    assert_df_parity(df_direct, df_campaign)
    assert_df_parity(df_direct, df_campaign_with_state)
    assert_df_parity(df_direct, df_run_oe)
    assert_df_parity(df_direct, df_run_vg)

    @test_throws ArgumentError aerobraking_campaign(Dict{Symbol, Any}(:foo => :bar))
end

@testset "Units/Normalization Consistency Audit" begin
    function strip_comments(src::String)
        no_block = replace(src, r"#=.*?=#"s => "")
        no_line = map(line -> first(split(line, '#'; limit=2)), split(no_block, '\n'; keepempty=true))
        return join(no_line, "\n")
    end

    run_src = strip_comments(read(joinpath(REPO_ROOT, "src", "simulation", "run_simulation.jl"), String))
    complete_src = strip_comments(read(joinpath(REPO_ROOT, "src", "simulation", "Complete_passage.jl"), String))

    # Typed path should stay SI-native; legacy path still carries DU/TU/MU normalization.
    @test !occursin("cnf.DU", run_src)
    @test !occursin("cnf.TU", run_src)
    @test !occursin("cnf.MU", run_src)
    @test occursin("cnf.DU", complete_src)
    @test occursin("cnf.TU", complete_src)
    @test occursin("cnf.MU", complete_src)

    sc = make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0)
    settings_norm_true = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        normalize=true
    )
    settings_norm_false = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        normalize=false
    )

    args_norm_true = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=700.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=settings_norm_true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    args_norm_false = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=700.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=settings_norm_false,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    # Build-state sanity: typed initial conditions are SI-scale, not O(1) normalized values.
    u0 = build_initial_conditions(args_norm_true)
    @test norm(SVector{3, Float64}(u0.sc[1].pos)) > 1e6
    @test norm(SVector{3, Float64}(u0.sc[1].vel)) > 1e3
    @test u0.sc[1].mass > 1.0

    df_norm_true = withenv("SPACEAGORA_WARN_NORMALIZE" => "0") do
        run_case_silent(args_norm_true)
    end
    df_norm_false = run_case_silent(args_norm_false)

    @test nrow(df_norm_true) == nrow(df_norm_false)
    @test nrow(df_norm_true) > 10
    sample_idxs = round.(Int, range(1, nrow(df_norm_true), length=8))
    for idx in sample_idxs
        t = Float64(df_norm_true.time[idx])
        p_true = SVector{3, Float64}(Float64(df_norm_true.sc1_pos_1[idx]), Float64(df_norm_true.sc1_pos_2[idx]), Float64(df_norm_true.sc1_pos_3[idx]))
        v_true = SVector{3, Float64}(Float64(df_norm_true.sc1_vel_1[idx]), Float64(df_norm_true.sc1_vel_2[idx]), Float64(df_norm_true.sc1_vel_3[idx]))
        p_false = SVector{3, Float64}(
            interp_linear(df_norm_false.time, df_norm_false.sc1_pos_1, t),
            interp_linear(df_norm_false.time, df_norm_false.sc1_pos_2, t),
            interp_linear(df_norm_false.time, df_norm_false.sc1_pos_3, t)
        )
        v_false = SVector{3, Float64}(
            interp_linear(df_norm_false.time, df_norm_false.sc1_vel_1, t),
            interp_linear(df_norm_false.time, df_norm_false.sc1_vel_2, t),
            interp_linear(df_norm_false.time, df_norm_false.sc1_vel_3, t)
        )
        @test norm(p_true - p_false) < 0.1
        @test norm(v_true - v_false) < 1e-4
    end
end

@testset "Normalize Flag Runtime Warning" begin
    args_warn = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=true),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    withenv("SPACEAGORA_WARN_NORMALIZE" => "1") do
        _normalize_warning_emitted[] = false
        @test_logs (:warn, r"normalize=true is legacy-only") run_simulation(args_warn)
        @test _normalize_warning_emitted[] == true
    end

    withenv("SPACEAGORA_WARN_NORMALIZE" => "0") do
        _normalize_warning_emitted[] = false
        @test_logs run_simulation(args_warn)
        @test _normalize_warning_emitted[] == false
    end
end

@testset "Verbose Gating for Callback/Runtime Logs" begin
    settings_quiet = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false
    )
    settings_verbose = SimulationSettings(
        results=true,
        verbose=true,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false
    )

    args_orbit_quiet = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=400e3, ν_deg=175.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=400.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=settings_quiet,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    _, out_orbit_quiet = run_case_capture_stdout(args_orbit_quiet)
    @test !occursin("Initial conditions:", out_orbit_quiet)
    @test !occursin("Orbit ", out_orbit_quiet)

    args_orbit_verbose = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=400e3, ν_deg=175.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=400.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=settings_verbose,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    _, out_orbit_verbose = run_case_capture_stdout(args_orbit_verbose)
    @test occursin("Initial conditions:", out_orbit_verbose)
    @test occursin("Orbit ", out_orbit_verbose)

    args_drag_quiet = build_config(
        spacecraft=make_spacecraft(ra_alt_m=220e3, rp_alt_m=100e3, ν_deg=180.0),
        density_model=ConstantDensityModel(1e-6, 240.0),
        orientation_sim=false,
        mission_time=700.0,
        EI_km=140.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        keplerian=false,
        simulation_settings=settings_quiet,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_orbit=5.0,
            dt_max_atmosphere=0.2
        )
    )
    _, out_drag_quiet = run_case_capture_stdout(args_drag_quiet)
    @test !occursin("Switching to atmosphere integration", out_drag_quiet)
    @test !occursin("Impact detected", out_drag_quiet)

    args_drag_verbose = build_config(
        spacecraft=make_spacecraft(ra_alt_m=220e3, rp_alt_m=100e3, ν_deg=180.0),
        density_model=ConstantDensityModel(1e-6, 240.0),
        orientation_sim=false,
        mission_time=700.0,
        EI_km=140.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        keplerian=false,
        simulation_settings=settings_verbose,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_orbit=5.0,
            dt_max_atmosphere=0.2
        )
    )
    _, out_drag_verbose = run_case_capture_stdout(args_drag_verbose)
    @test occursin("Switching to atmosphere integration", out_drag_verbose)
    @test occursin("Impact detected", out_drag_verbose)
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

    @testset "calcControlMassFlowRate" begin
        model = make_base_thruster_model(
            thrust=2.0,
            direction=0.0,
            start_burn_time=50.0,
            stop_burn_time=150.0,
            Isp=300.0
        )

        mdot_pre = calcControlMassFlowRate(model, u, p, 1, 49.9)
        @test mdot_pre == 0.0

        mdot_on = calcControlMassFlowRate(model, u, p, 1, 100.0)
        expected_mdot = -2.0 / (300.0 * 9.80665)
        @test isapprox(mdot_on, expected_mdot; atol=1e-14, rtol=0.0)

        mdot_post = calcControlMassFlowRate(model, u, p, 1, 150.1)
        @test mdot_post == 0.0

        model_bad_isp = make_base_thruster_model(
            thrust=2.0,
            direction=0.0,
            start_burn_time=50.0,
            stop_burn_time=150.0,
            Isp=0.0
        )
        @test calcControlMassFlowRate(model_bad_isp, u, p, 1, 100.0) == 0.0

        u_zero = ComponentVector(pos=[0.0, 0.0, 0.0], vel=[0.0, 0.0, 0.0], mass=1.0, heat_loads=[0.0])
        @test calcControlMassFlowRate(model, u_zero, p, 1, 100.0) == 0.0
        @test calcControlMassFlowRate(model, u, p, 2, 100.0) == 0.0

        model_subtyped = TimedTangentialThrusterModel(1.0, +1.0, 10.0, 20.0)
        @test calcControlMassFlowRate(model_subtyped, u, p, 1, 15.0) == 0.0

        struct UntypedControlEffector end
        @test calcControlMassFlowRate(UntypedControlEffector(), u, p, 1, 100.0) == 0.0
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

        sc_circular = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=330.0)
        args_circular = build_config(
            spacecraft=sc_circular,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=600.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true
        )
        p_circular = ODEParams{1}(args=args_circular)
        state_circular = build_initial_conditions(args_circular)
        model_circular = make_base_thruster_model(thrust=2.0, Δv=20.0, start_burn_time=-1.0, stop_burn_time=-1.0)
        calcControlEffect!(model_circular, state_circular, p_circular, 100.0, 1)
        @test model_circular.start_burn_time[1] != -1.0
        @test model_circular.stop_burn_time[1] != -1.0

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

        helper = SimulationModel.ControlEffectors._control_effector_exception_fallback
        p_dummy = (shared_buffers=(debug_control=Ref(false),),)
        err_eff = ErrorException("control-effector-fallback")
        @test helper(p_dummy, 1, err_eff, stacktrace()) === nothing
        withenv("SPACEAGORA_DEBUG_CONTROL" => "1", "SPACEAGORA_STRICT_CONTROL_EXCEPTIONS" => "0") do
            @test_logs (:warn, r"orbital-element conversion failed") helper(p_dummy, 1, err_eff, stacktrace())
        end
        withenv("SPACEAGORA_STRICT_CONTROL_EXCEPTIONS" => "1") do
            @test_throws ErrorException helper(p_dummy, 1, err_eff, stacktrace())
        end
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
                debug_key = "SPACEAGORA_DEBUG_THRUSTER"
                old_debug = get(ENV, debug_key, nothing)
                try
                    if haskey(ENV, debug_key)
                        delete!(ENV, debug_key)
                    end
                    thrust_calculation_schmitt_trigger!(link, thr, 1.0, 0.0)
                    @test thr.thrust == 0.0
                    @test !isfile("thruster_debug.csv")

                    thrust_calculation_schmitt_trigger!(link, thr, 9.0, 0.1)
                    @test thr.thrust > 0.0
                    @test thr.thrust <= thr.max_thrust + 1e-9
                    @test !isfile("thruster_debug.csv")

                    ENV[debug_key] = "1"
                    thrust_calculation_schmitt_trigger!(link, thr, 9.0, 0.2)
                    @test isfile("thruster_debug.csv")
                finally
                    if old_debug === nothing
                        if haskey(ENV, debug_key)
                            delete!(ENV, debug_key)
                        end
                    else
                        ENV[debug_key] = old_debug
                    end
                end

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

@testset "Control Mass Flow Coupling (End-to-End)" begin
    sc_no_control = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
    args_no_control = build_config(
        spacecraft=sc_no_control,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=180.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(),
        control_rates=Float64[],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    df_no_control = run_case_silent(args_no_control)
    @test "sc1_mass" in names(df_no_control)
    mass_no_control = Vector{Float64}(df_no_control.sc1_mass)
    @test maximum(abs.(mass_no_control .- first(mass_no_control))) < 1e-9

    function run_burn_case(isp::Float64)
        sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
        thruster = make_base_thruster_model(
            thrust=600.0,
            direction=0.0,
            Δv=0.0,
            start_burn_time=0.0,
            stop_burn_time=80.0,
            Isp=isp
        )
        args = build_config(
            spacecraft=sc,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=180.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            control_effectors=(thruster,),
            control_rates=[1.0],
            keplerian=true,
            tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
        )
        return run_case_silent(args)
    end

    df_low_isp = run_burn_case(200.0)
    df_high_isp = run_burn_case(400.0)

    mass_low = Vector{Float64}(df_low_isp.sc1_mass)
    mass_high = Vector{Float64}(df_high_isp.sc1_mass)
    Δm_low = first(mass_low) - last(mass_low)
    Δm_high = first(mass_high) - last(mass_high)

    @test Δm_low > 5.0
    @test Δm_high > 2.0
    @test all(diff(mass_low) .<= 1e-7)
    @test all(diff(mass_high) .<= 1e-7)
    @test isapprox(Δm_low / Δm_high, 2.0; atol=0.08, rtol=0.0)
end

@testset "Control Callback Multi-Spacecraft Mapping" begin
    sc1 = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=120.0)
    sc2 = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=120.0)
    shared_thruster = BaseThrusterModel(
        thrust=[800.0, 800.0],
        direction=[0.0, π],
        Δv=[20.0, 20.0],
        start_burn_time=[-1.0, -1.0],
        stop_burn_time=[-1.0, -1.0],
        Isp=[300.0, 300.0]
    )
    args = build_config_multi(
        spacecraft=[sc1, sc2],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=1000.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(shared_thruster,),
        control_rates=[1.0],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    df = run_case_silent(args; isolate_state=false)

    @test shared_thruster.start_burn_time[1] != -1.0
    @test shared_thruster.stop_burn_time[1] != -1.0
    @test shared_thruster.start_burn_time[2] != -1.0
    @test shared_thruster.stop_burn_time[2] != -1.0

    eps1 = 0.5 .* (df.sc1_vel_1.^2 .+ df.sc1_vel_2.^2 .+ df.sc1_vel_3.^2) .-
           EARTH.μ ./ sqrt.(df.sc1_pos_1.^2 .+ df.sc1_pos_2.^2 .+ df.sc1_pos_3.^2)
    eps2 = 0.5 .* (df.sc2_vel_1.^2 .+ df.sc2_vel_2.^2 .+ df.sc2_vel_3.^2) .-
           EARTH.μ ./ sqrt.(df.sc2_pos_1.^2 .+ df.sc2_pos_2.^2 .+ df.sc2_pos_3.^2)
    Δeps1 = last(eps1) - first(eps1)
    Δeps2 = last(eps2) - first(eps2)
    @test Δeps1 > 2e3
    @test Δeps2 < -2e3
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
