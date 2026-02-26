const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_OUTPUT_DIR = joinpath(REPO_ROOT, "output", "performance")
const SPICE_PATH = joinpath(REPO_ROOT, "GRAM_Data", "SPICE")
const EARTH_HARMONICS_FILE = joinpath(REPO_ROOT, "Gravity_harmonics_data", "EarthGGM05C.csv")

using CSV
using DataFrames
using Dates
using LinearAlgebra
using Random
using SPICE
using StaticArrays
using Statistics

include(joinpath(REPO_ROOT, "src", "simulation_model", "SimulationModel.jl"))
using .SimulationModel

# run_simulation.jl expects quat_mult in the including scope.
const quat_mult = SimulationModel.quat_mult
include(joinpath(REPO_ROOT, "src", "simulation", "run_simulation.jl"))

Base.@kwdef struct ProfileSpec
    name::String
    repeats::Int
    warmup::Int
    max_attempts::Int
    mission_short_s::Float64
    mission_long_s::Float64
    montecarlo_samples::Int
    montecarlo_mission_s::Float64
end

Base.@kwdef struct BenchmarkCase
    name::String
    category::String
    description::String
    args_template::SimulationConfiguration
    run_in_quick::Bool = true
end

@inline _safe_div(num::Float64, den::Float64) = den > 0.0 ? num / den : NaN

@inline function _alloc_calls(gcstats::Base.GC_Diff)::Int
    return gcstats.malloc + gcstats.realloc + gcstats.poolalloc + gcstats.bigalloc
end

@inline function _destat_int(sol, field::Symbol)
    if hasproperty(sol.destats, field)
        value = getproperty(sol.destats, field)
        return value isa Integer ? Int(value) : try
            Int(value)
        catch
            missing
        end
    end
    return missing
end

@inline function _solve_success(sol)::Bool
    retcode = string(sol.retcode)
    return retcode == "Success" || retcode == "Terminated"
end

@inline function _effector_signature(effectors::Tuple)
    isempty(effectors) && return "none"
    return join([string(nameof(typeof(e))) for e in effectors], "+")
end

@inline function _profile_from_name(name::String)::ProfileSpec
    if name == "full"
        return ProfileSpec(
            name="full",
            repeats=5,
            warmup=2,
            max_attempts=4,
            mission_short_s=3600.0,
            mission_long_s=14400.0,
            montecarlo_samples=20,
            montecarlo_mission_s=3600.0
        )
    elseif name == "quick"
        return ProfileSpec(
            name="quick",
            repeats=3,
            warmup=1,
            max_attempts=4,
            mission_short_s=1800.0,
            mission_long_s=7200.0,
            montecarlo_samples=8,
            montecarlo_mission_s=1800.0
        )
    else
        throw(ArgumentError("Unsupported profile '$name'. Valid values: quick, full."))
    end
end

function parse_cli()
    profile_name = lowercase(get(ENV, "SPACEAGORA_PERF_PROFILE", "quick"))
    outdir = get(ENV, "SPACEAGORA_PERF_OUTDIR", DEFAULT_OUTPUT_DIR)

    for arg in ARGS
        if arg in ("quick", "full")
            profile_name = arg
        elseif startswith(arg, "--profile=")
            profile_name = lowercase(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--outdir=")
            outdir = split(arg, "=", limit=2)[2]
        else
            throw(ArgumentError("Unknown argument '$arg'. Supported: [quick|full], --profile=..., --outdir=..."))
        end
    end

    return _profile_from_name(profile_name), abspath(outdir)
end

function make_spacecraft(
    planet::AbstractPlanet;
    id::Int=1,
    ra_alt_m::Float64=550e3,
    rp_alt_m::Float64=500e3,
    i_deg::Float64=35.0,
    ω_deg::Float64=40.0,
    Ω_deg::Float64=10.0,
    ν_deg::Float64=170.0,
    with_panel::Bool=true,
    orientation_state::Union{Nothing, Tuple{SVector{4, Float64}, SVector{3, Float64}}}=nothing,
    root_mass::Float64=500.0,
    root_area::Float64=12.0,
    panel_mass::Float64=30.0,
    panel_area::Float64=6.0,
    panel_offset_y::Float64=1.2,
    prop_mass::Float64=0.0
)::SpacecraftModel
    root = Link{0}(root=true, m=root_mass, ref_area=root_area)
    links = Link[root]
    if with_panel
        panel = Link{0}(root=false, m=panel_mass, ref_area=panel_area, r=MVector{3, Float64}(0.0, panel_offset_y, 0.0))
        push!(links, panel)
    end

    ra = planet.Rp_e + ra_alt_m
    rp = planet.Rp_e + rp_alt_m

    ic = if isnothing(orientation_state)
        InitialCondition(ra=ra, rp=rp, i=i_deg, ω=ω_deg, Ω=Ω_deg, ν=ν_deg)
    else
        q0, w0 = orientation_state
        a = (ra + rp) / 2.0
        e = (ra - rp) / (ra + rp)
        InitialCondition(a, e, i_deg, ω_deg, Ω_deg, ν_deg, q0, w0)
    end

    dry_mass = sum(link.m for link in links)
    return SpacecraftModel(
        Joint[],
        links,
        root,
        true,
        dry_mass,
        prop_mass,
        root.inertia,
        0,
        0,
        ic,
        id
    )
end

function make_constellation(planet::AbstractPlanet, n::Int; with_panel::Bool=false)::Vector{SpacecraftModel}
    sats = SpacecraftModel[]
    for i in 1:n
        ra_alt = 540e3 + 20e3 * (i - 1)
        rp_alt = 470e3 + 15e3 * (i - 1)
        if rp_alt >= ra_alt
            ra_alt = rp_alt + 8e3
        end
        ν = 140.0 + 180.0 * (i - 1) / n
        push!(sats, make_spacecraft(planet; id=i, ra_alt_m=ra_alt, rp_alt_m=rp_alt, ν_deg=ν, with_panel=with_panel))
    end
    return sats
end

function build_config(;
    planet::AbstractPlanet,
    spacecraft::Vector{SpacecraftModel},
    mission_time_s::Float64,
    orientation_sim::Bool,
    dynamic_effectors::Tuple,
    control_effectors::Tuple=(),
    control_rates::Vector{Float64}=Float64[],
    dt_max_orbit::Float64=1.0,
    reltol_orbit::Float64=1e-9,
    abstol_orbit::Float64=1e-9
)::SimulationConfiguration
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false,
            verbose=false,
            generate_plots=false,
            normalize=false,
            save_csv=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=mission_time_s,
            orientation_sim=orientation_sim,
            num_steps_to_save=400
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false
        ),
        dynamics_model=DynamicsModel(spacecraft, dynamic_effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=control_effectors, control_rates=control_rates),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=reltol_orbit,
            abstol_orbit=abstol_orbit,
            dt_max_orbit=dt_max_orbit
        )
    )
end

function make_montecarlo_config(seed::Int, planet::Earth, mission_time_s::Float64)::SimulationConfiguration
    rng = MersenneTwister(seed)
    ra_alt = 530e3 + randn(rng) * 20e3
    rp_alt = 490e3 + randn(rng) * 20e3
    rp_alt = max(rp_alt, 120e3)
    if rp_alt >= ra_alt
        ra_alt = rp_alt + 8e3
    end

    sc = make_spacecraft(
        planet;
        id=1,
        with_panel=false,
        root_mass=160.0,
        root_area=1.5,
        prop_mass=20.0,
        ra_alt_m=ra_alt,
        rp_alt_m=rp_alt,
        i_deg=28.0 + randn(rng) * 0.8,
        ω_deg=15.0 + randn(rng) * 2.0,
        Ω_deg=20.0 + randn(rng) * 2.0,
        ν_deg=160.0 + randn(rng) * 8.0
    )

    thruster = BaseThrusterModel(
        thrust=[0.6 + rand(rng) * 0.5],
        direction=[0.0],
        Δv=[0.0],
        start_burn_time=[0.0],
        stop_burn_time=[180.0 + rand(rng) * 40.0],
        Isp=[280.0 + rand(rng) * 40.0]
    )

    return build_config(
        planet=planet,
        spacecraft=[sc],
        mission_time_s=mission_time_s,
        orientation_sim=false,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(thruster,),
        control_rates=[1.0],
        dt_max_orbit=1.0,
        reltol_orbit=1e-9,
        abstol_orbit=1e-9
    )
end

function build_cases(spec::ProfileSpec, planet::Earth)::Vector{BenchmarkCase}
    harmonics20 = GravitationalHarmonicsModel(20, 20, EARTH_HARMONICS_FILE, planet)
    harmonics50 = GravitationalHarmonicsModel(50, 50, EARTH_HARMONICS_FILE, planet)
    nbody_sun_moon = NBodyGravityModel(body_names=("Sun", "Moon"), primary_body_name="Earth", planet=planet)

    q0 = normalize(SVector{4, Float64}(0.15, -0.05, 0.2, 0.96))
    w0 = SVector{3, Float64}(0.01, -0.02, 0.015)

    sc_baseline = [make_spacecraft(planet; id=1, with_panel=false)]
    sc_orientation = [make_spacecraft(planet; id=1, with_panel=true, orientation_state=(q0, w0))]
    sc_thruster = [make_spacecraft(planet; id=1, with_panel=false, prop_mass=30.0)]
    thruster = BaseThrusterModel(
        thrust=[0.8],
        direction=[0.0],
        Δv=[35.0],
        start_burn_time=[-1.0],
        stop_burn_time=[-1.0],
        Isp=[300.0]
    )

    cases = BenchmarkCase[
        BenchmarkCase(
            name="single_baseline_gravity",
            category="core",
            description="1 spacecraft, position-only, inverse-square gravity",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(),)
            )
        ),
        BenchmarkCase(
            name="single_orientation_aero",
            category="orientation",
            description="1 spacecraft, orientation dynamics on, aerodynamic model active",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_orientation,
                mission_time_s=spec.mission_short_s,
                orientation_sim=true,
                dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM())
            )
        ),
        BenchmarkCase(
            name="multi_2_gravity",
            category="satellite_scaling",
            description="2 spacecraft, position-only, inverse-square gravity",
            args_template=build_config(
                planet=planet,
                spacecraft=make_constellation(planet, 2; with_panel=false),
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(),)
            )
        ),
        BenchmarkCase(
            name="multi_4_gravity",
            category="satellite_scaling",
            description="4 spacecraft, position-only, inverse-square gravity",
            args_template=build_config(
                planet=planet,
                spacecraft=make_constellation(planet, 4; with_panel=false),
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(),)
            )
        ),
        BenchmarkCase(
            name="multi_8_gravity",
            category="satellite_scaling",
            description="8 spacecraft, position-only, inverse-square gravity",
            args_template=build_config(
                planet=planet,
                spacecraft=make_constellation(planet, 8; with_panel=false),
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(),)
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="single_j2",
            category="dynamics_fidelity",
            description="1 spacecraft, inverse-square + J2 gravity",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredJ2GravityModel(),)
            )
        ),
        BenchmarkCase(
            name="single_nbody_sun_moon",
            category="dynamics_fidelity",
            description="1 spacecraft, inverse-square gravity + N-body Sun/Moon perturbations",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), nbody_sun_moon)
            )
        ),
        BenchmarkCase(
            name="single_harmonics_l20",
            category="dynamics_fidelity",
            description="1 spacecraft, inverse-square gravity + spherical harmonics L=M=20",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), harmonics20)
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="single_harmonics_l50",
            category="dynamics_fidelity",
            description="1 spacecraft, inverse-square gravity + spherical harmonics L=M=50",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), harmonics50)
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="single_thruster_control",
            category="control",
            description="1 spacecraft, inverse-square gravity + BaseThrusterModel control callback",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_thruster,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(),),
                control_effectors=(thruster,),
                control_rates=[1.0]
            )
        ),
        BenchmarkCase(
            name="single_baseline_long_mission",
            category="mission_length",
            description="1 spacecraft baseline gravity with longer mission duration",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_long_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(),)
            )
        )
    ]

    return cases
end

function measure_case(case::BenchmarkCase, profile_name::String, repeat_idx::Int; seed::Union{Missing, Int}=missing, attempt::Int=1)
    timestamp_utc = string(now(UTC))
    args_meta = case.args_template
    n_sats = length(args_meta.dynamics_model.spacecraft)
    mission_time_s = args_meta.mission_configuration.mission_time

    GC.gc()
    copy_timed = @timed deepcopy(case.args_template)
    args_run = copy_timed.value

    GC.gc()
    solve_timed = @timed run_simulation(args_run; isolate_state=false, return_solution=true)
    sol = solve_timed.value
    solve_retcode = string(sol.retcode)
    solve_success = _solve_success(sol)

    copy_time_s = copy_timed.time
    solve_time_s = solve_timed.time
    total_time_s = copy_time_s + solve_time_s

    copy_bytes_mb = copy_timed.bytes / 1024^2
    solve_bytes_mb = solve_timed.bytes / 1024^2
    total_bytes_mb = copy_bytes_mb + solve_bytes_mb

    copy_alloc_calls = _alloc_calls(copy_timed.gcstats)
    solve_alloc_calls = _alloc_calls(solve_timed.gcstats)

    return (
        profile=profile_name,
        category=case.category,
        scenario=case.name,
        description=case.description,
        seed=seed,
        repeat=repeat_idx,
        attempt=attempt,
        satellites=n_sats,
        orientation=args_meta.mission_configuration.orientation_sim,
        mission_time_s=mission_time_s,
        dynamic_effectors=_effector_signature(args_meta.dynamics_model.dynamic_effectors),
        control_effectors=_effector_signature(args_meta.control_model.control_effectors),
        copy_time_s=copy_time_s,
        solve_time_s=solve_time_s,
        total_time_s=total_time_s,
        copy_compile_time_s=copy_timed.compile_time,
        solve_compile_time_s=solve_timed.compile_time,
        copy_gctime_s=copy_timed.gctime,
        solve_gctime_s=solve_timed.gctime,
        solve_retcode=solve_retcode,
        solve_success=solve_success,
        copy_bytes_mb=copy_bytes_mb,
        solve_bytes_mb=solve_bytes_mb,
        total_bytes_mb=total_bytes_mb,
        copy_alloc_calls=copy_alloc_calls,
        solve_alloc_calls=solve_alloc_calls,
        saved_points=length(sol.t),
        accepted_steps=_destat_int(sol, :naccept),
        rejected_steps=_destat_int(sol, :nreject),
        sim_seconds_per_wall_second=_safe_div(mission_time_s, total_time_s),
        satellite_sim_seconds_per_wall_second=_safe_div(mission_time_s * n_sats, total_time_s),
        timestamp_utc=timestamp_utc
    )
end

function run_warmup(case::BenchmarkCase, warmup::Int)
    for _ in 1:warmup
        args_run = deepcopy(case.args_template)
        run_simulation(args_run; isolate_state=false, return_solution=false)
    end
    return nothing
end

function run_case_batch!(rows::Vector{NamedTuple}, case::BenchmarkCase, spec::ProfileSpec, idx::Int, total::Int)
    println("[$idx/$total] $(case.name) :: warmup x$(spec.warmup), repeats x$(spec.repeats)")
    run_warmup(case, spec.warmup)
    for rep in 1:spec.repeats
        last_row = nothing
        for attempt in 1:spec.max_attempts
            row = measure_case(case, spec.name, rep; attempt=attempt)
            last_row = row
            if row.solve_success
                push!(rows, row)
                println("  repeat $(rep)/$(spec.repeats): total=$(round(row.total_time_s; digits=3)) s, solve=$(round(row.solve_time_s; digits=3)) s")
                break
            end
            println("  repeat $(rep)/$(spec.repeats) attempt $(attempt)/$(spec.max_attempts): failed retcode=$(row.solve_retcode), retrying")
        end
        if !(last_row === nothing) && !last_row.solve_success
            push!(rows, last_row)
            println("  repeat $(rep)/$(spec.repeats): failed after $(spec.max_attempts) attempts, retcode=$(last_row.solve_retcode)")
        end
    end
    return nothing
end

function run_montecarlo_batch!(rows::Vector{NamedTuple}, spec::ProfileSpec, planet::Earth)
    seeds = collect(1001:(1000 + spec.montecarlo_samples))
    warmup_case = BenchmarkCase(
        name="montecarlo_randomized",
        category="montecarlo",
        description="Randomized initial conditions + thruster, one run per seed",
        args_template=make_montecarlo_config(first(seeds), planet, spec.montecarlo_mission_s),
        run_in_quick=true
    )
    println("[montecarlo] warmup x$(spec.warmup), seeds=$(length(seeds))")
    run_warmup(warmup_case, spec.warmup)

    for (i, seed) in enumerate(seeds)
        case = BenchmarkCase(
            name="montecarlo_randomized",
            category="montecarlo",
            description="Randomized initial conditions + thruster, one run per seed",
            args_template=make_montecarlo_config(seed, planet, spec.montecarlo_mission_s),
            run_in_quick=true
        )
        last_row = nothing
        for attempt in 1:spec.max_attempts
            row = measure_case(case, spec.name, 1; seed=seed, attempt=attempt)
            last_row = row
            if row.solve_success
                push!(rows, row)
                println("  seed $(i)/$(length(seeds))=$(seed): total=$(round(row.total_time_s; digits=3)) s")
                break
            end
            println("  seed $(i)/$(length(seeds))=$(seed) attempt $(attempt)/$(spec.max_attempts): failed retcode=$(row.solve_retcode), retrying")
        end
        if !(last_row === nothing) && !last_row.solve_success
            push!(rows, last_row)
            println("  seed $(i)/$(length(seeds))=$(seed): failed after $(spec.max_attempts) attempts, retcode=$(last_row.solve_retcode)")
        end
    end
    return nothing
end

function run_benchmarks(spec::ProfileSpec, cases::Vector{BenchmarkCase}, planet::Earth)::DataFrame
    selected = spec.name == "full" ? cases : [c for c in cases if c.run_in_quick]
    rows = NamedTuple[]
    total = length(selected)

    for (idx, case) in enumerate(selected)
        run_case_batch!(rows, case, spec, idx, total)
    end

    run_montecarlo_batch!(rows, spec, planet)
    return DataFrame(rows)
end

function _safe_stat(values, op::Function)
    vec = collect(skipmissing(values))
    isempty(vec) && return missing
    return op(vec)
end

function summarize_results(raw_df::DataFrame)::DataFrame
    keys = [:category, :scenario, :description, :satellites, :orientation, :mission_time_s, :dynamic_effectors, :control_effectors]
    counts = combine(
        groupby(raw_df, keys),
        nrow => :samples_total,
        :solve_success => (v -> count(identity, v)) => :samples_success
    )
    counts[!, :samples_failed] = counts.samples_total .- counts.samples_success
    counts[!, :success_rate] = Float64.(counts.samples_success) ./ Float64.(counts.samples_total)

    success_df = raw_df[raw_df.solve_success .== true, :]
    metric_cols = [
        :samples,
        :copy_time_mean_s,
        :solve_time_mean_s,
        :total_time_mean_s,
        :total_time_median_s,
        :total_time_std_s,
        :total_time_min_s,
        :total_time_max_s,
        :total_time_p90_s,
        :total_bytes_mean_mb,
        :copy_alloc_mean,
        :solve_alloc_mean,
        :saved_points_mean,
        :accepted_steps_mean,
        :rejected_steps_mean,
        :sim_seconds_per_wall_second_mean,
        :satellite_sim_seconds_per_wall_second_mean
    ]

    summary = counts
    if nrow(success_df) > 0
        success_summary = combine(
            groupby(success_df, keys),
            nrow => :samples,
            :copy_time_s => (v -> _safe_stat(v, mean)) => :copy_time_mean_s,
            :solve_time_s => (v -> _safe_stat(v, mean)) => :solve_time_mean_s,
            :total_time_s => (v -> _safe_stat(v, mean)) => :total_time_mean_s,
            :total_time_s => (v -> _safe_stat(v, median)) => :total_time_median_s,
            :total_time_s => (v -> _safe_stat(v, x -> std(x; corrected=false))) => :total_time_std_s,
            :total_time_s => (v -> _safe_stat(v, minimum)) => :total_time_min_s,
            :total_time_s => (v -> _safe_stat(v, maximum)) => :total_time_max_s,
            :total_time_s => (v -> _safe_stat(v, x -> quantile(x, 0.9))) => :total_time_p90_s,
            :total_bytes_mb => (v -> _safe_stat(v, mean)) => :total_bytes_mean_mb,
            :copy_alloc_calls => (v -> _safe_stat(v, mean)) => :copy_alloc_mean,
            :solve_alloc_calls => (v -> _safe_stat(v, mean)) => :solve_alloc_mean,
            :saved_points => (v -> _safe_stat(v, mean)) => :saved_points_mean,
            :accepted_steps => (v -> _safe_stat(v, mean)) => :accepted_steps_mean,
            :rejected_steps => (v -> _safe_stat(v, mean)) => :rejected_steps_mean,
            :sim_seconds_per_wall_second => (v -> _safe_stat(v, mean)) => :sim_seconds_per_wall_second_mean,
            :satellite_sim_seconds_per_wall_second => (v -> _safe_stat(v, mean)) => :satellite_sim_seconds_per_wall_second_mean
        )
        summary = leftjoin(counts, success_summary, on=keys)
    else
        for col in metric_cols
            summary[!, col] = fill(missing, nrow(summary))
        end
    end

    baseline_idx = findfirst(==("single_baseline_gravity"), summary.scenario)
    if baseline_idx === nothing || ismissing(summary.total_time_mean_s[baseline_idx]) || summary.total_time_mean_s[baseline_idx] <= 0.0
        summary[!, :relative_to_baseline] = fill(missing, nrow(summary))
        summary[!, :speedup_vs_baseline] = fill(missing, nrow(summary))
    else
        baseline = summary.total_time_mean_s[baseline_idx]
        summary[!, :relative_to_baseline] = [
            ismissing(row_total) ? missing : row_total / baseline
            for row_total in summary.total_time_mean_s
        ]
        summary[!, :speedup_vs_baseline] = [
            ismissing(row_total) || row_total <= 0.0 ? missing : baseline / row_total
            for row_total in summary.total_time_mean_s
        ]
    end

    summary[!, :_sort_key] = [ismissing(x) ? -Inf : x for x in summary.total_time_mean_s]
    sort!(summary, :_sort_key, rev=true)
    select!(summary, Not(:_sort_key))
    return summary
end

@inline function _fmt(v; digits::Int=3)
    if v isa Missing
        return "n/a"
    elseif v isa AbstractFloat
        return isfinite(v) ? string(round(v; digits=digits)) : "n/a"
    else
        return string(v)
    end
end

function _scenario_metric(summary_df::DataFrame, scenario::String, metric::Symbol)
    idx = findfirst(==(scenario), summary_df.scenario)
    idx === nothing && return nothing
    return summary_df[idx, metric]
end

function write_report(path::String, spec::ProfileSpec, raw_df::DataFrame, summary_df::DataFrame)
    generated = string(now(UTC))
    julia_ver = string(VERSION)
    nthreads = Threads.nthreads()

    valid_rows = findall(x -> !ismissing(x), summary_df.total_time_mean_s)
    fastest = nothing
    slowest = nothing
    if !isempty(valid_rows)
        vals = summary_df.total_time_mean_s[valid_rows]
        fastest = summary_df[valid_rows[argmin(vals)], :]
        slowest = summary_df[valid_rows[argmax(vals)], :]
    end

    baseline = _scenario_metric(summary_df, "single_baseline_gravity", :total_time_mean_s)
    orientation = _scenario_metric(summary_df, "single_orientation_aero", :total_time_mean_s)
    multi2 = _scenario_metric(summary_df, "multi_2_gravity", :total_time_mean_s)
    multi4 = _scenario_metric(summary_df, "multi_4_gravity", :total_time_mean_s)
    harmonics20 = _scenario_metric(summary_df, "single_harmonics_l20", :total_time_mean_s)
    harmonics50 = _scenario_metric(summary_df, "single_harmonics_l50", :total_time_mean_s)

    mc_rows = raw_df[(raw_df.category .== "montecarlo") .& (raw_df.solve_success .== true), :]
    mc_mean = nrow(mc_rows) > 0 ? mean(mc_rows.total_time_s) : missing
    mc_p90 = nrow(mc_rows) > 0 ? quantile(mc_rows.total_time_s, 0.9) : missing
    mc_std = nrow(mc_rows) > 0 ? std(mc_rows.total_time_s; corrected=false) : missing

    open(path, "w") do io
        println(io, "# SpaceAGORA Computational Time Analysis (`$(spec.name)` profile)")
        println(io)
        println(io, "- Generated (UTC): $generated")
        println(io, "- Julia: `$julia_ver`")
        println(io, "- Threads: `$nthreads`")
        println(io, "- Repeats per deterministic scenario: `$(spec.repeats)`")
        println(io, "- Warmup runs per scenario: `$(spec.warmup)`")
        println(io, "- Monte Carlo seeds: `$(spec.montecarlo_samples)`")
        println(io)
        println(io, "## Key Findings")
        println(io)
        if fastest === nothing || slowest === nothing
            println(io, "- No successful runs were recorded.")
        else
            println(io, "- Slowest successful scenario: `$(slowest.scenario)` with mean total time `$(round(slowest.total_time_mean_s; digits=3)) s`.")
            println(io, "- Fastest successful scenario: `$(fastest.scenario)` with mean total time `$(round(fastest.total_time_mean_s; digits=3)) s`.")
        end
        if baseline !== nothing && orientation !== nothing && !ismissing(baseline) && !ismissing(orientation) && baseline > 0.0
            println(io, "- Orientation + aerodynamic run vs baseline: `$(round(orientation / baseline; digits=2))x` runtime.")
        end
        if baseline !== nothing && multi2 !== nothing && multi4 !== nothing && !ismissing(baseline) && !ismissing(multi2) && !ismissing(multi4) && baseline > 0.0
            println(io, "- Multi-satellite scaling: `2-sat=$(round(multi2 / baseline; digits=2))x`, `4-sat=$(round(multi4 / baseline; digits=2))x` relative to single-sat baseline.")
        end
        if harmonics20 !== nothing && harmonics50 !== nothing && !ismissing(harmonics20) && !ismissing(harmonics50) && harmonics20 > 0.0
            println(io, "- Harmonics scaling: `L=50` is `$(round(harmonics50 / harmonics20; digits=2))x` relative to `L=20`.")
        end
        if !(mc_mean isa Missing)
            println(io, "- Monte Carlo runtime spread: mean `$(round(mc_mean; digits=3)) s`, p90 `$(round(mc_p90; digits=3)) s`, std `$(round(mc_std; digits=3)) s`.")
        end
        failed_groups = summary_df[summary_df.samples_failed .> 0, :]
        if nrow(failed_groups) > 0
            println(io, "- Solver failures detected in `$(nrow(failed_groups))` scenario groups; timings only use successful runs.")
        end
        println(io)
        println(io, "## Scenario Summary")
        println(io)
        println(io, "| Scenario | Category | Success/Total | Mean Total (s) | P90 (s) | Mean Solve (s) | Mean Copy (s) | Sim sec / wall sec | Rel. Baseline |")
        println(io, "|---|---|---:|---:|---:|---:|---:|---:|---:|")
        for row in eachrow(summary_df)
            println(
                io,
                "| $(row.scenario) | $(row.category) | $(row.samples_success)/$(row.samples_total) | $(_fmt(row.total_time_mean_s)) | $(_fmt(row.total_time_p90_s)) | $(_fmt(row.solve_time_mean_s)) | $(_fmt(row.copy_time_mean_s)) | $(_fmt(row.sim_seconds_per_wall_second_mean)) | $(_fmt(row.relative_to_baseline)) |"
            )
        end
    end
end

function main()
    spec, outdir = parse_cli()
    mkpath(outdir)

    println("Performance runtime analysis profile: $(spec.name)")
    println("Output directory: $outdir")

    planet = Earth("", SPICE_PATH)
    cases = build_cases(spec, planet)
    raw_df = run_benchmarks(spec, cases, planet)
    summary_df = summarize_results(raw_df)

    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    raw_path = joinpath(outdir, "runtime_raw_$(spec.name)_$(stamp).csv")
    summary_path = joinpath(outdir, "runtime_summary_$(spec.name)_$(stamp).csv")
    report_path = joinpath(outdir, "runtime_report_$(spec.name)_$(stamp).md")

    CSV.write(raw_path, raw_df)
    CSV.write(summary_path, summary_df)
    write_report(report_path, spec, raw_df, summary_df)

    println("Analysis complete.")
    println("Raw results: $raw_path")
    println("Summary: $summary_path")
    println("Report: $report_path")
end

main()
