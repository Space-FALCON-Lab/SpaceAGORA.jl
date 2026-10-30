# Attribution runs: the same calibrated initial state under model variants (Sun/Moon/J3/drag/frozen ephemeris).
# Process route via run_monte_carlo. No src/ change: variants are built here by editing the effector tuple.
#   julia --project=. run_attribution.jl [--pilot] --out=DIR [--workers=N]
include(joinpath(@__DIR__, "run_static.jl"))

const ATTR_CONFIG_PATH = joinpath(@__DIR__, "attribution_config.toml")

"""Tidal (direct + indirect) point-mass acceleration of fixed Earth-centered body positions (same expression as the engine's NBodyGravityModel kernel)."""
struct FrozenTidalModel <: AbstractForceTorqueModel
    positions_ii::Vector{SVector{3, Float64}}
    mus::Vector{Float64}
end
SimulationModel.environment_requirements(::FrozenTidalModel) = SimulationModel.EffectorEnvironmentRequirements()
function SimulationModel.wrench(m::FrozenTidalModel, x::SimulationModel.StateSample, env::SimulationModel.EnvironmentSample, t::Float64)
    acc = SVector{3, Float64}(0.0, 0.0, 0.0)
    for k in eachindex(m.mus)
        pb = m.positions_ii[k]; d = pb - x.pos_ii
        acc += m.mus[k] * (d / norm(d)^3 - pb / norm(pb)^3)
    end
    return x.mass_kg * acc, SVector{3, Float64}(0.0, 0.0, 0.0)
end

variant_spec(v) = Dict(
    "v0" => (bodies=("Moon", "Sun"), frozen=false, L=3, drag=true), "v1" => (bodies=("Moon",), frozen=false, L=3, drag=true),
    "v2" => (bodies=("Sun",), frozen=false, L=3, drag=true), "v3" => (bodies=(), frozen=false, L=3, drag=true),
    "v4" => (bodies=("Moon", "Sun"), frozen=false, L=2, drag=true), "v5" => (bodies=("Moon", "Sun"), frozen=false, L=3, drag=false),
    "v6" => (bodies=("Moon",), frozen=true, L=3, drag=true), "v7" => (bodies=("Sun",), frozen=true, L=3, drag=true),
    "v8" => (bodies=("Moon",), frozen=false, L=3, drag=false), "v9" => (bodies=("Sun",), frozen=false, L=3, drag=false))[v]

function variant_args(cfg, study, epoch, horizon_s, dir, v)
    sp = variant_spec(v)
    base = build_args(cfg, study, epoch, horizon_s, dir; third_bodies=false)   # (harmonics L=3, drag); initial state from the calibrated cfg
    planet = base.environment_model.planet
    eff = Any[]
    for e in base.dynamics_model.dynamic_effectors
        if e isa GravitationalHarmonicsModel
            push!(eff, sp.L == 3 ? e : GravitationalHarmonicsModel(sp.L, 0, joinpath(REPO_ROOT, String(cfg.values["gravity"]["coefficient_file"])), planet))
        elseif e isa CorridorBallisticDrag
            sp.drag && push!(eff, e)
        else
            push!(eff, e)
        end
    end
    if !isempty(sp.bodies)
        if sp.frozen
            et0 = SimulationModel.ephemerides_time_seconds(base.initial_time, base.environment_model.ephemerides_model)
            pos = [SVector{3, Float64}(spkpos(uppercase(b), et0, "J2000", "NONE", "EARTH")[1]) .* 1e3 for b in sp.bodies]
            push!(eff, FrozenTidalModel(pos, [MU_BODY[b] for b in sp.bodies]))
        else
            push!(eff, NBodyGravityModel(body_names=Tuple(String.(sp.bodies)), primary_body_name="Earth", planet=nbody_planet()))
        end
    end
    return SimulationConfiguration(
        file_paths=base.file_paths, simulation_settings=base.simulation_settings, mission_configuration=base.mission_configuration,
        environment_model=base.environment_model, dynamics_model=DynamicsModel(base.dynamics_model.spacecraft, Tuple(eff)),
        guidance_model=base.guidance_model, navigation_model=base.navigation_model, control_model=base.control_model,
        initial_time=base.initial_time, integration_tolerances=base.integration_tolerances, solver_config=base.solver_config)
end

function attr_task(seed)
    study = load_study(seed.study_path); nbody_planet()
    sdir = joinpath(seed.ref_dir, "grid", seed.src_tag)
    cfg = cell_config(study, seed.inc, seed.raan, seed.omega; seed_alt=seed_alt_of(sdir))
    epoch = parse_epoch(study["grid"]["epochs_utc"][seed.epoch_idx])
    dir = joinpath(seed.out, "runs", seed.tag)
    args = variant_args(cfg, study, epoch, Float64(seed.horizon_days) * SECONDS_PER_DAY, dir, seed.variant)
    t0 = time(); sol = run_simulation(args; return_solution=true, isolate_state=false); wall = time() - t0
    rows = periapsis_table(sol, args, cfg, study)
    write_rows(joinpath(dir, "passes_$(seed.variant).csv"), rows)
    return (; tag=seed.tag, variant=seed.variant, wall_s=wall, retcode=string(sol.retcode), passes=length(rows), seed_alt_m=seed_alt_of(sdir),
        t_end_days=isempty(rows) ? NaN : rows[end].time_s / SECONDS_PER_DAY)
end

function attr_main(args=ARGS)
    o = parse_args(args); study = load_study(); ac = TOML.parsefile(ATTR_CONFIG_PATH); out = o["out"]; mkpath(out)
    workers = something(o["workers"], 4); hd = 60.0
    files = [joinpath(@__DIR__, "run_attribution.jl")]
    root = results_root(study); refdirs = Dict("border2" => joinpath(root, "border2"), "border" => joinpath(root, "border"))
    seeds = NamedTuple[]
    for c in ac["cells"]["list"]
        inc, raan, omega, ei, src = c
        tag = cell_tag(inc, raan, omega, ei)
        src_tag = src == "border" ? "R$(lpad(round(Int, raan), 3, '0'))_W$(lpad(round(Int, omega), 3, '0'))_E$(ei)" : tag
        vs = copy(ac["variants"]["all_cells"]); tag == ac["variants"]["frozen_cell"] && append!(vs, ac["variants"]["frozen"])
        for v in vs
            push!(seeds, (; study_path=STUDY_CONFIG_PATH, out, ref_dir=refdirs[src], src_tag, tag, inc, raan, omega, epoch_idx=Int(ei), horizon_days=hd, variant=v))
        end
    end
    if o["pilot"]
        t0 = time(); r = attr_task(first(filter(s -> s.variant == "v6" || s.tag == ac["variants"]["frozen_cell"] && s.variant == "v6", seeds)))
        println("PILOT v6 wall_s=$(round(time() - t0; digits=1)) $(r)"); return
    end
    res = run_process_campaign(attr_task, seeds, workers; mission_time_s=hd * SECONDS_PER_DAY, worker_files=files)
    report(res); write_rows(joinpath(out, "attribution_run_summary.csv"), [s.value for s in res.successful]); save_mc(res, out, "attribution")
end

if abspath(PROGRAM_FILE) == @__FILE__
    attr_main()
end
