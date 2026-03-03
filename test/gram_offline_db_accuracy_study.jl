const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

using Statistics
using Printf
using Dates
using LinearAlgebra
using StaticArrays
using CSV
using DataFrames
using Serialization

include(joinpath(REPO_ROOT, "src", "simulation_model", "SimulationModel.jl"))
using .SimulationModel

const EM = SimulationModel.EnvironmentModels
const CB = SimulationModel.SimulationCallbacks

function parse_cli(args::Vector{String})
    opts = Dict{String, String}()
    i = 1
    while i <= length(args)
        arg = args[i]
        startswith(arg, "--") || error("Unsupported argument '$arg'. Use --key=value.")
        body = arg[3:end]
        if occursin("=", body)
            key, value = split(body, "=", limit=2)
            opts[key] = value
        else
            key = body
            if i == length(args) || startswith(args[i + 1], "--")
                opts[key] = "true"
            else
                i += 1
                opts[key] = args[i]
            end
        end
        i += 1
    end
    return opts
end

@inline function _pctl(x::Vector{Float64}, p::Float64)
    isempty(x) && return 0.0
    xs = sort(x)
    idx = clamp(Int(ceil(p * length(xs))), 1, length(xs))
    return xs[idx]
end

function _planet_from_name(name::String, spice_path::String)
    key = lowercase(strip(name))
    if key == "earth"
        return Earth("", spice_path)
    elseif key == "mars"
        return Mars("", spice_path)
    elseif key == "venus"
        return Venus("", spice_path)
    elseif key == "titan"
        return Titan("", spice_path)
    end
    throw(ArgumentError("Unsupported planet '$name' for this analysis script. Supported: earth, mars, venus, titan"))
end

function _make_spacecraft(planet)::SpacecraftModel
    root = Link{0}(root=true, m=120.0, ref_area=1.2)
    ic = InitialCondition(
        ra=planet.Rp_e + 520e3,
        rp=planet.Rp_e + 500e3,
        i=25.0,
        ω=10.0,
        Ω=20.0,
        ν=170.0
    )
    return SpacecraftModel(
        joints=Joint[],
        links=Link[root],
        root=root,
        instant_actuation=true,
        prop_mass=0.0,
        inertia_tensor=root.inertia,
        n_reaction_wheels=0,
        n_thrusters=0,
        initial_condition=ic,
        id=1
    )
end

function _make_args(planet, density_model, initial_time)::SimulationConfiguration
    sc = _make_spacecraft(planet)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false,
            verbose=false,
            generate_plots=false,
            normalize=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=3600.0,
            orientation_sim=false,
            num_steps_to_save=1000
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=120.0,
            density_model=density_model,
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=true
        ),
        dynamics_model=DynamicsModel([sc], (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=initial_time,
        integration_tolerances=IntegrationTolerances()
    )
end

function _generate_track(
    planet;
    rp_alt_m::Float64,
    ra_alt_m::Float64,
    i_deg::Float64,
    Ω_deg::Float64,
    ω_deg::Float64,
    ν_start_deg::Float64,
    ν_end_deg::Float64,
    n_samples::Int
)
    rp = planet.Rp_e + rp_alt_m
    ra = planet.Rp_e + ra_alt_m
    a = 0.5 * (ra + rp)
    e = (ra - rp) / (ra + rp)
    i = deg2rad(i_deg)
    Ω = deg2rad(Ω_deg)
    ω = deg2rad(ω_deg)
    νs = collect(range(deg2rad(ν_start_deg), deg2rad(ν_end_deg), length=n_samples))

    pos = Vector{SVector{3, Float64}}(undef, n_samples)
    vel = Vector{SVector{3, Float64}}(undef, n_samples)
    alt = Vector{Float64}(undef, n_samples)
    lat = Vector{Float64}(undef, n_samples)
    lon = Vector{Float64}(undef, n_samples)
    t = Vector{Float64}(undef, n_samples)
    t[1] = 0.0

    for k in 1:n_samples
        oe = SVector{7, Float64}(a, e, i, Ω, ω, νs[k], 0.0)
        r_k, v_k = CB.orbitalelemtorv(oe, planet)
        pos[k] = SVector{3, Float64}(r_k)
        vel[k] = SVector{3, Float64}(v_k)
        rp_k, _ = CB.r_intor_p!(pos[k], vel[k], planet)
        alt[k], lat[k], lon[k] = CB.rtolatlong(rp_k, planet)

        if k > 1
            ds = norm(pos[k] - pos[k - 1])
            vbar = max(1.0, 0.5 * (norm(vel[k]) + norm(vel[k - 1])))
            t[k] = t[k - 1] + ds / vbar
        end
    end

    return (t=t, pos=pos, vel=vel, alt=alt, lat=lat, lon=lon)
end

function _evaluate_direct(track, density_model, p)
    n = length(track.t)
    rho = Vector{Float64}(undef, n)
    temp = Vector{Float64}(undef, n)
    wind = Vector{SVector{3, Float64}}(undef, n)

    elapsed_s = @elapsed begin
        @inbounds for k in 1:n
            rho[k], temp[k], wind[k] = Base.invokelatest(
                getDensity,
                density_model,
                track.alt[k],
                track.lat[k],
                track.lon[k],
                track.t[k],
                true,
                p
            )
        end
    end

    return (rho=rho, temp=temp, wind=wind, elapsed_s=elapsed_s)
end

function load_payload(path::String)
    open(path, "r") do io
        payload = deserialize(io)
        payload isa Dict || error("Expected Dict payload in '$path'.")
        return Dict{String, Any}(payload)
    end
end

@inline function _axis_segment(nodes::Vector{Float64}, x::Float64)
    n = length(nodes)
    n >= 2 || error("Axis must have at least two nodes.")
    xq = clamp(x, nodes[1], nodes[end])
    i0 = clamp(searchsortedlast(nodes, xq), 1, n - 1)
    i1 = i0 + 1
    x0 = nodes[i0]
    x1 = nodes[i1]
    w = x1 == x0 ? 0.0 : (xq - x0) / (x1 - x0)
    return i0, i1, clamp(w, 0.0, 1.0)
end

@inline function _lon_segment(nodes::Vector{Float64}, lon_deg::Float64)
    n = length(nodes)
    n >= 2 || error("Longitude axis must have at least two nodes.")
    period = 360.0
    xq = mod(lon_deg, period)
    xq < 0 && (xq += period)
    i0 = searchsortedlast(nodes, xq)
    i0 = i0 == 0 ? n : clamp(i0, 1, n)
    i1 = i0 == n ? 1 : i0 + 1
    x0 = nodes[i0]
    x1 = i1 == 1 ? nodes[1] + period : nodes[i1]
    xq_adj = i1 == 1 && xq < x0 ? xq + period : xq
    w = x1 == x0 ? 0.0 : (xq_adj - x0) / (x1 - x0)
    return i0, i1, clamp(w, 0.0, 1.0)
end

@inline _lerp(a::Float64, b::Float64, w::Float64) = a + w * (b - a)

@inline function _trilerp(
    c000::Float64, c100::Float64, c010::Float64, c110::Float64,
    c001::Float64, c101::Float64, c011::Float64, c111::Float64,
    wa::Float64, wb::Float64, wc::Float64
)
    c00 = _lerp(c000, c100, wa)
    c10 = _lerp(c010, c110, wa)
    c01 = _lerp(c001, c101, wa)
    c11 = _lerp(c011, c111, wa)
    c0 = _lerp(c00, c10, wb)
    c1 = _lerp(c01, c11, wb)
    return _lerp(c0, c1, wc)
end

function lookup_payload_state(payload::Dict{String, Any}, alt_km::Float64, lat_deg::Float64, lon_deg::Float64)
    get(payload, "status", "error") == "ok" || error("Payload is not usable (status=$(get(payload, "status", "unknown"))).")

    grid = payload["grid"]
    fields = payload["fields"]

    alt_nodes = Vector{Float64}(grid["alt_km"])
    lat_nodes = Vector{Float64}(grid["lat_deg"])
    lon_nodes = Vector{Float64}(grid["lon_deg"])

    ia0, ia1, wa = _axis_segment(alt_nodes, alt_km)
    ilat0, ilat1, wb = _axis_segment(lat_nodes, lat_deg)
    ilon0, ilon1, wc = _lon_segment(lon_nodes, lon_deg)

    rho_grid = fields["density_kgm3"]
    temp_grid = fields["temperature_K"]
    w_ew_grid = fields["wind_ew_ms"]
    w_ns_grid = fields["wind_ns_ms"]
    w_up_grid = fields["wind_up_ms"]

    rho = _trilerp(
        Float64(rho_grid[ia0, ilat0, ilon0]), Float64(rho_grid[ia1, ilat0, ilon0]),
        Float64(rho_grid[ia0, ilat1, ilon0]), Float64(rho_grid[ia1, ilat1, ilon0]),
        Float64(rho_grid[ia0, ilat0, ilon1]), Float64(rho_grid[ia1, ilat0, ilon1]),
        Float64(rho_grid[ia0, ilat1, ilon1]), Float64(rho_grid[ia1, ilat1, ilon1]),
        wa, wb, wc
    )
    temp = _trilerp(
        Float64(temp_grid[ia0, ilat0, ilon0]), Float64(temp_grid[ia1, ilat0, ilon0]),
        Float64(temp_grid[ia0, ilat1, ilon0]), Float64(temp_grid[ia1, ilat1, ilon0]),
        Float64(temp_grid[ia0, ilat0, ilon1]), Float64(temp_grid[ia1, ilat0, ilon1]),
        Float64(temp_grid[ia0, ilat1, ilon1]), Float64(temp_grid[ia1, ilat1, ilon1]),
        wa, wb, wc
    )
    w_ew = _trilerp(
        Float64(w_ew_grid[ia0, ilat0, ilon0]), Float64(w_ew_grid[ia1, ilat0, ilon0]),
        Float64(w_ew_grid[ia0, ilat1, ilon0]), Float64(w_ew_grid[ia1, ilat1, ilon0]),
        Float64(w_ew_grid[ia0, ilat0, ilon1]), Float64(w_ew_grid[ia1, ilat0, ilon1]),
        Float64(w_ew_grid[ia0, ilat1, ilon1]), Float64(w_ew_grid[ia1, ilat1, ilon1]),
        wa, wb, wc
    )
    w_ns = _trilerp(
        Float64(w_ns_grid[ia0, ilat0, ilon0]), Float64(w_ns_grid[ia1, ilat0, ilon0]),
        Float64(w_ns_grid[ia0, ilat1, ilon0]), Float64(w_ns_grid[ia1, ilat1, ilon0]),
        Float64(w_ns_grid[ia0, ilat0, ilon1]), Float64(w_ns_grid[ia1, ilat0, ilon1]),
        Float64(w_ns_grid[ia0, ilat1, ilon1]), Float64(w_ns_grid[ia1, ilat1, ilon1]),
        wa, wb, wc
    )
    w_up = _trilerp(
        Float64(w_up_grid[ia0, ilat0, ilon0]), Float64(w_up_grid[ia1, ilat0, ilon0]),
        Float64(w_up_grid[ia0, ilat1, ilon0]), Float64(w_up_grid[ia1, ilat1, ilon0]),
        Float64(w_up_grid[ia0, ilat0, ilon1]), Float64(w_up_grid[ia1, ilat0, ilon1]),
        Float64(w_up_grid[ia0, ilat1, ilon1]), Float64(w_up_grid[ia1, ilat1, ilon1]),
        wa, wb, wc
    )

    return rho, temp, w_ew, w_ns, w_up
end

function _evaluate_payload(track, payload::Dict{String, Any})
    n = length(track.t)
    rho = Vector{Float64}(undef, n)
    temp = Vector{Float64}(undef, n)
    wind = Vector{SVector{3, Float64}}(undef, n)

    elapsed_s = @elapsed begin
        @inbounds for k in 1:n
            r, t, we, wn, wu = lookup_payload_state(
                payload,
                track.alt[k] * 1e-3,
                rad2deg(track.lat[k]),
                rad2deg(track.lon[k])
            )
            rho[k] = r
            temp[k] = t
            wind[k] = SVector{3, Float64}(we, wn, wu)
        end
    end

    return (rho=rho, temp=temp, wind=wind, elapsed_s=elapsed_s)
end

function _summarize_case(scenario::String, mode::String, track, ref, eval)
    n = length(track.t)
    rho_abs = [abs(eval.rho[k] - ref.rho[k]) for k in 1:n]
    rho_rel = [rho_abs[k] / max(abs(ref.rho[k]), 1e-12) for k in 1:n]
    temp_abs = [abs(eval.temp[k] - ref.temp[k]) for k in 1:n]
    wind_abs = [norm(eval.wind[k] - ref.wind[k]) for k in 1:n]

    return (
        scenario=scenario,
        mode=mode,
        samples=n,
        t_span_s=track.t[end] - track.t[1],
        alt_min_km=minimum(track.alt) * 1e-3,
        alt_max_km=maximum(track.alt) * 1e-3,
        rho_abs_max=maximum(rho_abs),
        rho_abs_p95=_pctl(rho_abs, 0.95),
        rho_rel_max=maximum(rho_rel),
        rho_rel_p95=_pctl(rho_rel, 0.95),
        temp_abs_max=maximum(temp_abs),
        temp_abs_p95=_pctl(temp_abs, 0.95),
        wind_abs_max=maximum(wind_abs),
        wind_abs_p95=_pctl(wind_abs, 0.95),
        eval_s=eval.elapsed_s,
        speedup_vs_point_to_point=ref.elapsed_s / max(eval.elapsed_s, 1e-9)
    )
end

function _error_table(scenario::String, mode::String, track, ref, eval)
    n = length(track.t)
    rho_abs = [abs(eval.rho[k] - ref.rho[k]) for k in 1:n]
    rho_rel = [rho_abs[k] / max(abs(ref.rho[k]), 1e-12) for k in 1:n]
    temp_abs = [abs(eval.temp[k] - ref.temp[k]) for k in 1:n]
    wind_abs = [norm(eval.wind[k] - ref.wind[k]) for k in 1:n]

    return DataFrame(
        scenario=fill(scenario, n),
        mode=fill(mode, n),
        t_s=track.t,
        alt_m=track.alt,
        lat_deg=rad2deg.(track.lat),
        lon_deg=rad2deg.(track.lon),
        rho_ref=ref.rho,
        rho_eval=eval.rho,
        rho_abs_err=rho_abs,
        rho_rel_err=rho_rel,
        temp_ref=ref.temp,
        temp_eval=eval.temp,
        temp_abs_err=temp_abs,
        wind_ref=[norm(w) for w in ref.wind],
        wind_eval=[norm(w) for w in eval.wind],
        wind_abs_err=wind_abs
    )
end

function _scenario_defs(samples::Int)
    return [
        (
            name="drag_passage",
            rp_alt_m=110e3,
            ra_alt_m=4500e3,
            i_deg=25.0,
            Ω_deg=20.0,
            ω_deg=35.0,
            ν_start_deg=-20.0,
            ν_end_deg=30.0,
            n_samples=samples
        ),
        (
            name="entry",
            rp_alt_m=80e3,
            ra_alt_m=1500e3,
            i_deg=22.0,
            Ω_deg=35.0,
            ω_deg=10.0,
            ν_start_deg=-45.0,
            ν_end_deg=-2.0,
            n_samples=samples
        ),
        (
            name="orbit",
            rp_alt_m=450e3,
            ra_alt_m=450e3,
            i_deg=30.0,
            Ω_deg=15.0,
            ω_deg=0.0,
            ν_start_deg=0.0,
            ν_end_deg=360.0,
            n_samples=samples
        )
    ]
end

function _parse_scenarios(raw::String)
    txt = lowercase(strip(raw))
    if isempty(txt) || txt == "all"
        return Set(["drag_passage", "entry", "orbit"])
    end
    chosen = Set{String}()
    for token in split(txt, ",")
        v = lowercase(strip(token))
        v in ("drag_passage", "entry", "orbit") || error("Unsupported scenario '$v'.")
        push!(chosen, v)
    end
    return chosen
end

function run_study()
    opts = parse_cli(copy(ARGS))

    planet_name = lowercase(get(opts, "planet", "mars"))
    samples = parse(Int, get(opts, "samples", "220"))
    selected = _parse_scenarios(get(opts, "scenarios", "all"))

    static_root = get(opts, "static-root", joinpath(REPO_ROOT, "GRAMSuite.jl/GRAM Suite 2.0", "simulation", "GRAM", "static_grids"))
    grid_file = get(opts, "grid-file", joinpath(static_root, "$(planet_name)_grid.jls"))
    surrogate_file = get(opts, "surrogate-file", joinpath(static_root, "surrogates", "$(planet_name)_surrogate.jls"))

    out_summary = get(opts, "out-summary", joinpath(REPO_ROOT, "output", "gram_offline_db_accuracy_summary.csv"))
    out_errors = get(opts, "out-errors", joinpath(REPO_ROOT, "output", "gram_offline_db_accuracy_errors.csv"))

    spice_path = joinpath(REPO_ROOT, "GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
    planet = _planet_from_name(planet_name, spice_path)
    planet.L_PI .= [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]

    initial_time = InitialTime(year=2001, month=11, day=6, hour=19, minute=0, second=32.0)
    density_model = Base.invokelatest(EM.GRAMAtmosphereModel; planet_name=planet_name, initial_time=initial_time)
    args = _make_args(planet, density_model, initial_time)
    p = ODEParams{1}(args=args)

    grid_payload = isfile(grid_file) ? load_payload(grid_file) : nothing
    surrogate_payload = isfile(surrogate_file) ? load_payload(surrogate_file) : nothing

    if grid_payload === nothing
        @warn "Grid payload not found; grid comparison will be skipped." file=grid_file
    end
    if surrogate_payload === nothing
        @warn "Surrogate payload not found; surrogate comparison will be skipped." file=surrogate_file
    end

    rows = NamedTuple[]
    error_tables = DataFrame[]

    println("GRAM offline DB accuracy study")
    println(@sprintf("planet=%s samples=%d", planet_name, samples))
    println("grid_file=$(grid_file)")
    println("surrogate_file=$(surrogate_file)")

    for sc in _scenario_defs(samples)
        sc.name in selected || continue
        println("Running scenario: $(sc.name)")

        track = _generate_track(
            planet;
            rp_alt_m=sc.rp_alt_m,
            ra_alt_m=sc.ra_alt_m,
            i_deg=sc.i_deg,
            Ω_deg=sc.Ω_deg,
            ω_deg=sc.ω_deg,
            ν_start_deg=sc.ν_start_deg,
            ν_end_deg=sc.ν_end_deg,
            n_samples=sc.n_samples
        )

        direct = _evaluate_direct(track, density_model, p)
        push!(rows, merge(_summarize_case(sc.name, "point_to_point", track, direct, direct), (source="gram_live",)))

        if grid_payload !== nothing && get(grid_payload, "status", "error") == "ok"
            grid_eval = _evaluate_payload(track, grid_payload)
            push!(rows, merge(_summarize_case(sc.name, "offline_grid", track, direct, grid_eval), (source=basename(grid_file),)))
            push!(error_tables, _error_table(sc.name, "offline_grid", track, direct, grid_eval))
        end

        if surrogate_payload !== nothing && get(surrogate_payload, "status", "error") == "ok"
            s_eval = _evaluate_payload(track, surrogate_payload)
            push!(rows, merge(_summarize_case(sc.name, "offline_surrogate", track, direct, s_eval), (source=basename(surrogate_file),)))
            push!(error_tables, _error_table(sc.name, "offline_surrogate", track, direct, s_eval))
        end
    end

    summary_df = DataFrame(rows)
    errors_df = isempty(error_tables) ? DataFrame() : vcat(error_tables...)

    mkpath(dirname(out_summary))
    mkpath(dirname(out_errors))
    CSV.write(out_summary, summary_df)
    CSV.write(out_errors, errors_df)

    println("\nSummary:")
    show(summary_df, allrows=true, allcols=true)
    println("\n\nSaved CSVs:")
    println("  summary: $out_summary")
    println("  errors : $out_errors")

    return (summary=summary_df, errors=errors_df)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_study()
end
