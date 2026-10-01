# Full-arc cross-validation of the 4-body x 3-gravity x third-body-off/on parity
# matrix against one reference tool, keeping the per-epoch position errors.
#
# Usage (from the repository root):
#   julia --startup-file=no --project=. scripts/xval_fullarc.jl <target> <outdir> [variant]
#
#   target   gmat      data/telemetry/GMAT_Examples (real GMAT R2025a runs)
#            stk       data/telemetry/stk_results   (STK HPOP runs)
#            basilisk  Basilisk_Examples_Full; directory taken from XVAL_BASILISK_DIR
#   variant  committed (default): the force-model overrides in test/gmat_scenario_matrix.jl
#            as_basilisk: gmat reference, but the :basilisk overrides (sensitivity only)
#            earth_egm96: stk reference, Earth J2/J50 with EGM96_GMAT_L50.csv instead of GGM05C
#            moon_file_c20: stk reference, Moon J2 with the unmodified LP165P.csv C(2,0)
#
# Optional: XVAL_SCENARIOS=earth_j0_tbfalse,... restricts the case list.
#
# Every reference sample is compared (max_points raised past the file length), so
# the RMS covers the whole arc the reference file holds. Per-case output:
#   <outdir>/<target>_<variant>/<scenario>/manifest.toml   scenario as run
#   <outdir>/<target>_<variant>/<scenario>/summary.csv     TelemetryVerification summary
#   <outdir>/<target>_<variant>/<scenario>/series.arrow    t_s, reference and SpaceAGORA
#                                                          x/y/z (km), dx/dy/dz and |dr| (m)
# and one row per case in <outdir>/<target>_<variant>/results.csv.

using Arrow
using CSV
using DataFrames
using Dates
using Printf
using Statistics
using TOML

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

# Load only the definitions of the harness (everything before its first
# top-level testset guard), not the testsets themselves.
let
    path = joinpath(REPO_ROOT, "test", "gmat_scenario_matrix.jl")
    src = readlines(path)
    boundary = findfirst(l -> occursin("SPACEAGORA_SKIP_GMAT_MATRIX", l), src)
    boundary === nothing && error("boundary marker not found in $path")
    include_string(Main, join(src[1:boundary-1], "\n"), path)
end

const ALL_SCENARIOS = [
    "$(b)_$(g)_$(tb)" for b in ("earth", "mars", "venus", "moon")
    for g in ("j0", "j2", "j50") for tb in ("tbfalse", "tbtrue")
]

function _resolver(target::String)
    if target == "gmat"
        return _scenario_gmat_path
    elseif target == "stk"
        return _scenario_stk_path
    elseif target == "basilisk"
        dir = get(ENV, "XVAL_BASILISK_DIR", "")
        isdir(dir) || error("target basilisk needs XVAL_BASILISK_DIR pointing at Basilisk_Examples_Full")
        return name -> begin
            p = split(name, "_")
            body = p[1] == "moon" ? "Luna" : uppercasefirst(p[1])
            joinpath(dir, "Sim_$(body)_1M_$(uppercase(p[2]))_$(p[3] == "tbtrue" ? "TBTrue" : "TBFalse").feather")
        end
    end
    error("unknown target $target")
end

function _overrides(name::String, target::String, variant::String)
    planet, gtag, _ = split(name, "_")
    if variant == "committed"
        return _matrix_scenario_overrides(name, Symbol(target))
    elseif variant == "as_basilisk"
        target == "gmat" || error("variant as_basilisk applies to the gmat target only")
        return _matrix_scenario_overrides(name, :basilisk)
    elseif variant == "earth_egm96"
        target == "stk" || error("variant earth_egm96 applies to the stk target only")
        ov = _matrix_scenario_overrides(name, :stk)
        if planet == "earth" && gtag != "j0"
            ov["gravity_harmonics_file"] = _GMAT_HARMONICS_EARTH_EGM96_FILE
        end
        return ov
    elseif variant == "moon_file_c20"
        target == "stk" || error("variant moon_file_c20 applies to the stk target only")
        ov = _matrix_scenario_overrides(name, :stk)
        if planet == "moon" && gtag == "j2"
            ov["gravity_harmonics_file"] = _GMAT_HARMONICS_MOON_FILE
        end
        return ov
    end
    error("unknown variant $variant")
end

function _axis_rows(errors::DataFrame, name::String, event::String)
    rows = errors[(errors.scenario .== name) .& (errors.event .== event), :]
    return sort(rows, :idx)
end

function run_case(name::String, target::String, variant::String, ref_path::String, casedir::String)
    mkpath(casedir)
    ic = _matrix_initial_conditions(name)
    traj = _build_time_aligned_reference(
        ref_path, _scenario_planet_name(name), casedir, name;
        sma_km=ic.sma_km, ecc=ic.ecc, inc_deg=ic.inc_deg,
        aop_deg=ic.aop_deg, raan_deg=ic.raan_deg, ta_deg=ic.ta_deg
    )
    scenario = _base_scenario_dict(name, traj.telemetry_path)
    merge!(scenario, _overrides(name, target, variant))
    # Compare every reference sample: no truncation to the first N rows.
    scenario["max_points_quick"] = 1_000_000_000
    scenario["max_points_full"] = 1_000_000_000

    manifest_path = joinpath(casedir, "manifest.toml")
    open(manifest_path, "w") do io
        TOML.print(io, Dict{String, Any}("scenarios" => Any[scenario]))
    end
    tv_errors_path = joinpath(casedir, "errors_tv.csv")
    req = TV.VerificationRequest(
        profile=:quick,
        out_summary=joinpath(casedir, "summary.csv"),
        out_errors=tv_errors_path,
        manifest_path=manifest_path,
        enforce=false,
        generate_plots=false
    )
    env = Pair{String, String}[
        "SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH" => _gmat_planetary_kernel_relpath(),
        pairs(_telemetry_solver_env_overrides())...
    ]
    t_wall = @elapsed result = withenv(env...) do
        TV.run_verification(req)
    end

    xr = _axis_rows(result.errors, name, "state_x_time")
    yr = _axis_rows(result.errors, name, "state_y_time")
    zr = _axis_rows(result.errors, name, "state_z_time")
    t = Float64.(xr.telemetry_axis)
    (t == Float64.(yr.telemetry_axis) && t == Float64.(zr.telemetry_axis)) || error("axis epochs differ for $name")
    dx = Float64.(xr.error_km) .* 1e3
    dy = Float64.(yr.error_km) .* 1e3
    dz = Float64.(zr.error_km) .* 1e3
    dr = sqrt.(dx .^ 2 .+ dy .^ 2 .+ dz .^ 2)
    series = DataFrame(
        t_s=t,
        ref_x_km=Float64.(xr.telemetry_value_km), ref_y_km=Float64.(yr.telemetry_value_km), ref_z_km=Float64.(zr.telemetry_value_km),
        sa_x_km=Float64.(xr.sim_interp_value_km), sa_y_km=Float64.(yr.sim_interp_value_km), sa_z_km=Float64.(zr.sim_interp_value_km),
        dx_m=dx, dy_m=dy, dz_m=dz, err_m=dr
    )
    Arrow.write(joinpath(casedir, "series.arrow"), series)
    # The TV error table repeats the same data per axis and is large; the series
    # above holds it. The time-aligned reference copy is regenerated from ref_path.
    rm(tv_errors_path; force=true)
    rm(traj.telemetry_path; force=true)

    rms = sqrt(mean(dr .^ 2))
    # Cross-check against the harness's own per-axis RMS combination.
    s = result.summary
    axis_rmse = [Float64(s[(s.scenario .== name) .& (s.event .== ev), :rmse_km][1]) for ev in ("state_x_time", "state_y_time", "state_z_time")]
    rms_harness = sqrt(sum(axis_rmse .^ 2)) * 1e3
    isapprox(rms, rms_harness; rtol=1e-9) || error("RMS mismatch for $name: $rms vs $rms_harness")
    n10k = min(10_000, length(dr))
    retcode = "solver_retcode" in names(s) ? String(string(s[s.scenario .== name, :solver_retcode][1])) : ""
    ov = TOML.parsefile(manifest_path)["scenarios"][1]
    return (
        scenario=name,
        target=target,
        variant=variant,
        rms_m=rms,
        max_m=maximum(dr),
        t_at_max_s=t[argmax(dr)],
        end_err_m=dr[end],
        n_points=length(dr),
        arc_start_s=t[1],
        arc_end_s=t[end],
        rms_first10k_m=sqrt(mean(dr[1:n10k] .^ 2)),
        gravity_degree=ov["gravity_harmonics_degree"],
        gravity_order=ov["gravity_harmonics_order"],
        gravity_file=String(ov["gravity_harmonics_file"]),
        gm_override_m3s2=Float64(get(ov, "gravity_harmonics_gm_override_m3s2", NaN)),
        nbody_bodies=join(String.(ov["nbody_bodies"]), "+"),
        solver_retcode=retcode,
        wall_s=t_wall,
        reference_path=ref_path
    )
end

function main()
    length(ARGS) >= 2 || error("usage: scripts/xval_fullarc.jl <gmat|stk|basilisk> <outdir> [variant]")
    target = ARGS[1]
    outdir = abspath(ARGS[2])
    variant = length(ARGS) >= 3 ? ARGS[3] : "committed"
    selected = strip(get(ENV, "XVAL_SCENARIOS", ""))
    scenarios = isempty(selected) ? ALL_SCENARIOS : String.(strip.(split(selected, ",")))
    resolver = _resolver(target)
    rundir = joinpath(outdir, "$(target)_$(variant)")
    mkpath(rundir)
    commit = try readchomp(`git -C $REPO_ROOT rev-parse HEAD`) catch; "unknown" end
    dirty = try !isempty(readchomp(`git -C $REPO_ROOT status --porcelain --untracked-files=no -- src test scripts data/Gravity_harmonics_data`)) catch; true end
    println("xval_fullarc target=$target variant=$variant commit=$commit src/test/scripts dirty=$dirty")
    println("started $(now())  julia $(VERSION)  threads=$(Threads.nthreads())")
    rows = NamedTuple[]
    t_total = @elapsed for name in scenarios
        ref_path = resolver(name)
        isfile(ref_path) || error("missing reference for $name: $ref_path")
        row = run_case(name, target, variant, ref_path, joinpath(rundir, name))
        push!(rows, row)
        @printf("%-18s rms=%14.6f m  max=%14.6f m  n=%7d  arc_end=%11.3f s  wall=%7.1f s\n",
            name, row.rms_m, row.max_m, row.n_points, row.arc_end_s, row.wall_s)
        CSV.write(joinpath(rundir, "results.csv"), DataFrame(rows))
    end
    open(joinpath(rundir, "run_info.toml"), "w") do io
        TOML.print(io, Dict(
            "target" => target, "variant" => variant, "commit" => commit,
            "src_test_scripts_dirty" => dirty, "finished" => string(now()),
            "total_wall_s" => t_total, "julia" => string(VERSION),
            "hostname" => gethostname(), "cpu" => Sys.cpu_info()[1].model
        ))
    end
    @printf("done: %d cases in %.1f s\n", length(rows), t_total)
end

main()
