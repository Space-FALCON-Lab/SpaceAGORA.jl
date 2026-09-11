# Differential-correction fit of a time-aligned scenario's Cartesian initial
# condition against its telemetry, via the public manifest/run machinery.
#
# Method (the CYGNSS July-2026 campaign recipe): one baseline run plus six
# perturbation runs (+pos_step_m on each position axis, +vel_step_mps on each
# velocity axis, applied through the scenario's ic_offset_m/ic_offset_mps
# keys) build a finite-difference sensitivity matrix of the per-sample
# position residuals; a linear least-squares solve yields the 6-state IC
# correction, which is then verified by a validation propagation. Fitted
# state-vector components trade against un-modeled forces (an SMA offset
# absorbs along-track drag drift), so interpret the corrected drag/IC pair
# jointly — see the campaign's energy-consistency methodology.

const _IC_FIT_EVENTS = ("state_x_time", "state_y_time", "state_z_time")

function _ic_fit_series(errors_csv::String, scenario_name::String)
    df = CSV.read(errors_csv, DataFrame)
    df = df[(df.scenario .== scenario_name) .& in.(df.event, Ref(collect(_IC_FIT_EVENTS))), :]
    nrow(df) > 0 || throw(ArgumentError("fit_initial_state: no $( _IC_FIT_EVENTS ) rows for scenario '$scenario_name' in $errors_csv"))
    out = Dict{String, Dict{Float64, NTuple{2, Float64}}}()
    for ev in _IC_FIT_EVENTS
        sub = df[df.event .== ev, :]
        out[ev] = Dict(
            Float64(sub.telemetry_axis[i]) => (Float64(sub.sim_interp_value_km[i]), Float64(sub.error_km[i]))
            for i in 1:nrow(sub)
        )
    end
    return out
end

function _ic_fit_run(
    manifest::Dict{String, Any},
    scenario_index::Int,
    offsets_m::NTuple{3, Float64},
    offsets_mps::NTuple{3, Float64},
    label::String,
    workdir::String,
    profile::Symbol
)
    m = deepcopy(manifest)
    scen = m["scenarios"][scenario_index]
    scen["ic_offset_m"] = Any[offsets_m...]
    scen["ic_offset_mps"] = Any[offsets_mps...]
    manifest_path = joinpath(workdir, "icfit_manifest_$(label).toml")
    open(manifest_path, "w") do io
        TOML.print(io, m)
    end
    req = VerificationRequest(
        profile=profile,
        out_summary=joinpath(workdir, "icfit_summary_$(label).csv"),
        out_errors=joinpath(workdir, "icfit_errors_$(label).csv"),
        manifest_path=manifest_path,
        enforce=false,
        generate_plots=false
    )
    result = run_verification(req)
    return (
        series=_ic_fit_series(req.out_errors, String(scen["name"])),
        summary=result.summary
    )
end

@inline function _ic_fit_mean_axis_rmse(summary::DataFrame, scenario_name::String)::Float64
    rows = summary[(summary.scenario .== scenario_name) .& in.(summary.event, Ref(collect(_IC_FIT_EVENTS))), :]
    nrow(rows) == 3 || throw(ArgumentError("fit_initial_state: expected 3 position events in summary for '$scenario_name'."))
    return sum(Float64.(rows.rmse_km)) / 3.0
end

"""
    fit_initial_state(manifest_path, scenario_name;
                      profile=:quick, pos_step_m=100.0, vel_step_mps=0.005,
                      workdir=mktempdir(), validate=true)

Differential-correction fit of the named time-aligned scenario's Cartesian
initial condition: 1 baseline + 6 perturbation runs (finite-difference
sensitivities via the scenario's `ic_offset_m`/`ic_offset_mps` keys), a
linear least-squares solve over the stacked x/y/z per-sample residuals, and
(when `validate=true`) a validation propagation at the corrected offsets.

Perturbations are applied on top of any offsets already present in the
manifest, and the returned `offsets_total_*` are ready to paste back into it.
Returns a NamedTuple with `correction_m`, `correction_mps`,
`offsets_total_m`, `offsets_total_mps`, `rmse_before_km`,
`rmse_predicted_km` (linear model), and `rmse_after_km` (validated;
`NaN` when `validate=false`).
"""
function fit_initial_state(
    manifest_path::String,
    scenario_name::String;
    profile::Symbol=:quick,
    pos_step_m::Float64=100.0,
    vel_step_mps::Float64=0.005,
    workdir::String=mktempdir(),
    validate::Bool=true
)
    pos_step_m > 0.0 && vel_step_mps > 0.0 || throw(ArgumentError(
        "fit_initial_state: perturbation steps must be > 0."
    ))
    mkpath(workdir)
    manifest = TOML.parsefile(manifest_path)
    haskey(manifest, "scenarios") || throw(ArgumentError("fit_initial_state: no scenarios in $manifest_path"))
    idx = findfirst(s -> String(get(s, "name", "")) == scenario_name, manifest["scenarios"])
    idx === nothing && throw(ArgumentError("fit_initial_state: scenario '$scenario_name' not found in $manifest_path"))
    scen = manifest["scenarios"][idx]
    for key in ("x_ic", "y_ic", "z_ic", "vx_ic", "vy_ic", "vz_ic")
        haskey(get(scen, "telemetry_columns", Dict()), key) || throw(ArgumentError(
            "fit_initial_state: scenario '$scenario_name' must provide Cartesian IC telemetry columns (missing '$key')."
        ))
    end
    base_off_m = Tuple(Float64.(get(scen, "ic_offset_m", [0.0, 0.0, 0.0])))
    base_off_mps = Tuple(Float64.(get(scen, "ic_offset_mps", [0.0, 0.0, 0.0])))

    base = _ic_fit_run(manifest, idx, base_off_m, base_off_mps, "base", workdir, profile)
    axes_common = intersect(
        (Set(keys(base.series[ev])) for ev in _IC_FIT_EVENTS)...
    )

    perturbations = [
        ((pos_step_m, 0.0, 0.0), (0.0, 0.0, 0.0), pos_step_m),
        ((0.0, pos_step_m, 0.0), (0.0, 0.0, 0.0), pos_step_m),
        ((0.0, 0.0, pos_step_m), (0.0, 0.0, 0.0), pos_step_m),
        ((0.0, 0.0, 0.0), (vel_step_mps, 0.0, 0.0), vel_step_mps),
        ((0.0, 0.0, 0.0), (0.0, vel_step_mps, 0.0), vel_step_mps),
        ((0.0, 0.0, 0.0), (0.0, 0.0, vel_step_mps), vel_step_mps),
    ]
    pert_series = Vector{Any}(undef, 6)
    for (k, (dm, dv, _)) in enumerate(perturbations)
        run = _ic_fit_run(
            manifest, idx,
            base_off_m .+ dm, base_off_mps .+ dv,
            "p$(k)", workdir, profile
        )
        pert_series[k] = run.series
        for ev in _IC_FIT_EVENTS
            intersect!(axes_common, Set(keys(run.series[ev])))
        end
    end
    n_axes = length(axes_common)
    n_axes >= 12 || throw(ArgumentError(
        "fit_initial_state: only $n_axes common comparison samples across runs; need >= 12."
    ))
    axes_sorted = sort!(collect(axes_common))

    n_rows = 3 * n_axes
    b = Vector{Float64}(undef, n_rows)
    J = Matrix{Float64}(undef, n_rows, 6)
    for (ei, ev) in enumerate(_IC_FIT_EVENTS)
        base_ev = base.series[ev]
        for (ai, ax) in enumerate(axes_sorted)
            row = (ei - 1) * n_axes + ai
            sim_base, err_base = base_ev[ax]
            b[row] = -err_base * 1.0e3                       # residual to cancel [m]
            for k in 1:6
                step = perturbations[k][3]
                sim_pert = pert_series[k][ev][ax][1]
                J[row, k] = (sim_pert - sim_base) * 1.0e3 / step
            end
        end
    end
    delta = J \ b
    predicted_resid = b .- J * delta
    rmse_before = _ic_fit_mean_axis_rmse(base.summary, scenario_name)
    rmse_predicted = sqrt(sum(abs2, predicted_resid) / n_rows) * 1.0e-3

    correction_m = (delta[1], delta[2], delta[3])
    correction_mps = (delta[4], delta[5], delta[6])
    total_m = base_off_m .+ correction_m
    total_mps = base_off_mps .+ correction_mps

    rmse_after = NaN
    if validate
        validated = _ic_fit_run(manifest, idx, total_m, total_mps, "validated", workdir, profile)
        rmse_after = _ic_fit_mean_axis_rmse(validated.summary, scenario_name)
    end
    println("fit_initial_state[$scenario_name]: correction " *
            "d_r=$(round.(correction_m, digits=3)) m, d_v=$(round.(correction_mps .* 1.0e3, digits=4)) mm/s; " *
            "mean-axis RMSE $(round(rmse_before, digits=5)) -> " *
            "$(validate ? round(rmse_after, digits=5) : "(not validated)") km " *
            "(linear prediction $(round(rmse_predicted, digits=5)))")
    return (
        correction_m=correction_m,
        correction_mps=correction_mps,
        offsets_total_m=total_m,
        offsets_total_mps=total_mps,
        rmse_before_km=rmse_before,
        rmse_predicted_km=rmse_predicted,
        rmse_after_km=rmse_after,
        n_samples=n_axes,
        workdir=workdir
    )
end
