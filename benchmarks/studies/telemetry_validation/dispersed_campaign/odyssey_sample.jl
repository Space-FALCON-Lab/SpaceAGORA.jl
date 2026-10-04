# One member of the dispersed Odyssey MarsGRAM campaign. Included once per
# process (the campaign coordinator, and each process-pool worker on its first
# sample) by run_campaign.jl's sample closure.
#
# Thread-safe by construction: no withenv, no cd. Every run-wide setting is
# already in ENV (set by run_campaign.jl before the pool spawns, so workers
# inherit it); the only per-member inputs are the GRAM seed, whether the member
# is perturbed, and its own absolute results directory.
module OdysseyDispersedSample

using SpaceAGORA
using SpaceAGORA.TelemetryVerification
using DataFrames
using CSV
using TOML
using Statistics

const TV = SpaceAGORA.TelemetryVerification
const SE = SpaceAGORA.SimulationEngine
const SCB = SpaceAGORA.SimulationModel.SimulationCallbacks

const MANIFEST = normpath(joinpath(@__DIR__, "..", "manifests", "odyssey_marsgram.toml"))

function member_cfg(seed::Int, perturbed::Bool)
    doc = TOML.parsefile(MANIFEST)
    truth = doc["scenarios"][1]["atmosphere_truth"]
    if perturbed
        truth["gram_seed"] = seed
        truth["gram_perturbation_scales"] = [1.0, 0.0, 0.0, 0.0]
    end
    path = joinpath(mktempdir(), "manifest.toml")
    open(io -> TOML.print(io, doc), path, "w")
    return only(TV._load_scenarios_from_manifest(path))
end

save_fields() = vcat(TV._save_fields_for_study(), [
    SCB.SaveField(:heat_rate, (u, t, integ) -> SCB._save_heat_rate(1, u, t, integ); per_satellite=true),
    SCB.SaveField(:heat_load, (u, t, integ) -> SCB._save_heat_load(1, u, t, integ); per_satellite=true),
])

# Value at time tq from a sorted time series, nearest saved sample at or before tq.
at_time(t, y, tq) = y[clamp(searchsortedlast(t, tq), 1, length(t))]

"""
    run_member(seed, perturbed, orbits, out_dir) -> NamedTuple

Run one member and write its per-orbit and per-pass CSVs into `out_dir`.
"""
function run_member(seed::Int, perturbed::Bool, orbits::Int, out_dir::String)
    mkpath(out_dir)
    # A rerun into the same directory must not read an earlier attempt's files.
    for f in ("simulation_results.csv", "per_orbit.csv", "per_pass_r.csv", "gram_density_perturbation_log.csv")
        rm(joinpath(out_dir, f); force=true)
    end
    cfg = member_cfg(seed, perturbed)
    args = TV._with_study_settings(TV._make_orbit_args(cfg, orbits); quick=false)
    args = TV._with_configuration(args;
        simulation_settings=SpaceAGORA.SimulationModel.SimulationSettings(
            results=true, verbose=false, results_directory=out_dir, generate_plots=false,
            generate_filenames=false, normalize=false, save_csv=true),
        solver_config=nothing)
    result = nothing; err = ""
    solve_s = @elapsed try
        result = SE.run_simulation(args; isolate_state=false, save_fields=save_fields(),
                                   return_solution=true, return_solver_metadata=true,
                                   extra_callbacks=TV._scenario_extra_callbacks(cfg))
    catch e
        err = sprint(showerror, e)
    end

    csv = joinpath(out_dir, "simulation_results.csv")
    n_apo = 0; n_peri = 0; apo_first = NaN; apo_last = NaN; peri_min = NaN; heat_final = NaN; final_t = NaN
    if isfile(csv)
        df = CSV.read(csv, DataFrame)
        t = Float64.(df.time); final_t = t[end]
        ex = TV._extract_extrema_series(df, args.environment_model.planet, cfg.orbit_altitude_mode)
        hl_col = only(filter(c -> occursin("heat_load", c), names(df)))
        hr_col = only(filter(c -> occursin("heat_rate", c), names(df)))
        hl = Float64.(df[!, hl_col]); hr = Float64.(df[!, hr_col])
        apo_t = ex.apo.time_s; n_apo = length(ex.apo.altitude); n_peri = length(ex.peri.altitude)
        # Every apoapsis is kept: the decay over n passes needs the apoapsis after
        # the last of them, so the apoapses are not truncated to the periapsis
        # count. Rows run to the larger count; the shorter series is NaN-padded.
        # Pass k (periapsis k) lies between apoapsis k and apoapsis k+1, or the end
        # of the run when no later apoapsis was reached.
        n_rows = max(n_apo, n_peri)
        heat_pass = fill(NaN, n_rows); peak_rate = fill(NaN, n_rows)
        for k in 1:min(n_peri, n_apo)
            a = apo_t[k]; b = k + 1 <= n_apo ? apo_t[k + 1] : t[end]
            heat_pass[k] = at_time(t, hl, b) - at_time(t, hl, a)
            peak_rate[k] = maximum(hr[(t .>= a) .& (t .<= b)]; init=0.0)
        end
        pad(v) = vcat(Float64.(v), fill(NaN, n_rows - length(v)))
        off = something(cfg.epoch_orbit_offset, 0.0)
        CSV.write(joinpath(out_dir, "per_orbit.csv"), DataFrame(
            index=1:n_rows, flight_orbit=off .+ (0:n_rows-1),
            apo_time_s=pad(apo_t), apo_km=pad(ex.apo.altitude),
            peri_time_s=pad(ex.peri.time_s), peri_km=pad(ex.peri.altitude),
            pass_heat_load_Jcm2=heat_pass, pass_peak_heat_rate_Wcm2=peak_rate))
        n_apo > 0 && (apo_first = ex.apo.altitude[1]; apo_last = ex.apo.altitude[n_apo])
        n_peri > 0 && (peri_min = minimum(ex.peri.altitude))
        heat_final = hl[end]
        rm(csv)   # per-orbit series kept; the full state history is not needed
    end

    log = joinpath(out_dir, "gram_density_perturbation_log.csv")
    rdw_mean = NaN
    if isfile(log)
        lg = CSV.read(log, DataFrame)
        k = lg[lg.kind .== 2, :]   # B knots: walk samples with GRAM mean and sigma
        rows = map(collect(groupby(k, :pass))) do g
            w = g.mean_density_kgm3; low = g[g.alt_m .< 130e3, :]
            (pass=g.pass[1], knots=nrow(g),
             r_density_weighted=sum(w .* g.r) / sum(w),
             sigma_density_weighted=sum(w .* g.sigma_frac) / sum(w),
             r_std=std(g.r), r_std_below130=nrow(low) > 1 ? std(low.r) : NaN,
             sigma_median_below130=nrow(low) > 0 ? median(low.sigma_frac) : NaN)
        end
        pr = DataFrame(rows)
        CSV.write(joinpath(out_dir, "per_pass_r.csv"), pr)
        rdw_mean = mean(pr.r_density_weighted)
    end

    stats = result === nothing ? nothing : result.solution.stats
    retcode = result === nothing ? "ERROR" : string(result.solution.retcode)
    # The orbit-count stop and an impact both return Terminated; the run's own
    # orbit counter (one count per apoapsis crossing) tells them apart.
    completed = -1
    if result !== nothing
        try
            completed = Int(result.solution.prob.p.orbit_counter[1]) - 1
        catch
        end
    end
    cause = result === nothing ? "error" :
            retcode == "Success" ? "end_of_time_span" :
            retcode != "Terminated" ? "retcode_$(retcode)" :
            completed < 0 ? "terminated_unknown" :
            completed >= orbits ? "orbit_count" : "terminated_before_orbit_count"
    return (seed=seed, perturbed=perturbed, pid=getpid(), host=gethostname(),
            solve_s=solve_s,
            retcode=retcode, termination_cause=cause, completed_orbits=completed,
            error=err, final_t_s=final_t, requested_orbits=orbits, n_apo=n_apo, n_peri=n_peri, n_passes=n_peri,
            apo_first_km=apo_first, apo_final_km=apo_last, apo_decay_km=apo_first - apo_last,
            peri_min_km=peri_min, heat_load_final_Jcm2=heat_final, r_density_weighted_mean=rdw_mean,
            naccept=stats === nothing ? -1 : Int(stats.naccept),
            nreject=stats === nothing ? -1 : Int(stats.nreject),
            nf=stats === nothing ? -1 : Int(stats.nf))
end

end
