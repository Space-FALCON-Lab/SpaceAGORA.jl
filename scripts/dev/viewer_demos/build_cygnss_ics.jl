#!/usr/bin/env julia
#
# Build the CYGNSS constellation's flown state at 2025-06-06T00:00:00Z, the
# start of the 96-hour telemetry window, together with the flight track of all
# seven spacecraft across the window.
#
#   python3 scripts/dev/viewer_demos/fetch_cygnss_l1_states.py   # once, network
#   julia --project=. scripts/dev/run.jl viewer_demos/build_cygnss_ics.jl
#
# Outputs, all under the gitignored data/telemetry/CYGNSS/:
#
#   constellation_ics_20250606.json            seven flown states at the epoch
#   constellation_ics_20250606_catalog.json  the element-set states it replaced
#   catalog_vs_flown_20250606.json           the measured difference between them
#   cygnss_constellation_tracks_20250606_96hr.feather
#                                              all seven spacecraft, 1 Hz, J2000
#
# That directory is gitignored, which is why this script is committed: it is the
# auditable half of the files it produces.
#
# WHERE EACH STATE COMES FROM
#
#   all seven   flight. NASA CYGNSS Level 1 Science Data Record version 3.2,
#               the same product the reconstruction record uses, fetched by
#               variable subsetting over OPeNDAP (see the Python script above).
#               Earth-fixed in the product; converted here to J2000.
#   FM06        absent. CelesTrak's SATCAT records NORAD 41889 decayed on
#               2024-06-13, a year before this epoch, so the constellation at
#               this epoch is seven spacecraft and not eight. NASA's CMR lists
#               Level 1 granules for exactly the other seven on these dates.
#
# The historical element sets in cygnss_historical.tle are no longer the source
# of any state. They are kept, and still propagated here, for one reason: five
# of these seven spacecraft used to be catalog-derived, and now that every one
# of them has a flight state the catalog error can be measured directly instead
# of inferred from the two spacecraft that happened to have telemetry. That
# measurement is the point of `catalog_vs_flown_20250606.json`.

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

using Printf
using LinearAlgebra
using StaticArrays
using Statistics
using Dates
using JSON
using Arrow
using DataFrames
using SPICE
using SatelliteToolboxTle
using SatelliteToolboxSgp4
using SatelliteToolboxTransformations

include(joinpath(@__DIR__, "cygnss_ics.jl"))
include(joinpath(@__DIR__, "cygnss_tracks.jl"))
using .CygnssICs
using .CygnssTracks

const SPICE_PATH = get(ENV, "SPACEAGORA_SPICE_PATH", joinpath(REPO_ROOT, "data", "GRAMSuite.jl", "GRAM Suite 2.0", "SPICE"))
const TELEM_DIR = get(ENV, "SPACEAGORA_CYGNSS_DATA", joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS"))
const TLE_FILE = joinpath(@__DIR__, "cygnss_historical.tle")

const ECEF_FILE = joinpath(TELEM_DIR, "cygnss_l1_pvt_ecef_20250606_96hr.feather")
const OUT_FILE = joinpath(TELEM_DIR, "constellation_ics_20250606.json")
const CATALOG_FILE = joinpath(TELEM_DIR, "constellation_ics_20250606_catalog.json")
const COMPARISON_FILE = joinpath(TELEM_DIR, "catalog_vs_flown_20250606.json")
const TRACKS_FILE = joinpath(TELEM_DIR, "cygnss_constellation_tracks_20250606_96hr.feather")

# The eight flight models as launched, with the NORAD ids CelesTrak's SATCAT
# gives for international designator 2016-078 (retrieved 2026-09-16 from
# https://celestrak.org/satcat/records.php?INTDES=2016-078&FORMAT=json), and the
# `spacecraft_num` each carries in its Level 1 granules.
const ROSTER = [
    (fm = 1, name = "CYGNSS FM01", norad = 41887, catalog = "CYGFM01", intldes = "2016-078D"),
    (fm = 2, name = "CYGNSS FM02", norad = 41886, catalog = "CYGFM02", intldes = "2016-078C"),
    (fm = 3, name = "CYGNSS FM03", norad = 41891, catalog = "CYGFM03", intldes = "2016-078H"),
    (fm = 4, name = "CYGNSS FM04", norad = 41885, catalog = "CYGFM04", intldes = "2016-078B"),
    (fm = 5, name = "CYGNSS FM05", norad = 41884, catalog = "CYGFM05", intldes = "2016-078A"),
    (fm = 6, name = "CYGNSS FM06", norad = 41889, catalog = "CYGFM06", intldes = "2016-078F"),
    (fm = 7, name = "CYGNSS FM07", norad = 41890, catalog = "CYGFM07", intldes = "2016-078G"),
    (fm = 8, name = "CYGNSS FM08", norad = 41888, catalog = "CYGFM08", intldes = "2016-078E"),
]

# Position-arc fit: a degree-5 polynomial over FIT_SPAN_S of 1 Hz samples.
# The epoch is the first instant of the window, so the arc runs forward from it
# and the fit is evaluated at its own left endpoint. See `flown_epoch_state`.
const FIT_SPAN_S = 300.0
const FIT_ORDER = 5

# The archive's declared fill values for these variables, restated here so an
# unflagged fill would still be caught after the Feather round trip.
const POS_FILL_M = -99999999.0
const VEL_FILL_MPS = -9999.0

# ---------------------------------------------------------------------------
# SPICE
# ---------------------------------------------------------------------------

function furnish_kernels!()
    kernels = [
        joinpath(SPICE_PATH, "lsk", "naif0012.tls"),
        joinpath(SPICE_PATH, "pck", "earth_latest_high_prec.bpc"),
        joinpath(SPICE_PATH, "pck", "pck00011.tpc"),
    ]
    for k in kernels
        isfile(k) || error("missing SPICE kernel: $k")
        furnsh(k)
    end
    return kernels
end

"ITRF93 -> J2000 rotation and its time derivative, from SPICE's 6x6 state transform."
function itrf_state_transform(et::Float64)
    M = sxform("ITRF93", "J2000", et)     # maps [r;v]_ITRF93 to [r;v]_J2000
    # The upper-left block is (L') and the lower-left block is (dL/dt)', in the
    # sense `CygnssICs.earth_fixed_to_inertial` expects.
    l_pi = SMatrix{3, 3, Float64}(M[1:3, 1:3])'
    dl_pi = SMatrix{3, 3, Float64}(M[4:6, 1:3])'
    return (l_pi, dl_pi)
end

# ---------------------------------------------------------------------------
# Flight states
# ---------------------------------------------------------------------------

"""
    flown_track(sub, et_epoch) -> NamedTuple

One spacecraft's whole window in J2000.

The Earth-fixed to inertial step is `CygnssICs.earth_fixed_to_inertial` with
SPICE's exact rotation derivative, which is the repository's verified
conversion: the FM04 mirror carries an independently produced inertial position
and this transform reproduces it to under a millimeter (see the
"flight-data frame check" test). The alternative of adding a constant `omega x r`
formed in the inertial frame is wrong by about 1 m/s at this altitude, which is
why the derivative is carried rather than a spin rate.

The returned `et` is rebuilt from the UTC-based `pvt_unix_seconds` column plus a
single SPICE conversion of the epoch. No leap second falls inside 2025, so the
offset is constant across the window; the result is checked against the FM04
mirror's own `pvt_et_seconds` column in `main`.
"""
function flown_track(sub::AbstractDataFrame, et_epoch::Float64, unix_epoch::Float64)
    n = nrow(sub)
    et = Vector{Float64}(undef, n)
    r_ii = Matrix{Float64}(undef, 3, n)
    v_ii = Matrix{Float64}(undef, 3, n)
    for i in 1:n
        et[i] = et_epoch + (sub.pvt_unix_seconds[i] - unix_epoch)
        l_pi, dl_pi = itrf_state_transform(et[i])
        r, v = CygnssICs.earth_fixed_to_inertial(
            l_pi, dl_pi,
            SVector(sub.sc_pos_x_pvt_m[i], sub.sc_pos_y_pvt_m[i], sub.sc_pos_z_pvt_m[i]),
            SVector(sub.sc_vel_x_pvt_mps[i], sub.sc_vel_y_pvt_mps[i], sub.sc_vel_z_pvt_mps[i]),
        )
        r_ii[:, i] .= r
        v_ii[:, i] .= v
    end
    return (et = et, r_ii = r_ii, v_ii = v_ii)
end

"""
    arc_glitch_count(track; threshold_m) -> (count, worst_m, consistent, threshold_m)

Which samples' positions are inconsistent with their own neighbors, how many,
and the worst one.

Not every navigation fix in this product is good; the reconstruction record puts
the flagged fraction at 1.4% of the rows and traces an apparent 17.4 km outlier
to one of them. The Level 1 quality flags are not among the variables fetched,
so the check made here is dynamic instead of declared. Over a 1-second step the
second difference of position is the centripetal term plus the file's own
quantization: the measured median across this window is 9.1 m and the
99th percentile 20 m, so the 40 m default sits well above the noise and well
below anything that could be an orbit. Only samples that are actually one second
apart are tested, so a gap is never mistaken for a jump.

This flags the boundary of a bad stretch, not its interior, so it is a marker
and not a filter. The one gross case in this window, a six-second FM01 excursion
to a 927 km altitude, is caught by that boundary; everything else it flags is
tens to hundreds of meters.
"""
function arc_glitch_count(track; threshold_m::Float64 = 40.0)
    r = track.r_ii
    et = track.et
    n = size(r, 2)
    consistent = trues(n)
    count = 0
    worst = 0.0
    for i in 2:(n - 1)
        (abs(et[i] - et[i - 1] - 1.0) < 1e-3 && abs(et[i + 1] - et[i] - 1.0) < 1e-3) || continue
        d = norm(@views r[:, i + 1] .- 2 .* r[:, i] .+ r[:, i - 1])
        worst = max(worst, d)
        if d > threshold_m
            count += 1
            consistent[i] = false
        end
    end
    return (count = count, worst_m = worst, consistent = consistent, threshold_m = threshold_m)
end

"""
    flown_epoch_state(track, et_epoch; span_s, order) -> NamedTuple

The J2000 state at the common epoch, fitted from the inertial position arc that
starts there.

The arc is one-sided because the epoch is the first instant of the window and
the product carries nothing before it; the fit is therefore evaluated at its own
left endpoint, not extrapolated beyond it. `main` reports how far the answer
moves when the span and the polynomial order are changed, which is the honest
uncertainty on a one-sided fit.
"""
function flown_epoch_state(track, et_epoch::Float64;
                           span_s::Float64 = FIT_SPAN_S, order::Int = FIT_ORDER)
    idx = findall(t -> et_epoch - 1e-6 <= t <= et_epoch + span_s + 1e-6, track.et)
    length(idx) > order ||
        error("only $(length(idx)) samples within $span_s s of the epoch; need more than $order")
    positions = [SVector(track.r_ii[1, i], track.r_ii[2, i], track.r_ii[3, i]) for i in idx]
    r, v, rms = fit_state_from_arc(track.et[idx], positions, et_epoch; order = order)
    k = argmin(abs.(track.et .- et_epoch))
    return (
        r = r, v = v, resid_rms_m = rms,
        n_samples = length(idx),
        first_index = first(idx), last_index = last(idx),
        nearest_offset_s = track.et[k] - et_epoch,
        # The same instant read straight from the quantized velocity column, for
        # the record: the gap between the two is the quantization the fit avoids.
        quantized_dv = norm(SVector(track.v_ii[1, k], track.v_ii[2, k], track.v_ii[3, k]) - v),
    )
end

# ---------------------------------------------------------------------------
# Catalog states, kept only as the thing being measured
# ---------------------------------------------------------------------------

"""
    tle_state(record, jd_epoch_utc) -> NamedTuple

SGP4 an element set to `jd_epoch_utc` and rotate the result from TEME into
J2000. Returns meters and meters per second, plus the signed propagation span.

WGS-84 constants are used because that is the gravity model the element sets
that CelesTrak distributes are generated against.
"""
function tle_state(record::TleRecord, jd_epoch_utc::Float64)
    tle = read_tles(record.name * "\n" * record.line1 * "\n" * record.line2)[1]
    prop = sgp4_init(tle; sgp4c = sgp4c_wgs84)
    jd_tle = tle_epoch(tle)
    dt_min = (jd_epoch_utc - jd_tle) * 1440.0
    r_teme_km, v_teme_kms = sgp4!(prop, dt_min)
    R = r_eci_to_eci(TEME(), J2000(), jd_epoch_utc)
    r_j2000, v_j2000 = CygnssICs.teme_to_j2000(R, r_teme_km .* 1000.0, v_teme_kms .* 1000.0)
    return (
        r = r_j2000,
        v = v_j2000,
        jd_tle = jd_tle,
        offset_s = dt_min * 60.0,
        tle_epoch_utc = jd_to_iso(jd_tle),
    )
end

function jd_to_iso(jd::Float64)
    dt = julian2datetime(jd)
    return Dates.format(dt, "yyyy-mm-ddTHH:MM:SS") * "Z"
end

# ---------------------------------------------------------------------------
# Report helpers
# ---------------------------------------------------------------------------

function print_elements(label::AbstractString, r, v)
    el = classical_elements(r, v)
    @printf("  %-22s alt %7.2f km (%7.2f x %7.2f)  inc %7.4f deg  period %6.2f min  e %.6f  RAAN %8.3f  u %8.3f\n",
            label, el.altitude_m / 1e3, el.perigee_altitude_m / 1e3, el.apogee_altitude_m / 1e3,
            el.inclination_deg, el.period_s / 60, el.e, el.raan_deg, el.arg_latitude_deg)
    return el
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    println("CYGNSS constellation flown initial conditions at ", CYGNSS_EPOCH_UTC)
    println(repeat("=", 92))

    isfile(ECEF_FILE) || error(string(
        "missing ", relpath(ECEF_FILE, REPO_ROOT), ".\n",
        "Run: python3 scripts/dev/viewer_demos/fetch_cygnss_l1_states.py\n",
        "(one pass over NASA's OPeNDAP service; it needs a NASA Earthdata Login ",
        "bearer token in ~/.edl_token and caches every response so a second run is free)."))

    furnish_kernels!()
    et_epoch = utc2et("2025-06-06T00:00:00")
    jd_epoch = datetime2julian(DateTime(2025, 6, 6, 0, 0, 0))
    unix_epoch = 1_749_168_000.0          # 2025-06-06T00:00:00Z
    @printf("Epoch: ET %.6f s past J2000, JD(UTC) %.6f\n\n", et_epoch, jd_epoch)

    # --- flight states -----------------------------------------------------
    println("Flight states from NASA CYGNSS Level 1 v3.2, via OPeNDAP variable subsetting")
    println(repeat("-", 92))
    df = DataFrame(Arrow.Table(ECEF_FILE))
    @printf("  %s: %d rows, %d spacecraft\n", relpath(ECEF_FILE, REPO_ROOT),
            nrow(df), length(unique(df.spacecraft_num)))

    # The archive's fill values must not have survived into the table.
    any(<=(POS_FILL_M + 1.0), df.sc_pos_x_pvt_m) && error("a position fill value reached the table")
    any(<=(VEL_FILL_MPS + 1.0), df.sc_vel_x_pvt_mps) && error("a velocity fill value reached the table")
    # The mask is materialized into a BitVector first, deliberately. The Feather
    # table arrives as chunked Arrow columns, and indexing a DataFrame of those
    # with an `Arrow.BoolVector`-backed mask returns silently wrong rows on
    # Arrow 2.8.1 with DataFrames 1.8.1: every output row repeats row 1. The
    # failure is quiet, so the guard is kept even though it looks redundant.
    keep = BitVector(coalesce.(collect(df.valid_fix), false))
    all(keep) || @printf("  %d rows carry a fill value and are dropped\n", count(!, keep))
    df = df[keep, :]
    length(unique(df.spacecraft_num)) == 7 ||
        error("expected seven spacecraft after filtering, got $(unique(df.spacecraft_num))")

    tracks = Dict{Int, Any}()
    states = Dict{Int, Any}()
    for e in ROSTER
        sub = sort(df[BitVector(df.spacecraft_num .== e.fm), :], :pvt_unix_seconds)
        if isempty(sub)
            e.fm == 6 || error("no Level 1 samples for FM$(e.fm) in $(ECEF_FILE)")
            continue
        end
        track = flown_track(sub, et_epoch, unix_epoch)
        st = flown_epoch_state(track, et_epoch)
        glitch = arc_glitch_count(track)
        tracks[e.fm] = (track = track, sub = sub, glitch = glitch)
        states[e.fm] = st
        @printf("FM%02d  %d samples  %.3f h  nearest sample %+.6f s from the epoch\n",
                e.fm, nrow(sub), (track.et[end] - track.et[1]) / 3600, st.nearest_offset_s)
        @printf("      fit: degree %d over %d samples (%.0f s), position residual RMS %.3f m; ",
                FIT_ORDER, st.n_samples, FIT_SPAN_S, st.resid_rms_m)
        @printf("quantized velocity column is %.3f m/s away\n", st.quantized_dv)
        @printf("      neighbor consistency: worst second difference %.1f m, %d samples above 40 m\n",
                glitch.worst_m, glitch.count)
        print_elements("elements:", st.r, st.v)
    end
    println()

    # --- how stable is the one-sided fit -----------------------------------
    println("Fit robustness: how far the epoch state moves with the fit window and order")
    println(repeat("-", 92))
    # A variant is only evidence about the fit if it actually describes the arc.
    # A degree-3 polynomial over 600 s leaves a 160 m residual against a product
    # quantized to 1 m: that is a model too poor to represent 38 degrees of orbit,
    # not a measurement of the fit's sensitivity, so it is reported as rejected
    # rather than folded into the spread.
    worst_dr = 0.0
    worst_dv = 0.0
    worst_dv_production_order = 0.0
    rejected = 0
    accepted = 0
    for e in ROSTER
        haskey(tracks, e.fm) || continue
        base = states[e.fm]
        for span in (120.0, 300.0, 600.0), order in (3, 5, 7, 9)
            (span == FIT_SPAN_S && order == FIT_ORDER) && continue
            alt = flown_epoch_state(tracks[e.fm].track, et_epoch; span_s = span, order = order)
            if alt.resid_rms_m > 3 * base.resid_rms_m
                rejected += 1
                continue
            end
            accepted += 1
            worst_dr = max(worst_dr, norm(alt.r - base.r))
            worst_dv = max(worst_dv, norm(alt.v - base.v))
            order == FIT_ORDER &&
                (worst_dv_production_order = max(worst_dv_production_order, norm(alt.v - base.v)))
        end
    end
    @printf("  spans 120/300/600 s, orders 3/5/7/9: %d variants describe the arc, %d rejected\n",
            accepted, rejected)
    @printf("  worst position move %.3f m, worst velocity move %.4f m/s\n", worst_dr, worst_dv)
    @printf("  at the production order %d, worst velocity move across the spans %.4f m/s\n\n",
            FIT_ORDER, worst_dv_production_order)

    # --- against the two states the previous file already had from flight ---
    prev_ics_path = isfile(CATALOG_FILE) ? CATALOG_FILE : OUT_FILE
    if isfile(prev_ics_path)
        println("Recomputed flight states against the previous file's own two flight states")
        println(repeat("-", 92))
        for s_prev in load_constellation_ics(prev_ics_path).spacecraft
            s_prev.provenance == "telemetry" || continue
            idx = findfirst(e -> e.norad == s_prev.norad_id, ROSTER)
            idx === nothing && continue
            fm = ROSTER[idx].fm
            haskey(states, fm) || continue
            @printf("  %-12s dr %.3f m, dv %.4f m/s\n", s_prev.name,
                    norm(states[fm].r - s_prev.r_ii_m), norm(states[fm].v - s_prev.v_ii_m_s))
        end
        println()
    end

    # --- the catalog states these replace ----------------------------------
    println("The catalog states being replaced: ", relpath(TLE_FILE, REPO_ROOT))
    println(repeat("-", 92))
    blocks = parse_tle_blocks(TLE_FILE)
    by_block = Dict{String, Dict{Int, TleRecord}}()
    for (blk, recs) in blocks
        m = Dict{Int, TleRecord}()
        for rec in recs
            entry = findfirst(e -> e.catalog == rec.name, ROSTER)
            entry === nothing && error("element set \"$(rec.name)\" is not a known CYGNSS object")
            m[ROSTER[entry].fm] = rec
        end
        by_block[blk] = m
    end

    # Which spacecraft the previous file took from the catalog, read from the
    # file itself rather than assumed, so the comparison names the real thing.
    previously_catalog = Int[]
    previous_block = "B"
    if isfile(CATALOG_FILE) || isfile(OUT_FILE)
        prev_path = isfile(CATALOG_FILE) ? CATALOG_FILE : OUT_FILE
        prev = load_constellation_ics(prev_path)
        for s in prev.spacecraft
            # "catalogue" is the schema literal the initial-conditions file and
            # the viewer agreed on, so it stays as written even though prose and
            # file names in this repository use American English (CLAUDE.md).
            if s.provenance == "catalogue"
                idx = findfirst(e -> e.norad == s.norad_id, ROSTER)
                idx === nothing || push!(previously_catalog, ROSTER[idx].fm)
                m = match(r"Internet Archive capture (\w+)", s.source)
                m === nothing || (previous_block = m.captures[1])
            end
        end
        sort!(previously_catalog)
        @printf("  previous file %s: %d catalog states (%s), capture %s\n",
                basename(prev_path), length(previously_catalog),
                join((@sprintf("FM%02d", fm) for fm in previously_catalog), ", "), previous_block)
    end
    # A fresh machine has no earlier file to read. The five that were catalog
    # states there are named explicitly so the comparison still labels the right
    # spacecraft; FM01 and FM04 are the two that always had flight telemetry.
    # The generated comparison catalogue contains all seven catalogue states.
    # It is not the historical mixed file. On a repeated build, keep the same
    # historical five/two grouping used on a fresh machine.
    if isempty(previously_catalog) || length(previously_catalog) == length(states)
        previously_catalog = [2, 3, 5, 7, 8]
    end
    println()

    # --- the measurement ---------------------------------------------------
    println("Catalog state minus flown state, at the common epoch")
    println(repeat("-", 92))
    println("  block  sc     |dr| km   radial m   along-track km  cross-track km   |dv| m/s   plane deg   dt s")
    comparison = Dict{String, Any}()
    for blk in sort(collect(keys(by_block)))
        rows = Any[]
        for e in ROSTER
            haskey(tracks, e.fm) && haskey(by_block[blk], e.fm) || continue
            cs = tle_state(by_block[blk][e.fm], jd_epoch)
            fl = states[e.fm]
            dr = cs.r - fl.r
            dv = cs.v - fl.v
            rtn = rtn_offset(fl.r, fl.v, dr)
            plane = orbit_plane_angle_deg(fl.r, fl.v, cs.r, cs.v)
            dt = along_track_time_s(rtn.along_m, norm(fl.v))
            push!(rows, Dict(
                "name" => e.name, "norad_id" => e.norad,
                "was_catalog" => e.fm in previously_catalog,
                "tle_epoch_utc" => cs.tle_epoch_utc,
                "propagation_s" => cs.offset_s,
                "position_error_m" => norm(dr),
                "radial_m" => rtn.radial_m,
                "along_track_m" => rtn.along_m,
                "cross_track_m" => rtn.cross_m,
                "velocity_error_m_s" => norm(dv),
                "orbit_plane_deg" => plane,
                "along_track_time_s" => dt,
            ))
            @printf("  %-6s FM%02d %9.3f %10.1f %15.3f %15.3f %11.3f %11.5f %7.2f%s\n",
                    blk, e.fm, norm(dr) / 1e3, rtn.radial_m, rtn.along_m / 1e3,
                    rtn.cross_m / 1e3, norm(dv), plane, dt,
                    e.fm in previously_catalog ? "  *" : "")
        end
        errs = [r["position_error_m"] for r in rows]
        five = [r["position_error_m"] for r in rows if r["was_catalog"]]
        two = [r["position_error_m"] for r in rows if !r["was_catalog"]]
        comparison[blk] = Dict(
            "spacecraft" => rows,
            "mean_position_error_m" => mean(errs),
            "mean_position_error_m_previously_catalog" => isempty(five) ? NaN : mean(five),
            "mean_position_error_m_previously_telemetry" => isempty(two) ? NaN : mean(two),
        )
        @printf("  %-6s all seven mean %.3f km; the five that were catalog-derived %.3f km; ",
                blk, mean(errs) / 1e3, isempty(five) ? NaN : mean(five) / 1e3)
        @printf("the two that had telemetry %.3f km\n\n", isempty(two) ? NaN : mean(two) / 1e3)
    end
    println("  * marks a spacecraft whose state in the previous file came from the catalog.")
    println()

    # --- preserve the file being replaced ----------------------------------
    if !isfile(CATALOG_FILE) && isfile(OUT_FILE)
        cp(OUT_FILE, CATALOG_FILE)
        println("Preserved the previous file as ", relpath(CATALOG_FILE, REPO_ROOT))
    end
    if !isfile(CATALOG_FILE)
        # Nothing to preserve on a fresh machine: regenerate the catalog states
        # from the committed element sets, so the comparison above stays
        # reproducible without the earlier file.
        entries = Any[]
        for e in ROSTER
            haskey(by_block[previous_block], e.fm) || continue
            cs = tle_state(by_block[previous_block][e.fm], jd_epoch)
            push!(entries, Dict(
                "name" => e.name, "norad_id" => e.norad,
                "r_ii_m" => collect(cs.r), "v_ii_m_s" => collect(cs.v),
                "provenance" => "catalogue",
                "source" => string(
                    "CelesTrak GP element set for NORAD ", e.norad, " (", e.catalog, ", ",
                    e.intldes, "), epoch ", cs.tle_epoch_utc,
                    ", from the Internet Archive capture ", previous_block,
                    " recorded in scripts/dev/viewer_demos/cygnss_historical.tle",
                    "; SGP4 (WGS-84) to the common epoch, TEME -> J2000"),
                "epoch_offset_s" => cs.offset_s,
            ))
        end
        open(CATALOG_FILE, "w") do io
            JSON.print(io, Dict(
                "epoch_utc" => CYGNSS_EPOCH_UTC,
                "frame" => "J2000 Earth-centered inertial, meters and meters per second",
                "spacecraft" => entries,
                "notes" => string(
                    "Catalog states, regenerated from the committed element sets so that the ",
                    "catalog-versus-flight comparison in catalog_vs_flown_20250606.json can be ",
                    "reproduced on a machine that never held the earlier mixed file. These are not ",
                    "flown states; constellation_ics_20250606.json holds the flown ones."),
            ), 2)
        end
        println("Wrote ", relpath(CATALOG_FILE, REPO_ROOT), " (regenerated catalog states)")
    end

    open(COMPARISON_FILE, "w") do io
        JSON.print(io, Dict(
            "epoch_utc" => CYGNSS_EPOCH_UTC,
            "frame" => "J2000 Earth-centered inertial; differences resolved in the flown state's orbit frame",
            "description" => string(
                "Catalog state minus flown state at the common epoch, per archived capture. ",
                "The flown state is the NASA CYGNSS Level 1 v3.2 navigation solution; the catalog ",
                "state is the archived CelesTrak element set propagated to the epoch with SGP4. ",
                "Positive along-track means the catalog state runs ahead of the flown one."),
            "previous_file_capture" => previous_block,
            "previously_catalog_fm" => previously_catalog,
            "captures" => comparison,
        ), 2)
    end
    println("Wrote ", relpath(COMPARISON_FILE, REPO_ROOT))
    println()

    # --- assemble the flown file -------------------------------------------
    println("Constellation at the epoch, all flown")
    println(repeat("-", 92))
    entries = Any[]
    elements = NamedTuple[]
    names = String[]
    for e in ROSTER
        if !haskey(tracks, e.fm)
            e.fm == 6 || error("FM$(e.fm) has no flight track")
            println("FM06 (NORAD 41889) omitted: decayed 2024-06-13 per CelesTrak SATCAT.")
            continue
        end
        st = states[e.fm]
        sub = tracks[e.fm].sub
        src = string(
            "NASA CYGNSS Level 1 Science Data Record v3.2 (collection CYGNSS_L1_V3.2, PO.DAAC), ",
            "granule ", sub.source_file[1],
            "; variables sc_pos_{x,y,z}_pvt, sc_vel_{x,y,z}_pvt, pvt_timestamp_utc read by ",
            "OPeNDAP DAP4 variable subsetting from ",
            "https://opendap.earthdata.nasa.gov/collections/C2832195379-POCLOUD/granules/; ",
            "degree-", FIT_ORDER, " fit to the inertial position over ", st.n_samples,
            " samples spanning ", @sprintf("%.0f", FIT_SPAN_S), " s from the epoch, ",
            "position residual RMS ", @sprintf("%.3f", st.resid_rms_m), " m; ",
            "WGS84 ECEF (ITRF93) -> J2000 by SPICE sxform with earth_latest_high_prec.bpc",
        )
        push!(entries, Dict(
            "name" => e.name, "norad_id" => e.norad,
            "r_ii_m" => collect(st.r), "v_ii_m_s" => collect(st.v),
            "provenance" => "telemetry", "source" => src,
            # The product carries a navigation solution at exactly the common
            # epoch, so nothing is propagated to reach it.
            "epoch_offset_s" => st.nearest_offset_s,
        ))
        el = print_elements(@sprintf("%s [flight]", e.name), st.r, st.v)
        push!(elements, el); push!(names, e.name)
    end
    println()

    # --- sanity checks -----------------------------------------------------
    println("Sanity checks")
    println(repeat("-", 92))
    # Bounds are the constellation as flown in June 2025, not the published
    # design orbit. CYGNSS launched into a roughly 510 km orbit, but it carries
    # no propulsion and by this epoch drag had taken it to about 440 km and a
    # 93.4 minute period. A check written against 510 km would reject the real
    # constellation.
    ok = true
    for (nm, el) in zip(names, elements)
        bad = String[]
        (300e3 <= el.altitude_m <= 700e3) || push!(bad, @sprintf("altitude %.1f km", el.altitude_m / 1e3))
        (34.0 <= el.inclination_deg <= 36.0) || push!(bad, @sprintf("inclination %.3f deg", el.inclination_deg))
        (90 * 60 <= el.period_s <= 100 * 60) || push!(bad, @sprintf("period %.2f min", el.period_s / 60))
        (el.e < 0.01) || push!(bad, @sprintf("eccentricity %.5f", el.e))
        if isempty(bad)
            @printf("  %-14s OK\n", nm)
        else
            ok = false
            @printf("  %-14s FAIL: %s\n", nm, join(bad, ", "))
        end
    end
    println()

    raan = angle_spread_deg([el.raan_deg for el in elements])
    argu = angle_spread_deg([el.arg_latitude_deg for el in elements])
    println("Constellation geometry at the epoch")
    println(repeat("-", 92))
    @printf("  RAAN                 : %d values, span %.2f deg, largest gap %.2f deg\n",
            length(elements), raan.span_deg, raan.largest_gap_deg)
    @printf("                         %s\n", join((@sprintf("%.2f", el.raan_deg) for el in elements), ", "))
    @printf("  argument of latitude : span %.2f deg, largest gap %.2f deg\n",
            argu.span_deg, argu.largest_gap_deg)
    @printf("                         %s\n", join((@sprintf("%.2f", el.arg_latitude_deg) for el in elements), ", "))
    @printf("  mean altitude        : %.2f km   inclination %.4f .. %.4f deg\n",
            mean(el.altitude_m for el in elements) / 1e3,
            minimum(el.inclination_deg for el in elements),
            maximum(el.inclination_deg for el in elements))
    println()

    # --- the track table ---------------------------------------------------
    println("Flight tracks across the window")
    println(repeat("-", 92))
    parts = DataFrame[]
    for e in ROSTER
        haskey(tracks, e.fm) || continue
        sub = tracks[e.fm].sub
        tr = tracks[e.fm].track
        # The Earth-fixed columns keep the product's own Int32 type rather than
        # being widened to Float64: the archive really does quantize them to
        # whole meters and whole meters per second, and storing them as integers
        # says so and halves the file.
        push!(parts, DataFrame(
            name = fill(e.name, nrow(sub)),
            norad_id = fill(Int32(e.norad), nrow(sub)),
            spacecraft_num = fill(Int8(e.fm), nrow(sub)),
            time = tr.et .- et_epoch,
            pvt_unix_seconds = collect(sub.pvt_unix_seconds),
            sc_pos_x_pvt_m = Int32.(sub.sc_pos_x_pvt_m),
            sc_pos_y_pvt_m = Int32.(sub.sc_pos_y_pvt_m),
            sc_pos_z_pvt_m = Int32.(sub.sc_pos_z_pvt_m),
            sc_vel_x_pvt_mps = Int32.(sub.sc_vel_x_pvt_mps),
            sc_vel_y_pvt_mps = Int32.(sub.sc_vel_y_pvt_mps),
            sc_vel_z_pvt_mps = Int32.(sub.sc_vel_z_pvt_mps),
            pos_ii_1 = tr.r_ii[1, :], pos_ii_2 = tr.r_ii[2, :], pos_ii_3 = tr.r_ii[3, :],
            vel_ii_1 = tr.v_ii[1, :], vel_ii_2 = tr.v_ii[2, :], vel_ii_3 = tr.v_ii[3, :],
            arc_consistent = tracks[e.fm].glitch.consistent,
        ))
    end
    all_tracks = reduce(vcat, parts)
    Arrow.write(TRACKS_FILE, all_tracks; compress = :lz4)
    @printf("  %s: %d rows, %d spacecraft, %.1f MB\n", relpath(TRACKS_FILE, REPO_ROOT),
            nrow(all_tracks), length(parts), filesize(TRACKS_FILE) / 1e6)
    loaded = load_constellation_tracks(TRACKS_FILE)
    @printf("  loader check: %s\n", join(track_names(loaded), ", "))
    println()

    # --- write the flown file ----------------------------------------------
    b = comparison[previous_block]
    rows5 = [r for r in b["spacecraft"] if r["was_catalog"]]
    rows2 = [r for r in b["spacecraft"] if !r["was_catalog"]]
    fmt_km(rows, key) = join((@sprintf("%s %.2f km", replace(r["name"], "CYGNSS " => ""),
                                       abs(r[key]) / 1e3) for r in rows), ", ")
    notes = string(
        "Seven spacecraft, not eight, and every one of them flown. CYGNSS FM06 (NORAD 41889, ",
        "2016-078F) is omitted because it was no longer in orbit at this epoch: CelesTrak's SATCAT ",
        "gives it operational status code D, decayed, with decay date 2024-06-13 ",
        "(https://celestrak.org/satcat/records.php?CATNR=41889&FORMAT=json, retrieved 2026-09-16), ",
        "which is why the live GP service returns \"No GP data found\" for that id while the other ",
        "seven are present. NASA had already lost contact with FM06 in November 2022 ",
        "(https://science.nasa.gov/blogs/cygnss/2022/12/09/nasa-team-troubleshooting-out-of-contact-spacecraft-in-cygnss-constellation ",
        "and https://podaac.jpl.nasa.gov/announcements/2022-12-05-CYGNSS-Data-Outage). NASA's CMR ",
        "returns exactly seven Level 1 granules per day over this window, one per surviving ",
        "spacecraft: cyg01 to cyg05, cyg07 and cyg08.",
        "\n\n",
        "Every state is a flight state from the NASA CYGNSS Level 1 Science Data Record version 3.2 ",
        "(collection CYGNSS_L1_V3.2, PO.DAAC/POCLOUD). A whole granule is about 1.09 GB because it ",
        "carries the delay-Doppler map cube, so the granules were not downloaded: the eight ",
        "spacecraft-state variables were read by OPeNDAP DAP4 variable subsetting from ",
        "https://opendap.earthdata.nasa.gov/collections/C2832195379-POCLOUD/granules/. ",
        "scripts/dev/viewer_demos/fetch_cygnss_l1_states.py, which does that fetch and reports what ",
        "it moved, moved 153 MB across the whole window and all seven spacecraft against a ",
        "whole-granule equivalent of 29.84 GB, half a percent of it, and caches every response so a ",
        "second run moves nothing. Access needs a NASA Earthdata Login bearer token.",
        "\n\n",
        "Frame and units come from the granule's own metadata, not from an assumption. The ",
        "sc_pos_{x,y,z}_pvt variables declare units \"meter\" and the comment \"The X component of ",
        "the spacecraft WGS84 reference frame ECEF position, in meters, at pvt_timestamp_utc\"; the ",
        "sc_vel variables declare \"meter s-1\" and the matching ECEF velocity comment; ",
        "pvt_timestamp_utc is \"seconds since\" the granule's own time_coverage_start. Both position ",
        "and velocity are Int32, so the product is quantized to whole meters and whole meters per ",
        "second. The Earth-fixed to J2000 conversion is the repository's verified one ",
        "(scripts/dev/viewer_demos/cygnss_ics.jl, earth_fixed_to_inertial): the full state transform ",
        "from SPICE sxform with earth_latest_high_prec.bpc, carrying the rotation derivative rather ",
        "than adding a spin rate expressed in the inertial frame, which is wrong by about 1 m/s at ",
        "this altitude. It reproduces the independently derived pos_ii columns in the FM04 mirror to ",
        "better than a millimeter across the whole window.",
        "\n\n",
        "Velocity is fitted from the inertial position arc rather than read from the quantized ",
        "velocity column, as the reconstruction record does ",
        "(docs/spaceagora_cygnss_reconstruction_record.md, section 2.1). The archive carries a ",
        "navigation solution at exactly 2025-06-06T00:00:00.000Z for all seven spacecraft, so ",
        "epoch_offset_s is zero everywhere and nothing is propagated; the fit runs forward from the ",
        "epoch and is evaluated at its own left endpoint. Across every fit span of 120, 300 and 600 s ",
        "and polynomial order 3, 5, 7 and 9 whose residual stays at the product's 1 m quantization ",
        "floor, the epoch state moves by at most ", @sprintf("%.1f", worst_dr), " m in position and ",
        @sprintf("%.2f", worst_dv), " m/s in velocity; holding the order at ", FIT_ORDER,
        " and varying only the span, at most ", @sprintf("%.2f", worst_dv_production_order),
        " m/s. So the fitted velocity is good to a few tenths of a meter per second rather than to ",
        "the 0.05 m/s an earlier estimate quoted, but it is still several times better than the ",
        "product's own velocity column, which is quantized to 1 m/s and sits 0.2 to 1.4 m/s from the ",
        "fit here. The combinations that fail the residual check are the low-order fits over long ",
        "arcs, where a cubic cannot represent 38 degrees of orbit; they are excluded and counted, ",
        "not quietly averaged in.",
        "\n\n",
        "What this replaced, and what the replacement measured. The previous version of this file ",
        "had flight states only for FM01 and FM04; ",
        join((@sprintf("FM%02d", fm) for fm in previously_catalog), ", "),
        " were catalog states, archived CelesTrak element sets propagated to this epoch with SGP4. ",
        "That version put an error bar of about 12 km, almost all of it along-track, on those five, ",
        "inferred by scoring the same element sets against the two spacecraft that did have flight ",
        "states. With all seven now flown, the catalog error is measured directly rather than ",
        "inferred, and it is recorded in catalog_vs_flown_20250606.json. For the five that were ",
        "catalog-derived the total position errors are ", fmt_km(rows5, "position_error_m"),
        ", of which the along-track components are ", fmt_km(rows5, "along_track_m"),
        " and the cross-track components are ", fmt_km(rows5, "cross_track_m"),
        ". Their orbit planes differ from the flown planes by ",
        @sprintf("%.4f", minimum(r["orbit_plane_deg"] for r in rows5)), " to ",
        @sprintf("%.4f", maximum(r["orbit_plane_deg"] for r in rows5)), " degrees. ",
        "The two that already had flight states score ", fmt_km(rows2, "position_error_m"),
        " on the same test. So the shape of the earlier estimate held and its size did not. The ",
        "radial and cross-track parts together never exceed ",
        @sprintf("%.1f", 100 * maximum(sqrt(r["radial_m"]^2 + r["cross_track_m"]^2) /
                                       r["position_error_m"] for r in rows5)),
        " percent of the total, so the error is along-track as predicted, and the planes are good to ",
        "well under a hundredth of a degree, also as predicted. But ",
        "the five average ", @sprintf("%.2f", mean(r["position_error_m"] for r in rows5) / 1e3),
        " km and run as high as ", @sprintf("%.2f", maximum(r["position_error_m"] for r in rows5) / 1e3),
        " km, against the ", @sprintf("%.2f", mean(r["position_error_m"] for r in rows2) / 1e3),
        " km that the two flown spacecraft predicted. Two spacecraft were too small a sample to set ",
        "the error bar for five others.",
        "\n\n",
        "The flight track of all seven spacecraft across the full 96-hour window, 1 Hz, in the same ",
        "frame and units as the states here, is in ",
        "cygnss_constellation_tracks_20250606_96hr.feather beside this file; it is what a viewer ",
        "page should draw as each spacecraft's reference ghost.",
    )

    doc = Dict(
        "epoch_utc" => CYGNSS_EPOCH_UTC,
        "frame" => "J2000 Earth-centered inertial, meters and meters per second",
        "spacecraft" => entries,
        "notes" => notes,
    )
    mkpath(TELEM_DIR)
    open(OUT_FILE, "w") do io
        JSON.print(io, doc, 2)
    end
    println("Wrote ", relpath(OUT_FILE, REPO_ROOT), " (", filesize(OUT_FILE), " bytes)")

    # Read it back through the loader the viewer will use.
    ics = load_constellation_ics(OUT_FILE)
    @printf("Loader check: %d spacecraft, epoch %s, provenance %s\n",
            length(ics.spacecraft), ics.epoch_utc,
            join(unique(s.provenance for s in ics.spacecraft), "/"))
    ok || error("one or more states failed the CYGNSS sanity check")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
