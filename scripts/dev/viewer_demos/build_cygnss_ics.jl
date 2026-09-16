#!/usr/bin/env julia
#
# Build the CYGNSS constellation's state at 2025-06-06T00:00:00Z, the start of
# the 96-hour telemetry window, and write it to
# `data/telemetry/CYGNSS/constellation_ics_20250606.json`.
#
#   julia --project=. scripts/dev/viewer_demos/build_cygnss_ics.jl
#
# That directory is gitignored, which is why this script is committed: it is the
# auditable half of the file it produces. Everything it needs is either in the
# repository (the frozen element sets in `cygnss_historical.tle`, the SPICE
# kernels under `data/GRAMSuite.jl/`) or in the gitignored telemetry mirror
# (`cyg01_nasa_pvt_96hr.feather`, `cyg04_nasa_pvt_96hr.feather`). It needs no
# network access.
#
# WHERE EACH STATE COMES FROM
#
#   FM01, FM04  flight. NASA CYGNSS Level-1 PVT, the same product the
#               reconstruction record uses. Earth-fixed in the file; converted
#               here to J2000.
#   the rest    catalog. Historical element sets from an Internet Archive
#               capture of CelesTrak's General Perturbations service, propagated
#               with SGP4 to the common epoch and rotated TEME -> J2000.
#   FM06        absent. CelesTrak's SATCAT records NORAD 41889 decayed on
#               2024-06-13, a year before this epoch, so the constellation at
#               this epoch is seven spacecraft and not eight.
#
# Which archived capture the catalog states come from is decided by measurement,
# not by preference: FM01 and FM04 have both an element set and a flight state,
# so each capture is scored against them and the better one is used for the five
# spacecraft that have no flight state. The score is printed and is written into
# the file's notes as the honest error bar on those five.

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
using .CygnssICs

const SPICE_PATH = joinpath(REPO_ROOT, "data", "GRAMSuite.jl", "GRAM Suite 2.0", "SPICE")
const TELEM_DIR = joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS")
const TLE_FILE = joinpath(@__DIR__, "cygnss_historical.tle")
const OUT_FILE = joinpath(TELEM_DIR, "constellation_ics_20250606.json")

# The eight flight models as launched, with the NORAD ids CelesTrak's SATCAT
# gives for international designator 2016-078 (retrieved 2026-09-16 from
# https://celestrak.org/satcat/records.php?INTDES=2016-078&FORMAT=json).
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

# Flight models with telemetry mirrored locally, and the file each is in.
const TELEMETRY_FILES = Dict(
    1 => "cyg01_nasa_pvt_96hr.feather",
    4 => "cyg04_nasa_pvt_96hr.feather",
)

# Position-arc fit: a degree-5 polynomial over +/- FIT_HALF_WINDOW_S of 1 Hz
# samples. See `telemetry_state` for why velocity is fitted and not read.
const FIT_HALF_WINDOW_S = 150
const FIT_ORDER = 5

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
    # The upper-left block is (L')  and the lower-left block is (dL/dt)', in the
    # sense `CygnssICs.earth_fixed_to_inertial` expects.
    l_pi = SMatrix{3, 3, Float64}(M[1:3, 1:3])'
    dl_pi = SMatrix{3, 3, Float64}(M[4:6, 1:3])'
    return (l_pi, dl_pi)
end

# ---------------------------------------------------------------------------
# Telemetry states
# ---------------------------------------------------------------------------

"""
Ephemeris time of every row, and the frame-converted inertial positions.

`cyg04` carries a `pvt_et_seconds` column; `cyg01` does not. Both carry
`pvt_unix_seconds`, and no leap second falls inside the window, so one SPICE
conversion of the first timestamp plus the UTC offsets reproduces the column
exactly where the column exists (checked at run time below).
"""
function row_ephemeris_times(df::DataFrame)
    first_utc = replace(String(df.pvt_datetime_utc[1]), "+00:00" => "")
    et0 = utc2et(first_utc)
    u0 = df.pvt_unix_seconds[1]
    et = [et0 + (u - u0) for u in df.pvt_unix_seconds]
    if hasproperty(df, :pvt_et_seconds)
        # Tolerance is 1 ms: at 7.6 km/s that is 7.6 m of along-track time tag,
        # and the observed drift is 0.1 ms, below the file's own 1 m position
        # quantization. Anything larger would mean a leap second or a different
        # time system and must not pass silently.
        drift = maximum(abs.(et .- df.pvt_et_seconds))
        drift < 1e-3 || error("rebuilt ephemeris time disagrees with pvt_et_seconds by $drift s")
    end
    return et
end

"""
    telemetry_state(path, et_epoch) -> NamedTuple

The J2000 state at `et_epoch` from a CYGNSS Level-1 PVT file.

Position and velocity both come from a least-squares polynomial fit to the
inertial position over an arc centered on the epoch, rather than from the row
nearest the epoch, for two reasons that the data itself shows:

 1. The file's first sample is one second after the epoch, so *some*
    interpolation is needed regardless.
 2. Both the position and the velocity columns are quantized to integers, 1 m
    and 1 m/s. The reconstruction record says as much ("their difference is
    primarily file quantization: sub-meter in position and 1 m/s in velocity")
    and fits its own initial states "from position over an arc rather than
    taken from a single quantized velocity sample". Measured here, the fitted
    velocity is stable to about 0.05 m/s across half-windows from 30 s to 300 s
    and polynomial orders 3 and 5, while it sits about 1.2 m/s away from the
    quantized velocity column. The fit is the better number by more than an
    order of magnitude, and the residual of the fit is reported so the claim can
    be rechecked.

Rows flagged as invalid navigation fixes are not filtered here; the fit window
is 301 samples and the returned residual RMS would expose a glitch inside it.
"""
function telemetry_state(path::AbstractString, et_epoch::Float64)
    df = DataFrame(Arrow.Table(path))
    et = row_ephemeris_times(df)

    # Sample window centered as closely on the epoch as the file allows.
    # Keep the full sample count even when the epoch sits at the edge of the
    # file, sliding the window rather than truncating it: FM01's first sample is
    # one second after the epoch, so a centered window would be half as long.
    n = length(et)
    width = 2 * FIT_HALF_WINDOW_S + 1
    width <= n || error("telemetry file has only $n rows, fewer than the $width-sample fit window")
    kc = argmin(abs.(et .- et_epoch))
    lo = clamp(kc - FIT_HALF_WINDOW_S, 1, n - width + 1)
    hi = lo + width - 1
    idx = lo:hi

    t0 = et[kc]
    tau = [et[j] - t0 for j in idx]
    pos = Vector{SVector{3, Float64}}(undef, length(idx))
    for (m, j) in enumerate(idx)
        l_pi, dl_pi = itrf_state_transform(et[j])
        r_i, _ = CygnssICs.earth_fixed_to_inertial(
            l_pi, dl_pi,
            SVector{3, Float64}(df.sc_pos_x_pvt_m[j], df.sc_pos_y_pvt_m[j], df.sc_pos_z_pvt_m[j]),
            SVector{3, Float64}(0.0, 0.0, 0.0),
        )
        pos[m] = r_i
    end

    A = hcat((tau .^ p for p in 0:FIT_ORDER)...)
    dtau = et_epoch - t0
    r_fit = zeros(3)
    v_fit = zeros(3)
    fit_coeffs = Vector{Vector{Float64}}(undef, 3)
    sq = 0.0
    for c in 1:3
        y = [p[c] for p in pos]
        co = A \ y
        fit_coeffs[c] = co
        r_fit[c] = sum(co[p + 1] * dtau^p for p in 0:FIT_ORDER)
        v_fit[c] = sum(p * co[p + 1] * dtau^(p - 1) for p in 1:FIT_ORDER)
        sq += sum(abs2, A * co .- y)
    end
    resid_rms = sqrt(sq / (3 * length(idx)))

    # The same instant taken straight from the quantized columns, for the record.
    # Evaluated at the row's own time, not at the epoch: the two are one second
    # apart for FM01, which is 7.6 km of orbital motion and would swamp the
    # quantization the comparison is meant to expose.
    dt_row = et[kc] - t0
    r_at_row = SVector{3, Float64}(ntuple(c -> sum(fit_coeffs[c][p + 1] * dt_row^p for p in 0:FIT_ORDER), 3))
    v_at_row = SVector{3, Float64}(ntuple(c -> sum(p * fit_coeffs[c][p + 1] * dt_row^(p - 1) for p in 1:FIT_ORDER), 3))
    l_pi, dl_pi = itrf_state_transform(et[kc])
    r_q, v_q = CygnssICs.earth_fixed_to_inertial(
        l_pi, dl_pi,
        SVector{3, Float64}(df.sc_pos_x_pvt_m[kc], df.sc_pos_y_pvt_m[kc], df.sc_pos_z_pvt_m[kc]),
        SVector{3, Float64}(df.sc_vel_x_pvt_mps[kc], df.sc_vel_y_pvt_mps[kc], df.sc_vel_z_pvt_mps[kc]),
    )

    return (
        r = SVector{3, Float64}(r_fit),
        v = SVector{3, Float64}(v_fit),
        resid_rms_m = resid_rms,
        n_samples = length(idx),
        first_row = lo,
        last_row = hi,
        row = kc,
        row_utc = String(df.pvt_datetime_utc[kc]),
        row_offset_s = et[kc] - et_epoch,
        quantized_dv = norm(v_at_row - v_q),
        quantized_dr = norm(r_at_row - r_q),
        n_rows = nrow(df),
    )
end

# ---------------------------------------------------------------------------
# Catalog states
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
    @printf("  %-14s alt %7.2f km (%7.2f x %7.2f)  inc %7.4f deg  period %6.2f min  e %.6f  RAAN %8.3f  u %8.3f\n",
            label, el.altitude_m / 1e3, el.perigee_altitude_m / 1e3, el.apogee_altitude_m / 1e3,
            el.inclination_deg, el.period_s / 60, el.e, el.raan_deg, el.arg_latitude_deg)
    return el
end

function rtn_components(r_ref, v_ref, dr)
    rhat = normalize(SVector{3, Float64}(r_ref))
    h = cross(SVector{3, Float64}(r_ref), SVector{3, Float64}(v_ref))
    chat = normalize(h)
    that = cross(chat, rhat)
    d = SVector{3, Float64}(dr)
    return (radial = dot(d, rhat), along = dot(d, that), cross = dot(d, chat))
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    println("CYGNSS constellation initial conditions at ", CYGNSS_EPOCH_UTC)
    println(repeat("=", 78))

    furnish_kernels!()
    et_epoch = utc2et("2025-06-06T00:00:00")
    jd_epoch = datetime2julian(DateTime(2025, 6, 6, 0, 0, 0))
    @printf("Epoch: ET %.6f s past J2000, JD(UTC) %.6f\n\n", et_epoch, jd_epoch)

    # --- flight states -----------------------------------------------------
    println("Flight states from NASA CYGNSS Level-1 PVT")
    println(repeat("-", 78))
    flight = Dict{Int, Any}()
    for fm in sort(collect(keys(TELEMETRY_FILES)))
        path = joinpath(TELEM_DIR, TELEMETRY_FILES[fm])
        isfile(path) || error("missing telemetry file: $path")
        st = telemetry_state(path, et_epoch)
        flight[fm] = st
        @printf("FM%02d  %s  (%d rows)\n", fm, TELEMETRY_FILES[fm], st.n_rows)
        @printf("  nearest row %d at %s, %.6f s from the epoch\n", st.row, st.row_utc, st.row_offset_s)
        @printf("  fit: degree %d over %d samples, position residual RMS %.3f m\n",
                FIT_ORDER, st.n_samples, st.resid_rms_m)
        @printf("  fitted state vs the file's quantized columns: dr %.3f m, dv %.3f m/s\n",
                st.quantized_dr, st.quantized_dv)
        print_elements("elements:", st.r, st.v)
        println()
    end

    # --- element sets ------------------------------------------------------
    println("Historical element sets: ", relpath(TLE_FILE, REPO_ROOT))
    println(repeat("-", 78))
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
    for blk in sort(collect(keys(by_block)))
        @printf("  block %-10s %d element sets: FM%s\n", blk, length(by_block[blk]),
                join(sort(collect(keys(by_block[blk]))), ", FM"))
    end
    println()

    # --- score each capture against flight ---------------------------------
    println("Scoring each capture against the two flight states")
    println(repeat("-", 78))
    scores = Dict{String, Float64}()
    detail = Dict{String, Vector{NamedTuple}}()
    for blk in sort(collect(keys(by_block)))
        errs = Float64[]
        rows = NamedTuple[]
        for fm in sort(collect(keys(flight)))
            haskey(by_block[blk], fm) || continue
            cs = tle_state(by_block[blk][fm], jd_epoch)
            dr = cs.r - flight[fm].r
            dv = cs.v - flight[fm].v
            rtn = rtn_components(flight[fm].r, flight[fm].v, dr)
            push!(errs, norm(dr))
            push!(rows, (fm = fm, dr = norm(dr), dv = norm(dv), offset_d = cs.offset_s / 86400,
                         radial = rtn.radial, along = rtn.along, cross = rtn.cross,
                         tle_epoch = cs.tle_epoch_utc))
            @printf("  %-10s FM%02d  TLE epoch %s  propagated %+7.3f d  |dr| %8.3f km  |dv| %6.3f m/s  (R %+8.3f, T %+9.3f, N %+7.3f km)\n",
                    blk, fm, cs.tle_epoch_utc, cs.offset_s / 86400, norm(dr) / 1e3, norm(dv),
                    rtn.radial / 1e3, rtn.along / 1e3, rtn.cross / 1e3)
        end
        scores[blk] = mean(errs)
        detail[blk] = rows
        @printf("  %-10s mean position error against flight: %.3f km\n\n", blk, mean(errs) / 1e3)
    end

    best = reduce((a, b) -> scores[a] <= scores[b] ? a : b, sort(collect(keys(scores))))
    @printf("Chosen capture: %s (mean flight-referenced position error %.3f km)\n\n",
            best, scores[best] / 1e3)

    # --- assemble ----------------------------------------------------------
    println("Constellation at the epoch")
    println(repeat("-", 78))
    entries = Any[]
    elements = NamedTuple[]
    names = String[]
    for e in ROSTER
        if e.fm == 6
            println("FM06 (NORAD 41889) omitted: decayed 2024-06-13 per CelesTrak SATCAT.")
            continue
        end
        if haskey(flight, e.fm)
            st = flight[e.fm]
            src = string(
                "NASA CYGNSS Level-1 PVT, ", TELEMETRY_FILES[e.fm],
                "; degree-", FIT_ORDER, " fit to the inertial position over rows ",
                st.first_row, "..", st.last_row,
                " (nearest row ", st.row, " at ", st.row_utc,
                "), position residual RMS ", @sprintf("%.3f", st.resid_rms_m), " m; ",
                "ITRF93 -> J2000 by SPICE sxform with earth_latest_high_prec.bpc",
            )
            push!(entries, Dict(
                "name" => e.name, "norad_id" => e.norad,
                "r_ii_m" => collect(st.r), "v_ii_m_s" => collect(st.v),
                "provenance" => "telemetry", "source" => src,
                # Negative: the nearest sample is one second after the common
                # epoch, so the fit is evaluated one second back down its arc.
                "epoch_offset_s" => -st.row_offset_s,
            ))
            el = print_elements(@sprintf("%s [flight]", e.name), st.r, st.v)
            push!(elements, el); push!(names, e.name)
        else
            rec = by_block[best][e.fm]
            cs = tle_state(rec, jd_epoch)
            src = string(
                "CelesTrak GP element set for NORAD ", e.norad, " (", e.catalog, ", ", e.intldes,
                "), epoch ", cs.tle_epoch_utc,
                ", retrieved 2026-09-16 from the Internet Archive capture ", best,
                " recorded in scripts/dev/viewer_demos/cygnss_historical.tle",
                "; SGP4 (WGS-84) to the common epoch, TEME -> J2000",
            )
            push!(entries, Dict(
                "name" => e.name, "norad_id" => e.norad,
                "r_ii_m" => collect(cs.r), "v_ii_m_s" => collect(cs.v),
                "provenance" => "catalogue", "source" => src,
                "epoch_offset_s" => cs.offset_s,
            ))
            el = print_elements(@sprintf("%s [catalog]", e.name), cs.r, cs.v)
            push!(elements, el); push!(names, e.name)
        end
    end
    println()

    # --- sanity checks -----------------------------------------------------
    println("Sanity checks")
    println(repeat("-", 78))
    # Bounds are the constellation as flown in June 2025, not the published
    # design orbit. CYGNSS launched into a roughly 510 km orbit, but it carries
    # no propulsion and by this epoch drag had taken it to about 440 km and a
    # 93.4 minute period; the archived element sets agree (mean motion 15.41 to
    # 15.43 rev/day is a 93.3 to 93.5 minute period), and CelesTrak's SATCAT has
    # the survivors at 380 to 390 km by 2026. A check written against 510 km
    # would reject the real constellation.
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
    println(repeat("-", 78))
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

    # --- write -------------------------------------------------------------
    best_rows = detail[best]
    err_km = join((@sprintf("FM%02d %.2f km", r.fm, r.dr / 1e3) for r in best_rows), " and ")
    cross_km = join((@sprintf("%.2f km for FM%02d", abs(r.cross) / 1e3, r.fm) for r in best_rows), " and ")
    worst_block = reduce((a, b) -> scores[a] >= scores[b] ? a : b, sort(collect(keys(scores))))
    notes = string(
        "Seven spacecraft, not eight. CYGNSS FM06 (NORAD 41889, 2016-078F) is omitted because it was ",
        "no longer in orbit at this epoch: CelesTrak's SATCAT gives it operational status code D, ",
        "decayed, with decay date 2024-06-13 ",
        "(https://celestrak.org/satcat/records.php?CATNR=41889&FORMAT=json, retrieved ",
        "2026-09-16), which is why the live GP service returns \"No GP data found\" for that id while ",
        "the other seven are present. NASA had already lost contact with FM06 in November 2022 ",
        "(https://science.nasa.gov/blogs/cygnss/2022/12/09/nasa-team-troubleshooting-out-of-contact-spacecraft-in-cygnss-constellation ",
        "and https://podaac.jpl.nasa.gov/announcements/2022-12-05-CYGNSS-Data-Outage). Independently, ",
        "NASA's CMR lists Level-1 granules for exactly cyg01..cyg05, cyg07 and cyg08 on 2025-06-06. ",
        "\n\n",
        "FM01 and FM04 are flown states. Their source product reports position and velocity in a ",
        "WGS84 Earth-fixed frame; this was established from the reconstruction record ",
        "(docs/spaceagora_cygnss_reconstruction_record.md, section 2.1), from the column names, and ",
        "confirmed numerically: the median speed in the file is 7.240 km/s (FM01) and 7.245 km/s (FM04) ",
        "at a median radius of 6819 km and 6815 km, so reading those columns as inertial gives a ",
        "semimajor axis of 6185 km and 6181 km, both inside the Earth's 6378 km equatorial radius and ",
        "therefore impossible. Read as Earth-fixed and converted, the inertial speed is 7.6455 km/s ",
        "against a circular speed sqrt(mu/r) of 7.6459 km/s at the same radius, agreeing to better than ",
        "one part in ten thousand. The conversion used here reproduces the independently derived pos_ii columns in ",
        "the FM04 file to better than 1 mm across the whole window. Velocity is fitted from the ",
        "position arc rather than read, because both columns are quantized to integers (1 m and ",
        "1 m/s): the fitted velocity is stable to about 0.05 m/s across fit windows from 30 s to ",
        "300 s and polynomial orders 3 and 5, while the quantized column sits about 1 m/s away, and ",
        "the fit's position residual RMS of about 1.1 m is consistent with that quantization. ",
        "\n\n",
        "The other five are catalog states, not flown states: historical element sets propagated with ",
        "SGP4 across the gap given in epoch_offset_s and rotated from TEME to J2000. Their accuracy is ",
        "measured rather than assumed. Propagating the same capture's element sets for FM01 and FM04, ",
        "which do have flight states, to this epoch gives position errors of ", err_km, ". The error ",
        "is almost entirely along-track: cross-track is ", cross_km, ", so the orbit plane itself ",
        "(right ascension of the ascending node and inclination) is good to under 0.01 degrees and ",
        "only the phasing along the orbit carries the error, worth about 0.1 degrees of argument of ",
        "latitude or 1.6 seconds of orbital timing. That is the error bar to put on the five catalog ",
        "spacecraft. The other archived capture recorded in cygnss_historical.tle, block ",
        worst_block, ", scores ", @sprintf("%.0f", scores[worst_block] / scores[best]),
        " times worse on the same test (", @sprintf("%.1f", scores[worst_block] / 1e3),
        " km mean against ", @sprintf("%.1f", scores[best] / 1e3), " km) and is not used. ",
        "\n\n",
        "Why not better: element sets at the epoch itself are not publicly reachable from this machine. ",
        "CelesTrak's live GP service serves only the newest element set per object, and its own ",
        "historical archive (https://celestrak.org/NORAD/archives/) covers 1980 to 2004 and refers ",
        "any other date to a Special Data Request. Space-Track's gp_history class does serve ",
        "arbitrary past epochs, but it needs a registered Space-Track account (a username and ",
        "password for https://www.space-track.org), and there are no such credentials on this ",
        "machine. Flown states for all seven do exist publicly, in the same NASA CYGNSS ",
        "Level-1 product FM01 and FM04 come from, but downloading them needs a NASA Earthdata Login ",
        "account: the PO.DAAC granule URLs return HTTP 401 and redirect to urs.earthdata.nasa.gov. ",
        "To be exact about what was propagated and in which direction: the five catalog states come ",
        "from element sets whose epochs fall about six days after the common epoch, so SGP4 runs ",
        "backwards over that gap, and the resulting error is the 11 to 13 km measured above. No ",
        "present-day element set was used for anything; the current catalog epoch is 2026 day 259, ",
        "fifteen months after this window, and propagating that backwards was not attempted.",
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
