# Every convention `cygnss_slew.jl` depends on, re-measured from the data rather
# than asserted. Each check prints the numbers behind one claim in the
# `CygnssSlewTelemetry` docstring, and each ends in PASS or FAIL.
#
#   julia --project=. scripts/dev/run.jl viewer_demos/cygnss_slew_checks.jl
#
# It needs the gitignored FM01 telemetry under data/telemetry/CYGNSS/ and, for
# the epoch check, SPICE and the historical element sets in this directory.
include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "cygnss_slew_telemetry.jl"))
using .CygnssSlewTelemetry
using Arrow
using DataFrames
using Dates
using Printf
using Statistics
using SatelliteToolboxSgp4
using SatelliteToolboxTle
using SatelliteToolboxTransformations

const TELEMETRY_DIR = get(ENV, "SPACEAGORA_CYGNSS_DATA", joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS"))
const ADCS_PATH = joinpath(TELEMETRY_DIR, "cyg01_slew_adcs.feather")
const PV_PATH = joinpath(TELEMETRY_DIR, "cyg01_slew_pv_eci.feather")
const CONSTANTS_PATH = joinpath(TELEMETRY_DIR, "cyg01_adcs_constants.toml")
const TLE_PATH = joinpath(@__DIR__, "cygnss_historical.tle")
# CYGNSS FM01 is NORAD 41887 (international designator 2016-078D).
const FM01_NORAD = 41887

passed = Ref(true)
function report(name::AbstractString, ok::Bool)
    println(ok ? "  PASS  " : "  FAIL  ", name)
    ok || (passed[] = false)
    return ok
end

for path in (ADCS_PATH, PV_PATH, CONSTANTS_PATH)
    isfile(path) || error("missing $(path); these checks need the gitignored FM01 telemetry.")
end
tel = load_slew_telemetry(ADCS_PATH, PV_PATH)
constants = load_slew_constants(CONSTANTS_PATH)
n = length(tel.t_rel)

# --- 1. the quaternion convention ------------------------------------------
println("\n1. quaternion convention: the body-frame nadir direction must be nearly constant")
let
    spread_mapped = Float64[]
    spread_raw = Float64[]
    m = length(tel.pv_t_rel)
    mapped = Matrix{Float64}(undef, 3, m)
    raw = Matrix{Float64}(undef, 3, m)
    for (j, i) in enumerate(tel.pv_index)
        r = SVector{3, Float64}(tel.pos_m[:, j]); rhat = -r / norm(r)
        mapped[:, j] .= SM.rot(SVector{4, Float64}(tel.q[:, i])) * rhat
        raw[:, j] .= slew_scalar_first_attitude_matrix(tel.q_telemetry[:, i]) * rhat
    end
    spread_mapped = [std(mapped[k, :]) for k in 1:3]
    spread_raw = [std(raw[k, :]) for k in 1:3]
    @printf("   rot(mapped q):                       std %.4f %.4f %.4f\n", spread_mapped...)
    @printf("   scalar-first matrix, unmapped:       std %.4f %.4f %.4f\n", spread_raw...)
    report("the mapping (x,y,z,w) = (-t2,-t3,-t4,t1) is the nadir-holding one",
        maximum(spread_mapped) < 0.15 && minimum(spread_raw) > 0.2)
end

# --- 2. the body-rate sign --------------------------------------------------
println("\n2. body-rate sign: differentiate the mapped quaternion and solve the kinematics for omega")
let
    ratios = Float64[]
    for i in 2:200:(n - 1)
        h = tel.t_rel[i + 1] - tel.t_rel[i - 1]
        qd = (SVector{4, Float64}(tel.q[:, i + 1]) - SVector{4, Float64}(tel.q[:, i - 1])) / h
        q = SVector{4, Float64}(tel.q[:, i])
        qv = SVector{3, Float64}(q[1], q[2], q[3]); q4 = q[4]
        qdv = SVector{3, Float64}(qd[1], qd[2], qd[3]); qd4 = qd[4]
        # omega = 2 Xi(q)' qdot for Xi(q) = [q4 I + [qv x]; -qv']
        w_fd = 2.0 * (q4 * qdv - cross(qv, qdv) - qv * qd4)
        w_tel = SVector{3, Float64}(tel.omega[:, i])
        push!(ratios, dot(w_fd, w_tel) / dot(w_tel, w_tel))
    end
    @printf("   projection of the finite-difference rate onto the loaded rate: %.4f to %.4f over %d samples\n",
        minimum(ratios), maximum(ratios), length(ratios))
    report("omega_body = +w_eci (a projection of -1 would mean the opposite sign)",
        minimum(ratios) > 0.95 && maximum(ratios) < 1.05)
end

# --- 3. the wheel-axis sign -------------------------------------------------
println("\n3. wheel-axis sign: total angular momentum in inertial space, from measurements only")
let
    function closure_spread(axes)
        h = Matrix{Float64}(undef, 3, n)
        for i in 1:n
            om = SVector{3, Float64}(tel.omega[:, i])
            hw = SVector{3, Float64}(axes * (constants.wheel_inertia .* SVector{3, Float64}(tel.speeds_rad_s[i, :])))
            h[:, i] .= SM.rot(SVector{4, Float64}(tel.q[:, i]))' * (constants.inertia * om + hw)
        end
        window = findall(t -> 890.0 <= t <= 1250.0, tel.t_rel)
        return [std(h[k, window]) for k in 1:3]
    end
    used = closure_spread(constants.wheel_axes)
    unnegated = closure_spread(constants.wheel_axes_from_file)
    @printf("   axes as used (negated): per-axis std %.3e %.3e %.3e N m s\n", used...)
    @printf("   axes as in the file:    per-axis std %.3e %.3e %.3e N m s\n", unnegated...)
    report("the negated matrix holds the conserved momentum better", maximum(used) < 0.6 * maximum(unnegated))
end

# --- 4. the epoch -----------------------------------------------------------
println("\n4. epoch: which reading of the absolute time column puts FM01 on the orbit plane it flew")
let
    planet = Earth("", SPICE_PATH)
    # The GPS receiver's own stamps, when the export carries them, decide
    # between the two J2000 readings that the orbit plane cannot separate.
    adcs = DataFrame(Arrow.Table(ADCS_PATH))
    if all(c -> c in names(adcs), ("gps_week", "gps_sec", "gps_utc_offset_s"))
        gps_utc = DateTime(1980, 1, 6) + Millisecond(round(Int, 1000 * (Float64(adcs.gps_week[1]) * 604800 +
            Float64(adcs.gps_sec[1]) - Float64(adcs.gps_utc_offset_s[1]))))
        et_reading = slew_epoch_utc(tel.t_abs[1])
        utc_reading = SLEW_TIME_ORIGIN + Millisecond(round(Int, tel.t_abs[1] * 1000))
        d_et = (et_reading - gps_utc).value / 1000
        d_utc = (utc_reading - gps_utc).value / 1000
        @printf("   first row: GPS stamp %s UTC; counter as ET %s (%+.3f s); counter as UTC seconds %s (%+.3f s)\n",
            Dates.format(gps_utc, "yyyy-mm-ddTHH:MM:SS.sss"), Dates.format(et_reading, "yyyy-mm-ddTHH:MM:SS.sss"), d_et,
            Dates.format(utc_reading, "yyyy-mm-ddTHH:MM:SS.sss"), d_utc)
        report("the counter is ephemeris time: it lands within the fix age (0 to 2 s) of the GPS stamp, the UTC reading 69 s later",
            0.0 <= d_et <= 2.0 && 68.0 <= d_utc <= 71.0)
    else
        println("   (no GPS stamps in this export; the ET reading rests on the raw-export evidence in the module docstring)")
    end
    r = SVector{3, Float64}(tel.pos_m[:, 1]); v = SVector{3, Float64}(tel.vel_mps[:, 1])
    h = cross(r, v)
    inc_tel = rad2deg(acos(h[3] / norm(h)))
    raan_tel = mod(rad2deg(atan(h[1], -h[2])), 360.0)
    @printf("   telemetry at its first sample: i = %.4f deg, RAAN(J2000) = %.4f deg\n", inc_tel, raan_tel)

    blocks = readlines(TLE_PATH)
    idx = findlast(i -> startswith(blocks[i], "1 $(FM01_NORAD)"), eachindex(blocks))
    idx === nothing && error("no element set for NORAD $(FM01_NORAD) in $(TLE_PATH)")
    tle = read_tles("CYGFM01\n" * blocks[idx] * "\n" * blocks[idx + 1])[1]
    jd_tle = tle_epoch(tle)
    @printf("   element set epoch: JD %.5f\n", jd_tle)

    candidates = (
        ("UTC seconds from 2000-01-01T00:00:00", () -> et_of(Dates.format(DateTime(2000, 1, 1, 0, 0, 0) + Millisecond(round(Int, tel.t_abs[1] * 1000)), "yyyy-mm-ddTHH:MM:SS.sss"))),
        ("UTC seconds from 2000-01-01T12:00:00", () -> et_of(Dates.format(SLEW_TIME_ORIGIN + Millisecond(round(Int, tel.t_abs[1] * 1000)), "yyyy-mm-ddTHH:MM:SS.sss"))),
        ("TDB seconds past J2000 (ET, the reading used)", () -> slew_epoch_et(tel.t_abs[1])),
    )
    best = (name="", draan=Inf, dinc=Inf)
    for (name, et_fn) in candidates
        et = et_fn()
        stamp = slew_epoch_utc(et)
        jd = 2451545.0 + (et - deltet(et, "ET")) / 86400.0
        prop = sgp4_init(tle; sgp4c=sgp4c_wgs84)
        r_teme, v_teme = sgp4!(prop, (jd - jd_tle) * 1440.0)
        rot_teme = SMatrix{3, 3, Float64}(r_eci_to_eci(TEME(), J2000(), jd))
        rs = rot_teme * SVector{3, Float64}(r_teme .* 1000.0)
        vs = rot_teme * SVector{3, Float64}(v_teme .* 1000.0)
        hs = cross(rs, vs)
        inc = rad2deg(acos(hs[3] / norm(hs)))
        raan = mod(rad2deg(atan(hs[1], -hs[2])), 360.0)
        draan = abs(raan - raan_tel); dinc = abs(inc - inc_tel)
        @printf("   %-46s -> %s UTC, %.2f d propagation: dRAAN %.4f deg, di %.4f deg\n",
            name, Dates.format(stamp, "yyyy-mm-ddTHH:MM:SS"), (jd - jd_tle), draan, dinc)
        if draan < best.draan
            best = (name=name, draan=draan, dinc=dinc)
        end
    end
    # The plane separates midnight from noon (3.5 deg of node) but not the two
    # J2000 readings, which differ by 69 s; that is the GPS check above.
    report("a J2000 reading matches the flown plane and the midnight reading does not",
        best.name != "UTC seconds from 2000-01-01T00:00:00" && best.draan < 1.0)
end

# --- 5. the external torque, measured rather than modeled --------------------
println("\n5. external torque: the drift of the conserved momentum, against gravity gradient and the wheel exchange")
let
    adcs = DataFrame(Arrow.Table(ADCS_PATH))
    rods = all(c -> c in names(adcs), ("dpl_cmd_0", "dpl_cmd_1", "dpl_cmd_2")) ?
        hcat(Float64.(adcs.dpl_cmd_0), Float64.(adcs.dpl_cmd_1), Float64.(adcs.dpl_cmd_2)) : nothing
    ratio_window = NaN
    for (lo, hi) in ((890.0, 1100.0), (890.0, 1250.0), (890.0, tel.pv_t_rel[end]))
        L = momentum_ledger(tel, constants; window=(lo, hi), rod_duty=rods)
        @printf("   t_rel %5.0f-%5.0f s: measured dH/dt = [%+.3e %+.3e %+.3e] N m, |.| %.3e\n", lo, hi, L.drift_nm..., norm(L.drift_nm))
        @printf("                        mean gravity gradient = [%+.3e %+.3e %+.3e] N m, |.| %.3e\n", L.gravity_gradient_nm..., norm(L.gravity_gradient_nm))
        L.rod_duty === nothing || @printf("                        torque rods commanded, mean |duty| %.3f %.3f %.3f\n", L.rod_duty...)
        @printf("                        omitted external momentum %.3e N m s against wheel momentum exchanged %.3e N m s: ratio %.3f\n",
            L.omitted_nms, L.exchanged_nms, L.ratio)
        lo == 890.0 && hi == 1100.0 && (ratio_window = L.ratio)
    end
    # A wheels-only replay is meaningful over its window only if the external
    # momentum it omits is smaller than the wheel momentum it carries. The
    # ratio is printed so the reader sees the margin, or its absence.
    report("over the run window the external momentum the replay omits is smaller than the wheel momentum it exchanges",
        isfinite(ratio_window) && ratio_window < 1.0)
end

println()
passed[] || error("cygnss_slew_checks: at least one convention check failed")
println("cygnss_slew_checks: all conventions confirmed")
