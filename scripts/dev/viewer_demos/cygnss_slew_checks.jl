# Every convention `cygnss_slew.jl` depends on, re-measured from the data rather
# than asserted. Each check prints the numbers behind one claim in the
# `CygnssSlewTelemetry` docstring, and each ends in PASS or FAIL.
#
#   julia --project=. scripts/dev/viewer_demos/cygnss_slew_checks.jl
#
# It needs the gitignored FM01 telemetry under data/telemetry/CYGNSS/ and, for
# the epoch check, SPICE and the historical element sets in this directory.
include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "cygnss_slew_telemetry.jl"))
using .CygnssSlewTelemetry
using Dates
using Printf
using Statistics
using SatelliteToolboxSgp4
using SatelliteToolboxTle
using SatelliteToolboxTransformations

const TELEMETRY_DIR = joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS")
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
        ("seconds from 2000-01-01T00:00:00", DateTime(2000, 1, 1, 0, 0, 0)),
        ("seconds from 2000-01-01T12:00:00", SLEW_TIME_ORIGIN),
    )
    best = (name="", draan=Inf, dinc=Inf)
    for (name, origin) in candidates
        stamp = origin + Millisecond(round(Int, tel.t_abs[1] * 1000))
        et = et_of(Dates.format(stamp, "yyyy-mm-ddTHH:MM:SS.sss"))
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
        @printf("   %-34s -> %s UTC, %.2f d propagation: dRAAN %.4f deg, di %.4f deg\n",
            name, Dates.format(stamp, "yyyy-mm-ddTHH:MM:SS"), (jd - jd_tle), draan, dinc)
        if draan < best.draan
            best = (name=name, draan=draan, dinc=dinc)
        end
    end
    report("the J2000-epoch reading is the one that matches the flown plane",
        best.name == "seconds from 2000-01-01T12:00:00" && best.draan < 1.0)
end

# --- 5. the external torque, measured rather than modeled --------------------
println("\n5. external torque: the drift of the conserved momentum, against gravity gradient alone")
let
    mu = 3.986004418e14
    m = length(tel.pv_t_rel)
    h = Matrix{Float64}(undef, 3, m)      # total angular momentum, inertial
    tau = Matrix{Float64}(undef, 3, m)    # gravity-gradient torque, inertial
    for (j, i) in enumerate(tel.pv_index)
        q = SVector{4, Float64}(tel.q[:, i])
        c_bi = SM.rot(q)
        omega = SVector{3, Float64}(tel.omega[:, i])
        hw = SVector{3, Float64}(constants.wheel_axes * (constants.wheel_inertia .* SVector{3, Float64}(tel.speeds_rad_s[i, :])))
        h[:, j] .= c_bi' * (constants.inertia * omega + hw)
        r_body = c_bi * SVector{3, Float64}(tel.pos_m[:, j])
        rr = norm(r_body); rhat = r_body / rr
        tau[:, j] .= c_bi' * (3 * mu / rr^3 * cross(rhat, constants.inertia * rhat))
    end
    ok = true
    for (lo, hi) in ((890.0, 1250.0), (890.0, tel.pv_t_rel[end]))
        w = findall(x -> lo <= x <= hi, tel.pv_t_rel)
        tt = tel.pv_t_rel[w]
        drift = [sum((tt .- mean(tt)) .* (h[k, w] .- mean(h[k, w]))) / sum((tt .- mean(tt)) .^ 2) for k in 1:3]
        gg = [mean(tau[k, w]) for k in 1:3]
        @printf("   t_rel %5.0f-%5.0f s: measured dH/dt = [%+.3e %+.3e %+.3e] N m\n", lo, hi, drift...)
        @printf("                        mean gravity gradient = [%+.3e %+.3e %+.3e] N m  (ratio %+.2f %+.2f %+.2f)\n",
            gg..., (drift ./ gg)...)
        # The claim the run of record rests on: over the hour the NET external
        # torque is a small fraction of gravity gradient alone, so carrying
        # gravity gradient without the torque that cancels it is worse than
        # carrying neither.
        hi > 3000.0 && (ok = maximum(abs.(drift ./ gg)) < 0.2)
    end
    report("over the hour the net external torque is under a fifth of gravity gradient alone", ok)
end

println()
passed[] || error("cygnss_slew_checks: at least one convention check failed")
println("cygnss_slew_checks: all conventions confirmed")
