# Venus Express reference from ESA's reconstructed ephemeris (ORVM_T19), per
# docs/vex_orvm_reference_preregistration_2026-10-10.md.
#
#   julia --project=. benchmarks/studies/telemetry_validation/vex_orvm/build_reference.jl <ORVM.BSP> <out_dir>
#   julia --project=. benchmarks/studies/telemetry_validation/vex_orvm/build_reference.jl --selfcheck
#
# Writes to <out_dir>: orvm_apsides.csv (every apsis, refined), orvm_duplicates.csv, orvm_perpass_drag.csv,
# orvm_apoapsis_dv.csv, orvm_epoch.toml, convention_check.csv, and the scorer reference files
# vex_orvm_periapsis.feather (orbit = k, k-th periapsis after the epoch, 0-based) and
# vex_orvm_apoapsis.feather (orbit = N, N-th apoapsis after the epoch, 1-based).
using SPICE, SHA, Arrow, CSV, DataFrames, LinearAlgebra, Printf

const ORVM_SHA256 = "6b422c694f96bd33c2104de04433c2853e2aacd3d07205f10c6ac490ab512dc9"
const EPOCH_UTC = "2014-05-19T04:00:00"
const END_UTC = "2014-07-31T00:00:00"
const R_VENUS_KM = 6051.8            # SpaceAGORA Venus Rp_e (vacuum altitude convention)
const MU_KM3S2 = 3.24858599e5        # SpaceAGORA Venus mu
const N_PERI, N_APO = 55, 56         # compared span: periapsides 0..54, apoapsides 1..56
const LSK = joinpath(@__DIR__, "..", "..", "..", "..", "data", "GRAMSuite.jl", "GRAM Suite 2.0", "SPICE", "lsk", "naif0012.tls")

state(et) = spkezr("-248", et, "J2000", "NONE", "VENUS")[1]   # km, km/s
rdotv(s) = s[1] * s[4] + s[2] * s[5] + s[3] * s[6]
energy(s) = 0.5 * (s[4]^2 + s[5]^2 + s[6]^2) - MU_KM3S2 / norm(s[1:3])

# Bisection on f between a and b (f(a), f(b) of opposite sign) to tol seconds.
function refine(f, a, b; tol=1e-3)
    fa = f(a)
    while b - a > tol
        m = 0.5 * (a + b); fm = f(m)
        (fa < 0) == (fm < 0) ? (a = m; fa = fm) : (b = m)
    end
    return 0.5 * (a + b)
end

# Apsides as sign changes of r.v on a grid, refined; kind :peri for - to +, :apo for + to -.
function apsides(statefn, t0, t1, dt)
    ts = collect(t0:dt:t1); g = [rdotv(statefn(t)) for t in ts]
    out = Tuple{Symbol, Float64}[]
    for i in 1:length(ts)-1
        if g[i] < 0 <= g[i+1]
            push!(out, (:peri, refine(t -> rdotv(statefn(t)), ts[i], ts[i+1])))
        elseif g[i] > 0 >= g[i+1]
            push!(out, (:apo, refine(t -> rdotv(statefn(t)), ts[i], ts[i+1])))
        end
    end
    return out
end

# Analytic Kepler orbit check of `apsides` and `refine` (periapsis time and radius to 1 ms / 1 mm).
function selfcheck()
    a, e, mu = 39000.0, 0.84, MU_KM3S2; n = sqrt(mu / a^3); T = 2pi / n
    function kep(t)
        M = n * (t - 1000.0); E = M
        for _ in 1:50; E -= (E - e * sin(E) - M) / (1 - e * cos(E)); end
        r = a * (1 - e * cos(E)); x = a * (cos(E) - e); y = a * sqrt(1 - e^2) * sin(E)
        Ed = n / (1 - e * cos(E)); vx = -a * sin(E) * Ed; vy = a * sqrt(1 - e^2) * cos(E) * Ed
        return [x, y, 0.0, vx, vy, 0.0]
    end
    aps = apsides(kep, 2000.0, 2.2T, 5.0)
    peris = [t for (k, t) in aps if k === :peri]; apos = [t for (k, t) in aps if k === :apo]
    @assert length(peris) == 2 && length(apos) == 2
    @assert abs(peris[1] - (1000.0 + T)) < 1e-3 && abs(peris[2] - (1000.0 + 2T)) < 1e-3
    @assert abs(apos[1] - (1000.0 + T / 2)) < 1e-3
    @assert abs(norm(kep(peris[1])[1:3]) - a * (1 - e)) < 1e-6
    println("selfcheck ok: periapsis times to 1 ms, radius to 1 mm")
end

# The simulator's extractor on sampled states (linear interpolation at the r.v sign change), for the convention check.
function linear_extractor(t, S)
    r = [norm(s[1:3]) for s in S]; alt = r .- R_VENUS_KM; g = [rdotv(s) / norm(s[1:3]) for s in S]
    out = Tuple{Symbol, Float64, Float64}[]
    for i in 1:length(t)-1
        d = abs(g[i]) + abs(g[i+1]); w = d == 0 ? 0.5 : abs(g[i]) / d
        te = (1 - w) * t[i] + w * t[i+1]; ae = (1 - w) * alt[i] + w * alt[i+1]
        g[i] < 0 && g[i+1] >= 0 && push!(out, (:peri, te, ae))
        g[i] > 0 && g[i+1] <= 0 && push!(out, (:apo, te, ae))
    end
    return out
end

function main(bsp, outdir)
    h = open(io -> bytes2hex(sha256(io)), bsp)
    h == ORVM_SHA256 || error("ORVM SHA-256 mismatch: $h")
    furnsh(LSK); furnsh(bsp); mkpath(outdir)
    et0 = str2et(EPOCH_UTC); et1 = str2et(END_UTC)
    raw = apsides(state, et0, et1, 5.0)

    # Discontinuities can create a second, spurious r.v sign change pair: keep the extreme apsis of each kind within 2 h.
    rows = DataFrame(kind=String[], et=Float64[], alt_km=Float64[])
    for (k, t) in raw; push!(rows, (String(k), t, norm(state(t)[1:3]) - R_VENUS_KM)); end
    keep = trues(nrow(rows)); dups = DataFrame(kind=String[], kept_utc=String[], dropped_utc=String[], kept_alt_km=Float64[], dropped_alt_km=Float64[])
    for kind in ("peri", "apo")
        idx = findall(==(kind), rows.kind)
        for (i, j) in zip(idx[1:end-1], idx[2:end])
            keep[i] && keep[j] || continue
            if rows.et[j] - rows.et[i] < 7200
                better_j = kind == "peri" ? rows.alt_km[j] < rows.alt_km[i] : rows.alt_km[j] > rows.alt_km[i]
                kk, dd = better_j ? (j, i) : (i, j)
                keep[dd] = false
                push!(dups, (kind, et2utc(rows.et[kk], "ISOC", 3), et2utc(rows.et[dd], "ISOC", 3), rows.alt_km[kk], rows.alt_km[dd]))
            end
        end
    end
    rows = rows[keep, :]
    # The spurious member of a pair has the other kind's sign change next to it; drop opposite-kind events within 2 h of a dropped one.
    for d in eachrow(dups)
        et_d = str2et(d.dropped_utc)
        bad = findall(i -> rows.kind[i] != d.kind && abs(rows.et[i] - et_d) < 7200 &&
            (d.kind == "apo" ? rows.alt_km[i] > 60000 : rows.alt_km[i] < 2000), 1:nrow(rows))
        deleteat!(rows, bad)
    end
    rows.utc = [et2utc(t, "ISOC", 3) for t in rows.et]
    rows.t_since_epoch_s = rows.et .- et0
    rows.index = zeros(Int, nrow(rows))
    np = 0; na = 0
    for i in 1:nrow(rows)
        rows.kind[i] == "peri" ? (rows.index[i] = np; np += 1) : (na += 1; rows.index[i] = na)
    end
    S = [state(t) for t in rows.et]
    for (c, nm) in enumerate(("x_km", "y_km", "z_km", "vx_kmps", "vy_kmps", "vz_kmps"))
        rows[!, nm] = [s[c] for s in S]
    end
    rows.speed_kmps = [norm(s[4:6]) for s in S]
    CSV.write(joinpath(outdir, "orvm_apsides.csv"), rows)
    CSV.write(joinpath(outdir, "orvm_duplicates.csv"), dups)

    P = rows[rows.kind .== "peri", :]; A = rows[rows.kind .== "apo", :]
    Arrow.write(joinpath(outdir, "vex_orvm_periapsis.feather"), DataFrame(orbit=Float64.(P.index[1:N_PERI]), altitude=P.alt_km[1:N_PERI]))
    Arrow.write(joinpath(outdir, "vex_orvm_apoapsis.feather"), DataFrame(orbit=Float64.(A.index[1:N_APO]), altitude=A.alt_km[1:N_APO]))

    # Per-pass drag delta-V: two-body energy change across +-1800 s around each periapsis, over the periapsis speed.
    pd = DataFrame(k=P.index, utc=P.utc, alt_km=P.alt_km, v_p_kmps=P.speed_kmps)
    pd.dE_km2s2 = [energy(state(t + 1800)) - energy(state(t - 1800)) for t in P.et]
    pd.drag_dv_mps = -1000 .* pd.dE_km2s2 ./ pd.v_p_kmps
    CSV.write(joinpath(outdir, "orvm_perpass_drag.csv"), pd)

    # Apoapsis delta-V: energy change across +-3600 s around each apoapsis, over the apoapsis speed (burns at apoapsis).
    ad = DataFrame(N=A.index, utc=A.utc, alt_km=A.alt_km, v_a_kmps=A.speed_kmps)
    ad.dE_km2s2 = [energy(state(t + 3600)) - energy(state(t - 3600)) for t in A.et]
    ad.dv_mps = 1000 .* ad.dE_km2s2 ./ ad.v_a_kmps
    pk = P.alt_km
    ad.next_peri_step_km = [let i = findfirst(>(t), P.et); (i === nothing || i == 1) ? NaN : pk[i] - pk[i-1] end for t in A.et]
    CSV.write(joinpath(outdir, "orvm_apoapsis_dv.csv"), ad)

    # Epoch state and osculating elements (J2000).
    s0 = state(et0); el = oscltx(s0, et0, MU_KM3S2)
    open(joinpath(outdir, "orvm_epoch.toml"), "w") do io
        println(io, "epoch_utc = \"$EPOCH_UTC\"")
        println(io, "initial_state_j2000_m = [", join([@sprintf("%.6f", 1000 * x) for x in s0[1:3]], ", "), ", ",
            join([@sprintf("%.9f", 1000 * x) for x in s0[4:6]], ", "), "]")
        @printf(io, "rp_alt_km = %.6f\nra_km = %.6f\necc = %.9f\nta_deg = %.6f\nperiod_h = %.6f\n",
            el[1] - R_VENUS_KM, el[10] * (1 + el[2]), el[2], rad2deg(el[9]), el[11] / 3600)
    end

    # Convention check: the simulator's linear extractor on ORVM sampled at 60 s and 5 s against the refined apsides.
    # Grids are offset by half a step from the refined apsis (the worst case for the linear extractor) and by a quarter step.
    cc = DataFrame(kind=String[], index=Int[], refined_alt_km=Float64[], dt_s=Float64[], phase=Float64[], lin_alt_km=Float64[], lin_dt_s=Float64[])
    for r in eachrow(rows)
        r.kind == "peri" && r.index > 54 && continue
        r.kind == "apo" && r.index > 56 && continue
        for dt in (60.0, 5.0, 1.0), phase in (0.5, 0.25)
            ts = collect((r.et - 1800 + phase * dt):dt:(r.et + 1800)); ev = filter(e -> e[1] === Symbol(r.kind), linear_extractor(ts, [state(t) for t in ts]))
            push!(cc, (r.kind, r.index, r.alt_km, dt, phase, isempty(ev) ? NaN : ev[1][3], isempty(ev) ? NaN : ev[1][2] - r.et))
        end
    end
    CSV.write(joinpath(outdir, "convention_check.csv"), cc)
    println("apsides: ", count(==("peri"), rows.kind), " peri, ", count(==("apo"), rows.kind), " apo; duplicates dropped: ", nrow(dups))
end

if !isempty(ARGS) && ARGS[1] == "--selfcheck"
    selfcheck()
else
    length(ARGS) == 2 || error("usage: build_reference.jl <ORVM.BSP> <out_dir> | --selfcheck")
    main(ARGS[1], ARGS[2])
end
