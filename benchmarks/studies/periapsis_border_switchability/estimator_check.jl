# Part 3: feasibility check of a SINGLE-averaged third-body estimator (offline; no tuning).
# Per pass: orbit-average of the Gauss variational equations over the Kepler ellipse (eccentric-anomaly nodes, uniform in E,
# weight (1 - e cos E) = dM/dE) for a point-mass tidal perturbation from bodies held fixed at their position at the pass-interval
# midpoint, plus J2 secular rates; mean elements integrated pass to pass (midpoint rule, steps = truth pass times).
#   julia --project=. estimator_check.jl OUT [N]
include(joinpath(@__DIR__, "common.jl"))   # SPICE helpers, example support (GravitationalHarmonicsModel, Earth)
using LinearAlgebra, StaticArrays, CSV, DataFrames

const NQ = parse(Int, get(ARGS, 2, "256"))   # eccentric-anomaly nodes (ASSUMED; converged check below)

orbit_frame(i, Ω, ω) = begin
    h = SVector(sin(i) * sin(Ω), -sin(i) * cos(Ω), cos(i)); node = SVector(cos(Ω), sin(Ω), 0.0)
    P = cos(ω) * node + sin(ω) * cross(h, node); Q = cross(h, P)
    (P, Q, h)
end

"""Tidal acceleration of a body at geocentric `rb` on a spacecraft at `r`: exact point-mass (direct + indirect) or its quadrupole term."""
function tidal(r, rb, μb, quad::Bool)
    if quad
        rb3 = norm(rb)^3; u = rb / norm(rb)
        return μb / rb3 * (3 * dot(r, u) * u - r)
    end
    d = rb - r
    return μb * (d / norm(d)^3 - rb / norm(rb)^3)
end

"""Orbit-averaged rates (da, de, di, dΩ, dω) per second for bodies `bodies` = [(rb, μb), ...] (Gauss variational equations, numerical average)."""
function averaged_rates(a, e, i, Ω, ω, μ, bodies; quad=false, N=NQ)
    P, Q, W = orbit_frame(i, Ω, ω); p = a * (1 - e^2); h = sqrt(μ * p); acc = zeros(5)
    for k in 1:N
        E = 2π * (k - 0.5) / N; ce, se = cos(E), sin(E)
        r = a * (ce - e) * P + a * sqrt(1 - e^2) * se * Q; rm = norm(r)
        ν = atan(sqrt(1 - e^2) * se, ce - e); cν, sν = cos(ν), sin(ν)
        rh = r / rm; th = cross(W, rh)
        f = sum(tidal(r, rb, μb, quad) for (rb, μb) in bodies)
        R, S, Wc = dot(f, rh), dot(f, th), dot(f, W)
        g = SVector(2a^2 / h * (e * sν * R + p / rm * S),
            (p * sν * R + ((p + rm) * cν + rm * e) * S) / h,
            rm * cos(ω + ν) / h * Wc,
            rm * sin(ω + ν) / (h * sin(i)) * Wc,
            (-p * cν * R + (p + rm) * sν * S) / (h * e) - rm * sin(ω + ν) * cos(i) / (h * sin(i)) * Wc)
        acc .+= (1 - e * ce) .* g   # dM = (1 - e cos E) dE
    end
    return acc ./ N
end

"""Secular J2 rates of (Ω, ω) (Brouwer/Kozai; Vallado, 'Fundamentals of Astrodynamics and Applications', secular J2 section)."""
function j2_rates(a, e, i, μ, J2, Re)
    n = sqrt(μ / a^3); p = a * (1 - e^2); f = n * J2 * (Re / p)^2
    (dΩ=-1.5 * f * cos(i), dω=0.75 * f * (5 * cos(i)^2 - 1))
end

body_pos(name, et) = SVector{3, Float64}(spkpos(name, et, "J2000", "NONE", "EARTH")[1]) .* 1e3

"""Total rate vector (a, e, i, Ω, ω)/s at pass-midpoint epoch `et`."""
function total_rates(el, et, names, μ, J2, Re; quad=false, N=NQ, jfac=1.0)
    bodies = [(body_pos(uppercase(b), et), MU_BODY[b]) for b in names]
    r = isempty(bodies) ? zeros(5) : averaged_rates(el..., μ, bodies; quad, N)
    j = j2_rates(el[1], el[2], el[3], μ, J2, Re)
    r[4] += jfac * j.dΩ; r[5] += jfac * j.dω
    return r
end

"""Integrate mean elements over truth pass times (midpoint rule). `epoch_et0`: ephemeris time of t=0."""
function integrate_mean(rows, et0, names, μ, J2, Re; quad=false, N=NQ, epoch=:mid)
    el = [rows.a_m[1], rows.e[1], rows.i_rad[1], rows.raan_rad[1], rows.omega_rad[1]]
    out = [copy(el)]
    for k in 1:nrow(rows)-1
        t0, t1 = rows.time_s[k], rows.time_s[k+1]; dt = t1 - t0
        tb = epoch == :mid ? et0 + (t0 + t1) / 2 : et0 + t0
        r1 = total_rates(el, tb, names, μ, J2, Re; quad, N)
        elm = el .+ r1 .* (dt / 2)
        r2 = total_rates(elm, tb, names, μ, J2, Re; quad, N)
        el = el .+ r2 .* dt
        push!(out, copy(el))
    end
    return out
end

"""One-step prediction: start each step from the TRUTH osculating elements at pass k, predict pass k+1."""
function one_step(rows, et0, names, μ, J2, Re; quad=false, N=NQ)
    out = Vector{Vector{Float64}}()
    for k in 1:nrow(rows)-1
        el = [rows.a_m[k], rows.e[k], rows.i_rad[k], rows.raan_rad[k], rows.omega_rad[k]]
        t0, t1 = rows.time_s[k], rows.time_s[k+1]; dt = t1 - t0; tb = et0 + (t0 + t1) / 2
        r1 = total_rates(el, tb, names, μ, J2, Re; quad, N); elm = el .+ r1 .* (dt / 2)
        r2 = total_rates(elm, tb, names, μ, J2, Re; quad, N); push!(out, el .+ r2 .* dt)
    end
    return out
end

"""Previous studies' double-averaged hand formula, integrated along the TRUTH elements (as in the grid analysis)."""
function hand_history(rows, et0, names, μ)
    rp = [rows.radius_m[1]]
    for k in 1:nrow(rows)-1
        t0, t1 = rows.time_s[k], rows.time_s[k+1]; tm = et0 + (t0 + t1) / 2
        rate = 0.0; a, e = rows.a_m[k], rows.e[k]; nm = sqrt(μ / a^3)
        for b in names
            pa = perturber_angles(b, tm, rows.i_rad[k], rows.raan_rad[k], rows.omega_rad[k])
            rate += -a * (15 / 8) * (pa.mu / pa.r_m^3 / nm) * e * sqrt(1 - e^2) * sin(2pa.omega_p) * sin(pa.i_p)^2
        end
        push!(rp, rp[end] + rate * (t1 - t0))
    end
    return rp
end

# ---- convention check: double-averaged limit ----------------------------------------------------
function double_average_check(μ)
    a, e, i, Ω, ω = 36.4e6, 0.82, deg2rad(40.0), deg2rad(30.0), deg2rad(70.0)
    n = normalize(SVector(0.2, -0.5, 0.84)); q1 = normalize(cross(n, SVector(0.0, 0.0, 1.0))); q2 = cross(n, q1)
    rb, μb = 3.844e8, MU_BODY["Moon"]; M = 360
    de = 0.0; di = 0.0; dΩ = 0.0; dω = 0.0; da = 0.0
    for m in 1:M
        λ = 2π * (m - 0.5) / M; u = cos(λ) * q1 + sin(λ) * q2
        r = averaged_rates(a, e, i, Ω, ω, μ, [(rb * u, μb)]; quad=true, N=256); da += r[1]; de += r[2]; di += r[3]; dΩ += r[4]; dω += r[5]
    end
    de /= M; da /= M
    # hand formula (same angles definition as perturber_angles)
    h = SVector(sin(i) * sin(Ω), -sin(i) * cos(Ω), cos(i)); node = SVector(cos(Ω), sin(Ω), 0.0)
    ehat = cos(ω) * node + sin(ω) * cross(h, node); k = normalize(cross(n, h))
    ωp = atan(dot(h, cross(k, ehat)), dot(k, ehat)); ip = acos(dot(h, n)); nm = sqrt(μ / a^3)
    hand = (15 / 8) * (μb / rb^3 / nm) * e * sqrt(1 - e^2) * sin(2ωp) * sin(ip)^2
    return (de_numeric=de, de_hand=hand, ratio=de / hand, da_numeric=da, omega_p_deg=rad2deg(ωp), i_p_deg=rad2deg(ip))
end

function main()
    out = ARGS[1]; study = load_study(); nbody_planet()
    cfgv = load_corridor_config(CORRIDOR_CONFIG_PATH).values
    gpath = joinpath(REPO_ROOT, String(cfgv["gravity"]["coefficient_file"]))
    harm = GravitationalHarmonicsModel(3, 0, gpath, make_no_gram_planet(:earth))   # GM, Re and C20 of the sourced coefficient file
    μ, Re, J2 = harm.gm_m3s2, harm.reference_radius_m, -sqrt(5.0) * harm.C[3, 1]
    println("mu=$μ Re=$Re J2=$J2 NQ=$NQ")
    chk = double_average_check(μ); println("double_average_check: ", chk)
    open(joinpath(out, "estimator_convention_check.toml"), "w") do io
        for (k, v) in pairs(chk); println(io, "$k = $v"); end
        println(io, "mu = $μ\nRe = $Re\nJ2 = $J2\nNQ = $NQ")
    end
    # quadrature convergence on one real state
    rows0 = CSV.read(joinpath(out, "runs", "I0234_R000_W090_E1", "passes_v5.csv"), DataFrame)
    et0c = SimulationModel.ephemerides_time_seconds(parse_epoch(study["grid"]["epochs_utc"][1]), SpiceEphemeridesModel())
    el0 = [rows0.a_m[1], rows0.e[1], rows0.i_rad[1], rows0.raan_rad[1], rows0.omega_rad[1]]
    conv = [(N, total_rates(el0, et0c + 3e4, ["Moon", "Sun"], μ, J2, Re; N=N)) for N in (32, 64, 128, 256, 512, 1024)]
    open(joinpath(out, "estimator_quadrature_convergence.csv"), "w") do io
        println(io, "N,da_mps,de_per_s,di,draan,domega")
        for (N, r) in conv; println(io, "$N,$(join(r, ','))"); end
    end
    # J2 check against the neither-variant (v3) truth
    tags = unique(CSV.read(joinpath(out, "attribution_run_summary.csv"), DataFrame).tag)
    series = DataFrame(); stats = DataFrame()
    for tag in tags
        ei = parse(Int, last(split(tag, "_E"))); et0 = SimulationModel.ephemerides_time_seconds(parse_epoch(study["grid"]["epochs_utc"][ei]), SpiceEphemeridesModel())
        for (v, names) in (("v5", ["Moon", "Sun"]), ("v8", ["Moon"]), ("v9", ["Sun"]), ("v0", ["Moon", "Sun"]), ("v3", String[]))
            rows = CSV.read(joinpath(out, "runs", tag, "passes_$v.csv"), DataFrame)
            nrow(rows) < 8 && continue
            sa = integrate_mean(rows, et0, names, μ, J2, Re); sq = integrate_mean(rows, et0, names, μ, J2, Re; quad=true)
            sp = integrate_mean(rows, et0, names, μ, J2, Re; epoch=:start)
            os = one_step(rows, et0, names, μ, J2, Re); osq = one_step(rows, et0, names, μ, J2, Re; quad=true)
            hand = isempty(names) ? fill(rows.radius_m[1], nrow(rows)) : hand_history(rows, et0, names, μ)
            rp(el) = el[1] * (1 - el[2])
            df = DataFrame(tag=tag, variant=v, k=1:nrow(rows), time_s=rows.time_s, rp_truth=rows.radius_m,
                a_truth=rows.a_m, e_truth=rows.e, i_truth=rows.i_rad, raan_truth=rows.raan_rad, om_truth=rows.omega_rad,
                rp_sa=rp.(sa) .- rp(sa[1]) .+ rows.radius_m[1], rp_sq=rp.(sq) .- rp(sq[1]) .+ rows.radius_m[1], rp_sa_start=rp.(sp) .- rp(sp[1]) .+ rows.radius_m[1],
                rp_hand=hand, a_sa=[x[1] for x in sa], e_sa=[x[2] for x in sa], i_sa=[x[3] for x in sa], raan_sa=[x[4] for x in sa], om_sa=[x[5] for x in sa],
                rp_1step=[missing; [rp(x) - rp([rows.a_m[k], rows.e[k], 0, 0, 0]) + rows.radius_m[k] for (k, x) in enumerate(os)]],
                rp_1step_q=[missing; [rp(x) - rp([rows.a_m[k], rows.e[k], 0, 0, 0]) + rows.radius_m[k] for (k, x) in enumerate(osq)]])
            append!(series, df)
        end
        println("done $tag"); flush(stdout)
    end
    CSV.write(joinpath(out, "estimator_series.csv"), series)
end
main()
