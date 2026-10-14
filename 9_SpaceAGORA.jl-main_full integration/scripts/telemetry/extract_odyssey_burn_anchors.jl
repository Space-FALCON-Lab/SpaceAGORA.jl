# Extract the Odyssey NAV-kernel states that anchor the aerobraking replay after
# each trim burn, and print the manifest's [scenarios.state_anchors] block.
#
# Apoapsis passages are found as maxima of the Mars-centred radius of the
# reconstructed trajectory; periapses are numbered from P19 (the scenario epoch
# sits 16 min after P19, so the first periapsis after the epoch is P20). Burn B
# fires at the apoapsis after P(B); its anchor is that apoapsis plus the burn
# duration at the manifest's thrust and dry-plus-propellant mass, plus 600 s,
# so the replayed burn has finished. A burn whose anchor would fall within a
# quarter period of the next periapsis is skipped (the walk-out raise).
#
# Usage:
#   SPACEAGORA_SPICE_PATH=<dir with lsk/ pck/ spk/> julia --project=. scripts/telemetry/extract_odyssey_burn_anchors.jl
# The kernel m01_ab_v2.bsp (NAIF ody-m-spice-6-v1.0) must sit in spk/missions/.
using SPICE, LinearAlgebra, Printf, TOML
using Arrow, DataFrames
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SP = get(ENV, "SPACEAGORA_SPICE_PATH", joinpath(ROOT, "data", "GRAMSuite.jl", "GRAM Suite 2.0", "SPICE"))
furnsh(joinpath(SP, "lsk", "naif0012.tls"))
furnsh(joinpath(SP, "pck", "pck00011.tpc"))
for f in readdir(joinpath(SP, "spk", "planets"); join=true); endswith(lowercase(f), ".bsp") && furnsh(f); end
furnsh(joinpath(SP, "spk", "missions", "m01_ab_v2.bsp"))
manifest = TOML.parsefile(joinpath(ROOT, "test", "telemetry_benchmark_manifest.toml"))
ody = first(s for s in manifest["scenarios"] if s["name"] == "odyssey")
it = ody["initial_time"]
epoch = @sprintf("%04d-%02d-%02dT%02d:%02d:%06.3f", it["year"], it["month"], it["day"], it["hour"], it["minute"], it["second"])
et0 = str2et(epoch)
state(et) = (s = spkezr("-53", et, "J2000", "NONE", "MARS")[1]; (s[1:3] .* 1000.0, s[4:6] .* 1000.0))
radius(et) = 1000.0 * norm(spkezr("-53", et, "J2000", "NONE", "MARS")[1][1:3])
# sanity: initial state vs manifest
r0, v0 = state(et0)
println("epoch=", epoch, " et0=", et0)
println("kernel state at epoch (m, m/s): ", round.(r0; digits=3), " ", round.(v0; digits=6))
println("manifest initial_state_j2000_m:   ", ody["initial_state_j2000_m"])
# scan radius for extrema
span_s = 80 * 86400.0
dt = 30.0
ts = collect(et0:dt:(et0 + span_s))
rs = radius.(ts)
function refine(f, a, b; tol=1e-3, maximize=true)   # golden-section
    g = (sqrt(5) - 1) / 2
    c = b - g * (b - a); d = a + g * (b - a)
    while abs(b - a) > tol
        fc = f(c); fd = f(d)
        better = maximize ? (fc > fd) : (fc < fd)
        if better; b = d; else; a = c; end
        c = b - g * (b - a); d = a + g * (b - a)
    end
    return (a + b) / 2
end
peri = Float64[]; apo = Float64[]
for k in 2:length(ts)-1
    if rs[k] < rs[k-1] && rs[k] <= rs[k+1]; push!(peri, refine(radius, ts[k-1], ts[k+1]; maximize=false)); end
    if rs[k] > rs[k-1] && rs[k] >= rs[k+1]; push!(apo, refine(radius, ts[k-1], ts[k+1]; maximize=true)); end
end
println("periapses found: ", length(peri), "  apoapses found: ", length(apo))
Rp = 3396.19e3
# number periapses: first after epoch is P20
pnum(idx) = 19 + idx
println("first periapsis after epoch: P", pnum(1), " at ", et2utc(peri[1], "ISOC", 0), " alt_km=", round((radius(peri[1]) - Rp) / 1e3; digits=3))
burns = ody["maneuvers"]["orbit_numbers"]; dvs = ody["maneuvers"]["delta_v_mps"]
thrust = ody["maneuvers"]["thrust_n"]; m0 = ody["spacecraft"]["bus_mass_kg"] + 2 * ody["spacecraft"]["panel_mass_each_kg"] + ody["spacecraft"]["prop_mass_kg"]
println("burn mass estimate kg=", m0, " thrust N=", thrust)
rows = []
for (B, dv) in zip(burns, dvs)
    B >= 20 || continue                       # pre-epoch burns are in the initial condition
    ip = B - 19                               # index of P(B) in peri
    (ip >= 1 && ip + 1 <= length(peri)) || (println("skip burn P", B, ": periapsis not in span"); continue)
    tpB, tpB1 = peri[ip], peri[ip + 1]
    ia = findfirst(t -> tpB < t < tpB1, apo)
    ia === nothing && (println("skip burn P", B, ": no apoapsis between P", B, " and P", B + 1); continue)
    t_apo = apo[ia]
    t_burn = abs(dv) * m0 / thrust
    t_anchor = t_apo + t_burn + 600.0
    period = tpB1 - tpB
    (t_anchor < tpB1 - 0.25 * period) || (println("skip burn P", B, ": anchor too close to next periapsis (period ", round(period/3600; digits=2), " h)"); continue)
    r, v = state(t_anchor)
    push!(rows, (B=B, dv=dv, t_apo=t_apo, t_anchor=t_anchor, period_h=period/3600, apo_alt_km=(radius(t_apo)-Rp)/1e3, r=r, v=v))
end
apo_df = DataFrame(Arrow.Table(joinpath(ROOT, "data", "telemetry", "Odyssey", "True_Odyssey_apoapsis_alts_kernel.feather")))
feather_apo(B) = (i = findfirst(==(Float64(B)), apo_df.orbit); i === nothing ? NaN : apo_df.altitude[i])
println("\nanchors: ", length(rows), "   (apo_alt_km from kernel radius - Rp_e vs feather apoapsis altitude at that orbit)")
for row in rows
    @printf("P%-4d dv=%7.3f  apo=%s  apo_alt_km=%9.2f  feather=%9.2f  period_h=%5.2f  anchor_elapsed_s=%12.3f\n", row.B, row.dv, et2utc(row.t_apo, "ISOC", 0), row.apo_alt_km, feather_apo(row.B), row.period_h, row.t_anchor - et0)
end
# TOML block
out_path = get(ENV, "SPACEAGORA_ANCHORS_OUT", joinpath(ROOT, "output", "odyssey_state_anchors.toml"))
mkpath(dirname(out_path))
open(out_path, "w") do io
    println(io, "[scenarios.state_anchors]")
    println(io, "enabled = true")
    println(io, "burn_orbit_numbers = [", join(string.(r.B for r in rows), ", "), "]")
    println(io, "elapsed_s = [", join((@sprintf("%.3f", r.t_anchor - et0) for r in rows), ", "), "]")
    println(io, "states_j2000_m = [")
    for r in rows
        @printf(io, "    [%.6f, %.6f, %.6f, %.9f, %.9f, %.9f],\n", r.r[1], r.r[2], r.r[3], r.v[1], r.v[2], r.v[3])
    end
    println(io, "]")
end
println("wrote ", out_path)
