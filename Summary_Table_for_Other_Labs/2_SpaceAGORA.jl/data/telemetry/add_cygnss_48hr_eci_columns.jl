##
# Stage 2 of the CYGNSS 48hr data pipeline: convert the raw WGS84 ECEF GPS PVT
# (from convert_cygnss_48hr_xlsx.py) to J2000 and write the feather consumed by
# the verification harness (data/telemetry/CYGNSS/cygnss_data_48hr.feather).
#
# The rotation uses the repo's own SPICE kernels and the same ITRF93 frame the
# harness uses for planet-fixed comparisons, so telemetry and simulation share
# one frame convention. The output deliberately does NOT carry the raw
# OBS4.* position/velocity column names: the harness prefers those names over
# pos_ii_* and would interpret the ECEF values as inertial.
#
# Output columns: time [s from epoch], pos_ii_1..3 [m, J2000],
# vel_ii_1..3 [m/s, J2000], sc_pos_{x,y,z}_ecef_m, sc_vel_{x,y,z}_ecef_mps,
# numsats, gdop_x10, pvt_valid, utc.
#
# Usage: julia --project=. data/telemetry/add_cygnss_48hr_eci_columns.jl
##
using Arrow, DataFrames, SPICE, StaticArrays, Printf, Statistics

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
const SPICE_ROOT = joinpath(REPO, "data", "GRAMSuite.jl", "GRAM Suite 2.0", "SPICE")
const IN_FEATHER = joinpath(@__DIR__, "CYGNSS", "cygnss_48hr_raw_ecef.feather")
const OUT_FEATHER = joinpath(@__DIR__, "CYGNSS", "cygnss_data_48hr.feather")
const EPOCH_UTC = "2025-06-06T00:00:00"
const MU_EARTH = 3.98600436233e14  # m^3/s^2, matches src/environment/ephemerides/planets.jl

for k in ("lsk/naif0012.tls", "pck/earth_latest_high_prec.bpc", "pck/pck00011.tpc")
    p = joinpath(SPICE_ROOT, k)
    isfile(p) || error("missing SPICE kernel $p")
    furnsh(p)
end

df = DataFrame(Arrow.Table(IN_FEATHER))
const _PVT_COLS = [
    "OBS4.ENG_PVT.DDMI_PVT_SCPOS_X (m)", "OBS4.ENG_PVT.DDMI_PVT_SCPOS_Y (m)",
    "OBS4.ENG_PVT.DDMI_PVT_SCPOS_Z (m)", "OBS4.ENG_PVT.DDMI_PVT_SCVEL_X (m/s)",
    "OBS4.ENG_PVT.DDMI_PVT_SCVEL_Y (m/s)", "OBS4.ENG_PVT.DDMI_PVT_SCVEL_Z (m/s)",
]
n_raw = nrow(df)
df = dropmissing(df, _PVT_COLS)
println("dropped $(n_raw - nrow(df)) rows with missing PVT values ($(n_raw) -> $(nrow(df)))")
n = nrow(df)
t = Float64.(df[!, "TIME OFFSET"])
et0 = utc2et(EPOCH_UTC)

rx = Float64.(df[!, "OBS4.ENG_PVT.DDMI_PVT_SCPOS_X (m)"])
ry = Float64.(df[!, "OBS4.ENG_PVT.DDMI_PVT_SCPOS_Y (m)"])
rz = Float64.(df[!, "OBS4.ENG_PVT.DDMI_PVT_SCPOS_Z (m)"])
vx = Float64.(df[!, "OBS4.ENG_PVT.DDMI_PVT_SCVEL_X (m/s)"])
vy = Float64.(df[!, "OBS4.ENG_PVT.DDMI_PVT_SCVEL_Y (m/s)"])
vz = Float64.(df[!, "OBS4.ENG_PVT.DDMI_PVT_SCVEL_Z (m/s)"])

pos_ii = Matrix{Float64}(undef, n, 3)
vel_ii = Matrix{Float64}(undef, n, 3)
for i in 1:n
    m = sxform("ITRF93", "J2000", et0 + t[i])
    s = m * SVector{6, Float64}(rx[i], ry[i], rz[i], vx[i], vy[i], vz[i])
    pos_ii[i, :] .= s[1:3]
    vel_ii[i, :] .= s[4:6]
end

out = DataFrame(
    "time" => t,
    "pos_ii_1" => pos_ii[:, 1], "pos_ii_2" => pos_ii[:, 2], "pos_ii_3" => pos_ii[:, 3],
    "vel_ii_1" => vel_ii[:, 1], "vel_ii_2" => vel_ii[:, 2], "vel_ii_3" => vel_ii[:, 3],
    "sc_pos_x_ecef_m" => rx, "sc_pos_y_ecef_m" => ry, "sc_pos_z_ecef_m" => rz,
    "sc_vel_x_ecef_mps" => vx, "sc_vel_y_ecef_mps" => vy, "sc_vel_z_ecef_mps" => vz,
    "numsats" => df[!, "OBS4.ENG_PVT.DDMI_PVT_NUMSATS (numsats)"],
    "gdop_x10" => df[!, "OBS4.ENG_PVT.DDMI_PVT_GDOP (GDOP)"],
    "pvt_valid" => df[!, "OBS4.ENG_PVT.DDMI_PVT_VALID (null)"],
    "utc" => String.(df[!, "ENG_PVT"]),
)
Arrow.write(OUT_FEATHER, out)
println("wrote $OUT_FEATHER  rows=$n")

# --- validation ---------------------------------------------------------------
r_eci = sqrt.(out.pos_ii_1 .^ 2 .+ out.pos_ii_2 .^ 2 .+ out.pos_ii_3 .^ 2)
v_eci = sqrt.(out.vel_ii_1 .^ 2 .+ out.vel_ii_2 .^ 2 .+ out.vel_ii_3 .^ 2)
v_ecef = sqrt.(vx .^ 2 .+ vy .^ 2 .+ vz .^ 2)
sma_km = 1.0e-3 ./ (2.0 ./ r_eci .- v_eci .^ 2 ./ MU_EARTH)
valid = out.pvt_valid .== 2
@printf("|r|:      %.3f .. %.3f km\n", minimum(r_eci) / 1e3, maximum(r_eci) / 1e3)
@printf("|v| ECI:  %.4f .. %.4f km/s (ECEF was %.4f .. %.4f)\n",
    minimum(v_eci) / 1e3, maximum(v_eci) / 1e3, minimum(v_ecef) / 1e3, maximum(v_ecef) / 1e3)
@printf("vis-viva SMA (all):      median %.4f km  IQR [%.4f, %.4f]\n",
    median(sma_km), quantile(sma_km, 0.25), quantile(sma_km, 0.75))
@printf("vis-viva SMA (valid==2): median %.4f km  IQR [%.4f, %.4f]\n",
    median(sma_km[valid]), quantile(sma_km[valid], 0.25), quantile(sma_km[valid], 0.75))

# First-sample orbital elements — compare against the published Table 5 IC
# (SMA 6818.8611 km, ecc 4.7900e-4, inc 34.9357 deg, RAAN 177.3709,
# AOP 140.6318, TA 276.6553), which was derived from this same first sample.
r1 = SVector{3, Float64}(out.pos_ii_1[1], out.pos_ii_2[1], out.pos_ii_3[1])
v1 = SVector{3, Float64}(out.vel_ii_1[1], out.vel_ii_2[1], out.vel_ii_3[1])
h = cross(r1, v1)
evec = cross(v1, h) / MU_EARTH - r1 / norm(r1)
sma1 = 1.0 / (2.0 / norm(r1) - norm(v1)^2 / MU_EARTH)
inc1 = acosd(h[3] / norm(h))
nvec = SVector{3, Float64}(-h[2], h[1], 0.0)
raan1 = mod(atand(nvec[2], nvec[1]), 360.0)
ecc1 = norm(evec)
aop1 = mod(
    (evec[3] >= 0 ? 1 : -1) * acosd(clamp(dot(nvec, evec) / (norm(nvec) * ecc1), -1, 1)),
    360.0,
)
ta1 = mod(
    (dot(r1, v1) >= 0 ? 1 : -1) * acosd(clamp(dot(evec, r1) / (ecc1 * norm(r1)), -1, 1)),
    360.0,
)
@printf("first-sample elements: sma %.4f km  ecc %.4e  inc %.4f  raan %.4f  aop %.4f  ta %.4f\n",
    sma1 / 1e3, ecc1, inc1, raan1, aop1, ta1)
println("expected (Table 5):    sma 6818.8611 km  ecc 4.7900e-04  inc 34.9357  raan 177.3709  aop 140.6318  ta 276.6553")
