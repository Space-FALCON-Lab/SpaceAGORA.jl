##
# Build the along-track NOAA WAM-IPE density table consumed by
# run_drag_decay_calibration.jl: sample the assimilated density at the flight
# track's own (time, lat, lon, alt) at a fixed cadence and write
# data/telemetry/CYGNSS/wam_density_table.csv (time_s, rho columns).
#
# Requires WamIPEDensity.jl (github.com/Bourbon8464/WamIPEDensity.jl,
# develop branch with the fixes from PR #6: GridMetadata arity, fill_values
# field name, WFS cycle-folder flooring) available in the load path, e.g.
#     git clone -b develop https://github.com/Bourbon8464/WamIPEDensity.jl ~/.julia/dev/WamIPEDensity.jl
# and network access to the public NOAA bucket (noaa-nws-wam-ipe-pds).
#
# Usage: julia --project=<WamIPEDensity env> benchmarks/studies/cygnss_attitude_torque_reconstruction/build_wam_density_table.jl
##
using Dates

try
    @eval using WamIPEDensity
catch
    error("WamIPEDensity.jl not available in this environment — see the header for setup.")
end
using Arrow, DataFrames, CSV

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const TELEMETRY = joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "cygnss_data_48hr.feather")
const OUT = joinpath(REPO_ROOT, "data", "telemetry", "CYGNSS", "wam_density_table.csv")
const EPOCH = DateTime(2025, 6, 6, 0, 0, 0)
const CADENCE_S = 60.0
const RP_E = 6.3781366e6

isfile(TELEMETRY) || error("Missing $(TELEMETRY) — run scripts/dev/fetch_private_telemetry.sh")

df = DataFrame(Arrow.Table(TELEMETRY))
df = df[df.pvt_valid .== 2, :]
keep = [1]
for i in 2:nrow(df)
    if df.time[i] - df.time[keep[end]] >= CADENCE_S
        push!(keep, i)
    end
end
df = df[keep, :]
println("sampling WAM-IPE at $(nrow(df)) track points ($(CADENCE_S) s cadence)...")

itp = WAMInterpolator()
rhos = Vector{Float64}(undef, nrow(df))
for (i, row) in enumerate(eachrow(df))
    # ECEF -> geocentric lat/lon and geometric altitude (spherical, adequate
    # for sampling a ~2-degree density grid).
    x, y, z = row.sc_pos_x_ecef_m, row.sc_pos_y_ecef_m, row.sc_pos_z_ecef_m
    rr = sqrt(x^2 + y^2 + z^2)
    lat = asind(z / rr)
    lon = atand(y, x)
    h = rr - RP_E
    rhos[i] = get_density_at_point(itp, EPOCH + Millisecond(round(Int, row.time * 1000)),
        lat, lon, h; angles_in_deg=true)
    i % 200 == 0 && (println("  $(i)/$(nrow(df))"); flush(stdout))
end
bad = count(isnan, rhos)
CSV.write(OUT, DataFrame(time_s=df.time, rho=rhos))
println("wrote $(OUT): $(nrow(df)) rows, $(bad) NaN, rho range $(extrema(filter(!isnan, rhos)))")
