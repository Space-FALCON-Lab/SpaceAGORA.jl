#!/usr/bin/env julia
# Build a tabulated plume-surface gas field for a named engine and write it to
# data/psi/<engine>.json.
#
#   julia --project=. scripts/dev/psi/build_plume_field.jl                 # apollo_lmde
#   julia --project=. scripts/dev/psi/build_plume_field.jl --out /tmp/x.json
#
# The physics, the sources for every constant and the assumptions are documented
# on `PlumeNozzle` and `build_plume_field_table` in
# src/dynamics/coupled/force_torque_models/plume_gas_field.jl; data/psi/README.md
# records what the shipped tables were built from. The script prints a short
# report so a rebuild can be compared with the published reference points
# without opening the table.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."); io=devnull)

using SpaceAGORA
using Printf

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

# Engines the script knows how to build. Add one by adding a `PlumeNozzle`.
#
# `apollo_lmde` is the shipped default: the plain Simons source flow with its
# point source at the center of the nozzle exit plane. `apollo_lmde_virtual_source`
# is the same engine with the source moved 1.8 m downstream, the offset Morris
# 2012 section 4.6 measured for this nozzle; it reproduces that work's DSMC wall
# shear stress far better in the near field and is unusable within about 3 m of
# the ground, where the source would reach the surface and the field freezes.
# Both are written so a study can show the sensitivity; `data/psi/README.md`
# tabulates the difference.
const ENGINES = Dict{String, PlumeNozzle}(
    "apollo_lmde" => PlumeNozzle(),
    "apollo_lmde_virtual_source" => PlumeNozzle(name="apollo_lmde_virtual_source", virtual_source_offset_m=1.8),
)

function parse_args(args)
    engine = "apollo_lmde"
    out = ""
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--engine" && i < length(args)
            engine = args[i + 1]; i += 2
        elseif a == "--out" && i < length(args)
            out = args[i + 1]; i += 2
        elseif startswith(a, "--engine=")
            engine = split(a, "=", limit=2)[2]; i += 1
        elseif startswith(a, "--out=")
            out = split(a, "=", limit=2)[2]; i += 1
        else
            error("unknown argument $(a); usage: build_plume_field.jl [--engine NAME] [--out PATH]")
        end
    end
    haskey(ENGINES, engine) || error("unknown engine " * engine * "; known: " * join(sort(collect(keys(ENGINES))), ", "))
    isempty(out) && (out = joinpath(REPO_ROOT, "data", "psi", engine * ".json"))
    return engine, out
end

function main(args)
    engine, out = parse_args(args)
    nozzle = ENGINES[engine]
    table = build_plume_field_table(nozzle;
        source="scripts/dev/psi/build_plume_field.jl; Simons/Boynton source flow with the " *
               "nozzle-boundary-layer correction (EUCASS 2019-661 eqs 3-14) and classical " *
               "Newtonian impingement; engine $(engine)")
    cfg = PlumeSurfaceConfig(field=table)
    path = save_plume_field(out, table)
    @printf("wrote %s (%.0f kB)\n", path, filesize(path) / 1024)
    @printf("  engine            %s, reference thrust %.1f kN, exit diameter %.2f m\n",
            table.name, table.reference_thrust_n / 1e3, table.exit_diameter_m)
    @printf("  limit speed       %.0f m/s, chamber temperature %.0f K, gamma %.3f, R %.0f J/(kg K)\n",
            table.limit_speed_mps, table.chamber_temperature_k, table.gamma, table.gas_constant_j_kg_k)
    @printf("  grid              %d heights (h/D %.2f to %.0f) x %d radii (r/h 0 to %.2f)\n",
            length(table.height_over_diameter), first(table.height_over_diameter),
            last(table.height_over_diameter), length(table.radius_over_height),
            last(table.radius_over_height))

    # Reference points, printed so a rebuild can be checked by eye.
    println("\n  model versus the published reference points")
    # Morris 2012 section 4.4.2: the LMDE hovering 5 m above the surface at
    # 13.34 kN, peak laminar smooth-wall shear stress 92 Pa (DSMC).
    st = [plume_gas_state(table, cfg, 13_340.0, 5.0, r) for r in range(0.0, 12.0; length=1201)]
    peak_tau, k = findmax(s -> s.shear_pa, st)
    @printf("  h = 5 m, F = 13.34 kN: peak wall shear %.1f Pa at r = %.2f m (Morris 2012 s4.4.2 DSMC: 92 Pa)\n",
            peak_tau, (k - 1) * 12.0 / 1200)
    p0, R = plume_surface_footprint(cfg, 13_340.0, 5.0)
    @printf("                         stagnation pressure %.0f Pa, footprint radius %.2f m\n", p0, R)
    # Momentum: the integral of the wall pressure over the plane against the thrust.
    for h in (5.0, 20.0, 50.0)
        rs = range(0.0, 8.0 * h; length=4001)
        integral = 0.0
        for i in 2:length(rs)
            pa = plume_gas_state(table, cfg, 13_340.0, h, rs[i - 1]).pressure_pa * rs[i - 1]
            pb = plume_gas_state(table, cfg, 13_340.0, h, rs[i]).pressure_pa * rs[i]
            integral += pi * (pa + pb) * (rs[i] - rs[i - 1])
        end
        @printf("  h = %4.0f m: surface pressure integral / thrust = %.3f\n", h, integral / 13_340.0)
    end
    println()
    return path
end

main(ARGS)
