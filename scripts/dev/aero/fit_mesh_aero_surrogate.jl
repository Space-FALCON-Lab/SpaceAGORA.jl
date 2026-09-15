# Fit the free-molecular surrogate of a CAD model from the command line.
#
#   julia --project=. scripts/dev/aero/fit_mesh_aero_surrogate.jl --model=data/models/magellan_nasa_3d_resources.glb \
#         --scale=1.0 --rotation=-180,90,90 --out=output/aero/magellan.json [--degree=8] [--poly=2] \
#         [--directions=1200] [--speed-ratios=3,5,8,12,20] [--sigma-n=1] [--sigma-t=1] [--grid=256] \
#         [--ref-area=<m2>] [--ref-length=<m>] [--no-shadowing] [--outward-normals]
#
# Prints the fit residuals and a few reference coefficients (drag along each
# body axis at the middle speed ratio) so the normalisation can be sanity
# checked against the box model's numbers.
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if something(Base.active_project(), "") != joinpath(REPO_ROOT, "Project.toml")
    import Pkg
    Pkg.activate(REPO_ROOT; io=devnull)
end
using SpaceAGORA, StaticArrays

function parse_args(args)
    opts = Dict{String, String}()
    for a in args
        startswith(a, "--") || throw(ArgumentError("unexpected argument $(a)"))
        k, v = occursin("=", a) ? split(a[3:end], "=", limit=2) : (a[3:end], "true")
        opts[String(k)] = String(v)
    end
    return opts
end

opts = parse_args(ARGS)
haskey(opts, "model") || throw(ArgumentError("--model=<file> is required"))
model = opts["model"]
scale = parse(Float64, get(opts, "scale", "1.0"))
rotation = Tuple(parse.(Float64, split(get(opts, "rotation", "0,0,0"), ",")))
out = get(opts, "out", joinpath(REPO_ROOT, "output", "aero", splitext(basename(model))[1] * "_aero_surrogate.json"))
degree = parse(Int, get(opts, "degree", "8"))
poly = parse(Int, get(opts, "poly", "2"))
ndir = parse(Int, get(opts, "directions", "1200"))
srs = Tuple(parse.(Float64, split(get(opts, "speed-ratios", "3,5,8,12,20"), ",")))
sigma_n = parse(Float64, get(opts, "sigma-n", "1.0"))
sigma_t = parse(Float64, get(opts, "sigma-t", "1.0"))
grid = parse(Int, get(opts, "grid", "256"))
ref_area = haskey(opts, "ref-area") ? parse(Float64, opts["ref-area"]) : nothing
ref_length = haskey(opts, "ref-length") ? parse(Float64, opts["ref-length"]) : nothing
shadowing = !haskey(opts, "no-shadowing")
outward = haskey(opts, "outward-normals")

t0 = time()
panels = mesh_aero_panels(model; scale=scale, rotation_deg=rotation, reference_area_m2=ref_area, reference_length_m=ref_length, outward_normals=outward)
println("panels: ", length(panels), " facets, surface ", round(sum(panels.areas); digits=2), " m², extent ", round(panels.extent_m; digits=2),
    " m, reference area ", round(panels.reference_area_m2; digits=3), " m², reference length ", round(panels.reference_length_m; digits=3), " m")
for (name, v) in (("+x", SVector(1.0, 0.0, 0.0)), ("+y", SVector(0.0, 1.0, 0.0)), ("+z", SVector(0.0, 0.0, 1.0)))
    println("  projected area along ", name, ": ", round(panel_projected_area(panels, v; grid=grid); digits=3), " m²")
end
sur = fit_mesh_aero_surrogate(panels; degree=degree, poly_degree=poly, n_directions=ndir, speed_ratios=srs, sigma_n=sigma_n, sigma_t=sigma_t,
    shadowing=shadowing, grid=grid, verbose=true)
hold = sur.metadata["fit"]["holdout"]
println("holdout at s=", hold["speed_ratio"], ":")
for c in ("CFx", "CFy", "CFz", "CMx", "CMy", "CMz")
    println("  ", c, ": max error ", round(hold[c]["max_error"]; sigdigits=3), " of scale ", round(hold[c]["scale"]; sigdigits=3))
end
s_mid = srs[(length(srs) + 1) ÷ 2]
for (name, v) in (("+x", SVector(1.0, 0.0, 0.0)), ("-x", SVector(-1.0, 0.0, 0.0)), ("+y", SVector(0.0, 1.0, 0.0)), ("+z", SVector(0.0, 0.0, 1.0)))
    cf, cm = mesh_aero_coefficients(sur, v, s_mid)
    println("  flow along ", name, " at s=", s_mid, ": CD=", round(-sum(cf .* v); digits=4), "  CF=", round.(cf; digits=4), "  CM=", round.(cm; digits=4))
end
write_mesh_aero_surrogate(out, sur)
println("wrote ", out, " in ", round(time() - t0; digits=1), " s")
