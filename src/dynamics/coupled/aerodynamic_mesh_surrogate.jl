# Free-molecular aerodynamics of an arbitrary triangle mesh (a CAD model)
# and its analytic surrogate.
#
# Offline: `mesh_aero_panels` reads a model (through the structure layer's
# mesh readers, with the same scale/rotation/centring the viewer applies),
# `panel_aero_coefficients` integrates the Schaaf-Chambre pressure and shear
# over the facets for one flow direction and speed ratio with depth-buffer
# self-shadowing, and `fit_mesh_aero_surrogate` tabulates those over the
# sphere of flow directions and a range of speed ratios and fits real
# spherical harmonics (in direction) times a polynomial (in 1/s) to the six
# body-frame force and moment coefficients. The diffuse re-emission term is
# the only place the wall temperature enters, and it enters linearly through
# sqrt(Tw/T), so each coefficient is fitted as A(v, s) + sqrt(Tw/T) B(v, s)
# and the wall temperature stays a free runtime parameter.
#
# Runtime: `AerodynamicCoefficientMeshSurrogate` holds one surrogate per
# link and evaluates it from the link-frame airspeed direction, the speed
# ratio and the wall temperature, so a fitted surrogate costs a few hundred
# multiply-adds per link per RHS call, no mesh in the loop.
#
# Conventions: coefficients are normalised by the dynamic pressure
# 0.5 rho V^2 of the wind-relative airspeed, `reference_area_m2` and (for
# moments) `reference_length_m`; the force coefficient is in the mesh's
# (link) frame, the moment about `moment_reference_m` in that frame. `vhat`
# is the unit vector of the body's velocity relative to the gas, expressed
# in the link frame. A facet whose outward normal has a positive projection
# on `vhat` faces the flow.

using JSON

# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------

"""
    MeshAeroPanels

A triangle mesh prepared for the panel method: vertices (3 x 3N, metres),
per-facet centroids, unit normals (from the vertex winding) and areas, plus
the normalisation the coefficients use. Build with [`mesh_aero_panels`](@ref).
"""
struct MeshAeroPanels
    vertices::Matrix{Float64}
    centroids::Matrix{Float64}
    normals::Matrix{Float64}
    areas::Vector{Float64}
    reference_area_m2::Float64
    reference_length_m::Float64
    moment_reference_m::SVector{3, Float64}
    extent_m::Float64
    source::String
end

Base.length(p::MeshAeroPanels) = length(p.areas)

"""
    mesh_aero_panels(path; scale=1.0, rotation_deg=(0, 0, 0), center=true,
                     reference_area_m2=nothing, reference_length_m=nothing,
                     moment_reference_m=nothing, outward_normals=false,
                     min_area_m2=0.0) -> MeshAeroPanels
    mesh_aero_panels(triangles::AbstractMatrix; kwargs...)

Prepare a model file (STL, OBJ, glTF/GLB) or a 3 x 3N triangle matrix for the
panel method. `scale`, `rotation_deg`, `center` and `articulations` are the
viewer's model transform (metres per model unit, XYZ Euler angles in
degrees, bounding-box centre to the origin, and parts rotated about an axis
in model units, see `articulate_triangles`), so the panels sit exactly where
the page draws the model and the link frame is the model frame after that
transform.

Normals come from the vertex winding (counter-clockwise seen from outside,
as glTF and STL require); `outward_normals=true` instead flips every normal
that points toward the mesh's centroid, a heuristic for meshes with mixed
winding that is exact for convex bodies and wrong on concave parts such as
the inside of a dish. Facets below `min_area_m2` are dropped.

`reference_area_m2` defaults to the mean projected area, total surface area
over four (exact for convex bodies, Cauchy's theorem); `reference_length_m`
to the bounding-box diagonal; `moment_reference_m` to the origin (the link's
centre of mass once the model is centred).
"""
function mesh_aero_panels(
    path::AbstractString;
    scale::Real=1.0,
    rotation_deg=(0.0, 0.0, 0.0),
    center::Bool=true,
    articulations=(),
    kwargs...
)::MeshAeroPanels
    tris = Structure.load_model_triangles(path; scale=scale, rotation_deg=rotation_deg, articulations=articulations)
    if center
        c = SVector{3, Float64}(
            0.5 * (minimum(tris[1, :]) + maximum(tris[1, :])),
            0.5 * (minimum(tris[2, :]) + maximum(tris[2, :])),
            0.5 * (minimum(tris[3, :]) + maximum(tris[3, :])))
        tris = tris .- c
    end
    return mesh_aero_panels(tris; source=basename(String(path)), kwargs...)
end

function mesh_aero_panels(
    triangles::AbstractMatrix{<:Real};
    reference_area_m2::Union{Nothing, Real}=nothing,
    reference_length_m::Union{Nothing, Real}=nothing,
    moment_reference_m=nothing,
    outward_normals::Bool=false,
    min_area_m2::Real=0.0,
    source::AbstractString="triangles"
)::MeshAeroPanels
    size(triangles, 1) == 3 || throw(ArgumentError("triangles must be 3 x 3N, got $(size(triangles))."))
    size(triangles, 2) % 3 == 0 || throw(ArgumentError("triangles must hold whole triangles (3 columns each), got $(size(triangles, 2)) columns."))
    ntri = size(triangles, 2) ÷ 3
    ntri >= 1 || throw(ArgumentError("the mesh has no triangles."))
    V = Matrix{Float64}(triangles)
    keep = Int[]
    cents = SVector{3, Float64}[]
    norms = SVector{3, Float64}[]
    areas = Float64[]
    @inbounds for t in 1:ntri
        a = SVector{3, Float64}(V[1, 3t - 2], V[2, 3t - 2], V[3, 3t - 2])
        b = SVector{3, Float64}(V[1, 3t - 1], V[2, 3t - 1], V[3, 3t - 1])
        c = SVector{3, Float64}(V[1, 3t], V[2, 3t], V[3, 3t])
        n2 = cross(b - a, c - a)
        twice = norm(n2)
        area = 0.5 * twice
        (isfinite(area) && area > max(Float64(min_area_m2), 0.0)) || continue
        push!(keep, t)
        push!(cents, (a + b + c) / 3)
        push!(norms, n2 / twice)
        push!(areas, area)
    end
    isempty(keep) && throw(ArgumentError("the mesh has no facets above min_area_m2 = $(min_area_m2)."))
    if outward_normals
        mesh_center = sum(cents[k] * areas[k] for k in eachindex(areas)) / sum(areas)
        @inbounds for k in eachindex(norms)
            dot(norms[k], cents[k] - mesh_center) < 0.0 && (norms[k] = -norms[k])
        end
    end
    n = length(keep)
    vertices = Matrix{Float64}(undef, 3, 3n)
    centroids = Matrix{Float64}(undef, 3, n)
    normals = Matrix{Float64}(undef, 3, n)
    @inbounds for (k, t) in enumerate(keep)
        for v in 1:3, c in 1:3
            vertices[c, 3k - 3 + v] = V[c, 3t - 3 + v]
        end
        for c in 1:3
            centroids[c, k] = cents[k][c]
            normals[c, k] = norms[k][c]
        end
    end
    lo = SVector{3, Float64}(minimum(vertices[1, :]), minimum(vertices[2, :]), minimum(vertices[3, :]))
    hi = SVector{3, Float64}(maximum(vertices[1, :]), maximum(vertices[2, :]), maximum(vertices[3, :]))
    extent = max(norm(hi - lo), eps(Float64))
    total_area = sum(areas)
    a_ref = reference_area_m2 === nothing ? total_area / 4 : Float64(reference_area_m2)
    a_ref > 0.0 || throw(ArgumentError("reference_area_m2 must be positive."))
    l_ref = reference_length_m === nothing ? extent : Float64(reference_length_m)
    l_ref > 0.0 || throw(ArgumentError("reference_length_m must be positive."))
    m_ref = moment_reference_m === nothing ? SVector{3, Float64}(0.0, 0.0, 0.0) : SVector{3, Float64}(moment_reference_m)
    return MeshAeroPanels(vertices, centroids, normals, areas, a_ref, l_ref, m_ref, extent, String(source))
end

# Orthonormal pair spanning the plane perpendicular to a unit vector.
@inline function _plane_basis(v::SVector{3, Float64})
    helper = abs(v[3]) < 0.9 ? SVector{3, Float64}(0.0, 0.0, 1.0) : SVector{3, Float64}(1.0, 0.0, 0.0)
    e1 = normalize(cross(helper, v))
    e2 = cross(v, e1)
    return e1, e2
end

"""
    panel_shadow_mask(panels, vhat; grid=256) -> BitVector

`true` for every facet whose centroid lies behind another facet along the
flow direction `vhat` (unit, link frame): the mesh is rasterised into a
`grid x grid` depth buffer on the plane perpendicular to `vhat`, each pixel
keeping the most upstream facet, and a facet is shadowed when the pixel under
its centroid belongs to a different facet that sits upstream by more than
the pixel's own depth uncertainty. Facets facing away from the flow are
tested the same way; the free-molecular formulas already give them only the
thermal (exponentially small) flux. Facets exactly edge-on to the flow are
never shadowed: they carry no pressure, and the flow runs along both sides
of a thin body (so a box in axial flow gets shear on all four sides, where
the Hart closed forms count one face of each pair).
"""
function panel_shadow_mask(panels::MeshAeroPanels, vhat::SVector{3, Float64}; grid::Integer=256)::BitVector
    n = length(panels)
    mask = falses(n)
    n <= 1 && return mask
    g = max(Int(grid), 8)
    e1, e2 = _plane_basis(vhat)
    V = panels.vertices
    nv = size(V, 2)
    u = Vector{Float64}(undef, nv); w = Vector{Float64}(undef, nv); d = Vector{Float64}(undef, nv)
    umin = Inf; umax = -Inf; wmin = Inf; wmax = -Inf
    @inbounds for i in 1:nv
        p = SVector{3, Float64}(V[1, i], V[2, i], V[3, i])
        u[i] = dot(p, e1); w[i] = dot(p, e2); d[i] = dot(p, vhat)
        umin = min(umin, u[i]); umax = max(umax, u[i]); wmin = min(wmin, w[i]); wmax = max(wmax, w[i])
    end
    span = max(umax - umin, wmax - wmin, eps(Float64))
    pixel = span / (g - 2)
    # centre the footprint in the buffer with a one-pixel margin
    u0 = umin - 0.5 * (pixel * (g - 2) - (umax - umin)) - pixel
    w0 = wmin - 0.5 * (pixel * (g - 2) - (wmax - wmin)) - pixel
    depth = fill(-Inf, g, g)
    owner = zeros(Int32, g, g)
    @inbounds for t in 1:n
        i1, i2, i3 = 3t - 2, 3t - 1, 3t
        ua, wa, da = u[i1], w[i1], d[i1]
        ub, wb, db = u[i2], w[i2], d[i2]
        uc, wc, dc = u[i3], w[i3], d[i3]
        det = (ub - ua) * (wc - wa) - (uc - ua) * (wb - wa)
        abs(det) <= eps(Float64) * span^2 && continue   # edge-on to the flow
        inv_det = 1.0 / det
        px_lo = max(1, floor(Int, (min(ua, ub, uc) - u0) / pixel) + 1)
        px_hi = min(g, floor(Int, (max(ua, ub, uc) - u0) / pixel) + 1)
        py_lo = max(1, floor(Int, (min(wa, wb, wc) - w0) / pixel) + 1)
        py_hi = min(g, floor(Int, (max(wa, wb, wc) - w0) / pixel) + 1)
        for py in py_lo:py_hi
            wp = w0 + (py - 0.5) * pixel
            for px in px_lo:px_hi
                up = u0 + (px - 0.5) * pixel
                # barycentric coordinates of the pixel centre
                l2 = ((up - ua) * (wc - wa) - (uc - ua) * (wp - wa)) * inv_det
                l3 = ((ub - ua) * (wp - wa) - (up - ua) * (wb - wa)) * inv_det
                l1 = 1.0 - l2 - l3
                (l1 >= -1e-9 && l2 >= -1e-9 && l3 >= -1e-9) || continue
                dp = l1 * da + l2 * db + l3 * dc
                if dp > depth[px, py]
                    depth[px, py] = dp
                    owner[px, py] = Int32(t)
                end
            end
        end
    end
    C = panels.centroids
    N = panels.normals
    @inbounds for t in 1:n
        c = SVector{3, Float64}(C[1, t], C[2, t], C[3, t])
        px = clamp(floor(Int, (dot(c, e1) - u0) / pixel) + 1, 1, g)
        py = clamp(floor(Int, (dot(c, e2) - w0) / pixel) + 1, 1, g)
        o = owner[px, py]
        (o == 0 || o == t) && continue
        gamma = abs(N[1, t] * vhat[1] + N[2, t] * vhat[2] + N[3, t] * vhat[3])
        # An edge-on facet carries no pressure and its shear acts on a surface
        # the flow runs along, not one it has to reach: leave it unshadowed.
        gamma <= 1e-9 && continue
        # depth uncertainty of the sample under this facet: one pixel times its slope along the flow
        slope = sqrt(max(0.0, 1.0 - gamma^2)) / max(gamma, 0.1)
        tol = pixel * (1.0 + slope) + 1e-9 * panels.extent_m
        depth[px, py] > dot(c, vhat) + tol && (mask[t] = true)
    end
    return mask
end

# Schaaf-Chambre pressure and shear on one facet: gamma = n . vhat (cosine of
# the angle between the outward normal and the body's velocity relative to
# the gas), s the speed ratio. Returns (Cp_A, Cp_B, Ct) with
# Cp = Cp_A + sqrt(Tw/T) Cp_B and the shear coefficient Ct multiplying the
# unnormalised tangential direction (vhat - gamma n), all per unit area and
# per dynamic pressure. Every term is finite for gamma in [-1, 1].
@inline function _schaaf_chambre(gamma::Float64, s::Float64, sigma_n::Float64, sigma_t::Float64)
    sg = s * gamma
    E = exp(-sg * sg)
    F = 1.0 + erf(sg)
    inv_s2 = 1.0 / (s * s)
    cp_a = inv_s2 * ((2.0 - sigma_n) * inv_sqrt_π * sg * E + (2.0 - sigma_n) * (sg * sg + 0.5) * F)
    cp_b = inv_s2 * (0.5 * sigma_n * E + 0.5 * sigma_n * sqrt_π * sg * F)
    # shear: sigma_t sin(theta)/(s sqrt(pi)) [E + sqrt(pi) s gamma F], and the
    # tangential unit vector (vhat - gamma n)/sin(theta); the sines cancel.
    ct = sigma_t * inv_sqrt_π / s * (E + sqrt_π * sg * F)
    return cp_a, cp_b, ct
end

"""
    panel_aero_coefficients_split(panels, vhat, s; sigma_n=1.0, sigma_t=1.0, shadowing=true, grid=256)
        -> (CF_A, CM_A, CF_B, CM_B)

Body-frame force and moment coefficients of the mesh for the airspeed
direction `vhat` (unit, link frame) and speed ratio `s`, split into the part
independent of the wall temperature and the part multiplying
sqrt(Tw / T_inf): `CF = CF_A + sqrt(Tw/T) CF_B`, likewise for `CM`. The force
is normalised by `reference_area_m2`, the moment (about
`moment_reference_m`) by `reference_area_m2 * reference_length_m`.
`sigma_n`, `sigma_t` are the normal and tangential momentum accommodation
coefficients (1 = fully diffuse). With `shadowing` on, facets behind others
along the flow contribute nothing.
"""
function panel_aero_coefficients_split(
    panels::MeshAeroPanels,
    vhat::SVector{3, Float64},
    s::Real;
    sigma_n::Real=1.0,
    sigma_t::Real=1.0,
    shadowing::Bool=true,
    grid::Integer=256
)
    sval = Float64(s)
    sval > 0.0 || throw(ArgumentError("the speed ratio must be positive, got $(s)."))
    v = normalize(vhat)
    sn, st = Float64(sigma_n), Float64(sigma_t)
    mask = shadowing ? panel_shadow_mask(panels, v; grid=grid) : falses(length(panels))
    N = panels.normals
    C = panels.centroids
    fa = MVector{3, Float64}(0.0, 0.0, 0.0); fb = MVector{3, Float64}(0.0, 0.0, 0.0)
    ma = MVector{3, Float64}(0.0, 0.0, 0.0); mb = MVector{3, Float64}(0.0, 0.0, 0.0)
    r0 = panels.moment_reference_m
    @inbounds for t in 1:length(panels)
        mask[t] && continue
        n = SVector{3, Float64}(N[1, t], N[2, t], N[3, t])
        gamma = dot(n, v)
        cp_a, cp_b, ct = _schaaf_chambre(gamma, sval, sn, st)
        area = panels.areas[t]
        tangent = v - gamma * n
        # pressure pushes along -n, shear along the tangential direction of the body's motion reversed
        f_a = area * (-cp_a * n - ct * tangent)
        f_b = area * (-cp_b * n)
        lever = SVector{3, Float64}(C[1, t], C[2, t], C[3, t]) - r0
        fa .+= f_a; fb .+= f_b
        ma .+= cross(lever, f_a); mb .+= cross(lever, f_b)
    end
    inv_a = 1.0 / panels.reference_area_m2
    inv_al = inv_a / panels.reference_length_m
    return SVector{3, Float64}(fa) * inv_a, SVector{3, Float64}(ma) * inv_al, SVector{3, Float64}(fb) * inv_a, SVector{3, Float64}(mb) * inv_al
end

"""
    panel_aero_coefficients(panels, vhat, s; tw_ratio=1.0, kwargs...) -> (CF, CM)

Force and moment coefficients for the wall-to-freestream temperature ratio
`tw_ratio`; see [`panel_aero_coefficients_split`](@ref) for the conventions.
"""
function panel_aero_coefficients(panels::MeshAeroPanels, vhat::SVector{3, Float64}, s::Real; tw_ratio::Real=1.0, kwargs...)
    tw_ratio >= 0.0 || throw(ArgumentError("tw_ratio must be non-negative."))
    cf_a, cm_a, cf_b, cm_b = panel_aero_coefficients_split(panels, vhat, s; kwargs...)
    k = sqrt(Float64(tw_ratio))
    return cf_a + k * cf_b, cm_a + k * cm_b
end

"""
    panel_projected_area(panels, vhat; grid=512) -> Float64

Projected (silhouette) area of the mesh along `vhat`, from the depth buffer.
"""
function panel_projected_area(panels::MeshAeroPanels, vhat::SVector{3, Float64}; grid::Integer=512)::Float64
    v = normalize(vhat)
    mask = panel_shadow_mask(panels, v; grid=grid)
    total = 0.0
    N = panels.normals
    @inbounds for t in 1:length(panels)
        mask[t] && continue
        gamma = N[1, t] * v[1] + N[2, t] * v[2] + N[3, t] * v[3]
        gamma > 0.0 && (total += gamma * panels.areas[t])
    end
    return total
end

# ---------------------------------------------------------------------------
# Real spherical harmonics (fully normalised, no Condon-Shortley phase)
# ---------------------------------------------------------------------------

@inline _sh_count(degree::Int)::Int = (degree + 1)^2
@inline _sh_index(l::Int, m::Int)::Int = l * l + l + m + 1

# Fills `Y[1:(L+1)^2]` with Y_lm(v) for l = 0..L, m = -l..l, indexed by
# `_sh_index`. Standard three-term recursion on the normalised associated
# Legendre functions; exact for any unit vector, including the poles.
function real_sh_basis!(Y::AbstractVector{Float64}, degree::Int, v::SVector{3, Float64})
    L = degree
    x, y, z = v
    ct = clamp(z, -1.0, 1.0)
    st = sqrt(max(0.0, 1.0 - ct * ct))
    phi = atan(y, x)
    # P[m+1, l+1] holds the normalised P_l^m
    P = zeros(Float64, L + 1, L + 1)
    P[1, 1] = 1.0 / sqrt(4pi)
    @inbounds for m in 1:L
        P[m + 1, m + 1] = -sqrt((2m + 1) / (2m)) * st * P[m, m]
    end
    @inbounds for m in 0:(L - 1)
        P[m + 1, m + 2] = sqrt(2m + 3.0) * ct * P[m + 1, m + 1]
    end
    @inbounds for m in 0:L, l in (m + 2):L
        a = sqrt((4.0 * l * l - 1.0) / (l * l - m * m))
        b = sqrt(((l - 1.0)^2 - m * m) / (4.0 * (l - 1.0)^2 - 1.0))
        P[m + 1, l + 1] = a * (ct * P[m + 1, l] - b * P[m + 1, l - 1])
    end
    @inbounds for l in 0:L
        Y[_sh_index(l, 0)] = P[1, l + 1]
        for m in 1:l
            c = sqrt(2.0) * P[m + 1, l + 1]
            Y[_sh_index(l, m)] = c * cos(m * phi)
            Y[_sh_index(l, -m)] = c * sin(m * phi)
        end
    end
    return Y
end

# Basis of the full fit: Y_lm(v) * x^k with x = 1/s, k = 0..P, laid out as
# [k = 0 block][k = 1 block]... so the direction basis is reused.
@inline _basis_count(degree::Int, poly_degree::Int)::Int = _sh_count(degree) * (poly_degree + 1)

function _fit_basis!(row::AbstractVector{Float64}, Y::AbstractVector{Float64}, degree::Int, poly_degree::Int, v::SVector{3, Float64}, s::Float64)
    nsh = _sh_count(degree)
    real_sh_basis!(Y, degree, v)
    x = 1.0 / s
    xk = 1.0
    @inbounds for k in 0:poly_degree
        off = k * nsh
        for j in 1:nsh
            row[off + j] = Y[j] * xk
        end
        xk *= x
    end
    return row
end

# Fibonacci sphere: n well-spread unit vectors.
function fibonacci_directions(n::Integer)::Vector{SVector{3, Float64}}
    n >= 1 || throw(ArgumentError("n must be >= 1."))
    golden = pi * (3.0 - sqrt(5.0))
    out = Vector{SVector{3, Float64}}(undef, Int(n))
    @inbounds for i in 1:Int(n)
        z = 1.0 - 2.0 * (i - 0.5) / n
        r = sqrt(max(0.0, 1.0 - z * z))
        phi = golden * (i - 1)
        out[i] = SVector{3, Float64}(r * cos(phi), r * sin(phi), z)
    end
    return out
end

# ---------------------------------------------------------------------------
# Surrogate
# ---------------------------------------------------------------------------

"""
    MeshAeroSurrogate

Analytic approximation of a mesh's free-molecular force and moment
coefficients: for each of the six components, real spherical harmonics of
degree `degree` in the link-frame airspeed direction times a polynomial of
degree `poly_degree` in the inverse speed ratio, one set (`coeff_a`) for the
wall-temperature-independent part and one (`coeff_b`) multiplying
sqrt(Tw / T_inf). Evaluate with [`mesh_aero_coefficients`](@ref); build with
[`fit_mesh_aero_surrogate`](@ref); persist with
[`write_mesh_aero_surrogate`](@ref) / [`read_mesh_aero_surrogate`](@ref).
`metadata` records the source mesh, the accommodation coefficients, the
sampling and the fit residuals.
"""
struct MeshAeroSurrogate
    degree::Int
    poly_degree::Int
    coeff_a::Matrix{Float64}   # 6 x nbasis, rows CFx CFy CFz CMx CMy CMz
    coeff_b::Matrix{Float64}
    reference_area_m2::Float64
    reference_length_m::Float64
    moment_reference_m::SVector{3, Float64}
    sigma_n::Float64
    sigma_t::Float64
    speed_ratio_min::Float64
    speed_ratio_max::Float64
    metadata::Dict{String, Any}
end

function MeshAeroSurrogate(degree::Integer, poly_degree::Integer, coeff_a::AbstractMatrix, coeff_b::AbstractMatrix,
                           reference_area_m2::Real, reference_length_m::Real, moment_reference_m,
                           sigma_n::Real, sigma_t::Real, speed_ratio_min::Real, speed_ratio_max::Real,
                           metadata::AbstractDict=Dict{String, Any}())
    nb = _basis_count(Int(degree), Int(poly_degree))
    size(coeff_a) == (6, nb) || throw(ArgumentError("coeff_a must be 6 x $(nb) for degree $(degree), poly_degree $(poly_degree); got $(size(coeff_a))."))
    size(coeff_b) == (6, nb) || throw(ArgumentError("coeff_b must be 6 x $(nb); got $(size(coeff_b))."))
    reference_area_m2 > 0 || throw(ArgumentError("reference_area_m2 must be positive."))
    reference_length_m > 0 || throw(ArgumentError("reference_length_m must be positive."))
    0 < speed_ratio_min <= speed_ratio_max || throw(ArgumentError("need 0 < speed_ratio_min <= speed_ratio_max."))
    return MeshAeroSurrogate(Int(degree), Int(poly_degree), Matrix{Float64}(coeff_a), Matrix{Float64}(coeff_b),
        Float64(reference_area_m2), Float64(reference_length_m), SVector{3, Float64}(moment_reference_m),
        Float64(sigma_n), Float64(sigma_t), Float64(speed_ratio_min), Float64(speed_ratio_max),
        Dict{String, Any}(String(k) => v for (k, v) in metadata))
end

"""
    mesh_aero_coefficients(surrogate, vhat, s; tw_ratio=1.0) -> (CF, CM)

Force and moment coefficients (link frame; see
[`panel_aero_coefficients_split`](@ref)) from the fitted surrogate for the
airspeed direction `vhat`, speed ratio `s` (clamped to the fitted range) and
wall-to-freestream temperature ratio `tw_ratio`.
"""
function mesh_aero_coefficients(sur::MeshAeroSurrogate, vhat::SVector{3, Float64}, s::Real; tw_ratio::Real=1.0)
    sc = clamp(Float64(s), sur.speed_ratio_min, sur.speed_ratio_max)
    v = normalize(vhat)
    nsh = _sh_count(sur.degree)
    Y = MVector{441, Float64}(undef)   # degree <= 20
    real_sh_basis!(Y, sur.degree, v)
    x = 1.0 / sc
    k_tw = sqrt(max(Float64(tw_ratio), 0.0))
    out = MVector{6, Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    A = sur.coeff_a; B = sur.coeff_b
    xk = 1.0
    @inbounds for k in 0:sur.poly_degree
        off = k * nsh
        for j in 1:nsh
            w = Y[j] * xk
            for c in 1:6
                out[c] += w * (A[c, off + j] + k_tw * B[c, off + j])
            end
        end
        xk *= x
    end
    return SVector{3, Float64}(out[1], out[2], out[3]), SVector{3, Float64}(out[4], out[5], out[6])
end

"""
    fit_mesh_aero_surrogate(panels; degree=8, poly_degree=2, n_directions=1200,
                            speed_ratios=(3, 5, 8, 12, 20), sigma_n=1.0, sigma_t=1.0,
                            shadowing=true, grid=256, holdout_directions=200, verbose=false)
        -> MeshAeroSurrogate

Tabulate the panel-method coefficients of `panels` over `n_directions`
Fibonacci-spread airspeed directions and the given speed ratios, and fit the
surrogate by least squares. The fit residual (rms and max absolute error per
coefficient, over the training samples and over `holdout_directions` fresh
directions at the middle speed ratio) is stored in the metadata under
`"fit"`; compare it with the coefficient magnitudes there before trusting a
low degree. Cost is one shadow buffer per direction and speed ratio pair.
"""
function fit_mesh_aero_surrogate(
    panels::MeshAeroPanels;
    degree::Integer=8,
    poly_degree::Integer=2,
    n_directions::Integer=1200,
    speed_ratios=(3.0, 5.0, 8.0, 12.0, 20.0),
    sigma_n::Real=1.0,
    sigma_t::Real=1.0,
    shadowing::Bool=true,
    grid::Integer=256,
    holdout_directions::Integer=200,
    verbose::Bool=false
)::MeshAeroSurrogate
    L = Int(degree); P = Int(poly_degree)
    0 <= L <= 20 || throw(ArgumentError("degree must be within 0..20."))
    P >= 0 || throw(ArgumentError("poly_degree must be >= 0."))
    srs = sort!(Float64[Float64(s) for s in speed_ratios])
    isempty(srs) && throw(ArgumentError("speed_ratios is empty."))
    all(>(0.0), srs) || throw(ArgumentError("speed ratios must be positive."))
    dirs = fibonacci_directions(n_directions)
    nb = _basis_count(L, P)
    nsamp = length(dirs) * length(srs)
    nsamp >= nb || throw(ArgumentError("$(nsamp) samples cannot determine $(nb) basis functions; raise n_directions or lower degree."))
    X = Matrix{Float64}(undef, nsamp, nb)
    YA = Matrix{Float64}(undef, nsamp, 6)
    YB = Matrix{Float64}(undef, nsamp, 6)
    Ybuf = Vector{Float64}(undef, _sh_count(L))
    row = Vector{Float64}(undef, nb)
    started = time_ns()
    k = 0
    for s in srs, v in dirs
        k += 1
        cf_a, cm_a, cf_b, cm_b = panel_aero_coefficients_split(panels, v, s; sigma_n=sigma_n, sigma_t=sigma_t, shadowing=shadowing, grid=grid)
        _fit_basis!(row, Ybuf, L, P, v, s)
        @inbounds X[k, :] .= row
        @inbounds YA[k, :] .= (cf_a[1], cf_a[2], cf_a[3], cm_a[1], cm_a[2], cm_a[3])
        @inbounds YB[k, :] .= (cf_b[1], cf_b[2], cf_b[3], cm_b[1], cm_b[2], cm_b[3])
        if verbose && k % 500 == 0
            println("  tabulated ", k, " / ", nsamp, " samples (", round((time_ns() - started) / 1e9; digits=1), " s)")
        end
    end
    tab_s = (time_ns() - started) / 1e9
    CA = Matrix{Float64}((X \ YA)')
    CB = Matrix{Float64}((X \ YB)')
    # residuals
    RA = X * CA' .- YA
    RB = X * CB' .- YB
    names = ("CFx", "CFy", "CFz", "CMx", "CMy", "CMz")
    fit = Dict{String, Any}(
        "samples" => nsamp,
        "basis" => nb,
        "tabulation_s" => tab_s,
        "train" => Dict{String, Any}(names[c] => Dict("rms_a" => sqrt(sum(abs2, RA[:, c]) / nsamp), "max_a" => maximum(abs, RA[:, c]),
                                                    "rms_b" => sqrt(sum(abs2, RB[:, c]) / nsamp), "max_b" => maximum(abs, RB[:, c]),
                                                    "scale_a" => maximum(abs, YA[:, c]), "scale_b" => maximum(abs, YB[:, c])) for c in 1:6),
    )
    sur = MeshAeroSurrogate(L, P, CA, CB, panels.reference_area_m2, panels.reference_length_m, panels.moment_reference_m,
        sigma_n, sigma_t, srs[1], srs[end], Dict{String, Any}())
    if holdout_directions > 0
        s_mid = srs[(length(srs) + 1) ÷ 2]
        hold = fibonacci_directions(Int(holdout_directions) + 7)[8:end]   # offset so they differ from the training set
        err = zeros(6); scale = zeros(6)
        for v in hold
            cf_a, cm_a, cf_b, cm_b = panel_aero_coefficients_split(panels, v, s_mid; sigma_n=sigma_n, sigma_t=sigma_t, shadowing=shadowing, grid=grid)
            truth = (cf_a + cf_b, cm_a + cm_b)   # Tw = T
            cf, cm = mesh_aero_coefficients(sur, v, s_mid; tw_ratio=1.0)
            for c in 1:3
                err[c] = max(err[c], abs(cf[c] - truth[1][c])); scale[c] = max(scale[c], abs(truth[1][c]))
                err[c + 3] = max(err[c + 3], abs(cm[c] - truth[2][c])); scale[c + 3] = max(scale[c + 3], abs(truth[2][c]))
            end
        end
        fit["holdout"] = Dict{String, Any}("speed_ratio" => s_mid, "directions" => length(hold),
            (names[c] => Dict("max_error" => err[c], "scale" => scale[c]) for c in 1:6)...)
    end
    meta = Dict{String, Any}(
        "source" => panels.source, "facets" => length(panels), "surface_area_m2" => sum(panels.areas), "extent_m" => panels.extent_m,
        "n_directions" => length(dirs), "speed_ratios" => srs, "shadowing" => shadowing, "grid" => Int(grid),
        "sigma_n" => Float64(sigma_n), "sigma_t" => Float64(sigma_t), "fit" => fit,
    )
    if verbose
        println("fit: ", nb, " basis functions from ", nsamp, " samples in ", round(tab_s; digits=1), " s")
        for c in 1:6
            t = fit["train"][names[c]]
            println("  ", names[c], ": rms ", round(t["rms_a"]; sigdigits=3), " of scale ", round(t["scale_a"]; sigdigits=3))
        end
    end
    return MeshAeroSurrogate(sur.degree, sur.poly_degree, sur.coeff_a, sur.coeff_b, sur.reference_area_m2, sur.reference_length_m,
        sur.moment_reference_m, sur.sigma_n, sur.sigma_t, sur.speed_ratio_min, sur.speed_ratio_max, meta)
end

"""
    write_mesh_aero_surrogate(path, surrogate) -> path

Save a surrogate as JSON (coefficients, normalisation, fit range, metadata).
"""
function write_mesh_aero_surrogate(path::AbstractString, sur::MeshAeroSurrogate)::String
    doc = Dict{String, Any}(
        "schema" => "spaceagora_mesh_aero_surrogate_v1",
        "degree" => sur.degree, "poly_degree" => sur.poly_degree,
        "coeff_a" => [sur.coeff_a[c, :] for c in 1:6], "coeff_b" => [sur.coeff_b[c, :] for c in 1:6],
        "reference_area_m2" => sur.reference_area_m2, "reference_length_m" => sur.reference_length_m,
        "moment_reference_m" => collect(sur.moment_reference_m),
        "sigma_n" => sur.sigma_n, "sigma_t" => sur.sigma_t,
        "speed_ratio_min" => sur.speed_ratio_min, "speed_ratio_max" => sur.speed_ratio_max,
        "metadata" => sur.metadata,
    )
    mkpath(dirname(abspath(String(path))))
    open(String(path), "w") do io
        JSON.print(io, doc)
    end
    return String(path)
end

"""
    read_mesh_aero_surrogate(path) -> MeshAeroSurrogate
"""
function read_mesh_aero_surrogate(path::AbstractString)::MeshAeroSurrogate
    doc = try
        JSON.parsefile(String(path))
    catch err
        throw(ArgumentError("$(path) is not a mesh aero surrogate file: $(sprint(showerror, err))"))
    end
    doc isa AbstractDict || throw(ArgumentError("$(path) is not a mesh aero surrogate file."))
    get(doc, "schema", "") == "spaceagora_mesh_aero_surrogate_v1" || throw(ArgumentError("$(path) is not a mesh aero surrogate file (schema $(repr(get(doc, "schema", nothing))))."))
    ca = Matrix{Float64}(reduce(vcat, [Float64.(r)' for r in doc["coeff_a"]]))
    cb = Matrix{Float64}(reduce(vcat, [Float64.(r)' for r in doc["coeff_b"]]))
    return MeshAeroSurrogate(Int(doc["degree"]), Int(doc["poly_degree"]), ca, cb,
        Float64(doc["reference_area_m2"]), Float64(doc["reference_length_m"]), Float64.(doc["moment_reference_m"]),
        Float64(doc["sigma_n"]), Float64(doc["sigma_t"]), Float64(doc["speed_ratio_min"]), Float64(doc["speed_ratio_max"]),
        Dict{String, Any}(get(doc, "metadata", Dict{String, Any}())))
end

# ---------------------------------------------------------------------------
# Runtime effector
# ---------------------------------------------------------------------------

"""
    AerodynamicCoefficientMeshSurrogate(surrogates::Dict{Int, MeshAeroSurrogate}; wall_temperature_k=300.0)
    AerodynamicCoefficientMeshSurrogate(surrogate::MeshAeroSurrogate; kwargs...)

Free-molecular aerodynamics from fitted mesh surrogates, one per link keyed
by the link's index in `spacecraft.links` (1 is the root); the single-argument
form puts one whole-vehicle surrogate on the root. Links without a surrogate
carry no aerodynamic load. Each link's coefficients are evaluated from the
airspeed direction in that link's frame, the speed ratio
V / sqrt(2 R T) of the wind-relative airspeed and the ratio of
`wall_temperature_k` to the local temperature, then scaled by the dynamic
pressure and the surrogate's reference area (and length for the moment).

Attitude: with `orientation_sim` the propagated root attitude and the stored
child attitudes place every link; without it the root is held in the
velocity-aligned frame (body x along the airspeed, body z toward nadir, the
attitude the viewer draws), child links relative to it, so a link's
configured quaternion selects its incidence as the fM model's `:attitude`
mode does. Torque about the root centre of mass is returned only when the
attitude is propagated, like the other aerodynamic effectors.
"""
struct AerodynamicCoefficientMeshSurrogate <: AbstractForceTorqueModel
    surrogates::Dict{Int, MeshAeroSurrogate}
    wall_temperature_k::Float64
    function AerodynamicCoefficientMeshSurrogate(surrogates::AbstractDict; wall_temperature_k::Real=300.0)
        isempty(surrogates) && throw(ArgumentError("AerodynamicCoefficientMeshSurrogate needs at least one surrogate."))
        wall_temperature_k > 0 || throw(ArgumentError("wall_temperature_k must be positive."))
        d = Dict{Int, MeshAeroSurrogate}()
        for (k, v) in surrogates
            Int(k) >= 1 || throw(ArgumentError("surrogate keys are 1-based link indices, got $(k)."))
            v isa MeshAeroSurrogate || throw(ArgumentError("surrogate for link $(k) must be a MeshAeroSurrogate."))
            d[Int(k)] = v
        end
        return new(d, Float64(wall_temperature_k))
    end
end

AerodynamicCoefficientMeshSurrogate(sur::MeshAeroSurrogate; kwargs...) = AerodynamicCoefficientMeshSurrogate(Dict(1 => sur); kwargs...)

@inline environment_requirements(::AerodynamicCoefficientMeshSurrogate) = EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)
@inline solver_partition(::AerodynamicCoefficientMeshSurrogate) = :implicit

# Velocity-aligned reference frame in planet-fixed coordinates: columns are
# the frame's x (along the airspeed), y and z (toward nadir) axes.
@inline function _velocity_aligned_columns(pos_pp::SVector{3, Float64}, vhat_pp::SVector{3, Float64})::SMatrix{3, 3, Float64, 9}
    nadir = -pos_pp
    z = nadir - dot(nadir, vhat_pp) * vhat_pp
    if norm(z) <= eps(Float64) * max(1.0, norm(pos_pp))
        helper = abs(vhat_pp[3]) < 0.9 ? SVector{3, Float64}(0.0, 0.0, 1.0) : SVector{3, Float64}(1.0, 0.0, 0.0)
        z = helper - dot(helper, vhat_pp) * vhat_pp
    end
    z = z / norm(z)
    y = cross(z, vhat_pp)
    return SMatrix{3, 3, Float64, 9}(vhat_pp[1], vhat_pp[2], vhat_pp[3], y[1], y[2], y[3], z[1], z[2], z[3])
end

# Shared evaluation. Returns (force_ii, torque_root, drag_ii, lift_ii, cross_ii)
# with the same contract as `_aero_pure_wrench`.
function _mesh_aero_wrench(
    model::AerodynamicCoefficientMeshSurrogate,
    spacecraft,
    q_ib::Union{Nothing, SVector{4, Float64}},
    planet,
    pos_pp::SVector{3, Float64},
    vel_pp::SVector{3, Float64},
    alt::Float64,
    lat::Float64,
    lon::Float64,
    l_pi::SMatrix{3, 3, Float64, 9},
    rho::Float64,
    T::Float64,
    wind,
)::NTuple{5, SVector{3, Float64}}
    if !isfinite(rho) || rho <= eps(Float64) || !isfinite(T) || T <= 0.0
        return _AERO_ZERO5
    end
    uD, uN, uE = latlongtoNED((alt, lat, lon))
    wE, wN, wU = wind
    wind_pp = wN * uN + wE * uE - wU * uD
    vel_rw = vel_pp - wind_pp
    speed = norm(vel_rw)
    speed <= eps(Float64) && return _AERO_ZERO5
    vhat_pp = vel_rw / speed
    q_dyn = 0.5 * rho * speed * speed
    s = speed / sqrt(2.0 * planet.R * T)
    tw_ratio = model.wall_temperature_k / T
    h_pp = cross(pos_pp, vel_pp)
    h_mag = norm(h_pp)
    (isfinite(h_mag) && h_mag > eps(Float64)) || return _AERO_ZERO5
    lift_pp_hat = normalize(cross(h_pp / h_mag, vhat_pp))
    drag_pp_hat = -vhat_pp
    cross_pp_hat = cross(drag_pp_hat, lift_pp_hat)
    l_pi_t = l_pi'
    orientation_sim = q_ib !== nothing
    vhat_pi = l_pi_t * vhat_pp
    R_root = orientation_sim ? SMatrix{3, 3, Float64, 9}(rot(q_ib)') : SMatrix{3, 3, Float64, 9}(I)
    C_pp_from_F = orientation_sim ? SMatrix{3, 3, Float64, 9}(I) : _velocity_aligned_columns(pos_pp, vhat_pp)
    e1 = SVector{3, Float64}(1.0, 0.0, 0.0)
    root = spacecraft.root
    rot_root_from_F = rot(SVector{4, Float64}(root.q...))

    force_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    torque_root = MVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for (k, body) in enumerate(spacecraft.links)
        sur = get(model.surrogates, k, nothing)
        sur === nothing && continue
        q_child = SVector{4, Float64}(body.q...)
        if orientation_sim
            R_link = body.root ? R_root : SMatrix{3, 3, Float64, 9}(R_root * rot(q_child)')   # link -> inertial
            vhat_link = R_link' * vhat_pi
        else
            rot_link_from_F = body.root ? rot_root_from_F : SMatrix{3, 3, Float64, 9}(rot(q_child) * rot_root_from_F)
            vhat_link = rot_link_from_F * e1
        end
        cf, cm = mesh_aero_coefficients(sur, vhat_link, s; tw_ratio=tw_ratio)
        f_link = (q_dyn * sur.reference_area_m2) * cf
        if orientation_sim
            f_ii = R_link * f_link
            m_link = (q_dyn * sur.reference_area_m2 * sur.reference_length_m) * cm
            f_root = rot(q_ib) * f_ii
            m_root = body.root ? m_link : rot(q_child)' * m_link
            lever_root = body.root ? SVector{3, Float64}(0.0, 0.0, 0.0) : SVector{3, Float64}(body.r)
            torque_root .+= m_root + cross(lever_root, f_root)
        else
            f_F = rot_link_from_F' * f_link
            f_ii = l_pi_t * (C_pp_from_F * f_F)
        end
        force_ii .+= f_ii
    end
    f = SVector{3, Float64}(force_ii)
    drag_ii_hat = l_pi_t * drag_pp_hat
    lift_ii_hat = l_pi_t * lift_pp_hat
    cross_ii_hat = l_pi_t * cross_pp_hat
    drag_ii = dot(f, drag_ii_hat) * drag_ii_hat
    lift_ii = dot(f, lift_ii_hat) * lift_ii_hat
    cross_ii = dot(f, cross_ii_hat) * cross_ii_hat
    return (f, SVector{3, Float64}(torque_root), drag_ii, lift_ii, cross_ii)
end

@inline function _mesh_aero_wrench_sampled(model::AerodynamicCoefficientMeshSurrogate, x::StateSample, env::EnvironmentSample)
    x.spacecraft === nothing && throw(ArgumentError("Mesh aerodynamic wrench evaluation requires StateSample.spacecraft."))
    pf = env.planet_frame
    at = env.atmosphere
    pf === nothing && throw(ArgumentError("Mesh aerodynamic wrench evaluation requires env.planet_frame."))
    at === nothing && throw(ArgumentError("Mesh aerodynamic wrench evaluation requires env.atmosphere."))
    return _mesh_aero_wrench(model, x.spacecraft, x.q_ib, env.planet, pf.pos_pp, pf.vel_pp, pf.alt_m, pf.lat_rad, pf.lon_rad,
        SMatrix{3, 3, Float64, 9}(pf.l_pi), at.rho_kg_m3, at.temperature_k, at.wind_pp)
end

@inline function wrench(model::AerodynamicCoefficientMeshSurrogate, x::StateSample, env::EnvironmentSample, t::Float64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    force, torque, _, _, _ = _mesh_aero_wrench_sampled(model, x, env)
    return force, torque
end

@inline function wrench_caching!(model::AerodynamicCoefficientMeshSurrogate, x::StateSample, env::EnvironmentSample, t::Float64, p::ODEParams, sat_idx::Int)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    force, torque, drag_ii, lift_ii, cross_ii = _mesh_aero_wrench_sampled(model, x, env)
    _store_aero_caches!(p, sat_idx, drag_ii, lift_ii, cross_ii)
    return force, torque
end

function calcForceTorque(model::AerodynamicCoefficientMeshSurrogate, x::AbstractVector{Float64}, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    planet = param.args.environment_model.planet
    orientation_sim = param.args.mission_configuration.orientation_sim
    spacecraft = param.args.dynamics_model.spacecraft[i]
    env_state = SimulationModel.SimulationCallbacks._stage_environment_state(x, param, i, param.shared_buffers.current_time[]; write_buffers=true)
    q_ib = if orientation_sim
        hasproperty(x, :q) ? SVector{4, Float64}(x.q) : SVector{4, Float64}(spacecraft.root.q...)
    else
        nothing
    end
    force, torque, drag_ii, lift_ii, cross_ii = _mesh_aero_wrench(model, spacecraft, q_ib, planet, env_state.pos_pp, env_state.vel_pp,
        env_state.alt, env_state.lat, env_state.lon, SMatrix{3, 3, Float64, 9}(env_state.l_pi), env_state.rho, env_state.T, env_state.wind)
    _store_aero_caches!(param, i, drag_ii, lift_ii, cross_ii)
    return force, torque
end
