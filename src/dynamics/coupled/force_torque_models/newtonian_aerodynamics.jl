"""
    modified_newtonian_cp_max(mach, gamma)

Return the stagnation-point pressure coefficient ``C_{p,max}`` used by
modified Newtonian theory. For supersonic flow, the calculation combines a
normal shock with isentropic deceleration to the stagnation point.

The `mach <= 1` branch is an isentropic, continuous numerical extension. It is
useful when a trajectory crosses Mach 1, but it does not extend the physical
validity of Newtonian surface aerodynamics outside hypersonic flow.
"""
@inline function modified_newtonian_cp_max(mach::Real, gamma::Real)
    mach_float = Float64(mach)
    gamma_float = Float64(gamma)
    isfinite(mach_float) && mach_float >= 0.0 || throw(
        DomainError(mach, "mach must be finite and nonnegative."),
    )
    isfinite(gamma_float) && gamma_float > 1.0 || throw(
        DomainError(gamma, "gamma must be finite and greater than 1."),
    )

    if mach_float == 0.0
        return 1.0
    elseif mach_float <= 1.0
        log_pressure_ratio = (gamma_float / (gamma_float - 1.0)) *
                             log1p(0.5 * (gamma_float - 1.0) * mach_float^2)
        return (2.0 / (gamma_float * mach_float^2)) * expm1(log_pressure_ratio)
    end

    mach_sq = mach_float^2
    normal_shock_total_pressure_factor = (
        ((gamma_float + 1.0)^2 * mach_sq) /
        (4.0 * gamma_float * mach_sq - 2.0 * (gamma_float - 1.0))
    )^(gamma_float / (gamma_float - 1.0))
    post_shock_static_pressure_ratio =
        (1.0 - gamma_float + 2.0 * gamma_float * mach_sq) / (gamma_float + 1.0)
    return (2.0 / (gamma_float * mach_sq)) *
           (normal_shock_total_pressure_factor * post_shock_static_pressure_ratio - 1.0)
end

"""
    AerodynamicSurfaceQuadrature

Regime-neutral continuous-surface quadrature samples. `positions_body_m` are
measured from the requested moment reference, `inward_normals_body` are unit
vectors, and `area_weights_m2` integrate the exact parametrized surface rather
than a faceted approximation. The same samples can be used by continuum,
free-molecular, and transitional aerodynamic laws.
"""
struct AerodynamicSurfaceQuadrature
    positions_body_m::Vector{SVector{3, Float64}}
    inward_normals_body::Vector{SVector{3, Float64}}
    area_weights_m2::Vector{Float64}

    function AerodynamicSurfaceQuadrature(
        positions_body_m::AbstractVector,
        inward_normals_body::AbstractVector,
        area_weights_m2::AbstractVector,
    )
        sample_count = length(positions_body_m)
        sample_count > 0 || throw(ArgumentError("aerodynamic surface quadrature must contain samples."))
        length(inward_normals_body) == sample_count || throw(DimensionMismatch(
            "surface position and normal counts must match.",
        ))
        length(area_weights_m2) == sample_count || throw(DimensionMismatch(
            "surface position and area-weight counts must match.",
        ))

        positions = Vector{SVector{3, Float64}}(undef, sample_count)
        normals = Vector{SVector{3, Float64}}(undef, sample_count)
        weights = Vector{Float64}(undef, sample_count)
        @inbounds for index in 1:sample_count
            position = SVector{3, Float64}(positions_body_m[index])
            normal = SVector{3, Float64}(inward_normals_body[index])
            weight = Float64(area_weights_m2[index])
            all(isfinite, position) || throw(DomainError(position, "surface positions must be finite."))
            all(isfinite, normal) || throw(DomainError(normal, "surface normals must be finite."))
            normal_magnitude = norm(normal)
            normal_magnitude > eps(Float64) || throw(DomainError(normal, "surface normals must be nonzero."))
            isfinite(weight) && weight > 0.0 || throw(DomainError(
                weight,
                "surface quadrature weights must be finite and positive.",
            ))
            positions[index] = position
            normals[index] = normal / normal_magnitude
            weights[index] = weight
        end
        return new(positions, normals, weights)
    end
end

const NewtonianSurfaceQuadrature = AerodynamicSurfaceQuadrature

"""
    AerodynamicGeometry(surface, reference_area_m2, reference_length_m; label="")

Geometry and normalization data for complete surface force and moment
integrals. The reference area and length use the coefficient normalization of
Grant and Braun Eqs. (4) and (5), independent of the surface pressure law.
"""
struct AerodynamicGeometry
    surface::AerodynamicSurfaceQuadrature
    reference_area_m2::Float64
    reference_length_m::Float64
    label::String

    function AerodynamicGeometry(
        surface::AerodynamicSurfaceQuadrature,
        reference_area_m2::Real,
        reference_length_m::Real;
        label::AbstractString="",
    )
        area = Float64(reference_area_m2)
        length = Float64(reference_length_m)
        isfinite(area) && area > 0.0 || throw(DomainError(
            reference_area_m2,
            "reference area must be finite and positive.",
        ))
        isfinite(length) && length > 0.0 || throw(DomainError(
            reference_length_m,
            "reference length must be finite and positive.",
        ))
        return new(surface, area, length, String(label))
    end
end

const NewtonianAerodynamicGeometry = AerodynamicGeometry

"""
    SurfaceAerodynamicCoefficients

Complete surface-integrated aerodynamic result. `force_body` stores
`(CX, CY, CZ)` and `moment_body` stores `(Cl, Cm, Cn)` in the Grant-Braun body
axes. `CA`, `CS`, and `CN` are the entry-vehicle body coefficients, while
`CD`, `CY_wind`, and `CL` are the corresponding wind-axis coefficients.
"""
struct SurfaceAerodynamicCoefficients
    force_body::SVector{3, Float64}
    moment_body::SVector{3, Float64}
    CA::Float64
    CS::Float64
    CN::Float64
    CD::Float64
    CY_wind::Float64
    CL::Float64
    exposed_area_fraction::Float64
end


const NewtonianAerodynamicCoefficients = SurfaceAerodynamicCoefficients

@inline function _surface_aerodynamic_coefficients(
    force_body::SVector{3, Float64},
    moment_body::SVector{3, Float64},
    freestream::SVector{3, Float64},
    active_area_fraction::Real,
)
    CA = -force_body[1]
    CS = force_body[2]
    CN = -force_body[3]
    alpha = atan(-freestream[3], -freestream[1])
    lift_unit_body = SVector{3, Float64}(sin(alpha), 0.0, -cos(alpha))
    vehicle_velocity_unit_body = -freestream
    side_unit_body = cross(vehicle_velocity_unit_body, lift_unit_body)
    return SurfaceAerodynamicCoefficients(
        force_body,
        moment_body,
        CA,
        CS,
        CN,
        dot(force_body, freestream),
        dot(force_body, side_unit_body),
        dot(force_body, lift_unit_body),
        Float64(active_area_fraction),
    )
end

@inline function _newtonian_cp_scale(
    pressure_model::Symbol,
    mach::Union{Nothing, Real},
    gamma::Union{Nothing, Real},
)
    if pressure_model === :regular_newtonian || pressure_model === :newtonian
        return 2.0
    elseif pressure_model === :modified_newtonian
        mach === nothing && throw(ArgumentError("modified Newtonian theory requires `mach`."))
        gamma === nothing && throw(ArgumentError("modified Newtonian theory requires `gamma`."))
        return modified_newtonian_cp_max(mach, gamma)
    end
    throw(ArgumentError(
        "Unsupported hypersonic pressure model $(repr(pressure_model)); " *
        "expected :regular_newtonian, :newtonian, or :modified_newtonian.",
    ))
end

function _gauss_legendre_interval(order::Int, lower::Float64, upper::Float64)
    order >= 1 || throw(ArgumentError("Gauss-Legendre quadrature order must be positive."))
    isfinite(lower) && isfinite(upper) && upper > lower || throw(ArgumentError(
        "quadrature bounds must be finite and strictly increasing.",
    ))
    if order == 1
        return [0.5 * (lower + upper)], [upper - lower]
    end

    indices = collect(1:(order - 1))
    off_diagonal = indices ./ sqrt.(4.0 .* indices.^2 .- 1.0)
    decomposition = eigen(SymTridiagonal(zeros(Float64, order), off_diagonal))
    midpoint = 0.5 * (lower + upper)
    half_width = 0.5 * (upper - lower)
    nodes = midpoint .+ half_width .* decomposition.values
    weights = 2.0 .* half_width .* decomposition.vectors[1, :].^2
    return nodes, weights
end

"""
    aerodynamic_surface_quadrature(
        position, inward_area_vector, u_bounds, v_bounds;
        quadrature_order_u=12, quadrature_order_v=48,
    )

Discretize the *integration*, not the shape, of a continuously parametrized
surface. `position(u, v)` implements Eq. (9), while
`inward_area_vector(u, v)` returns the inward `r_u x r_v` vector whose norm is
the differential-area Jacobian in Eqs. (10) and (12).
"""
function aerodynamic_surface_quadrature(
    position,
    inward_area_vector,
    u_bounds::Tuple{<:Real, <:Real},
    v_bounds::Tuple{<:Real, <:Real};
    quadrature_order_u::Int=12,
    quadrature_order_v::Int=48,
)
    u_nodes, u_weights = _gauss_legendre_interval(
        quadrature_order_u,
        Float64(u_bounds[1]),
        Float64(u_bounds[2]),
    )
    v_nodes, v_weights = _gauss_legendre_interval(
        quadrature_order_v,
        Float64(v_bounds[1]),
        Float64(v_bounds[2]),
    )
    sample_count = quadrature_order_u * quadrature_order_v
    positions = Vector{SVector{3, Float64}}(undef, sample_count)
    normals = Vector{SVector{3, Float64}}(undef, sample_count)
    area_weights = Vector{Float64}(undef, sample_count)

    sample_index = 1
    @inbounds for u_index in eachindex(u_nodes), v_index in eachindex(v_nodes)
        u = u_nodes[u_index]
        v = v_nodes[v_index]
        area_vector = SVector{3, Float64}(inward_area_vector(u, v))
        area_jacobian = norm(area_vector)
        area_jacobian > eps(Float64) || throw(DomainError(
            area_vector,
            "parametric surface has a singular quadrature point.",
        ))
        positions[sample_index] = SVector{3, Float64}(position(u, v))
        normals[sample_index] = area_vector / area_jacobian
        area_weights[sample_index] =
            u_weights[u_index] * v_weights[v_index] * area_jacobian
        sample_index += 1
    end
    return AerodynamicSurfaceQuadrature(positions, normals, area_weights)
end

"""
    newtonian_surface_quadrature(args...; kwargs...)

Backward-compatible name for [`aerodynamic_surface_quadrature`](@ref).
"""
newtonian_surface_quadrature(args...; kwargs...) =
    aerodynamic_surface_quadrature(args...; kwargs...)

"""
    aerodynamic_plate_surface(
        center_body_m, span_u_body_m, span_v_body_m;
        inward_normal_body=cross(span_u_body_m, span_v_body_m),
        moment_reference_body_m=(0, 0, 0),
    )

Construct one gas-facing rectangular plate surface. The two span vectors are
the plate's full edge vectors, so its area is `norm(cross(span_u, span_v))`.
Because both supported local pressure laws are uniform over a planar plate,
the centroid sample integrates its force and moment exactly.

For a zero-thickness plate exposed on both sides, create two surfaces with
opposite `inward_normal_body` values and pass both to
[`combine_aerodynamic_surfaces`](@ref).
"""
function aerodynamic_plate_surface(
    center_body_m::AbstractVector{<:Real},
    span_u_body_m::AbstractVector{<:Real},
    span_v_body_m::AbstractVector{<:Real};
    inward_normal_body::AbstractVector{<:Real}=cross(
        SVector{3, Float64}(span_u_body_m),
        SVector{3, Float64}(span_v_body_m),
    ),
    moment_reference_body_m::AbstractVector{<:Real}=SVector{3, Float64}(0.0, 0.0, 0.0),
)
    center = SVector{3, Float64}(center_body_m)
    span_u = SVector{3, Float64}(span_u_body_m)
    span_v = SVector{3, Float64}(span_v_body_m)
    inward_normal = SVector{3, Float64}(inward_normal_body)
    reference_point = SVector{3, Float64}(moment_reference_body_m)
    all(isfinite, center) || throw(DomainError(center_body_m, "plate center must be finite."))
    all(isfinite, span_u) || throw(DomainError(span_u_body_m, "plate span vectors must be finite."))
    all(isfinite, span_v) || throw(DomainError(span_v_body_m, "plate span vectors must be finite."))
    all(isfinite, inward_normal) || throw(DomainError(
        inward_normal_body,
        "plate inward normal must be finite.",
    ))
    all(isfinite, reference_point) || throw(DomainError(
        moment_reference_body_m,
        "moment reference coordinates must be finite.",
    ))
    area_vector = cross(span_u, span_v)
    area = norm(area_vector)
    area > eps(Float64) || throw(DomainError(
        (span_u_body_m, span_v_body_m),
        "plate span vectors must define a nonzero area.",
    ))
    normal_magnitude = norm(inward_normal)
    normal_magnitude > eps(Float64) || throw(DomainError(
        inward_normal_body,
        "plate inward normal must be nonzero.",
    ))
    unit_normal = inward_normal / normal_magnitude
    alignment = abs(dot(unit_normal, area_vector / area))
    alignment >= 1.0 - 1.0e-10 || throw(DomainError(
        inward_normal_body,
        "plate inward normal must be perpendicular to both span vectors.",
    ))
    return AerodynamicSurfaceQuadrature(
        [center - reference_point],
        [unit_normal],
        [area],
    )
end

"""
    combine_aerodynamic_surfaces(
        surfaces; reference_area_m2, reference_length_m, label="",
    )

Combine continuously integrated basic surfaces using common reference values.
This implements the normalization principle in Grant-Braun Eq. (23), but the
resulting geometry is not tied to a particular flow regime.
"""
function combine_aerodynamic_surfaces(
    surfaces::AbstractVector{AerodynamicSurfaceQuadrature};
    reference_area_m2::Real,
    reference_length_m::Real,
    label::AbstractString="",
)
    isempty(surfaces) && throw(ArgumentError("at least one aerodynamic surface is required."))
    positions = SVector{3, Float64}[]
    normals = SVector{3, Float64}[]
    weights = Float64[]
    for surface in surfaces
        append!(positions, surface.positions_body_m)
        append!(normals, surface.inward_normals_body)
        append!(weights, surface.area_weights_m2)
    end
    return AerodynamicGeometry(
        AerodynamicSurfaceQuadrature(positions, normals, weights),
        reference_area_m2,
        reference_length_m;
        label=label,
    )
end

"""
    combine_newtonian_surfaces(args...; kwargs...)

Backward-compatible name for [`combine_aerodynamic_surfaces`](@ref).
"""
combine_newtonian_surfaces(args...; kwargs...) =
    combine_aerodynamic_surfaces(args...; kwargs...)

"""
    sphere_cone_aerodynamic_geometry(
        nose_radius, base_radius, cone_half_angle_rad;
        moment_reference_body_m=(0, 0, 0),
        quadrature_order_meridional=12, quadrature_order_azimuth=48,
        include_base=true,
    )

Construct the exact spherical-segment and conical-frustum parametrizations in
Grant-Braun Eqs. (15), (16), and (23). Coordinates follow the paper: `+x`
points forward, `+y` right, and `+z` down. By default the rear base disk is
included to form the complete exterior surface required by free-molecular
aerodynamics. Set `include_base=false` for the Grant-Braun forebody-only
Newtonian model.
"""
function sphere_cone_aerodynamic_geometry(
    nose_radius::Real,
    base_radius::Real,
    cone_half_angle_rad::Real;
    moment_reference_body_m::AbstractVector{<:Real}=SVector{3, Float64}(0.0, 0.0, 0.0),
    quadrature_order_meridional::Int=12,
    quadrature_order_azimuth::Int=48,
    include_base::Bool=true,
)
    rn = Float64(nose_radius)
    rb = Float64(base_radius)
    delta = Float64(cone_half_angle_rad)
    isfinite(rn) && rn >= 0.0 || throw(DomainError(
        nose_radius,
        "nose radius must be finite and nonnegative.",
    ))
    isfinite(rb) && rb > 0.0 || throw(DomainError(
        base_radius,
        "base radius must be finite and positive.",
    ))
    isfinite(delta) && 0.0 < delta < π / 2.0 || throw(DomainError(
        cone_half_angle_rad,
        "cone half-angle must lie strictly between 0 and pi/2 radians.",
    ))
    quadrature_order_meridional >= 2 || throw(ArgumentError(
        "sphere-cone meridional quadrature order must be at least 2.",
    ))
    quadrature_order_azimuth >= 8 || throw(ArgumentError(
        "sphere-cone azimuth quadrature order must be at least 8.",
    ))

    reference_point = SVector{3, Float64}(moment_reference_body_m)
    all(isfinite, reference_point) || throw(DomainError(
        moment_reference_body_m,
        "moment reference coordinates must be finite.",
    ))
    sin_delta = sin(delta)
    cos_delta = cos(delta)
    junction_radius = rn * cos_delta
    junction_radius <= rb || throw(DomainError(
        rn / rb,
        "nose/base-radius ratio is incompatible with the tangent sphere-cone geometry.",
    ))
    junction_x = rn * sin_delta
    virtual_apex_x = rn / sin_delta
    base_x = virtual_apex_x - rb / tan(delta)
    nose_x = rn
    reference_length = nose_x - base_x
    reference_length > 0.0 || throw(DomainError(
        reference_length,
        "sphere-cone reference length must be positive.",
    ))

    surfaces = AerodynamicSurfaceQuadrature[]
    if rn > 0.0 && nose_x - junction_x > eps(Float64)
        sphere_position = function (u, v)
            radial = sqrt(max(0.0, rn^2 - u^2))
            return SVector{3, Float64}(u, radial * cos(v), -radial * sin(v)) - reference_point
        end
        sphere_inward_area = function (u, v)
            radial = sqrt(max(0.0, rn^2 - u^2))
            # Eq. (16) produces an outward r_u x r_v for the paper's axis
            # orientation, so use its negative to satisfy Eq. (12)'s inward normal.
            return -SVector{3, Float64}(u, radial * cos(v), -radial * sin(v))
        end
        push!(surfaces, aerodynamic_surface_quadrature(
            sphere_position,
            sphere_inward_area,
            (junction_x, nose_x),
            (0.0, 2.0 * π);
            quadrature_order_u=quadrature_order_meridional,
            quadrature_order_v=quadrature_order_azimuth,
        ))
    end

    if rb - junction_radius > eps(Float64)
        cone_position = function (radius, v)
            x = virtual_apex_x - radius / tan(delta)
            return SVector{3, Float64}(x, radius * cos(v), -radius * sin(v)) - reference_point
        end
        cone_inward_area = function (radius, v)
            return SVector{3, Float64}(
                -radius,
                -radius * cos(v) / tan(delta),
                radius * sin(v) / tan(delta),
            )
        end
        push!(surfaces, aerodynamic_surface_quadrature(
            cone_position,
            cone_inward_area,
            (junction_radius, rb),
            (0.0, 2.0 * π);
            quadrature_order_u=quadrature_order_meridional,
            quadrature_order_v=quadrature_order_azimuth,
        ))
    end

    if include_base
        base_position = function (radius, v)
            return SVector{3, Float64}(
                base_x,
                radius * cos(v),
                -radius * sin(v),
            ) - reference_point
        end
        base_inward_area = function (radius, _)
            # The exterior base normal points rearward (-x), so the inward
            # normal used by the pressure-force convention points forward.
            return SVector{3, Float64}(radius, 0.0, 0.0)
        end
        push!(surfaces, aerodynamic_surface_quadrature(
            base_position,
            base_inward_area,
            (0.0, rb),
            (0.0, 2.0 * π);
            quadrature_order_u=quadrature_order_meridional,
            quadrature_order_v=quadrature_order_azimuth,
        ))
    end

    return combine_aerodynamic_surfaces(
        surfaces;
        reference_area_m2=π * rb^2,
        reference_length_m=reference_length,
        label=include_base ? "closed sphere-cone" : "Grant-Braun sphere-cone forebody",
    )
end

"""
    sphere_cone_newtonian_geometry(args...; kwargs...)

Construct the Grant-Braun forebody surface. This compatibility constructor
deliberately excludes the rear base, preserving the published Newtonian model.
Use [`sphere_cone_aerodynamic_geometry`](@ref) for a closed geometry shared by
continuum and free-molecular laws.
"""
function sphere_cone_newtonian_geometry(args...; kwargs...)
    haskey(kwargs, :include_base) && throw(ArgumentError(
        "sphere_cone_newtonian_geometry fixes include_base=false; use sphere_cone_aerodynamic_geometry to select it.",
    ))
    return sphere_cone_aerodynamic_geometry(args...; kwargs..., include_base=false)
end

"""
    newtonian_aerodynamic_coefficients(
        geometry, freestream_unit_body;
        pressure_model=:regular_newtonian, mach=nothing, gamma=nothing,
    )

Evaluate Grant-Braun Eqs. (1)-(5) over the complete continuously parametrized
surface. At every quadrature point, `Cp = scale*(Vinf dot n)^2` when
`Vinf dot n > 0`; otherwise `Cp = 0` exactly as required by the paper's shadow
rule. Arbitrary angle of attack, sideslip, partial shadowing, force, and moment
coefficients are supported.
"""
function newtonian_aerodynamic_coefficients(
    geometry::NewtonianAerodynamicGeometry,
    freestream_unit_body::AbstractVector{<:Real};
    pressure_model::Symbol=:regular_newtonian,
    mach::Union{Nothing, Real}=nothing,
    gamma::Union{Nothing, Real}=nothing,
)
    freestream = SVector{3, Float64}(freestream_unit_body)
    all(isfinite, freestream) || throw(DomainError(
        freestream_unit_body,
        "freestream direction must be finite.",
    ))
    freestream_magnitude = norm(freestream)
    freestream_magnitude > eps(Float64) || throw(DomainError(
        freestream_unit_body,
        "freestream direction must be nonzero.",
    ))
    freestream /= freestream_magnitude
    cp_scale = _newtonian_cp_scale(pressure_model, mach, gamma)

    force_integral = MVector{3, Float64}(0.0, 0.0, 0.0)
    moment_integral = MVector{3, Float64}(0.0, 0.0, 0.0)
    exposed_area = 0.0
    total_area = 0.0
    surface = geometry.surface
    @inbounds for index in eachindex(surface.area_weights_m2)
        normal = surface.inward_normals_body[index]
        area_weight = surface.area_weights_m2[index]
        incidence = dot(freestream, normal)
        total_area += area_weight
        if incidence > 0.0
            cp_area = cp_scale * incidence^2 * area_weight
            force_integral .+= cp_area * normal
            moment_integral .+= cp_area * cross(surface.positions_body_m[index], normal)
            exposed_area += area_weight
        end
    end

    force_body = SVector{3, Float64}(force_integral) / geometry.reference_area_m2
    moment_body = SVector{3, Float64}(moment_integral) /
                  (geometry.reference_area_m2 * geometry.reference_length_m)
    return _surface_aerodynamic_coefficients(
        force_body,
        moment_body,
        freestream,
        exposed_area / total_area,
    )
end

function newtonian_aerodynamic_coefficients(
    geometry::NewtonianAerodynamicGeometry,
    alpha_rad::Real,
    beta_rad::Real;
    kwargs...,
)
    alpha = Float64(alpha_rad)
    beta = Float64(beta_rad)
    isfinite(alpha) || throw(DomainError(alpha_rad, "angle of attack must be finite."))
    isfinite(beta) || throw(DomainError(beta_rad, "sideslip must be finite."))
    freestream = SVector{3, Float64}(
        -cos(alpha) * cos(beta),
        -sin(beta),
        -sin(alpha) * cos(beta),
    )
    return newtonian_aerodynamic_coefficients(geometry, freestream; kwargs...)
end

"""
    newtonian_stability_derivatives(
        geometry, alpha_rad, beta_rad;
        pressure_model=:regular_newtonian, mach=nothing, gamma=nothing,
    )

Evaluate the pitch and yaw stability derivatives defined by Grant-Braun Eqs.
(6) and (7) by differentiating the complete surface integrand. No coefficient
finite differencing is used. The moving shadow boundary contributes no boundary
term because `Cp` is zero at that boundary.
"""
function newtonian_stability_derivatives(
    geometry::NewtonianAerodynamicGeometry,
    alpha_rad::Real,
    beta_rad::Real;
    pressure_model::Symbol=:regular_newtonian,
    mach::Union{Nothing, Real}=nothing,
    gamma::Union{Nothing, Real}=nothing,
)
    alpha = Float64(alpha_rad)
    beta = Float64(beta_rad)
    isfinite(alpha) || throw(DomainError(alpha_rad, "angle of attack must be finite."))
    isfinite(beta) || throw(DomainError(beta_rad, "sideslip must be finite."))
    cp_scale = _newtonian_cp_scale(pressure_model, mach, gamma)
    sin_alpha, cos_alpha = sincos(alpha)
    sin_beta, cos_beta = sincos(beta)
    freestream = SVector{3, Float64}(
        -cos_alpha * cos_beta,
        -sin_beta,
        -sin_alpha * cos_beta,
    )
    freestream_alpha = SVector{3, Float64}(
        sin_alpha * cos_beta,
        0.0,
        -cos_alpha * cos_beta,
    )
    freestream_beta = SVector{3, Float64}(
        cos_alpha * sin_beta,
        -cos_beta,
        sin_alpha * sin_beta,
    )

    moment_alpha_integral = MVector{3, Float64}(0.0, 0.0, 0.0)
    moment_beta_integral = MVector{3, Float64}(0.0, 0.0, 0.0)
    surface = geometry.surface
    @inbounds for index in eachindex(surface.area_weights_m2)
        normal = surface.inward_normals_body[index]
        incidence = dot(freestream, normal)
        if incidence > 0.0
            moment_arm = cross(surface.positions_body_m[index], normal)
            common = 2.0 * cp_scale * incidence * surface.area_weights_m2[index]
            moment_alpha_integral .+= common * dot(freestream_alpha, normal) * moment_arm
            moment_beta_integral .+= common * dot(freestream_beta, normal) * moment_arm
        end
    end
    normalization = geometry.reference_area_m2 * geometry.reference_length_m
    return (
        Cm_alpha=moment_alpha_integral[2] / normalization,
        Cn_beta=moment_beta_integral[3] / normalization,
    )
end

@inline function _validate_sphere_cone_newtonian_inputs(
    alpha_rad::Float64,
    nose_radius::Float64,
    base_radius::Float64,
    cone_half_angle_rad::Float64,
)
    isfinite(alpha_rad) || throw(DomainError(alpha_rad, "angle of attack must be finite."))
    isfinite(nose_radius) && nose_radius >= 0.0 || throw(
        DomainError(nose_radius, "nose radius must be finite and nonnegative."),
    )
    isfinite(base_radius) && base_radius > 0.0 || throw(
        DomainError(base_radius, "base radius must be finite and positive."),
    )
    isfinite(cone_half_angle_rad) && 0.0 < cone_half_angle_rad < π / 2.0 || throw(
        DomainError(
            cone_half_angle_rad,
            "cone half-angle must lie strictly between 0 and pi/2 radians.",
        ),
    )
    abs(alpha_rad) <= cone_half_angle_rad || throw(
        DomainError(
            alpha_rad,
            "the unshadowed sphere-cone relation requires abs(alpha) <= cone half-angle.",
        ),
    )

    radius_ratio = nose_radius / base_radius
    frustum_factor = 1.0 - (radius_ratio * cos(cone_half_angle_rad))^2
    frustum_factor >= 0.0 || throw(
        DomainError(
            radius_ratio,
            "nose/base-radius ratio is incompatible with the tangent sphere-cone geometry.",
        ),
    )
    return radius_ratio, frustum_factor
end

"""
    _legacy_unshadowed_sphere_cone_cn_ca(
        alpha_rad, nose_radius, base_radius, cone_half_angle_rad,
    ) -> (CN, CA)

Preserve the preliminary study's compact zero-sideslip relation only as an
internal regression comparison. It is not the complete Grant-Braun spherical
segment solution because it omits part of the cap's angle-dependent force.
"""
@inline function _legacy_unshadowed_sphere_cone_cn_ca(
    alpha_rad::Real,
    nose_radius::Real,
    base_radius::Real,
    cone_half_angle_rad::Real,
)
    alpha_float = Float64(alpha_rad)
    nose_radius_float = Float64(nose_radius)
    base_radius_float = Float64(base_radius)
    half_angle_float = Float64(cone_half_angle_rad)
    radius_ratio, frustum_factor = _validate_sphere_cone_newtonian_inputs(
        alpha_float,
        nose_radius_float,
        base_radius_float,
        half_angle_float,
    )

    sin_delta = sin(half_angle_float)
    cos_delta = cos(half_angle_float)
    ratio_sq = radius_ratio^2
    CN = frustum_factor * cos_delta^2 * sin(2.0 * alpha_float)
    CA = (1.0 - sin_delta^4) * ratio_sq +
         (
             2.0 * sin_delta^2 * cos(alpha_float)^2 +
             cos_delta^2 * sin(alpha_float)^2
         ) * frustum_factor
    return CN, CA
end

"""
    regular_newtonian_sphere_cone_cn_ca(
        alpha_rad, nose_radius, base_radius, cone_half_angle_rad,
    ) -> (CN, CA)

Evaluate the complete regular-Newtonian sphere-cone surface model at zero
sideslip, including spherical-cap force and surface shadowing.
"""
function regular_newtonian_sphere_cone_cn_ca(
    alpha_rad::Real,
    nose_radius::Real,
    base_radius::Real,
    cone_half_angle_rad::Real,
)
    geometry = sphere_cone_newtonian_geometry(
        nose_radius,
        base_radius,
        cone_half_angle_rad,
    )
    coefficients = newtonian_aerodynamic_coefficients(
        geometry,
        alpha_rad,
        0.0;
        pressure_model=:regular_newtonian,
    )
    return coefficients.CN, coefficients.CA
end

"""
    modified_newtonian_sphere_cone_cn_ca(
        alpha_rad, nose_radius, base_radius, cone_half_angle_rad, mach, gamma,
    ) -> (CN, CA)

Evaluate the complete modified-Newtonian sphere-cone surface model. Grant and
Braun identify the stagnation-pressure multiplier as the conversion from their
regular Newtonian analytic results to modified Newtonian theory.
"""
@inline function modified_newtonian_sphere_cone_cn_ca(
    alpha_rad::Real,
    nose_radius::Real,
    base_radius::Real,
    cone_half_angle_rad::Real,
    mach::Real,
    gamma::Real,
)
    geometry = sphere_cone_newtonian_geometry(
        nose_radius,
        base_radius,
        cone_half_angle_rad,
    )
    coefficients = newtonian_aerodynamic_coefficients(
        geometry,
        alpha_rad,
        0.0;
        pressure_model=:modified_newtonian,
        mach=mach,
        gamma=gamma,
    )
    return coefficients.CN, coefficients.CA
end

"""
    sphere_cone_cn_ca(
        alpha_rad, nose_radius, base_radius, cone_half_angle_rad;
        pressure_model=:regular_newtonian, mach=nothing, gamma=nothing,
    ) -> (CN, CA)

Select the regular or modified Newtonian sphere-cone model. Accepted regular
selectors are `:regular_newtonian` and `:newtonian`. The
`:modified_newtonian` selector requires both `mach` and `gamma`.
"""
function sphere_cone_cn_ca(
    alpha_rad::Real,
    nose_radius::Real,
    base_radius::Real,
    cone_half_angle_rad::Real;
    pressure_model::Symbol=:regular_newtonian,
    mach::Union{Nothing, Real}=nothing,
    gamma::Union{Nothing, Real}=nothing,
)
    if pressure_model === :regular_newtonian || pressure_model === :newtonian
        return regular_newtonian_sphere_cone_cn_ca(
            alpha_rad,
            nose_radius,
            base_radius,
            cone_half_angle_rad,
        )
    elseif pressure_model === :modified_newtonian
        mach === nothing && throw(ArgumentError("modified Newtonian theory requires `mach`."))
        gamma === nothing && throw(ArgumentError("modified Newtonian theory requires `gamma`."))
        return modified_newtonian_sphere_cone_cn_ca(
            alpha_rad,
            nose_radius,
            base_radius,
            cone_half_angle_rad,
            mach,
            gamma,
        )
    end
    throw(ArgumentError(
        "Unsupported hypersonic pressure model $(repr(pressure_model)); " *
        "expected :regular_newtonian, :newtonian, or :modified_newtonian.",
    ))
end

"""
    cn_ca_to_cl_cd(alpha_rad, CN, CA) -> (CL, CD)

Transform zero-sideslip body-axis normal and axial coefficients into wind-axis
lift and drag coefficients using the Grant-Braun axis convention.
"""
@inline function cn_ca_to_cl_cd(alpha_rad::Real, CN::Real, CA::Real)
    alpha_float = Float64(alpha_rad)
    isfinite(alpha_float) || throw(DomainError(alpha_rad, "angle of attack must be finite."))
    CN_float = Float64(CN)
    CA_float = Float64(CA)
    isfinite(CN_float) || throw(DomainError(CN, "CN must be finite."))
    isfinite(CA_float) || throw(DomainError(CA, "CA must be finite."))
    CL = CN_float * cos(alpha_float) - CA_float * sin(alpha_float)
    CD = CA_float * cos(alpha_float) + CN_float * sin(alpha_float)
    return CL, CD
end

"""
    sphere_cone_cl_cd(args...; kwargs...) -> (CL, CD)

Evaluate sphere-cone `CN` and `CA` with [`sphere_cone_cn_ca`](@ref), then
transform them to `CL` and `CD`.
"""
function sphere_cone_cl_cd(
    alpha_rad::Real,
    nose_radius::Real,
    base_radius::Real,
    cone_half_angle_rad::Real;
    kwargs...,
)
    CN, CA = sphere_cone_cn_ca(
        alpha_rad,
        nose_radius,
        base_radius,
        cone_half_angle_rad;
        kwargs...,
    )
    return cn_ca_to_cl_cd(alpha_rad, CN, CA)
end
