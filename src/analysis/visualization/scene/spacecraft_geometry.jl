# Spacecraft model -> renderable geometry.
#
# Non-root link poses are not integrated state on the standard dynamics path:
# they live in the mutable `Link` objects and are changed by control code
# through `rotate_link`. The `link_pose` save field therefore snapshots
# `Link.r`/`Link.q` at every saved step, and this file fixes the link order
# that both the geometry and those columns share.

@inline _is_root_link(model::SpacecraftModel, link::Link)::Bool = link === model.root

"""
    link_pose_link_indices(model) -> Vector{Int}

Indices into `model.links` of every non-root link, in the order the
`link_pose` save field and `SpacecraftGeometry.links[2:end]` use.
"""
function link_pose_link_indices(model::SpacecraftModel)::Vector{Int}
    return [i for (i, link) in enumerate(model.links) if !_is_root_link(model, link)]
end

"""
    link_pose_vector(model) -> Vector{Float64}

Concatenated `[rx, ry, rz, qx, qy, qz, qw]` for every non-root link, in
`link_pose_link_indices` order; empty when the spacecraft is a single body.
Written to columns `sc{i}_link_pose_{1..7n}`.
"""
function link_pose_vector(model::SpacecraftModel)::Vector{Float64}
    indices = link_pose_link_indices(model)
    out = Vector{Float64}(undef, LINK_POSE_STRIDE * length(indices))
    k = 0
    @inbounds for i in indices
        link = model.links[i]
        out[k + 1] = link.r[1]
        out[k + 2] = link.r[2]
        out[k + 3] = link.r[3]
        out[k + 4] = link.q[1]
        out[k + 5] = link.q[2]
        out[k + 6] = link.q[3]
        out[k + 7] = link.q[4]
        k += LINK_POSE_STRIDE
    end
    return out
end

@inline function _link_box(name::String, link::Link, is_root::Bool)::LinkBox
    r = is_root ? SVector{3, Float64}(0.0, 0.0, 0.0) : _svec3(link.r)
    q = is_root ? SVector{4, Float64}(0.0, 0.0, 0.0, 1.0) : _svec4(link.q)
    return LinkBox(name, is_root, _svec3(link.dims), r, q, Float64(link.m))
end

# Half the box diagonal bounds every corner regardless of the link's rotation,
# so |r| plus that bound is a rotation-invariant radius for the whole assembly.
@inline _link_extent_m(box::LinkBox)::Float64 = norm(box.r_m) + 0.5 * norm(box.dims_m)

function _geometry_index_map(model::SpacecraftModel)::Dict{Int, Int}
    # model.links index -> SpacecraftGeometry.links index (root is 1)
    mapping = Dict{Int, Int}()
    root_idx = findfirst(link -> _is_root_link(model, link), model.links)
    root_idx === nothing || (mapping[root_idx] = 1)
    for (k, i) in enumerate(link_pose_link_indices(model))
        mapping[i] = k + 1
    end
    return mapping
end

"""
    arm_geometry(plan::RobotArmPlan) -> ArmGeometry

Renderable arm from a planned motion's `ClothArmModel`.
"""
function arm_geometry(plan::RobotArmPlan)::ArmGeometry
    model = plan.model
    links = ArmLinkGeometry[
        ArmLinkGeometry(String(link.name), _svec3(link.vector_parent), _svec3(link.com_offset_parent), Float64(link.radius_m), Float64(link.mass_kg))
        for link in model.links
    ]
    return ArmGeometry(links, _svec3(model.mount_offset_body), Float64(cloth_total_reach(model)))
end

"""
    robot_arm_plan_for(args, sat_idx) -> Union{Nothing, RobotArmPlan}

The `RobotArmPlan` carried by a control effector for spacecraft `sat_idx`
(the same lookup the engine uses to couple the arm state), or `nothing`.
"""
function robot_arm_plan_for(args::SimulationConfiguration, sat_idx::Integer)
    hasproperty(args, :control_model) && hasproperty(args.control_model, :control_effectors) || return nothing
    for effector in args.control_model.control_effectors
        hasproperty(effector, :plan) && hasproperty(effector, :spacecraft_idx) || continue
        getproperty(effector, :spacecraft_idx) == sat_idx || continue
        plan = getproperty(effector, :plan)
        plan isa RobotArmPlan && return plan
    end
    return nothing
end

"""
    arm_pose_vector(sc_view, r_ii_m) -> Vector{Float64}

`[rx, ry, rz, qx, qy, qz, qw]` per arm link from a spacecraft state view
carrying `arm_r`/`arm_q` (inertial COM position in metres, inertial
quaternion), with positions relative to the spacecraft position `r_ii_m`.
Empty when the view has no arm state.
"""
function arm_pose_vector(sc_view, r_ii_m)::Vector{Float64}
    hasproperty(sc_view, :arm_r) && hasproperty(sc_view, :arm_q) || return Float64[]
    arm_r = getproperty(sc_view, :arm_r)
    arm_q = getproperty(sc_view, :arm_q)
    n = size(arm_r, 2)
    out = Vector{Float64}(undef, LINK_POSE_STRIDE * n)
    @inbounds for i in 1:n
        k = LINK_POSE_STRIDE * (i - 1)
        out[k + 1] = Float64(arm_r[1, i]) - Float64(r_ii_m[1])
        out[k + 2] = Float64(arm_r[2, i]) - Float64(r_ii_m[2])
        out[k + 3] = Float64(arm_r[3, i]) - Float64(r_ii_m[3])
        out[k + 4] = Float64(arm_q[1, i])
        out[k + 5] = Float64(arm_q[2, i])
        out[k + 6] = Float64(arm_q[3, i])
        out[k + 7] = Float64(arm_q[4, i])
    end
    return out
end

"""
    spacecraft_geometry(model; name="sc<id>", stl_path=nothing, arm=nothing) -> SpacecraftGeometry

Convert a `SpacecraftModel` into boxes and glyphs. Link order is root first,
then `link_pose_link_indices(model)`. Thruster, facet and joint glyphs refer to
links by that geometry index. `stl_path` records a CAD override the viewer may
use instead of the boxes; it is stored verbatim.
"""
function spacecraft_geometry(
    model::SpacecraftModel;
    name::AbstractString="sc$(model.id)",
    stl_path::Union{Nothing, AbstractString}=nothing,
    arm::Union{Nothing, ArmGeometry, RobotArmPlan}=nothing
)::SpacecraftGeometry
    arm_geom = arm isa RobotArmPlan ? arm_geometry(arm) : arm
    boxes = LinkBox[_link_box("root", model.root, true)]
    for (k, i) in enumerate(link_pose_link_indices(model))
        push!(boxes, _link_box("link$(k)", model.links[i], false))
    end

    index_map = _geometry_index_map(model)
    thrusters = ThrusterGlyph[]
    facets = FacetGlyph[]
    for (i, link) in enumerate(model.links)
        geometry_idx = get(index_map, i, 0)
        geometry_idx == 0 && continue
        for thruster in link.thrusters
            push!(thrusters, ThrusterGlyph(geometry_idx, _svec3(thruster.location), _svec3(thruster.direction), Float64(thruster.max_thrust)))
        end
        for facet in link.SRP_facets
            push!(facets, FacetGlyph(geometry_idx, String(facet.name), _svec3(facet.normal_vector), _svec3(facet.cp), Float64(facet.area)))
        end
    end

    joints = JointGlyph[]
    for joint in model.joints
        i1 = findfirst(link -> link === joint.link1, model.links)
        i2 = findfirst(link -> link === joint.link2, model.links)
        g1 = i1 === nothing ? 0 : get(index_map, i1, 0)
        g2 = i2 === nothing ? 0 : get(index_map, i2, 0)
        push!(joints, JointGlyph(g1, g2, _svec3(joint.p1ᵇ), _svec3(joint.p2ᵇ)))
    end

    bounding_radius = maximum(_link_extent_m, boxes)
    if arm_geom !== nothing
        bounding_radius = max(bounding_radius, norm(arm_geom.mount_offset_body_m) + arm_geom.reach_m)
    end
    return SpacecraftGeometry(
        Int(model.id),
        String(name),
        boxes,
        thrusters,
        facets,
        joints,
        bounding_radius,
        stl_path === nothing ? nothing : String(stl_path),
        arm_geom
    )
end

"""
    velocity_aligned_quaternion(r_i, v_i) -> SVector{4, Float64}

Attitude for runs without orientation simulation: body +x along the inertial
velocity and body +z along the nadir direction with its velocity component
removed. Returned in the state convention (scalar-last, `rot(q)` maps
inertial to body), so it can stand in for the saved `sc{i}_q` columns.
"""
function velocity_aligned_quaternion(r_i::AbstractVector{<:Real}, v_i::AbstractVector{<:Real})::SVector{4, Float64}
    r = _svec3(r_i)
    v = _svec3(v_i)
    speed = norm(v)
    speed > 0.0 || throw(ArgumentError("velocity_aligned_quaternion requires a non-zero velocity."))
    x_b = v / speed
    nadir = -r
    z_b = nadir - dot(nadir, x_b) * x_b
    if norm(z_b) <= eps(Float64) * max(1.0, norm(r))
        # Radial flight: any axis orthogonal to x_b will do; pick a stable one.
        helper = abs(x_b[3]) < 0.9 ? SVector{3, Float64}(0.0, 0.0, 1.0) : SVector{3, Float64}(1.0, 0.0, 0.0)
        z_b = helper - dot(helper, x_b) * x_b
    end
    z_b = z_b / norm(z_b)
    y_b = cross(z_b, x_b)
    # Rows of the inertial->body DCM are the body axes expressed in inertial.
    dcm = SMatrix{3, 3, Float64}(
        x_b[1], y_b[1], z_b[1],
        x_b[2], y_b[2], z_b[2],
        x_b[3], y_b[3], z_b[3]
    )
    return _svec4(dcm_to_quaternion(dcm))
end
