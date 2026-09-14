module SceneVisualization

# Post-hoc visualization scene layer. Builds a renderer-independent description
# of a run (planet, spacecraft geometry, epoch, frame rotation samples) from the
# `SimulationConfiguration`, and writes it as a JSON sidecar next to the
# results bundle. Nothing here touches the integrator; the trajectory itself
# stays in the Arrow results file and is joined to this sidecar by the viewer.
# Design record: docs/architecture/interactive_visualization_plan.md.

using StaticArrays
using LinearAlgebra
using Dates
using JSON
using TOML
using Base64
using Arrow
using DataFrames

using ..AbstractTypes: AbstractPlanet, AbstractEphemeridesModel
using ..SpacecraftModels: SpacecraftModel, Link
using ..RobotArmPlanning: RobotArmPlan
using ..Robotics: ClothArmModel, cloth_total_reach
using ..EphemeridesModels: planet_frame_lpi, ephemerides_time_seconds
using ..EphemeridesModels: _initial_time_datetime
using ..SimConfig: SimulationConfiguration, SimulationSettings
using ..IOConfig
using ..IOSerialization
import ..SimulationModel: rot, dcm_to_quaternion

export VisualizationScene, PlanetSpec, SpacecraftGeometry, LinkBox, AtmosphereSpec, atmosphere_spec
export ThrusterGlyph, FacetGlyph, JointGlyph, ArmLinkGeometry, ArmGeometry
export arm_geometry, robot_arm_plan_for, arm_pose_vector
export spacecraft_geometry, link_pose_link_indices, link_pose_vector
export planet_spec, planet_rotation_table, rotation_sample_times
export build_visualization_scene, visualization_scene_path, with_visualization_scene
export write_visualization_scene, read_visualization_scene, write_visualization_scene!
export velocity_aligned_quaternion, visualization_frame_budget
export export_visualization, write_viewer_dev_payload
export texture_manifest, texture_entry, texture_payload, build_viewer_frames, viewer_payload, model_payloads
export render_viewer_html, viewer_import_map, kept_row_indices
export EnsembleSample, sample_results_directory, with_results_directory, default_sample_scalar
export write_ensemble_manifest, read_ensemble_manifest, discover_ensemble_samples
export build_ensemble_frames, ensemble_time_axis, export_ensemble_visualization

const SCENE_SCHEMA_VERSION = 1
# Save field that snapshots every non-root link pose (see `link_pose_vector`).
const LINK_POSE_FIELD = :link_pose
const LINK_POSE_STRIDE = 7
const LINK_POSE_LAYOUT = ("rx", "ry", "rz", "qx", "qy", "qz", "qw")
# Save field for the integrated cloth robot-arm chain: per arm link the COM
# position relative to the spacecraft (inertial, metres) and its inertial
# quaternion, same 7-float layout.
const ARM_POSE_FIELD = :arm_pose

"""
    LinkBox

One rigid link rendered as an axis-aligned box in its own frame. `r_m` and `q`
place the link relative to the root bus (root frame, scalar-last quaternion);
both are zero/identity for the root itself.
"""
struct LinkBox
    name::String
    root::Bool
    dims_m::SVector{3, Float64}
    r_m::SVector{3, Float64}
    q::SVector{4, Float64}
    mass_kg::Float64
end

struct ThrusterGlyph
    link::Int
    location_m::SVector{3, Float64}
    direction::SVector{3, Float64}
    max_thrust_n::Float64
end

struct FacetGlyph
    link::Int
    name::String
    normal::SVector{3, Float64}
    cp_m::SVector{3, Float64}
    area_m2::Float64
end

struct JointGlyph
    link1::Int
    link2::Int
    p1_m::SVector{3, Float64}
    p2_m::SVector{3, Float64}
end

"""
    ArmLinkGeometry

One link of a cloth robot arm: the link vector and centre-of-mass offset in
the link's own frame (metres), its radius and mass. Drawn as a cylinder from
the joint along `vector_m`.
"""
struct ArmLinkGeometry
    name::String
    vector_m::SVector{3, Float64}
    com_offset_m::SVector{3, Float64}
    radius_m::Float64
    mass_kg::Float64
end

"""
    ArmGeometry

A robot arm mounted on a spacecraft: its links, the mount offset in the body
frame and the total reach. Poses come from the `arm_pose` save field.
"""
struct ArmGeometry
    links::Vector{ArmLinkGeometry}
    mount_offset_body_m::SVector{3, Float64}
    reach_m::Float64
end

"""
    SpacecraftGeometry

Renderable description of one spacecraft: its links as boxes (root first, then
the non-root links in `link_pose_link_indices` order), thruster and facet
glyphs, joint attachment points, a bounding radius for level-of-detail
switching, an optional STL override path, and an optional robot arm.
"""
struct SpacecraftGeometry
    id::Int
    name::String
    links::Vector{LinkBox}
    thrusters::Vector{ThrusterGlyph}
    facets::Vector{FacetGlyph}
    joints::Vector{JointGlyph}
    bounding_radius_m::Float64
    stl_path::Union{Nothing, String}
    arm::Union{Nothing, ArmGeometry}
end

"""
    PlanetSpec

Central body for the viewer: shape, spin, texture key, and a sampled table of
J2000-to-body-fixed quaternions so the viewer never needs SPICE.
"""
struct PlanetSpec
    name::String
    equatorial_radius_m::Float64
    polar_radius_m::Float64
    spin_rad_s::SVector{3, Float64}
    inertial_frame::String
    texture::String
    rotation_t_s::Vector{Float64}
    rotation_q_pi::Vector{SVector{4, Float64}}
end

include(joinpath(@__DIR__, "atmosphere_spec.jl"))

"""
    VisualizationScene

Everything the viewer needs besides the trajectory rows: epoch, planet,
spacecraft geometry, the atmosphere (when the run had one), and where the
trajectory and link-pose columns live.
"""
struct VisualizationScene
    schema::Int
    epoch_et_start_s::Float64
    epoch_utc::String
    planet::PlanetSpec
    spacecraft::Vector{SpacecraftGeometry}
    orientation_sim::Bool
    results_feather::String
    link_pose_field::String
    link_pose_stride::Int
    atmosphere::Union{Nothing, AtmosphereSpec}
end

for T in (LinkBox, ThrusterGlyph, FacetGlyph, JointGlyph, ArmLinkGeometry, ArmGeometry, SpacecraftGeometry, PlanetSpec, VisualizationScene)
    @eval Base.:(==)(a::$T, b::$T) = all(getfield(a, f) == getfield(b, f) for f in fieldnames($T))
end

@inline _svec3(v) = SVector{3, Float64}(Float64(v[1]), Float64(v[2]), Float64(v[3]))
@inline _svec4(v) = SVector{4, Float64}(Float64(v[1]), Float64(v[2]), Float64(v[3]), Float64(v[4]))

include(joinpath(@__DIR__, "spacecraft_geometry.jl"))
include(joinpath(@__DIR__, "planet_spec.jl"))
include(joinpath(@__DIR__, "scene_export.jl"))
include(joinpath(@__DIR__, "viewer_bundle.jl"))
include(joinpath(@__DIR__, "ensemble_export.jl"))

end # module SceneVisualization
