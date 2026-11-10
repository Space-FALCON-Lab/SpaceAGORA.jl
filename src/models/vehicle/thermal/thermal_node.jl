"""
    ThermalNodes

A spacecraft's lumped thermal node, mounted on its root link with
[`add_thermal!`](@ref). The record is `ThermalNode` (`Components`).
"""
module ThermalNodes

using ..Components: ThermalNode
using ..SpacecraftModels: Link, SpacecraftModel

export add_thermal!, has_thermal_node, thermal_node

"""
    add_thermal!(model::SpacecraftModel, link::Link, node::ThermalNode) -> node

Mount the spacecraft's thermal node on its root link (at most one). Returns the node.
"""
function add_thermal!(model::SpacecraftModel, link::Link, node::ThermalNode)
    any(l -> l === link, model.links) || throw(ArgumentError(
        "The thermal node mounts on one of the spacecraft's own links; this link is not in its link list."))
    link === model.root || throw(ArgumentError("The thermal node mounts on the spacecraft's root link."))
    isempty(link.thermal_nodes) || throw(ArgumentError("The spacecraft already has a thermal node; mount at most one."))
    push!(link.thermal_nodes, node)
    return node
end

"""Whether `model` carries a thermal node."""
has_thermal_node(model::SpacecraftModel)::Bool = !isempty(model.root.thermal_nodes)

"""The spacecraft's thermal node, or `nothing`."""
thermal_node(model::SpacecraftModel) = has_thermal_node(model) ? model.root.thermal_nodes[1] : nothing

end
