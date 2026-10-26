module LaserForceEffectors

# The forces between linked spacecraft come from each link's `link_force`
# (src/models/vehicle/interlinks/links.jl). These are the names they had when
# every link was a laser link: deprecated, they forward to that contract.

using StaticArrays
using ...InterLinkModels: InterlinkModel, AbstractInterlink, LinkContext, link_force, link_endpoint,
    interlink_force_on_spacecraft, _interlink_deprecation

export force_on_endpoint, laser_force_on_spacecraft

"""
    force_on_endpoint(current, other, link)

Deprecated, use `link_force`: the force `link` (a `LaserMomentumLink` or another
kind) puts on the spacecraft state `current` from `other`. The old name will be
removed in SpaceAGORA 1.0.
"""
function force_on_endpoint(current, other, link::AbstractInterlink)
    _interlink_deprecation((InterlinkModel, :force_on_endpoint),
        "force_on_endpoint is deprecated and will be removed in SpaceAGORA 1.0; use " *
        "link_force(link, a, b, ctx), which returns the forces on both endpoints.")
    return _endpoint_force(current, other, link)
end

# The force on `current`, without the deprecation: for callers that still pass
# a `candidate_force` of the old form.
_endpoint_force(current, other, link) =
    first(link_force(link, link_endpoint(current, nothing, 0), link_endpoint(other, nothing, 0), LinkContext()))

"""
    laser_force_on_spacecraft(model, state, satellite; is_active=nothing)

Deprecated, use `interlink_force_on_spacecraft`, which sums every selected link
of any kind. The old name will be removed in SpaceAGORA 1.0.
"""
function laser_force_on_spacecraft(model::InterlinkModel, state, satellite::Int; is_active=nothing)
    _interlink_deprecation((InterlinkModel, :laser_force_on_spacecraft),
        "laser_force_on_spacecraft is deprecated and will be removed in SpaceAGORA 1.0; use " *
        "interlink_force_on_spacecraft(model, state, satellite; is_active).")
    return interlink_force_on_spacecraft(model, state, satellite; is_active)
end

end
