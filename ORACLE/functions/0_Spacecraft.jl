# Builds a rigid single-body SpacecraftModel from an id, mass, and initial condition.
function _spacecraft(id::Int, mass_kg::Float64, ic::InitialCondition)
    bus = Link(root=true, m=mass_kg, ref_area=1.0)
    return SpacecraftModel(
        Joint[],
        [bus],
        bus,
        true,
        mass_kg,
        0.0,
        bus.inertia,
        0,
        0,
        ic,
        id,
    )
end
