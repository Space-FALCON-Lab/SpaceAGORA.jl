"""Compatibility entry point for HYPR's optional SpaceAGORA adapter."""
module SpaceAGORAHYPR
import SpaceAGORA
import HYPR
const Adapter = Base.get_extension(HYPR, :HYPRSpaceAGORAExt)
Adapter === nothing && error("HYPR's SpaceAGORA extension did not load.")
const SwarmPolicy = Adapter.SwarmPolicy
const RPO = Adapter.RPO
const RobotArm = Adapter.RobotArm
const PlannerAdapter = Adapter.PlannerAdapter
end
