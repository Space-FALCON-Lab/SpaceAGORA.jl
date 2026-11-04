"""Compatibility entry point for HYPR's optional SpaceAGORA adapter."""
module SpaceAGORAHYPR
import SpaceAGORA
import HYPR
const REQUIRED_BINDINGS = (:SwarmPolicy, :RPO, :RobotArm, :PlannerAdapter, :initialized)
function checked_adapter(adapter; require_initialized::Bool)
    complete = adapter isa Module && all(name -> isdefined(adapter, name), REQUIRED_BINDINGS)
    complete || throw(SpaceAGORA.HYPRServices.CompatibilityError(
        "HYPR's SpaceAGORA extension is missing or incomplete. Start a fresh process with the supported package pair."))
    if require_initialized && !(adapter.initialized() && SpaceAGORA.hypr_available())
        throw(SpaceAGORA.HYPRServices.CompatibilityError(
            "HYPR's SpaceAGORA extension did not initialize successfully. Discard this process and install the supported package pair."))
    end
    return adapter
end
# Bindings are checked before aliases; runtime initialization is checked again
# when this package image is restored in a fresh process.
const Adapter = checked_adapter(Base.get_extension(HYPR, :HYPRSpaceAGORAExt); require_initialized=false)
const SwarmPolicy = Adapter.SwarmPolicy
const RPO = Adapter.RPO
const RobotArm = Adapter.RobotArm
const PlannerAdapter = Adapter.PlannerAdapter
function __init__()
    checked_adapter(Adapter; require_initialized=true)
end
end
