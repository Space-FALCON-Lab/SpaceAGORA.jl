const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function _active_source(path::String)
    src = read(path, String)
    return join(
        (
            strip(first(split(line, '#', limit=2)))
            for line in split(src, '\n')
            if !isempty(strip(first(split(line, '#', limit=2))))
        ),
        '\n',
    )
end

const COMMAND_TYPES_PATH = joinpath(REPO_ROOT, "src", "models", "gnc", "contracts", "command_types.jl")
const GUIDANCE_HOOKS_PATH = joinpath(REPO_ROOT, "src", "models", "gnc", "guidance", "hooks.jl")
const GUIDANCE_THRUSTER_PATH = joinpath(REPO_ROOT, "src", "models", "gnc", "guidance", "maneuvers", "thruster_guidance_functions.jl")
const CONTROL_HOOKS_PATH = joinpath(REPO_ROOT, "src", "models", "gnc", "control", "hooks.jl")
const CONTROL_PROPULSIVE_PATH = joinpath(REPO_ROOT, "src", "models", "gnc", "control", "maneuvers", "propulsive_maneuvers.jl")

for path in (
    COMMAND_TYPES_PATH,
    GUIDANCE_HOOKS_PATH,
    GUIDANCE_THRUSTER_PATH,
    CONTROL_HOOKS_PATH,
    CONTROL_PROPULSIVE_PATH,
)
    isfile(path) || error("Missing typed-command boundary file: $(relpath(path, REPO_ROOT))")
end

command_types_src = _active_source(COMMAND_TYPES_PATH)
guidance_hooks_src = _active_source(GUIDANCE_HOOKS_PATH)
guidance_thruster_src = _active_source(GUIDANCE_THRUSTER_PATH)
control_hooks_src = _active_source(CONTROL_HOOKS_PATH)
control_propulsive_src = _active_source(CONTROL_PROPULSIVE_PATH)

occursin("struct PropulsiveManeuverCommand", command_types_src) ||
    error("Typed-command boundary violation: missing PropulsiveManeuverCommand owner.")
occursin("struct AerobrakingControlCommand", command_types_src) ||
    error("Typed-command boundary violation: missing AerobrakingControlCommand owner.")
occursin("using ..CommandTypes: PropulsiveManeuverCommand, AerobrakingControlCommand", guidance_hooks_src) ||
    error("Typed-command boundary violation: GuidanceHooks is not importing typed commands from CommandTypes.")
occursin("using ..CommandTypes: PropulsiveManeuverCommand", control_hooks_src) ||
    error("Typed-command boundary violation: ControlHooks is not importing PropulsiveManeuverCommand from CommandTypes.")
occursin("PropulsiveManeuverCommand(", guidance_thruster_src) ||
    error("Typed-command boundary violation: guidance thruster logic is not constructing PropulsiveManeuverCommand values.")
occursin("p.shared_buffers.maneuver_commands[i] = command", guidance_thruster_src) ||
    error("Typed-command boundary violation: guidance thruster logic is not writing maneuver commands into shared_buffers.")
occursin("commands = p.shared_buffers.maneuver_commands", control_propulsive_src) ||
    error("Typed-command boundary violation: control propulsive maneuvers are not reading the typed maneuver command buffer.")
occursin("PropulsiveBurnPlan(", control_propulsive_src) ||
    error("Typed-command boundary violation: control propulsive maneuvers are not constructing a separate typed burn plan.")
occursin("maneuver_burn_plans", control_propulsive_src) ||
    error("Typed-command boundary violation: control propulsive maneuvers are not storing actuation state in the maneuver_burn_plans buffer.")
occursin("command.delta_v_mps", control_propulsive_src) ||
    error("Typed-command boundary violation: control propulsive maneuvers are not consuming command.delta_v_mps.")
occursin("command.direction_rad", control_propulsive_src) ||
    error("Typed-command boundary violation: control propulsive maneuvers are not consuming command.direction_rad.")
occursin("controlModel.Δv[i] = command.delta_v_mps", control_propulsive_src) &&
    error("Typed-command boundary violation: control propulsive maneuvers are still copying guidance delta-v into BaseThrusterModel state.")
occursin("controlModel.direction[i] = command.direction_rad", control_propulsive_src) &&
    error("Typed-command boundary violation: control propulsive maneuvers are still copying guidance direction into BaseThrusterModel state.")

const GUIDANCE_ROOT = joinpath(REPO_ROOT, "src", "models", "gnc", "guidance")
const CONTROL_ROOT = joinpath(REPO_ROOT, "src", "models", "gnc", "control")
const GUIDANCE_ROOTS = (GUIDANCE_ROOT,
    (joinpath(REPO_ROOT, "src", "models", "gnc", owner) for owner in ("planning", "companions"))...,
    joinpath(REPO_ROOT, "src", "interfaces", "analysis", "planning"))
# The general contract serves guidance and control alike: the GNC-wide checks
# apply to it, the guidance-only patterns do not.
const CONTRACTS_ROOT = joinpath(REPO_ROOT, "src", "models", "gnc", "contracts")
const GNC_ROOTS = (GUIDANCE_ROOTS..., CONTROL_ROOT, CONTRACTS_ROOT)
const ALLOWED_MANEUVER_COMMAND_OWNERS = Set([
    joinpath("src", "models", "gnc", "guidance", "maneuvers", "thruster_guidance_functions.jl"),
    joinpath("src", "models", "gnc", "control", "maneuvers", "propulsive_maneuvers.jl"),
])
const GUIDANCE_FORBIDDEN_PATTERNS = (
    r"\bBaseThrusterModel\b" => "guidance must not depend on BaseThrusterModel directly",
    r"\bAbstractThrusterModel\b" => "guidance must not depend on AbstractThrusterModel directly",
    r"\bThrusterModels\b" => "guidance must not depend on ThrusterModels directly",
    r"\.\s*Δv\s*\[" => "guidance must not mutate thruster delta-v state directly",
    r"\.\s*direction\s*\[" => "guidance must not mutate thruster direction state directly",
    r"\.\s*start_burn_time\s*\[" => "guidance must not mutate thruster burn scheduling directly",
    r"\.\s*stop_burn_time\s*\[" => "guidance must not mutate thruster burn scheduling directly",
    r"\.\s*thrust\s*\[" => "guidance must not mutate thruster thrust state directly",
)

# Comparison execution owns the concrete six-axis actuator used by its plant.
# Permit that one binding only; every other thruster dependency/state access and
# the GNC-wide checks below still apply to this file.
function _guidance_scan_source(rel::String, src::String)
    rel == joinpath("src", "interfaces", "analysis", "planning", "rpo_comparison_tracking.jl") || return src
    return replace(src, r"(?m)^using \.\.ThrusterModels: SixAxisThrusterModel$" => ""; count=1)
end

violations = String[]

for root in GNC_ROOTS
    isdir(root) || error("Missing typed GNC scan root: $(relpath(root, REPO_ROOT))")
    for (dir, _, files) in walkdir(root)
        for file in files
            endswith(file, ".jl") || continue
            path = joinpath(dir, file)
            rel = relpath(path, REPO_ROOT)
            src = _active_source(path)

            occursin("DynamicEffectors", src) &&
                push!(violations, "$rel: GNC source still depends on DynamicEffectors directly")

            if occursin("maneuver_commands", src) && !(rel in ALLOWED_MANEUVER_COMMAND_OWNERS)
                push!(violations, "$rel: maneuver_commands ownership is restricted to typed guidance writer/control consumer files")
            end

            if root in GUIDANCE_ROOTS
                guidance_src = _guidance_scan_source(rel, src)
                for (rx, message) in GUIDANCE_FORBIDDEN_PATTERNS
                    occursin(rx, guidance_src) || continue
                    push!(violations, "$rel: $message")
                end
            end
        end
    end
end

if !isempty(violations)
    println("gnc_typed_command_boundary_violations:")
    for violation in sort(unique(violations))
        println("  - $violation")
    end
    error("GNC typed-command boundary gate failed")
end

println("gnc_typed_command_boundary_gate_ok")
