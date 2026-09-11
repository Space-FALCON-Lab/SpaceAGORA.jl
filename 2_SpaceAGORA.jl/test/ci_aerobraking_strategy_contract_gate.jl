const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

function _read(rel)
    path = joinpath(REPO_ROOT, rel)
    isfile(path) || error("Missing required file: $rel")
    return read(path, String)
end

interfaces = _read(joinpath("src", "gnc", "guidance", "aerobraking", "interfaces.jl"))
dispatcher = _read(joinpath("src", "gnc", "guidance", "aerobraking", "dispatcher.jl"))
command_types = _read(joinpath("src", "gnc", "command_types.jl"))
policy_types = _read(joinpath("src", "mission", "operations", "aerobraking_policy", "policy_types.jl"))
selector_stub = _read(joinpath("src", "mission", "operations", "aerobraking_policy", "selector_stub.jl"))

for token in (
    "abstract type AbstractAerobrakingStrategy",
    "struct AerobrakingGuidanceInput",
    "struct AerobrakingGuidanceOutput",
)
    occursin(token, interfaces) || error("Missing strategy interface token: $token")
end

for token in (
    "struct AerobrakingControlCommand",
    "struct PropulsiveManeuverCommand",
)
    occursin(token, command_types) || error("Missing command-type token: $token")
end

occursin("dispatch_aerobraking_guidance", dispatcher) ||
    error("Missing dispatcher entrypoint dispatch_aerobraking_guidance")
occursin("strategy_from_kind", dispatcher) ||
    error("Missing strategy_from_kind mapping")

occursin("@enum AerobrakingStrategyKind", policy_types) ||
    error("Missing AerobrakingStrategyKind enum")
occursin("AbstractAerobrakingPolicySelector", policy_types) ||
    error("Missing AbstractAerobrakingPolicySelector contract")
occursin("DRLPolicyAdapterStub", selector_stub) ||
    error("Missing DRLPolicyAdapterStub selector contract")
occursin("Not implemented: DRLPolicyAdapterStub.select_strategy", selector_stub) ||
    error("DRLPolicyAdapterStub must have explicit Not implemented behavior")
!occursin("AbstractDict", selector_stub) ||
    error("DefaultAerobrakingPolicySelector must not inspect AbstractDict fallback inputs")
!occursin(":aerobraking_strategy", selector_stub) ||
    error("DefaultAerobrakingPolicySelector must not read strategy tokens from input.args")

println("aerobraking_strategy_contract_gate_ok")
