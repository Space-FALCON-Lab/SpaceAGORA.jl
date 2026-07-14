module IncludeOrderSandbox
end
const INCLUDE_ORDER_SANDBOX = IncludeOrderSandbox

module ExportImportSandbox
end
const EXPORT_IMPORT_SANDBOX = ExportImportSandbox

const GUIDANCE_SANDBOX_LOADED = Ref(false)
module GuidanceSandbox
using ..SimulationModel
using ..SimulationModel.AbstractTypes: AbstractGuidanceModel
using ComponentArrays
using LinearAlgebra
using StaticArrays
end
const GUIDANCE_SANDBOX = GuidanceSandbox

function ensure_guidance_sandbox_loaded!()
    if GUIDANCE_SANDBOX_LOADED[]
        return
    end

    Core.eval(GUIDANCE_SANDBOX, :(include(joinpath(Main.REPO_ROOT, "src", "gnc", "guidance", "thruster_guidance", "thruster_guidance_models.jl"))))
    Core.eval(GUIDANCE_SANDBOX, :(include(joinpath(Main.REPO_ROOT, "src", "gnc", "guidance", "thruster_guidance", "thruster_guidance_functions.jl"))))
    GUIDANCE_SANDBOX_LOADED[] = true
end
