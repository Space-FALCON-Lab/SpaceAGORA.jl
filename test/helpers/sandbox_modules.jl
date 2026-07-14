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
