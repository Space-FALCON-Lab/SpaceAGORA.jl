module ThrusterModels

using ..AbstractTypes: AbstractThrusterModel
using StaticArrays
using LinearAlgebra

export BaseThrusterModel, SixAxisThrusterModel

include(joinpath(@__DIR__, "thruster_models.jl"))

end # module ThrusterModels
