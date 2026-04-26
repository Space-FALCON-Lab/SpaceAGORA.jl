module ThrusterModels

using ..AbstractTypes: AbstractThrusterModel

export BaseThrusterModel

include(joinpath(@__DIR__, "thruster_models.jl"))

end # module ThrusterModels
