# The two EDG decisions SpaceAGORA calls at run time. The standalone companion
# owns them and every numerical kernel behind them.
module EDGAlgorithms
using ..EDGServices: _invoke_kernel

control_decision!(args...; kwargs...) = _invoke_kernel(:control_decision!, args...; kwargs...)
guidance_decision!(args...; kwargs...) = _invoke_kernel(:guidance_decision!, args...; kwargs...)

end # module EDGAlgorithms
