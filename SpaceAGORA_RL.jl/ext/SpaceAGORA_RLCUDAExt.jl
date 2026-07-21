module SpaceAGORA_RLCUDAExt

using CUDA
using SpaceAGORA_RL

cuda_functional() = CUDA.functional()

to_cuda_device_array(x::AbstractArray{Float32}) = CUDA.cu(x)
to_cuda_device_array(x::AbstractArray{<:Real}) = CUDA.cu(Float32.(x))
to_cuda_device_array(x::AbstractArray{Bool}) = CUDA.cu(x)
maybe_to_cpu_array(x) = nothing
maybe_to_cpu_array(x::CUDA.CuArray) = Base.invokelatest(Array, x)

function SpaceAGORA_RL.to_device_array(::SpaceAGORA_RL.CUDATrainingDevice,
                                       x::AbstractArray{Float32})
    return to_cuda_device_array(x)
end

function SpaceAGORA_RL.to_device_array(::SpaceAGORA_RL.CUDATrainingDevice,
                                       x::AbstractArray{<:Real})
    return to_cuda_device_array(x)
end

function SpaceAGORA_RL.to_device_array(::SpaceAGORA_RL.CUDATrainingDevice,
                                       x::AbstractArray{Bool})
    return to_cuda_device_array(x)
end

SpaceAGORA_RL.to_cpu_array(x::CUDA.CuArray) = maybe_to_cpu_array(x)
SpaceAGORA_RL.cpu_scalar(x::CUDA.CuArray) = only(maybe_to_cpu_array(x))

end
