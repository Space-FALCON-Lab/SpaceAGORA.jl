module SpaceAGORA_RLCUDAExt

using CUDA
using SpaceAGORA_RL

SpaceAGORA_RL.cuda_functional() = CUDA.functional()

function SpaceAGORA_RL.to_device_array(::SpaceAGORA_RL.CUDATrainingDevice,
                                       x::AbstractArray{Float32})
    return CUDA.cu(x)
end

function SpaceAGORA_RL.to_device_array(::SpaceAGORA_RL.CUDATrainingDevice,
                                       x::AbstractArray{<:Real})
    return CUDA.cu(Float32.(x))
end

function SpaceAGORA_RL.to_device_array(::SpaceAGORA_RL.CUDATrainingDevice,
                                       x::AbstractArray{Bool})
    return CUDA.cu(x)
end

SpaceAGORA_RL.to_cpu_array(x::CUDA.CuArray) = Array(x)
SpaceAGORA_RL.cpu_scalar(x::CUDA.CuArray) = only(Array(x))

end
