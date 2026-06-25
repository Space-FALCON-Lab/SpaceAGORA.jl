abstract type AbstractTrainingDevice end

struct CPUTrainingDevice <: AbstractTrainingDevice end

struct CUDATrainingDevice <: AbstractTrainingDevice
    device_id::Int
end

training_device_name(::CPUTrainingDevice) = "cpu"
training_device_name(device::CUDATrainingDevice) = "cuda:$(device.device_id)"

cuda_functional() = false

function try_load_cuda()
    Base.find_package("CUDA") === nothing && return false
    try
        Base.require(Base.PkgId(Base.UUID("052768ef-5323-5732-b1bb-66c8b64840ba"), "CUDA"))
    catch
        return false
    end
    return cuda_functional()
end

function resolve_training_device(name::Union{AbstractString,Symbol})
    key = Symbol(lowercase(String(name)))
    if key == :cpu
        return CPUTrainingDevice()
    elseif key == :auto
        return try_load_cuda() ? CUDATrainingDevice(0) : CPUTrainingDevice()
    elseif key == :cuda || key == :gpu
        try_load_cuda() || throw(ArgumentError("CUDA training requested, but CUDA.jl is not available or not functional"))
        return CUDATrainingDevice(0)
    end
    throw(ArgumentError("unknown training device: $name"))
end

to_device_array(::CPUTrainingDevice, x::AbstractArray{Float32}) = Array(x)
to_device_array(::CPUTrainingDevice, x::AbstractArray{<:Real}) = Float32.(Array(x))
to_device_array(::CPUTrainingDevice, x::AbstractArray{Bool}) = Array(x)

function to_device_array(::CUDATrainingDevice, x)
    throw(ArgumentError("CUDA extension is not loaded; run with CUDA.jl available or use device=\"cpu\""))
end

to_cpu_array(x::AbstractArray) = Array(x)
cpu_scalar(x::Number) = x
cpu_scalar(x) = only(to_cpu_array(x))

zero_array_like(x::AbstractArray{Float32}) = fill!(similar(x), 0f0)
