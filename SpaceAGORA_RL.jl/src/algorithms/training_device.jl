abstract type AbstractTrainingDevice end

struct CPUTrainingDevice <: AbstractTrainingDevice end

struct CUDATrainingDevice <: AbstractTrainingDevice
    device_id::Int
end

training_device_name(::CPUTrainingDevice) = "cpu"
training_device_name(device::CUDATrainingDevice) = "cuda:$(device.device_id)"

function _cuda_extension()
    return Base.get_extension(@__MODULE__, :SpaceAGORA_RLCUDAExt)
end

function cuda_functional()
    ext = _cuda_extension()
    ext === nothing && return false
    return Bool(Base.invokelatest(getproperty(ext, :cuda_functional)))
end

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

function to_cpu_array(x::AbstractArray)
    ext = _cuda_extension()
    if ext !== nothing && isdefined(ext, :maybe_to_cpu_array)
        converted = Base.invokelatest(getproperty(ext, :maybe_to_cpu_array), x)
        converted === nothing || return converted
    end
    return Base.invokelatest(Array, x)
end

to_device_array(::CPUTrainingDevice, x::AbstractArray{Float32}) = to_cpu_array(x)
to_device_array(::CPUTrainingDevice, x::AbstractArray{<:Real}) = Float32.(to_cpu_array(x))
to_device_array(::CPUTrainingDevice, x::AbstractArray{Bool}) = to_cpu_array(x)

function to_device_array(::CUDATrainingDevice, x)
    ext = _cuda_extension()
    ext === nothing && throw(ArgumentError("CUDA extension is not loaded; run with CUDA.jl available or use device=\"cpu\""))
    return Base.invokelatest(getproperty(ext, :to_cuda_device_array), x)
end

cpu_scalar(x::Number) = x
cpu_scalar(x) = only(to_cpu_array(x))

zero_array_like(x::AbstractArray{Float32}) = fill!(similar(x), 0f0)
