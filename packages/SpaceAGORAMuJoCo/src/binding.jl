"""
The only file in this package that talks to libmujoco.

It loads MuJoCo 3.11.0 from the package's lazy artifact (official release tarball, checksum pinned in
`Artifacts.toml`), resolves the C entry points it needs, and exposes raw-pointer field access at byte
offsets taken from the 3.11.0 headers with `offsetof` (`mujoco/mjmodel.h`, `mjdata.h`, `mjoption.h`).
Those offsets are valid for exactly that version; `load!` refuses any other `mj_version()`. Moving to
another MuJoCo release means regenerating `OFFSETS` and nothing else.

Pointers are plain `Ptr{Cvoid}`; ownership lives in `MjModel` and `MjData`, whose finalizers free the
native memory.
"""
module Binding

using Libdl
using Artifacts
using LazyArtifacts

const MUJOCO_VERSION = 3_011_000

# offsetof() values from the 3.11.0 headers (x86_64, mjtNum = double, mjtBool = unsigned char).
const OFFSETS = (
    model = (nq = 0, nv = 8, nu = 16, na = 40, nbody = 48, njnt = 88, neq = 472,
        body_rootid = 1816, body_jntnum = 1840, body_jntadr = 1848, body_dofnum = 1856,
        jnt_type = 2104, jnt_qposadr = 2112, jnt_dofadr = 2120,
        body_mass = 1944, body_inertia = 1960, opt = 792),
    opt = (timestep = 0, gravity = 56, integrator = 248),
    data = (warning = 160264, time = 160640, qpos = 160680, qvel = 160688, act = 160696, ctrl = 160728,
        xfrc_applied = 160744, eq_active = 160752, xquat = 160840, xipos = 160856, ximat = 160864),
)
# mjWarningStat is {int lastinfo; int number}; mjtWarning enum values from mujoco/mjtype.h, mjNWARNING = 7.
const WARNING_STRIDE = 8
const WARNING_NUMBER = 4
const WARN_BADQPOS = 3
const WARN_BADQVEL = 4
const WARN_BADQACC = 5
const SIZEOF_MJMODEL = 5672
const SIZEOF_MJDATA = 161976
const SIZEOF_MJOPTION = 304

# enums (mjtIntegrator, mjtObj, mjtJoint, mjtState)
const INT_EULER = Cint(0)
const INT_RK4 = Cint(1)
const INT_IMPLICIT = Cint(2)
const INT_IMPLICITFAST = Cint(3)
const OBJ_BODY = Cint(1)
const JNT_FREE = Cint(0)
const STATE_INTEGRATION = Cint(1 << 0 | 1 << 1 | 1 << 2 | 1 << 3 | 1 << 4 | 1 << 13 |            # FULLPHYSICS
    1 << 6 | 1 << 7 | 1 << 8 | 1 << 9 | 1 << 10 | 1 << 11 | 1 << 12 | 1 << 5)                    # USER + WARMSTART

const FUNCTIONS = (:mj_version, :mj_loadXML, :mj_parseXMLString, :mj_compile, :mjs_getError, :mj_deleteSpec,
    :mj_makeData, :mj_deleteData, :mj_deleteModel, :mj_copyModel, :mj_step1, :mj_step2, :mj_forward,
    :mj_resetData, :mj_stateSize, :mj_getState, :mj_setState, :mj_name2id, :mj_id2name, :mj_objectVelocity, :mj_getTotalmass)

struct Api
    handle::Ptr{Cvoid}
    fn::NamedTuple{FUNCTIONS, NTuple{length(FUNCTIONS), Ptr{Cvoid}}}
end

const _API = Ref{Union{Nothing, Api}}(nothing)   # the loaded library, process-wide by nature
const _API_LOCK = ReentrantLock()

"""Path of `libmujoco`; throws on platforms without a pinned official asset."""
function library_path()
    (Sys.islinux() && Sys.ARCH === :x86_64) || throw(ErrorException(
        "SpaceAGORAMuJoCo: unsupported platform $(Sys.KERNEL)/$(Sys.ARCH). MuJoCo 3.11.0 is pinned for Linux x86_64 only."))
    root = artifact"mujoco_3_11_0_linux_x86_64"
    for dir in (joinpath(root, "lib"), joinpath(root, "mujoco-3.11.0", "lib"))
        isdir(dir) || continue
        for name in readdir(dir)
            startswith(name, "libmujoco.so.") && return joinpath(dir, name)
        end
    end
    error("libmujoco.so.* not found under the artifact at $root")
end

function _load()
    handle = Libdl.dlopen(library_path(), Libdl.RTLD_LOCAL | Libdl.RTLD_NOW)
    fn = NamedTuple{FUNCTIONS}(Tuple(Libdl.dlsym(handle, f) for f in FUNCTIONS))
    version = ccall(fn.mj_version, Cint, ())
    version == MUJOCO_VERSION || error("SpaceAGORAMuJoCo needs MuJoCo $MUJOCO_VERSION, loaded $version")
    return Api(handle, fn)
end

function api()::Api
    a = _API[]
    a === nothing || return a
    lock(_API_LOCK) do
        _API[] === nothing && (_API[] = _load())
        return _API[]::Api
    end
end

version() = Int(ccall(api().fn.mj_version, Cint, ()))

# --- native memory owners ---------------------------------------------------------------------------

mutable struct MjModel
    ptr::Ptr{Cvoid}
    function MjModel(ptr::Ptr{Cvoid})
        ptr == C_NULL && error("MuJoCo returned a null model")
        self = new(ptr)
        finalizer(self) do m
            m.ptr == C_NULL || ccall(api().fn.mj_deleteModel, Cvoid, (Ptr{Cvoid},), m.ptr)
            m.ptr = C_NULL
        end
        return self
    end
end

mutable struct MjData
    ptr::Ptr{Cvoid}
    function MjData(model::MjModel)
        ptr = ccall(api().fn.mj_makeData, Ptr{Cvoid}, (Ptr{Cvoid},), model.ptr)
        ptr == C_NULL && error("mj_makeData failed")
        self = new(ptr)
        finalizer(self) do d
            d.ptr == C_NULL || ccall(api().fn.mj_deleteData, Cvoid, (Ptr{Cvoid},), d.ptr)
            d.ptr = C_NULL
        end
        return self
    end
end

# --- loading ----------------------------------------------------------------------------------------

function load_xml_file(path::AbstractString)
    err = zeros(UInt8, 1024)
    ptr = GC.@preserve err ccall(api().fn.mj_loadXML, Ptr{Cvoid}, (Cstring, Ptr{Cvoid}, Ptr{UInt8}, Cint), path, C_NULL, err, 1024)
    ptr == C_NULL && error("mj_loadXML failed: ", unsafe_string(pointer(err)))
    return MjModel(ptr)
end

function load_xml_string(xml::AbstractString)
    f = api().fn
    err = zeros(UInt8, 1024)
    spec = GC.@preserve err ccall(f.mj_parseXMLString, Ptr{Cvoid}, (Cstring, Ptr{Cvoid}, Ptr{UInt8}, Cint), xml, C_NULL, err, 1024)
    spec == C_NULL && error("mj_parseXMLString failed: ", unsafe_string(pointer(err)))
    try
        ptr = ccall(f.mj_compile, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}), spec, C_NULL)
        ptr == C_NULL && error("mj_compile failed: ", unsafe_string(ccall(f.mjs_getError, Cstring, (Ptr{Cvoid},), spec)))
        return MjModel(ptr)
    finally
        ccall(f.mj_deleteSpec, Cvoid, (Ptr{Cvoid},), spec)
    end
end

copy_model(model::MjModel) = MjModel(ccall(api().fn.mj_copyModel, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}), C_NULL, model.ptr))

# --- stepping and state -----------------------------------------------------------------------------

step1!(m::MjModel, d::MjData) = ccall(api().fn.mj_step1, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), m.ptr, d.ptr)
step2!(m::MjModel, d::MjData) = ccall(api().fn.mj_step2, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), m.ptr, d.ptr)
forward!(m::MjModel, d::MjData) = ccall(api().fn.mj_forward, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), m.ptr, d.ptr)
reset_data!(m::MjModel, d::MjData) = ccall(api().fn.mj_resetData, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), m.ptr, d.ptr)
state_size(m::MjModel, sig::Integer=STATE_INTEGRATION) = Int(ccall(api().fn.mj_stateSize, Cint, (Ptr{Cvoid}, Cint), m.ptr, sig))
get_state!(buf::Vector{Float64}, m::MjModel, d::MjData, sig::Integer=STATE_INTEGRATION) =
    (ccall(api().fn.mj_getState, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Float64}, Cint), m.ptr, d.ptr, buf, sig); buf)
set_state!(m::MjModel, d::MjData, buf::Vector{Float64}, sig::Integer=STATE_INTEGRATION) =
    (ccall(api().fn.mj_setState, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Float64}, Cint), m.ptr, d.ptr, buf, sig); nothing)
name2id(m::MjModel, objtype::Integer, name::AbstractString) =
    Int(ccall(api().fn.mj_name2id, Cint, (Ptr{Cvoid}, Cint, Cstring), m.ptr, objtype, name))
function id2name(m::MjModel, objtype::Integer, id::Integer)
    p = ccall(api().fn.mj_id2name, Cstring, (Ptr{Cvoid}, Cint, Cint), m.ptr, objtype, id)
    return p == C_NULL ? "" : unsafe_string(p)
end
total_mass(m::MjModel) = ccall(api().fn.mj_getTotalmass, Float64, (Ptr{Cvoid},), m.ptr)

"""6D velocity (angular; linear) of a body about its center of mass, in world axes. Needs fresh `cvel`."""
function body_velocity!(res::Vector{Float64}, m::MjModel, d::MjData, body::Integer)
    ccall(api().fn.mj_objectVelocity, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Cint, Ptr{Float64}, Cint), m.ptr, d.ptr, OBJ_BODY, body, res, 0)
    return res
end

# --- field access -----------------------------------------------------------------------------------

@inline _int(m::MjModel, off) = Int(unsafe_load(Ptr{Cint}(m.ptr + off)))
@inline _arr(base::Ptr{Cvoid}, off, ::Type{T}, n) where {T} =
    n == 0 ? T[] : unsafe_wrap(Array, unsafe_load(Ptr{Ptr{T}}(base + off)), n; own = false)

nq(m::MjModel) = _int(m, OFFSETS.model.nq)
nv(m::MjModel) = _int(m, OFFSETS.model.nv)
nu(m::MjModel) = _int(m, OFFSETS.model.nu)
na(m::MjModel) = _int(m, OFFSETS.model.na)
nbody(m::MjModel) = _int(m, OFFSETS.model.nbody)
njnt(m::MjModel) = _int(m, OFFSETS.model.njnt)
neq(m::MjModel) = _int(m, OFFSETS.model.neq)

body_mass(m::MjModel) = _arr(m.ptr, OFFSETS.model.body_mass, Float64, nbody(m))
body_inertia(m::MjModel) = _arr(m.ptr, OFFSETS.model.body_inertia, Float64, 3 * nbody(m))
body_rootid(m::MjModel) = _arr(m.ptr, OFFSETS.model.body_rootid, Cint, nbody(m))
body_jntnum(m::MjModel) = _arr(m.ptr, OFFSETS.model.body_jntnum, Cint, nbody(m))
body_jntadr(m::MjModel) = _arr(m.ptr, OFFSETS.model.body_jntadr, Cint, nbody(m))
body_dofnum(m::MjModel) = _arr(m.ptr, OFFSETS.model.body_dofnum, Cint, nbody(m))
jnt_type(m::MjModel) = _arr(m.ptr, OFFSETS.model.jnt_type, Cint, njnt(m))
jnt_qposadr(m::MjModel) = _arr(m.ptr, OFFSETS.model.jnt_qposadr, Cint, njnt(m))
jnt_dofadr(m::MjModel) = _arr(m.ptr, OFFSETS.model.jnt_dofadr, Cint, njnt(m))

timestep(m::MjModel) = unsafe_load(Ptr{Float64}(m.ptr + OFFSETS.model.opt + OFFSETS.opt.timestep))
set_timestep!(m::MjModel, dt::Float64) = unsafe_store!(Ptr{Float64}(m.ptr + OFFSETS.model.opt + OFFSETS.opt.timestep), dt)
integrator(m::MjModel) = Int(unsafe_load(Ptr{Cint}(m.ptr + OFFSETS.model.opt + OFFSETS.opt.integrator)))
set_integrator!(m::MjModel, i::Integer) = unsafe_store!(Ptr{Cint}(m.ptr + OFFSETS.model.opt + OFFSETS.opt.integrator), Cint(i))
gravity(m::MjModel) = ntuple(i -> unsafe_load(Ptr{Float64}(m.ptr + OFFSETS.model.opt + OFFSETS.opt.gravity), i), 3)
set_gravity!(m::MjModel, g) = for i in 1:3
    unsafe_store!(Ptr{Float64}(m.ptr + OFFSETS.model.opt + OFFSETS.opt.gravity), Float64(g[i]), i)
end

data_time(d::MjData) = unsafe_load(Ptr{Float64}(d.ptr + OFFSETS.data.time))
qpos(m::MjModel, d::MjData) = _arr(d.ptr, OFFSETS.data.qpos, Float64, nq(m))
qvel(m::MjModel, d::MjData) = _arr(d.ptr, OFFSETS.data.qvel, Float64, nv(m))
ctrl(m::MjModel, d::MjData) = _arr(d.ptr, OFFSETS.data.ctrl, Float64, nu(m))
xfrc_applied(m::MjModel, d::MjData) = unsafe_wrap(Array, unsafe_load(Ptr{Ptr{Float64}}(d.ptr + OFFSETS.data.xfrc_applied)), (6, nbody(m)); own = false)
xipos(m::MjModel, d::MjData) = unsafe_wrap(Array, unsafe_load(Ptr{Ptr{Float64}}(d.ptr + OFFSETS.data.xipos)), (3, nbody(m)); own = false)
xquat(m::MjModel, d::MjData) = unsafe_wrap(Array, unsafe_load(Ptr{Ptr{Float64}}(d.ptr + OFFSETS.data.xquat)), (4, nbody(m)); own = false)
"""Times MuJoCo raised warning `w` (a `WARN_*` index) since the counter was last cleared."""
warning_count(d::MjData, w::Integer) = Int(unsafe_load(Ptr{Cint}(d.ptr + OFFSETS.data.warning + WARNING_STRIDE * w + WARNING_NUMBER)))
eq_active(m::MjModel, d::MjData) = _arr(d.ptr, OFFSETS.data.eq_active, UInt8, neq(m))

end # module Binding
