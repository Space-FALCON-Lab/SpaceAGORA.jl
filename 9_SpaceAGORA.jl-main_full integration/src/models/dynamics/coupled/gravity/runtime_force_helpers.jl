@inline function _gravity_runtime_field(args, name::Symbol, default)
    if args !== nothing && hasproperty(args, name)
        return getproperty(args, name)
    end
    if args !== nothing && applicable(get, args, name, default)
        return get(args, name, default)
    end
    if args !== nothing && applicable(getindex, args, name)
        return getindex(args, name)
    end
    return default
end

# Compiled once, never inlined, called by every route: the flat pre-pass and
# the per-satellite wrench must not each get their own codegen of the same
# arithmetic (StaticArrays products use `muladd`, whose fusion is decided per
# compilation context).
@noinline function _inverse_squared_force_ii(pos_ii::SVector{3, Float64}, mass::Float64, planet)::SVector{3, Float64}
    return mass * _inverse_squared_gravity_accel(pos_ii, planet)
end

@noinline function _inverse_squared_j2_force_ii(pos_ii::SVector{3, Float64}, mass::Float64, l_pi::SMatrix{3, 3, Float64, 9}, planet)::SVector{3, Float64}
    pos_pp = SVector{3, Float64}(l_pi * pos_ii)
    gravity_pp = _inverse_squared_j2_gravity_accel(pos_pp, planet)
    return mass * (l_pi' * gravity_pp)
end
