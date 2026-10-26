"""Shared once-per-key warnings: deprecated names and pre-run configuration checks."""
module HookDeprecation

"""
    _hook_deprecation(key::Tuple, msg::String)

Warn once per process about a deprecated name, with `key` naming its owner and
the old name, for example `(StateAnchor, :sat_idx)`. Under
`julia --depwarn=error` it throws an `ErrorException` with `msg` instead, so a
test run can refuse deprecated names.

Its callers are deprecated names in use: an old keyword passed, or an old
property read or written (the renamed index names below, for example).
`Base.depwarn` is not used: its once-only filter is keyed on the call site, so
one call site shared by several owners would warn about the first owner only,
and it is silent unless Julia runs with `--depwarn=yes`.
"""
function _hook_deprecation(key::Tuple, msg::String)
    return _preflight_warning(key, msg; group=:depwarn)
end

"""
    _preflight_warning(key::Tuple, msg::String; group::Symbol=:preflight)

Emit a pre-run warning once per key through the active logger. Refuse the
configuration under `--depwarn=error`, even if the logger has already seen the
key. Non-deprecation configuration checks use the default `:preflight` group.
"""
function _preflight_warning(key::Tuple, msg::String; group::Symbol=:preflight)
    Base.JLOptions().depwarn == 2 && throw(ErrorException(msg))
    @warn msg _group=group _id=key maxlog=1
    return nothing
end

# Old-style index names deprecated in 0.5.0 (`StateAnchor.sat_idx` and the
# `sat_idx` keyword of `get_periapsis_pulse_callback`; `SchedulingPolicyModel`'s
# `target_idx` and `target_idxs` were removed in 0.6.0): an old keyword or property
# warns through `_hook_deprecation`, keyed by `(owner, old name)`, so the warning
# shows once per owner and name. The owner is the type or the function that had
# the name.
_renamed_name(owner::Type, old::Symbol) = "$(nameof(owner)).$(old)"
_renamed_name(owner::Function, old::Symbol) = "The $(old) keyword of $(nameof(owner))"
@noinline function _renamed_field_warning(owner, old::Symbol, new::Symbol)
    _hook_deprecation((owner, old),
        "$(_renamed_name(owner, old)) is deprecated and will be removed in SpaceAGORA 1.0; use $(new) instead: " *
        "only the name changes.")
    return new
end

# The default of a renamed keyword and of its old name: "not passed". It is not
# `nothing`, so an explicit `nothing` reaches the field.
struct _NotPassed end
const _NOT_PASSED = _NotPassed()

"""
    _renamed_keyword(owner, old, old_value, new, new_value, default)

The value of a keyword that was renamed from `old` to `new`, where `_NOT_PASSED`
means "not passed" (any other value, `nothing` included, counts as passed).
Passing `old` warns once per owner and name (see `_renamed_field_warning`);
passing both names throws an `ArgumentError`; passing neither gives `default`.
"""
function _renamed_keyword(owner, old::Symbol, old_value, new::Symbol, new_value, default)
    old_value === _NOT_PASSED && return new_value === _NOT_PASSED ? default : new_value
    new_value === _NOT_PASSED || throw(ArgumentError(
        "$(nameof(owner)): pass $(new) only; $(old) is its deprecated name."))
    _renamed_field_warning(owner, old, new)
    return old_value
end

end # module HookDeprecation
