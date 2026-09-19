"""
    with_density_model_epoch(model::AbstractDensityModel, initial_time)

Return a density model whose native elapsed-time origin matches `initial_time`.
The default implementation returns `model` unchanged; custom epoch-dependent
models may extend this hook. `run_simulation` calls it before state isolation.

For a keyword-built `GRAMAtmosphereModel`, an unchanged epoch returns the same
object. A changed epoch reconstructs a fresh native model from its owned
construction recipe, preserving the other recorded options. Equality uses the
stored integer calendar fields and `Float32` seconds, without a time tolerance.
A raw-core wrapper with an unknown recipe cannot be realigned safely and throws.

A `GRAMAtmosphereModelSurrogate` backed by GRAM rejects a changed epoch: rebuilding
its native fallback cannot validate or re-epoch its fixed table. Construct and
validate a table/model for the requested epoch explicitly. An unchanged epoch
preserves the surrogate object, file and fallback setting. Surrogates with a
non-GRAM base retain the default no-op behavior.

This hook changes the construction epoch only. It does not convert time systems,
validate atmospheric coordinate/datum conventions, certify cached or surrogate
data, or preserve a previously advanced native random stream or manual handle
mutations. Existing native cache and environment policies still apply.
"""
with_density_model_epoch(model::AbstractDensityModel, initial_time) = model

function _density_epoch_key(epoch)
    fields = (:year, :month, :day, :hour, :minute, :second)
    all(name -> hasproperty(epoch, name), fields) ||
        throw(ArgumentError("A density model epoch must provide year, month, day, hour, minute and second."))
    values = try
        (Int32(epoch.year), Int16(epoch.month), Int16(epoch.day),
         Int16(epoch.hour), Int16(epoch.minute), Float32(epoch.second))
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("Density model epoch fields must be representable in InitialTime: $(sprint(showerror, err))"))
    end
    isfinite(values[6]) || throw(ArgumentError("Density model epoch seconds must be finite."))
    return values
end

function _gram_epoch_constructor_kwargs(model::GRAMAtmosphereModel, initial_time)
    _density_epoch_key(model.core.initial_time) == _density_epoch_key(initial_time) && return nothing
    model.constructor_kwargs === nothing && throw(ArgumentError(
        "Cannot realign a raw-core GRAMAtmosphereModel with unknown construction settings. " *
        "Construct a keyword-built GRAMAtmosphereModel at the requested initial_time explicitly."))
    recipe = deepcopy(model.constructor_kwargs)
    recipe[:initial_time] = deepcopy(initial_time)
    return recipe
end

# The extension provides native construction under the existing setup lock.
_rebuild_gram_epoch_model(recipe) = _gram_not_loaded_error("with_density_model_epoch")

function with_density_model_epoch(model::GRAMAtmosphereModel, initial_time)
    recipe = _gram_epoch_constructor_kwargs(model, initial_time)
    return recipe === nothing ? model : _rebuild_gram_epoch_model(recipe)
end

function with_density_model_epoch(model::GRAMAtmosphereModelSurrogate, initial_time)
    model.base_model isa GRAMAtmosphereModel || return model
    _density_epoch_key(model.base_model.core.initial_time) == _density_epoch_key(initial_time) && return model
    throw(ArgumentError(
        "Cannot change the epoch of a fixed GRAMAtmosphereModelSurrogate table. " *
        "Construct and validate a matching surrogate table/model for the requested epoch explicitly; " *
        "rebuilding the native fallback does not re-epoch the table."))
end
