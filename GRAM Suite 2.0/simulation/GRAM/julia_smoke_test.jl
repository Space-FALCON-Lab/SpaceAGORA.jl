function resolve_gram_root()
    gram_root = get(ENV, "GRAM_ROOT", "")
    if isempty(gram_root)
        candidates = [
            normpath(joinpath(@__DIR__, "..", "..")),
            normpath(joinpath(@__DIR__, "..", "..", "GRAM Suite 2.0")),
            normpath(joinpath(@__DIR__, "..", "..", "GRAM")),
        ]
        for c in candidates
            if isdir(joinpath(c, "Build")) && isdir(joinpath(c, "Julia"))
                gram_root = c
                break
            end
        end
    end

    isempty(gram_root) && error("Set GRAM_ROOT environment variable to your GRAM Suite root path.")
    return gram_root
end

gram_root = resolve_gram_root()
include(joinpath(gram_root, "Julia", "GRAM.jl"))
using .GRAM

function main(gram_root::AbstractString)
    libext = Sys.iswindows() ? "dll" : (Sys.isapple() ? "dylib" : "so")
    libpath = get(ENV, "GRAM_LIB", joinpath(gram_root, "Build", "lib", "libGRAM.$libext"))
    set_library!(libpath)
    initialize!(joinpath(gram_root, "SPICE"))

    atmos = create_atmosphere(BODY_MARS; data_path = joinpath(gram_root, "Mars", "data"))
    set_start_time!(atmos; year = 2020, month = 3, day = 15, hour = 0, minute = 0, seconds = 0.0, scale = 1, frame = 1)
    set_position!(atmos; height = 50.0, latitude = 22.0, longitude = 48.0, elapsed_time = 100.0)

    err = update!(atmos)
    if err != 0
        error("GRAM update failed: $(get_error_message())")
    end

    pos = get_position(atmos)
    dyn = get_dynamics_state(atmos)
    println("GRAM smoke test passed")
    println("height_km=$(pos.height)")
    println("temperature_K=$(dyn.temperature)")
    println("density_kg_m3=$(dyn.density)")

    close!(atmos)
end

main(gram_root)
