include(joinpath(@__DIR__, "GRAM.jl"))
using .GRAM

libext = Sys.iswindows() ? "dll" : (Sys.isapple() ? "dylib" : "so")
set_library!(joinpath(@__DIR__, "..", "Build", "lib", "libGRAM.$libext"))
initialize!(joinpath(@__DIR__, "..", "SPICE"))

atmos = create_atmosphere(BODY_MARS; data_path = joinpath(@__DIR__, "..", "Mars", "data"))
set_start_time!(atmos; year = 2020, month = 3, day = 15, hour = 0, minute = 0, seconds = 0.0, scale = 1, frame = 1)
set_position!(atmos; height = 50.0, latitude = 22.0, longitude = 48.0, elapsed_time = 100.0)

err = update!(atmos)
if err != 0
    error("GRAM update failed: $(get_error_message())")
end

pos = get_position(atmos)
dyn = get_dynamics_state(atmos)

println("Height (km): ", pos.height)
println("Temperature (K): ", dyn.temperature)
println("Density (kg/m^3): ", dyn.density)

close!(atmos)
