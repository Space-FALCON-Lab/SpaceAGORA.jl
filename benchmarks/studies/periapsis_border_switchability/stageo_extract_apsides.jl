# Stage O, WP-O1: the one short Arrow read (no pyarrow on the host).
# usage: julia --project=. stageo_extract_apsides.jl <dir with the two kernel feather files> <out dir>
# Writes True_Odyssey_{periapsis,apoapsis}_alts_kernel.csv (columns orbit, altitude; km).
# `stageo_reference.py build --dumps <out dir>` merges them into stageo_data/kernel_apsides.csv.
using Arrow, DataFrames
D, O = ARGS[1], ARGS[2]
for n in ("True_Odyssey_periapsis_alts_kernel", "True_Odyssey_apoapsis_alts_kernel")
    df = DataFrame(Arrow.Table(joinpath(D, n * ".feather")))
    open(joinpath(O, n * ".csv"), "w") do io
        println(io, join(names(df), ","))
        for r in eachrow(df); println(io, join(string.(collect(r)), ",")); end
    end
end
