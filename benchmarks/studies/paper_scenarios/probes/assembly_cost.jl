# Prices the results-table assembly: `_direct_assembly_columns!` against the
# generic comprehension path it bypasses, on synthetic save data of the two
# shapes the built-in catalogue produces (a 3-vector field and a scalar field).
#
# Assembly only -- no solve -- so the numbers are not diluted by the run around
# them. Prints one `ASM ...` line per (size, rows, field, implementation) with
# wall time and allocation.
#
#   julia --project=. <this file>
using SpaceAGORA, DataFrames, StaticArrays, Statistics
const IOO = SpaceAGORA.SimulationModel.IOOutputs
const SaveData = SpaceAGORA.SimulationModel.SaveData
struct F; name::Symbol; per_satellite::Bool; column_prefix::String; end

function generic!(df, field, saved_data, num_sats)
    field_series = [snapshot[field.name] for snapshot in saved_data]
    for sat_idx in 1:num_sats
        sat_series = [value[sat_idx] for value in field_series]
        IOO._append_series_columns!(df, "sc$(sat_idx)_$(field.column_prefix)", sat_series)
    end
end
fast!(df, field, saved_data, num_sats) = IOO._direct_assembly_columns!(df, field, saved_data, num_sats)

function bench(num_sats, n_rows)
    # One vector field and one scalar field, the two shapes the catalogue has.
    vec_data = [ (d=SaveData(); d[:position]=[SVector{3,Float64}(t+i,t-i,t*i) for i in 1:num_sats]; d) for t in 1:n_rows ]
    sca_data = [ (d=SaveData(); d[:altitude]=[Float64(t*i) for i in 1:num_sats]; d) for t in 1:n_rows ]
    for (label, data, field) in (("vector(3)", vec_data, F(:position,true,"pos")),
                                 ("scalar",    sca_data, F(:altitude,true,"alt")))
        for (impl, fn) in (("generic", generic!), ("direct", fast!))
            fn(DataFrame(), field, data, num_sats)   # warm
            GC.gc(); GC.gc()
            s = @timed fn(DataFrame(), field, data, num_sats)
            println("ASM sats=$num_sats rows=$n_rows field=$label impl=$impl " *
                    "s=$(round(s.time,digits=4)) alloc_mib=$(round(s.bytes/2^20,digits=2))")
        end
    end
end
for (n, r) in ((256, 121), (1024, 121), (4096, 121), (1024, 601))
    bench(n, r)
end
