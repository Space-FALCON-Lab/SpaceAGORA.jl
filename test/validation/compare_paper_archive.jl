# Optional regression against a fresh, untouched ZIP run; not a runtime dependency.
using CSV, DataFrames, Test
length(ARGS)==1 || error("Pass the original ZIP run's exact/ output directory.")
original=ARGS[1]
integrated=normpath(joinpath(@__DIR__,"..","..","outputs","paper_constellation_interlinks"))
rows=NamedTuple[]
@testset "Untouched paper ZIP versus SpaceAGORA integration" begin
    for objective in ("a","i","omega","combined")
        folder="maximize_$objective"
        old=CSV.read(joinpath(original,folder,"history.csv"),DataFrame)
        new=CSV.read(joinpath(integrated,folder,"history.csv"),DataFrame)
        os=CSV.read(joinpath(original,folder,"schedule.csv"),DataFrame)
        ns=CSV.read(joinpath(integrated,folder,"schedule.csv"),DataFrame)
        @test old.time_s ≈ new.time_s atol=1e-8 rtol=0
        differences=sum(Matrix(os[:,2:end]).!=Matrix(ns[:,3:end]))
        @test differences==0
        # SI versus kilometre KS scaling changes adaptive integration and dense
        # interpolation. Require identical schedules and tighter linear parity;
        # nonlinear history limits are 1 cm, 1e-9 in e and 1e-6 degrees.
        for (element,suffix,factor,linear_tolerance,nonlinear_tolerance) in
                (("a","a_m",1000.,1e-6,1e-2),("e","e",1.,1e-12,1e-9),
                 ("i","i_deg",1.,1e-10,1e-6),("raan","RAAN_deg",1.,1e-10,1e-6),
                 ("omega","omega_deg",1.,1e-10,1e-6))
            prediction_error=maximum(abs.(factor*old[!,"predicted_$element"]-new[!,"predicted_$suffix"]))
            nonlinear_error=maximum(abs.(factor*old[!,"actual_$element"]-new[!,"actual_$suffix"]))
            @test prediction_error<linear_tolerance
            @test nonlinear_error<nonlinear_tolerance
            push!(rows,(;objective,element=suffix,changed_schedule_cells=differences,
                max_linear_history_difference=prediction_error,max_nonlinear_history_difference=nonlinear_error,
                linear_terminal=new[end,"predicted_$suffix"],nonlinear_terminal=new[end,"actual_$suffix"]))
        end
    end
end
CSV.write(joinpath(integrated,"archive_regression.csv"),DataFrame(rows))
show(stdout,DataFrame(rows);allrows=true,allcols=true)
println()
