function interp_linear(times::AbstractVector{<:Real}, values::AbstractVector{<:Real}, t::Float64)
    idx = searchsortedlast(times, t)
    if idx <= 0
        return Float64(values[1])
    elseif idx >= length(times)
        return Float64(values[end])
    end
    t1 = Float64(times[idx])
    t2 = Float64(times[idx + 1])
    y1 = Float64(values[idx])
    y2 = Float64(values[idx + 1])
    α = (t - t1) / (t2 - t1)
    return (1.0 - α) * y1 + α * y2
end

function run_case(args::SimulationConfiguration; isolate_state::Bool=true)
    return mktempdir() do tmp
        cd(tmp) do
            results_csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
            run_simulation(args; isolate_state=isolate_state)
            @test isfile(results_csv_path)
            return CSV.read(results_csv_path, DataFrame)
        end
    end
end

function run_case_silent(args::SimulationConfiguration; isolate_state::Bool=true)
    return mktempdir() do tmp
        cd(tmp) do
            results_csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
            redirect_stdout(devnull) do
                run_simulation(args; isolate_state=isolate_state)
            end
            @test isfile(results_csv_path)
            return CSV.read(results_csv_path, DataFrame)
        end
    end
end

function run_case_capture_stdout(args::SimulationConfiguration; expect_results_csv::Bool=true, isolate_state::Bool=true)
    return mktempdir() do tmp
        cd(tmp) do
            results_csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
            output = ""
            mktemp() do path, io
                redirect_stdout(io) do
                    run_simulation(args; isolate_state=isolate_state)
                end
                flush(io)
                seekstart(io)
                output = read(io, String)
            end
            if expect_results_csv
                @test isfile(results_csv_path)
                return CSV.read(results_csv_path, DataFrame), output
            else
                @test !isfile(results_csv_path)
                return DataFrame(), output
            end
        end
    end
end
