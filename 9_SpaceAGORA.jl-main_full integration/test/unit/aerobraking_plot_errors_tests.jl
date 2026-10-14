module AerobrakingPlotErrorTests
using Test

const PROJECT_BEFORE_INCLUDE = Base.active_project()
const DIRECTORY_BEFORE_INCLUDE = pwd()
include(joinpath(@__DIR__, "..", "..", "examples", "support", "aerobraking_plot_errors.jl"))

@testset "standalone RTN plot-error preparation" begin
    @test Base.active_project() == PROJECT_BEFORE_INCLUDE
    @test pwd() == DIRECTORY_BEFORE_INCLUDE
    for binding in (:SpaceAGORA, :Plots, :PlotlyJS, :SPICE, :SimulationConfiguration,
                    :CSV, :DataFrames, :RuntimeServices)
        @test !isdefined(@__MODULE__, binding)
    end

    # The second reference frame has R=+y, T=+z and N=+x. Distinct,
    # signed errors expose swapped axes, row reordering and signed outputs.
    # A radial velocity component in row 1 must not tilt the transverse axis.
    samples = ([-2, -5], [3, -6], [-4, 7],
               [2, 0], [0, 3], [0, 0],
               [3, 0], [4, 0], [0, 2])
    original = deepcopy(samples)
    @test _rtn_error_components(samples...) == ([2.0, 6.0], [3.0, 7.0], [4.0, 5.0])
    @test samples == original
    @test _rtn_error_components((view(v, 2:2) for v in samples)...) ==
          ([6.0], [7.0], [5.0])

    # A 3-4-5 rotation gives R=(3/5,4/5,0), T=(-4/5,3/5,0).
    oblique = ([1.0], [-7.0], [-4.0],
               [3.0], [4.0], [0.0], [-4.0], [3.0], [0.0])
    r, t, n = _rtn_error_components(oblique...)
    @test r ≈ [5.0]
    @test t ≈ [5.0]
    @test n ≈ [4.0]
    # The same helper prepares position and velocity errors without converting
    # their units. Rescale errors alone while holding the reference frame fixed.
    scaled = (map(v -> v ./ 1000, oblique[1:3])..., oblique[4:9]...)
    r, t, n = _rtn_error_components(scaled...)
    @test r ≈ [0.005]
    @test t ≈ [0.005]
    @test n ≈ [0.004]

    floor = eps(Float64)
    @test _rtn_error_components([0.0], [-floor/2], [floor],
                               [1.0], [0.0], [0.0], [0.0], [1.0], [0.0]) ==
          ([floor], [floor], [floor])
    empty_result = _rtn_error_components((Float64[] for _ in 1:9)...)
    @test empty_result == (Float64[], Float64[], Float64[])
    @test all(v -> v isa Vector{Float64}, empty_result)
    @test length(unique(objectid.(empty_result))) == 3

    # Both zero position and parallel position/velocity fail on the offending
    # row. Keep the existing error text and sample index for the plot callers.
    for bad_ref in (([1.0, 0.0], [0.0, 0.0], [0.0, 0.0],
                     [0.0, 0.0], [1.0, 1.0], [0.0, 0.0]),
                    ([1.0, 1.0], [0.0, 0.0], [0.0, 0.0],
                     [0.0, 2.0], [1.0, 0.0], [0.0, 0.0]))
        err = try
            _rtn_error_components(ones(2), ones(2), ones(2), bad_ref...)
            nothing
        catch caught
            caught
        end
        @test err isa ArgumentError
        @test sprint(showerror, err) ==
              "ArgumentError: Cannot construct RTN frame at sample 2: degenerate reference position/angular momentum."
    end
end
end # module AerobrakingPlotErrorTests
