using CSV
using Plots
using DataFrames
using LinearAlgebra
gr()
# plotly()

δ(i, j) = ==(i, j)
"""
    alf_IDR(x::Real, N::Integer, M::Integer)

Fully normalized Associate Legendre Functions (fnALFs) from degree 0 to N and order 0 to M,
evaluate at x ∈ [0, 1]. The ALFs are computed through Increasing Degree Recursion (IDR)
formulas that keep either the degree or the order fixed. All ALFs for M > N are zero by
definition.
"""
function fnALF_IDR(x, N, M)
    # Allocate the output array.
    # TODO: this should be done in-place to avoid allocating a new array at each run.
    local A = zeros(N+1,M+1)

    # Precompute sqrt
    ξ = √(1 - x^2)
    # Initialize top left element
    A[1,1] = 1;

    # Sectorials
    if M > 0
        for i=2:min(N+1,M+1)
            # For ease of notation
            n = i - 1

            # TODO: preallocate fn's for speed.
            fn = √( ((1 + δ(1,n))*(2n + 1)) / 2n )
            A[i,i] = fn * ξ * A[i-1,i-1]

        end
    end

    # Zonals and tesserals
    for j = 1:M+1
        for i = j+1:N+1
            # For ease of notation
            m = j - 1
            n = i - 1
            #TODO: preallocate gnm, hnm for speed
            gnm = √(((2n + 1)*(2n-1)) / ((n+m)*(n-m)))

            if i == j + 1
                A[i,j] = gnm * x * A[i-1,j]
            else
                hnm = √(((2n + 1)*(n-m-1)*(n+m-1))/((2n-3)*(n+m)*(n-m)))
                A[i,j] = gnm * x * A[i-1,j] - hnm * A[i-2,j]
            end
            
        end      
    end

    return A
end

data = CSV.read("/workspaces/ABTS.jl/Topography_harmonics_data/MOLA.csv", DataFrame)
data_true = CSV.read("/workspaces/ABTS.jl/Topography_harmonics_data/topogrd.csv", DataFrame, header=false)
println(size(data_true))
data_true = reshape(Matrix(data_true), (180, 360))
data_true = circshift(circshift(data_true, (120, 0)), (-180, 0))
data_true = reshape(data_true, 360*180)
lats_grid = range(-90, 90, length=180)' .* ones(360, 1)
lons_grid = range(-180, 180, length=360) .* ones(1, 180)
lats_grid = reshape(lats_grid, 360*180)
lons_grid = reshape(lons_grid, 360*180)
println(size(lats_grid))
println(size(lons_grid))
println(size(data_true))
degree = maximum(data[:, 1]) + 1

Clm = zeros(degree,degree)

Slm = zeros(degree,degree)

for i = 1:degree
    local l = data[i, 1] + 1
    local m = data[i, 2] + 1
    Clm[l, m] = data[i, 3]
    Slm[l, m] = data[i, 4]
end

# Compute the fully normalized ALFs

num_points = 360
lats = range(-π/2, π/2, length=Int64(round(num_points/2)))
lons = range(-π, π, length=num_points)
elevations = zeros(size(lats,1)*size(lons,1), 3)
count = 1
l = 20
m = 20
println("Starting to calculate elevations")
for lon in lons
    λ = lon
    if m == 0
        sinmλ = [0.0]
        cosmλ = [1.0]
    elseif m == 1
        sinmλ = [0.0; sin(λ)]
        cosmλ = [1.0; cos(λ)]
    elseif m > 1
        sinmλ = [0.0; sin(λ); zeros(m-1)]
        cosmλ = [1.0; cos(λ); zeros(m-1)]
    end

    for j = 3:m+1
        sinmλ[j] = 2cosmλ[2] * sinmλ[j-1] - sinmλ[j-2]
        cosmλ[j] = 2cosmλ[2] * cosmλ[j-1] - cosmλ[j-2]
    end
    for lat in lats
        global count
        x = sin(lat)
        # A = zeros(l+1,m+1)
        local A = fnALF_IDR(x, l, m)
        # println("Calculating elevation at lat: ", lat, " lon: ", lon)
        for i = 1:l+1
            for j = 1:m+1
                # println("i: ", i, " j: ", j, " Clm: ", Clm[i,j], " Slm: ", Slm[i,j])

                elevations[count,1] += A[i,j] * (Clm[i,j]*cosmλ[j] + Slm[i,j]*sinmλ[j])
            end
        end
        elevations[count,2] = lat
        elevations[count,3] = lon
        count += 1
    end
end
println("Calculated elevations")
mean_radius = 6051.848
p = plot(rad2deg.(elevations[:,3]), rad2deg.(elevations[:,2]), elevations[:,1], st = :surface)
# plot!(p, lons_grid, lats_grid, data_true, st = :surface)
println("Made plot")
display(p)
println("Displayed plot")