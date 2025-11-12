using LoopVectorization
δ(i, j) = ==(i, j)
"""
    alf_IDR(x::Real, N::Integer, M::Integer)

Fully normalized Associate Legendre Functions (fnALFs) from degree 0 to N and order 0 to M,
evaluate at x ∈ [0, 1]. The ALFs are computed through Increasing Degree Recursion (IDR)
formulas that keep either the degree or the order fixed. All ALFs for M > N are zero by
definition.
"""
function fnALF_IDR!(A, x, N, M)
    # Allocate the output array.
    # TODO: this should be done in-place to avoid allocating a new array at each run.
    # A = zeros(N+1,M+1)

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
        @inbounds for i = j+1:N+1
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

function calculate_topography_harmonics!(args, Clm, Slm, latitude, longitude, A)
    """
        Calculate the topography harmonics for a given planet. This function is used with the
        planet-specific conversion functions because all the planetary data gives the output in
        different formats. (ノಠ益ಠ)ノ彡┻━┻

        Parameters
        ----------
        args : Dict
            Dictionary with the ABTS input arguments.
        Clm : Array{Float64}
            Array with the cosine harmonics.
        Slm : Array{Float64}
            Array with the sine harmonics.
        latitude : Float64
            Latitude in radians.
        longitude : Float64
            Longitude in radians.
        
        
        Returns
        -------
        harmonic : Float64
            Value of the topography harmonic at the given latitude and longitude.
    """
    # Precompute the fnALFs
    l = args[:topo_degree] + 1
    m = args[:topo_order] + 1
    Plm = fnALF_IDR!(A, sin(latitude), l, m)

    # Initialize the harmonic
    harmonic = 0.0
    λ = longitude
    # Precompute sin/cos
    if m == 0
        sinmλ = [0.0]
        cosmλ = [1.0]
    elseif m == 1
        sinmλ = [0.0; sin(λ)]
        cosmλ = [1.0; cos(λ)]
        constant_cosmλ = 2cosmλ[2]
    elseif m > 1
        sinmλ = [0.0; sin(λ); zeros(m-1)]
        cosmλ = [1.0; cos(λ); zeros(m-1)]
        constant_cosmλ = 2cosmλ[2]
    end

    @inbounds for j = 3:m+1
        sinmλ[j] = constant_cosmλ * sinmλ[j-1] - sinmλ[j-2]
        cosmλ[j] = constant_cosmλ * cosmλ[j-1] - cosmλ[j-2]
    end
    harmonic_list = zeros(l*m)
    @inbounds @inline @turbo for i = 1:l
        for j = 1:m
            harmonic += Plm[i, j] * (Clm[i,j] * cosmλ[j] + Slm[i,j] * sinmλ[j])
        end
    end

    return harmonic
end


function Mars_elevation!(args, Clm, Slm, latitude, longitude, A)
    """
        Calculate the topography harmonics for Mars. This function is required because
        Mars data is given as the distance from the center of mass of Mars, and the Venus data is given as a 
        multiple of the planet's mean radius. (ノಠ益ಠ)ノ彡┻━┻

        Parameters
        ----------
        args : Dict
            Dictionary with the ABTS input arguments.
        Clm : Array{Float64}
            Array with the cosine harmonics.
        Slm : Array{Float64}
            Array with the sine harmonics.
        latitude : Float64
            Latitude in radians.
        longitude : Float64
            Longitude in radians.
        
        
        Returns
        -------
        elevation : Float64
            Radius of Mars at the given latitude and longitude.
    """

    harmonic = calculate_topography_harmonics!(args, Clm, Slm, latitude, longitude, A)
    # mean_radius = 3389.5e3
    elevation = harmonic
    return elevation
end

function Venus_elevation!(args, Clm, Slm, latitude, longitude, A)
    """
        Calculate the topography harmonics for Venus. This function is required because
        Mars data is given as the distance from the center of mass of Mars, and the Venus data is given as a 
        multiple of the planet's mean radius. (ノಠ益ಠ)ノ彡┻━┻

        Parameters
        ----------
        args : Dict
            Dictionary with the ABTS input arguments.
        Clm : Array{Float64}
            Array with the cosine harmonics.
        Slm : Array{Float64}
            Array with the sine harmonics.
        latitude : Float64
            Latitude in radians.
        longitude : Float64
            Longitude in radians.
        
        
        Returns
        -------
        elevation : Float64
            Radius of Venus at the given latitude and longitude.
    """

    harmonic = calculate_topography_harmonics!(args, Clm, Slm, latitude, longitude, A)
    mean_radius = 6051.8e3
    elevation = mean_radius + harmonic
    return elevation
end

function Earth_elevation!(args, Clm, Slm, latitude, longitude, A)
    """
        Calculate the topography harmonics for Earth. This function is required because
        Earth data is given as the distance from the center of mass of Earth, and the Venus data is given as a 
        multiple of the planet's mean radius. (ノಠ益ಠ)ノ彡┻━┻

        Parameters
        ----------
        args : Dict
            Dictionary with the ABTS input arguments.
        Clm : Array{Float64}
            Array with the cosine harmonics.
        Slm : Array{Float64}
            Array with the sine harmonics.
        latitude : Float64
            Latitude in radians.
        longitude : Float64
            Longitude in radians.
        
        
        Returns
        -------
        elevation : Float64
            Radius of Earth at the given latitude and longitude.
    """
    harmonic = calculate_topography_harmonics!(args, Clm, Slm, latitude, longitude, A)
    elevation = harmonic
    return elevation
end