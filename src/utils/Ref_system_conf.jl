module ref_sys
    export OE, cartesian, R_RA_DEC, H_LAN_LON, uDuNuE, clock

    mutable struct OE
        a::Float64
        e::Float64
        i::Float64
        Ω::Float64
        ω::Float64
        vi::Float64
        m::Float64
    end

    mutable struct cartesian
        x::Float64
        y::Float64
        z::Float64
    end

    mutable struct R_RA_DEC
        r::Float64
        RA::Float64
        dec::Float64
    end

    mutable struct H_LAN_LON
        h::Float64
        LAT::Float64
        LON::Float64
    end

    mutable struct uDuNuE
        uD::Float64
        uN::Float64
        uE::Float64
    end

    mutable struct clock
        year::Int64
        month::Int64
        day::Int64
        hour::Int64
        minute::Int64
        second::Float64
    end

end