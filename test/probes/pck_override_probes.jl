# Real CSPICE pool/rotation witnesses with tiny synthetic kernels, no ephemeris
# datasets, native GRAM, trajectory propagation or shared-process kernel changes.
using Test, SPICE, SpaceAGORA
const P = SpaceAGORA.SimulationModel.Planets

function text_kernel(root, relative, body)
    path = joinpath(root, relative)
    mkpath(dirname(path))
    write(path, "KPL/PCK\n\\begindata\n" * body * "\n\\begintext\n")
    path
end

function binary_kernel(root, relative, angle)
    path = joinpath(root, relative)
    mkpath(dirname(path))
    handle = SPICE.pckopn(path, "Synthetic override precedence fixture", 0)
    try
        for id in (3000, 31006)
            # One degree-zero record with three Euler coefficients. Call the
            # CSPICE writer directly: SPICE.jl's pckw02 wrapper derives degree
            # from the whole coefficient row, not its three angle components.
            coefficients = Float64[angle, 0.5, 0.25]
            ccall((:pckw02_c, SPICE.libcspice), Cvoid,
                (Cint, Cint, Cstring, Cdouble, Cdouble, Cstring, Cdouble,
                 Cint, Cint, Ptr{Cdouble}, Cdouble),
                handle, id, "J2000", 0.0, 1.0, "test segment", 1.0,
                1, 0, coefficients, 0.0)
            SPICE.handleerror()
        end
    finally
        SPICE.pckcls(handle)
    end
    path
end

function clear_pool()
    SPICE.kclear()
    P._reset_furnished_kernels!()
end

@testset "PCK override precedence and cache policy" begin
    mktempdir() do root
        constants = join(["BODY$(id)_RADII = ( 10 10 9 )\nBODY$(id)_GM = ( 100 )\n" *
                          "BODY$(id)_POLE_RA = ( 0 0 0 )\nBODY$(id)_POLE_DEC = ( 90 0 0 )\n" *
                          "BODY$(id)_PM = ( 0 1 0 )" for id in (399, 301, 499, 299, 606)], "\n")
        text_kernel(root, "pck/pck00011.tpc", constants * "\nXVAL_TEST_PRIORITY = ( 1 )\nBODY4_NUT_PREC_ANGLES = ( " * join(fill("0", 52), " ") * " )")
        for relative in ("lsk/naif0012.tls", "spk/planets/de430.bsp",
                         "spk/satellites/mar097.bsp", "spk/satellites/sat441.bsp")
            text_kernel(root, relative, "XVAL_DUMMY = ( 1 )")
        end
        text_kernel(root, "pck/gm_de440.tpc", "BODY301_GM = ( 101 )")
        fk = text_kernel(root, "tf/SPICELunaFrameKernel.tf", """
            FRAME_MOON_PA_DE421 = 31006
            FRAME_31006_NAME = 'MOON_PA_DE421'
            FRAME_31006_CLASS = 2
            FRAME_31006_CLASS_ID = 31006
            FRAME_31006_CENTER = 301
            XVAL_TEST_PRIORITY = ( 2 )
            """)
        binary_kernel(root, "pck/earth_latest_high_prec.bpc", 0.1)
        binary_kernel(root, "spk/satellites/SPICELunaCurrentKernel.bpc", 0.2)
        override = binary_kernel(root, "override.bpc", 0.7)
        textoverride = text_kernel(root, "override.tpc", "XVAL_TEST_PRIORITY = ( 9 )\nBODY301_GM = ( 999 )")
        clear_pool()
        SPICE.furnsh(fk)
        SPICE.furnsh(override)
        expected_earth = SPICE.pxform("J2000", "ITRF93", 0.5)
        expected_moon = SPICE.pxform("J2000", "MOON_PA_DE421", 0.5)
        try
            for order in ((P.Earth, P.Moon), (P.Moon, P.Earth))
                clear_pool()
                # Relative and absolute paths are both part of the supported API.
                withenv("SPACEAGORA_SPICE_PCK_OVERRIDES" => "override.bpc,$textoverride",
                        "SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH" => nothing) do
                    first_body = order[1]("", root)
                    order[2]("", root)
                    @test SPICE.pxform("J2000", "ITRF93", 0.5) == expected_earth
                    @test SPICE.pxform("J2000", "MOON_PA_DE421", 0.5) == expected_moon
                    @test SPICE.gdpool("XVAL_TEST_PRIORITY"; start=1, room=1) == [9.0]
                    @test P.Moon("", root).μ == 999e9
                    @test order[1]("", root) === first_body
                    count = SPICE.ktotal("ALL")
                    for _ in 1:20
                        P.Earth("", root); P.Moon("", root)
                    end
                    @test SPICE.ktotal("ALL") == count
                    # Loading other standard kernels cannot erase the override.
                    for body in (P.Mars, P.Venus, P.Titan)
                        body("", root)
                        @test SPICE.pxform("J2000", "ITRF93", 0.5) == expected_earth
                        @test SPICE.gdpool("XVAL_TEST_PRIORITY"; start=1, room=1) == [9.0]
                    end
                    withenv("SPACEAGORA_SPICE_PCK_OVERRIDES" => nothing) do
                        @test_throws ArgumentError P.Earth("", root)
                    end
                    withenv("SPACEAGORA_SPICE_PCK_OVERRIDES" => textoverride) do
                        @test_throws ArgumentError P.Moon("", root)
                    end
                end
            end
            clear_pool()
            withenv("SPACEAGORA_SPICE_PCK_OVERRIDES" => "missing.tpc") do
                @test_throws ArgumentError P.Earth("", root)
                @test SPICE.ktotal("ALL") == 0
            end
            withenv("SPACEAGORA_SPICE_PCK_OVERRIDES" => nothing,
                    "SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH" => nothing) do
                P.Earth("", root)
                @test SPICE.pxform("J2000", "ITRF93", 0.5) != expected_earth
                count = SPICE.ktotal("ALL")
                P.Earth("", root)
                @test SPICE.ktotal("ALL") == count
                withenv("SPACEAGORA_SPICE_PCK_OVERRIDES" => override) do
                    @test_throws ArgumentError P.Earth("", root)
                end
            end
            clear_pool()
            withenv("SPACEAGORA_SPICE_PCK_OVERRIDES" => override,
                    "SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH" => nothing) do
                P.Earth("", root)
                @test SPICE.pxform("J2000", "ITRF93", 0.5) == expected_earth
            end
            # The shipped Mars text override retains the same rotation as a
            # direct CSPICE load, including after other bodies are constructed.
            mars_override = normpath(joinpath(@__DIR__, "..", "..", "data", "spice_overrides", "mars_iau2009_pole.tpc"))
            clear_pool()
            SPICE.furnsh(joinpath(root, "pck/pck00011.tpc"))
            SPICE.furnsh(mars_override)
            epochs = (0.0, 1000.0)
            mars_expected = [SPICE.pxform("J2000", "IAU_MARS", t) for t in epochs]
            clear_pool()
            withenv("SPACEAGORA_SPICE_PCK_OVERRIDES" => mars_override,
                    "SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH" => nothing) do
                P.Mars("", root)
                P.Moon("", root)
                P.Earth("", root)
                for (t, expected) in zip(epochs, mars_expected)
                    @test SPICE.pxform("J2000", "IAU_MARS", t) == expected
                end
                @test SPICE.gdpool("BODY499_POLE_RA") == [317.68143, -0.1061, 0.0]
                @test all(iszero, SPICE.gdpool("BODY499_NUT_PREC_PM"))
            end
        finally
            clear_pool()
        end
    end
end
