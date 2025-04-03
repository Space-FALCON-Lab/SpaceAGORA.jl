# SPICE Data

This folder contains ephemeris data for the NAIF SPICE library utilized within the GRAM Suite.  The GRAM team considers this a SPICE kernels "starter pack" so that one can quickly set up and run the GRAM Suite.  **For analytical purposes, it is important for each user of the GRAM Suite to choose and use the appropriate SPICE kernels for their application from the NAIF website.** 

To that end, it is now possible to override the default starter pack kernels using NAMELIST file entries.  These can be specified in each NAMELIST file, or per folder in a file named "spice.txt".  All overrides are paths relative to the specified SpicePath.

+ SpicePath    = Path to NAIF SPICE data.
+ SpiceLsk     = Optional override of the SPICE leapseconds LSK file.
+ SpicePck     = Optional override of the SPICE planetary constants PCK file.
+ SpiceVenus   = Optional override of the SPICE Venus kernel.
+ SpiceEarth   = Optional override of the SPICE Earth kernel.
+ SpiceMars    = Optional override of the SPICE Mars kernel.
+ SpiceJupiter = Optional override of the SPICE Jupiter kernel.
+ SpiceSaturn  = Optional override of the SPICE Saturn kernel.
+ SpiceUranus  = Optional override of the SPICE Uranus kernel.
+ SpiceNeptune = Optional override of the SPICE Neptune kernel.
+ SpiceTitan   = Optional override of the SPICE Saturn kernel (used for Titan).

Note that each GRAM program requires the LSK file, the PCK file, and one planetary kernel.  If one is running VenusGRAM, say, then one only need override SpiceVenus.

The SPICE kernel starter pack provided with the GRAM Suite originated from the NAIF SPICE website.  The planetary kernels were reduced in size using the "spkmerge" tool from the cspice toolkit.  The data for moons was stripped from the original files and the time frames were reduced.  These kernels are restricted to the dates below.

+ BEGIN_TIME = 1 JAN 2000 00:00:00.000
+ END_TIME = 1 JAN 2100 00:00:00.000

If the range of these kernels does not meet your needs or if you need the latest ephemeris data, then you should visit the NAIF SPICE website and download the appropriate kernels for your project.  We also suggest that you remove the starter pack files so they are not accidentally utilized.

The NAIF SPICE kernels are available via File Transfer Protocol (FTP) from ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels. Information about the SPICE data is available from https://naif.jpl.nasa.gov/naif/data.html and help downloading is available from https://naif.jpl.nasa.gov/naif/download_tip.html.