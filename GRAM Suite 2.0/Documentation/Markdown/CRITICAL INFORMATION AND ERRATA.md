# Critical Information

This document points out significant modifications to the GRAM Suite documentation.

## Earth-GRAM MERRA-2 Data Download
**Multiple downloads are required to run Earth-GRAM with the MERRA-2 model.**  The GRAM Suite zip bundle contains programs, code, and some data required by the GRAM models.  Earth-GRAM also requires lower atmosphere data files for either the MERRA-2 model or the NCEP model.

The MERRA-2 data 2.0-degree files are divided into 9 separate data links: *MERRA2data_00Z.zip, MERRAdata_03Z.zip, MERR2data_06Z.zip, MERRA2data_09Z.zip, MERRA2data_12Z.zip, MERRA2data_15Z.zip, MERRA2data_18Z.zip, MERRA2data_21Z.zip*, and *MERRA2data_All_Mean.zip*.  The 9 files, at 2.0-degree resolution, are for each time of day (00Z, 03Z, 06Z, 09Z, 12Z, 15Z, 18Z, 21Z) and the file for all times averaged (All Mean). These zip files contain binary files for each month. 
   + Extract the data into the GRAM Suite folders GRAM/Earth/data/MERRA2data.  This is the default location.  Afterwards, the MERRA2data folder should contain nine folders (00Z, 03Z, 06Z, 09Z, 12Z, 15Z, 18Z, 21Z, All Mean), the MERRA2info.txt file, and a README.md markdown file.
   + Extract the data into a location of your choice, say /mypath/MERRA2data, and set MERRA2Path = '/mypath/MERRA2data' in the NAMELIST file.  Afterwards, the MERRA2data folder should contain nine folders, the MERRA2info.txt file, and a README.md markdown file.

## Earth-GRAM NCEP Data Download
**Two downloads are required to run Earth-GRAM with the NCEP model.**  The GRAM Suite zip bundle contains programs, code, and most data required by the GRAM models.  The EarthGRAM_NCEP_Data.zip bundle contains some large data files that are critical to the proper operation of the Earth-GRAM program.  The Earth-GRAM user has two options for integrating the data files.

+ Extract the data into the GRAM Suite folders GRAM/Earth/data/NCEPdata.  This is the default location.  Afterwards, the NCEPData folder should contain two folders, ACSCII and FixedBin, and a README.md markdown file.
+ Extract the data into a location of your choice, say /mypath/NCEPdata, and set NCEPPath = '/mypath/NCEPdata' in the NAMELIST file.  Afterwards, the NCEPData folder should contain two folders, ACSCII and FixedBin, and a README.md markdown file.

## SPICE Overrides
**For analytical purposes, it is important for each user of the GRAM Suite to choose and use the appropriate SPICE kernels for their application from the NAIF website.** 

To that end, it is now possible to override the default kernels using NAMELIST file entries.  These can be specified in each NAMELIST file, or per folder in a file named "spice.txt".  All overrides are paths relative to the specified SpicePath.

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

## SPICE Kernels Starter Pack
The GRAM Suite distribution now contains a folder of ephemeris data for the NAIF SPICE library utilized within the GRAM Suite.  The GRAM team considers this folder as a SPICE kernels "starter pack" so that one can quickly set up and run the GRAM Suite.

The SPICE kernel starter pack provided with the GRAM Suite originated from the NAIF SPICE website.  The planetary kernels were reduced in size using the "spkmerge" tool from the cspice toolkit.  The data for moons was stripped from the original files and the time frames were reduced.  These kernels are restricted to the dates below.

+ BEGIN_TIME = 1 JAN 2000 00:00:00.000
+ END_TIME = 1 JAN 2100 00:00:00.000

The new defaults for the SPICE kernels are as shown below.  The original NAIF SPICE file name is contained within each file name.

+ SpiceLsk     = /lsk/naif0012.tls
+ SpicePck     = /pck/pck00011.tpc
+ SpiceVenus   = /spk/planets/de440_GRAM.bsp
+ SpiceEarth   = /spk/planets/de440_GRAM.bsp
+ SpiceMars    = /spk/satellites/mar097_GRAM.bsp
+ SpiceJupiter = /spk/satellites/jup365_GRAM.bsp
+ SpiceSaturn  = /spk/satellites/sat441_GRAM.bsp
+ SpiceUranus  = /spk/satellites/ura116_GRAM.bsp
+ SpiceNeptune = /spk/satellites/nep101_GRAM.bsp
+ SpiceTitan   = /spk/satellites/sat441_GRAM.bsp

---------------------------------------------------------------------

# Venus-GRAM User Guide New Information and Errata

+ The Venus-GRAM 1 degree by 1 degree topography grid is derived from Ford and Pettengill (Journal of Geophysical Research, 97, E8, 13103-13114, 1992) Magellan topography data.  The data was acquired from the PDS Geosciences Node.
+ Appendix A: Two output fields have been added.
     + SurfaceHeight_km: Height of the surface relative to the reference ellipsoid.
     + HeightAboveSurface_km: Current height relative to the surface.

---------------------------------------------------------------------

# All User Guides New Information and Errata

+ Section 3.2 (Earth-GRAM Section 4.2): The SPICE files are now provided in a reduced form as part of the release bundle.
+ Appendix A: Corrections.
     + Density perturbation (kg/m3) *should be* Density perturbation (%)
     + DensityStandardDeviation_kgm3 *should be* DensityStandardDeviation_pct
     + Standard deviation of the density (kg/m3) *should be* Standard deviation of the density (%)
