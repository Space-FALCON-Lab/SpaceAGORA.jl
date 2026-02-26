# The GRAM Suite Change Log

This document tracks significant modifications to the GRAM Suite.
+ GRAM suite versions are numbered using a major.minor.patch scheme.
  Version 2 marks the model upgrade phase of the GRAM Suite.
+ The framework versions are numbered by the year of the update.  If 
  multiple updates are made in the same year, then a letter is appended.
+ The model versions are numbered by the year of each major release. If
  minor modifications are made to the model to support the framework, 
  then a letter is appended.

---------------------------------------------------------------------

## Version 2.0 (October 2023)

+ Common framework upgraded to GRAMlib 2023b.
    + Bug fix: Missing units added to output header.
+ Earth-GRAM upgraded to 2023.
    + The MERRA-2 model has been added for the lower atmosphere. The NCEP model has been retained as an option.
    + Bug fix to NCEP calculation of vertical winds.
    + Bug fix to NCEP calculation of initial standard deviations.
    + Bug fix to EarthCorrelator.
    + Updated User Guide.
+ Adding Python interface.
    + Wraps the C++ interface.
    + Working on Windows and Linux OS.
    + Updated Programmer's Manual.
    
---------------------------------------------------------------------

## Version 1.5.1 (June 2023)

+ Common framework upgraded to GRAMlib 2023a.
    + SPICE file overrides have been added to the C and Fortran interfaces.
    + Bug fix: Limit on random number seeds increased to 2^24.
+ Venus-GRAM upgraded to 2023a.
    + SPICE file overrides have been added to the C and Fortran interfaces.
+ Earth-GRAM upgraded to 2021c.
    + SPICE file overrides have been added to the C and Fortran interfaces.
    + Bug fix to correctly override NCEPPath.
    + Bug fix to make sure data file reader handles both Linux and PC text file formats.
+ Mars-GRAM upgraded to 2021b.
    + SPICE file overrides have been added to the C and Fortran interfaces.
+ Jupiter-GRAM upgraded to 2021a.
    + SPICE file overrides have been added to the C and Fortran interfaces.
+ Titan-GRAM upgraded to 2020c.
    + SPICE file overrides have been added to the C and Fortran interfaces.
+ Uranus-GRAM upgraded to 2021b.
    + SPICE file overrides have been added to the C and Fortran interfaces.
+ Neptune-GRAM upgraded to 2019d.
    + SPICE file overrides have been added to the C and Fortran interfaces.
    
---------------------------------------------------------------------

## Version 1.5 (April 2023)

+ Common framework upgraded to GRAMlib 2023.
    + SPICE file overrides have been added (details below).
    + A SPICE kernel "starter pack" has been added to the release bundle (details below).
    + The cspice toolkit has been upgraded to version N0067.
    + The GoogleTest unit testing framework has been upgraded to version 1.12.1.
    + A new Topography class has been added (currently used by Venus-GRAM).
    + Improvements have been made to the makefile system.
+ Venus-GRAM upgraded to 2023.
    + Added topography model (details below).
+ Earth-GRAM upgraded to 2021b.
    + Bug fix to low/high density outputs.
    + Bug fix to memory management.
    + Speed improvements.  dayofyear, surface computations
    + Bug fix to error condition.
    + For distribution, the NCEP data files are now in a separate zip bundle.

### SPICE File Overrides
It is now possible to override the default SPICE kernels using NAMELIST file entries.  These can be specified in each NAMELIST file, or per folder in a file named "spice.txt".  All overrides are paths relative to the specified SpicePath.

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

### SPICE Kernels Starter Pack
**Please read "CRITICAL INFORMATION AND ERRATA.pdf" for detailed instructions.**  A SPICE folder has been added to the release bundle containing SPICE kernels needed to run the GRAM Suite.  This SPICE kernel starter pack originated from the NAIF SPICE website.  The planetary kernels were reduced in size using the "spkmerge" tool from the cspice toolkit.  The data for moons was stripped from the original files and the time frames were reduced.  The starter pack kernels are restricted to the dates below.

+ BEGIN_TIME = 1 JAN 2000 00:00:00.000
+ END_TIME = 1 JAN 2100 00:00:00.000

### Venus Topography
The Venus-GRAM 1 degree by 1 degree topography grid is derived from Ford and Pettengill (Journal of Geophysical Research, 97, E8, 13103-13114, 1992) Magellan topography data.  The data was acquired from the PDS Geosciences Node.  Venus-GRAM now outputs two new fields, SurfaceHeight and HeightAboveSurface.

---------------------------------------------------------------------

## Version 1.4.1 (July 2022)

+ Common framework upgraded to GRAMlib 2021c.
    + Bug fix to support vertical wind perturbation scales.
+ Mars-GRAM upgraded to 2021a.
    + Added setHeightAboveSurface_M to the C and Fortran interfaces.
+ Earth-GRAM upgraded to 2021a.
    + Bug fix to implementation of patchy turbulence.
    + Bug fix to the C and Fortran interfaces. The setRRAParameters 
      did not handle 2013 and 2019 cases.
    + Bug fix to RRA path initialization.
    + Bug fixes to EarthCorrMonte.
+ Bug fix to error message routine in all GRAM C interfaces.
 
---------------------------------------------------------------------

## Version 1.4 (October 2021)

+ Mars-GRAM 2021 added to the suite.
+ Earth-GRAM 2021 added to the suite.

---------------------------------------------------------------------

## Version 1.3 (September 2021)

+ Common framework upgraded to GRAMlib 2021b.
    + All GRAM models use the upgraded framework.
    + Heat capacity ratio uses improved formulation (details below).
    + Performance improvements in output.
    + Bug fixes.
+ Neptune-GRAM upgraded to 2019c.
    + Modified to support framework upgrades.
+ Titan-GRAM upgraded to 2020b.
    + Modified to support framework upgrades.
+ Uranus-GRAM upgraded to 2021a.
    + Modified to support framework upgrades.
+ Jupiter-GRAM 2021 added to the suite.
+ Venus-GRAM 2021 added to the suite.

### Heat Capacity Ratio
Prior to version 1.3, all GRAM models used a formulation of heat capacity 
ratio that was based on temperature and constituent gases using Shomate's 
equation.  This method is replaced by a formulation which takes into 
account both temperature and pressure to compute the Cv and Cp heat 
capacities of each constituent gas. These values are summed and used to 
compute the heat capacity ratio. 

---------------------------------------------------------------------

## Version 1.2.1 (July 2021)

+ Bug fixed: Name list reader was not handling file paths with spaces.
+ Bug fixed: Distance routine for auxiliary profiles used flat distances.
+ Bug fixed: Generic C interface for auxiliary profiles was not working.

---------------------------------------------------------------------

## Version 1.2 (June 2021)

+ Common framework upgraded to GRAMlib 2021.
    + All GRAM models use the upgraded framework.
    + The random number generator has been improved (details below).
    + Bug fixes and performance improvements.
+ Neptune-GRAM upgraded to 2019b.
    + Minor improvements and bug fixes.
    + Modified to support framework upgrades.
+ Titan-GRAM upgraded to 2020a.
    + Minor improvements and bug fixes.
    + Modified to support framework upgrades.
+ Uranus-GRAM 2021 added to the suite.

### Random Number Generator
Prior to version 1.2, many GRAM models used a random number generator 
with only 30,000 seeds and a potentially short period.  The framework 
has adopted the random number generator used in EarthGRAM.  This RNG
has 2^24^ = 16,777,216 seeds and a period of about 2^570^ = 10^171^.
More details can be found in the GRAM Programmer's Manual within the
RandomNumberGenerator class.

---------------------------------------------------------------------

## Version 1.1 (September 2020)

+ Common framework upgraded to GRAMlib 2020.
    + All GRAM models use the upgraded framework.
    + Heat capacity ratio uses improved formulation (details below).
    + Upgrades to the framework in order to support new models.
    + Improved documentation.
+ Neptune-GRAM upgraded to 2019a.
    + Improved documentation.
+ Titan-GRAM 2020 added to the suite.

### Heat Capacity Ratio
Prior to version 1.1, all GRAM models used a constant heat capacity 
ratio for the speed of sound computation. The heat capacity ratio is 
now based on temperature and constituent gases.  Shomate's equation is
used to compute the specific heat capacity with constant pressure for 
each constituent gas. 

---------------------------------------------------------------------

## Version 1.0 (May 2020)

+ Common framework, GRAMlib 2019, added to the suite.
+ Neptune-GRAM 2019 added to the suite.

