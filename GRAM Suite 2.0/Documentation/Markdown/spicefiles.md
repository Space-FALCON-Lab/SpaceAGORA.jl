NAIF SPICE Files Needed to Run GRAMs
====================================

This document details the NAIF SPICE files that are required to run particular GRAM models.  These files have been be acquired via FTP from ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels and reduced in size by removing unused moons and restricting start and end dates (1/1/2000 to 1/1/2100).  Each of these files can be overridden with the original NAIF kernel or a kernel of your choosing.


| %GRAM            | File            | Path                 | Original File   |
|------------------|-----------------|----------------------|-----------------|
| All              | naif0012.tls    | spice/lsk            | naif0012.tls    |
| All              | pck0011.tcp     | spice/pck            | pck0011.tcp     |
| Venus            | de440_GRAM.bsp  | spice/spk/planets    | de440.bsp       |
| Earth            | de440_GRAM.bsp  | spice/spk/planets    | de440.bsp       |
| Mars             | mar097_GRAM.bsp | spice/spk/satellites | mar097.bsp      |
| Jupiter          | jup365_GRAM.bsp | spice/spk/satellites | jup365.bsp      |
| Saturn/Titan     | sat441_GRAM.bsp | spice/spk/satellites | sat441.bsp      |
| Uranus           | ura116_GRAM.bsp | spice/spk/satellites | ura116.bsp      |
| Neptune          | nep101_GRAM.bsp | spice/spk/satellites | nep101.bsp      |

--------------------------------------------------------------------------------------------

SPICE Initialization
====================

One can optionally override the default SPICE kernels.  This should be done before calling setSpiceDataPath.
~~~~~~~~~~
  SpiceLoader::setSpiceLsk("/lsk/naif0012.tls");
  SpiceLoader::setSpicePck("/pck/pck00011.tpc");
  SpiceLoader::setSpiceKernel(NEPTUNE, "/spk/satellites/nep101.bsp");
  SpiceLoader::setSpiceKernel(URANUS, "/spk/satellites/ura116.bsp");
~~~~~~~~~~
The SPICE library requires that the files shown above be loaded into memory before calls are made into the library.  Since the Ephemeris and GramTime classes depend on SPICE libraries, care must be taken to initialize SPICE before these classes are invoked.  In particular, the SPICE data path must be set so that the location of these files is known. The SPICE data path can be set using the SpiceLoader class.  The SpiceLoader class is responsible for loading SPICE data files into memory.
~~~~~~~~~~
  SpiceLoader spiceLoader;
  spiceLoader.setSpiceDataPath(pathToSpiceData);
~~~~~~~~~~
The SPICE data path can also be set using a BodyAtmosphere class.
~~~~~~~~~~
  MarsAtmosphere mars;
  mars.setSpiceDataPath(pathToSpiceData);
~~~~~~~~~~
Or the input parameters can be used to set the SPICE data path.
~~~~~~~~~~
  MarsAtmosphere mars;
  MarsInputParameters params;
  params.spicePath = pathToSpiceData;
  mars.setInputParameters(params);
~~~~~~~~~~
The NamelistReader class has a convenient method which looks for a "spice.txt" name list file.  If found, it will parse the SPICE data path from the name list entries.
~~~~~~~~~~
  InputParameters params;
  NamelistReader reader;
  reader.tryGetSpicePath(params);
  SpiceLoader spiceLoader;
  spiceLoader.setSpiceDataPath(params.spicePath);
~~~~~~~~~~
Individual kernels can be overridden using the InputParameters class.  Note that all files are assumed to reside in the spicePath.
~~~~~~~~~~
  MarsAtmosphere mars;
  MarsInputParameters params;
  params.spicePath = pathToSpiceData;
  params.spiceMars = "/spk/satellites/mar097.bsp"; // Use the full sized NAIF kernel.
  mars.setInputParameters(params);
~~~~~~~~~~
