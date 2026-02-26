# GRAM Suite 2.0.0

Please read "CRITICAL INFORMATION AND ERRATA.pdf" located in the GRAM\Documentation folder.

If you'd like to run one of the pre-built GRAM programs, please refer to the User Guides found in GRAM\Documentation folder.

If you'd like to integrate the GRAM library into your program, please refer to the GRAM Programmer's Manual found in GRAM\Documentation folder.

For developers working with GCC or MinGW, check out the README file in the Build subfolder.

For developers working with MSVS 2017, check out the README file in the MSVS subfolder.

For Python developers, check out the README in the Python subfolder.

For Julia developers, check out the README in the Julia subfolder.

If you want to regenerate the latest version of the GRAM Programmer's Manual:
 1. Download and install doxygen (http://www.doxygen.nl/download.html)
 2. Run the doxywizard
 3. Open the file "Doxyfile" in this same folder.
 4. Click on the Run tab. Click on the "Run doxygen" button.
 5. Click on the "Show HTML Output" button.

Files and Subfolders:
- Documentation: User Guides, GRAM Programmer's Manual, a report detailing the updated methodology to compute the speed of sound in the GRAM Suite, the GRAM Suite change log, and subfolders of additional Mars-GRAM and Earth-GRAM documentation.
- Windows: Prebuilt Windows executables.
- Linux: Prebuilt Linux executables.
- Build: A makefile system for building the GRAM Suite.
- MSVS: A Visual Studio solution for building the GRAM suite (no FORTRAN).
- Python: A Python interface into the GRAM Suite.
- Julia: A minimal Julia wrapper around the GRAM C interface.
- common: A framework shared by all GRAM models.
- Planet folders: The model specific code, examples, and tests for each planet.
- GRAM: Source files for examples that combine all GRAM models.
- SPICE: Data kernels for ephemeris computations using SPICE.
- Doxyfile and DoxygenLayout.html: Configuration files used to generate the Programmer's Manual.

Please note that the GRAM Suite uses the NAIF SPICE libraries. The files provided with this bundle have been reduced in size for easy distribution.  If the full NAIF data is required, that data must be downloaded from the NAIF web site.  More information on this proceedure is available in the Documentation folder.
