Building the GRAM Interfaces
============================

Quick Links
===========

- [Files Needed for Building GRAMs](@ref BuildFiles)
- [Build using Microsoft Visual Studio](@ref BuildMSVS)
- [Build using GNU Make on Linux](@ref BuildGNU)

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

The GRAM suite is a collection of static libraries and programs.  The libraries can be easily linked into any simulation, or the provided source can be compiled directly.  This section details the build process and compiler requirements.

Compiling the GRAM suite requires a compiler that supports the following language standards:
 - C++11 
 - C99 (for the C interface)
 - FORTRAN 2003 (for the FORTRAN interface)

 The GRAM suite has been successfully built with the following compilers:
 - Windows
     + Microsoft Visual Studio 2017 (without the FORTRAN interface)
	 + Microsoft Visual Studio 2015 (without the FORTRAN interface)
	 + MinGW-w64 7.1 (under MSYS2 2.8.1)
	 + MinGW-w64 8.2 (under MSYS2 2.11.2)
 - Linux
     + GCC 6.3, 7.3, 8.2

- - - - - - - - - - - 	
Files Needed for Building GRAMs   {#BuildFiles}
===============================

The source code for the GRAM suite is organized in a modular fashion.  The *common* folder contains the common framework on which all GRAM models are built.  The planetary models are in folders orgranized by the body name.  The *MSVS* folder contains the Microsoft Visual Studio projects.  And the *Build* folder contains a GNU make system.  The location of the *MSVS* and *Build* folders relative to the source code should not be modified.
The source folders are organized as follows:
 - Planetary models and common
     + source: C, C++ and FORTRAN source code.
     + include: C++ and C header files.
     + examples: main programs that utilize the library.
     + unittest: Unit tests using the gtest framework.
 - Common
     + cspice: Header files for the cspice library (version N0066).
	 + googletest: Header and source files for the gtest framework (version 1.8.0).

As noted in the table above, the GRAM suite depends on two third party libraries.  The [NAIF CSPICE toolkit](https://naif.jpl.nasa.gov/naif/toolkit_C.html) provides ephemeris values for the atmosphere models.  This library should be built using the instructions provided with the NAIF CSPICE toolkit.  Some pre-built cspice libraries are found in the GRAM/common/cspice/lib folder.

Unit testing in the GRAM suite utilizes the GoogleTest framework.  The build systems will compile this library as part of the build process.

So what files are actually needed to build a GRAM progam or integrate with a simulation?  As noted, the code is organized so that one can create libraries to be linked with a main program.  Alternatively, one can compile directly by simply pointing to select folders.


###The Common Library
 + Source Folders
	  - Compile all files in GRAM/common/source
		   + These are the base classes used by all GRAM models.
	  - Optionally, compile all files in GRAM/common/examples
		   + These are generic examples that may be called by any GRAM model.
		   + These files are not required when linking GRAM models to a simulation.
		   + The examples require the GoogleTest framework.
 + Include Folders
	  - Add the folder GRAM/common/include to the include path  (gcc: -I GRAM/common/include)
	  - Add the folder GRAM/common/cspice/include to the include path 
	  - Optionally, add the folder GRAM/common/examples to the include path  
	       + If the examples are included, then also include the GoogleTest framework.
	  - Optionally, add the folder GRAM/common/googletest to the include path  
	  - Optionally, add the folder GRAM/common/googletest/include to the include path  
		  
###The Planetary Model Library (assuming the name Body)
 + Source Folders
	  - Compile all cpp files in GRAM/Body/source
		   + These are the classes that define the Body GRAM model.
	  - Optionally, compile the f90 file in GRAM/Body/source
		   + This is the FORTRAN interface for the Body GRAM model.
		   + This file is not required if the main program is not FORTRAN.
 + Include Folders
	  - Add the folder GRAM/Body/include to the include path 
		  
###Multiple Planetary Models
 + Simply repeat the section above for each desired body.
	 
###The SPICE Library
 + Instructions for building the SPICE library are contained in the toolkit distribution. 
 + The cspice library is required.  The csupport library is not required.
	 
###The GoogleTest Library
 + Source Files
	  - Compile the file GRAM/common/googletest/src/gtest-all.cc
 + Include Folders
	  - Add the folder GRAM/common/googletest to the include path 
	  - Add the folder GRAM/common/googletest/include to the include path  

###The Main Program
 + Either link in or directly compile in all of the above sections.
 + Source Folders
	  - Optionally, compile all unit test files.
		   + Common unit tests in GRAM/common/unittest
		   + Body model unit tests in GRAM/Body/unittest
		   + Unit tests cannot be added to the libraries, they must be compiled with the main program.
		   + Unit tests are not required when integrating GRAM models with a simulation.
 + Include Folders
	  - Add the folder GRAM/Body/include to the include path  
	  - Add the folder GRAM/common/include to the include path 
	  - Add the folder GRAM/common/cspice/include to the include path  
	  - Optionally, add the folder GRAM/common/examples to the include path 
	  - Optionally, add the folder GRAM/common/googletest to the include path  
	  - Optionally, add the folder GRAM/common/googletest/include to the include path 
	 
###Summary (all together)
 + Source Folders
	  - GRAM/Body/source (for each body model)
	  - GRAM/common/source
	  - GRAM/common/examples (optional)
	  - GRAM/common/googletest/src/gtest-all.cc (optional)
	  - GRAM/common/unittest (optional)
	  - GRAM/Body/unittest (optional)
 + Include Folders
	  - GRAM/Body/include (for each body model)
	  - GRAM/common/include 
	  - GRAM/common/cspice/include
	  - GRAM/common/examples (optional)
	  - GRAM/common/googletest (optional)
	  - GRAM/common/googletest/include (optional)

- - - - - - - - - - - 	
Build using Microsoft Visual Studio   {#BuildMSVS}
===================================

Steps for building the GRAM suite using MSVS 2017:

 1. Open the solution GRAM\\MSVS\\GRAMs.sln.  
 2. Set the Solution Configuration to "Release".
 3. Set the Solution Platform to "x64".
 4. Select "Rebuild Solution" from the Build menu.
 
These steps will compile all GRAM libraries, programs, and examples.  The binaries can be found in the GRAM\\MSVS\\x64\\Release folder.  At this point, only the x64 platform is supported.

The projects are organized to build static libraries that are then combined to create the executables.  All executables require, at a minimum, the common library created in GramCommonLib and at least one body specific library, such as the one created in NeptuneLib.  Alternatively, an executable can link with the all inclusive GRAM suite library created in the GRAMLib project.

For programs that wish to run the GRAM unit tests, the googletest library should be compiled using the googletest project.  The code required to run the unit tests appears throughout the examples.

Most of the project settings are the default values.  Settings of particular interest are discussed in the table below.
 - C/C++ : General : Additional Include Directories
     + The include directories should correspond with the linked libraries.
     + Required: common\\include for the common library.
	 + Required: at least one body library, such as Neptune\\include. 
	 + Optional: common\\googletest\\include for the googletest library.
 - Linker : Input : Additional Dependencies
     + The GRAM suite option
	     - Required: GRAM.lib
		 - Optional: googletest.lib
     + The body specific option
	     - Required: GRAMCommon.lib
	     - Required: At least one body specific library, such as Neptune.lib. 
		 - Optional: Multiple body libraries may be used together.
		 - Optional: googletest.lib
 - C/C++ : Preprocesor : Preprocesor Definitions
	 + _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING  (moving to GTest 1.8.1 may alleviate this need)
 - Librarian : General : Additional Dependencies (and Additional Library Directories)
     + Includes the NAIF CSPICE library (cspice_x64.lib or cspice_x86.lib).
	 + This library should be built following the CSPICE Toolkit instructions.
	     - Start up the appropriate Visual Studio command prompt.
		 - Change to the cspice\\src\\cspice folder.
		 - Issue the command "mkprodct.bat".
		 - The library will be found in the cspice\\lib folder.
	 + Some pre-built binaries are provided in the GRAM\\common\\cspice\\lib folder.
 - Other settings
     + All other choices are non-critical and can be set to the preferences of the developer.

For developers running MSVS 2015, the same solution and projects can be used with some modifications.  Follow these steps:
 1. Start MSVS 2015 and open GRAM\\MSVS\\GRAMs.sln.
 2. In the Solution Explorer use ctrl-click to select all projects excluding the Shared Items projects.
 3. Right-click and select *Properties*.
 4. In the dropdowns select "All Configurations" and "All Platforms".
 5. Under *Configuration Properties : General* set the *Platform Toolset* to "Visual Studio 2015 (v140)".
 6. Set the *Target Platform Version* to "8.1".
 7. Apply the changes and save the solution.
 
- - - - - - - - - - - 	
Build using GNU Make on Linux (or MinGW)    {#BuildGNU}
========================================

Steps for building the GRAM suite using GCC and GNU Make:

 1. Set the build environment in makefile.defs.
 2. Open a command shell and set the directory to GRAM/Build.
 3. Enter the command “make clean”.
 4. Enter the command “make –j”.

These steps will compile all GRAM libraries, programs, and examples.  The executables can be found in the GRAM/Build/bin folder.  The libraries can be found in the GRAM/Build/lib folder.

The projects are organized to build static libraries that are then combined to create the executables.  All executables require, at a minimum, the common library created in Build/GramCommon and at least one body specific library, such as the one created in Build/Neptune.  Alternatively, an executable can link with the all inclusive GRAM suite library created in the Build/GRAM folder.

For programs that wish to run the GRAM unit tests, the googletest library should be compiled using the makefile in the Build/googletest folder.  The code required to run the unit tests appears throughout the examples.

The Build/makefile will default to building all targets with "make -j".  To build libraries and programs for a specific body, use the body name as the target as in "make Neptune -j".

The settings used in the makefiles can be modified in Build/makefile.defs.  The table below summarizes some of the key settings.  Other settings are detailed in the makefile.defs file.

 - BUILD
     + "release" is optimized for speed
	 + "debug" is not optimized and contains information for debuggers
 - CXX, CC, FF, LNK
     + The command that invokes the C++ compiler, C compiler, FORTRAN compiler, and the linker.
 - CXX_FLAGS
     + Must be set to use the C++11 standard.
 - C_FLAGS
     + Must be set to use the C99 standard.
 - F_FLAGS
     + Must be set to use the FORTRAN 2003 standard.
 - BINARY_EXTENSION
     + Appended to the names of the executables (typically blank on Linux and ".exe" under MinGW).
 - SPICE_LIB
     + Path to the NAIF CSPICE library. 
	 + Some pre-built binaries are provided in the GRAM/common/cspice/lib folder.
	 + This library should be built following the CSPICE Toolkit instructions.
		 - Change to the cspice/src/cspice folder.
		 - Issue the command "mkprodct.csh".
		 - The resulting cspice library will be found in the cspice/lib folder.
