# EarthGRAM

**Multiple downloads are required to run Earth-GRAM.**  The GRAM Suite zip bundle contains programs, code, and some data required by the GRAM models.  Earth-GRAM also requires lower atmosphere data files for either the MERRA-2 model or the NCEP model.

The MERRA-2 data 2.0-degree files are divided into 9 separate data links: *MERRA2data_00Z.zip, MERRAdata_03Z.zip, MERR2data_06Z.zip, MERRA2data_09Z.zip, MERRA2data_12Z.zip, MERRA2data_15Z.zip, MERRA2data_18Z.zip, MERRA2data_21Z.zip*, and *MERRA2data_All_Mean.zip*.  The 9 files, at 2.0-degree resolution, are for each time of day (00Z, 03Z, 06Z, 09Z, 12Z, 15Z, 18Z, 21Z) and the file for all times averaged (All Mean). These zip files contain binary files for each month. 
   + Extract the data into the GRAM Suite folders GRAM/Earth/data/MERRA2data.  This is the default location.  Afterwards, the MERRA2data folder should contain nine folders (00Z, 03Z, 06Z, 09Z, 12Z, 15Z, 18Z, 21Z, All Mean), the MERRA2info.txt file, and a README.md markdown file.
   + Extract the data into a location of your choice, say /mypath/MERRA2data, and set MERRA2Path = '/mypath/MERRA2data' in the NAMELIST file.  Afterwards, the MERRA2data folder should contain nine folders, the MERRA2info.txt file, and a README.md markdown file.

The *EarthGRAM_NCEPdata.zip* bundle contains some large data files that are critical to the proper operation of the Earth-GRAM program.  The Earth-GRAM user has two options for integrating the data files.
   + Extract the data into the GRAM Suite folders GRAM/Earth/data/NCEPdata.  This is the default location.  Afterwards, the NCEPData folder should contain two folders, ACSCII and FixedBin, and a README.md markdown file.
   + Extract the data into a location of your choice, say /mypath/NCEPdata, and set NCEPPath = '/mypath/NCEPdata' in the NAMELIST file.  Afterwards, the NCEPData folder should contain two folders, ACSCII and FixedBin, and a README.md markdown file.

To build EarthGRAM using a Linux/MinGW makefile, navigate to the folder /GRAM/Builds and follow the README instructions.

To build EarthGRAM using MSVC++ 2017, follow the instructions in the README file in \GRAM\MSVS.

Files and Subfolders:
- data: Contains model data for the EarthGRAM program.
- examples: Examples and the GRAM program for this model.
- include: Header files for the model.
- sample_inputs: Sample input and output files.
- source: Souce code for the model.
- unittest: Source code for unit tests.
- md files: Markdown files used to build the Programmer's Guide.

