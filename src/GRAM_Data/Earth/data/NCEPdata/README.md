# EarthGRAM NCEP data

**Two downloads are required to run Earth-GRAM.**  The GRAM Suite zip bundle contains programs, code, and most data required by the GRAM models.  The EarthGRAM_NCEP_Data.zip bundle contains some large data files that are critical to the proper operation of the Earth-GRAM program.  The Earth-GRAM user has two options for integrating the data files.

+ Extract the data into the folder GRAM/Earth/data.  This is the default location.  Afterwards, the NCEPData folder should contain two folders, ACSCII and FixedBin, and this README.md markdown file.  During the extraction process, you may be warned about overwriting files.
+ Extract the data into a location of your choice, say /mypath/NCEPdata, and set NCEPPath = '/mypath/NCEPdata' in the NAMELIST file.  Afterwards, the NCEPData folder should contain two folders, ACSCII and FixedBin, and a README.md markdown file.

Files and Subfolders (after download):
- ASCII: Readable text form of the NCEP data. Conversion utility programs.
- FixedBin: Binary form of the NCEP data.
