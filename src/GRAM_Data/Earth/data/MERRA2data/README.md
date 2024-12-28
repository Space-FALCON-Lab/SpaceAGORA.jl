# EarthGRAM MERRA-2 data

**Multiple downloads are required to run Earth-GRAM.**  The GRAM Suite zip bundle contains programs, code, and some data required by the GRAM models.  Earth-GRAM also requires lower atmosphere data files for either the MERRA-2 model or the NCEP model.

The MERRA-2 data 2.0-degree files are divided into 9 separate data links: *MERRA2data_00Z.zip, MERRAdata_03Z.zip, MERR2data_06Z.zip, MERRA2data_09Z.zip, MERRA2data_12Z.zip, MERRA2data_15Z.zip, MERRA2data_18Z.zip, MERRA2data_21Z.zip*, and *MERRA2data_All_Mean.zip*.  The 9 files, at 2.0-degree resolution, are for each time of day (00Z, 03Z, 06Z, 09Z, 12Z, 15Z, 18Z, 21Z) and the file for all times averaged (All Mean). These zip files contain binary files for each month. 
   + Extract the data into the GRAM Suite folders GRAM/Earth/data/MERRA2data.  This is the default location.  Afterwards, the MERRA2data folder should contain nine folders (00Z, 03Z, 06Z, 09Z, 12Z, 15Z, 18Z, 21Z, All mean), the MERRA2info.txt file, and a README.md markdown file.
   + Extract the data into a location of your choice, say /mypath/MERRA2data, and set MERRA2Path = '/mypath/MERRA2data' in the NAMELIST file.  Afterwards, the MERRA2data folder should contain nine folders, the MERRA2info.txt file, and a README.md markdown file.

Files and Subfolders (after download):
- 00Z:  Binary data folder.
- 03Z:  Binary data folder.
- 06Z:  Binary data folder.
- 09Z:  Binary data folder.
- 12Z:  Binary data folder.
- 15Z:  Binary data folder.
- 18Z:  Binary data folder.
- 21Z:  Binary data folder.
- All Mean:  Binary data folder.
- MERRA2info.txt: Contains a description of the binary data.
- README.md: This file.
