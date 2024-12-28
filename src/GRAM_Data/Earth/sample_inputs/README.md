# EarthGRAM

Each of the input files contains a SPICE data path that should be edited to match your system.

Files and Subfolders:
- ref_input.txt: Reference namelist inputs for EarthGRAM using the small (2 degree) MERRA-2 data set.
     + ref_LIST.md: List style output produced from ref_input.txt.
     + ref_OUTPUT.csv: CSV style output produced from ref_input.txt.
- merra2_input.txt: Reference namelist inputs for EarthGRAM using the large (0.5 degree) MERRA-2 data set.
     + merra2_LIST.md: List style output produced from merra2_input.txt.
     + merra2_OUTPUT.csv: CSV style output produced from merra2_input.txt.
- ncep_input.txt: Reference namelist inputs for EarthGRAM using the NCEP data set.
     + ncep_LIST.md: List style output produced from ncep_input.txt.
     + ncep_OUTPUT.csv: CSV style output produced from ncep_input.txt.
- traj_input.txt: Reference namelist inputs for EarthGRAM with trajectory inputs using the small (2 degree) MERRA-2 data set.
     + traj_data.txt: Input trajectory data for traj_input.txt.
     + traj_LIST.md: List style output produced from traj_input.txt.
     + traj_OUTPUT.csv: CSV style output produced from traj_input.txt.
- aux_input.txt: Reference namelist inputs for EarthGRAM with auxiliary atmosphere using the small (2 degree) MERRA-2 data set.
     + aux_data.txt: Input auxiliary atmosphere data for aux_input.txt.
     + aux_LIST.md: List style output produced from aux_input.txt.
     + aux_OUTPUT.csv: CSV style output produced from aux_input.txt.
- rra_input.txt: Reference namelist inputs for EarthGRAM using a range reference atmosphere using the small (2 degree) MERRA-2 data set.
     + rra_LIST.md: List style output produced from rra_input.txt.
     + rra_OUTPUT.csv: CSV style output produced from rra_input.txt.
- cm_input.txt: Reference namelist inputs for EarthCorrMonte using the small (2 degree) MERRA-2 data set.
     + cm_LIST.md: List style output produced from cm_input.txt.
     + cm_OUTPUT.csv: CSV style output produced from cm_input.txt.
     + Nominal_cm_LIST.md: List style nominal output produced from cm_input.txt.
     + Nominal_cm_OUTPUT.csv: CSV style nominal output produced from cm_input.txt.
- find_dates.txt: Sample namelist inputs for finding LS/LTST dates using EarthGRAM.
- spice.txt: Namelist input file for overriding the default SPICE data path.