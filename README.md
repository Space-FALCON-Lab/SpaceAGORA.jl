# ABTS.jl
The Julia Aerobraking Trajectory Simulator is the successor to the Python-based Aerobraking Trajectory Simulator ([link](https://github.com/Space-FALCON-Lab/Aerobraking-Trajectory-Simulator)). This simulator can be used to model a wide range of atmospheric and non-atmospheric missions with high-fidelity models of perturbing effects such as third-body gravity, solar radiation pressure, and gravitational harmonics. In addition, this simulator includes the capability to use high-fidelity atmospheric models through the GRAM Suite, with the option to use spherical harmonic topography models to determine altitude. This iteration of the simulator provides order-of-magnitude speed improvements compared to the Python version, as well as the inclusion of a wider range of perturbations and planets. Finally, the simulator provides estimates of several important physical parameters in atmospheric flight, including dynamic pressure, heat rate, heat load, and energy depletion rate. This mean it is not only an aerobraking mission simulator, but is also easily extensible to other mission profiles, such as entry and aerocapture. 

# Setting up the Docker Environment (Currently not working for Apple Silicon Macs)
To ensure that all the package versions are consistent, a Docker environment has been configured for use with ABTS.jl. This guarantees consistent results between computers. This is needed for running on Apple Silicon Macs because of how GRAM was compiled. For other virtualization options, see [GRAM Setup](https://github.com/Space-FALCON-Lab/Aerobraking-Trajectory-Simulator/tree/GRAM-updates?tab=readme-ov-file#gram-setup). The process of setting it up is as follows:
1. Download and sign in to Docker Desktop ([link](https://www.docker.com/get-started/)). If on Linux, follow these instructions to sign in ([link](https://docs.docker.com/desktop/setup/sign-in/))
2. In VSCode, download the Docker and Dev Container extensions
3. Open the ABTS.jl directory in VSCode and click on the remote window icon in the bottom left corner
4. Select "Reopen in container". This will open a window with the ABTS.jl repository open in the preconfigured Docker environment.
5. Open a terminal window in VSCode. Start Julia, enter "]" to activate the package environment, and enter ```activate .ABTS```, then ```instantiate```. This will activate the Julia project environment that has been set up and install all the packages needed to run ABTS.jl
6. Now, you should be able to run any existing scenario or create and run a new scenario

# GRAM access
If part of the Space-FALCON Lab, access through the lab's Drive in ABTS/GRAM. Download GRAM_Data.zip and GRAMpy.zip and extract them to the root folder of your local repository (i.e., you should have two new folders: ABTS.jl/GRAM_Data/ and ABTS.jl/GRAMpy/).

If not part of Space-FALCON lab, GRAM may be requested [here](https://software.nasa.gov/software/MFS-33888-1). Follow the instructions to build the Python version and follow the same steps as above to get all the data.

# Getting Started
Start by cloning this repository:

```cmd
git clone https://github.com/Space-FALCON-Lab/ABTS.jl
```

Next, follow the instructions above to download the GRAM Suite, which is used for high-fidelity atmospheric modeling. This will include all the GRAM and SPICE files required to properly run the simulations. Example code is provided in ```src/ABTS_*.jl```.

If GRAM is installed properly and the included Docker environment is being used, this code will run as-is. If for some reason you are unable to get access to GRAM, make sure that the following SPICE files are present in the folder designated in the ```directory_Spice``` argument, organized into subfolders ```pck```, ```lsk```, ```spk/planets```, etc.:
> PCK files
> * ```pck00011.tpc```: Orientation, size, and shape data

> SPK files
>> Planets
>> * ```de440s.bsp```: Planetary barycenters w.r.t. solar system barycenter

>> Satellites
>> * ```sat441.bsp```: Saturn/Titan positions

> LSK files
>> * ```naif0012.tls```: Timing data

The following described necessary modifications to the example code in several cases:
> GRAM is not installed properly
* Change the ```density_model``` argument to ```Constant``` or ```Exponential```. This will change the model used to calculate the atmospheric density at each step from GRAM to either an exponential model or a constant density model.
* Change the ```directory_GRAM``` argument to ```""```.
> Docker is not properly set up
* Change the file path arguments (```directory_Gram```, ```directory_Spice```, etc.) to match your system. Curently they are configured for the Docker environment assuming all GRAM and SPICE files are present in the correct locations.
