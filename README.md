# ABTS.jl
The Julia Aerobraking Trajectory Simulator is the successor to the Python-based Aerobraking Trajectory Simulator ([link](https://github.com/Space-FALCON-Lab/Aerobraking-Trajectory-Simulator)). This simulator can be used to model a wide range of atmospheric and non-atmospheric missions with high-fidelity models of perturbing effects such as third-body gravity, solar radiation pressure, and gravitational harmonics. In addition, this simulator includes the capability to use high-fidelity atmospheric models through the GRAM Suite, with the option to use spherical harmonic topography models to determine altitude. This iteration of the simulator provides order-of-magnitude speed improvements compared to the Python version, as well as the inclusion of a wider range of perturbations and planets. Finally, the simulator provides estimates of several important physical parameters in atmospheric flight, including dynamic pressure, heat rate, heat load, and energy depletion rate. This means it is not only an aerobraking mission simulator, but is also easily extensible to other mission profiles, such as entry and aerocapture. 

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

# Usage Guide
This section provides a brief overview of the structure and usage of this simulator. Generally, user interaction with the simulator is through a program similar to those given in the ```ABTS_*.jl``` files. This is where simulation settings and variables are defined and where the run function is called. The user may use this format to easily modify settings or perform a Monte Carlo analysis.
## Initial Conditions
The initial conditions are defined either through initial orbital elements or through initial velocity and flight-path angle. If running a full aerobraking campaign, the terminal condition is defined through the desired final apoapsis radius.

## Vehicle Definition
The spacecraft is modeled either as a box with two flat plates attached, representing the spacecraft main bus with a solar panel on each side, or as a blunted cone, referred to as "Spacecraft" and "Blunted Cone", respectively. In the "Spacecraft" case, the main bus is defined with a length, width, and height, and the solar panels are defined with a length and height. The overall spacecraft is defined with a dry mass, propellant mass, reflection coefficient, and thermal accomodation factor. The engine thrust, used for propulsive maneuvers to change the periapsis altitude, must also be defined. The "Blunted Cone" case is defined similarly, with the length, width, and height parameters being replaced by cone angle, nose radius, and base radius. 

For propulsive maneuvers, a function must be defined in ```utils/maneuver_plans.jl``` with the following input arguments: 
```julia
planet=nothing, ra=0.0, rp=0.0, numberofpassage=0.0, args=nothing
```
It must also modify the ```:delta_v``` and ```:phi``` fields of ```args``` and return ```args```. The ```:delta_v``` field specifies the $\Delta V$ of the propulsive maneuver, and the ```:phi``` field specifies the direction in which the thruster is oriented, with $\phi=\pi$ being a raise maneuver, and $\phi=0$ being a lower maneuver.

## Planet
The planet model is mostly predefined, with the option to specify the models used for different physical characteristics such as atmospheric density and gravity. The planet shape and size are predefined in ```src/physical_models/Planet_data.jl```. The orientation is determined using the SPICE system for accurate latitude and longitude calculations.

### Density Models
The atmospheric density calculation can be done using several built-in methods. The most accurate is to use the GRAM Suite, the set up for which is discussed previously. Other models include an exponential atmosphere model and a constant density model.

### Gravity Models
The nominal gravitational acceleration can also be computed using different built-in functions. First, a constant law sets the gravity to be constant at all altitudes, latitudes, and longitudes. Next, an inverse square law calculates the gravity as a function of altitude. Finally, an inverse square law accounting for the J2 effect determines the gravity as a function of both altitude and latitude.

### Aerodynamic and Thermal Models
The aerodynamic coefficients are Mach-dependent and use the [Flow of Rarefied Gases](https://books.google.com/books?hl=en&lr=&id=DIIrDgAAQBAJ&oi=fnd&pg=PP1&dq=rarefied+flow+schaaf+and+chambre&ots=PWLd04BJmj&sig=DaKV6gVakAuvKRgQDM3ZE9uFrdQ#v=onepage&q=rarefied%20flow%20schaaf%20and%20chambre&f=false) theory from Shaaf and Chambre, and the thermal model follows the Maxwellian heat transfer theory also from Shaaf and Chambre.

## Perturbing Forces
The most significant perturbing forces&mdash;third-body gravity, solar radiation pressure (SRP), and gravitational harmonics&mdash;are built into the simulator.

### Third-body Gravity
The gravitational effect of massive third-bodies is calculated using SPICE to determine the relative positions of each body. The third-bodies to consider are defined using the ```:n_bodies``` field, which is a vector of strings containing the SPICE names (often just the standard names, e.g., "Sun", "Moon") of the bodies to be considered. Additional SPICE files may be required if performing an analysis at destinations other than Venus, Earth, Mars, or Titan.

### Solar Radiation Pressure
Solar radiation pressure is also calculated using SPICE, this time to determine the relative position of the Sun and the central body. This is important both for determining the direction and magnitude of the SRP, but also in the case of a partial eclipse of the spacecraft, the area receiving sunlight. This is all done behind the scenes, and the only user input required is a flag that indicates whether to consider SRP.

### Gravitational Harmonics
This is the most complex perturbation, computationally. Very high degree and order models can cause significant performance loss with negligible benefit to the accuracy of the simulation. In the examples, 50th degree and order models were found to be sufficient for most purposes. This requires the user to specify the file from which the spherical harmonic coefficients are drawn and the degree and order (```:L``` and ```:M```, respectively).

## Control
The simulator currently enables two forms of control: propulsive control, discussed previously, and atmospheric-drag control. As described in [Energy Depletion Guidance for Aerobraking Atmospheric Passes](https://arc.aiaa.org/doi/abs/10.2514/1.G006171), drag-modulation trajectory control is possible through the use of rotating solar panels. This form of control is only available for the "Spacecraft" body type. The goal of this control method is to maximize the energy depletion during the drag passage without exceeding the spacecraft's thermal limits. There are thus three configurations for this type of control: 1) heat rate limited, 2) heat load limited, and 3) heat rate and heat load limited.
