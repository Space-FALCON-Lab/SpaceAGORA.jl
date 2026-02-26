The GRAM C++ Interface
======================

This documents the usage of the primary interface classes in the GRAM C++ framework.  Since all GRAM models are built around this framework, the contents are apropos for each model.  When refering to a model specific class or item, simply replace <em>Body</em> with the model of choice (e.g. NeptuneAtmosphere for <em>Body</em>Atmosphere).  For details on a specific class or method, please refer to the class documentation.

Quick Links
===========

- [The Interface Classes](@ref a0)
- [Usage of the BodyAtmosphere and Data Classes](@ref a1)
- [Usage of the Data Generation Classes](@ref a2)
- [Usage of the Utility Classes](@ref a3)

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	
The Interface Classes    {#a0}
=====================
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	

Simulations which control time and position will probably only need the Atmosphere Model and Data classes.  Data Generation and Utility classes are provided for simulations which require a full trajectory.

- [Atmosphere Model Classes](#a1)
	+ <em>Body</em>Atmosphere
		- [VenusAtmosphere](Venus/Venus.md)
		- [EarthAtmosphere](Earth/Earth.md)
		- [MarsAtmosphere](Mars/Mars.md)
		- [JupiterAtmosphere](Jupiter/Jupiter.md)
		- [TitanAtmosphere](Titan/Titan.md)
		- [UranusAtmosphere](Uranus/Uranus.md)
		- [NeptuneAtmosphere](Neptune/Neptune.md)
- [Data Classes](#a1)
	+ <em>Body</em>InputParameters
		- [InputParameters](@ref GRAM::InputParameters)
	+ [GramTime](@ref GRAM::GramTime)
	+ [Position](@ref GRAM::Position)
	+ [AtmosphereState](@ref GRAM::AtmosphereState)
	+ [EphemerisState](@ref GRAM::EphemerisState)
	+ [PerturbationState](@ref GRAM::PerturbationState)
	+ [ConstituentGas](@ref GRAM::ConstituentGas)
- [Data Generation Classes](#a2)
	+ [Trajectory](@ref GRAM::Trajectory)
	+ [VerticalProfile](@ref GRAM::VerticalProfile)
	+ [EphemerisProfile](@ref GRAM::EphemerisProfile)
	+ [MonteCarlo](@ref GRAM::MonteCarlo)
- [Utility Classes](#a3)
	+ [SpiceLoader](@ref GRAM::SpiceLoader)
	+ <em>Body</em>NamelistReader
		- [NamelistReader](@ref GRAM::NamelistReader)
	+ [ProfilePrinter](@ref GRAM::ProfilePrinter)

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	
Usage of the BodyAtmosphere and Data Classes  {#a1}
============================================
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	

A simple example of using an atmosphere class can be seen [here](@ref Neptune/Examples/Neptune.cpp).

Of course, headers for the appropriate classes must be included.  Header files typically have the same name as the class.  Including the atmosphere model will also include all of the data classes.

~~~~~~~~~~
#include "SpiceLoader.h"
#include "BodyAtmosphere.h"
~~~~~~~~~~

All of the atmosphere models are in the GRAM namespace.
~~~~~~~~~~
using namespace GRAM;
~~~~~~~~~~

### Setup - Call these methods once

*First and Foremost:* The GRAM models use the NAIF Spice library for time and ephemeris computations.  The location of the Spice data must be set before any other call into the GRAM models.  This can be done in one of two ways. (Note that using both methods will result in the latter overriding the former.)

One can optionally override the default SPICE kernels.  This should be done before calling setSpiceDataPath.
~~~~~~~~~~
  SpiceLoader::setSpiceLsk("/lsk/naif0012.tls");
  SpiceLoader::setSpicePck("/pck/pck00011.tpc");
  SpiceLoader::setSpiceKernel(NEPTUNE, "/spk/satellites/nep101.bsp");
  SpiceLoader::setSpiceKernel(URANUS, "/spk/satellites/ura116.bsp");
~~~~~~~~~~
Use the SpiceLoader static initialization method:
~~~~~~~~~~
  SpiceLoader::setSpiceDataPath("/my/path/to/Spice/data");
~~~~~~~~~~

*or*

Use the <em>Body</em>InputParameters class. This method is preferred if other input parameters are modified.
~~~~~~~~~~
  BodyInputParameters inputParameters;
  inputParameters.spicePath = "/my/path/to/Spice/data";
  inputParameters.spiceLsk = "/lsk/naif0012.tls";  // optional
  inputParameters.spicePck = "/pck/pck00011.tpc";  // optional
  inputParameters.spiceNeptune = "/spk/satellites/nep101.bsp";  // optional
  inputParameters.spiceUranus = "/spk/satellites/ura116.bsp";   // optional
  BodyAtmosphere body;
  body.setInputParameters(inputParameters);
~~~~~~~~~~

Next, create an atmosphere model. (The second option above has already done this.)
~~~~~~~~~~
  BodyAtmosphere body;
~~~~~~~~~~

The atmosphere model needs to be initialized.  First, set the start time.  In GRAM models, time is computed using the number of elapsed seconds past the start time.
~~~~~~~~~~
  GramTime time;
  time.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT);
  body.setStartTime(time);
~~~~~~~~~~

Choose a seed and set parameters for the random perturbations.
~~~~~~~~~~
  body.setSeed(1234);
  body.setPerturbationScaleFactor(1.5);
  body.setMinRelativeStepSize(0.5);
~~~~~~~~~~

Some models may require further setup.  See the model specific documentation for details.

- - - - - - - - - - - 	
### In the loop - Call these methods for each time step of the simulation.

Set the position and elapsed time.
~~~~~~~~~~
  Position position;
  position.height = 50.0;       // km
  position.latitude = 22.0;     // degrees
  position.longitude = 48.0;    // degrees, east positive
  position.elapsedTime = 100.0; // seconds
  position.isGeocentric = true; // optional, true is the default, set to false for geodetic
  body.setPosition(position);
~~~~~~~~~~

Perform the model computations with an update.
~~~~~~~~~~
  body.update();
~~~~~~~~~~

Retrieve the updated values
~~~~~~~~~~
 position = body.getPosition(); // contains position related metrics
 const AtmosphereState& atmos = body.getAtmosphereState();    // atmosphere model values
 const EphemerisState ephem& = body.getEphemerisState();      // planetary position values
 const PerturbationState pert& = body.getPerturbationState(); // (optional) random numbers used in perturbations
~~~~~~~~~~

Process the information.  Update the position and elapsed time.  And repeat.

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	
Usage of the Data Generation Classes   {#a2}
====================================
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	

The GRAM::Profile class is a base class for building one-dimensional atmosphere profiles. The interface classes for profiles are:
- GRAM::VerticalProfile: Uses a starting position and a change in height to generate a vertical column of data.
- GRAM::Trajectory: Two options available to generate data along a curve.
    + Use a starting point and a change in position.
	+ Read positions from a file.
- GRAM::EphemerisProfile: Set a starting ephemeris and a change in ephemeris to generate a sweep of values over the course of a day or season.

Always start by setting the Spice data path:
~~~~~~~~~~
  SpiceLoader::setSpiceDataPath("/my/path/to/Spice/data");
~~~~~~~~~~

Next, create an atmosphere model. And set any appropriate parameters.
~~~~~~~~~~
  BodyAtmosphere body;
  body.setSeed(1234);
  body.setPerturbationScaleFactor(1.5);
  body.setMinRelativeStepSize(0.5);
  GramTime time;
  time.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT);
  body.setStartTime(time);
~~~~~~~~~~

Create a profile (here we use the VerticalProfile) and set the atmosphere.
~~~~~~~~~~
  VerticalProfile profile;
  profile.setAtmosphere(body);
~~~~~~~~~~

Now set the profile specific parameters.  For a VerticalProfile:
~~~~~~~~~~
  // Set the initial position.
  Position position;
  position.height = 50.0;       // km
  position.latitude = 22.0;     // degrees
  position.longitude = 48.0;    // degrees, east positive
  position.elapsedTime = 0.0;   // seconds
  position.isGeocentric = true; // optional, true is the default, set to false for geodetic
  profile.setInitialPosition(position);
  
  // Set the change in height and number of points to generate.
  profile.setHeightConditions(50.0, 20);
~~~~~~~~~~

For a Trajectory:
~~~~~~~~~~
  // Set the initial position.
  Position position;
  position.height = 50.0;       // km
  position.latitude = 30.0;     // degrees
  position.longitude = 48.0;    // degrees, east positive
  position.elapsedTime = 0.0;   // seconds
  position.isGeocentric = true; // optional, true is the default, set to false for geodetic
  profile.setInitialPosition(position);

  // Set the change in position (including time).
  Position delta;
  delta.height = 5.0;        // km
  delta.latitude = -1.0;     // degrees
  delta.longitude = 0.0;     // degrees, east positive
  delta.elapsedTime = 10.0;  // seconds
  profile.setDeltaPosition(position);

  profile.setNumberOfPoints(61);
~~~~~~~~~~
  
For an EphemerisProfile:
~~~~~~~~~~
  // Set the position.
  Position position;
  position.height = 50.0;       // km
  position.latitude = 22.0;     // degrees
  position.longitude = 48.0;    // degrees, east positive
  position.elapsedTime = 0.0;   // seconds
  profile.setPosition(position);

  // Set the initial ephemeris values.
  EphemerisState initialEphemeris;
  initialEphemeris.longitudeSun = 0.0_deg;
  initialEphemeris.solarTime = 12.0_hr;
  profile.setInitialEphemeris(initialEphemeris);

  // Set the change in ephemeris values.
  EphemerisState deltaEphemeris;
  deltaEphemeris.longitudeSun = 10.0;
  deltaEphemeris.solarTime = 0.0;
  profile.setDeltaEphemeris(deltaEphemeris);

  // Choose the number of datasteps to generate
  profile.setNumberOfPoints(37);
~~~~~~~~~~

Generate the profile data.
~~~~~~~~~~
  profile.generate();
~~~~~~~~~~

The data is stored in a vector of ProfileData objects. Get the data and process as desired.  Here the data is passed to a ProfilePrinter. 
~~~~~~~~~~
  const std::vector<ProfileData>& myData = profile.getProfile();
  ProfilePrinter printer;
  printer.setStyle(printer.GRAM_COLUMN_STYLE | printer.GRAM_LIST_STYLE);
  printer.print(body, myData);
~~~~~~~~~~

The GRAM::MonteCarlo class will generate any number of profiles (VerticalProfile, Trajectory, EphemerisProfile, or custom subclass of Profile) with a new pertubation seed for each profile.

A MonteCarlo object will need a Profile and a ProfilePrinter object preconfigured for use.  Start by creating the objects that will be needed.
~~~~~~~~~~
  MonteCarlo monte;
  Trajectory trajectory;
  ProfilePrinter printer;
  BodyAtmosphere body;
~~~~~~~~~~

Initialize each of the objects and pass them to the MonteCarlo.
~~~~~~~~~~
  // A start time is needed
  GramTime ttime;
  ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT);
  body.setStartTime(ttime);

  // Also a starting position
  Position initialPosition;
  initialPosition.height = 0.0;
  initialPosition.latitude = 22.0;
  initialPosition.longitude = 88.0;
  initialPosition.elapsedTime = 0.0;
  trajectory.setInitialPosition(initialPosition);

  // And a change in position
  Position deltaPosition;
  deltaPosition.height = 10.0;
  deltaPosition.latitude = 0.30;
  deltaPosition.longitude = -0.50;
  deltaPosition.elapsedTime = 500.00;
  trajectory.setDeltaPosition(deltaPosition);

  // The number of data steps in each trajectory
  trajectory.setNumberOfPoints(201);

  // The Monte Carlo also needs a profile printer
  printer.setStyle(printer.GRAM_MONTE_CARLO_STYLE);
  
  // Assign the objects to the Monte Carlo object
  monte.setProfile(trajectory);
  monte.setProfilePrinter(printer);
~~~~~~~~~~

The initial seed for random perturbations can be set or a list of seeds can be provided.  If the initial seed is set, then subsequent seeds will be generated automatically.
~~~~~~~~~~
  // Set the initial seed.
  monte.setInitialSeed(1234);
  
  // OR, provide a std::vector<int> of seeds.
  monte.setSeedList(mySeedList);
~~~~~~~~~~

Finally, generate the trajectories.
~~~~~~~~~~
  monte.generate();
~~~~~~~~~~

The seeds used to generate the trajectories can be retreived.
~~~~~~~~~~
  const std::vector<int>& seeds = monte.getSeedList();
~~~~~~~~~~

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	
Usage of the Utility Classes   {#a3}
============================
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	

As seen in the example above, the GRAM::ProfilePrinter is used to write profile data to a file.  The ProfilePrinter class can be used as is or it can be subclassed to provide custom output styles.  Basic usage of the ProfilePrinter is fairly simple.
~~~~~~~~~~
  // Instantiate a printer object.
  ProfilePrinter printer;
  
  // Set the desired output styles by createing a bitmask with bitwise "or".
  printer.setStyle(printer.GRAM_COLUMN_STYLE | printer.GRAM_LIST_STYLE);
  
  // Print the data to a file.
  printer.print(body, myData);
~~~~~~~~~~

Each style is output to a separate file. The filenames are automatic.  Collision of filenames can be avoided by setting a unique prefix before printing.
~~~~~~~~~~
 printer.setFileNamePrefix("Example1_");
~~~~~~~~~~

Changing the longitude convention for output is possible.
~~~~~~~~~~
 printer.setWestLongitudePositive();
~~~~~~~~~~

The primary output styles are now in CVS (comma seperated value) format or MD (markdown) format. Legacy output styles have been maintained for comparison of data between versions. Legacy output styles are in a space delimited column format.  Here are the available print styles:

| Style                  | Description                                                             |
|------------------------|-------------------------------------------------------------------------|
| GRAM_CSV_STYLE         | All data output in a comma separated value file.                        |
| GRAM_MD_STYLE          | A long list of parameters formatted using markdown syntax.              |
| GRAM_COLUMN_STYLE      | Density, temperature, wind, and gas data.                               |
| GRAM_LIST_STYLE        | A long list of parameters in a paragraph format.                        |
| GRAM_EPHEMERIS_STYLE   | Ephemeris data.                                                         |
| GRAM_DENSITY_STYLE     | Density ranges, mean and perturbed density.                             |
| GRAM_PERTURB_STYLE     | Standard deviations for density and winds.                              |
| GRAM_WINDS_STYLE       | Mean and perturbed wind data.                                           |
| GRAM_TPRESHGT_STYLE    | Temperature, pressure, scale heights, and gas mole fractions.           |
| GRAM_SOUND_STYLE       | Temperature, pressure, density, speed of sound, and gas mass fractions. |
| GRAM_MONTE_CARLO_STYLE | A combination of DENSITY, PERTURB, WINDS, and TPRESHGT styles.          |
| GRAM_ALL_STYLES        | A combination of all styles.                                            |

The GRAM::SpiceLoader class is a C++ interface to the Spice furnsh function.  Spice data is loaded into memory by passing a file name to the furnsh function.  Multiple calls to furnsh function with the identical file name will result in memory bloat.  The SpiceLoader class ensures the furnsh calls are unique.

The GRAM codebase automates the loading of the necessary Spice data.  However, a simulation that is also using Spice may want to take advantage of the SpiceLoader to avoid overloading the Spice memory.

The first call into the GRAM library should be to set the path to the Spice data files.
~~~~~~~~~~
 // Using a static function.
 SpiceLoader::setSpiceDataPath("/my/path/to/Spice/data");
 
 // Or instantiating a SpiceLoader.
 SpiceLoader spiceLoader;
 spiceLoader.setSpiceDataPath("/my/path/to/Spice/data");
~~~~~~~~~~

Loading a Spice file will ignore calls that duplicate an existing file.
~~~~~~~~~~
 // Using a static function.
 SpiceLoader::loadFile("/lsk/naif0012.tls");
 
 // Or instantiating a SpiceLoader.
 SpiceLoader spiceLoader;
 spiceLoader.loadFile("/lsk/naif0012.tls");
~~~~~~~~~~ 

The GRAM::NamelistReader class supports the legacy capability of reading FORTRAN namelist files as a means of loading GRAM parameters.  Each <em>Body</em>Atmosphere class will subclass the NamelistReader as <em>Body</em>NamelistReader to handle body specific parameters.  The NamelistReaders support the legacy input field names as well as the new variable names.

~~~~~~~~~~ 
  // Create a namelist reader.
  BodyNamelistReader reader;
  
  // Get the input parameters.
  BodyInputParameters params;
  reader.getParameters("myInputs.txt", params);
  
  // Set the parameters.
  BodyAtmosphere body;
  body.setInputParameters(params);
~~~~~~~~~~ 
