The GRAM C Interface
====================

This documents the C interface for the GRAM framework.  Since all GRAM models are built around this framework, the contents are apropos for each model.  When refering to a model specific function or item, simply replace <em>Body</em> with the model of choice (e.g. NeptuneAtmosphere_C for <em>Body</em>Atmosphere_C). Also, functions ending with "_B" (B for Body) should replace with the "B" with the first letter of the model of choice (e.g Neptune would use update_N instead of update_B).  This naming convention is to avoid name collision. All structures will end in "_C" to distinguish them from the C++ structures and the FORTRAN structures.  The C interface also includes a generic (not planet specific) interface in which all functions end with "_C".

- - - - - - - - - - - 	
The Interface Functions and Structures
---------------------
- - - - - - - - - - - 	


- C Interfaces for Atmosphere Models
	- [Venus (_V)](Venus/Venus_C.md)
	- [Earth (_E)](Earth/Earth_C.md)
	- [Mars (_M)](Mars/Mars_C.md)
	- [Jupiter (_J)](Jupiter/Jupiter_C.md)
	- [Titan (_T)](Titan/Titan_C.md)
	- [Uranus (_U)](Uranus/Uranus_C.md)
	- [Neptune (_N)](Neptune/Neptune_C.md)
	- [The Generic C Interface (_C)](Documentation/Markdown/generic_C.md)
- Atmosphere Model Functions
	+ setSpiceLsk_B
	+ setSpicePck_B
	+ setSpiceKernel_B
	+ initialize_B
	+ createAtmosphere_B
	+ deleteAtmosphere_B
	+ setStartTime_B
	+ setPerturbationFactors_B
	+ addAuxiliaryAtmosphere_B
	+ setPosition_B
	+ setDelta_B
	+ setPerturbationAction_B
	+ setEphemerisState_B
	+ update_B
	+ getStartTime_B
	+ getPerturbationFactors_B
	+ getPosition_B
	+ getDynamicsState_B
	+ getDensityState_B
	+ getWindsState_B
	+ getGasesState_B
	+ getEphemerisState_B
	+ getPerturbationState_B
- Utility Functions
	+ getVersionString_B
	+ loadSpiceFile_B
- Data Structures
	+ <em>Body</em>Atmosphere_C
		- [VenusAtmosphere_C](@ref GRAM::VenusAtmosphere_C)
		- [EarthAtmosphere_C](@ref GRAM::EarthAtmosphere_C)
		- [MarsAtmosphere_C](@ref GRAM::MarsAtmosphere_C)
		- [JupiterAtmosphere_C](@ref GRAM::JupiterAtmosphere_C)
		- [TitanAtmosphere_C](@ref GRAM::TitanAtmosphere_C)
		- [UranusAtmosphere_C](@ref GRAM::UranusAtmosphere_C)
		- [NeptuneAtmosphere_C](@ref GRAM::NeptuneAtmosphere_C)
	+ [GramTime_C](@ref GRAM::GramTime_C)
	+ [Position_C](@ref GRAM::Position_C)
	+ [DynamicsState_C](@ref GRAM::DynamicsState_C)
	+ [DensityState_C](@ref GRAM::DensityState_C)
	+ [WindsState_C](@ref GRAM::WindsState_C)
	+ [GasesState_C](@ref GRAM::GasesState_C)
	+ [EphemerisState_C](@ref GRAM::EphemerisState_C)
	+ [PerturbationState_C](@ref GRAM::PerturbationState_C)
	+ [ConstituentGas_C](@ref GRAM::ConstituentGas_C)

- - - - - - - - - - - 	
Usage of the C Interface
------------------------
- - - - - - - - - - - 	

A simple example of using the C interface can be seen [here](@ref Neptune/Examples/Neptune_C.c).

Include the header for the C interface.  This will declare all functions and data structures.

~~~~~~~~~~
#include "BodyAtmosphere_C.h"
~~~~~~~~~~

### Setup - Call these functions once

*First and Foremost:* The GRAM models use the NAIF SPICE library for time and ephemeris computations.  The location of the Spice data must be set before any other call into the GRAM models.  

One can optionally override the default SPICE kernels.  This should be done before initializing.
~~~~~~~~~~
  setSpiceLsk_B("/lsk/naif0012.tls");
  setSpicePck_B("/pck/pck00011.tpc");
  setSpiceKernel_B("/spk/satellites/bod123.bsp");
~~~~~~~~~~

Call the initialize function:
~~~~~~~~~~
  initialize_B("/my/path/to/spice/data");
~~~~~~~~~~

Next, create an atmosphere model. Do this by declaring a pointer to the <em>Body</em>Atmosphere_C structure and calling the createAtmosphere_B function. The function returns a [handle](https://en.wikipedia.org/wiki/Handle_(computing)) to the atmosphere model that will be needed by all successive calls.  It is okay to create more than one atmosphere handle.
~~~~~~~~~~
  BodyAtmosphere_C* body = createAtmosphere_B();
~~~~~~~~~~

The atmosphere model needs to be initialized.  First, set the start time.  In GRAM models, time is computed using the number of elapsed seconds past the start time.
~~~~~~~~~~
  GramTime_C time;
  time.year = 2020;
  time.month = 3;
  time.day = 15;
  time.hour = 0;
  time.minute = 0;
  time.seconds = 0.0;
  time.frame = 1;   // ERT
  time.scale = 1;   // UTC
  setStartTime_B(body, &time);
~~~~~~~~~~

Set the factors that control random perturbations.
~~~~~~~~~~
  setSeed_B(body, 1234);
  setMinRelativeStepSize_B(body, 0.05);
  setPerturbationScales_B(body, 1.5, 1.5, 1.5, 1.5);
~~~~~~~~~~

Some models may require further setup.  See the model specific documentation for details.

- - - - - - - - - - - 	
### In the loop - Call these methods for each time step of the simulation.

Set the position and elapsed time.
~~~~~~~~~~
  Position_C position;
  position.height = 50.0;       // km
  position.latitude = 22.0;     // degrees
  position.longitude = 48.0;    // degrees, east positive
  position.elapsedTime = 100.0; // seconds
  position.isGeocentric = 1;    // 1 is for geocentric, set to 0 for geodetic
  setPosition_B(body, position);
~~~~~~~~~~

Perform the model computations with an update.
~~~~~~~~~~
  update_B(body);
~~~~~~~~~~

Retrieve the updated values
~~~~~~~~~~
  getPosition_B(body, &position);         // contains position related metrics
  DynamicsState_C dynamics;
  getDynamicsState_B(body, &dynamics);    // atmosphere model values
  DensityState_C densities;
  getDensityState_B(body, &densities);    // atmosphere density values
  WindsState_C winds;
  getWindsState_B(body, &winds);          // atmosphere wind values
  GasesState_C gases;
  getGasesState_B(body, &gases);          // atmosphere constituent gas values
  EphemerisState_C ephem;
  getEphemerisState_B(body, &ephem);      // planetary position values
  PerturbationState_C perts;
  getPerturbationState_B(body, &perts);   // (optional) random numbers used in perturbations
~~~~~~~~~~
The call to getGasesState_B may be diffenent for each model.  See the model specific documentation for the precise call signature.

Process the information.  Update the position and elapsed time.  And repeat.

- - - - - - - - - - - 	
### On Shutdown

The createAtmosphere_B() function allocated memory which needs to be freed.  Once deleteAtmosphere_B() has been called, the body handle becomes invalid.  It is a good practice to nullify the handle.
~~~~~~~~~~
  deleteAtmosphere_B(body);
  body = NULL;
~~~~~~~~~~
