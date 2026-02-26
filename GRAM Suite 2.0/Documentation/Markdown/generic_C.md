The Generic GRAM C Interface
============================

This documents the generic (not planet specific) C interface for the GRAM Suite.  This interface is provided to enable simulations that require multiple planet atmospheres operated via a single interface.  Documentation of all functions and structures appearing in this interface can be found here: \ref C_GRAM.

- - - - - - - - - - - 	
A Programmer's Warning
---------------------
This interface is *not type safe*.  The interface utilizes void pointers as handles for the atmosphere objects.  The void pointers are cast to the appropriate objects in the C++ layer of the interface.  It is the responsibility of the developer to ensure the validity of the void pointers in order to avoid catastrophe.

- - - - - - - - - - - 	
Usage of the Generic C Interface
--------------------------------
- - - - - - - - - - - 	

Include GRAM_C.h, the header for the C interface.  This will declare all functions and data structures.

~~~~~~~~~~
#include "GRAM_C.h"
~~~~~~~~~~

### Setup - Call these functions once

*First and Foremost:* The GRAM models use the NAIF Spice library for time and ephemeris computations.  The location of the Spice data must be set before any other call into the GRAM models.  

One can optionally override the default SPICE kernels.  This should be done before initializing.
~~~~~~~~~~
  setSpiceLsk_C("/lsk/naif0012.tls");
  setSpicePck_C("/pck/pck00011.tpc");
  setSpiceKernel_C(NEPTUNE_C, "/spk/satellites/nep101.bsp");
  setSpiceKernel_C(URANUS_C, "/spk/satellites/ura116.bsp");
~~~~~~~~~~

Call the initialize function:
~~~~~~~~~~
  initialize_C("/my/path/to/Spice/data");
~~~~~~~~~~

Next, create an atmosphere model. Do this by declaring a void pointer and calling the createAtmosphere_C function. The function returns a [handle](https://en.wikipedia.org/wiki/Handle_(computing)) to the atmosphere model that will be needed by all successive calls.  It is okay to create more than one atmosphere handle for all available atmospheres.  The type of the atmosphere must be provided from the \ref GRAM_BODY_C enumeration.
~~~~~~~~~~
  void* neptune = createAtmosphere_C(NEPTUNE_C);
  void* uranus = createAtmosphere_C(URANUS_C);
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
  setStartTime_C(neptune, &time);
  setStartTime_C(uranus, &time);
~~~~~~~~~~

Set the factors that control random perturbations.
~~~~~~~~~~
  setSeed_C(neptune, 1234);
  setMinRelativeStepSize_C(neptune, 0.05);
  setPerturbationScales_C(neptune, 1.5, 1.5, 1.5, 1.5);
  setSeed_C(uranus, 1234);
  setMinRelativeStepSize_C(uranus, 0.05);
  setPerturbationScales_C(uranus, 1.5, 1.5, 1.5, 1.5);
~~~~~~~~~~

Some models may require further setup.  See the model specific documentation for details.  The model specific setup functions are currently not available via the generic interface.
~~~~~~~~~~
  setMinMaxFactor_N((NeptuneAtmosphere_C*)neptune, -0.2, 1);
~~~~~~~~~~

After setup, the atmosphere models are placed into a more genereic data structure for loop processing.
~~~~~~~~~~
  void *body[2];
  body[0] = neptune;
  body[1] = uranus;
~~~~~~~~~~

- - - - - - - - - - - 	

### In the loop - Call these methods for each time step of the simulation. 

Process all bodies in the simulation.
~~~~~~~~~~
  for (int i=0; i<2; ++i) {
~~~~~~~~~~

Set the position and elapsed time.
~~~~~~~~~~
  Position_C position;
  position.height = 50.0;       // km
  position.latitude = 22.0;     // degrees
  position.longitude = 48.0;    // degrees, east positive
  position.elapsedTime = 100.0; // seconds
  position.isGeocentric = 1;    // 1 is for geocentric, set to 0 for geodetic
  setPosition_C(body[i], position);
~~~~~~~~~~

Perform the model computations with an update.
~~~~~~~~~~
  update_C(body[i]);
~~~~~~~~~~

Retrieve the updated values
~~~~~~~~~~
  getPosition_C(body[i], &position);         // contains position related metrics
  DynamicsState_C dynamics;
  getDynamicsState_C(body[i], &dynamics);    // atmosphere model values
  DensityState_C densities;
  getDensityState_C(body[i], &densities);    // atmosphere density values
  WindsState_C winds;
  getWindsState_C(body[i], &winds);          // atmosphere wind values
  GasesState_C gases;
  getGasesState_C(body[i], &gases);          // atmosphere constituent gas values
  EphemerisState_C ephem;
  getEphemerisState_C(body[i], &ephem);      // planetary position values
  PerturbationState_C perts;
  getPerturbationState_C(body[i], &perts);   // (optional) random numbers used in perturbations
~~~~~~~~~~

The set of constituent gases may be diffenent for each model.  See the model specific documentation for the set of available gases.  If a particular gas is not part of a model, it will be returned with all zeroes.  The desired gas is identified using the \ref GasType_C enumeration.
~~~~~~~~~~
  ConstituentGas_C n2, ch4, h2;
  getConstituentGas_C(body[i], DINITROGEN_C, &n2);
  getConstituentGas_C(body[i], METHANE_C, &ch4);
  getConstituentGas_C(body[i], DIHYDROGEN_C, &h2);
~~~~~~~~~~

Process the information for that body.  Then close the loop for processing all bodies.
~~~~~~~~~~
  } // end of the body loop
~~~~~~~~~~

Process the information.  Update the simulation position and elapsed time.  And repeat.

- - - - - - - - - - - 	
### On Shutdown

The createAtmosphere_C() function allocated memory which needs to be freed.  Once deleteAtmosphere_C() has been called, the body handle becomes invalid.  It is a good practice to nullify the handle.
~~~~~~~~~~
  deleteAtmosphere_C(neptune);
  neptune = NULL;
  deleteAtmosphere_C(uranus);
  neptune = uranus;
~~~~~~~~~~
