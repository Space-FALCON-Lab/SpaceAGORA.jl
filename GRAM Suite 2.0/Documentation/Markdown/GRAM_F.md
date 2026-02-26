The GRAM FORTRAN Interface
====================

This documents the FORTRAN interface for the GRAM framework.  Since all GRAM models are built around this framework, the contents are apropos for each model.  When refering to a model specific function or item, simply replace <em>Body</em> with the model of choice (e.g. NeptuneAtmosphere_F for <em>Body</em>Atmosphere_F). Also, functions ending with "_B" (B for Body) should replace with the "B" with the first letter of the model of choice (e.g Neptune would use update_N instead of update_B).  This naming convention is to avoid name collision. All structures will end in "_F" to distinguish them from the C++ structures and the C structures.

- - - - - - - - - - - 	
The Interface Functions and Structures
---------------------
- - - - - - - - - - - 	


- FORTRAN Interfaces for Atmosphere Models
	- [Venus (_V)](Venus/Venus_F.md)
	- [Earth (_E)](Earth/Earth_F.md)
	- [Mars (_M)](Mars/Mars_F.md)
	- [Jupiter (_J)](Jupiter/Jupiter_F.md)
	- [Titan (_T)](Titan/Titan_F.md)
	- [Uranus (_U)](Uranus/Uranus_F.md)
	- [Neptune (_N)](Neptune/Neptune_F.md)
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
	+ [GramTime_F](@ref gramstructs::gramtime_f)
	+ [Position_F](@ref gramstructs::position_f)
	+ [DynamicsState_F](@ref gramstructs::dynamicsstate_f)
	+ [DensityState_F](@ref gramstructs::densitystate_f)
	+ [WindsState_F](@ref gramstructs::windsstate_f)
	+ [GasesState_F](@ref gramstructs::gasesstate_f)
	+ [EphemerisState_F](@ref gramstructs::ephemerisstate_f)
	+ [PerturbationState_F](@ref gramstructs::perturbationstate_f)
	+ [ConstituentGas_F](@ref gramstructs::constituentgas_f)

- - - - - - - - - - - 	
Usage of the FORTRAN Interface
------------------------
- - - - - - - - - - - 	

A simple example of using the FORTRAN interface can be seen [here](@ref Neptune/Examples/Neptune_F.f90).

The FORTRAN interface uses the ISO C binding.  The interface functions and data structures are declared in two modules.

~~~~~~~~~~
USE, INTRINSIC :: iso_c_binding
USE GramStructs
USE bodyGRAM
~~~~~~~~~~

### Setup - Call these functions once

*First and Foremost:* The GRAM models use the NAIF Spice library for time and ephemeris computations.  The location of the Spice data must be set before any other call into the GRAM models.  

One can optionally override the default SPICE kernels.  This should be done before initializing.
~~~~~~~~~~
call setSpiceLsk_B("/lsk/naif0012.tls");
call setSpicePck_B("/pck/pck00011.tpc");
call setSpiceKernel_B("/spk/satellites/bod123.bsp");
~~~~~~~~~~

Call the initialize function:
~~~~~~~~~~
call initialize_B("/my/path/to/Spice/data")
~~~~~~~~~~

Next, create an atmosphere model. Do this by declaring a pointer to the <em>Body</em>Atmosphere_F structure and calling the createAtmosphere_B function. The function returns a [handle](https://en.wikipedia.org/wiki/Handle_(computing)) to the atmosphere model that will be needed by all successive calls.  It is okay to create more than one atmosphere handle.
~~~~~~~~~~
  TYPE(C_PTR) :: body          ! Move to declaration section
  body = createAtmosphere_B()
~~~~~~~~~~

The atmosphere model needs to be initialized.  First, set the start time.  In GRAM models, time is computed using the number of elapsed seconds past the start time.
~~~~~~~~~~
  TYPE(GramTime_F) :: time          ! Move to declaration section
  time%year = 2020
  time%month = 3
  time%day = 15
  time%hour = 0
  time%minute = 0
  time%seconds = 0.0
  time%frame = 1   ! PET
  time%scale = 1   ! UTC
  call setStartTime_B(body, time)
~~~~~~~~~~

Set the factors that control random perturbations.
~~~~~~~~~~
  call setPerturbationFactors_B(body, 1.5, 0.5, 1234)
~~~~~~~~~~

Some models may require further setup.  See the model specific documentation for details.

- - - - - - - - - - - 	
### In the loop - Call these methods for each time step of the simulation.

Set the position and elapsed time.
~~~~~~~~~~
  TYPE(Position_F) :: position           ! Move to declaration section
  position%height = 50.0        ! km
  position%latitude = 22.0      ! degrees
  position%longitude = 48.0     ! degrees, east positive
  position%elapsedTime = 100.0  ! seconds
  position%isPlanetoCentric = 1 ! 1 is for planeto-centric, set to 0 for planeto-graphic
  call setPosition_B(body, position)
~~~~~~~~~~

Perform the model computations with an update.
~~~~~~~~~~
  call update_B(body)
~~~~~~~~~~

Retrieve the updated values
~~~~~~~~~~
  TYPE(DynamicsState_F) :: dynamics    ! Move to declaration section
  TYPE(DensityState_F) ::  densities
  TYPE(WindsState_F) :: winds
  TYPE(GasesState_F) :: gases
  TYPE(EphemerisState_F) :: ephem
  TYPE(PerturbationState_F) :: perts

  call getPosition_B(body, position)         ! contains position related metrics
  call getDynamicsState_B(body, dynamics)    ! atmosphere model values
  call getDensityState_B(body, densities)    ! atmosphere density values
  call getWindsState_B(body, winds)          ! atmosphere wind values
  call getGasesState_B(body, gases)          ! atmosphere constituent gas values
  call getEphemerisState_B(body, ephem)      ! planetary position values
  call getPerturbationState_B(body, perts)   ! (optional) random numbers used in perturbations
~~~~~~~~~~
The call to getGasesState_B may be diffenent for each model.  See the model specific documentation for the precise call signature.

Process the information.  Update the position and elapsed time.  And repeat.

- - - - - - - - - - - 	
### On Shutdown

The createAtmosphere_B() function allocated memory which needs to be freed.  Once deleteAtmosphere_B() has been called, the body handle becomes invalid.  It is a good practice to nullify the handle.
~~~~~~~~~~
  call deleteAtmosphere_B(body);
  body = 0;
~~~~~~~~~~
