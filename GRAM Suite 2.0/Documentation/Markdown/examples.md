Example Programs
================

Quick Links
===========

- [Examples for Documentation](@ref ExampleDocs)
- [The C++ Example Suite](@ref ExampleCpp)

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

The GRAM suite is a collection of static libraries and programs.  This section describes the example programs that are provided in the distribution. Similar examples are provided for each body atmosphere.


- - - - - - - - - - - 	
Examples for Documentation   {#ExampleDocs}
==========================

Each body atmosphere is distributed with three examples used in model documentation of the C++, C, and Fortran interfaces.  The three examples produce equivalent results, differing only in the choice of interface.

- BodyExample_Cpp.exe is the executable for GRAM/Body/examples/Body.cpp.
- BodyExample_C.exe is the executable for GRAM/Body/examples/Body_C.c.
- BodyExample_F.exe is the executable for GRAM/Body/examples/Body_F.f90.

Each of these examples is a simple demonstration of using the interface to produce atmospheric values.  Control logic and looping are purposely omitted to avoid obfuscating the examples.  The examples show the following concepts:

- Reading and setting (via input parameters) the Spice data path.
- Creating multiple (2) atmosphere models.
- Getting and printing the version string.
- Setting model specific parameters (if applicable).
- Setting perturbation parameters.
- Setting the start time (epoch).
- Setting a position.
- Updating the model.
- Getting and printing various metrics.
- Freeing memory (C and Fortran).

For the developer that is building a simulation, the steps of setting the position, updating the model, and processing the metrics would need to be appropriately placed within the simulation logic.

- - - - - - - - - - - 	
The C++ Example Suite   {#ExampleCpp}
=====================

Each body atmosphere is distributed with a program, BodyExamples.exe, that will run a number of C++ interface examples.  This section describes the intent of each example.

The source code for the C++ examples can be found in GRAM/common/examples.  The main program which calls these examples is GRAM/Body/examples/ExampleRunner.cpp.


### Unit Testing  
[common/examples/UnitTestRunner.cpp](@ref common/examples/UnitTestRunner.cpp)

This is an extremely simple example of how to use the GoogleTest framework.  Note that the Spice data path must be correctly set before running the unit tests.  The Spice data path is set in ExampleRunner.cpp.

### Atmosphere example
[common/examples/AtmosphereExample.cpp](@ref common/examples/AtmosphereExample.cpp)

The Atmosphere example is the "Hello World" of the GRAM examples.  A BodyAtmosphere object is instantiated in ExampleRunner.cpp.  This object is passed to the Atmosphere example.  The start time and a position are set.  The model is updated.  And a few metrics are retrieved and printed.

### Trajectory example
[common/examples/TrajectoryExample.cpp](@ref common/examples/TrajectoryExample.cpp)

This example demonstrates the use of a Trajectory object to generate a profile.  In this example, the following occurs:

- A Trajectory object is created.
- The PerturbedAtmosphere object is initialized and given to the Trajectory.
    - The Trajectory object uses the PerturbedAtmosphere to generate metrics at each point of the trajectory.
- The trajectory is defined.
    - The initial position is set.
	- A constant change in position is defined.
	- The number of points to generate is set.
- The trajectory is generated.  Within this function, the following occurs:
    - The positions are automatically generated from the initial position and delta.
	- The atmosphere is updated at each position.
	- The results are stored in the Trajectory object.
- A ProfilePrinter is created and output styles are chosen.
- The data is sent to the ProfilePrinter for output.

### Monte Carlo example
[common/examples/MonteCarloExample.cpp](@ref common/examples/MonteCarloExample.cpp)

Building on the Trajectory example, the use of the MonteCarlo object is demonstrated.  The MonteCarlo object will run multiple profiles using different random seeds for the density and wind perturbations.  Usage of the MonteCarlo object requires a Profile object, such as a trajectory, a ProfilePrinter object, and a PerturbedAtmosphere object.  This example demonstrates how the objects are initialized and how the Monte Carlo run is generated.

### Namelist example
[common/examples/NamelistExample.cpp](@ref common/examples/NamelistExample.cpp)

The Namelist example demonstrates how to generate a Monte Carlo using a namelist file for inputs.  This example closely resembles the actual BodyGRAM program.  Since some input parameters are body specific, the reading of the namelist file occurs in the ExampleRunner.cpp file using a BodyNamelistReader object.  The BodyInputParameters are passed into the Namelist example function via the BodyAtmosphere object.

The Namelist example is very similar to the MonteCarlo example.  The primary difference is that the objects are initialized via the InputParameters object as opposed to the explicit initialization in the Monte Carlo example.

### Ephemeris example
[common/examples/EphemerisExample.cpp](@ref common/examples/EphemerisExample.cpp)

This example demonstrates the usage of the EphemerisProfile class.  The example generates a profile similar to the Trajectory example.  However, the actual ephemeris values are overridden by auto-generated ephemerides.  This allows one to observe the behaviour of a metric in response to a changing ephemeride.  This particular example alters longitude of the sun from 0 to 360 by 10 degree increments in order to observe the seasonal effects.  

### Perturbation example
[common/examples/PerturbationExample.cpp](@ref common/examples/PerturbationExample.cpp)

The use of the PerturbationState is demonstrated in the Perturbation example.  The first portion of this example generates atmospheric metrics for five positions.  The example saves the positions, the atmospheric states, and the perturbation states for these five points.

In the second portion of the example, the five positions and perturbation states are used to regenerate the atmosphere states.  The perturbed values of the two runs are compared in the output.  They should be identical.  In this second run, the user supplied perturbation states override the random number generators with the values computed in the first run.

This method of overriding the PerturbationState with saved states is useful for developers that need greater control on a trajectory sequence and its effects on the perturbations.


 
