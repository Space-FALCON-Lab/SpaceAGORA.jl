How to Create a New Atmosphere Model
====================================

This documents the method suggested for creating a new atmosphere model which utilizes the common framework.

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	
The Model Creation Process (Quick Links)
========================================
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	

1. [Identify the data model.](#c0)
2. [Duplicate an existing atmosphere code base with a similar data model.](#c1)
3. [Identify and modify model specific code and data.](#c2)
4. [Test the new model.](#c3)

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	
Identify the Data Model   {#c0}
=======================
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	

The HeightModel is the simplest data model based on a single height based profile of temperature, pressure, and density.  This model is used by JupiterGRAM and UranusGRAM.

The MinMaxModel is a two-dimensional model based on three height based profiles of temperature, pressure, and density.  A minMaxFactor must be formulated to appropriately interpolate between the min, mean, and max profiles.  This model is used by NeptuneGRAM and the TitanGRAM Yelle model.

MarsGRAM, VenusGRAM, and the TitanGRAM GCM models are not part of the common framework.  (They will be described in greater detail here eventually.)

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	
Duplicating a Code Base  {#c1}
=======================
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	

In the GRAM folder, copy the folder of the planet with the desired data model.  Rename the folder to the new planet.  Rename each of the source and header files to use the new planet name.

In multi-file editer (like that in MSVS), perform some multi-file edits to replace the old planet name with the new planet name.  It is best to do this with three case sensitive replaces (replace "Jupiter" with "Saturn", replace "jupiter" with "saturn", and replace "JUPITER" with "SATURN").

Duplicate and edit the build systems in GRAM/Build and GRAM/MSVS.  Test the copy/rename process by building the examples.

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	
Identify and Modify Model Specific Code and Data  {#c2}
================================================
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	

There are parts of each model that are peculiar to that model alone.  This section will list some of the functions that are model specific, and will help identify others.

### Parameters

Update the common parameters class named BodyCommon (replace "Body" with the appropriate body name).  These include planetary constants, gas molecular weights, and a specific heat ratio.  The constructor of this class should set the planetary constants.  It should also set the molecular weight for each gas of concern.  Gases will be ignored if the molecular weight is not set.

### BodyAtmosphere

Next, update code in the BodyAtmosphere class.  This class is the primary interface to the user for this model.  All parameters required of the user should have a "set" method in this class (and possibly a "get" method).  The update() method should be fine, but peruse it nonetheless.  The three perturbation methods will need to be updated as appropriate for the model.  These are the getPerturbationFactors(), getScaleParameters(), and getWindDeviations() methods. If no perturbation model is known, set the perturbation factors and scale parameters to one and the wind deviations to zero. Finally, make certain that the ephemeris body name was properly changed (search for setBody).

### Models

The BodyAtmosphere class will own one or more atmosphere model classes.  Each of these classes will need to be vigourously updated.  Many of the model classes are coded in two .cpp files, an implementation file and a data file.

The data for a model is typically found in a file ending in Data.cpp (such as NeptuneMinMaxData.cpp).  Most good editors have a column editing mode which is probably needed to accomplish this task in a timely manner.  Sizes may need to be modified in a header file and/or the initialData() function.  If the nature of the new data is different (e.g. units), one may need to modify the initializeData() function to standardize the model data.

In the implementation file, be sure to look at each method for model specific code.  It is likely that all code, except for the constructors, in this file is model specific.  In particular, the updateWinds() method is typically coded in the model class.  If no wind model is known, set the winds to zero.  

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	
Test the New Model  {#c3}
==================
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 	

The unit tests should have been copied over from the source project.  While these test will fail, they are a good starting point for creating new tests.  If a data model from the common framework is being used (HeightModel or MinMaxModel), then unit tests already exist to verify the framework. In some of the model specific unit tests, data initialization is verified by chosing a height that appears in the input data.  Since no interpolation is required, the density, temperature, and pressure ought to update to the values shown in the input data. For model specific functions, such as getPerturbationFactors(), the test values can be independently verified using a tool such as Excel or Matlab.

One simple test of the conversion process is to build the new BodyGRAM program and run a vertical profile from a namelist file.  Make sure the heights extend over the entire data set.  If the outputs appear to be appropriate, then the program is functional.  Proceed to system level testing.
