The GRAM Python Interface
====================

This documents the Python interface for the GRAM framework.  Since all GRAM models are built around this framework, the contents are apropos for each model.  When refering to a model specific function or item, simply replace <em>Body</em> with the model of choice (e.g. NeptuneAtmosphere for <em>Body</em>Atmosphere). For details on a specific class or method, please refer to the class documentation. 

The provided Python interface simply wraps GRAM's core 
C++ code, allowing the user to call C++ functions from a program written in Python. All exposed classes and methods should have the same behavior as they would in C++, meaning that Python programmers can reference any of the existing C++ documentation in addition to Python-specific documentation. 

- - - - - - - - - - - 	
The Interface Classes
---------------------
- - - - - - - - - - - 	

The Python interface provides wrappers for the following classes:
- Atmosphere Model Classes
	+ <em>Body</em>Atmosphere
		- VenusAtmosphere
		- EarthAtmosphere
		- MarsAtmosphere
		- JupiterAtmosphere
		- TitanAtmosphere
		- UranusAtmosphere
		- NeptuneAtmosphere	
- Data Classes
	+ <em>Body</em>InputParameters
	+ GramTime
	+ Position
	+ AtmosphereState
	+ EphemerisState
	+ Perturbation State
	+ ConstituentGas
- Utility Classes
	+ SpiceLoader
	+ <em>Body</em>NamelistReader

The classes above are wrappers around the C++ classes.  Please refer to the C++ class documentation for details about the member functions.

- - - - - - - - - - - 	
Usage of the Python Interface
------------------------
- - - - - - - - - - - 	

The instructions below assume that you have already built the Python interface on your computer following the instructions in **Python/README.md**.

A simple example using the Python interface can be seen [here](@ref Python/JupiterExample.py).

Python installation via package managers (pip/conda/etc.) is not currently supported. To use the `gram` module, users must import it using a relative import in their programs. 

The simplest way to do this would be to copy the entire `GRAMpy` directory into the folder where you are writing your Python program, then import `gram` like this in your program: 
~~~~~~~~~~{.py}
from GRAMpy import gram
~~~~~~~~~~

GRAM classes (such as JupiterAtmosphere) can now be accessed like this:
~~~~~~~~~~{.py}
jupiter = gram.JupiterAtmosphere()
~~~~~~~~~~

Advanced users who wish to set their import up differently certainly can- just note that  `_gram.so` must be in the same directory as `gram.py`. 

### Setup - Call these functions once

*First and Foremost:* The GRAM models use the NAIF SPICE library for time and ephemeris computations.  The location of the Spice data must be set before any other call into the GRAM models.  

One can optionally override the default SPICE kernels.  This should be done before initializing.
Use the <em>Body</em>InputParameters class. 
~~~~~~~~~~{.py}
  # Set Spice path
  inputParameters = gram.BodyInputParameters()
  reader = gram.BodyNamelistReader()
  reader.tryGetSpicePath(inputParameters)

  # Create an atmosphere model and set input parameters
  body = gram.BodyAtmosphere()
  body.setInputParameters(inputParameters)
~~~~~~~~~~

The atmosphere model needs to be initialized.  First, set the start time.  In GRAM models, time is computed using the number of elapsed seconds past the start time.
~~~~~~~~~~{.py}
  time = gram.GramTime()
  time.setStartTime(2020, 3, 15, 0, 0, 0, gram.UTC, gram.ERT)
  body.setStartTime(time)
~~~~~~~~~~

Choose a seed and set parameters for the random perturbations.
~~~~~~~~~~{.py}
  body.setSeed(1234);
  body.setPerturbationScaleFactor(1.5);
  body.setMinRelativeStepSize(0.5);
~~~~~~~~~~

Some models may require further setup.  See the model specific documentation for details.

- - - - - - - - - - - 	
### In the loop - Call these methods for each time step of the simulation.

Set the position and elapsed time.
~~~~~~~~~~{.py}
  position = gram.Position()
  position.height = 50.0       # km
  position.latitude = 22.0     # degrees
  position.longitude = 48.0    # degrees, east positive
  position.elapsedTime = 100.0 # seconds
  position.isGeocentric = 1    # 1 is for geocentric, set to 0 for geodetic
  body.setPosition(position)
~~~~~~~~~~

Perform the model computations with an update.
~~~~~~~~~~{.py}
  body.update()
~~~~~~~~~~

Retrieve the updated values
~~~~~~~~~~{.py}
 position = body.getPosition() # contains position related metrics
 atmos = body.getAtmosphereState()    # atmosphere model values
 ephem = body.getEphemerisState()      # planetary position values
 pert = body.getPerturbationState() # (optional) random numbers used in perturbations
~~~~~~~~~~

Process the information.  Update the position and elapsed time.  And repeat.

- - - - - - - - - - - 	

The **gramExample.ipynb** Jupyter notebook contains several more involved examples. Those familiar with using a Jupyter notebook can follow the instructions in **Python/README.md** to run them interactively. 

\example Python/JupiterExample.py
