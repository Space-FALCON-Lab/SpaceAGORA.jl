The MarsGRAM C++ Interface
=============================

This documents the \ref GRAM::MarsAtmosphere interface class for MarsGRAM.  Since the GRAM framework is common to all models, only the Mars specific information is documented on this page.  See the [C++ interface](Documentation/Markdown/GRAM.md) page for documentation on the common framework.

A simple example of using the Mars atmosphere class can be seen [here](@ref Mars/examples/Mars.cpp).

- - - - - - - 

Path to Mars Data
-----------------

The location of the Mars data files is specified with the method \ref GRAM::MarsAtmosphere::setDataPath.  

~~~~~~~~~~
mars.setDataPath("\path\tomy\marsdata\");
~~~~~~~~~~

Min Max Conmputations
---------------------

The method \ref GRAM::MarsAtmosphere::setMinMax flags the computation of daily minimum and maximum values.  

~~~~~~~~~~
mars.setMinMax(true);
~~~~~~~~~~

Planetary Radii
---------------

The default equatorial and polar radii can by overridden with the method \ref GRAM::MarsAtmosphere::setPlanetaryRadii.  

~~~~~~~~~~
mars.setPlanetaryRadii(3397.0, 3375.0);
~~~~~~~~~~

Model and Map Year
------------------

Use the method \ref GRAM::MarsAtmosphere::setMapYear to select the MGCM model (0) or the TES model year (1 or 2).  

~~~~~~~~~~
mars.setMapYear(2);
~~~~~~~~~~

The Height Offset Model
-----------------------

The height offset model is chosen with the method \ref GRAM::MarsAtmosphere::setHeightOffsetModel.  The models are enumerated with GRAM::MarsOffsetModel.  The \p MARS_CONSTANT and \p MARS_SEASONAL models require an additional height offset arguement.

~~~~~~~~~~
mars.setHeightOffsetModel(MARS_CONSTANT, 0.45);

// or

mars.setHeightOffsetModel(MARS_GLOBAL_MEAN);
~~~~~~~~~~

MGCM Dust Levels
---------------

Within the MGCM model, dust optical depths may be controlled using the method \ref GRAM::MarsAtmosphere::setMGCMDustLevels.  
~~~~~~~~~~
mars.setMGCMDustLevels(0.3, 0.3, 1.0);
~~~~~~~~~~

Solar Flux
----------

Use the method \ref GRAM::MarsAtmosphere::setF107 to set the 10.7 cm solar flux.

~~~~~~~~~~
mars.setF107(68.5);
~~~~~~~~~~

Perturbation Wave Length Scale
------------------------------

The scale factor for perturbation wave lengths can by set with the method \ref GRAM::MarsAtmosphere::setPerturbationWaveLengthScale.  

~~~~~~~~~~
mars.setPerturbationWaveLengthScale(2.0);
~~~~~~~~~~

MOLA Heights
------------

The method \ref GRAM::MarsAtmosphere::setMOLAHeights is used to designate that heights are relative to the MOLA topographic heights (true) or the reference ellipsoid (false).  

~~~~~~~~~~
mars.setMOLAHeights(true);
~~~~~~~~~~

Dust Storm Parameters
---------------------

A dust storm can be simulated with the method \ref GRAM::MarsAtmosphere::setDustStorm.  

~~~~~~~~~~
mars.setDustStorm(90.0, 5.0, 2.0, 1000.0, 35.5, 140.0);
~~~~~~~~~~

Dust Density Parameters
-----------------------

The parameters of the method \ref GRAM::MarsAtmosphere::setDustDensity are used to control the dust density output.  

~~~~~~~~~~
mars.setDustDensity(0.003, 4.0, 2500.0);
~~~~~~~~~~

Exospheric Temperature Parameters
---------------------------------

Use the method \ref GRAM::MarsAtmosphere::setExosphericTemperature to set the exospheric temperature parameters for the Stewart model 

~~~~~~~~~~
mars.setExosphericTemperature(50.0, -2.0);
~~~~~~~~~~

Wave Perturbation Parameters
----------------------------

Set the default values for the wave perturbation parameters with the method \ref GRAM::MarsAtmosphere::setWaveDefaults.  

~~~~~~~~~~
mars.setWaveDefaults(date, scale, mean, a1, p1, r1, a2, p2, r2, a3, p3, r3);
~~~~~~~~~~

Wave Perturbation File Input
----------------------------

Set the wave perturbation parameter file with the method \ref GRAM::MarsAtmosphere::setWaveFile.  

~~~~~~~~~~
mars.setWaveFile("/path/tomy/waveFile.txt");
~~~~~~~~~~

Wind Scales
-----------

The method \ref GRAM::MarsAtmosphere::setWindScales will set the wind scale factors.  

~~~~~~~~~~
mars.setWindScales(1.5, 3.0);
~~~~~~~~~~

