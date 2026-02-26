The MarsGRAM C Interface
===========================

This documents the C interface for MarsGRAM.  Since the GRAM framework is common to all models, only the Mars specific information is documented on this page.  See the [C interface](Documentation/Markdown/GRAM_C.md) page for documentation on the common framework.

A simple example of using the Mars atmosphere can be seen [here](@ref Mars/examples/Mars_C.c).

- - - - - - - 

Path to Mars Data
-----------------

The location of the Mars data files are specified during the creation of the \p MarsAtmosphere_C handle with the function \ref createAtmosphere_M.  

~~~~~~~~~~
MarsAtmosphere_C* mars = createAtmosphere_M("\path\tomy\marsdata\");
~~~~~~~~~~

Min Max Conmputations
---------------------

The function \ref setMinMax_M flags the computation of daily minimum and maximum values.  

~~~~~~~~~~
setMinMax_M(mars, 1);
~~~~~~~~~~

Planetary Radii
---------------

The default equatorial and polar radii can by overridden with the function \ref setPlanetaryRadii_M.  

~~~~~~~~~~
setPlanetaryRadii_M(mars, 3397.0, 3375.0);
~~~~~~~~~~

Model and Map Year
------------------

Use the function \ref setMapYear_M to select the MGCM model (0) or the TES model year (1 or 2).  

~~~~~~~~~~
setMapYear_M(mars, 2);
~~~~~~~~~~

The Height Offset Model
-----------------------

The height offset model is chosen with the function \ref setHeightOffsetModel_M.  The \p MARS_CONSTANT (0) and \p MARS_SEASONAL (1) models require an additional height offset arguement. For other models, this parameter is ignored.

~~~~~~~~~~
setHeightOffsetModel_M(mars, 0, 0.45);

// or

setHeightOffsetModel_M(mars, 2, 0.0);
~~~~~~~~~~

MGCM Dust Levels
---------------

Within the MGCM model, dust optical depths may be controlled using the function \ref setMGCMDustLevels_M.  
~~~~~~~~~~
setMGCMDustLevels_M(mars, 0.3, 0.3, 1.0);
~~~~~~~~~~

Solar Flux
----------

Use the function \ref setF107_M to set the 10.7 cm solar flux.

~~~~~~~~~~
setF107_M(mars, 68.5);
~~~~~~~~~~

Perturbation Wave Length Scale
------------------------------

The scale factor for perturbation wave lengths can by set with the function \ref setPerturbationWaveLengthScale_M.  

~~~~~~~~~~
setPerturbationWaveLengthScale_M(mars, 2.0);
~~~~~~~~~~

MOLA Heights
------------

The function \ref setMOLAHeights_M is used to designate that heights are relative to the MOLA topographic heights (true) or the reference ellipsoid (false).  

~~~~~~~~~~
setMOLAHeights_M(mars, 1);
~~~~~~~~~~

Dust Storm Parameters
---------------------

A dust storm can be simulated with the function \ref setDustStorm_M.  

~~~~~~~~~~
setDustStorm_M(mars, 90.0, 5.0, 2.0, 1000.0, 35.5, 140.0);
~~~~~~~~~~

Dust Density Parameters
-----------------------

The parameters of the function \ref setDustDensity_M are used to control the dust density output.  

~~~~~~~~~~
setDustDensity_M(mars, 0.003, 4.0, 2500.0);
~~~~~~~~~~

Exospheric Temperature Parameters
---------------------------------

Use the function \ref setExosphericTemperature_M to set the exospheric temperature parameters for the Stewart model 

~~~~~~~~~~
setExosphericTemperature_M(mars, 50.0, -2.0);
~~~~~~~~~~

Wave Perturbation Parameters
----------------------------

Set the default values for the wave perturbation parameters with the function \ref setWaveDefaults_M.  

~~~~~~~~~~
setWaveDefaults_M(mars, date, scale, mean, a1, p1, r1, a2, p2, r2, a3, p3, r3);
~~~~~~~~~~

Wave Perturbation File Input
----------------------------

Set the wave perturbation parameter file with the function \ref setWaveFile_M.  

~~~~~~~~~~
setWaveFile_M(mars, "/path/tomy/waveFile.txt");
~~~~~~~~~~

Wind Scales
-----------

The function \ref setWindScales_M will set the wind scale factors.  

~~~~~~~~~~
setWindScales_M(mars, 1.5, 3.0);
~~~~~~~~~~

