The MarsGRAM F Interface
===========================

This documents the FORTRAN interface for MarsGRAM.  Since the GRAM framework is common to all models, only the Mars specific information is documented on this page.  See the [FORTRAN interface](Documentation/Markdown/GRAM_F.md) page for documentation on the common framework.

A simple example of using the Mars atmosphere can be seen [here](@ref Mars/examples/Mars_F.f90).

- - - - - - - 

Path to Mars Data
-----------------

The location of the Mars data files are specified during the creation of the Mars atmosphere handle with the function \ref marsgram::createAtmosphere_M.

~~~~~~~~~~
TYPE(C_PTR) :: mars
mars = createAtmosphere_M("\path\tomy\marsdata\")
~~~~~~~~~~

Min Max Conmputations
---------------------

The function \ref marsgram::setMinMax_M flags the computation of daily minimum and maximum values.  

~~~~~~~~~~
call setMinMax_M(mars, 1)
~~~~~~~~~~

Planetary Radii
---------------

The default equatorial and polar radii can by overridden with the function \ref marsgram::setPlanetaryRadii_M.  

~~~~~~~~~~
call setPlanetaryRadii_M(mars, 3397.0d0, 3375.0d0)
~~~~~~~~~~

Model and Map Year
------------------

Use the function \ref marsgram::setMapYear_M to select the MGCM model (0) or the TES model year (1 or 2).  

~~~~~~~~~~
call setMapYear_M(mars, 2)
~~~~~~~~~~

The Height Offset Model
-----------------------

The height offset model is chosen with the function \ref marsgram::setHeightOffsetModel_M.  The \p MARS_CONSTANT (0) and \p MARS_SEASONAL (1) models require an additional height offset arguement. For other models, this parameter is ignored.

~~~~~~~~~~
call setHeightOffsetModel_M(mars, 0, 0.45d0)

! or

call setHeightOffsetModel_M(mars, 2, 0.0d0)
~~~~~~~~~~

MGCM Dust Levels
---------------

Within the MGCM model, dust optical depths may be controlled using the function \ref marsgram::setMGCMDustLevels_M.  

~~~~~~~~~~
call setMGCMDustLevels_M(mars, 0.3d0, 0.3d0, 1.0d0)
~~~~~~~~~~

Solar Flux
----------

Use the function \ref marsgram::setF107_M to set the 10.7 cm solar flux.

~~~~~~~~~~
call setF107_M(mars, 68.5d0)
~~~~~~~~~~

Perturbation Wave Length Scale
------------------------------

The scale factor for perturbation wave lengths can by set with the function \ref marsgram::setPerturbationWaveLengthScale_M.  

~~~~~~~~~~
call setPerturbationWaveLengthScale_M(mars, 2.0d0)
~~~~~~~~~~

MOLA Heights
------------

The function \ref marsgram::setMOLAHeights_M is used to designate that heights are relative to the MOLA topographic heights (true) or the reference ellipsoid (false).  

~~~~~~~~~~
call setMOLAHeights_M(mars, 1)
~~~~~~~~~~

Dust Storm Parameters
---------------------

A dust storm can be simulated with the function \ref marsgram::setDustStorm_M.  

~~~~~~~~~~
call setDustStorm_M(mars, 90.0d0, 5.0d0, 2.0d0, 1000.0d0, 35.5d0, 140.0d0)
~~~~~~~~~~

Dust Density Parameters
-----------------------

The parameters of the function \ref marsgram::setDustDensity_M are used to control the dust density output.  

~~~~~~~~~~
call setDustDensity_M(mars, 0.003d0, 4.0d0, 2500.0d0)
~~~~~~~~~~

Exospheric Temperature Parameters
---------------------------------

Use the function \ref marsgram::setExosphericTemperature_M to set the exospheric temperature parameters for the Stewart model 

~~~~~~~~~~
call setExosphericTemperature_M(mars, 50.0d0, -2.0d0)
~~~~~~~~~~

Wave Perturbation Parameters
----------------------------

Set the default values for the wave perturbation parameters with the function \ref marsgram::setWaveDefaults_M.  

~~~~~~~~~~
call setWaveDefaults_M(mars, date, scale, mean, a1, p1, r1, a2, p2, r2, a3, p3, r3)
~~~~~~~~~~

Wave Perturbation File Input
----------------------------

Set the wave perturbation parameter file with the function \ref marsgram::setWaveFile_M.  

~~~~~~~~~~
call setWaveFile_M(mars, "/path/tomy/waveFile.txt")
~~~~~~~~~~

Wind Scales
-----------

The function \ref marsgram::setWindScales_M will set the wind scale factors.  

~~~~~~~~~~
call setWindScales_M(mars, 1.5d0, 3.0d0)
~~~~~~~~~~

