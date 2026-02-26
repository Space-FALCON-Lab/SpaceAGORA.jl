Variable Name Mapping (Legacy versus New GRAM)
==============================================

This document contains a list of new variable names corresponding to the legacy FORTRAN names.  Since the legacy GRAMs did not always follow a naming convention, this list may not be exhaustive.  The units of the legacy and new GRAMs may differ.  EarthGRAM names tend to differ from the other GRAMs and are presented in [this table](@ref vnm).

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Variable Name Mapping  
=====================

| Legacy Variable Name | New Name                   | Classes          |
|----------------------|----------------------------|------------------|
| chgt                 | height                     | Position |
| clat                 | latitude                   | Position |
| clon                 | longitude                  | Position |
| csec                 | elapsedTime                | Position, GramTime |
| radtotal, rsc        | totalRadius                | Position |
| rref                 | latitudeRadius             | Position |
| ----                 | surfaceHeight              | Position |
| gz                   | gravity                    | Position |
| lst, tlocal          | solarTime                  | EphemerisState |
| als                  | longitudeSun               | EphemerisState |
| sunlat               | subsolarLatitude           | EphemerisState |
| sunlon               | subsolarLongitude          | EphemerisState |
| (body)au             | orbitalRadius              | EphemerisState |
| owlt                 | oneWayLightTime            | EphemerisState, GramTime |
| sza                  | solarZenithAngle           | EphemerisState |
| onesol, day(body)    | secondsPerSol              | EphemerisState |
| temp                 | temperature                | AtmosphereState, DynamicsState |
| pres                 | pressure                   | AtmosphereState, DynamicsState |
| dens                 | density                    | AtmosphereState, DynamicsState |
| hscale               | pressureScaleHeight        | AtmosphereState, DynamicsState |
| hrho                 | densityScaleHeight         | AtmosphereState, DynamicsState |
| ssp2  (squared)      | speedOfSound               | AtmosphereState, DynamicsState |
| dnpn, dyel, dvra     | referenceDensity           | AtmosphereState, DynamicsState |
| tnpn, tyel, tvra     | referenceTemperature       | AtmosphereState, DynamicsState |
| pnpn, pyel, pvra     | referencePressure          | AtmosphereState, DynamicsState |
| sigmalevel           | sigmaLevel                 | AtmosphereState, DynamicsState |
| patsurf              | pressureAtSurface          | AtmosphereState, DynamicsState |
| preshgt              | pressureAltitude           | AtmosphereState, DynamicsState |
| denslo               | lowDensity                 | AtmosphereState, DensityState |
| denshi               | highDensity                | AtmosphereState, DensityState |
| densp                | densityPerturbation        | AtmosphereState, DensityState |
| denrand              | densityPerturbation        | AtmosphereState, DensityState |
| denstot              | perturbedDensity           | AtmosphereState, DensityState |
| sigd                 | densityStandardDeviation   | AtmosphereState, DensityState |
| ----                 | perturbedSpeedOfSound      | AtmosphereState, DensityState |
| corlim               | relativeStepSize           | AtmosphereState, DensityState |
| devav                | densityDeviation           | AtmosphereState, DensityState |
| devlo                | lowDensityDeviation        | AtmosphereState, DensityState |
| devhi                | highDensityDeviation       | AtmosphereState, DensityState |
| devtot               | perturbedDensityDeviation  | AtmosphereState, DensityState |
| ewwind, u            | ewWind                     | AtmosphereState, WindsState |
| nswind, v            | nsWind                     | AtmosphereState, WindsState |
| ----                 | verticalWind               | AtmosphereState, WindsState |
| ewpert               | ewWindPerturbation         | AtmosphereState, WindsState |
| nspert               | nsWindPerturbation         | AtmosphereState, WindsState |
| vwpert               | verticalWindPerturbation   | AtmosphereState, WindsState |
| ewtot                | perturbedEWWind            | AtmosphereState, WindsState |
| nstot                | perturbedNSWind            | AtmosphereState, WindsState |
| ----                 | perturbedVerticalWind      | AtmosphereState, WindsState |
| sigu                 | ewStandardDeviation        | AtmosphereState, WindsState |
| sigv                 | nsStandardDeviation        | AtmosphereState, WindsState |
| sigw                 | verticalStandardDeviation  | AtmosphereState, WindsState |
| ytotnd, totnd        | totalNumberDensity         | AtmosphereState, GasesState |
| amz                  | averageMolecularWeight     | AtmosphereState, GasesState |
| zeta                 | compressibilityFactor      | AtmosphereState, GasesState |
| rgas, r              | specificGasConstant        | AtmosphereState, GasesState |
| y(gas)nd             | numberDensity              | AtmosphereState, ConstituentGas |
| fmol(gas)            | moleFraction               | AtmosphereState, ConstituentGas |
| y(gas)p              | massFraction               | AtmosphereState, ConstituentGas |
| am(gas)              | avgMolecularWeight         | AtmosphereState, ConstituentGas |
| fminmax              | minMaxFactor               | AtmosphereState |
| ifmm                 | computeMinMaxFactor        | MinMaxModel |
| profwgt              | profileWeight              | AuxiliaryAtmosphere |
| profile              | auxFileName                | AuxiliaryAtmosphere |
| profnear             | innerRadius                | AuxiliaryAtmosphere |
| proffar              | outerRadius                | AuxiliaryAtmosphere |
| loneast              | eastLongitudePositive      | AuxiliaryAtmosphere, ProfilePrinter |
| nr1                  | initialSeed                | MonteCarlo |
| nmonte               | sampleSize                 | MonteCarlo |
| trajfl               | trajectoryFileName         | Trajectory |
| fhgt,flat,flon,fsec  | initialPostion             | Trajectory |
| delhgt,dellat,dellon,delsec  | deltaPostion       | Trajectory |
| npos                 | numPoints                  | Trajectory |
| iutc                 | timeScale                  | GramTime |
| iert                 | timeFrame                  | GramTime |
| day0                 | spiceEphemerisTime (closest) | GramTime |
| z1 (reused local)    | densityRandomNumber        | PerturbationState |
| z1 (reused local)    | ewWindRandomNumber         | PerturbationState |
| z1 (reused local)    | nsWindRandomNumber         | PerturbationState |
| z1 (reused local)    | verticalWindRandomNumber   | PerturbationState |
| iupdate              | perturbationAction         | PerturbedAtmosphere |
| gm                   | mu                         | Atmosphere |
| requ                 | equatorialRadius           | Atmosphere |
| rpol                 | PolarRadius                | Atmosphere |
| J2                   | J2                         | Atmosphere |
| per                  | period                     | Atmosphere |
| ----                 | specificHeatRatio          | Atmosphere |
| logscale             | densityPrintScale          | ProfilePrinter |


EarthGRAM Variable Name Mapping   {#vnm}
===============================
The list below is not exhaustive.  Other variable changes are documented within the EarthGRAM code.

| Legacy Variable Name | New Name                   | Classes          |
|----------------------|----------------------------|------------------|
| ch, h                | height                     | Position |
| phi                  | latitude                   | Position |
| thet                 | longitude                  | Position |
| elt                  | elapsedTime                | Position, GramTime |
| ri                   | totalRadius                | Position |
| ----                 | latitudeRadius             | Position |
| hsrf1                | surfaceHeight              | Position |
| g                    | gravity                    | Position |
| tgh                  | temperature                | AtmosphereState, DynamicsState |
| pgh                  | pressure                   | AtmosphereState, DynamicsState |
| dgh                  | density                    | AtmosphereState, DynamicsState |
| hgtp                 | pressureScaleHeight        | AtmosphereState, DynamicsState |
| hgtd                 | densityScaleHeight         | AtmosphereState, DynamicsState |
| csp0                 | speedOfSound               | AtmosphereState, DynamicsState |
| d1                   | referenceDensity           | AtmosphereState, DynamicsState |
| t1                   | referenceTemperature       | AtmosphereState, DynamicsState |
| p1                   | referencePressure          | AtmosphereState, DynamicsState |
| ----                 | sigmaLevel                 | AtmosphereState, DynamicsState |
| psrf                 | pressureAtSurface          | AtmosphereState, DynamicsState |
| ----                 | pressureAltitude           | AtmosphereState, DynamicsState |
| ----                 | lowDensity                 | AtmosphereState, DensityState |
| ----                 | highDensity                | AtmosphereState, DensityState |
| drh                  | densityPerturbation        | AtmosphereState, DensityState |
| dh                   | perturbedDensity           | AtmosphereState, DensityState |
| sdh                  | densityStandardDeviation   | AtmosphereState, DensityState |
| sph                  | pressureStandardDeviation  | AtmosphereState, DensityState |
| sth                  | temperatureStandardDeviation | AtmosphereState, DensityState |
| csp                  | perturbedSpeedOfSound      | AtmosphereState, DensityState |
| ----                 | relativeStepSize           | AtmosphereState, DensityState |
| dghp                 | densityDeviation           | AtmosphereState, DensityState |
| ----                 | lowDensityDeviation        | AtmosphereState, DensityState |
| ----                 | highDensityDeviation       | AtmosphereState, DensityState |
| dhp                  | perturbedDensityDeviation  | AtmosphereState, DensityState |
| ugh                  | ewWind                     | AtmosphereState, WindsState |
| vgh                  | nsWind                     | AtmosphereState, WindsState |
| wgh                  | verticalWind               | AtmosphereState, WindsState |
| urh                  | ewWindPerturbation         | AtmosphereState, WindsState |
| vrh                  | nsWindPerturbation         | AtmosphereState, WindsState |
| wrh                  | verticalWindPerturbation   | AtmosphereState, WindsState |
| uh                   | perturbedEWWind            | AtmosphereState, WindsState |
| vh                   | perturbedNSWind            | AtmosphereState, WindsState |
| wh                   | perturbedVerticalWind      | AtmosphereState, WindsState |
| suh                  | ewStandardDeviation        | AtmosphereState, WindsState |
| svh                  | nsStandardDeviation        | AtmosphereState, WindsState |
| swh                  | verticalStandardDeviation  | AtmosphereState, WindsState |
| totnd                | totalNumberDensity         | AtmosphereState, GasesState |
| mwnd                 | averageMolecularWeight     | AtmosphereState, GasesState |
| ----                 | compressibilityFactor      | AtmosphereState, GasesState |
| gasconst             | specificGasConstant        | AtmosphereState, GasesState |
| (gas)nd              | numberDensity              | AtmosphereState, ConstituentGas |
| ----                 | moleFraction               | AtmosphereState, ConstituentGas |
| ----                 | massFraction               | AtmosphereState, ConstituentGas |
| wt(gas)              | avgMolecularWeight         | AtmosphereState, ConstituentGas |
| profwgt              | profileWeight              | AuxiliaryAtmosphere |
| profile              | auxFileName                | AuxiliaryAtmosphere |
| profnear, sitenear   | innerRadius                | AuxiliaryAtmosphere |
| proffar, sitelim     | outerRadius                | AuxiliaryAtmosphere |
| nr1                  | initialSeed                | MonteCarlo |
| mc                   | sampleSize                 | MonteCarlo |
| trapath              | trajectoryFileName         | Trajectory |
| h1, phi1, thet1      | initialPostion             | Trajectory |
| dz, dphi, dthet, delt | deltaPostion              | Trajectory |
| nmax                 | numPoints                  | Trajectory |
| iupdate              | perturbationAction         | PerturbedAtmosphere |
| prh                  | pressurePerturbation       | EarthAtmosphereState |
| ph                   | perturbedPressure          | EarthAtmosphereState |
| trh                  | temperaturePerturbation    | EarthAtmosphereState |
| th                   | perturbedTemperature       | EarthAtmosphereState |
| prhs                 | presPertSmall              | EarthAtmosphereState |
| sphs                 | presSDSmall                | EarthAtmosphereState |
| prhl                 | presPertLarge              | EarthAtmosphereState |
| sphl                 | presSDLarge                | EarthAtmosphereState |
| eofT                 | vaporPressure              | EarthAtmosphereState |
| seofT                | vaporPressureSD            | EarthAtmosphereState |
| drhs                 | densPertSmall              | EarthAtmosphereState |
| sdhs                 | densSDSmall                | EarthAtmosphereState |
| drhl                 | densPertLarge              | EarthAtmosphereState |
| sdhl                 | densSDLarge                | EarthAtmosphereState |
| rhov                 | vaporDensity               | EarthAtmosphereState |
| srhov                | vaporDensitySD             | EarthAtmosphereState |
| trhs                 | tempPertSmall              | EarthAtmosphereState |
| sths                 | tempSDSmall                | EarthAtmosphereState |
| trhl                 | tempPertLarge              | EarthAtmosphereState |
| sthl                 | tempSDLarge                | EarthAtmosphereState |
| tdgh                 | dewPoint                   | EarthAtmosphereState |
| stdgh                | dewPointSD                 | EarthAtmosphereState |
| urhs                 | ewWindPertSmall            | EarthAtmosphereState |
| suhs                 | ewWindSDSmall              | EarthAtmosphereState |
| urhl                 | ewWindPertLarge            | EarthAtmosphereState |
| suhl                 | ewWindSDLarge              | EarthAtmosphereState |
| vrhs                 | nsWindPertSmall            | EarthAtmosphereState |
| svhs                 | nsWindSDSmall              | EarthAtmosphereState |
| vrhl                 | nsWindPertLarge            | EarthAtmosphereState |
| svhl                 | nsWindSDLarge              | EarthAtmosphereState |
| rhp                  | relativeHumidity           | EarthAtmosphereState |
| srhp                 | relativeHumiditySD         | EarthAtmosphereState |
| dsrf                 | densityAtSurface           | EarthAtmosphereState |
| tsrf                 | temperatureAtSurface       | EarthAtmosphereState |
| usrf                 | ewWindAtSurface            | EarthAtmosphereState |
| vsrf                 | nsWindAtSurface            | EarthAtmosphereState |
| hsrf                 | surfaceHeight              | EarthAtmosphereState |
| spsrf                | pressureSDAtSurface        | EarthAtmosphereState |
| sdsrf                | densitySDAtSurface         | EarthAtmosphereState |
| stsrf                | temperatureSDAtSurface     | EarthAtmosphereState |
| susrf                | ewWindSDAtSurface          | EarthAtmosphereState |
| svsrf                | nsWindSDAtSurface          | EarthAtmosphereState |
| spdavsrf             | windSpeedAtSurface         | EarthAtmosphereState |
| spdsdsrf             | windSpeedSDAtSurface       | EarthAtmosphereState |
| uvtsrf               | windCorrelationAtSurface   | EarthAtmosphereState |
| isev                 | severityLevel              | EarthAtmosphereState |
| z0                   | surfaceRoughness           | EarthAtmosphereState |
| rday                 | solarDays                  | EarthAtmosphereState |
| el                   | solarElevation             | EarthAtmosphereState |
| elmn                 | elevationAtMidnight        | EarthAtmosphereState |
| elmd                 | elevationAtNoon            | EarthAtmosphereState |
| nri                  | netRadiationIndex          | EarthAtmosphereState |
| s                    | stability                  | EarthAtmosphereState |
| ool                  | inverseLength              | EarthAtmosphereState |
| ustar                | frictionVelocity           | EarthAtmosphereState |
| bvfsq                | BVFrequencySquare          | EarthAtmosphereState |
| hbl                  | boundaryLayerDepth         | EarthAtmosphereState |
| hn                   | neutralBoundaryLayerDepth  | EarthAtmosphereState |
| blfact               | unstableBLFactor           | EarthAtmosphereState |
| sigrat               | sigmaRatio                 | EarthAtmosphereState |
| swb                  | sigmaW                     | EarthAtmosphereState |
| chb                  | metersAboveSurface         | EarthAtmosphereState |
| spdsrf               | perturbedWindSpeedAtSurface | EarthAtmosphereState |