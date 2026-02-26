GRAM Units - Standard Units used in Interfaces
==============================================

This document details the standardized units used in all GRAM programming interfaces.  This does not include print or file I/O which may be customized.

- All longitudes use an east positive convention.
- All percentages are decimal percentages between 0 and 1.

| Variable                   | Units           | Classes/Structs  |
|----------------------------|:---------------:|------------------|
| height                     | \f$ km      \f$ | Position |
| latitude                   | \f$ deg     \f$ | Position |
| longitude                  | \f$ deg     \f$ | Position |
| elapsedTime                | \f$ s       \f$ | Position |
| totalRadius                | \f$ km      \f$ | Position |
| latitudeRadius             | \f$ km      \f$ | Position |
| surfaceHeight              | \f$ km      \f$ | Position |
| gravity                    | \f$ m/s^2   \f$ | Position |
| solarTime                  | \f$ hr      \f$ | EphemerisState |
| longitudeSun               | \f$ deg     \f$ | EphemerisState |
| subsolarLatitude           | \f$ deg     \f$ | EphemerisState |
| subsolarLongitude          | \f$ deg     \f$ | EphemerisState |
| orbitalRadius              | \f$ AU      \f$ | EphemerisState |
| oneWayLightTime            | \f$ min     \f$ | EphemerisState |
| solarZenithAngle           | \f$ deg     \f$ | EphemerisState |
| secondsPerSol              | \f$ s       \f$ | EphemerisState |
| temperature                | \f$ K       \f$ | AtmosphereState, DynamicsState |
| pressure                   | \f$ Pa      \f$ | AtmosphereState, DynamicsState |
| density                    | \f$ kg/m^3  \f$ | AtmosphereState, DynamicsState |
| pressureScaleHeight        | \f$ km      \f$ | AtmosphereState, DynamicsState |
| densityScaleHeight         | \f$ km      \f$ | AtmosphereState, DynamicsState |
| speedOfSound               | \f$ m/s     \f$ | AtmosphereState, DynamicsState |
| referenceDensity           | \f$ kg/m^3  \f$ | AtmosphereState, DynamicsState |
| referenceTemperature       | \f$ K       \f$ | AtmosphereState, DynamicsState |
| referencePressure          | \f$ Pa      \f$ | AtmosphereState, DynamicsState |
| sigmaLevel                 | \f$ --      \f$ | AtmosphereState, DynamicsState |
| pressureAtSurface          | \f$ Pa      \f$ | AtmosphereState, DynamicsState |
| pressureAltitude           | \f$ km      \f$ | AtmosphereState, DynamicsState |
| lowDensity                 | \f$ kg/m^3  \f$ | AtmosphereState, DensityState |
| highDensity                | \f$ kg/m^3  \f$ | AtmosphereState, DensityState |
| densityPerturbationPercent | \f$ \%      \f$ | AtmosphereState, DensityState |
| perturbedDensity           | \f$ kg/m^3  \f$ | AtmosphereState, DensityState |
| densityStandardDeviation   | \f$ kg/m^3  \f$ | AtmosphereState, DensityState |
| perturbedSpeedOfSound      | \f$ m/s     \f$ | AtmosphereState, DensityState |
| relativeStepSize           | \f$ --      \f$ | AtmosphereState, DensityState |
| densityDeviation           | \f$ \%      \f$ | AtmosphereState, DensityState |
| lowDensityDeviation        | \f$ \%      \f$ | AtmosphereState, DensityState |
| highDensityDeviation       | \f$ \%      \f$ | AtmosphereState, DensityState |
| perturbedDensityDeviation  | \f$ \%      \f$ | AtmosphereState, DensityState |
| ewWind                     | \f$ m/s     \f$ | AtmosphereState, WindsState |
| nsWind                     | \f$ m/s     \f$ | AtmosphereState, WindsState |
| verticalWind               | \f$ m/s     \f$ | AtmosphereState, WindsState |
| ewWindPerturbation         | \f$ m/s     \f$ | AtmosphereState, WindsState |
| nsWindPerturbation         | \f$ m/s     \f$ | AtmosphereState, WindsState |
| verticalWindPerturbation   | \f$ m/s     \f$ | AtmosphereState, WindsState |
| perturbedEWWind            | \f$ m/s     \f$ | AtmosphereState, WindsState |
| perturbedNSWind            | \f$ m/s     \f$ | AtmosphereState, WindsState |
| perturbedVerticalWind      | \f$ m/s     \f$ | AtmosphereState, WindsState |
| ewStandardDeviation        | \f$ m/s     \f$ | AtmosphereState, WindsState |
| nsStandardDeviation        | \f$ m/s     \f$ | AtmosphereState, WindsState |
| verticalStandardDeviation  | \f$ m/s     \f$ | AtmosphereState, WindsState |
| totalNumberDensity         | \f$ \#/m^3  \f$ | AtmosphereState, GasesState |
| averageMolecularWeight     | \f$ --      \f$ | AtmosphereState, GasesState |
| compressibilityFactor      | \f$ --      \f$ | AtmosphereState, GasesState |
| specificGasConstant        | \f$ J/(kg K)\f$ | AtmosphereState, GasesState |
| numberDensity              | \f$ \#/m^3  \f$ | AtmosphereState, ConstituentGas |
| moleFraction               | \f$ \%      \f$ | AtmosphereState, ConstituentGas |
| massFraction               | \f$ \%      \f$ | AtmosphereState, ConstituentGas |
| avgMolecularWeight         | \f$ --      \f$ | AtmosphereState, ConstituentGas |

