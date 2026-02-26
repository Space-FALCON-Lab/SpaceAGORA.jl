//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"

struct Position_C
{
  greal height;          //!< \copydoc GRAM::Position::height
  greal latitude;        //!< \copydoc GRAM::Position::latitude
  greal longitude;       //!< \copydoc GRAM::Position::longitude
  greal elapsedTime;     //!< \copydoc GRAM::Position::elapsedTime
  int isPlanetoCentric;  //!< 1 if planeto-centric, 0 if planeto-graphic

  greal totalRadius;     //!< \copydoc GRAM::Position::totalRadius
  greal latitudeRadius;  //!< \copydoc GRAM::Position::latitudeRadius
  greal surfaceHeight;   //!< \copydoc GRAM::Position::surfaceHeight
  greal gravity;         //!< \copydoc GRAM::Position::gravity
};
typedef struct Position_C Position_C;

struct GramTime_C
{
  int year;        //!< A four digit year
  int month;       //!< 1 through 12
  int day;         //!< 1 through 31
  int hour;        //!< 0 through 23
  int minute;      //!< 0 through 59
  double seconds;  //!< 0 through 60, may contain fractions
  int scale;       //!< 0 = TDT, 1 = UTC, 2 = TDB
  int frame;       //!< 0 = PET, 1 = ERT
};
typedef struct GramTime_C GramTime_C;

struct DynamicsState_C
{
  greal temperature;              //!< \copydoc GRAM::AtmosphereState::temperature
  greal pressure;                 //!< \copydoc GRAM::AtmosphereState::pressure
  greal density;                  //!< \copydoc GRAM::AtmosphereState::density
  greal pressureScaleHeight;      //!< \copydoc GRAM::AtmosphereState::pressureScaleHeight
  greal densityScaleHeight;       //!< \copydoc GRAM::AtmosphereState::densityScaleHeight
  greal speedOfSound;             //!< \copydoc GRAM::AtmosphereState::speedOfSound
  greal referenceDensity;         //!< \copydoc GRAM::AtmosphereState::referenceDensity
  greal referenceTemperature;     //!< \copydoc GRAM::AtmosphereState::referenceTemperature
  greal referencePressure;        //!< \copydoc GRAM::AtmosphereState::referencePressure
  greal sigmaLevel;               //!< \copydoc GRAM::AtmosphereState::sigmaLevel
  greal pressureAtSurface;        //!< \copydoc GRAM::AtmosphereState::pressureAtSurface
  greal pressureAltitude;         //!< \copydoc GRAM::AtmosphereState::pressureAltitude
};
typedef struct DynamicsState_C DynamicsState_C;

struct DensityState_C
{
  greal density;                     //!< \copydoc GRAM::AtmosphereState::density
  greal lowDensity;                  //!< \copydoc GRAM::AtmosphereState::lowDensity
  greal highDensity;                 //!< \copydoc GRAM::AtmosphereState::highDensity
  greal densityPerturbation;         //!< \copydoc GRAM::AtmosphereState::densityPerturbation
  greal perturbedDensity;            //!< \copydoc GRAM::AtmosphereState::perturbedDensity
  greal densityStandardDeviation;    //!< \copydoc GRAM::AtmosphereState::densityStandardDeviation
  greal perturbedSpeedOfSound;       //!< \copydoc GRAM::AtmosphereState::perturbedSpeedOfSound
  greal relativeStepSize;            //!< \copydoc GRAM::AtmosphereState::relativeStepSize
  greal referenceDensity;            //!< \copydoc GRAM::AtmosphereState::referenceDensity
  greal densityDeviation;            //!< \copydoc GRAM::AtmosphereState::densityDeviation
  greal lowDensityDeviation;         //!< \copydoc GRAM::AtmosphereState::lowDensityDeviation
  greal highDensityDeviation;        //!< \copydoc GRAM::AtmosphereState::highDensityDeviation
  greal perturbedDensityDeviation;   //!< \copydoc GRAM::AtmosphereState::perturbedDensityDeviation
  int updateStatus;                  //!< \copydoc GRAM::AtmosphereState::updateStatus
};
typedef struct DensityState_C DensityState_C;

struct WindsState_C
{
  greal ewWind;                    //!< \copydoc GRAM::AtmosphereState::ewWind
  greal nsWind;                    //!< \copydoc GRAM::AtmosphereState::nsWind
  greal verticalWind;              //!< \copydoc GRAM::AtmosphereState::verticalWind
  greal ewWindPerturbation;        //!< \copydoc GRAM::AtmosphereState::ewWindPerturbation
  greal nsWindPerturbation;        //!< \copydoc GRAM::AtmosphereState::nsWindPerturbation
  greal verticalWindPerturbation;  //!< \copydoc GRAM::AtmosphereState::verticalWindPerturbation
  greal perturbedEWWind;           //!< \copydoc GRAM::AtmosphereState::perturbedEWWind
  greal perturbedNSWind;           //!< \copydoc GRAM::AtmosphereState::perturbedNSWind
  greal perturbedVerticalWind;     //!< \copydoc GRAM::AtmosphereState::perturbedVerticalWind
  greal ewStandardDeviation;       //!< \copydoc GRAM::AtmosphereState::ewStandardDeviation
  greal nsStandardDeviation;       //!< \copydoc GRAM::AtmosphereState::nsStandardDeviation
  greal verticalStandardDeviation; //!< \copydoc GRAM::AtmosphereState::verticalStandardDeviation
};
typedef struct WindsState_C WindsState_C;

struct GasesState_C
{
  greal totalNumberDensity;       //!< \copydoc GRAM::AtmosphereState::totalNumberDensity
  greal averageMolecularWeight;   //!< \copydoc GRAM::AtmosphereState::averageMolecularWeight
  greal compressibilityFactor;    //!< \copydoc GRAM::AtmosphereState::compressibilityFactor
  greal specificGasConstant;      //!< \copydoc GRAM::AtmosphereState::specificGasConstant
  greal specificHeatRatio;        //!< \copydoc GRAM::AtmosphereState::specificHeatRatio
};
typedef struct GasesState_C GasesState_C;

struct ConstituentGas_C
{
  greal moleFraction;             //!< \copydoc GRAM::ConstituentGas::moleFraction
  greal massFraction;             //!< \copydoc GRAM::ConstituentGas::massFraction
  greal numberDensity;            //!< \copydoc GRAM::ConstituentGas::numberDensity
  greal averageMolecularWeight;   //!< \copydoc GRAM::ConstituentGas::averageMolecularWeight
  greal specificHeatCapacity;     //!< \copydoc GRAM::ConstituentGas::specificHeatCapacity
};
typedef struct ConstituentGas_C ConstituentGas_C;

struct EphemerisState_C
{
  greal solarTime;         //!< \copydoc GRAM::EphemerisState::solarTime
  greal longitudeSun;      //!< \copydoc GRAM::EphemerisState::longitudeSun
  greal subsolarLatitude;  //!< \copydoc GRAM::EphemerisState::subsolarLatitude
  greal subsolarLongitude; //!< \copydoc GRAM::EphemerisState::subsolarLongitude
  greal orbitalRadius;     //!< \copydoc GRAM::EphemerisState::orbitalRadius
  greal oneWayLightTime;   //!< \copydoc GRAM::EphemerisState::oneWayLightTime
  greal solarZenithAngle;  //!< \copydoc GRAM::EphemerisState::solarZenithAngle
  greal secondsPerSol;     //!< \copydoc GRAM::EphemerisState::secondsPerSol
};
typedef struct EphemerisState_C EphemerisState_C;

struct PerturbationState_C
{
  greal densityRandomNumber;      //!< \copydoc GRAM::PerturbationState::densityRandomNumber
  greal ewWindRandomNumber;       //!< \copydoc GRAM::PerturbationState::ewWindRandomNumber
  greal nsWindRandomNumber;       //!< \copydoc GRAM::PerturbationState::nsWindRandomNumber
  greal verticalWindRandomNumber; //!< \copydoc GRAM::PerturbationState::verticalWindRandomNumber
};
typedef struct PerturbationState_C PerturbationState_C;

struct TrajectoryPoint_C
{
  Position_C position;            //!< The output position at this sample.
  DynamicsState_C dynamics;       //!< The dynamics output at this sample.
  WindsState_C winds;             //!< The wind output at this sample.
  EphemerisState_C ephemeris;     //!< The ephemeris output at this sample.
};
typedef struct TrajectoryPoint_C TrajectoryPoint_C;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The Generic C Interface for all GRAM models.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef __cplusplus
#include "Gram_C_wrapper.h"
extern "C" {
#endif

  enum GasType_C {
    ARGON_C, HELIUM_C, HYDROGEN_C, DIHYDROGEN_C, NITROGEN_C, DINITROGEN_C, OXYGEN_C, DIOXYGEN_C,
    METHANE_C, CARBON_MONOXIDE_C, CARBON_DIOXIDE_C, OZONE_C, NITROUS_OXIDE_C, WATER_C
  };

  enum GRAM_BODY_C { NO_BODY_C, VENUS_C, EARTH_C, MARS_C, JUPITER_C, SATURN_C, URANUS_C, NEPTUNE_C, TITAN_C};

  void enableConsoleOutput_C(int flag);
  void tryGetSpicePath_C(char* spicePath, int bufferSize);
  void tryGetDataPaths_C(char* spicePath, char* dataPath, int bufferSize);
  void setSpiceLsk_C(const char* lsk);
  void setSpicePck_C(const char* pck);
  void setSpiceKernel_C(enum GRAM_BODY_C body, const char* bsp);
  void initialize_C(const char* spicePath);
  void loadSpiceFile_C(const char* fileName);

  void* createAtmosphere_C(enum GRAM_BODY_C body);
  void* copyAtmosphere_C(void* atmos);
  void deleteAtmosphere_C(void* atmos);

  void setStartTime_C(void* atmos, const GramTime_C* time);
  void setPosition_C(void* atmos, const Position_C* pos);
  void setDelta_C(void* atmos, const Position_C* delta);
  void setSeed_C(void* atmos, int seed);
  void setMinRelativeStepSize_C(void* atmos, greal minRelativeStepSize);
  void setPerturbationScales_C(void* atmos, greal densityScale, greal ewWindScale, greal nsWindScale, greal verticalWindScale);
  void addAuxiliaryAtmosphere_C(void* atmos, const char* fileName, greal innerRadius, greal outerRadius, int isEastLongitudePositive);
  void setAuxiliaryValues_C(void* atmos, greal dens, greal pres, greal temp, greal ew, greal ns);
  void setPerturbationAction_C(void* atmos, int action);
  void setEphemerisState_C(void* atmos, const EphemerisState_C* state);
  void setEphemerisFastModeOn_C(void* atmos, int flag);
  void setSubsolarUpdateTime_C(void* atmos, greal utime);
  int generateTrajectory_C(void* atmos, const Position_C* initial, const Position_C* delta, int numberOfPoints, int updateInitialPerturbations, TrajectoryPoint_C* trajectory);

  int update_C(void* atmos);
  size_t getErrorMessage_C(char *message, size_t bufferSize);

  void getStartTime_C(void* atmos, GramTime_C* time);
  void getPosition_C(void* atmos, Position_C* position);
  void getDynamicsState_C(void* atmos, DynamicsState_C* state);
  void getDensityState_C(void* atmos, DensityState_C* state);
  void getWindsState_C(void* atmos, WindsState_C* state);
  void getGasesState_C(void* atmos, GasesState_C* state);
  void getConstituentGas_C(void* atmos, enum GasType_C gasTypeC, ConstituentGas_C* gas);
  void getEphemerisState_C(void* atmos, EphemerisState_C* state);
  void getPerturbationState_C(void* atmos, PerturbationState_C* state);

#ifdef __cplusplus
} // extern "C"
#endif
