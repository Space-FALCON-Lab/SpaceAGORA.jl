//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once
#ifdef _MSC_VER
#pragma warning( disable : 4250 )
#endif

#include "Gram_C.h"

struct MarsState_C
{
  greal planetoGraphicHeight;         //!< \copydoc GRAM::MarsAtmosphereState::planetoGraphicHeight
  greal planetoGraphicLatitude;       //!< \copydoc GRAM::MarsAtmosphereState::planetoGraphicLatitude
  greal referenceHeight;              //!< \copydoc GRAM::MarsAtmosphereState::referenceHeight
  greal referenceRadius;              //!< \copydoc GRAM::MarsAtmosphereState::referenceRadius
  greal groundTemperature;            //!< \copydoc GRAM::MarsAtmosphereState::groundTemperature
  greal thermosphereBaseHeight;       //!< \copydoc GRAM::MarsAtmosphereState::thermosphereBaseHeight
  greal thermosphereBaseTemperature;  //!< \copydoc GRAM::MarsAtmosphereState::thermosphereBaseTemperature
  greal exosphericTemperature;        //!< \copydoc GRAM::MarsAtmosphereState::exosphericTemperature
  greal f1PeakHeight;                 //!< \copydoc GRAM::MarsAtmosphereState::f1PeakHeight
  greal albedo;                       //!< \copydoc GRAM::MarsAtmosphereState::albedo
  greal heightOffset;                 //!< \copydoc GRAM::MarsAtmosphereState::heightOffset
  greal localHeightOffset;            //!< \copydoc GRAM::MarsAtmosphereState::localHeightOffset
  greal dustOpticalDepth;             //!< \copydoc GRAM::MarsAtmosphereState::dustOpticalDepth
  greal dustColumnArealDensity;       //!< \copydoc GRAM::MarsAtmosphereState::dustColumnArealDensity
  greal dustMixingRatio;              //!< \copydoc GRAM::MarsAtmosphereState::dustMixingRatio
  greal dustMassDensity;              //!< \copydoc GRAM::MarsAtmosphereState::dustMassDensity
  greal dustNumberDensity;            //!< \copydoc GRAM::MarsAtmosphereState::dustNumberDensity
  int   iceIsPresent;                 //!< \copydoc GRAM::MarsAtmosphereState::iceIsPresent
};
typedef struct MarsState_C MarsState_C;

struct DailyDynamicsState_C
{
  greal temperatureDaily;         //!< \copydoc GRAM::MarsAtmosphereState::temperatureDaily
  greal pressureDaily;            //!< \copydoc GRAM::MarsAtmosphereState::pressureDaily
  greal densityDaily;             //!< \copydoc GRAM::MarsAtmosphereState::densityDaily
  greal ewWindDaily;              //!< \copydoc GRAM::MarsAtmosphereState::ewWindDaily
  greal nsWindDaily;              //!< \copydoc GRAM::MarsAtmosphereState::nsWindDaily
  greal densityMin;               //!< \copydoc GRAM::MarsAtmosphereState::densityMin
  greal densityMax;               //!< \copydoc GRAM::MarsAtmosphereState::densityMax
  greal temperatureMin;           //!< \copydoc GRAM::MarsAtmosphereState::temperatureMin
  greal temperatureMax;           //!< \copydoc GRAM::MarsAtmosphereState::temperatureMax
};
typedef struct DailyDynamicsState_C DailyDynamicsState_C;

#ifdef __cplusplus

#include "MarsAtmosphere.h"
using namespace GRAM;
typedef MarsAtmosphere MarsAtmosphere_C;

#else

typedef struct MarsAtmosphere MarsAtmosphere_C;
typedef double(*TopoCallback)(double, double, void*);

#endif


#ifdef __cplusplus
extern "C" {
#endif

  void tryGetSpicePath_M(char* spicePath, int bufferSize);
  void tryGetDataPaths_M(char* spicePath, char* dataPath, int bufferSize);
  void setSpiceLsk_M(const char* lsk);
  void setSpicePck_M(const char* pck);
  void setSpiceKernel_M(const char* bsp);
  void initialize_M(const char* spicePath);
  void loadSpiceFile_M(const char* fileName);

  MarsAtmosphere_C* createAtmosphere_M(const char* dataPath);
  MarsAtmosphere_C* copyAtmosphere_M(MarsAtmosphere_C* atmos);
  void deleteAtmosphere_M(MarsAtmosphere_C* atmos);

  void setPlanetaryRadii_M(MarsAtmosphere_C* atmos, greal equatorial, greal polar);
  void setMapYear_M(MarsAtmosphere_C* atmos, int year);
  void setHeightOffsetModel_M(MarsAtmosphere_C* atmos, int model, greal hgtOffset);
  void setHeightAboveSurface_M(MarsAtmosphere_C* atmos, greal heightAboveSurface);
  void setMGCMDustLevels_M(MarsAtmosphere_C* atmos, greal constantDustLevel, greal minDustLevel, greal maxDustLevel);
  void setF107_M(MarsAtmosphere_C* atmos, greal f107);
  void setPerturbationWaveLengthScale_M(MarsAtmosphere_C* atmos, greal scale);
  void setMOLAHeights_M(MarsAtmosphere_C* atmos, int isMola);
  void setMinMax_M(MarsAtmosphere_C* atmos, int minMax);
  void setDustStorm_M(MarsAtmosphere_C* atmos, greal longitudeSun, greal duration, greal intensity, greal maxRadius, greal lat, greal lon);
  void setDustDensity_M(MarsAtmosphere_C* atmos, greal nu, greal diameter, greal density);
  void setExosphericTemperature_M(MarsAtmosphere_C* atmos, greal offset, greal factor);
  void setWaveDefaults_M(MarsAtmosphere_C* atmos, greal date, greal scale, greal mean, greal a1, greal p1, greal r1, greal a2, greal p2, greal r2, greal a3, greal p3, greal r3);
  void setWaveFile_M(MarsAtmosphere_C* atmos, const char* waveFile);
  void setWindScales_M(MarsAtmosphere_C* atmos, greal meanWinds, greal boundaryLayerWinds);

  void setStartTime_M(MarsAtmosphere_C* atmos, const GramTime_C* time);
  void setSeed_M(MarsAtmosphere_C* atmos, int seed);
  void setMinRelativeStepSize_M(MarsAtmosphere_C* atmos, greal minRelativeStepSize);
  void setPerturbationScales_M(MarsAtmosphere_C* atmos, greal densityScale, greal ewWindScale, greal nsWindScale, greal verticalWindScale);
  void addAuxiliaryAtmosphere_M(MarsAtmosphere_C* atmos, const char* fileName, greal innerRadius, greal outerRadius, int isEastLongitudePositive);
  void setAuxiliaryValues_M(MarsAtmosphere_C* atmos, greal dens, greal pres, greal temp, greal ew, greal ns);

  void setPosition_M(MarsAtmosphere_C* atmos, const Position_C* pos);
  void setDelta_M(MarsAtmosphere_C* atmos, const Position_C* delta);
  void setPerturbationAction_M(MarsAtmosphere_C* atmos, int action);
  void setEphemerisState_M(MarsAtmosphere_C* atmos, const EphemerisState_C* state);
  void setEphemerisFastModeOn_M(MarsAtmosphere_C* atmos, int onFlag);
  void setSubsolarUpdateTime_M(MarsAtmosphere_C* atmos, greal utime);

  void setAreoidRadiusCallback_M(MarsAtmosphere_C* atmos, TopoCallback callback);
  void setTopographicHeightCallback_M(MarsAtmosphere_C* atmos, TopoCallback callback);
  void setCallbackData_M(MarsAtmosphere_C* atmos, void* pointer);

  int update_M(MarsAtmosphere_C* atmos);
  size_t getErrorMessage_M(char *message, size_t bufferSize);

  void getStartTime_M(MarsAtmosphere_C* atmos, GramTime_C* time);
  size_t getVersionString_M(MarsAtmosphere_C* atmos, char* versionString, size_t bufferSize);
  void getPosition_M(MarsAtmosphere_C* atmos, Position_C* position);
  void getDynamicsState_M(MarsAtmosphere_C* atmos, DynamicsState_C* dstate);
  void getDensityState_M(MarsAtmosphere_C* atmos, DensityState_C* state);
  void getWindsState_M(MarsAtmosphere_C* atmos, WindsState_C* state);
  void getGasesState_M(MarsAtmosphere_C* atmos, GasesState_C* state,
    ConstituentGas_C* argon, ConstituentGas_C* carbonDioxide, ConstituentGas_C* carbonMonoxide,
    ConstituentGas_C* dihydrogen, ConstituentGas_C* dinitrogen, ConstituentGas_C* dioxygen,
    ConstituentGas_C* helium, ConstituentGas_C* hydrogen, ConstituentGas_C* oxygen, ConstituentGas_C* water);
  void getEphemerisState_M(MarsAtmosphere_C* atmos, EphemerisState_C* state);
  void getPerturbationState_M(MarsAtmosphere_C* atmos, PerturbationState_C* state);
  void getDailyDynamicsState_M(MarsAtmosphere_C* atmos, DailyDynamicsState_C* state);
  void getMarsState_M(MarsAtmosphere_C* atmos, MarsState_C* state);

#ifdef __cplusplus
} // extern "C"
#endif
