//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Earth-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once
#ifdef _MSC_VER
#pragma warning( disable : 4250 )
#endif

#include "Gram_C.h"

struct EarthState_C
{
  greal perturbedTemperature;       //!< \copydoc GRAM::EarthAtmosphereState::perturbedTemperature
  greal temperaturePerturbation;    //!< \copydoc GRAM::EarthAtmosphereState::temperaturePerturbation
  greal temperatureStandardDeviation; //!< \copydoc GRAM::AtmosphereState::temperatureStandardDeviation
  greal perturbedPressure;          //!< \copydoc GRAM::EarthAtmosphereState::perturbedPressure
  greal pressurePerturbation;       //!< \copydoc GRAM::EarthAtmosphereState::pressurePerturbation
  greal pressureStandardDeviation;  //!< \copydoc GRAM::AtmosphereState::pressureStandardDeviation

  greal vaporPressure;              //!< \copydoc GRAM::EarthAtmosphereState::vaporPressure
  greal vaporPressureSD;            //!< \copydoc GRAM::EarthAtmosphereState::vaporPressureSD
  greal vaporDensity;               //!< \copydoc GRAM::EarthAtmosphereState::vaporDensity
  greal vaporDensitySD;             //!< \copydoc GRAM::EarthAtmosphereState::vaporDensitySD
  greal dewPoint;                   //!< \copydoc GRAM::EarthAtmosphereState::dewPoint
  greal dewPointSD;                 //!< \copydoc GRAM::EarthAtmosphereState::dewPointSD
  greal relativeHumidity;           //!< \copydoc GRAM::EarthAtmosphereState::relativeHumidity
  greal relativeHumiditySD;         //!< \copydoc GRAM::EarthAtmosphereState::relativeHumiditySD

  greal geodeticLatitude;           //!< \copydoc GRAM::EarthAtmosphereState::geodeticLatitude
  greal rraWeight;                  //!< \copydoc GRAM::EarthAtmosphereState::rraWeight
  char rraSiteName[6];              //!< \copydoc GRAM::EarthAtmosphereState::rraSiteName

  greal windSpeed;                  //!< \copydoc GRAM::EarthAtmosphereState::windSpeed
  greal windSpeedStandardDeviation; //!< \copydoc GRAM::EarthAtmosphereState::windSpeedStandardDeviation
  greal windCorrelation;            //!< \copydoc GRAM::EarthAtmosphereState::windCorrelation
  int  severityLevel;               //!< \copydoc GRAM::EarthAtmosphereState::severityLevel
};
typedef struct EarthState_C EarthState_C;

struct EarthPerts_C
{
  greal presPertSmall;              //!< \copydoc GRAM::EarthAtmosphereState::presPertSmall
  greal densPertSmall;              //!< \copydoc GRAM::EarthAtmosphereState::densPertSmall
  greal tempPertSmall;              //!< \copydoc GRAM::EarthAtmosphereState::tempPertSmall
  greal ewWindPertSmall;            //!< \copydoc GRAM::EarthAtmosphereState::ewWindPertSmall
  greal nsWindPertSmall;            //!< \copydoc GRAM::EarthAtmosphereState::nsWindPertSmall
                                    
  greal presSDSmall;                //!< \copydoc GRAM::EarthAtmosphereState::presSDSmall
  greal densSDSmall;                //!< \copydoc GRAM::EarthAtmosphereState::densSDSmall
  greal tempSDSmall;                //!< \copydoc GRAM::EarthAtmosphereState::tempSDSmall
  greal ewWindSDSmall;              //!< \copydoc GRAM::EarthAtmosphereState::ewWindSDSmall
  greal nsWindSDSmall;              //!< \copydoc GRAM::EarthAtmosphereState::nsWindSDSmall
                                    
  greal presPertLarge;              //!< \copydoc GRAM::EarthAtmosphereState::presPertLarge
  greal densPertLarge;              //!< \copydoc GRAM::EarthAtmosphereState::densPertLarge
  greal tempPertLarge;              //!< \copydoc GRAM::EarthAtmosphereState::tempPertLarge
  greal ewWindPertLarge;            //!< \copydoc GRAM::EarthAtmosphereState::ewWindPertLarge
  greal nsWindPertLarge;            //!< \copydoc GRAM::EarthAtmosphereState::nsWindPertLarge
                                    
  greal presSDLarge;                //!< \copydoc GRAM::EarthAtmosphereState::presSDLarge
  greal densSDLarge;                //!< \copydoc GRAM::EarthAtmosphereState::densSDLarge
  greal tempSDLarge;                //!< \copydoc GRAM::EarthAtmosphereState::tempSDLarge
  greal ewWindSDLarge;              //!< \copydoc GRAM::EarthAtmosphereState::ewWindSDLarge
  greal nsWindSDLarge;              //!< \copydoc GRAM::EarthAtmosphereState::nsWindSDLarge
};
typedef struct EarthPerts_C EarthPerts_C;

struct EarthSurface_C
{
  greal windSpeedAtSurface;         //!< \copydoc GRAM::EarthAtmosphereState::windSpeedAtSurface
  greal windSpeedSDAtSurface;       //!< \copydoc GRAM::EarthAtmosphereState::windSpeedSDAtSurface
  greal temperatureAtSurface;       //!< \copydoc GRAM::EarthAtmosphereState::temperatureAtSurface
  greal temperatureSDAtSurface;     //!< \copydoc GRAM::EarthAtmosphereState::temperatureSDAtSurface
  greal pressureSDAtSurface;        //!< \copydoc GRAM::EarthAtmosphereState::pressureSDAtSurface
  greal densitySDAtSurface;         //!< \copydoc GRAM::EarthAtmosphereState::densitySDAtSurface
  greal densityAtSurface;           //!< \copydoc GRAM::EarthAtmosphereState::densityAtSurface
  greal ewWindAtSurface;            //!< \copydoc GRAM::EarthAtmosphereState::ewWindAtSurface
  greal nsWindAtSurface;            //!< \copydoc GRAM::EarthAtmosphereState::nsWindAtSurface
  greal ewWindSDAtSurface;          //!< \copydoc GRAM::EarthAtmosphereState::ewWindSDAtSurface
  greal nsWindSDAtSurface;          //!< \copydoc GRAM::EarthAtmosphereState::nsWindSDAtSurface
  greal windCorrelationAtSurface;   //!< \copydoc GRAM::EarthAtmosphereState::windCorrelationAtSurface
};
typedef struct EarthSurface_C EarthSurface_C;

struct EarthBoundaryLayer_C
{
  int   landCode;                      //!< \copydoc GRAM::EarthAtmosphereState::landCode
  greal surfaceRoughness;              //!< \copydoc GRAM::EarthAtmosphereState::surfaceRoughness
  greal netRadiationIndex;             //!< \copydoc GRAM::EarthAtmosphereState::netRadiationIndex
  greal stability;                     //!< \copydoc GRAM::EarthAtmosphereState::stability
  greal inverseLength;                 //!< \copydoc GRAM::EarthAtmosphereState::inverseLength
  greal frictionVelocity;              //!< \copydoc GRAM::EarthAtmosphereState::frictionVelocity
  greal BVFrequencySquare;             //!< \copydoc GRAM::EarthAtmosphereState::BVFrequencySquare
  greal boundaryLayerDepth;            //!< \copydoc GRAM::EarthAtmosphereState::boundaryLayerDepth
  greal neutralBoundaryLayerDepth;     //!< \copydoc GRAM::EarthAtmosphereState::neutralBoundaryLayerDepth
  greal unstableBLFactor;              //!< \copydoc GRAM::EarthAtmosphereState::unstableBLFactor
  greal sigmaRatio;                    //!< \copydoc GRAM::EarthAtmosphereState::sigmaRatio
  greal sigmaW;                        //!< \copydoc GRAM::EarthAtmosphereState::sigmaW
  greal metersAboveSurface;            //!< \copydoc GRAM::EarthAtmosphereState::metersAboveSurface
  greal perturbedWindSpeedAtSurface;   //!< \copydoc GRAM::EarthAtmosphereState::perturbedWindSpeedAtSurface
};
typedef struct EarthBoundaryLayer_C EarthBoundaryLayer_C;

#ifdef __cplusplus

#include "EarthAtmosphere.h"
using namespace GRAM;
typedef EarthAtmosphere EarthAtmosphere_C;

#else

typedef struct EarthAtmosphere EarthAtmosphere_C;

#endif


#ifdef __cplusplus
extern "C" {
#endif

  void tryGetSpicePath_E(char* spicePath, int bufferSize);
  void tryGetDataPaths_E(char* spicePath, char* dataPath, int bufferSize);
  void setSpiceLsk_E(const char* lsk);
  void setSpicePck_E(const char* pck);
  void setSpiceKernel_E(const char* bsp);
  void initialize_E(const char* spicePath);
  void loadSpiceFile_E(const char* fileName);

  EarthAtmosphere_C* createAtmosphere_E(const char* dataPath);
  EarthAtmosphere_C* copyAtmosphere_E(EarthAtmosphere_C* atmos);
  void deleteAtmosphere_E(EarthAtmosphere_C* atmos);

  void setThermosphereModel_E(EarthAtmosphere_C* atmos, int model);
  void setSurfaceRoughness_E(EarthAtmosphere_C* atmos, greal z0);
  void setAtmosPath_E(EarthAtmosphere_C* atmos, const char* path);
  void setMERRA2Path_E(EarthAtmosphere_C* atmos, const char* path);
  void setNCEPPath_E(EarthAtmosphere_C* atmos, const char* path);
  void setMERRA2Parameters_E(EarthAtmosphere_C* atmos, int M2Hour, greal latMin, greal latMax, greal lonMin, greal lonMax);
  void setUseNCEP_E(EarthAtmosphere_C* atmos, int useNCEP);
  void setNCEPParameters_E(EarthAtmosphere_C* atmos, int NCEPYear, int NCEPHour);
  void setRRAPath_E(EarthAtmosphere_C* atmos, const char* path);
  void setRRASiteList_E(EarthAtmosphere_C* atmos, const char* fileName);
  void setRRAParameters_E(EarthAtmosphere_C* atmos, int year, greal innerRadius, greal outerRadius);
  void setUseRRA_E(EarthAtmosphere_C* atmos, int useFlag);
  void setSolarParameters_E(EarthAtmosphere_C* atmos, greal dailyF10, greal meanF10, greal ap);
  void setJB2008Parameters_E(EarthAtmosphere_C* atmos, greal dailyS10, greal meanS10, greal dailyXM10, greal meanXM10, greal dailyY10, greal meanY10, greal dstdtc);
  void setInitialPerturbations_E(EarthAtmosphere_C* atmos, greal densPert, greal tempPert, greal ewPert, greal nsPert, greal verticalPert);
  void setPatchiness_E(EarthAtmosphere_C* atmos, int usePatchiness);

  void setStartTime_E(EarthAtmosphere_C* atmos, const GramTime_C* time);
  void setSeed_E(EarthAtmosphere_C* atmos, int seed);
  void setPerturbationScales_E(EarthAtmosphere_C* atmos, greal randomScale, greal horizontalWindScale, greal verticalWindScale);
  void addAuxiliaryAtmosphere_E(EarthAtmosphere_C* atmos, const char* fileName, greal innerRadius, greal outerRadius, int isEastLongitudePositive);
  void setAuxiliaryValues_E(EarthAtmosphere_C* atmos, greal dens, greal pres, greal temp, greal ew, greal ns);

  void setPosition_E(EarthAtmosphere_C* atmos, const Position_C* pos);
  void setDelta_E(EarthAtmosphere_C* atmos, const Position_C* delta);
  void setPerturbationAction_E(EarthAtmosphere_C* atmos, int action);
  void setEphemerisState_E(EarthAtmosphere_C* atmos, const EphemerisState_C* state);
  void setEphemerisFastModeOn_E(EarthAtmosphere_C* atmos, int onFlag);
  void setSubsolarUpdateTime_E(EarthAtmosphere_C* atmos, greal utime);

  int update_E(EarthAtmosphere_C* atmos);
  size_t getErrorMessage_E(char *message, size_t bufferSize);

  void getStartTime_E(EarthAtmosphere_C* atmos, GramTime_C* time);
  size_t getVersionString_E(EarthAtmosphere_C* atmos, char* versionString, size_t bufferSize);
  void getPosition_E(EarthAtmosphere_C* atmos, Position_C* position);
  void getDynamicsState_E(EarthAtmosphere_C* atmos, DynamicsState_C* dstate);
  void getDensityState_E(EarthAtmosphere_C* atmos, DensityState_C* state);
  void getWindsState_E(EarthAtmosphere_C* atmos, WindsState_C* state);
  void getGasesState_E(EarthAtmosphere_C* atmos, GasesState_C* state,
    ConstituentGas_C* argon, ConstituentGas_C* carbonDioxide, ConstituentGas_C* carbonMonoxide, ConstituentGas_C* dinitrogen,
    ConstituentGas_C* dioxygen, ConstituentGas_C* helium, ConstituentGas_C* hydrogen, ConstituentGas_C* methane,
    ConstituentGas_C* nitrogen, ConstituentGas_C* nitrousOxide, ConstituentGas_C* oxygen, ConstituentGas_C* ozone, ConstituentGas_C* water);
  void getEphemerisState_E(EarthAtmosphere_C* atmos, EphemerisState_C* state);
  void getPerturbationState_E(EarthAtmosphere_C* atmos, PerturbationState_C* state);
  void getEarthState_E(EarthAtmosphere_C* atmos, EarthState_C* state);
  void getEarthPerts_E(EarthAtmosphere_C* atmos, EarthPerts_C* state);
  void getEarthSurface_E(EarthAtmosphere_C* atmos, EarthSurface_C* state);
  void getEarthBoundaryLayer_E(EarthAtmosphere_C* atmos, EarthBoundaryLayer_C* state);

#ifdef __cplusplus
} // extern "C"
#endif
