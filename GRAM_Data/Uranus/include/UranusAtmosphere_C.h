//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Uranus-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once
#ifdef _MSC_VER
#pragma warning( disable : 4250 )
#endif

#include "Gram_C.h"

#ifdef __cplusplus

#include "UranusAtmosphere.h"
using namespace GRAM;
typedef UranusAtmosphere UranusAtmosphere_C;

#else

typedef struct UranusAtmosphere UranusAtmosphere_C;

#endif


#ifdef __cplusplus
extern "C" {
#endif

  void tryGetSpicePath_U(char* spicePath, int bufferSize);
  void setSpiceLsk_U(const char* lsk);
  void setSpicePck_U(const char* pck);
  void setSpiceKernel_U(const char* bsp);
  void initialize_U(const char* spicePath);
  void loadSpiceFile_U(const char* fileName);

  UranusAtmosphere_C* createAtmosphere_U();
  UranusAtmosphere_C* copyAtmosphere_U(UranusAtmosphere_C* atmos);
  void deleteAtmosphere_U(UranusAtmosphere_C* atmos);

  void setStartTime_U(UranusAtmosphere_C* atmos, const GramTime_C* time);
  void setSeed_U(UranusAtmosphere_C* atmos, int seed);
  void setMinRelativeStepSize_U(UranusAtmosphere_C* atmos, greal minRelativeStepSize);
  void setPerturbationScales_U(UranusAtmosphere_C* atmos, greal densityScale, greal ewWindScale, greal nsWindScale);
  void addAuxiliaryAtmosphere_U(UranusAtmosphere_C* atmos, const char* fileName, greal innerRadius, greal outerRadius, int isEastLongitudePositive);
  void setAuxiliaryValues_U(UranusAtmosphere_C* atmos, greal dens, greal pres, greal temp, greal ew, greal ns);

  void setPosition_U(UranusAtmosphere_C* atmos, const Position_C* pos);
  void setDelta_U(UranusAtmosphere_C* atmos, const Position_C* delta);
  void setPerturbationAction_U(UranusAtmosphere_C* atmos, int action);
  void setEphemerisState_U(UranusAtmosphere_C* atmos, const EphemerisState_C* state);
  void setEphemerisFastModeOn_U(UranusAtmosphere_C* atmos, int onFlag);
  void setSubsolarUpdateTime_U(UranusAtmosphere_C* atmos, greal utime);

  int update_U(UranusAtmosphere_C* atmos);
  size_t getErrorMessage_U(char *message, size_t bufferSize);

  void getStartTime_U(UranusAtmosphere_C* atmos, GramTime_C* time);
  size_t getVersionString_U(UranusAtmosphere_C* atmos, char* versionString, size_t bufferSize);
  void getPosition_U(UranusAtmosphere_C* atmos, Position_C* position);
  void getDynamicsState_U(UranusAtmosphere_C* atmos, DynamicsState_C* dstate);
  void getDensityState_U(UranusAtmosphere_C* atmos, DensityState_C* state);
  void getWindsState_U(UranusAtmosphere_C* atmos, WindsState_C* state);
  void getGasesState_U(UranusAtmosphere_C* atmos, GasesState_C* state,
    ConstituentGas_C* dihydrogen, ConstituentGas_C* methane, ConstituentGas_C* helium);
  void getEphemerisState_U(UranusAtmosphere_C* atmos, EphemerisState_C* state);
  void getPerturbationState_U(UranusAtmosphere_C* atmos, PerturbationState_C* state);

#ifdef __cplusplus
} // extern "C"
#endif
