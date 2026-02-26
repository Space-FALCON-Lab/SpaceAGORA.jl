//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Jupiter-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once
#ifdef _MSC_VER
#pragma warning( disable : 4250 )
#endif

#include "Gram_C.h"

#ifdef __cplusplus

#include "JupiterAtmosphere.h"
using namespace GRAM;
typedef JupiterAtmosphere JupiterAtmosphere_C;

#else

typedef struct JupiterAtmosphere JupiterAtmosphere_C;

#endif


#ifdef __cplusplus
extern "C" {
#endif

  void tryGetSpicePath_J(char* spicePath, int bufferSize);
  void setSpiceLsk_J(const char* lsk);
  void setSpicePck_J(const char* pck);
  void setSpiceKernel_J(const char* bsp);
  void initialize_J(const char* spicePath);
  void loadSpiceFile_J(const char* fileName);

  JupiterAtmosphere_C* createAtmosphere_J();
  JupiterAtmosphere_C* copyAtmosphere_J(JupiterAtmosphere_C* atmos);
  void deleteAtmosphere_J(JupiterAtmosphere_C* atmos);

  void setStartTime_J(JupiterAtmosphere_C* atmos, const GramTime_C* time);
  void setSeed_J(JupiterAtmosphere_C* atmos, int seed);
  void setMinRelativeStepSize_J(JupiterAtmosphere_C* atmos, greal minRelativeStepSize);
  void setPerturbationScales_J(JupiterAtmosphere_C* atmos, greal densityScale, greal ewWindScale, greal nsWindScale);
  void addAuxiliaryAtmosphere_J(JupiterAtmosphere_C* atmos, const char* fileName, greal innerRadius, greal outerRadius, int isEastLongitudePositive);
  void setAuxiliaryValues_J(JupiterAtmosphere_C* atmos, greal dens, greal pres, greal temp, greal ew, greal ns);

  void setPosition_J(JupiterAtmosphere_C* atmos, const Position_C* pos);
  void setDelta_J(JupiterAtmosphere_C* atmos, const Position_C* delta);
  void setPerturbationAction_J(JupiterAtmosphere_C* atmos, int action);
  void setEphemerisState_J(JupiterAtmosphere_C* atmos, const EphemerisState_C* state);
  void setEphemerisFastModeOn_J(JupiterAtmosphere_C* atmos, int onFlag);
  void setSubsolarUpdateTime_J(JupiterAtmosphere_C* atmos, greal utime);

  int update_J(JupiterAtmosphere_C* atmos);
  size_t getErrorMessage_J(char *message, size_t bufferSize);

  void getStartTime_J(JupiterAtmosphere_C* atmos, GramTime_C* time);
  size_t getVersionString_J(JupiterAtmosphere_C* atmos, char* versionString, size_t bufferSize);
  void getPosition_J(JupiterAtmosphere_C* atmos, Position_C* position);
  void getDynamicsState_J(JupiterAtmosphere_C* atmos, DynamicsState_C* dstate);
  void getDensityState_J(JupiterAtmosphere_C* atmos, DensityState_C* state);
  void getWindsState_J(JupiterAtmosphere_C* atmos, WindsState_C* state);
  void getGasesState_J(JupiterAtmosphere_C* atmos, GasesState_C* state);
  void getEphemerisState_J(JupiterAtmosphere_C* atmos, EphemerisState_C* state);
  void getPerturbationState_J(JupiterAtmosphere_C* atmos, PerturbationState_C* state);

#ifdef __cplusplus
} // extern "C"
#endif
