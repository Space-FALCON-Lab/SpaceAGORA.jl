//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Neptune-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once
#ifdef _MSC_VER
#pragma warning( disable : 4250 )
#endif

#include "Gram_C.h"

#ifdef __cplusplus

#include "NeptuneAtmosphere.h"
using namespace GRAM;
typedef NeptuneAtmosphere NeptuneAtmosphere_C;

#else

typedef struct NeptuneAtmosphere NeptuneAtmosphere_C;

#endif


#ifdef __cplusplus
extern "C" {
#endif

  void tryGetSpicePath_N(char* spicePath, int bufferSize);
  void setSpiceLsk_N(const char* lsk);
  void setSpicePck_N(const char* pck);
  void setSpiceKernel_N(const char* bsp);
  void initialize_N(const char* spicePath);
  void loadSpiceFile_N(const char* fileName);

  NeptuneAtmosphere_C* createAtmosphere_N();
  NeptuneAtmosphere_C* copyAtmosphere_N(NeptuneAtmosphere_C* atmos);
  void deleteAtmosphere_N(NeptuneAtmosphere_C* atmos);

  void setStartTime_N(NeptuneAtmosphere_C* atmos, const GramTime_C* time);
  void setDinitrogenMoleFraction_N(NeptuneAtmosphere_C* atmos, greal n2mf);
  void setMinMaxFactor_N(NeptuneAtmosphere_C* atmos, greal factor, int computeFlag);
  void setSeed_N(NeptuneAtmosphere_C* atmos, int seed);
  void setMinRelativeStepSize_N(NeptuneAtmosphere_C* atmos, greal minRelativeStepSize);
  void setPerturbationScales_N(NeptuneAtmosphere_C* atmos, greal densityScale, greal ewWindScale, greal nsWindScale);
  void addAuxiliaryAtmosphere_N(NeptuneAtmosphere_C* atmos, const char* fileName, greal innerRadius, greal outerRadius, int isEastLongitudePositive);
  void setAuxiliaryValues_N(NeptuneAtmosphere_C* atmos, greal dens, greal pres, greal temp, greal ew, greal ns);

  void setPosition_N(NeptuneAtmosphere_C* atmos, const Position_C* pos);
  void setDelta_N(NeptuneAtmosphere_C* atmos, const Position_C* delta);
  void setPerturbationAction_N(NeptuneAtmosphere_C* atmos, int action);
  void setEphemerisState_N(NeptuneAtmosphere_C* atmos, const EphemerisState_C* state);
  void setEphemerisFastModeOn_N(NeptuneAtmosphere_C* atmos, int onFlag);
  void setSubsolarUpdateTime_N(NeptuneAtmosphere_C* atmos, greal utime);

  int update_N(NeptuneAtmosphere_C* atmos);
  size_t getErrorMessage_N(char *message, size_t bufferSize);

  void getStartTime_N(NeptuneAtmosphere_C* atmos, GramTime_C* time);
  void getMinMaxFactor_N(NeptuneAtmosphere_C* atmos, greal *factor);
  size_t getVersionString_N(NeptuneAtmosphere_C* atmos, char* versionString, size_t bufferSize);
  void getPosition_N(NeptuneAtmosphere_C* atmos, Position_C* position);
  void getDynamicsState_N(NeptuneAtmosphere_C* atmos, DynamicsState_C* dstate);
  void getDensityState_N(NeptuneAtmosphere_C* atmos, DensityState_C* state);
  void getWindsState_N(NeptuneAtmosphere_C* atmos, WindsState_C* state);
  void getGasesState_N(NeptuneAtmosphere_C* atmos, GasesState_C* state,
    ConstituentGas_C* dihydrogen, ConstituentGas_C* methane,
    ConstituentGas_C* helium, ConstituentGas_C* dinitrogen);
  void getEphemerisState_N(NeptuneAtmosphere_C* atmos, EphemerisState_C* state);
  void getPerturbationState_N(NeptuneAtmosphere_C* atmos, PerturbationState_C* state);

#ifdef __cplusplus
} // extern "C"
#endif
