//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Titan-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once
#ifdef _MSC_VER
#pragma warning( disable : 4250 )
#endif

#include "Gram_C.h"

#ifdef __cplusplus

#include "TitanAtmosphere.h"
using namespace GRAM;
typedef TitanAtmosphere TitanAtmosphere_C;

#else

typedef struct TitanAtmosphere TitanAtmosphere_C;

#endif


#ifdef __cplusplus
extern "C" {
#endif

  void tryGetSpicePath_T(char* spicePath, int bufferSize);
  void setSpiceLsk_T(const char* lsk);
  void setSpicePck_T(const char* pck);
  void setSpiceKernel_T(const char* bsp);
  void initialize_T(const char* spicePath);
  void loadSpiceFile_T(const char* fileName);

  TitanAtmosphere_C* createAtmosphere_T();
  TitanAtmosphere_C* copyAtmosphere_T(TitanAtmosphere_C* atmos);
  void deleteAtmosphere_T(TitanAtmosphere_C* atmos);

  void setStartTime_T(TitanAtmosphere_C* atmos, const GramTime_C* time);
  void setMinMaxFactor_T(TitanAtmosphere_C* atmos, greal factor, int computeFlag);
  void setMethaneMoleFraction_T(TitanAtmosphere_C *atmos, greal mmf);
  void setModelType_T(TitanAtmosphere_C *atmos, int type);
  void setSeed_T(TitanAtmosphere_C* atmos, int seed);
  void setMinRelativeStepSize_T(TitanAtmosphere_C* atmos, greal minRelativeStepSize);
  void setPerturbationScales_T(TitanAtmosphere_C* atmos, greal densityScale, greal ewWindScale, greal nsWindScale);
  void addAuxiliaryAtmosphere_T(TitanAtmosphere_C* atmos, const char* fileName, greal innerRadius, greal outerRadius, int isEastLongitudePositive);
  void setAuxiliaryValues_T(TitanAtmosphere_C* atmos, greal dens, greal pres, greal temp, greal ew, greal ns);

  void setPosition_T(TitanAtmosphere_C* atmos, const Position_C* pos);
  void setDelta_T(TitanAtmosphere_C* atmos, const Position_C* delta);
  void setPerturbationAction_T(TitanAtmosphere_C* atmos, int action);
  void setEphemerisState_T(TitanAtmosphere_C* atmos, const EphemerisState_C* state);
  void setEphemerisFastModeOn_T(TitanAtmosphere_C* atmos, int onFlag);
  void setSubsolarUpdateTime_T(TitanAtmosphere_C* atmos, greal utime);

  int update_T(TitanAtmosphere_C* atmos);
  size_t getErrorMessage_T(char *message, size_t bufferSize);

  void getStartTime_T(TitanAtmosphere_C* atmos, GramTime_C* time);
  void getMinMaxFactor_T(TitanAtmosphere_C* atmos, greal *factor);
  size_t getVersionString_T(TitanAtmosphere_C* atmos, char* versionString, size_t bufferSize);
  void getPosition_T(TitanAtmosphere_C* atmos, Position_C* position);
  void getDynamicsState_T(TitanAtmosphere_C* atmos, DynamicsState_C* dstate);
  void getDensityState_T(TitanAtmosphere_C* atmos, DensityState_C* state);
  void getWindsState_T(TitanAtmosphere_C* atmos, WindsState_C* state);
  void getGasesState_T(TitanAtmosphere_C* atmos, GasesState_C* state,
    ConstituentGas_C* argon, ConstituentGas_C* methane, ConstituentGas_C* dinitrogen);
  void getEphemerisState_T(TitanAtmosphere_C* atmos, EphemerisState_C* state);
  void getPerturbationState_T(TitanAtmosphere_C* atmos, PerturbationState_C* state);

#ifdef __cplusplus
} // extern "C"
#endif
