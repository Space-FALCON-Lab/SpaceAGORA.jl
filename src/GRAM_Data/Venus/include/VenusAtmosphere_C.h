//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once
#ifdef _MSC_VER
#pragma warning( disable : 4250 )
#endif

#include "Gram_C.h"

#ifdef __cplusplus

#include "VenusAtmosphere.h"
using namespace GRAM;
typedef VenusAtmosphere VenusAtmosphere_C;

#else

typedef struct VenusAtmosphere VenusAtmosphere_C;

#endif


#ifdef __cplusplus
extern "C" {
#endif

  void tryGetSpicePath_V(char* spicePath, int bufferSize);
  void setSpiceLsk_V(const char* lsk);
  void setSpicePck_V(const char* pck);
  void setSpiceKernel_V(const char* bsp);
  void initialize_V(const char* spicePath);
  void loadSpiceFile_V(const char* fileName);

  VenusAtmosphere_C* createAtmosphere_V();
  VenusAtmosphere_C* copyAtmosphere_V(VenusAtmosphere_C* atmos);
  void deleteAtmosphere_V(VenusAtmosphere_C* atmos);

  void setStartTime_V(VenusAtmosphere_C* atmos, const GramTime_C* time);
  void setSeed_V(VenusAtmosphere_C* atmos, int seed);
  void setMinRelativeStepSize_V(VenusAtmosphere_C* atmos, greal minRelativeStepSize);
  void setPerturbationScales_V(VenusAtmosphere_C* atmos, greal densityScale, greal ewWindScale, greal nsWindScale);
  void addAuxiliaryAtmosphere_V(VenusAtmosphere_C* atmos, const char* fileName, greal innerRadius, greal outerRadius, int isEastLongitudePositive);
  void setAuxiliaryValues_V(VenusAtmosphere_C* atmos, greal dens, greal pres, greal temp, greal ew, greal ns);

  void setPosition_V(VenusAtmosphere_C* atmos, const Position_C* pos);
  void setDelta_V(VenusAtmosphere_C* atmos, const Position_C* delta);
  void setPerturbationAction_V(VenusAtmosphere_C* atmos, int action);
  void setEphemerisState_V(VenusAtmosphere_C* atmos, const EphemerisState_C* state);
  void setEphemerisFastModeOn_V(VenusAtmosphere_C* atmos, int onFlag);
  void setSubsolarUpdateTime_V(VenusAtmosphere_C* atmos, greal utime);

  int update_V(VenusAtmosphere_C* atmos);
  size_t getErrorMessage_V(char *message, size_t bufferSize);

  void getStartTime_V(VenusAtmosphere_C* atmos, GramTime_C* time);
  size_t getVersionString_V(VenusAtmosphere_C* atmos, char* versionString, size_t bufferSize);
  void getPosition_V(VenusAtmosphere_C* atmos, Position_C* position);
  void getDynamicsState_V(VenusAtmosphere_C* atmos, DynamicsState_C* dstate);
  void getDensityState_V(VenusAtmosphere_C* atmos, DensityState_C* state);
  void getWindsState_V(VenusAtmosphere_C* atmos, WindsState_C* state);
  void getGasesState_V(VenusAtmosphere_C* atmos, GasesState_C* state,
    ConstituentGas_C* hydrogen, ConstituentGas_C* helium, ConstituentGas_C* oxygen,
    ConstituentGas_C* nitrogen, ConstituentGas_C* dinitrogen,
    ConstituentGas_C* carbonMonoxide, ConstituentGas_C* carbonDioxide);
  void getEphemerisState_V(VenusAtmosphere_C* atmos, EphemerisState_C* state);
  void getPerturbationState_V(VenusAtmosphere_C* atmos, PerturbationState_C* state);

#ifdef __cplusplus
} // extern "C"
#endif
