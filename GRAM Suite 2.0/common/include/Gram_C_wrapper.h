//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

// This is the C to C++ implementation layer for the GRAM C interface.
// It is separated from the GRAM_C.h file for documentation purposes.

#ifdef __cplusplus
#include "PerturbedAtmosphere.h"

namespace GRAM {

typedef void* (*CreateFunction)();
typedef void* (*CopyFunction)(void*);
typedef int (*UpdateFunction)(void*);
typedef void (*DeleteFunction)(void*);
void registerBody_C(GRAM_BODY gramBody, CreateFunction createFunction, CopyFunction copyFunction, UpdateFunction updateFunction, DeleteFunction deleteFunction);

void setStartTime_C(PerturbedAtmosphere* atmos, const GramTime_C* time);
void setPosition_C(PerturbedAtmosphere* atmos, const Position_C* position);
void setDelta_C(PerturbedAtmosphere* atmos, const Position_C* position);
void setSeed_C(PerturbedAtmosphere* atmos, int seed);
void setMinRelativeStepSize_C(PerturbedAtmosphere* atmos, greal minRelativeStepSize);
void setPerturbationScales_C(PerturbedAtmosphere* atmos, greal densityScale, greal ewWindScale, greal nsWindScale, greal verticalWindScale);
void addAuxiliaryAtmosphere_C(PerturbedAtmosphere* adapter, const char* fileName, greal innerRadius, greal outerRadius, int isEastLongitudePositive);
void setAuxiliaryValues_C(PerturbedAtmosphere* adapter, greal dens, greal pres, greal temp, greal ew, greal ns);
void setPerturbationAction_C(PerturbedAtmosphere* atmos, int action);
void setEphemerisState_C(PerturbedAtmosphere* atmos, const EphemerisState_C* state);
void setEphemerisFastModeOn_C(PerturbedAtmosphere* atmos, int flag);
void setSubsolarUpdateTime_C(PerturbedAtmosphere* atmos, greal utime);

void getStartTime_C(PerturbedAtmosphere* atmos, GramTime_C* time);
void getPosition_C(PerturbedAtmosphere* atmos, Position_C* position);
void getDynamicsState_C(PerturbedAtmosphere* atmos, DynamicsState_C* state);
void getDensityState_C(PerturbedAtmosphere* atmos, DensityState_C* state);
void getWindsState_C(PerturbedAtmosphere* atmos, WindsState_C* state);
void getGasesState_C(PerturbedAtmosphere* atmos, GasesState_C* state);
void getConstituentGas_C(PerturbedAtmosphere* atmos, ConstituentGas_C* gas, GasType type);
void getEphemerisState_C(PerturbedAtmosphere* atmos, EphemerisState_C* state);
void getPerturbationState_C(PerturbedAtmosphere* atmos, PerturbationState_C* state);

extern std::string errorMessage;
}
#endif

