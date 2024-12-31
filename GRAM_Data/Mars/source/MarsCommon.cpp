//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "MarsCommon.h"

namespace GRAM {

std::string MarsCommon::dataPath;

//! \brief The default constructor
//!
//! This basic constructor is called when objects are created without parameters.
MarsCommon::MarsCommon(Atmosphere* atm)
  : atmosphere(atm)
{
  if (atmosphere) {
    atmosphere->setPlanetaryConstants(muDefault, equatorialRadiusDefault, polarRadiusDefault,
		J2Default, periodDefault, specificHeatRatioDefault);
    atmosphere->setGasConstant(DINITROGEN, amwDinitrogen);
    atmosphere->setGasConstant(ARGON, amwArgon);
    atmosphere->setGasConstant(DIOXYGEN, amwDioxygen);
    atmosphere->setGasConstant(CARBON_DIOXIDE, amwCarbonDioxide);
    atmosphere->setGasConstant(CARBON_MONOXIDE, amwCarbonMonoxide);
    atmosphere->setGasConstant(OXYGEN, amwOxygen);
    atmosphere->setGasConstant(HELIUM, amwHelium);
    atmosphere->setGasConstant(DIHYDROGEN, amwDihydrogen);
    atmosphere->setGasConstant(HYDROGEN, amwHydrogen);
    atmosphere->setGasConstant(WATER, amwWater);
  }
}

//! \fn  MarsCommon::MarsCommon(const MarsCommon& orig)
//! \brief The copy constructor
//!
//! This constructor enables the copying of objects.  This can be done by using
//! the constructor: Object newobject(oldobject).  Or by assignment: newobject = oldobject.

//! \fn  MarsCommon::~MarsCommon()
//! \brief Destructor.
//!
//! The destructor is called when an object is deleted or goes out of scope.

} // namespace
