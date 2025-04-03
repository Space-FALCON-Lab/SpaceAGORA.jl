//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Uranus-GRAM
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "UranusCommon.h"

namespace GRAM {

//! \brief The default constructor
//!
//! This basic constructor is called when objects are created without parameters.
UranusCommon::UranusCommon(Atmosphere* atmos)
  : atmosphere(*atmos)
{
  atmosphere.setPlanetaryConstants(muDefault, equatorialRadiusDefault, polarRadiusDefault,
	  J2Default, periodDefault, specificHeatRatioDefault);
  atmosphere.setGasConstant(HELIUM, amwHelium);
  atmosphere.setGasConstant(DIHYDROGEN, amwDihydrogen);
  atmosphere.setGasConstant(METHANE, amwMethane);
}

//! \fn  UranusCommon::UranusCommon(const UranusCommon& orig)
//! \brief The copy constructor
//!
//! This constructor enables the copying of objects.  This can be done by using
//! the constructor: Object newobject(oldobject).  Or by assignment: newobject = oldobject.

//! \fn  UranusCommon::~UranusCommon()
//! \brief Destructor.
//!
//! The destructor is called when an object is deleted or goes out of scope.

} // namespace
