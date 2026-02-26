//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Titan-GRAM
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "TitanCommon.h"

namespace GRAM {

//! \brief The default constructor
//!
//! This basic constructor is called when objects are created without parameters.
TitanCommon::TitanCommon(Atmosphere* atmos)
  : atmosphere(*atmos)
{
  atmosphere.setPlanetaryConstants(muDefault, equatorialRadiusDefault, polarRadiusDefault,
	  J2Default, periodDefault, specificHeatRatioDefault);
  atmosphere.setGasConstant(ARGON, amwArgon);
  atmosphere.setGasConstant(DINITROGEN, amwDinitrogen);
  atmosphere.setGasConstant(METHANE, amwMethane);
}

//! \fn  TitanCommon::TitanCommon(const TitanCommon& orig)
//! \brief The copy constructor
//!
//! This constructor enables the copying of objects.  This can be done by using
//! the constructor: Object newobject(oldobject).  Or by assignment: newobject = oldobject.

//! \fn  TitanCommon::~TitanCommon()
//! \brief Destructor.
//!
//! The destructor is called when an object is deleted or goes out of scope.

} // namespace
