//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Jupiter-GRAM
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "JupiterCommon.h"

namespace GRAM {

//! \brief The default constructor
//!
//! This basic constructor is called when objects are created without parameters.
JupiterCommon::JupiterCommon(Atmosphere* atmos)
  : atmosphere(*atmos)
{
  atmosphere.setPlanetaryConstants(muDefault, equatorialRadiusDefault, polarRadiusDefault,
	  J2Default, periodDefault, specificHeatRatioDefault);
}

//! \fn  JupiterCommon::JupiterCommon(const JupiterCommon& orig)
//! \brief The copy constructor
//!
//! This constructor enables the copying of objects.  This can be done by using
//! the constructor: Object newobject(oldobject).  Or by assignment: newobject = oldobject.

//! \fn  JupiterCommon::~JupiterCommon()
//! \brief Destructor.
//!
//! The destructor is called when an object is deleted or goes out of scope.

} // namespace
