//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "VenusCommon.h"

namespace GRAM {

//! \brief The default constructor
//!
//! This basic constructor is called when objects are created without parameters.
VenusCommon::VenusCommon(Atmosphere* atmos)
  : atmosphere(*atmos)
{
  atmosphere.setPlanetaryConstants(muDefault, equatorialRadiusDefault, polarRadiusDefault,
	  J2Default, periodDefault, specificHeatRatioDefault);
  atmosphere.setGasConstant(HELIUM, amwHelium);
  atmosphere.setGasConstant(HYDROGEN, amwHydrogen);
  atmosphere.setGasConstant(NITROGEN, amwNitrogen);
  atmosphere.setGasConstant(DINITROGEN, amwDinitrogen);
  atmosphere.setGasConstant(OXYGEN, amwOxygen);
  atmosphere.setGasConstant(CARBON_MONOXIDE, amwCarbonMonoxide);
  atmosphere.setGasConstant(CARBON_DIOXIDE, amwCarbonDioxide);
}

//! \fn  VenusCommon::VenusCommon(const VenusCommon& orig)
//! \brief The copy constructor
//!
//! This constructor enables the copying of objects.  This can be done by using
//! the constructor: Object newobject(oldobject).  Or by assignment: newobject = oldobject.

//! \fn  VenusCommon::~VenusCommon()
//! \brief Destructor.
//!
//! The destructor is called when an object is deleted or goes out of scope.

} // namespace
