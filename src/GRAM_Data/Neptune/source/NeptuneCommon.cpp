//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Neptune-GRAM
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "NeptuneCommon.h"

namespace GRAM {

//! \brief The default constructor
//!
//! This basic constructor is called when objects are created without parameters.
NeptuneCommon::NeptuneCommon(Atmosphere* atmos)
  : atmosphere(*atmos)
{
  atmosphere.setPlanetaryConstants(muDefault, equatorialRadiusDefault, polarRadiusDefault,
	  J2Default, periodDefault, specificHeatRatioDefault);
  atmosphere.setGasConstant(HELIUM, amwHelium);
  atmosphere.setGasConstant(DIHYDROGEN, amwDihydrogen);
  atmosphere.setGasConstant(METHANE, amwMethane);
}

//! \fn  NeptuneCommon::NeptuneCommon(const NeptuneCommon& orig)
//! \brief The copy constructor
//!
//! This constructor enables the copying of objects.  This can be done by using
//! the constructor: Object newobject(oldobject).  Or by assignment: newobject = oldobject.

//! \fn  NeptuneCommon::~NeptuneCommon()
//! \brief Destructor.
//!
//! The destructor is called when an object is deleted or goes out of scope.

//! \brief Activates the computations for molecular nitrogen.
void NeptuneCommon::useDinitrogen()
{
  atmosphere.setGasConstant(DINITROGEN, amwDinitrogen);
}

} // namespace
