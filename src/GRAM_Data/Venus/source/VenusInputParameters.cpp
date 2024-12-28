//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "VenusInputParameters.h"

namespace GRAM {

//! \brief The default constructor
//!
//! This basic constructor is called when objects are created without parameters.
VenusInputParameters::VenusInputParameters()
{
  // Default time = Pioneer probe encounter date, UTC & ERT
  year = 1978;
  month = 12;
  day = 9;
  hour = 19;
  minute = 47;
  seconds = 59.0;
  timeScale = UTC;
  timeFrame = ERT;

}

//! \fn  VenusInputParameters::VenusInputParameters(const VenusInputParameters& orig)
//! \brief The copy constructor
//!
//! This constructor enables the copying of objects.  This can be done by using
//! the constructor: Object newobject(oldobject).  Or by assignment: newobject = oldobject.

//! \fn  VenusInputParameters::~VenusInputParameters()
//! \brief Destructor.
//!
//! The destructor is called when an object is deleted or goes out of scope.

} // namespace
