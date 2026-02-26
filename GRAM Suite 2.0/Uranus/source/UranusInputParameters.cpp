//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Uranus-GRAM
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "UranusInputParameters.h"

namespace GRAM {

//! \brief The default constructor
//!
//! This basic constructor is called when objects are created without parameters.
UranusInputParameters::UranusInputParameters()
{
  // Default time = Voyager 2 encounter date, UTC & ERT
  year = 1986;
  month = 1;
  day = 22;
  hour = 0;
  minute = 0;
  seconds = 0.0;
  timeScale = UTC;
  timeFrame = ERT;
}

//! \fn  UranusInputParameters::UranusInputParameters(const UranusInputParameters& orig)
//! \brief The copy constructor
//!
//! This constructor enables the copying of objects.  This can be done by using
//! the constructor: Object newobject(oldobject).  Or by assignment: newobject = oldobject.

//! \fn  UranusInputParameters::~UranusInputParameters()
//! \brief Destructor.
//!
//! The destructor is called when an object is deleted or goes out of scope.

} // namespace
