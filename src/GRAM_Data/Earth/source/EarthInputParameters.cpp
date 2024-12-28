//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Earth-GRAM
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "EarthInputParameters.h"

namespace GRAM {

//! \brief The default constructor
//!
//! This basic constructor is called when objects are created without parameters.
EarthInputParameters::EarthInputParameters()
{
  // Default time = Cassini encounter date, UTC & ERT
  year = 2000;
  month = 12;
  day = 30;
  hour = 0;
  minute = 0;
  seconds = 0.0;
  timeScale = UTC;
  timeFrame = ERT;
}

//! \fn  EarthInputParameters::EarthInputParameters(const EarthInputParameters& orig)
//! \brief The copy constructor
//!
//! This constructor enables the copying of objects.  This can be done by using
//! the constructor: Object newobject(oldobject).  Or by assignment: newobject = oldobject.

//! \fn  EarthInputParameters::operator=(const EarthInputParameters& orig)
//! \brief The assignment operator
//!
//! This operator enables the copying of objects by assignment: newobject = oldobject.

//! \fn  EarthInputParameters::~EarthInputParameters()
//! \brief Destructor.
//!
//! The destructor is called when an object is deleted or goes out of scope.

} // namespace
