//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "MarsInputParameters.h"

namespace GRAM {

//! \brief The default constructor
//!
//! This basic constructor is called when objects are created without parameters.
MarsInputParameters::MarsInputParameters()
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

//! \fn  MarsInputParameters::MarsInputParameters(const MarsInputParameters& orig)
//! \brief The copy constructor
//!
//! This constructor enables the copying of objects.  This can be done by using
//! the constructor: Object newobject(oldobject).  Or by assignment: newobject = oldobject.

//! \fn  MarsInputParameters::~MarsInputParameters()
//! \brief Destructor.
//!
//! The destructor is called when an object is deleted or goes out of scope.

} // namespace
