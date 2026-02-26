//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Titan-GRAM
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "TitanInputParameters.h"

namespace GRAM {

//! \brief The default constructor
//!
//! This basic constructor is called when objects are created without parameters.
TitanInputParameters::TitanInputParameters()
{
  // Default time = Voyager 1 encounter date, UTC & ERT
  year = 1980;
  month = 11;
  day = 12;
  hour = 0;
  minute = 0;
  seconds = 0.0;
  timeScale = UTC;
  timeFrame = ERT;
}

//! \fn  TitanInputParameters::TitanInputParameters(const TitanInputParameters& orig)
//! \brief The copy constructor
//!
//! This constructor enables the copying of objects.  This can be done by using
//! the constructor: Object newobject(oldobject).  Or by assignment: newobject = oldobject.

//! \fn  TitanInputParameters::~TitanInputParameters()
//! \brief Destructor.
//!
//! The destructor is called when an object is deleted or goes out of scope.

} // namespace
