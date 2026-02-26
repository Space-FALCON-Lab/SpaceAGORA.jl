//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include "EphemerisState.h"

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
EphemerisState::EphemerisState()
{
}

//! \fn EphemerisState::EphemerisState(const EphemerisState& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn EphemerisState::~EphemerisState()
//! \copydoc Atmosphere::~Atmosphere()

//! \brief Resets values to zero.
void EphemerisState::reset()
{
  solarTime = 0.0;        
  longitudeSun = 0.0;     
  subsolarLatitude = 0.0; 
  subsolarLongitude = 0.0;
  orbitalRadius = 0.0;    
  oneWayLightTime = 0.0;  
  solarZenithAngle = 0.0;
  secondsPerSol = 0.0;
}

//! \fn EphemerisState::getSubsolarLongitude(bool eastPositive) const
//! \brief Get the subsolar longitude in an east/west positive convention.
//!
//! \param eastPositive If true value is east postive, else it is west positive.
//! \returns The subsolar longitude in degrees.

} // namespace
