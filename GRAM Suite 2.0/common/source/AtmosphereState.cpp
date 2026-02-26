//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include "AtmosphereState.h"

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
AtmosphereState::AtmosphereState()
  : AtmosphereStateBase()
{
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
AtmosphereState::AtmosphereState(const AtmosphereState& orig)
  : AtmosphereStateBase(orig)
{
  planetSpecificMetrics = orig.planetSpecificMetrics;
}

//! \copydoc Atmosphere::~Atmosphere()
AtmosphereState::~AtmosphereState()
{
}

AtmosphereState& AtmosphereState::operator=(const AtmosphereState &value)
{
  AtmosphereStateBase::operator=(value);
  planetSpecificMetrics = value.planetSpecificMetrics;
  return *this;
}

} // namespace