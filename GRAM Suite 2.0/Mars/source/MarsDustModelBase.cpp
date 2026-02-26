//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United
// States without explicit approval by NASA Marshall Space Flight Center.
//
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cmath>
#include "MarsDustModelBase.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
MarsDustModelBase::MarsDustModelBase()
  : MarsCommon(NULL)
{
}

//! \fn  MarsDustModelBase::MarsDustModelBase(const MarsDustModelBase& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  MarsDustModelBase::~MarsDustModelBase()
//! \copydoc Atmosphere::~Atmosphere()

//! \fn  MarsDustModelBase::setPosition()
//! \copydoc PerturbedAtmosphere::setPosition()

//! \fn  MarsDustModelBase::setLongitudeSun()
//! \brief Set the longitude of the sun.
//! \param lonSun The longitude of the sun (degrees east). 

//! \fn  MarsDustModelBase::getDustOpticalDepth()
//! \brief Get the dust optical depth.
//! \returns The dust optical depth.

//! \fn  MarsDustModelBase::getDustOffset()
//! \brief Get the dust height offset.
//! \returns The dust height offset \units{km}.


//! \copydoc PerturbedAtmosphere::setInputParameters()
void MarsDustModelBase::setInputParameters(const MarsInputParameters& params)
{
  stormLongitudeSun = params.stormLongitudeSun;
  stormDuration = params.stormDuration;
  stormIntensity = params.stormIntensity;
  stormMaxRadius = params.stormMaxRadius;
  stormLatitude = params.stormLatitude;
  stormLongitude = params.stormLongitude;
  if (!params.isEastLongitudePositiveOnInput) {
    stormLongitude = 360.0_deg - stormLongitude;
  }
}

//! \brief Interface for the primary dust model computations.
//!
//! Evaluates the dust optical depth and intensity factor.
//! The model must implement updateDustOpticalDepth() called by this method.
void MarsDustModelBase::update()
{
  //  Evaluate dust OD
  updateDustOpticalDepth();

  // add storm intensity to dust OD
  updateDustIntensity();
  dustOpticalDepth += intensity;
}

//! \brief Computes dust optical depth offset for dust storms.
//!
//! Computes dust optical depth offset for dust storms defined by user input parameters.
//! The dust storm starts at stormLongitudeSun and continues for stormDuration degrees.
//! The center of the dust storm is at stormLatitude, stormLongitude with a size of
//! stormMaxRadius (zero radius implies a global storm).  The intensity of the storm is
//! specified by stormIntensity.
//!
//! \b Inputs
//! \arg #position (latitude, longitude, latitudeRadius)
//! \arg #longitudeSun
//! \arg #stormLatitude
//! \arg #stormLongitude
//! \arg #stormLongitudeSun
//! \arg #stormDuration
//! \arg #stormIntensity
//! \arg #stormMaxRadius
//!
//! \retval #intensity
void MarsDustModelBase::updateDustIntensity()
{
  // Current duration (in angle) of the dust storm
  greal dlsmax = clamp(stormDuration, 12.0_deg, 48.0_deg);

  // Delta from the beginning of the storm.
  greal dls = longitudeSun - stormLongitudeSun;
  // Handle the case when Ls is near 360 degrees.
  if (stormLongitudeSun > 360.0_deg - dlsmax && longitudeSun < dlsmax) {
    dls += 360.0_deg;
  }

  // if outside storm regime (not started or ended or no storm), no intensity
  if (dls <= 0.0 || dls > dlsmax || stormIntensity <= 0.0) {
    intensity = 0.0;
    return;
  }

  // ramp up intensity (0-1)from start of storm until 1/8 of duration
  if (dls <= dlsmax / 8.0) {
    intensity = 8.0 * dls / dlsmax;
  }
  // ramp down intensity (1-0) from storm midpoint to end
  else if (dls >= dlsmax / 2.0) {
    intensity = 2.0 * (1.0 - dls / dlsmax);
  }
  // intensity=1 between 1/8 and 1/2 duration
  else {
    intensity = 1.0;
  }

  // For global storms, the default size factor is 1.
  greal sizeFactor = 1.0;

  // If regional storm and not global storm, then evaluate a size factor
  // based on distance from the storm center.
  if (abs(stormMaxRadius) > 0.0) {
    // N-S distance (km) between current position and dust storm center
    greal dns = position.latitudeRadius * toRadians(position.latitude - stormLatitude);
    // longitude delta (degrees) between current position and dust center
    greal dlon = abs(position.longitude - stormLongitude);
    // force dlon to range 0-180
    if (dlon > 180.0_deg) {
      dlon = 360.0_deg - dlon;
    }
    // E-W distance (km) between current position and dust storm center
    greal dew = position.latitudeRadius * cos(toRadians(position.latitude)) * toRadians(dlon);
    // radial distance (km) between current position and dust stormcenter
    greal rad = sqrt(pow(dns, 2) + pow(dew, 2));
    // current storm radius
    greal raddust = intensity * stormMaxRadius;
    // if radial distance < 2*current storm radius, compute size factor
    if (rad < 2.0 * raddust) {
      sizeFactor = 0.5 * (1.0 + cos(toRadians(90.0_deg * rad / raddust)));
    }
    // Too far away, zero out size factor
    else {
      sizeFactor = 0.0;
    }
  } // end if local storm

  // intensity factor is multiplied by size_factor and by input stormIntensity
  intensity *= sizeFactor * stormIntensity;
}

//! \brief Computes the dust storm offset.
//!
//! Computes the dust storm offset.
//!
//! \b Inputs
//! \arg #intensity
//! \arg #longitudeSun
//! \arg #stormLongitudeSun
//! \arg #stormDuration
//! \arg #stormIntensity
//!
//! \returns The dust storm offset.
greal MarsDustModelBase::getDustOffset() const
{
  // Protect against zero division.
  if (stormIntensity == 0.0) {
    return 0.0;
  }

  // dust storm offset
  return 5.0 * intensity * sqrt(stormIntensity) / stormIntensity;
}

} // namespace
