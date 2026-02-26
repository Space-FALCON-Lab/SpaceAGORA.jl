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
#include "MGCMInterpolator.h"

using namespace std;

namespace GRAM {

const greal MGCMInterpolator::dustOD[OD_SIZE] = { greal(0.3),  greal(1.0),  greal(3.0) };

//! \copydoc Atmosphere::Atmosphere()
MGCMInterpolator::MGCMInterpolator()
: MarsInterpolatorBase()
{
}

//! \fn  MGCMInterpolator::MGCMInterpolator(const MGCMInterpolator& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  MGCMInterpolator::~MGCMInterpolator()
//! \copydoc Atmosphere::~Atmosphere()

//! \fn  MGCMInterpolator::setDustOpticalDepth(greal depth)
//! \brief Set the dust optical depth parameter.
//! \param depth The dust optical depth.

//! \brief Update latitude, ls, and od base indices and displacements.
//!
//! This method update the base indices and displacements for latitude, ls, and optical depth.
//! Call this method from each subclass implementation of updateIndices().
//!
//! \param latitudeArraySize The size of the latitude dimension in the GCM data.
//! \param latitudeOffset Offset of latitude levels (from 0 degrees).
//! \b Inputs
//! \arg #latitude
//! \arg #longitudeSun
//! \arg #dustOpticalDepth
//!
//! \returns #baseIndex and #displacements
void MGCMInterpolator::updateBaseIndices(size_t latitudeArraySize, greal latitudeOffset)
{
  // latitude index
  greal latitudeStepSize = (180.0_deg - 2.0 * latitudeOffset) / greal(latitudeArraySize - 1);
  greal latitudeShift = 90.0_deg - latitudeOffset;
  // Need to use int as this may go negative due to the latitude offset.
  int indx = int(floor((latitude + latitudeShift) / latitudeStepSize));
  if (indx < 0) {
    baseIndex.lat = 0;
  }
  else {
    baseIndex.lat = clampSize(size_t(indx), latitudeArraySize - 1);
  }
  // Displacement for winds
  greal latdisp = (latitude + latitudeShift - latitudeStepSize * greal(baseIndex.lat)) / latitudeStepSize;
  displacements.lat = clamp(latdisp, 0.0, 1.0);

  // Compute a pole factor to dampen tidal amplitudes near the poles (used on MTGCM data).
  // Away from the poles, the factor has no effect.
  displacements.tpolefac = 1.0;
  // Near the south pole  (lowest lat layer)
  if (baseIndex.lat == 0) {
    // Polefac goes to from 0.5 to 0 at the pole.
    displacements.tpolefac = 0.5;
    if (latdisp <= 0.0) {
      displacements.tpolefac = 1.0 + (85.0 - abs(latitude)) / 5.0;
    }
  }
  // Near the north pole (uppermost lat layer)
  else if (baseIndex.lat >= latitudeArraySize - 2) {
    // Polefac goes to from 0.5 to 0 at the pole.
    displacements.tpolefac = 0.5;
    if (latdisp >= 1.0) {
      displacements.tpolefac = 1.0 - (abs(latitude) - 85.0) / 5.0;
    }
  }


  // latitude index for winds
  // Note: the MGCM wind data is offset in latitude by a half step.  
  constexpr greal wlatitudeShift = 86.25_deg; // That's 90.0_deg - 7.5_deg / 2.0

  size_t idx = int(1 + floor((latitude + wlatitudeShift) / latitudeStepSize));
  baseIndex.wlat = min(latitudeArraySize - 2, max(size_t(1), idx));
  // Displacement for winds
  greal wdisp = (latitude + wlatitudeShift - latitudeStepSize * greal(baseIndex.wlat - 1)) / latitudeStepSize;
  displacements.wlat = clamp(wdisp, 0.0, 1.0);

  // Pole factor logic for winds (for MGCM data)
  // Away from the poles, the factor has no effect.
  displacements.wpolefac = 1.0;
  // Near the south pole  (lowest wlat layer)
  if (baseIndex.wlat == 0) {
    displacements.wpolefac = 0.75;
    if (wdisp <= 0.0) {
      displacements.wpolefac = 1.0 - (85.0 - abs(latitude)) / 5.0;
    }
  }
  // Near the north pole (uppermost wlat layer)
  else if (baseIndex.wlat >= latitudeArraySize - 2) {
    displacements.wpolefac = 0.75;
    if (wdisp >= 1.0) {
      displacements.wpolefac = 1.0 - (abs(latitude) - 85.0) / 5.0;
    }
  }

  // Ls index for surface, MGCM, and MTGCM data
  idx = size_t(floor(longitudeSun / 30.0_deg));
  baseIndex.ls = clampSize(idx, LS_SIZE - 1);
  greal lsdisp = (longitudeSun - 30.0_deg * float(baseIndex.ls)) / 30.0_deg;
  displacements.ls = clamp(lsdisp, 0.0, 1.0);

  // data array index for dust optical depth
  baseIndex.od = (dustOpticalDepth < dustOD[1]) ? 0 : 1;
  greal oddisp = log(dustOpticalDepth / dustOD[baseIndex.od]) / log(dustOD[baseIndex.od + 1] / dustOD[baseIndex.od]);
  if (oddisp < 0.0 && baseIndex.od == 0) {
    // adjust displacements.od for negative values
    oddisp = (dustOpticalDepth - dustOD[0]) / (dustOD[1] - dustOD[0]);
  }
  // should not extrapolate. limit displacements.od at 1.0
  displacements.od = min(1.0, oddisp);
}

} // namespace