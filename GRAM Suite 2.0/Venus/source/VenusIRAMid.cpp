//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//
// Adapted from Venus-GRAM 2005 developed by Dr. C. G. (Jere) Justus.
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include "VenusIRAMid.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
VenusIRAMid::VenusIRAMid()
: Atmosphere(), VenusCommon(this)
{
  initializeData();
}

//! \fn  VenusIRAMid::VenusIRAMid(const VenusIRAMid& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  VenusIRAMid::~VenusIRAMid()
//! \copydoc Atmosphere::~Atmosphere()

//! \brief Returns heights near the lower and upper limits of the data.
//!
//! \retval a Height at first data index (lowest).
//! \retval b Height at second data index.
//! \retval y Height at penultimate data index.
//! \retval z Height at last data index (highest).
void VenusIRAMid::getInterpolationLimits(greal& a, greal& b, greal& y, greal& z)
{
  a = zmd[0];
  b = zmd[1];
  y = zmd[heightSize - 2];
  z = zmd[heightSize - 1];
}

//! \brief Finds and sets the lower bound index for the current height.
//!
//! This routine finds the interpolation segment in the VIRA Mid data based
//! on the current height (so setPostion() should be called first).  It then
//! sets the internal height index to the lower bound of that segment.
//!
//! \b Inputs
//! \arg position The position must be set before calling this function.
//!
//! \returns The height at the lower bound index.
greal VenusIRAMid::setLowerHeightIndex()
{
  // Find height index for interpolation
  size_t j;
  for (j = 1; j < heightSize; ++j) {
    if (position.height <= zmd[j]) {
      break;
    }
  }
  // Use the lower bound index 
  heightIndex = j - 1;

  // Return the height at that index
  return zmd[heightIndex];
}

//! \brief Finds and sets the upper bound index for the current height.
//!
//! This routine finds the interpolation segment in the VIRA Mid data based
//! on the current height (so setPostion() should be called first).  It then
//! sets the internal height index to the upper bound of that segment.
//!
//! \b Inputs
//! \arg position The position must be set before calling this function.
//!
//! \returns The height at the upper bound index.
greal VenusIRAMid::setUpperHeightIndex()
{
  // If the current heightIndex is the lower bound index ...
  if (heightIndex < zmd.size() - 1 && height > zmd[heightIndex] && height <= zmd[heightIndex + 1]) {
    // then just increment the index to the upper bound.
    heightIndex += 1;
  }
  else {
    // Find height index for interpolation
    size_t j;
    for (j = 1; j < heightSize; ++j) {
      // Find height index for interpolation
      if (position.height <= zmd[j]) {
        break;
      }
    }

    // Use the upper bound index
    heightIndex = min(j, heightSize - 1);
  }

  // Return the height at that index
  return zmd[heightIndex];
}

//! \brief Interface for the primary atmosphere computations.
//!
//! This routine computes the atmospheric state for the current position
//! but at the height indicated by the current height index. Either setLowerHeightIndex() or
//! setUpperHeightIndex() must be called prior to calling this function in order
//! to set the desired height index.
//!
//! \b Inputs
//! \arg position The position must be set before calling this function.
//! \arg heightIndex The height index must be set proior to calling this function.
//!
//! \returns The AtmosphereState is populated.
void VenusIRAMid::update()
{
  size_t i = heightIndex;

  // Diurnal time variation parameter
  greal ftime = sin(toRadians(45.0_deg * (solarTime - 6.0) / 3.0));

  // Latitude factor (suppresses diurnal variation near poles)
  greal flat = 1.0;
  greal alat = fabs(latitude);
  if (alat >= 70.0_deg) {
    flat = 1.0 - pow((alat - 70.0_deg) / 20.0_deg, 2.0);
  }

  // Combine factors for computational efficiency.
  greal factor = flat * ftime;

  // Time-dependent temperature
  temperature = getLinearMean(tmd[1][i], tmd[0][i], factor);

  // Time-dependent pressure
  pressure = getLogMean(pmd[1][i], pmd[0][i], factor);

  // Gas constant at LST = 0 and 12 hr
  greal R0 = pmd[0][i] / (dmd[0][i] * tmd[0][i]);
  greal R1 = pmd[1][i] / (dmd[1][i] * tmd[1][i]);
  // Time-dependent gas constant
  specificGasConstant = getLinearMean(R1, R0, factor);
  compressibilityFactor = 1.0;

  // Density from perfect gas law
  density = pressure / (specificGasConstant * temperature);

  // Time-dependent CO2 number density
  carbonDioxide.numberDensity = getLogMean(yco2md[1][i], yco2md[0][i], factor);

  // Time-dependent N2 number density
  dinitrogen.numberDensity = getLogMean(yn2md[1][i], yn2md[0][i], factor);

  // Time-dependent O number density
  oxygen.numberDensity = getLogMean(yomd[1][i], yomd[0][i], factor);

  // Time-dependent CO number density
  carbonMonoxide.numberDensity = getLogMean(ycomd[1][i], ycomd[0][i], factor);

  // Time-dependent He number density
  helium.numberDensity = getLogMean(yhemd[1][i], yhemd[0][i], factor);

  // Time-dependent N number density
  nitrogen.numberDensity = getLogMean(ynmd[1][i], ynmd[0][i], factor);

  // No hydrogen.
  hydrogen.numberDensity = 0.0;
}

//! \brief Computes diurnal mean (linear) with variation due to time
//!
//! \param noonValue The value at local solar time = 12.
//! \param midnightValue The value at local solar time = 0.
//! \param amplitudeFactor Variational factor applied to the amplitude.
//! \returns the Mean value with diurnal variation.
greal VenusIRAMid::getLinearMean(greal noonValue, greal midnightValue, greal amplitudeFactor)
{
  // Mean and diurnal amplitude.
  greal mean = 0.5 * (noonValue + midnightValue);
  greal amplitude = 0.5 * (noonValue - midnightValue);
  // Apply amplitude factor.
  return mean + amplitude * amplitudeFactor;
}

//! \brief Computes diurnal mean (log) with variation due to time
//!
//! \param noonValue The value at local solar time = 12.
//! \param midnightValue The value at local solar time = 0.
//! \param amplitudeFactor Variational factor applied to the amplitude.
//! \returns the Mean value with diurnal variation.
greal VenusIRAMid::getLogMean(greal noonValue, greal midnightValue, greal amplitudeFactor)
{
  // Mean and diurnal amplitude.
  greal mean = 0.5 * (log(noonValue * midnightValue));
  greal amplitude = 0.5 * (log(noonValue / midnightValue));
  // Apply amplitude factor.
  return exp(mean + amplitude * amplitudeFactor);
}

} // namespace
