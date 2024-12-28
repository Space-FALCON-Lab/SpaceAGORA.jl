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
#include "VenusIRAHigh.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
VenusIRAHigh::VenusIRAHigh()
: Atmosphere(), VenusCommon(this)
{
  initializeData();
}

//! \fn  VenusIRAHigh::VenusIRAHigh(const VenusIRAHigh& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  VenusIRAHigh::~VenusIRAHigh()
//! \copydoc Atmosphere::~Atmosphere()

//! \brief Returns heights near the lower and upper limits of the data.
//!
//! \retval a Height at first data index (lowest).
//! \retval b Height at second data index.
//! \retval y Height at penultimate data index.
//! \retval z Height at last data index (highest).
void VenusIRAHigh::getInterpolationLimits(greal& a, greal& b, greal& y, greal& z)
{
  a = zhi[0];
  b = zhi[1];
  y = zhi[heightSize - 2];
  z = zhi[heightSize - 1];
}

//! \brief Finds and sets the lower bound index for the current height.
//!
//! This routine finds the interpolation segment in the VIRA High data based
//! on the current height (so setPostion() should be called first).  It then
//! sets the internal height index to the lower bound of that segment.
//!
//! \b Inputs
//! \arg position The position must be set before calling this function.
//!
//! \returns The height at the lower bound index.
greal VenusIRAHigh::setLowerHeightIndex()
{
  // Find height index for interpolation
  size_t j;
  for (j = 1; j < heightSize; ++j) {
    if (position.height <= zhi[j]) {
      break;
    }
  }
  // Use the lower bound index 
  heightIndex = j - 1;

  // Return the height at that index
  return zhi[heightIndex];
}

//! \brief Finds and sets the upper bound index for the current height.
//!
//! This routine finds the interpolation segment in the VIRA High data based
//! on the current height (so setPostion() should be called first).  It then
//! sets the internal height index to the upper bound of that segment.
//!
//! \b Inputs
//! \arg position The position must be set before calling this function.
//!
//! \returns The height at the upper bound index.
greal VenusIRAHigh::setUpperHeightIndex()
{
  if (height > zhi[heightSize - 1]) {
    heightIndex = heightSize - 1;
  }
  // If the current heightIndex is the lower bound index ...
  else if (height > zhi[heightIndex] && height <= zhi[heightIndex + 1]) {
    // then just increment the index to the upper bound.
    heightIndex += 1;
  }
  else {
    // Find height index for interpolation
    size_t j;
    for (j = 1; j < heightSize; ++j) {
      // Find height index for interpolation
      if (position.height <= zhi[j]) {
        break;
      }
    }

    // Use the upper bound index
    heightIndex = min(j, heightSize - 1);
  }

  // Return the height at that index
  return zhi[heightIndex];
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
void VenusIRAHigh::update()
{
  size_t i = heightIndex;
  size_t j, jm;
  // Solar zenith angle high-altitude VIRA data
  greal xsza = fabs(ephem.solarZenithAngle);
  // Use VIRA data for lat=0-30 if abs(lat)<30
  if (xsza <= vsza[0]) {
    pressure = phi[0][i];
    density = dhi[0][i];
    temperature = thi[0][i];
    carbonDioxide.numberDensity = yco2hi[0][i];
    dinitrogen.numberDensity = yn2hi[0][i];
    oxygen.numberDensity = yohi[0][i];
    carbonMonoxide.numberDensity = ycohi[0][i];
    helium.numberDensity = yhehi[0][i];
    nitrogen.numberDensity = ynhi[0][i];
    hydrogen.numberDensity = yhhi[0][i];
  }    // Use VIRA data for lat=85-90 if abs(lat)>85
  else if (xsza >= vsza[6]) {
    pressure = phi[6][i];
    density = dhi[6][i];
    temperature = thi[6][i];
    carbonDioxide.numberDensity = yco2hi[6][i];
    dinitrogen.numberDensity = yn2hi[6][i];
    oxygen.numberDensity = yohi[6][i];
    carbonMonoxide.numberDensity = ycohi[6][i];
    helium.numberDensity = yhehi[6][i];
    nitrogen.numberDensity = ynhi[6][i];
    hydrogen.numberDensity = yhhi[6][i];
  } else {
    // Find latitude interpolation index j
    for (j = 1; j < szaSize; ++j) {
      if (xsza <= vsza[j]) {
        break;
      }
    }
    jm = j - 1;
    // Get the solar zenith angle interpolation factor
    Interpolator szaTerp;
    szaTerp.makeFraction(vsza[jm], vsza[j], xsza);

    // Logarithmic interpolation on pressure
    pressure = szaTerp.log(phi[jm][i], phi[j][i]);

    // Linear interpolation on temperature, zeta, and gas constant
    temperature = szaTerp.linear(thi[jm][i], thi[j][i]);

    // Gas constant 
    greal R1 = phi[jm][i] / (dhi[jm][i] * thi[jm][i]);
    greal R2 = phi[j][i] / (dhi[j][i] * thi[j][i]);
    specificGasConstant = szaTerp.linear(R1, R2);
    compressibilityFactor = 1.0;

    // Density from perfect gas law
    density = pressure / (specificGasConstant * temperature);

    // Logarithmic interpolation for number densities 
    carbonDioxide.numberDensity = szaTerp.log(yco2hi[jm][i], yco2hi[j][i]);
    dinitrogen.numberDensity = szaTerp.log(yn2hi[jm][i], yn2hi[j][i]);
    carbonMonoxide.numberDensity = szaTerp.log(ycohi[jm][i], ycohi[j][i]);
    oxygen.numberDensity = szaTerp.log(yohi[jm][i], yohi[j][i]);
    helium.numberDensity = szaTerp.log(yhehi[jm][i], yhehi[j][i]);
    nitrogen.numberDensity = szaTerp.log(ynhi[jm][i], ynhi[j][i]);
    hydrogen.numberDensity = szaTerp.log(yhhi[jm][i], yhhi[j][i]);
  }
}

} // namespace
