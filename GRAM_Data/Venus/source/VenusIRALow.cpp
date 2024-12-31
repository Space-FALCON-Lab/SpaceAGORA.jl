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
#include "VenusIRALow.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
VenusIRALow::VenusIRALow()
: Atmosphere(), VenusCommon(this)
{
  initializeData();
}

//! \fn  VenusIRALow::VenusIRALow(const VenusIRALow& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  VenusIRALow::~VenusIRALow()
//! \copydoc Atmosphere::~Atmosphere()


//! \brief Returns heights near the lower and upper limits of the data.
//!
//! \retval a Height at first data index (lowest).
//! \retval b Height at second data index.
//! \retval y Height at penultimate data index.
//! \retval z Height at last data index (highest).
void VenusIRALow::getInterpolationLimits(greal& a, greal& b, greal& y, greal& z)
{
  a = zlo[0];
  b = zlo[1];
  y = zlo[heightSize - 2];
  z = zlo[heightSize - 1];
}

//! \brief Finds and sets the lower bound index for the current height.
//!
//! This routine finds the interpolation segment in the VIRA Low data based
//! on the current height (so setPostion() should be called first).  It then
//! sets the internal height index to the lower bound of that segment.
//!
//! \b Inputs
//! \arg position The position must be set before calling this function.
//!
//! \returns The height at the lower bound index.
greal VenusIRALow::setLowerHeightIndex()
{
  // Find height index for interpolation
  size_t j;
  for (j = 1; j < heightSize; ++j) {
    if (height <= zlo[j]) {
      break;
    }
  }

  // Use the lower bound index 
  heightIndex = j - 1;

  // Return the height at that index
  return zlo[heightIndex];
}

//! \brief Finds and sets the upper bound index for the current height.
//!
//! This routine finds the interpolation segment in the VIRA Low data based
//! on the current height (so setPostion() should be called first).  It then
//! sets the internal height index to the upper bound of that segment.
//!
//! \b Inputs
//! \arg position The position must be set before calling this function.
//!
//! \returns The height at the upper bound index.
greal VenusIRALow::setUpperHeightIndex()
{
  // If the current heightIndex is the lower bound index ...
  if (height > zlo[heightIndex] && height <= zlo[heightIndex + 1]) {
    // then just increment the index to the upper bound.
    heightIndex += 1;
  }
  else {
    // Find height index for interpolation
    size_t j;
    for (j = 1; j < heightSize; ++j) {
      if (height <= zlo[j]) {
        break;
      }
    }

    // Use the upper bound index
    heightIndex = min(j, heightSize - 1);
  }

  // Return the height at that index
  return zlo[heightIndex];
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
void VenusIRALow::update()
{
  size_t i = heightIndex;
  size_t j,jm;
  // Latitudes for low-altitude VIRA data

  greal xlat = fabs(latitude);
  // TROPICAL: Use VIRA data for lat=0-30 if abs(lat)<30
  if (xlat <= vlat[0]) {
    pressure = plo[0][i];
    density = dlo[0][i];
    temperature = tlo[0][i];
    compressibilityFactor = zetalo[0][i];
    specificGasConstant = Rlo[0][i];
  }    
  // POLAR: Use VIRA data for lat=85-90 if abs(lat)>85
  else if (xlat >= vlat[4]) {
    pressure = plo[4][i];
    density = dlo[4][i];
    temperature = tlo[4][i];
    compressibilityFactor = zetalo[4][i];
    specificGasConstant = Rlo[4][i];
  } 
  // TEMPERATE: Latitude is betwen 30 and 85.
  else {
    // Find latitude interpolation index j
    for (j = 1; j < 5; ++j) {
      if (xlat <= vlat[j]) {
        break;
      }
    }
    jm = j - 1;
    // Get latitude interpolation factor
    Interpolator latTerp;
    latTerp.makeFraction(vlat[jm], vlat[j], xlat);
    // Logarithmic interpolation on pressure
    pressure = latTerp.log(plo[jm][i], plo[j][i]);

    // Linear interpolation on temperature, zeta, and gas constant
    temperature = latTerp.linear(tlo[jm][i], tlo[j][i]);
    compressibilityFactor = latTerp.linear(zetalo[jm][i], zetalo[j][i]);
    specificGasConstant = latTerp.linear(Rlo[jm][i], Rlo[j][i]);

    // Density from perfect gas law (with compressibility zeta)
    density = pressure / (compressibilityFactor * specificGasConstant * temperature);
  }

  greal AM = UNIVERSAL_GAS * 1000.0 / specificGasConstant;
  // Use uniform mixing ratios of CO2&N2 for heights up to 82 km
  greal fco2, fn2, fo, fco;
  if (i < 72) { // was 73 in fortran
    fco2 = 0.965;
    fn2 = 0.035;
    fco = 0.0;
    fo = 0.0;
  } else {
    // Use height-variable CO2/N2/O/CO mixing ratios from 82 to 100 km
    greal AMx = min(AM, (greal)43.44);
    // Get mixing ratios from mean molecular weight
    fco2 = -1.723936 + 6.19e-2 * AMx;
    fn2 = 2.515424 - 5.71e-2 * AMx;
    fo = max(0.0, 3.4752e-2 - 8.0e-4 * AMx);
    fco = max(0.0, 0.17376 - 4.0e-3 * AMx);
  }
  // Avogadro's number
  greal AVn = AVOGADRO * 1000.0;
  // Get number densities from mixing ratios and mass density
  carbonDioxide.numberDensity = fco2 * density * AVn / AM;
  dinitrogen.numberDensity = fn2 * density * AVn / AM;
  carbonMonoxide.numberDensity = fco * density * AVn / AM;
  oxygen.numberDensity = fo * density * AVn / AM;
  helium.numberDensity = 0.0;
  nitrogen.numberDensity = 0.0;
  hydrogen.numberDensity = 0.0;
}

} // namespace
