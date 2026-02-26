//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Neptune-GRAM
//
// Adapted from Neptune-GRAM 2004 developed by Dr. C. G. (Jere) Justus.
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include "NeptuneMinMax.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
NeptuneMinMax::NeptuneMinMax()
: MinMaxModel(), NeptuneCommon(this)
{
  initializeGases();
  initializeData();
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
NeptuneMinMax::NeptuneMinMax(const NeptuneMinMax& orig)
  : MinMaxModel(orig), NeptuneCommon(this)
{
  initializeGases();
  fmolnitro = orig.fmolnitro;
}

//! \fn  NeptuneMinMax::~NeptuneMinMax()
//! \copydoc Atmosphere::~Atmosphere()

//! \copydoc Atmosphere::update()
void NeptuneMinMax::update()
{
  updateMinMaxFactor();
  updateAtmosphereState();
  updateWinds();
  atmos.minMaxFactor = minMaxFactor;
}

//!  \brief Calculate the effects of season, latitude, and local solar time.
//!
//! If the #minMaxFactor has not been set by the user, then the factor is
//! updated based on its current value and input parameters.
//!
//! \b Inputs
//! \arg #latitude     
//! \arg #longitudeSun 
//! \arg #solarTime    
//! \arg #computeMinMaxFactor 
//! \arg #userMinMaxFactor   
//!
//! \retval #minMaxFactor 
void NeptuneMinMax::updateMinMaxFactor()
{
  // Inputs are computeMinMaxFactor = true to compute longitudeSun, latitude, and solarTime effect,
  // or computeMinMaxFactor = false to use input value (userMinMaxFactor)
  // If evaluated, the input value massDensityFactor
  // is re-scaled by a factor of 1 minus magnitude of Ls-Lat-LST effect,
  // so that final calculated value cannot exceed range +/- 1
  if (computeMinMaxFactor) {
    greal alat = 0.44;  // *** SOURCE UNKNOWN ***
    greal als = 0.2;   // *** SOURCE UNKNOWN ***
    greal alst = 0.2;   // *** SOURCE UNKNOWN ***
    greal latr = toRadians(latitude);
    minMaxFactor = 
      alat * cos(4.0 * latr) * (1.0 + als * sin(toRadians(longitudeSun)) * sin(latr))
      - alst * cos(toRadians(15.0_deg * solarTime))
      + userMinMaxFactor * (1.0 - alat * (1.0 + als) - alst);
  }
  else {
    minMaxFactor = userMinMaxFactor;
  }
}

//! \brief Compute zonal and meridional winds at height and latitude.
//!
//! Zonal wind model from ESA SP-1177, page 140.
//! Latitude variation of zonal wind from equation (2) page 143 of ESA SP-1177.
//!
//! \b Inputs
//! \arg #height  
//! \arg #latitude  
//!
//! \retval #ewWind 
//! \retval #nsWind 
void NeptuneMinMax::updateWinds()
{
  // Compute zonal and meridional winds (u and v) at height and latitude

  // Zonal wind model from "Neptune and Triton" pages 626-628, and 636
  // zonal velocity at the equator (U sub e)
  greal U_e = -400.0;
  if (height > 0.0) {
    U_e += 0.7 * height;
  }
  if (U_e > 0.0) {
    U_e = 0.0;
  }

  // Latitude variation of zonal wind from equation (2) page 143 of 
  // "Huygens Science, Payload and Mission", ESA SP - 1177, fit to Neptune data
  greal coslat = cos(toRadians(latitude));
  // Latitude variation exponent (2 / R_i - 1)   *** SOURCE UNKNOWN ***
  const greal xpon = 0.63;
  // Omega_a in m/s
  const greal Omega_a = 2683.0;
  ewWind = (U_e + Omega_a) * pow(coslat, xpon) - Omega_a * coslat;
  nsWind = 0.0;
}

//! \brief Set the molecular nitrogen (N2) mole fraction.
//!
//! Diatomic nitrogen (N2) quantities will not be computed unless the 
//! mole fraction is set with this function prior to an update.
//! \param n2mf The mole fraction.
void NeptuneMinMax::setDinitrogenMoleFraction(greal n2mf)
{
  // Activate dinitrogen
  if (n2mf > 0.0) {
    useDinitrogen();
  }
  // Set the mole fraction.
  fmolnitro = clamp(n2mf, 0.0, 0.006);
}

//! \brief Compute the mole fraction of each constituent gases.
//!
//! The routine computes the mole fraction of each of the active gases.
//!
//! \b Inputs
//! \arg #methane 
//! \arg #dinitrogen 
//! \arg #dihydrogen 
//! \arg #helium 
//! \arg #fmolnitro 
//!
//! \retval #gases 
void NeptuneMinMax::updateMoleFractions()
{
  // Assumptions (checking occurs in debug mode only)
  assert(totalNumberDensity != 0.0);
  assert(helium.averageMolecularWeight != 0.0);
  assert(dihydrogen.averageMolecularWeight != 0.0);

  // Methane number density will be held constant
  greal fmolCH4 = methane.numberDensity / totalNumberDensity;
  greal fmolN2 = 0.0;
  // Recompute number densities if input fmolnitro > 0
  if (dinitrogen.isPresent) {
    // Adjust N2 mole fraction for height
    fmolN2 = n2mixr(height, fmolnitro);

    // Get H2 and HE mole fractions from number densities.
    greal fmolH2 = 100.0 * dihydrogen.numberDensity / totalNumberDensity;
    greal fmolHe = 100.0 * helium.numberDensity / totalNumberDensity;

    // Given N2 and CH4, the rest must be H2 and HE
    greal fmolH2HE = 1.0 - fmolN2 - fmolCH4;
    greal amwH2HE = averageMolecularWeight 
      - dinitrogen.averageMolecularWeight * fmolN2 
      - methane.averageMolecularWeight * fmolCH4;

    // Preserve mean molecular weight
    fmolH2 = (fmolH2HE * helium.averageMolecularWeight - amwH2HE) / (helium.averageMolecularWeight - dihydrogen.averageMolecularWeight);
    fmolHe = (amwH2HE - dihydrogen.averageMolecularWeight * fmolH2HE) / (helium.averageMolecularWeight - dihydrogen.averageMolecularWeight);

    // Get number densities from recomputed mole fractions
    dinitrogen.numberDensity = fmolN2 * totalNumberDensity;
    dihydrogen.numberDensity = fmolH2 * totalNumberDensity;
    helium.numberDensity = fmolHe * totalNumberDensity;
  }
  // Update mole fractions based on number densities.
  MinMaxModel::updateMoleFractions();
}

//! \brief Adjust the dinitrogen mole fraction based on height.
//!
//! \param ht Height in km.
//! \param fmolnitro Nominal N2 mole fraction.
//! \returns Adjusted N2 mole fraction.
greal NeptuneMinMax::n2mixr(greal ht, greal fmolnitro)
{
  // Interpolates between fmolnitro at surface ...
  greal N2mr0 = fmolnitro;
  // ... and 1e-10 at z = 580.
  greal zm = 580.0;
  greal N2mrzm = 1.0e-10;

  greal al0 = log(N2mr0);
  greal alzm = log(N2mrzm);
  greal A = 2.0*alzm - al0;
  greal B = 2.0*(al0 - alzm);

  // Interpolation factor is height
  greal x = (ht - zm) / 59.0;
  greal expx = exp(-x);

  // Interpolate 
  greal alf = A + B * expx / (1.0 + expx);
  return exp(alf);
}

} // namespace