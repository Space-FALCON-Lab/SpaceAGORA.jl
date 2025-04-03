//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Titan-GRAM
//
// Adapted from Titan-GRAM 2004 developed by Dr. C. G. (Jere) Justus.
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "TitanGCMMid.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
TitanGCMMid::TitanGCMMid(TitanGCMTerp& terp)
: TitanCommon(this), gcmTerp(terp)
{
}

//! \fn  TitanGCMMid::TitanGCMMid(const TitanGCMMid& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  TitanGCMMid::~TitanGCMMid()
//! \copydoc Atmosphere::~Atmosphere()

//! \brief Function for temperature versus altitude in region between top
//! of GCM data (at 0.3 mbar) and bottom of exosphere (at 925 km).
//!
//! Function satisfies constraints that:
//!    -# T = T03mb at height z03mb, with dT/dz = 0 at height z03mb,
//!    -# T = Texos at 925 km, with dT/dz = 0 at this altitude, and
//!    -# T has a minimum value of Texos-40 at altitude 565 km.
//! These conditions make the middle-atmosphere profiles similar to those
//! of the Yelle et al. profiles in ESA SP-1177.
//! \param z      Height in km. Should be between z03mb and 925 km.
//! \param z03mb  Height with pressure of 0.3 mbar.
//! \param T03mb  Temperature at z03mb.
//! \param Texos
//! \returns Temperature
greal TitanGCMMid::getMiddleAtmosphereTemperature(greal z, greal z03mb, greal T03mb, greal Texos)
{
  // Wavelength for sinusoidal variation between z03mb and 565 km
  greal waveLength = 2.0 * (565.0 - z03mb);

  // Average temperature & average height between 0.3 mbar and 565 km
  greal Tb = 0.5 * (T03mb + Texos - 40.0);
  greal zb = 0.5 * (z03mb + 565.0);

  // Amplitude for sinusoidal term between 0.3 mbar level and 565 km
  greal A = 0.5 * (T03mb - Texos + 40.0);

  greal T = 0;
  if (z >= 565.0) {
    // Sinusoidal term for layer between 565 and 925 km
    T = Texos + 20.0 * (sin(PI * (z - 745.0) / 360.0) - 1.0);
  }
  else {
    // Different sinusoidal term for layer between 0.3 mbar & 565 km
    T = Tb - A * sin(2.0 * PI * (z - zb) / waveLength);
  }
  return T;
}

//! \brief Compute methane and argon mole fractions at given height.
//!
//! Yelle et al. functions for methane and argon mole fractions (ESA
//! SP 1177, page 245), adjusted to give provided values of mole
//! fractions at 0.3 mbar level (fch403mb and far03mb)
//! \param z              Height.
//! \param z03mb          Height at pressure of 0.3 mbar.
//! \param mfMethane03mb  Methane mole fraction at pressure of 0.3 mbar.
//! \param mfArgon03mb    Argon mole fraction at pressure of 0.3 mbar.
//! \param[out] mfMethane A greal.
//! \param[out] mfArgon   A greal.
//! \retval mfMethane     Methane mole fraction at height z.
//! \retval mfArgon       Argon mole fraction at height z.
void TitanGCMMid::getMoleFractions(greal z, greal z03mb, greal mfMethane03mb, greal mfArgon03mb, greal& mfMethane, greal& mfArgon)
{
  // Variables are named to match equations (4) and (5) of ESA SP 1177.
  // Mole fractions at 0.3 mbar
  greal A1 = mfMethane03mb;
  greal A3 = mfArgon03mb;

  // Yelle et al. recommended height parameter zh (km)
  const greal zh = 1050.0;

  // Titan radius (km)
  const greal Rt = 2575.5;

  // Yelle et al. recommended parameter
  const greal kappa = 0.625;
  greal omKappa = 1.0 - kappa;

  // Exponents in mole fraction functions
  greal expMethane = 3.0 / (7.0 * omKappa);
  greal expArgon = -0.3 / omKappa;

  // Height function at desired height z  (given below equation (5)).
  greal x = 1.76e5 * (z - zh) / ((Rt + zh) * (Rt + z));
  // Height function at height of 0.3 mbar level
  greal x03 = 1.76e5 * (z03mb - zh) / ((Rt + zh) * (Rt + z03mb));

  // Exponential terms for use in methane and argon functions
  greal zBase = 1.0 + exp(omKappa * x);
  greal z03Base = 1.0 + exp(omKappa * x03);

  // Terms that ensure mole fractions = desired values at 0.3 mbar
  greal A2 = A1 * (1.0 - pow(z03Base, expMethane));
  greal A4 = A3 * (1.0 - pow(z03Base, expArgon));

  // Compute methane and argon mole fractions at height z
  mfMethane = A1 * pow(zBase, expMethane) + A2;
  mfArgon = A3 * pow(zBase, expArgon) + A4;
}

//! \copydoc Atmosphere::update()
void TitanGCMMid::update()
{
  // Make sure height is below bottom of exosphere (925 km)
  if (height > maxHeight) {
    height = maxHeight;
  }

  // Get mole fractions at 0.3 mbar level (depend on input fmolmeth)
  greal z03mb = gcmTerp.getHeightAt03mb();
  greal t03mb = gcmTerp.getTemperatureAt03mb();

  // Assumptions (checking occurs in debug mode only)
  assert(height >= z03mb);

  // Use height at 0.3 mbar.
  Position p = position;
  p.height = z03mb;
  gcmTerp.setPosition(p);
  gcmTerp.setEphemerisState(ephem);
  gcmTerp.update();

  AtmosphereState atmos03mb = gcmTerp.getAtmosphereState();
  // Assumptions (checking occurs in debug mode only)
  assert(atmos03mb.totalNumberDensity != 0.0);
  greal mfMethane03mb = atmos03mb.methane.numberDensity / atmos03mb.totalNumberDensity;
  greal mfArgon03mb = atmos03mb.argon.numberDensity / atmos03mb.totalNumberDensity;
  ewWind = atmos03mb.ewWind;
  nsWind = atmos03mb.nsWind;

  // Start height steps at 0.3 mbar level for integration and
  // interpolation; initialize temperature and pressure there
  greal z1 = z03mb;
  greal pres1 = 30.0;
  greal temp1 = t03mb;

  // Get gravity at latitude and 0.3 mbar level
  greal radius = getRadius(latitude);
  greal grav1 = getGravity(latitude, radius, z1);

  // Use Yelle et al. mixing ratio functions (ESA SP-1177, page 245)
  // to get mixing ratios at height z1
  greal mfArgon, mfMethane;
  getMoleFractions(z1, z03mb, mfMethane03mb, mfArgon03mb, mfMethane, mfArgon);

  // Get average molecular weight at height z1
  greal amw1 = (1.0 - mfMethane - mfArgon) * dinitrogen.averageMolecularWeight
    + mfMethane * methane.averageMolecularWeight
    + mfArgon * argon.averageMolecularWeight;

  // Assumptions (checking occurs in debug mode only)
  assert(amw1 != 0.0);
  assert(grav1 != 0.0);
  assert(temp1 != 0.0);

  // Get pressure scale height and density at height z1
  greal psh1 = UNIVERSAL_GAS * 1000.0 * temp1 / (1000.0 * amw1 * grav1);
  greal dens1 = pres1 * amw1 / (UNIVERSAL_GAS * 1.0e3 * temp1);

  // Start loop on height in 5-km steps
  // If height z > z2=z1 + 5 km, then do hydrostatic integration to
  // compute conditions at height z2
  greal amw2 = 0, temp2 = 0, pres2 = 0, dens2 = 0, grav2 = 0, psh2 = 0, mfNitrogen = 0;
  greal z2 = z1 + 5.0;
  while (height > z2) {
    // Compute temperature and gravity at height z2
    temp2 = getMiddleAtmosphereTemperature(z2, z03mb, t03mb, exosphericTemperature);
    grav2 = getGravity(latitude, radius, z2);

    // Compute Yelle et al. mole fractions at height z2
    getMoleFractions(z2, z03mb, mfMethane03mb, mfArgon03mb, mfMethane, mfArgon);
    mfNitrogen = 1.0 - mfMethane - mfArgon;

    // Compute average molecular weight at height z2
    amw2 = mfNitrogen * dinitrogen.averageMolecularWeight 
      + mfMethane * methane.averageMolecularWeight
      + mfArgon * argon.averageMolecularWeight;

    // Assumptions (checking occurs in debug mode only)
    assert(amw2 != 0.0);
    assert(grav2 != 0.0);

    // Compute pressure scale height at height z2
    psh2 = UNIVERSAL_GAS * temp2 / (amw2 * grav2);

    // Assumptions (checking occurs in debug mode only)
    assert(temp2 != 0.0);
    assert(psh1 + psh2 != 0.0);

    // Compute pressure and density at height z2
    pres2 = pres1 * exp(-2.0 * (z2 - z1) / (psh1 + psh2));
    dens2 = pres2 * amw2 / (UNIVERSAL_GAS * 1.0e3 * temp2);

    // Store conditions at height z2 as conditions at next height z1
    // for next height step of 5 km
    z1 = z2;
    temp1 = temp2;
    grav1 = grav2;
    amw1 = amw2;
    psh1 = psh2;
    pres1 = pres2;
    dens1 = dens2;
    z2 = z1 + 5.0;
  }

  // At desired height between z1 and z2=z1+5km, do hydrostatics for
  // interpolation between z1 and z2
  // Get temperature and gravity at desired height
  temperature = getMiddleAtmosphereTemperature(height, z03mb, t03mb, exosphericTemperature);
  gravity = getGravity(latitude, radius, height);

  // Compute Yelle et al. mole fractions at height
  getMoleFractions(z2, z03mb, mfMethane03mb, mfArgon03mb, mfMethane, mfArgon);
  mfNitrogen = 1.0 - mfMethane - mfArgon;

  // Compute average molecular weight at height
  averageMolecularWeight = mfNitrogen * dinitrogen.averageMolecularWeight
    + mfMethane * methane.averageMolecularWeight 
    + mfArgon * argon.averageMolecularWeight;

  // Assumptions (checking occurs in debug mode only)
  assert(averageMolecularWeight != 0.0);
  assert(gravity != 0.0);
  assert(temperature != 0.0);

  // Compute pressure scale height at height
  psh2 = UNIVERSAL_GAS * temperature / (averageMolecularWeight * gravity);

  // Get average scale height over level z1 to height
  pressureScaleHeight = (psh1 + psh2) / 2.0;

  // Assumptions (checking occurs in debug mode only)
  assert(pressureScaleHeight != 0.0);

  // Compute pressure and density at height
  pressure = pres1 * exp(-1.0 * (height - z1) / pressureScaleHeight);
  density = pressure * averageMolecularWeight / (UNIVERSAL_GAS * 1.0e3 * temperature);

  // Assumptions (checking occurs in debug mode only)
  assert(density != 0.0);

  // Compute average density scale height over level z1 to height
  if (height == z1) {
    densityScaleHeight = pressureScaleHeight;
  } 
  else {
    densityScaleHeight = (height - z1) / log(dens1 / density);
  }

  // Specific Gas Constant at height
  specificGasConstant = pressure / (density * temperature);

  // Compute species number densities, total number density, and
  // average molecular weight
  dinitrogen.numberDensity = mfNitrogen * density * (AVOGADRO * 1.0e3) / averageMolecularWeight;
  methane.numberDensity = mfMethane * density * (AVOGADRO * 1.0e3) / averageMolecularWeight;
  argon.numberDensity = mfArgon * density * (AVOGADRO * 1.0e3) / averageMolecularWeight;
  totalNumberDensity = dinitrogen.numberDensity + methane.numberDensity + argon.numberDensity;
 }

} // namespace
