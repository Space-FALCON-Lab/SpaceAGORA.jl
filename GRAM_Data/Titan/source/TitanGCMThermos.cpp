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
#include "TitanGCMThermos.h"

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
TitanGCMThermos::TitanGCMThermos(TitanGCMTerp& gcmTerp)
: TitanCommon(this), gcmMid(gcmTerp)
{
}

//! \fn  TitanGCMThermos::TitanGCMThermos(const TitanGCMThermos& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  TitanGCMThermos::~TitanGCMThermos()
//! \copydoc Atmosphere::~Atmosphere()

//! \copydoc Atmosphere::update()
void TitanGCMThermos::update()
{
  // Get gravity (gref), radius (Rf) at surface and given latitude
  greal radius= getRadius(latitude);
  greal gravityAtSurface = getGravity(latitude, radius, 0.0);

  // Compute number densities at height zf
  Position p = position;
  p.height = baseHeight;
  gcmMid.setPosition(p);
  gcmMid.setEphemerisState(ephem);
  gcmMid.setExosphericTemperature(exosphericTemperature);
  gcmMid.update();
  AtmosphereState atmosAtBaseHeight = gcmMid.getAtmosphereState();
  greal n2ndBase = atmosAtBaseHeight.dinitrogen.numberDensity;
  greal ch4ndBase = atmosAtBaseHeight.methane.numberDensity;
  greal arndBase = atmosAtBaseHeight.argon.numberDensity;

  ewWind = atmosAtBaseHeight.ewWind;
  nsWind = atmosAtBaseHeight.nsWind;

  // Temperature = exospheric temperature (constant with height,
  // dependent on latitude, local time, and seasonal Ls angle)
  temperature = exosphericTemperature;

  // Constituent partial pressures at height zf
  greal ppn2Base = BOLTZMANN * n2ndBase * temperature;
  greal ppch4Base = BOLTZMANN * ch4ndBase * temperature;
  greal pparBase = BOLTZMANN * arndBase * temperature;

  // Assumptions (checking occurs in debug mode only)
  assert(radius != 0.0);
  assert(height != 1.0);
  assert(baseHeight != 1.0);
  assert(dinitrogen.averageMolecularWeight != 0.0);
  assert(methane.averageMolecularWeight != 0.0);
  assert(argon.averageMolecularWeight != 0.0);
  assert(gravityAtSurface != 0.0);
  assert(temperature != 0.0);

  // Height variation parameter y
  greal y = (height / (1.0 + height / radius))-(baseHeight / (1.0 + baseHeight / radius));

  // Scale height, partial pressure, and number density for N2
  greal Hn2 = UNIVERSAL_GAS * temperature / (dinitrogen.averageMolecularWeight * gravityAtSurface);
  assert(Hn2 != 0.0);
  greal ppn2 = ppn2Base * exp(-y / Hn2);
  dinitrogen.numberDensity = ppn2 / (BOLTZMANN * temperature);

  // Scale height, partial pressure, and number density for CH4
  greal Hch4 = UNIVERSAL_GAS * temperature / (methane.averageMolecularWeight * gravityAtSurface);
  assert(Hch4 != 0.0);
  greal ppch4 = ppch4Base * exp(-y / Hch4);
  methane.numberDensity = ppch4 / (BOLTZMANN * temperature);

  // Scale height, partial pressure, and number density for argon
  greal Har = UNIVERSAL_GAS * temperature / (argon.averageMolecularWeight * gravityAtSurface);
  assert(Har != 0.0);
  greal ppar = pparBase * exp(-y / Har);
  argon.numberDensity = ppar / (BOLTZMANN * temperature);

  // Total pressure = sum of partial pressures
  pressure = ppn2 + ppch4 + ppar;

  // Total number density = sum of number densities
  updateTotalNumberDensity();

  // Assumptions (checking occurs in debug mode only)
  assert(totalNumberDensity != 0.0);

  // Average molecular weight
  averageMolecularWeight = getTotalMass() / totalNumberDensity;

  // Assumptions (checking occurs in debug mode only)
  assert(averageMolecularWeight != 0.0);

  // Density from gas law
  specificGasConstant = (UNIVERSAL_GAS * 1.0e3) / averageMolecularWeight;
  density = pressure / (specificGasConstant * temperature);

  // Gravity at local height z
  greal gravity = getGravity(latitude, radius, height);

  // Assumptions (checking occurs in debug mode only)
  assert(gravity != 0.0);

  // Pressure scale height at local height z
  pressureScaleHeight = UNIVERSAL_GAS * temperature / (averageMolecularWeight * gravity);

  // Height variation parameter evaluated at height z+1 km
  greal yPlus = ((height + 1.0) / (1.0 + (height + 1.0) / radius))-(baseHeight / (1.0 + baseHeight / radius));

  // Partial pressures and number densities at height z+1 km
  greal ppn2Plus = ppn2Base * exp(-yPlus / Hn2);
  greal ppch4Plus = ppch4Base * exp(-yPlus / Hch4);
  greal pparPlus = pparBase * exp(-yPlus / Har);
  greal n2ndPlus = ppn2Plus / (BOLTZMANN * temperature);
  greal ch4ndPlus = ppch4Plus / (BOLTZMANN * temperature);
  greal arndPlus = pparPlus / (BOLTZMANN * temperature);

  // Total number density an molecular weight at height z+1 km
  greal totndPlus = n2ndPlus + ch4ndPlus + arndPlus;
  assert(totndPlus != 0.0);
  greal amwPlus = (n2ndPlus * dinitrogen.averageMolecularWeight
    + ch4ndPlus * methane.averageMolecularWeight
    + arndPlus * argon.averageMolecularWeight) / totndPlus;

  // Density scale height, including effect of height variation of
  // molecular weight (effect of height variation of temperature is
  // zero, since temperature is constant exospheric value)
  densityScaleHeight = pressureScaleHeight / (1.0 - pressureScaleHeight * (amwPlus / averageMolecularWeight - 1.0));
}

} // namespace
