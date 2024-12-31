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
#include <fstream>
#include <iostream>
#include <iomanip>
#include "TitanGCMTerp.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
TitanGCMTerp::TitanGCMTerp()
  : TitanCommon(this)
{
  initializeData();
}

//! \fn  TitanGCMTerp::TitanGCMTerp(const TitanGCMTerp& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  TitanGCMTerp::~TitanGCMTerp()
//! \copydoc Atmosphere::~Atmosphere()

//! \brief Updates GCM data from model of Hourdin et al.
//!
//! Uses Titan GCM data from model of Hourdin et al., Icarus, vol.
//! 117, pp. 358-374 (1995) - Figures 7 and 8, for temperature and
//! zonal wind versus pressure at Ls = 0 and Ls = 270.  Estimates
//! other seasons by 180 degree Ls offset and hemispheric reversal.
//! Generates GCM profiles of temperature (Tgcm), density (dgcm),
//! and zonal wind (Ugcm), versus height (zgcm) and pressure (pNm2)
//! at 31 pressure levels of GCM data.
//!
//! \b Inputs
//! \arg #latitude      
//! \arg #solarTime      
//! \arg #longitudeSun      
//!
//! \retval #temperature03mb 
//! \retval #height03mb 
//! \retval #gcmHeight
//! \retval #gcmTemperature 
//! \retval #gcmDensity
//! \retval #gcmEwWind 
//! \retval #gcmPressureScaleHeight 
void TitanGCMTerp::updateGCMTables() {
  // The following tables have been populated:
  // gcmPressure[31], temperature000[31][19], temperature270[31][19],
  // ewWind000[31][19], ewWind270[31][19],
  // zetaPressure[6], Zeta[16][6]

  // Compute indexes for latitude interpolation of GCM data
  // Indices correspond to latitudes as follows:
  //   0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18
  // -90,-80,-70,-60,-50,-40,-30,-20,-10,  0, 10, 20, 30, 40, 50, 60, 70, 80, 90
  int lowLat = clampSize(int((latitude + 90.0_deg) / 10.0), latSize);
  int highLat = lowLat + 1;
  int oppLowLat = 18 - lowLat;
  int oppHighLat = oppLowLat - 1;

  // Interpolation factor for latitude interpolation
  greal lowerLatitude = -90.0_deg + 10.0_deg * lowLat;
  greal upperLatitude = lowerLatitude + 10.0_deg;
  Interpolator interp;
  interp.makeFraction(lowerLatitude, upperLatitude, latitude);

  // Pressure factor for first-guess altitudes
  greal totalDeltaLogPressure = log10(gcmPressure[0] / gcmPressure[presSize - 1]);

  // Sines and cosines for Ls evaluation
  greal sinLs = sin(toRadians(longitudeSun));
  greal cosLs = cos(toRadians(longitudeSun));
  greal cos2Ls = cos(toRadians(2.0 * longitudeSun));

  // Latitude factor to force diurnal variation to 0 at poles
  greal afLat = 1.0 - pow((latitude / 90.0_deg), 6);

  greal radius = getRadius(latitude);

  // Loop for hydrostatics calculations using first-guess altitudes
  greal deltaLogPressure[presSize] = {0.0};
  greal gravity[presSize] = {0};
  for (size_t i = 0; i < presSize; ++i) {
    // Assumptions (checking occurs in debug mode only)
    assert(gcmPressure[i] != 0.0);

    // Inter-level pressure changes for hydrostatic integration
    greal deltaLogPressure0 = log10(gcmPressure[0] / gcmPressure[i]);
    if (i > 0) {
      deltaLogPressure[i] = log(gcmPressure[i - 1] / gcmPressure[i]);
    }

    // First-guess altitudes from pressures
    gcmHeight[i] = 240.0_km * deltaLogPressure0 / totalDeltaLogPressure;

    // Daily average temperature at pressure level, interpolated to
    // latitude, at four seasonal Ls values
    greal Td000 = interp.linear(temperature000[i][lowLat], temperature000[i][highLat]);
    greal Td270 = interp.linear(temperature270[i][lowLat], temperature270[i][highLat]);
    greal Td180 = interp.linear(temperature000[i][oppLowLat], temperature000[i][oppHighLat]);
    greal Td090 = interp.linear(temperature270[i][oppLowLat], temperature270[i][oppHighLat]);

    // Coefficients for seasonal temperature variation versus Ls
    greal Acoef = 0.25 * (Td000 + Td090 + Td180 + Td270);
    greal Bcoef = 0.5  * (Td090 - Td270);
    greal Ccoef = 0.5  * (Td000 - Td180);
    greal Dcoef = 0.25 * (Td000 + Td180 - Td090 - Td270);

    // Daily average temperature at given latitude and Ls, with
    //  added diurnal variation, using assumed amplitude factor
    //  versus pressure level and latitude
    gcmTemperature[i] = Acoef + Bcoef * sinLs + Ccoef * cosLs + Dcoef * cos2Ls
            + afLat * (2.0 + 1.8 * deltaLogPressure0) * cos(toRadians(15.0_deg * (solarTime - 13.0_hr)));

    // Compressibility factor (zeta = p/rho*R*T)
    greal zetap = getNitrogenCompressibiltyFactor(gcmTemperature[i], gcmPressure[i]);

    // Daily average zonal wind at pressure level, interpolated to
    // latitude, at four seasonal Ls values
    greal Ud000 = interp.linear(ewWind000[i][lowLat], ewWind000[i][highLat]);
    greal Ud270 = interp.linear(ewWind270[i][lowLat], ewWind270[i][highLat]);
    greal Ud180 = interp.linear(ewWind000[i][oppLowLat], ewWind000[i][oppHighLat]);
    greal Ud090 = interp.linear(ewWind270[i][oppLowLat], ewWind270[i][oppHighLat]);

    // Coefficients for seasonal wind variation versus Ls
    Acoef = 0.25 * (Ud000 + Ud090 + Ud180 + Ud270);
    Bcoef = 0.5 * (Ud090 - Ud270);
    Ccoef = 0.5 * (Ud000 - Ud180);
    Dcoef = 0.25 * (Ud000 + Ud180 - Ud090 - Ud270);

    // Daily average zonal wind at given latitude and Ls , with no
    // diurnal variation added
    gcmEwWind[i] = Acoef + Bcoef * sinLs + Ccoef * cosLs + Dcoef*cos2Ls;

    // Gravity at latitude and first-guess altitude
    gravity[i] = getGravity(latitude, radius, gcmHeight[i]);

    // Assumptions (checking occurs in debug mode only)
    assert(zetap != 0.0);
    assert(gravity[i] != 0.0);
    assert(gcmTemperature[i] != 0.0);

    // Get pressure scale height and density
    gcmPressureScaleHeight[i] = zetap * UNIVERSAL_GAS * gcmTemperature[i] / (avgMolWt * gravity[i]);
    gcmDensity[i] = gcmPressure[i] * avgMolWt / (zetap * UNIVERSAL_GAS * 1.0e3 * gcmTemperature[i]);
  }

  // Get second-guess altitudes, refine gravity calculation, and
  // revise scale height values
  for (int i = 1; i < presSize; ++i) {
    gcmHeight[i] = gcmHeight[i - 1] + 0.5 * (gcmPressureScaleHeight[i - 1] + gcmPressureScaleHeight[i]) * deltaLogPressure[i];
    greal grav = getGravity(latitude, radius, gcmHeight[i]);
    // Assumptions (checking occurs in debug mode only)
    assert(grav != 0.0);
    gcmPressureScaleHeight[i] = gcmPressureScaleHeight[i] * gravity[i] / grav;
  }

  // Get final altitude values from revised scale heights
  for (int i = 1; i < presSize; ++i) {
    gcmHeight[i] = gcmHeight[i - 1] + 0.5 * (gcmPressureScaleHeight[i - 1] + gcmPressureScaleHeight[i]) * deltaLogPressure[i];
  }

  // Store height and temperature at top GCM level (0.3 mb)
  height03mb = gcmHeight[presSize - 1];
  temperature03mb = gcmTemperature[presSize - 1];
}

//! \brief Get Nitrogen (N2) compressibility factor.
//!
//! Compressibility factor zeta = p/(rho*R*T) for nitrogen (N2),
//! as function of temperature and pressure.
//! \param temp Temperature in Kelvin.
//! \param pres Pressure in N / m^2.
//! \returns Compressibility factor zeta = p/(rho*R*T) for nitrogen.
greal TitanGCMTerp::getNitrogenCompressibiltyFactor(greal temp, greal pres) {
  // Data base of Zeta factors have been read in and loaded into array zeta

  // Get p in atmospheres to interpolate over data
  greal patm = pres / STANDARD_ATMOSPHERE;

  // Get pressure index values k1,k2,k3 across which to interpolate
  size_t k1 = 3;
  size_t k2 = 4;
  size_t k3 = 5;
  for (size_t k = 1; k < zPresSize - 1; ++k) {
    if (patm < zetaPressure[k]) {
      k1 = k - 1;
      k2 = k;
      k3 = k + 1;
      break;
    }
  }

  // Interpolation factor is based on current temperature.
  greal dT = fmod(temp, 10.0) / 10.0;
  Interpolator interp(dT);

  // Get temperature index values j1 and j2 across which to
  // interpolate
  int j1 = clampSize(int((temp - 90.0) / 10.0) - 1, zTempSize);
  int j2 = clampSize(j1 + 1, zTempSize);

  // Linear interpolation on T to get zeta=Zeta1 at p(k1)
  greal zeta1 = interp.linear(zeta[j1][k1], zeta[j2][k1]);

  // Linear interpolation on T to get zeta=Zeta2 at p(k2)
  greal zeta2 = interp.linear(zeta[j1][k2], zeta[j2][k2]);

  // Linear interpolation on T to get zeta=Zeta3 at p(k3)
  greal zeta3 = interp.linear(zeta[j1][k3], zeta[j2][k3]);

  // Quadratic interpolation on pressure to get zeta at input pressure
  greal zetaN2;
  LagrangeQuad(zetaPressure[k1], zetaPressure[k2], zetaPressure[k3],
                             zeta1, zeta2, zeta3, patm, zetaN2);
  return zetaN2;
}

//! \copydoc Atmosphere::update()
void TitanGCMTerp::update()
{
  // Uses profile of temperature, density, pressure, wind, versus
  // height, at latitude, local time, and Ls value, generated by
  // subroutine GCMptduz_T04 from Titan GCM data from model of
  // Hourdin et al., Icarus, vol.117, pp. 358-374 (1995) - Figures 7
  // and 8, and interpolates to desired height z.

  // Loop to find height index values i1 and i2=i1+1 for height interpolation
  size_t i;
  for ( i = 0; i < presSize - 2; ++i) {
    if (height <= gcmHeight[i + 1])
      break;
  }
  size_t lowerHt = i;
  size_t upperHt = i + 1;

  // Height interpolation factor
  Interpolator interp;
  interp.makeFraction(gcmHeight[lowerHt], gcmHeight[upperHt], height);
  greal deltaHeight = (gcmHeight[upperHt] - gcmHeight[lowerHt]);

  // Linear interpolation on temperature
  temperature = interp.linear(gcmTemperature[lowerHt], gcmTemperature[upperHt]);

  // Logarithmic interpolation on pressure
  pressure = interp.log(gcmPressure[lowerHt], gcmPressure[upperHt]);

  // Assumptions (checking occurs in debug mode only)
  assert(gcmDensity[lowerHt] != 0.0);
  assert(gcmDensity[upperHt] != 0.0);
  assert(gcmTemperature[lowerHt] != 0.0);
  assert(gcmTemperature[upperHt] != 0.0);
  assert(gcmPressure[lowerHt] != 0.0);
  assert(gcmPressure[upperHt] != 0.0);

  // Linear interpolation between gas constant at two heights
  greal R1 = gcmPressure[lowerHt] / (gcmDensity[lowerHt] * gcmTemperature[lowerHt]);
  greal R2 = gcmPressure[upperHt] / (gcmDensity[upperHt] * gcmTemperature[upperHt]);
  specificGasConstant = interp.linear(R1, R2);

  // Assumptions (checking occurs in debug mode only)
  assert(specificGasConstant != 0.0);
  assert(temperature != 0.0);

  // Density from pressure, temperature and interpolated gas constant
  density = pressure / (specificGasConstant * temperature);

  // Pressure scale height between levels i1 and i2
  pressureScaleHeight = deltaHeight / log(gcmPressure[lowerHt] / gcmPressure[upperHt]);

  // Density scale height between levels i1 and i2
  densityScaleHeight = deltaHeight / log(gcmDensity[lowerHt] / gcmDensity[upperHt]);

  // Molecular weight for GCM height range
  averageMolecularWeight = avgMolWt;

  // Nominal mole fractions in GCM height range
  greal mfNitrogen = mfNitrogenDefault;
  greal mfMethane = mfMethaneDefault;
  greal mfArgon = mfArgonDefault;

  // Adjust mole fractions if user supplies methane mole fraction from NAMELIST input file
  if (methane.moleFraction > 0.0) {
    mfMethane = methane.moleFraction;
    mfNitrogen = 1.022857 - 2.428571 * mfMethane;
    mfArgon = -0.022857 + 1.428571 * mfMethane;
    // Limit Argon mole fraction to >= 0 (limits CH4 mole fraction for molecular weight to remain
    // unchanged)
    if (mfArgon <= 0.0) {
      mfArgon = 0.0;
      mfMethane = 0.016;
      mfNitrogen = 0.984;
    }
  }

  // Use nominal mole fractions in GCM height range
  // Compute species number densities and total number density
  dinitrogen.numberDensity = 1.0e3 * mfNitrogen * density * AVOGADRO / averageMolecularWeight;
  methane.numberDensity = 1.0e3 * mfMethane * density * AVOGADRO / averageMolecularWeight;
  argon.numberDensity = 1.0e3 * mfArgon * density * AVOGADRO / averageMolecularWeight;
  totalNumberDensity = dinitrogen.numberDensity + methane.numberDensity + argon.numberDensity;

  // Linear interpolation on zonal wind (daily value - no diurnal model applied)
  ewWind = interp.linear(gcmEwWind[lowerHt], gcmEwWind[upperHt]);

  // Set meridional wind to zero
  nsWind = 0.0;
}

//! \brief Quadratic interpolation routine.
//!
//! Quadratic Lagrange interpolation over x1,y1,  x2,y2,  and 
//! x3,y3  to get y at x.  Assumes x values are in order      
//! x1 <= x2 <= x3
//! \param x1  Smallest x value.
//! \param x2  Middle x value.
//! \param x3  Largest x value.
//! \param y1  y value paired with x1.
//! \param y2  y value paired with x2.
//! \param y3  y value paired with x3.
//! \param x   Current x value.
//! \param[out] y   A greal.
//! \retval y  Current y value.
void TitanGCMTerp::LagrangeQuad(greal x1, greal x2, greal x3,
  greal y1, greal y2, greal y3, greal x, greal& y)
{
  // Make sure x is >= x1
  if (x < x1) x = x1;
  // Make sure x <= x3
  if (x > x3) x = x3;
  // Use linear Lagrange interpolation if x1=x2 or x1=x3 or x2=x3
  if (x1 == x2) {
    LagrangeLin(x2, x3, y2, y3, x, y);
  }
  else if (x1 == x3) {
    LagrangeLin(x1, x2, y1, y2, x, y);
  }
  else if (x2 == x3) {
    LagrangeLin(x1, x2, y1, y2, x, y);
  }
  else {
    // Canonical Lagrange interpolation formula
    y = y1 * (x - x2)*(x - x3) / ((x1 - x2)*(x1 - x3))
      + y2 * (x - x1)*(x - x3) / ((x2 - x1)*(x2 - x3))
      + y3 * (x - x1)*(x - x2) / ((x3 - x1)*(x3 - x2));
  }
}

//! \brief Linear interpolation routine.
//!
//! Linear Lagrange interpolation between x1,y1 and x2,y2 to get y at x.
//! \param x1  Smallest x value.
//! \param x2  Largest x value.
//! \param y1  y value paired with x1.
//! \param y2  y value paired with x2.
//! \param x   Current x value.
//! \param[out] y   A greal.
//! \retval y  Current y value.
void TitanGCMTerp::LagrangeLin(greal x1, greal x2, greal y1, greal y2, greal x, greal& y)
{
  // Assume y = y1 (= y2) if y1 = y2
  if (x1 == x2) {
    y = y1;
  }
  else {
    // Canonical Lagrange interpolation formula
    y = y1 * (x - x2) / (x1 - x2) + y2 * (x - x1) / (x2 - x1);
  }
}

} // namespace
