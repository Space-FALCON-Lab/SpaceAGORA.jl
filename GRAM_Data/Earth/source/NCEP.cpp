//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United
// States without explicit approval by NASA Marshall Space Flight Center.
//
// Module: Earth-GRAM
// Point of Contact: P. White
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <algorithm>
#include <cmath>

#include "NCEP.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {


//! \copydoc Atmosphere::Atmosphere()
NCEP::NCEP()
  : Atmosphere(), EarthCommon(this)
{
  atmos.setPlanetSpecificMetrics(earthAtmos);

  temp = new greal[arraySize3D];
  dens = new greal[arraySize3D];
  dewp = new greal[arraySize3D];
  uwnd = new greal[arraySize3D];
  vwnd = new greal[arraySize3D];
  geop = new greal[arraySize3D];
  stmp = new greal[arraySize3D];
  sden = new greal[arraySize3D];
  sdewp = new greal[arraySize3D];
  suwd = new double[arraySize3D];
  svwd = new greal[arraySize3D];
  sprs = new greal[arraySize3D];
  rhum = new greal[arraySize3D];
  srhum = new greal[arraySize3D];
  vprs = new greal[arraySize3D];
  svprs = new greal[arraySize3D];
  spdav = new greal[arraySize3D];
  spdsd = new greal[arraySize3D];
  uvcor = new greal[arraySize3D];

  slpn = new greal[arraySize2D];
  sfcp = new greal[arraySize2D];
  sslp = new greal[arraySize2D];
  ssfcp = new greal[arraySize2D];
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
NCEP::NCEP(const NCEP& orig)
  : Atmosphere(orig), EarthCommon(this)
{
  initialized = false;
  earthAtmos = orig.earthAtmos;
  NCEPYear = orig.NCEPYear;
  NCEPHour = orig.NCEPHour;
  hour = orig.hour;
  month = orig.month;
  densityPerturbationScale = orig.densityPerturbationScale;
  pressurePerturbationScale = orig.pressurePerturbationScale;
  temperaturePerturbationScale = orig.temperaturePerturbationScale;
  ewWindPerturbationScale = orig.ewWindPerturbationScale;
  nsWindPerturbationScale = orig.nsWindPerturbationScale;
  NCEPPath = orig.NCEPPath;

  atmos.setPlanetSpecificMetrics(earthAtmos);

  temp = new greal[arraySize3D];
  dens = new greal[arraySize3D];
  dewp = new greal[arraySize3D];
  uwnd = new greal[arraySize3D];
  vwnd = new greal[arraySize3D];
  geop = new greal[arraySize3D];
  stmp = new greal[arraySize3D];
  sden = new greal[arraySize3D];
  sdewp = new greal[arraySize3D];
  suwd = new double[arraySize3D];
  svwd = new greal[arraySize3D];
  sprs = new greal[arraySize3D];
  rhum = new greal[arraySize3D];
  srhum = new greal[arraySize3D];
  vprs = new greal[arraySize3D];
  svprs = new greal[arraySize3D];
  spdav = new greal[arraySize3D];
  spdsd = new greal[arraySize3D];
  uvcor = new greal[arraySize3D];

  slpn = new greal[arraySize2D];
  sfcp = new greal[arraySize2D];
  sslp = new greal[arraySize2D];
  ssfcp = new greal[arraySize2D];
}

//! \copydoc Atmosphere::~Atmosphere()
NCEP::~NCEP()
{
  delete[] temp;
  delete[] dens;
  delete[] dewp;
  delete[] uwnd;
  delete[] vwnd;
  delete[] geop;
  delete[] stmp;
  delete[] sden;
  delete[] sdewp;
  delete[] suwd;
  delete[] svwd;
  delete[] sprs;
  delete[] rhum;
  delete[] srhum;
  delete[] vprs;
  delete[] svprs;
  delete[] spdav;
  delete[] spdsd;
  delete[] uvcor;

  delete[] slpn;
  delete[] sfcp;
  delete[] sslp;
  delete[] ssfcp;
}

//! \copydoc PerturbedAtmosphere::setInputParameters()
void NCEP::setInputParameters(const EarthInputParameters& params)
{
  NCEPYear = params.NCEPYear;
  NCEPHour = params.NCEPHour;
  hour = params.hour;
  month = params.month;
  densityPerturbationScale = params.randomPerturbationScale;
  pressurePerturbationScale = params.randomPerturbationScale;
  temperaturePerturbationScale = params.randomPerturbationScale;
  ewWindPerturbationScale = params.horizontalWindPerturbationScale;
  nsWindPerturbationScale = params.horizontalWindPerturbationScale;
}


void NCEP::updateSurface() {
  surfaceUpdate = true;
  greal saveHeight = height;
  height = position.surfaceHeight;
  update();
  height = saveHeight;
  surfaceUpdate = false;
}

//! \fn NCEP::setYearAndHour(int NCEPyr, int NCEPhr)
//! \brief Identify the desired NCEP data.
//!
//! The NCEP data is identified by specifying two codes, a year and hour code.  NCEP monthly 
//! climatology is determined by the value of month in the initial time input.
//! 
//! \param NCEPyr   The NCEPYear code y1y2 specifies NCEP climatology for period-of-record (POR) from year
//!                   y1 through year y2 (e.g. NCEPYear = 9008 for POR = 1990 through 2008). 
//! \param NCEPhr   The NCEPHour code specifies the UT hour of day for NCEP climatology:
//!                   1 = 00 UT, 2 = 06 UT, 3 = 12 UT, 4 = 18 UT, 5 = all times of day combined,
//!                   or 0 to use NCEP time-of-day based on input UTC hour.

//! \fn NCEP::getTemperatureGradient()
//! \brief Gets the temperature gradient with respect to height.
//!
//! \returns Temperature gradient with respect to height \units{K}.

//! \brief Compute mean pressure, density, temperature, wind components, 
//! dewpoint temperature and their standard deviations, from the NCEP data
//!
//! This routine will read the NCEP data, if necessary, and perform interpolations based on the
//! current position.
//!
//! \b Inputs
//! \arg #position          
//! \arg NCEP data          
//!
//! \retval #pressure
//! \retval #temperature
//! \retval #density
//! \retval #ewWind
//! \retval #nsWind
//! \retval #verticalWind
//! \retval #dewPoint
//! \retval #pressureStandardDeviation
//! \retval #densityStandardDeviation
//! \retval #temperatureStandardDeviation
//! \retval #ewStandardDeviation
//! \retval #nsStandardDeviation
//! \retval #dewPointSD
//! \retval #relativeHumidity
//! \retval #relativeHumiditySD
//! \retval #vaporPressure
//! \retval #vaporPressureSD
//! \retval #windSpeed
//! \retval #windSpeedStandardDeviation
//! \retval #windCorrelation
void NCEP::update()
{
  constexpr greal grid = 0.4;
  constexpr greal onek = 1000.0;

  readNCEPFile();

  // Terminate if phi or thet too large
  if ((abs(latitude) > 90.0) & (abs(longitude) > 360.0)) {
    return;
  }

  // Find i, j indexes for next lower lat-lon grid point
  size_t latIndex = int((latitude + 90.0) * grid);
  if (latIndex == latSize - 1) {
    latIndex = latSize - 2;
  }

  greal elon = longitude;
  if (elon < 0.0) {
    elon = 360.0_deg;
  }

  size_t lonIndex = int(elon * grid);
  if (lonIndex > lonSize - 1) {
    lonIndex = lonIndex - (lonSize - 1);
  }

  greal p[2][2][presSize], t[2][2][presSize], rho[2][2][presSize], u[2][2][presSize], v[2][2][presSize], w[2][2][presSize];
  greal sp[2][2][presSize], st[2][2][presSize], srho[2][2][presSize], su[2][2][presSize], sv[2][2][presSize];
  greal td[2][2][presSize], vp[2][2][presSize], rh[2][2][presSize], spdavl[2][2][presSize];
  greal std[2][2][presSize], svp[2][2][presSize], srh[2][2][presSize], spdsdl[2][2][presSize];
  greal uvcorrl[2][2][presSize], h[2][2][presSize];
  greal lonFactor;
  greal latFactor;

  lonFactor = elon * grid - int(elon * grid);
  latFactor = (latitude + 90.0 - (latIndex - 1.0) / grid) * grid - 1.0;

  // Get profiles at "corners" of lat-lon "square" for interpolation
  for (size_t iLat = 0; iLat < 2; iLat++) {
    size_t jLat = latIndex + iLat;

    for (size_t iLon = 0; iLon < 2; iLon++) {
      size_t jLon = lonIndex + iLon;

      // Get pressure, sigma-pressure, and vertical wind at sea level and surface
      size_t index2D = getIndex(jLat, jLon);
      p[iLat][iLon][0] = slpn[index2D];
      p[iLat][iLon][1] = sfcp[index2D];
      sp[iLat][iLon][0] = sslp[index2D];
      sp[iLat][iLon][1] = ssfcp[index2D];
      w[iLat][iLon][0] = 0.0;
      w[iLat][iLon][1] = 0.0;

      // Get "corner" profiles from full NCEP arrays
      for (size_t ijPres = 0; ijPres < presSize; ijPres++) {
        size_t index3D = getIndex(ijPres, jLat, jLon);
        if (ijPres > 1) {
          p[iLat][iLon][ijPres] = pn[ijPres];
          sp[iLat][iLon][ijPres] = sprs[index3D] * pressurePerturbationScale;
        }
        t[iLat][iLon][ijPres] = temp[index3D];
        td[iLat][iLon][ijPres] = dewp[index3D];
        rho[iLat][iLon][ijPres] = dens[index3D];
        u[iLat][iLon][ijPres] = uwnd[index3D];
        v[iLat][iLon][ijPres] = vwnd[index3D];
        h[iLat][iLon][ijPres] = geop[index3D];
        st[iLat][iLon][ijPres] = stmp[index3D] * temperaturePerturbationScale;
        std[iLat][iLon][ijPres] = sdewp[index3D];
        srho[iLat][iLon][ijPres] = sden[index3D] * densityPerturbationScale;
        su[iLat][iLon][ijPres] = suwd[index3D] * ewWindPerturbationScale;
        sv[iLat][iLon][ijPres] = svwd[index3D] * nsWindPerturbationScale;
        rh[iLat][iLon][ijPres] = rhum[index3D];
        srh[iLat][iLon][ijPres] = srhum[index3D];
        vp[iLat][iLon][ijPres] = vprs[index3D];
        svp[iLat][iLon][ijPres] = svprs[index3D];
        spdavl[iLat][iLon][ijPres] = spdav[index3D];
        spdsdl[iLat][iLon][ijPres] = spdsd[index3D];
        uvcorrl[iLat][iLon][ijPres] = uvcor[index3D];
      }
    }
  }

  // Compute surface gravity and radius
  greal g0, re;
  getEffectiveRadius(latitude, g0, re);

  if (!surfaceUpdate) {

    // Compute vertical winds from Montgomery stream function
    // Get distances dx, dy across lat-lon "square"
    greal dy = onek * TO_RADIANS * totalRadius / grid;
    greal dx = 0.01 * dy;
    if (abs(latitude) < 89.0) {
      dx = dy * cos(toRadians(latitude));
    }

    // Get horizontal temperature gradients
    for (size_t iPres = 2; iPres < presSize; iPres++) {
      dtdx[iPres] = (temp[getIndex(iPres, latIndex, lonIndex + 1)]
        - temp[getIndex(iPres, latIndex, lonIndex)])
        / dx;
      dtdy[iPres] = (temp[getIndex(iPres, latIndex + 1, lonIndex)]
        - temp[getIndex(iPres, latIndex, lonIndex)])
        / dy;
    }

    // Get profiles at "corners" of lat-lon "square" for interpolation
    for (size_t iLat = 0; iLat < 2; iLat++) {
      for (size_t iLon = 0; iLon < 2; iLon++) {
        // Get vertical wind at "corners" of lat-lon "square"
        for (size_t iPres = 2; iPres < presSize; iPres++) {
          // This clampSize ensures k1 a max of presSize-2.
          size_t k1 = clampSize(iPres, presSize - 1);

          greal dtdz = (t[iLat][iLon][k1 + 1] - t[iLat][iLon][k1])
                       / (onek * (h[iLat][iLon][k1 + 1] - h[iLat][iLon][k1]));
          greal cp = 7.0 * getGasConstant(td[iLat][iLon][iPres], pn[iPres]) / 2.0;
          greal dtdz1 = max(1.0e-4, dtdz + g0 / cp);

          w[iLat][iLon][iPres] = -(u[iLat][iLon][iPres] * dtdx[iPres] + v[iLat][iLon][iPres] * dtdy[iPres]) / dtdz1;
        }
      }
    }
  }

  // Convert input geometric height to geopotential height, for interpolation
  greal hz = getGeopotentialHeight(g0, 1000.0 * re, 1000.0 * height) / 1000.0;

  greal pz[2][2], tz[2][2], rhoz[2][2], uz[2][2], vz[2][2], wz[2][2];
  greal spz[2][2], stz[2][2], srhoz[2][2], suz[2][2], svz[2][2];
  greal tdz[2][2], vpn[2][2], rhn[2][2], spdavz[2][2];
  greal stdz[2][2], svpn[2][2], srhn[2][2], spdsdz[2][2];
  greal uvcorrz[2][2], dtdz[2][2];

  // Interpolate over height using the pressure level dimension.
  // Step through 4 "corners" of lat-lon "square"
  for (size_t iLat = 0; iLat < 2; iLat++) {

    for (size_t iLon = 0; iLon < 2; iLon++) {

      greal zsort[2][2][presSize];
      zsort[iLat][iLon][0] = 0.0;
      zsort[iLat][iLon][1] = h[iLat][iLon][1];
      for (size_t iPres = 2; iPres < presSize; iPres++) {
        zsort[iLat][iLon][iPres] = h[iLat][iLon][iPres];
      }

      // Sort levels to put surface and sea level in right order
      size_t lsort[2][2][presSize];
      sortLevel(lsort, zsort);

      // Find index values lowerHeight and upperHeight for interpolation
      size_t lowerHeight, upperHeight;
      if (hz < 0.0) {
        lowerHeight = lsort[iLat][iLon][0];
        upperHeight = lsort[iLat][iLon][1];
      }
      else {
        for (size_t iPres = 0; iPres < presSize - 1; iPres++) {
          lowerHeight = lsort[iLat][iLon][iPres];
          upperHeight = lsort[iLat][iLon][iPres + 1];
          if ((hz < h[iLat][iLon][lowerHeight]) || (hz >= h[iLat][iLon][upperHeight]))
            continue;
          break;
        }
      }

      // Density scale height, km
      greal dh;
      if (abs(h[iLat][iLon][lowerHeight] - h[iLat][iLon][upperHeight]) == 0.0) {
        dtdz[iLat][iLon] = 0.0;
        dh = 0.0;
      }
      else {
        dtdz[iLat][iLon] = (t[iLat][iLon][upperHeight] - t[iLat][iLon][lowerHeight])
          / (h[iLat][iLon][upperHeight] - h[iLat][iLon][lowerHeight]);
        dh = (hz - h[iLat][iLon][lowerHeight])
          / (h[iLat][iLon][upperHeight] - h[iLat][iLon][lowerHeight]);
      }

      greal hden;
      if (abs(rho[iLat][iLon][lowerHeight] - rho[iLat][iLon][upperHeight]) == 0.0
        || abs(h[iLat][iLon][lowerHeight] - h[iLat][iLon][upperHeight]) == 0.0) {
        hden = 1.0;
        rhoz[iLat][iLon] = rho[iLat][iLon][lowerHeight];
      }
      else {
        hden = (h[iLat][iLon][upperHeight] - h[iLat][iLon][lowerHeight])
          / log(rho[iLat][iLon][lowerHeight] / rho[iLat][iLon][upperHeight]);
        // Logarithmic interpolation for density
        rhoz[iLat][iLon] = rho[iLat][iLon][lowerHeight] * exp((h[iLat][iLon][lowerHeight] - hz) / hden);
      }

      Interpolator interpHeight(dh);

      // Linear interpolation on temperature, dewpoint, and gas constant
      tz[iLat][iLon] = interpHeight.linear(t[iLat][iLon][lowerHeight], t[iLat][iLon][upperHeight]);
      tdz[iLat][iLon] = interpHeight.linear(td[iLat][iLon][lowerHeight], td[iLat][iLon][upperHeight]);

      greal r1 = p[iLat][iLon][lowerHeight]
        / (rho[iLat][iLon][lowerHeight] * t[iLat][iLon][lowerHeight]);
      greal r2 = p[iLat][iLon][upperHeight]
        / (rho[iLat][iLon][upperHeight] * t[iLat][iLon][upperHeight]);
      greal rz = interpHeight.linear(r1, r2);

      // Pressure (mb) from gas law
      pz[iLat][iLon] = rhoz[iLat][iLon] * rz * tz[iLat][iLon];

      // Interpolate winds and standard deviations linearly
      uz[iLat][iLon] = interpHeight.linear(u[iLat][iLon][lowerHeight], u[iLat][iLon][upperHeight]);
      vz[iLat][iLon] = interpHeight.linear(v[iLat][iLon][lowerHeight], v[iLat][iLon][upperHeight]);
      spz[iLat][iLon] = interpHeight.linear(sp[iLat][iLon][lowerHeight], sp[iLat][iLon][upperHeight]);
      stz[iLat][iLon] = interpHeight.linear(st[iLat][iLon][lowerHeight], st[iLat][iLon][upperHeight]);
      suz[iLat][iLon] = interpHeight.linear(su[iLat][iLon][lowerHeight], su[iLat][iLon][upperHeight]);
      svz[iLat][iLon] = interpHeight.linear(sv[iLat][iLon][lowerHeight], sv[iLat][iLon][upperHeight]);
      uvcorrz[iLat][iLon] = interpHeight.linear(uvcorrl[iLat][iLon][lowerHeight], uvcorrl[iLat][iLon][upperHeight]);
      spdavz[iLat][iLon] = interpHeight.linear(spdavl[iLat][iLon][lowerHeight], spdavl[iLat][iLon][upperHeight]);
      spdsdz[iLat][iLon] = interpHeight.linear(spdsdl[iLat][iLon][lowerHeight], spdsdl[iLat][iLon][upperHeight]);

      if (!surfaceUpdate) {
        wz[iLat][iLon] = interpHeight.linear(w[iLat][iLon][lowerHeight], w[iLat][iLon][upperHeight]);
        stdz[iLat][iLon] = interpHeight.linear(std[iLat][iLon][lowerHeight], std[iLat][iLon][upperHeight]);
        srhoz[iLat][iLon] = interpHeight.linear(srho[iLat][iLon][lowerHeight], srho[iLat][iLon][upperHeight]);
        rhn[iLat][iLon] = interpHeight.linear(rh[iLat][iLon][lowerHeight], rh[iLat][iLon][upperHeight]);
        srhn[iLat][iLon] = interpHeight.linear(srh[iLat][iLon][lowerHeight], srh[iLat][iLon][upperHeight]);
        vpn[iLat][iLon] = interpHeight.log(vp[iLat][iLon][lowerHeight], vp[iLat][iLon][upperHeight]);
        svpn[iLat][iLon] = interpHeight.linear(svp[iLat][iLon][lowerHeight], svp[iLat][iLon][upperHeight]);
      }
    }
  }

  // Lat-lon interpolation to get output values at input height
  Interpolator interpLatLon(latFactor, lonFactor);

  pressure = interpLatLon.linear(pz);
  temperature = interpLatLon.linear(tz);
  density = interpLatLon.linear(rhoz);
  ewWind = interpLatLon.linear(uz);
  nsWind = interpLatLon.linear(vz);
  pressureStandardDeviation = square(interpLatLon.linear(spz) / pressure);
  temperatureStandardDeviation = square(interpLatLon.linear(stz) / temperature);
  densityStandardDeviation = square(interpLatLon.linear(srhoz) / density);
  ewStandardDeviation = interpLatLon.linear(suz);
  nsStandardDeviation = interpLatLon.linear(svz);
  windCorrelation = interpLatLon.linear(uvcorrz);
  windSpeed = interpLatLon.linear(spdavz);
  windSpeedStandardDeviation = interpLatLon.linear(spdsdz);

  if (!surfaceUpdate) {
    verticalWind = interpLatLon.linear(wz);
    dewPoint = interpLatLon.linear(tdz);
    dewPointSD = interpLatLon.linear(stdz);
    temperatureGradient = interpLatLon.linear(dtdz);
    relativeHumidity = interpLatLon.linear(rhn);
    relativeHumiditySD = interpLatLon.linear(srhn);
    vaporPressure = interpLatLon.linear(vpn);
    vaporPressureSD = interpLatLon.linear(svpn);
  }
}

//! \brief Compute gas constant from dewpoint and pressure
//!
//! The function computes the gas constant (R) based on dewpont temperature and pressure.
//!
//! \param tdc Dewpoint temperature \units{\text{degrees K}}.
//! \param pres Pressure \units{Pa}.
//!
//! \retval The gas constant (R).
greal NCEP::getGasConstant(greal tdc, greal pres)
{
  constexpr greal omeps = 0.37803;  // one minus epsilon
//  constexpr greal ckf = 273.15;     // C to K
  constexpr greal onec = 100.0;

//  greal e = wexler(tdc + ckf, onec * pres) / onec;
  greal e = wexler(tdc, pres) / onec;

  greal gascon_out = dryAirGasConstant / (1.0 - omeps * e / pres);

  return gascon_out;
}

//! \brief Converts input geometric height into output geopotential height. 
//! 
//! Note: r, z, and H must all be in the same units (e.g. all in meters
//! or all in km)
//!
//! \param g0  local surface gravity (m / s^2),
//! \param r   local effective Earth radius
//! \param z   gemoetric height
//!
//! \returns Geopotential height
greal NCEP::getGeopotentialHeight(greal g0, greal r, greal z)
{
  greal geopotentialHeight = r * z * g0 / (referenceGravity * (r + z));

  return geopotentialHeight;
}

//! \brief Calculate heights of 10 mb and 20mb pressure levels.
//!
//! Compute geometric height of lowest 10 mb level and highest 20 mb
//! level from four 2.5 by 2.5 degree NCEP profiles surrounding given lat-lon.
//!
//! \param lat  Latitude.
//! \param lon  Longitude.
//! \param[out] highest20mb  A greal.
//! \param[out] lowest10mb   A greal.
//!
//! \retval highest20mb  Height \units{km} of the highest 20 mb pressure level.
//! \retval lowest10mb   Height \units{km} of the lowest 10 mb pressure level.
void NCEP::getHeights(greal lat, greal lon, greal &highest20mb, greal &lowest10mb)
{
  readNCEPFile();

  // Do we need this?
  greal elon = lon;
  if (elon < 0.0) elon = elon + 360.0;

  constexpr greal grid = 0.4;

  // Find lat and lon indexes for next lower lat-lon grid point
  // Max out at the penultimate index.
  size_t latIndex = clampSize(size_t((lat + 90.0) * grid), latSize - 1);
  size_t lonIndex = clampSize(size_t(elon * grid), lonSize - 1);

  // Convert latitude to radians and get local gravity and effective radius
  greal gg0, effectiveRadius=0;
  getEffectiveRadius(lat, gg0, effectiveRadius);
  greal grav = getGravity(lat, getRadius(lat), 0.0);

  highest20mb = 0.0;
  lowest10mb = 99.9;

  //Get geometric heights at corners of NCEP lat-lon "sqaure"
  for (size_t m = 0; m < 2; m++) {
    size_t iLat = latIndex + m;
    for (size_t n = 0; n < 2; n++) {
      size_t iLon = lonIndex + n;

      // Get geopotential heights of 10mb and 20mb levels
      greal h10 = greal(geop[getIndex(presSize - 1, iLat, iLon)]);
      greal h20 = greal(geop[getIndex(presSize - 2, iLat, iLon)]);

      // Get geometric heights at 10mb and 20mb levels
      greal z10 = (h10 * referenceGravity * effectiveRadius
                   / (effectiveRadius * grav - h10 * referenceGravity));
      greal z20 = (h20 * referenceGravity * effectiveRadius
                   / (effectiveRadius * grav - h20 * referenceGravity));

      // Find lowest10mb = min z10 and highest20mb = max z20
      if (z10 < lowest10mb) {
        lowest10mb = z10;
      }
      if (z20 > highest20mb) {
        highest20mb = z20;
      }
    }
  }
}

//! \brief Sorts level numbers, to put surface and sea level in right order.
//!
//! \param[out] lsort   An array of type size_t.
//! \param zsort        A data array.
//! 
//! \retval lsort  An indirect reference array. That is, an array of indices.
void NCEP::sortLevel(size_t lsort[2][2][presSize], greal zsort[2][2][presSize])
{
  // Step through "corners" of lat-lon "square"
  for (size_t iLat = 0; iLat < 2; iLat++) {
    for (size_t iLon = 0; iLon < 2; iLon++) {
      // Initialize level numbers
      for (size_t i = 0; i < presSize; i++) {
        lsort[iLat][iLon][i] = i;
      }

      // Store surface (level 1) and sea level (level 0) altitudes
      greal zsort1 = zsort[iLat][iLon][1];
      greal zsort0 = zsort[iLat][iLon][0];

      // Sort to put surface level number into right order
      for (size_t iPres = 2; iPres < presSize; iPres++) {
        if (zsort1 > zsort[iLat][iLon][iPres]) {
          lsort[iLat][iLon][iPres - 1] = lsort[iLat][iLon][iPres];
          lsort[iLat][iLon][iPres] = 1;
        }
        else {
          break;
        }
      }

      // Store altitudes with surface altitude in right order
      greal qqsort[presSize];
      for (size_t iPres = 0; iPres < presSize; iPres++) {
        size_t lvl = lsort[iLat][iLon][iPres];
        qqsort[iPres] = zsort[iLat][iLon][lvl];
      }

      // Sort to put sea level number into right order
      for (size_t iPres = 1; iPres < presSize; iPres++) {
        if (zsort0 > qqsort[iPres]) {
          lsort[iLat][iLon][iPres - 1] = lsort[iLat][iLon][iPres];
          lsort[iLat][iLon][iPres] = 0;
        }
        else {
          break;
        }
      }
    }
  }
}

//! \fn NCEP::getIndex(size_t ipres, size_t ilat, size_t ilon)
//! \brief Generates a 1D index from 3 indices.
//! \param ipres A pressure index.
//! \param ilat  A latitude index.
//! \param ilon  A longitude index.
//! \returns A 1D index.

//! \fn NCEP::getIndex(size_t ilat, size_t ilon)
//! \brief Generates a 1D index from 2 indices.
//! \param ilat  A latitude index.
//! \param ilon  A longitude index.
//! \returns A 1D index.


} // namespace
