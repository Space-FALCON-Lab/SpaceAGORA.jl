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

#include "MERRA2.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {


//! \copydoc Atmosphere::Atmosphere()
MERRA2::MERRA2() : Atmosphere(), EarthCommon(this)
{
  atmos.setPlanetSpecificMetrics(earthAtmos);
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
MERRA2::MERRA2(const MERRA2& orig) : Atmosphere(orig), EarthCommon(this)
{
  initialized = false;
  earthAtmos = orig.earthAtmos;
  atmos.setPlanetSpecificMetrics(earthAtmos);
  hour = orig.hour;
  month = orig.month;
  densityPerturbationScale = orig.densityPerturbationScale;
  pressurePerturbationScale = orig.pressurePerturbationScale;
  temperaturePerturbationScale = orig.temperaturePerturbationScale;
  ewWindPerturbationScale = orig.ewWindPerturbationScale;
  nsWindPerturbationScale = orig.nsWindPerturbationScale;
  M2Path = orig.M2Path;
  M2Hour = orig.M2Hour;
  minimumLatitude = orig.minimumLatitude;
  maximumLatitude = orig.maximumLatitude;
  minimumLongitude = orig.minimumLongitude;
  maximumLongitude = orig.maximumLongitude;
  blockOfDoubles = nullptr;
}

//! \copydoc Atmosphere::~Atmosphere()
MERRA2::~MERRA2()
{
  deleteArrays();
}

//! \copydoc PerturbedAtmosphere::setInputParameters()
void MERRA2::setInputParameters(const EarthInputParameters& params)
{

  M2Hour = params.M2Hour;
  minimumLatitude = params.minimumLatitude;
  maximumLatitude = params.maximumLatitude;
  minimumLongitude = params.minimumLongitude;
  maximumLongitude = params.maximumLongitude;
  hour = params.hour;
  month = params.month;
  densityPerturbationScale = params.randomPerturbationScale;
  pressurePerturbationScale = params.randomPerturbationScale;
  temperaturePerturbationScale = params.randomPerturbationScale;
  ewWindPerturbationScale = params.horizontalWindPerturbationScale;
  nsWindPerturbationScale = params.horizontalWindPerturbationScale;
}

//! \fn MERRA2::setHour(int M2hr)
//! \brief Identify the desired MERRA-2 data.
//!
//! The MERRA-2 data is identified by specifying an hour code.  MERRA-2 monthly 
//! climatology is determined by the value of month in the initial time input.
//! 
//! \param M2hr   The M2hr code specifies the UT hour of day for MERRA-2 climatology:
//!                   1 = 00 UT, 2 = 03 UT, 3 = 06 UT, 4 = 09 UT,
//!                   5 = 12 UT, 6 = 15 UT, 7 = 18 UT, 8 = 21 UT,
//!                   9 = all times of day combined,
//!                   or 0 to use NCEP time-of-day based on input UTC hour.

//! \fn MERRA2::getTemperatureGradient()
//! \brief Gets the temperature gradient with respect to height.
//!
//! \returns Temperature gradient with respect to height \units{K}.


//! \brief Perform an update using surface height.
//!
//! This routine temporarily sets height to 0 and performs an update.
//! The height is returned to its previous value.
void MERRA2::updateSurface()
{
  surfaceUpdate = true;
  greal saveHeight = height;
  //P. White - Set to 0.0 to default to MERRA-2 surface data.  Would like to set the MERRA-2 surface height as the EG surface height instead of topographic map.
  height = 0.0;
  update();
  height = saveHeight;
  surfaceUpdate = false;
}

//! \brief Derive latitude and longitude interpolation indices and factors.
//!
//! This routine will use the sizes of the loaded MERRA2 data set to determine
//! interpolation indices and factors for latitude and longitude.
//!
//! \b Inputs
//! \arg #latitude
//! \arg #longitude
//! \arg #minimumLatitude
//! \arg #maximumLatitude
//! \arg #minimumLongitude
//! \arg #maximumLongitude
//! \arg MERRA2 data sizes
//!
//! \retval #latIndex
//! \retval #lonIndex
//! \retval #lonSurfaceIndex
//! \retval #latFactor
//! \retval #lonFactor
//! \retval #lonSurfaceFactor
void MERRA2::updateIndices()
{
  // Inverses of the grid size
  gridLat =(latFullSize - 1) / 180.0_deg;
  gridLon = lonFullSize / 360.0_deg;
  const greal gridSurfaceLon = lonSurfaceFullSize / 360.0_deg;

  //cout << "lat/lon = " << latitude << " " << longitude << '\n';
  greal clampLat = clamp(latitude, minimumLatitude, maximumLatitude);
  if (latitude != clampLat) {
    throw string("Error: Latitude is out of bounds.\n"
      "       Adjust MinimumLatitude and MaximumLatitude.\n");
  }
  greal shiftLat = clampLat + 90.0_deg;
  latIndex = size_t(shiftLat * gridLat) - latStartIndex;
  if (latIndex == latSize - 1) {
    latIndex = latSize - 2;
  }

  greal clampLon = clamp(longitude, minimumLongitude, maximumLongitude);
  if (minimumLongitude >= 0.0_deg && longitude != clampLon) {
    throw string("Error: Longitude is out of bounds.\n"
      "       Adjust MinimumLongitude and MaximumLongitude.\n");
  }
  if (minimumLongitude < 0.0_deg) {
    greal adjusted = longitude;
    if (longitude > maximumLongitude) {
      adjusted = longitude - 360.0_deg;
      clampLon = clamp(adjusted, minimumLongitude, maximumLongitude);
    }
    if (adjusted != clampLon) {
      throw string("Error: Longitude is out of bounds.\n"
        "       Adjust MinimumLongitude and MaximumLongitude.\n");
    }
  }
  if (longitude == 0.0 && maximumLongitude == 360.0_deg && minimumLongitude > 0.0_deg) {
    clampLon = 360.0_deg;
  }
  greal shiftLon = clampLon - lonShift;
  lonIndex = size_t(shiftLon * gridLon);
  greal shiftSurfaceLon = clampLon - lonSurfaceShift;
  lonSurfaceIndex = size_t(shiftSurfaceLon * gridSurfaceLon);
  //cout << "Index " << lonIndex << "  " << latIndex << '\n';

  lonFactor = shiftLon * gridLon - lonIndex ;
  lonSurfaceFactor = shiftSurfaceLon * gridSurfaceLon - lonSurfaceIndex;
  latFactor = shiftLat * gridLat - (latIndex + latStartIndex);
}

//! \brief Compute mean pressure, density, temperature, wind components, 
//! dewpoint temperature and their standard deviations, from the MERRA2 data
//!
//! This routine will read the MERRA2 data, if necessary, and perform interpolations based on the
//! current position.
//!
//! \b Inputs
//! \arg #position          
//! \arg MERRA2 data          
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
void MERRA2::update()
{
  // Initialize data, if necessary
  readM2File();

  constexpr greal onek = 1000.0;

  greal t[2][2][presSize], rho[2][2][presSize];
  greal u[2][2][presSize], v[2][2][presSize], w[2][2][presSize];
  greal sp[2][2][presSize], st[2][2][presSize], srho[2][2][presSize];
  greal su[2][2][presSize], sv[2][2][presSize];
  greal td[2][2][presSize], rh[2][2][presSize], spdavl[2][2][presSize]; // vp[2][2][presSize], 
  greal std[2][2][presSize], svp[2][2][presSize], srh[2][2][presSize], spdsdl[2][2][presSize];
  greal uvcorrl[2][2][presSize], h[2][2][presSize];

  greal usf[2][2], susf[2][2], vsf[2][2], svsf[2][2], tsf[2][2], stsf[2][2], dsf[2][2], sdsf[2][2],
      psf[2][2], spsf[2][2], tdsf[2][2], stdsf[2][2], svpsf[2][2], ghsf[2][2],
      hsf[2][2], spdavsf[2][2], spdsdsf[2][2], uvcorsf[2][2], rhsf[2][2], srhsf[2][2], wsf[2][2]; // vpsf[2][2], 

  greal pz[2][2], tz[2][2], rhoz[2][2], uz[2][2], vz[2][2], wz[2][2];
  greal spz[2][2], stz[2][2], srhoz[2][2], suz[2][2], svz[2][2];
  greal tdz[2][2], rhn[2][2], spdavz[2][2];  //, vpn[2][2]
  greal stdz[2][2], svpn[2][2], srhn[2][2], spdsdz[2][2], h0p3[2][2], h0p1[2][2];
  greal uvcorrz[2][2];

  // Compute indices and interpolation factors.
  updateIndices();

  // Compute surface gravity and radius
  greal g0, re;
  getEffectiveRadius(latitude, g0, re);
  greal re1 = re * onek;

  greal latrad = getRadius(latitude);
  greal totrad = latrad + height;

  // Compute vertical winds from Montgomery stream function
  // Get distances dx, dy across lat-lon "square"
  greal dy = onek * TO_RADIANS * totrad / gridLat;
  greal dx = onek * TO_RADIANS * totrad / gridLon;
  //greal dx1 = onek * TO_RADIANS * totrad / gridSurfaceLon;

  // Get horizontal temperature gradients
  greal dtdx[presSize] = { 0.0 };
  greal dtdy[presSize] = { 0.0 };
  for (size_t iPres = 0; iPres < presSize; iPres++) {

    dtdx[iPres] = (temp[getIndex(iPres, latIndex, lonIndex + 1)]
                   - temp[getIndex(iPres, latIndex, lonIndex)])
                  / dx;
    dtdy[iPres] = (temp[getIndex(iPres, latIndex + 1, lonIndex)]
                   - temp[getIndex(iPres, latIndex, lonIndex)])
                  / dy;
  }
  for (size_t iLat = 0; iLat < 2; iLat++) {
    size_t jLat = latIndex + iLat;
    for (size_t iLon = 0; iLon < 2; iLon++) {
      size_t jLon = lonIndex + iLon;
      for (size_t iPres = 0; iPres < presSize; iPres++) {
        size_t index3d = getIndex(iPres, jLat, jLon);

        rho[iLon][iLat][iPres] = dens[index3d];
        h[iLon][iLat][iPres] = hgt[index3d] * 1.0e-3;
        u[iLon][iLat][iPres] = uwnd[index3d];
        v[iLon][iLat][iPres] = vwnd[index3d];
        t[iLon][iLat][iPres] = temp[index3d];
        td[iLon][iLat][iPres] = dewp[index3d];
        rh[iLon][iLat][iPres] = rhum[index3d];
        //vp[iLon][iLat][iPres] = vprs[index3d];
        srho[iLon][iLat][iPres] = sden[index3d] * densityPerturbationScale;
        sp[iLon][iLat][iPres] = spres[index3d] * pressurePerturbationScale * 100.0;
        st[iLon][iLat][iPres] = stmp[index3d] * temperaturePerturbationScale;
        std[iLon][iLat][iPres] = sdewp[index3d];
        su[iLon][iLat][iPres] = suwd[index3d] * ewWindPerturbationScale;
        sv[iLon][iLat][iPres] = svwd[index3d] * nsWindPerturbationScale;
        uvcorrl[iLon][iLat][iPres] = uvcor[index3d];
        svp[iLon][iLat][iPres] = svprs[index3d];
        spdavl[iLon][iLat][iPres] = spdav[index3d];
        spdsdl[iLon][iLat][iPres] = spdsd[index3d];
        srh[iLon][iLat][iPres] = srhum[index3d];
      }

      // Get vertical wind at "corners" of lat-lon "square"
      for (size_t iPres = 0; iPres < presSize; iPres++) {
	    // This clampSize ensures cPres a max of presSize-2.
        size_t cPres = clampSize(iPres, presSize - 1);
        greal dtdz1 = (t[iLon][iLat][cPres + 1] - t[iLon][iLat][cPres])
                      / (onek * (h[iLon][iLat][cPres + 1] - h[iLon][iLat][cPres]));
        greal cp = 7.0 * getGasConstant(td[iLon][iLat][iPres], pn[iPres]) / 2.0;
        greal dtdz2 = max(1.0e-4, dtdz1 + g0 / cp);
        w[iLon][iLat][iPres] = -(u[iLon][iLat][iPres] * dtdx[iPres] 
                              + v[iLon][iLat][iPres] * dtdy[iPres]) / dtdz2;
        //Check for nan vertical wind and set to 0.0.  Bug propagated for 2.0 resolution data. - P. White
        if (isnan(w[iLon][iLat][iPres])) {
          w[iLon][iLat][iPres] = 0.0;
        }
      }
    }
  }

  for (size_t iLat = 0; iLat < 2; iLat++) {
    size_t jLat = latIndex + iLat;
    for (size_t iLon = 0; iLon < 2; iLon++) {
      size_t jLon = lonIndex + iLon;
      size_t jSurfaceLon = lonSurfaceIndex + iLon;

      size_t index2d = getIndex(jLat, jLon);
      size_t index2dSurface = getSurfaceIndex(jLat, jSurfaceLon);

      usf[iLon][iLat] = usfc[index2dSurface];
      susf[iLon][iLat] = susfc[index2dSurface];
      vsf[iLon][iLat] = vsfc[index2dSurface];
      svsf[iLon][iLat] = svsfc[index2dSurface];
      tsf[iLon][iLat] = tsfc[index2dSurface];
      stsf[iLon][iLat] = stsfc[index2dSurface];
      dsf[iLon][iLat] = dsfc[index2dSurface];
      sdsf[iLon][iLat] = sdsfc[index2dSurface];
      psf[iLon][iLat] = psfc[index2d];
      spsf[iLon][iLat] = spsfc[index2d];
      tdsf[iLon][iLat] = tdsfc[index2d];
      stdsf[iLon][iLat] = stdsfc[index2d];
      //vpsf[iLon][iLat] = vpsfc[index2d];
      svpsf[iLon][iLat] = svpsfc[index2d];
      ghsf[iLon][iLat] = ghsfc[index2d] * 1.0e-3;
      hsf[iLon][iLat] = ghsf[iLon][iLat] * re1 / (g0 * re1 - ghsf[iLon][iLat]);
      spdavsf[iLon][iLat] = spdavsfc[index2dSurface];
      spdsdsf[iLon][iLat] = spdsdsfc[index2dSurface];
      uvcorsf[iLon][iLat] = uvcorsfc[index2dSurface];
      rhsf[iLon][iLat] = rhsfc[index2d];
      srhsf[iLon][iLat] = srhsfc[index2d];
      wsf[iLon][iLat] = 0.0;
    }
  }

  greal dtdz[2][2];
  for (size_t iLon = 0; iLon < 2; iLon++) {
    for (size_t iLat = 0; iLat < 2; iLat++) {
      size_t lowerHeight = 0, upperHeight = 0;
      // Below the surface
      // height <= hsf
      if (height <= hsf[iLon][iLat]) {
        rhoz[iLon][iLat] = dsf[iLon][iLat];
        pz[iLon][iLat] = psf[iLon][iLat];
        uz[iLon][iLat] = usf[iLon][iLat];
        vz[iLon][iLat] = vsf[iLon][iLat];
        tz[iLon][iLat] = tsf[iLon][iLat];
        tdz[iLon][iLat] = tdsf[iLon][iLat];
        //vpn[iLon][iLat] = vpsf[iLon][iLat];
        rhn[iLon][iLat] = rhsf[iLon][iLat];
        srhoz[iLon][iLat] = sdsf[iLon][iLat];
        spz[iLon][iLat] = spsf[iLon][iLat];
        stz[iLon][iLat] = stsf[iLon][iLat];
        stdz[iLon][iLat] = stdsf[iLon][iLat];
        suz[iLon][iLat] = susf[iLon][iLat];
        svz[iLon][iLat] = svsf[iLon][iLat];
        svpn[iLon][iLat] = svpsf[iLon][iLat];
        srhn[iLon][iLat] = srhsf[iLon][iLat];
        spdavz[iLon][iLat] = spdavsf[iLon][iLat];
        spdsdz[iLon][iLat] = spdsdsf[iLon][iLat];
        uvcorrz[iLon][iLat] = uvcorsf[iLon][iLat];
        wz[iLon][iLat] = 0.0;
        dtdz[iLon][iLat] = 0.0;
      }
      // Above the surface, but below all h values
      // hsf < height < h[0]
      else if (hsf[iLon][iLat] < height && height < h[iLon][iLat][lowerHeight]) {
        greal dh;
        // Density scale height, km
        if (abs(hsf[iLon][iLat] - h[iLon][iLat][lowerHeight]) == 0.0) {
          dtdz[iLon][iLat] = 0.0;
          dh = 0.0;
        }
        else {
          dh = (height - hsf[iLon][iLat]) / (h[iLon][iLat][lowerHeight] - hsf[iLon][iLat]);
          dtdz[iLon][iLat] = (t[iLon][iLat][lowerHeight] - tsf[iLon][iLat])
            / (h[iLon][iLat][lowerHeight] - hsf[iLon][iLat]);
        }
        greal hden;
        if (abs(dsf[iLon][iLat] - rho[iLon][iLat][lowerHeight]) == 0.0
          || abs(hsf[iLon][iLat] - h[iLon][iLat][lowerHeight]) == 0.0) {
          hden = 1.0;
          rhoz[iLon][iLat] = dsf[iLon][iLat];
        }
        else {
          hden = (h[iLon][iLat][lowerHeight] - hsf[iLon][iLat])
            / log(dsf[iLon][iLat] / rho[iLon][iLat][lowerHeight]);
          rhoz[iLon][iLat] = dsf[iLon][iLat] * exp((hsf[iLon][iLat] - height) / hden);
        }
        Interpolator interpHeight(dh);

        tz[iLon][iLat] = interpHeight.linear(tsf[iLon][iLat], t[iLon][iLat][lowerHeight]);
        greal r1 = psf[iLon][iLat] / (dsf[iLon][iLat] * tsf[iLon][iLat]);
        greal r2 = pn[lowerHeight] / (rho[iLon][iLat][lowerHeight] * t[iLon][iLat][lowerHeight]);
        greal rz = interpHeight.linear(r1, r2);

        // Pressure (mb) from gas law
        pz[iLon][iLat] = rhoz[iLon][iLat] * rz * tz[iLon][iLat];
        uz[iLon][iLat] = interpHeight.linear(usf[iLon][iLat], u[iLon][iLat][lowerHeight]);
        vz[iLon][iLat] = interpHeight.linear(vsf[iLon][iLat], v[iLon][iLat][lowerHeight]);
        wz[iLon][iLat] = interpHeight.linear(wsf[iLon][iLat], w[iLon][iLat][lowerHeight]);
        tdz[iLon][iLat] = interpHeight.linear(tdsf[iLon][iLat], td[iLon][iLat][lowerHeight]);
        //vpn[iLon][iLat] = interpHeight.log(vpsf[iLon][iLat], vp[iLon][iLat][lowerHeight]);
        rhn[iLon][iLat] = interpHeight.linear(rhsf[iLon][iLat], rh[iLon][iLat][lowerHeight]);
        srhoz[iLon][iLat] = interpHeight.linear(sdsf[iLon][iLat], srho[iLon][iLat][lowerHeight]);
        spz[iLon][iLat] = interpHeight.linear(spsf[iLon][iLat], sp[iLon][iLat][lowerHeight]);
        stz[iLon][iLat] = interpHeight.linear(stsf[iLon][iLat], st[iLon][iLat][lowerHeight]);
        stdz[iLon][iLat] = interpHeight.linear(stdsf[iLon][iLat], std[iLon][iLat][lowerHeight]);
        suz[iLon][iLat] = interpHeight.linear(susf[iLon][iLat], su[iLon][iLat][lowerHeight]);
        svz[iLon][iLat] = interpHeight.linear(svsf[iLon][iLat], sv[iLon][iLat][lowerHeight]);
        svpn[iLon][iLat] = interpHeight.linear(svpsf[iLon][iLat], svp[iLon][iLat][lowerHeight]);
        srhn[iLon][iLat] = interpHeight.linear(srhsf[iLon][iLat], srh[iLon][iLat][lowerHeight]);
        spdavz[iLon][iLat] = interpHeight.linear(spdavsf[iLon][iLat], spdavl[iLon][iLat][lowerHeight]);
        spdsdz[iLon][iLat] = interpHeight.linear(spdsdsf[iLon][iLat], spdsdl[iLon][iLat][lowerHeight]);
        uvcorrz[iLon][iLat] = interpHeight.linear(uvcorsf[iLon][iLat], uvcorrl[iLon][iLat][lowerHeight]);
      }
      // Above the surface and in the h value range
      // h[0] <= hsf < height
      else {
        // Find such that h[lower] <= height <= h[upper]
        for (size_t iPres = 0; iPres < presSize - 1; iPres++) {
          lowerHeight = iPres;
          upperHeight = iPres + 1;

          if (height <= h[iLon][iLat][upperHeight])
            break;
        }

        //P. White - Case to handle surface level between a lower and upper pressure level
        // h[lower] < hsf < height < h[upper]
        if (hsf[iLon][iLat] > h[iLon][iLat][lowerHeight] && hsf[iLon][iLat] < h[iLon][iLat][upperHeight])
        {
          greal dh;
          if (abs(hsf[iLon][iLat] - h[iLon][iLat][upperHeight]) == 0.0) {
            dh = 0.0;
            dtdz[iLon][iLat] = 0.0;
          }
          else {
            dh = (height - hsf[iLon][iLat]) / (h[iLon][iLat][upperHeight] - hsf[iLon][iLat]);
            dtdz[iLon][iLat] = (t[iLon][iLat][upperHeight] - tsf[iLon][iLat])
              / (h[iLon][iLat][upperHeight] - hsf[iLon][iLat]);
          }
          greal hden;
          if (abs(dsf[iLon][iLat] - rho[iLon][iLat][upperHeight]) == 0.0
            || abs(hsf[iLon][iLat] - h[iLon][iLat][upperHeight]) == 0.0) {
            hden = 1.0;
            rhoz[iLon][iLat] = dsf[iLon][iLat];
          }
          else {
            hden = (h[iLon][iLat][upperHeight] - hsf[iLon][iLat])
              / log(dsf[iLon][iLat] / rho[iLon][iLat][upperHeight]);
            // Logarithmic interpolation for density
            rhoz[iLon][iLat] = dsf[iLon][iLat] * exp((hsf[iLon][iLat] - height) / hden);
          }

          Interpolator interpHeight(dh);
          // Linear interpolation on temperature and gas constant
          tz[iLon][iLat] = interpHeight.linear(tsf[iLon][iLat], t[iLon][iLat][upperHeight]);
          greal r1 = psf[iLon][iLat] / (dsf[iLon][iLat] * tsf[iLon][iLat]);
          greal r2 = pn[upperHeight] / (rho[iLon][iLat][upperHeight] * t[iLon][iLat][upperHeight]);
          greal rz = interpHeight.linear(r1, r2);
          // Pressure (mb) from gas law
          pz[iLon][iLat] = rhoz[iLon][iLat] * rz * tz[iLon][iLat];
          uz[iLon][iLat] = interpHeight.linear(usf[iLon][iLat], u[iLon][iLat][upperHeight]);
          vz[iLon][iLat] = interpHeight.linear(vsf[iLon][iLat], v[iLon][iLat][upperHeight]);
          wz[iLon][iLat] = interpHeight.linear(0.0, w[iLon][iLat][upperHeight]);
          tdz[iLon][iLat] = interpHeight.linear(tdsf[iLon][iLat], td[iLon][iLat][upperHeight]);
          //vpn[iLon][iLat] = interpHeight.log(vpsf[iLon][iLat], vp[iLon][iLat][upperHeight]);
          rhn[iLon][iLat] = interpHeight.linear(rhsf[iLon][iLat], rh[iLon][iLat][upperHeight]);
          srhoz[iLon][iLat] = interpHeight.linear(sdsf[iLon][iLat], srho[iLon][iLat][upperHeight]);
          spz[iLon][iLat] = interpHeight.linear(spsf[iLon][iLat], sp[iLon][iLat][upperHeight]);
          stz[iLon][iLat] = interpHeight.linear(stsf[iLon][iLat], st[iLon][iLat][upperHeight]);
          stdz[iLon][iLat] = interpHeight.linear(stdsf[iLon][iLat], std[iLon][iLat][upperHeight]);
          suz[iLon][iLat] = interpHeight.linear(susf[iLon][iLat], su[iLon][iLat][upperHeight]);
          svz[iLon][iLat] = interpHeight.linear(svsf[iLon][iLat], sv[iLon][iLat][upperHeight]);
          svpn[iLon][iLat] = interpHeight.linear(svpsf[iLon][iLat], svp[iLon][iLat][upperHeight]);
          srhn[iLon][iLat] = interpHeight.linear(srhsf[iLon][iLat], srh[iLon][iLat][upperHeight]);
          spdavz[iLon][iLat]
            = interpHeight.linear(spdavsf[iLon][iLat], spdavl[iLon][iLat][upperHeight]);
          spdsdz[iLon][iLat]
            = interpHeight.linear(spdsdsf[iLon][iLat], spdsdl[iLon][iLat][upperHeight]);
          uvcorrz[iLon][iLat]
            = interpHeight.linear(uvcorsf[iLon][iLat], uvcorrl[iLon][iLat][upperHeight]);

          //          cout << height << "  " << usf[iLon][iLat] << "  " << u[iLon][iLat][upperHeight] << "  "
          //               << hsf[iLon][iLat] << "  " << h[iLon][iLat][upperHeight] << " " << upperHeight
          //               << '\n';
        }
        //  hsf < h[lower] < height < h[upper]
        else {
            greal dh;
            if (abs(h[iLon][iLat][lowerHeight] - h[iLon][iLat][upperHeight]) == 0.0) {
              dh = 0.0;
              dtdz[iLon][iLat] = 0.0;
            }
            else {
              dh = (height - h[iLon][iLat][lowerHeight])
                / (h[iLon][iLat][upperHeight] - h[iLon][iLat][lowerHeight]);
              dtdz[iLon][iLat] = (t[iLon][iLat][upperHeight] - t[iLon][iLat][lowerHeight])
                / (h[iLon][iLat][upperHeight] - h[iLon][iLat][lowerHeight]);
            }
            greal hden;
            if (abs(rho[iLon][iLat][lowerHeight] - rho[iLon][iLat][upperHeight]) == 0.0
              || abs(h[iLon][iLat][lowerHeight] - h[iLon][iLat][upperHeight]) == 0.0) {
              hden = 1.0;
              rhoz[iLon][iLat] = rho[iLon][iLat][lowerHeight];
            }
            else {
              hden = (h[iLon][iLat][upperHeight] - h[iLon][iLat][lowerHeight])
                / log(rho[iLon][iLat][lowerHeight] / rho[iLon][iLat][upperHeight]);
              // Logarithmic interpolation for density
              rhoz[iLon][iLat] = rho[iLon][iLat][lowerHeight]
                * exp((h[iLon][iLat][lowerHeight] - height) / hden);
            }

            Interpolator interpHeight(dh);
            // Linear interpolation on temperature and gas constant
            tz[iLon][iLat] = interpHeight.linear(t[iLon][iLat][lowerHeight], t[iLon][iLat][upperHeight]);
            greal r1 = pn[lowerHeight] / (rho[iLon][iLat][lowerHeight] * t[iLon][iLat][lowerHeight]);
            greal r2 = pn[upperHeight] / (rho[iLon][iLat][upperHeight] * t[iLon][iLat][upperHeight]);
            greal rz = interpHeight.linear(r1, r2);
            // Pressure (mb) from gas law
            pz[iLon][iLat] = rhoz[iLon][iLat] * rz * tz[iLon][iLat];
            uz[iLon][iLat] = interpHeight.linear(u[iLon][iLat][lowerHeight], u[iLon][iLat][upperHeight]);
            vz[iLon][iLat] = interpHeight.linear(v[iLon][iLat][lowerHeight], v[iLon][iLat][upperHeight]);
            wz[iLon][iLat] = interpHeight.linear(w[iLon][iLat][lowerHeight], w[iLon][iLat][upperHeight]);
            tdz[iLon][iLat] = interpHeight.linear(td[iLon][iLat][lowerHeight], td[iLon][iLat][upperHeight]);
            //vpn[iLon][iLat] = interpHeight.log(vp[iLon][iLat][lowerHeight], vp[iLon][iLat][upperHeight]);
            rhn[iLon][iLat] = interpHeight.linear(rh[iLon][iLat][lowerHeight], rh[iLon][iLat][upperHeight]);
            srhoz[iLon][iLat] = interpHeight.linear(srho[iLon][iLat][lowerHeight], srho[iLon][iLat][upperHeight]);
            spz[iLon][iLat] = interpHeight.linear(sp[iLon][iLat][lowerHeight], sp[iLon][iLat][upperHeight]);
            stz[iLon][iLat] = interpHeight.linear(st[iLon][iLat][lowerHeight], st[iLon][iLat][upperHeight]);
            stdz[iLon][iLat] = interpHeight.linear(std[iLon][iLat][lowerHeight], std[iLon][iLat][upperHeight]);
            suz[iLon][iLat] = interpHeight.linear(su[iLon][iLat][lowerHeight], su[iLon][iLat][upperHeight]);
            svz[iLon][iLat] = interpHeight.linear(sv[iLon][iLat][lowerHeight], sv[iLon][iLat][upperHeight]);
            svpn[iLon][iLat] = interpHeight.linear(svp[iLon][iLat][lowerHeight], svp[iLon][iLat][upperHeight]);
            srhn[iLon][iLat] = interpHeight.linear(srh[iLon][iLat][lowerHeight], srh[iLon][iLat][upperHeight]);
            spdavz[iLon][iLat] = interpHeight.linear(spdavl[iLon][iLat][lowerHeight], spdavl[iLon][iLat][upperHeight]);
            spdsdz[iLon][iLat] = interpHeight.linear(spdsdl[iLon][iLat][lowerHeight], spdsdl[iLon][iLat][upperHeight]);
            uvcorrz[iLon][iLat] = interpHeight.linear(uvcorrl[iLon][iLat][lowerHeight], uvcorrl[iLon][iLat][upperHeight]);

            //For uvcorr nan cases interpolate between surface to upper height.  Propagated from 2.0 resolution data. - P. White 
            if (isnan(uvcorrl[iLon][iLat][lowerHeight])) {
              uvcorrz[iLon][iLat]
                = interpHeight.linear(uvcorsf[iLon][iLat], uvcorrl[iLon][iLat][upperHeight]);
            }
            //        cout << height << "  " << u[iLon][iLat][lowerHeight] << "  " << u[iLon][iLat][upperHeight] << "  "
            //             << h[iLon][iLat][lowerHeight] << "  " << h[iLon][iLat][upperHeight] << " " << lowerHeight << "  " << upperHeight << '\n';
          }
      }

      h0p3[iLon][iLat] = h[iLon][iLat][40];
      h0p1[iLon][iLat] = h[iLon][iLat][41];
    }
  }

//  cout << "Lat/Lon Factors " << latFactor << "  " << lonFactor << '\n';

  Interpolator interpLatLon(lonFactor, latFactor);
  pressure = interpLatLon.linear(pz);
  temperature = interpLatLon.linear(tz);
  dewPoint = interpLatLon.linear(tdz);
  density = interpLatLon.linear(rhoz);
  ewWind = interpLatLon.linear(uz);
  nsWind = interpLatLon.linear(vz);
  verticalWind = interpLatLon.linear(wz);
  pressureStandardDeviation = square( (interpLatLon.linear(spz)) / pressure);
  temperatureStandardDeviation = square(interpLatLon.linear(stz) / temperature);
  densityStandardDeviation = square(interpLatLon.linear(srhoz) / density);
  dewPointSD = interpLatLon.linear(stdz);
  ewStandardDeviation = interpLatLon.linear(suz);
  nsStandardDeviation = interpLatLon.linear(svz);
  temperatureGradient = interpLatLon.linear(dtdz);
  relativeHumidity = interpLatLon.linear(rhn);
  relativeHumiditySD = interpLatLon.linear(srhn);
//  vaporPressure = interpLatLon.linear(vpn);
  vaporPressure = 610.94 * exp((17.625 * (dewPoint - 273.15)) / (243.04 + (dewPoint - 273.15)));
  vaporPressureSD = interpLatLon.linear(svpn);
  windSpeed = interpLatLon.linear(spdavz);
  windSpeedStandardDeviation = interpLatLon.linear(spdsdz);
  windCorrelation = interpLatLon.linear(uvcorrz);
  hp1 = interpLatLon.linear(h0p1);
  hp3 = interpLatLon.linear(h0p3);

//  cout << "latitude " << latitude << "  " << height << "  " << temperature << "  " << density
//       << "  " << pressure << "  " << vaporPressure << "  " << relativeHumidity << " " << dewPoint << '\n';
}

//! \brief Compute gas constant from dewpoint and pressure
//!
//! The function computes the gas constant (R) based on dewpont temperature and pressure.
//!
//! \param tdc Dewpoint temperature \units{\text{degrees K}}.
//! \param pres Pressure \units{Pa}.
//!
//! \retval The gas constant (R).
greal MERRA2::getGasConstant(greal tdc, greal pres)
{
  constexpr greal omeps = 0.37803; // one minus epsilon
  constexpr greal onec = 100.0;

  greal e = wexler(tdc, pres) / onec;

  greal gascon_out = dryAirGasConstant / (1.0 - omeps * e / pres);

  return gascon_out;
}

//! \brief Calculate heights of 0.3 mb and 0.1 mb pressure levels.
//!
//! Compute geometric height of lowest 0.3 mb level and highest 0.1 mb
//! level from four profiles surrounding given lat-lon.
//!
//! \param lat  Latitude.
//! \param lon  Longitude.
//! \param[out] lowest0p3mb  A greal.
//! \param[out] highest0p1mb   A greal.
//!
//! \retval highest0p1mb  Height \units{km} of the highest 0.1 mb pressure level.
//! \retval lowest0p3mb   Height \units{km} of the lowest 0.3 mb pressure level.
void MERRA2::getHeights(greal lat, greal lon, greal& lowest0p3mb, greal& highest0p1mb)
{
  // Initialize data, if necessary
  readM2File();

  latitude = lat;
  longitude = lon;

  // Compute indices and interpolation factors.
  updateIndices();

  greal hLow[2][2], hHigh[2][2];
  for (size_t iLat = 0; iLat < 2; iLat++) {
    size_t jLat = latIndex + iLat;
    for (size_t iLon = 0; iLon < 2; iLon++) {
      size_t jLon = lonIndex + iLon;
      size_t index3dLow = getIndex(presSize - 2, jLat, jLon);
      size_t index3dHigh = getIndex(presSize - 1, jLat, jLon);
      hLow[iLon][iLat] = hgt[index3dLow] * 1.0e-3;
      hHigh[iLon][iLat] = hgt[index3dHigh] * 1.0e-3;
    }
  }

  Interpolator interpLatLon(lonFactor, latFactor);

  highest0p1mb = interpLatLon.linear(hHigh);
  lowest0p3mb = interpLatLon.linear(hLow);
}

//! \fn MERRA2::getIndex(size_t ipres, size_t ilat, size_t ilon)
//! \brief Generates a 1D index from 3 indices.
//! \param ipres A pressure index.
//! \param ilat  A latitude index.
//! \param ilon  A longitude index.
//! \returns A 1D index.

//! \fn MERRA2::getIndex(size_t ilat, size_t ilon)
//! \brief Generates a 1D index from 2 indices.
//! \param ilat  A latitude index.
//! \param ilon  A longitude index.
//! \returns A 1D index.

//! \fn MERRA2::getSurfaceIndex(size_t ilat, size_t ilon)
//! \brief Generates a 1D index from 2 indices for surface data.
//! \param ilat  A latitude index.
//! \param ilon  A longitude index.
//! \returns A 1D index.
} // namespace GRAM