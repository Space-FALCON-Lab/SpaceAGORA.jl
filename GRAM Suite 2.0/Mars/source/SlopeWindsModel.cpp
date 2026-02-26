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
#include "SlopeWindsModel.h"
#include "MOLATopography.h"

using namespace std;

namespace GRAM {


//! \copydoc Atmosphere::Atmosphere()
SlopeWindsModel::SlopeWindsModel()
{
}

//! \fn  SlopeWindsModel::SlopeWindsModel(const SlopeWindsModel& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  SlopeWindsModel::~SlopeWindsModel()
//! \copydoc Atmosphere::~Atmosphere()

//! \fn  SlopeWindsModel::setPosition()
//! \copydoc Atmosphere::setPosition()

//! \fn  SlopeWindsModel::setSolarTime()
//! \brief Set the solar time.

//! \fn  SlopeWindsModel::setTemperature()
//! \brief Set the temperature.

//! \fn  SlopeWindsModel::setMeanWindsScale()
//! \brief Set the mean winds scale.

//! \fn  SlopeWindsModel::setBoundaryLayerWindsScale()
//! \brief Set the boundary layer winds scale.

//! \fn  SlopeWindsModel::setWinds()
//! \brief Set the winds (before boundary layer effects).

//! \fn  SlopeWindsModel::setDailyAverageWinds()
//! \brief Set the average winds (before boundary layer effects).

//! \fn  SlopeWindsModel::getWinds()
//! \brief Get the winds with boundary layer effects.

//! \fn  SlopeWindsModel::getDailyAverageWinds()
//! \brief Get the average winds with boundary layer effects.

//! \fn  SlopeWindsModel::setAreoidRadiusCallback()
//! \copydoc MOLATopography::setAreoidRadiusCallback()

//! \fn  SlopeWindsModel::setTopographicHeightCallback()
//! \copydoc MOLATopography::setTopographicHeightCallback()

//! \fn  SlopeWindsModel::setCallbackData()
//! \copydoc MOLATopography::setCallbackData()


//! \brief Interface for the slope winds computations.
//!
//! This defines the interface for the slope winds computations.
//! Precede this function by setting all parameters and current winds.
//! The EW, NS, and vertical winds are updated.
void SlopeWindsModel::update()
{
  updateSlopes();
  updateWinds();
}

//! \brief Computes the EW and NS slopes.
//!
//! Computes the EW and NS slopes as a delta of the topographic height over a
//! delta of the surface distance.  The deltas are taken over half degree change
//! in latitude or longitude (as appropriate) centered at the current position.
//!
//! \b Inputs
//! \arg #latitude
//! \arg #longitude
//! \arg #latitudeRadius
//!
//! \retval #ewSlope
//! \retval #nsSlope
void SlopeWindsModel::updateSlopes()
{
  // NORTH-SOUTH Slope computations
  // current latitude +- 0.25 degrees to be used to compute northward delta
  greal latPlus = latitude + 0.25_deg;
  greal latMinus = latitude - 0.25_deg;

  // initially, no change in longitude needed to compute northward delta
  greal lonPlus = longitude;
  greal lonMinus = longitude;

  //=================================================================================================//
  // adjust latPlus and lonPlus if latPlus crosses N pole //
  //-------------------------------------------------------------------------------------------------//
  if (latPlus > 90.0_deg) {
    // adjust latPlus (= 90 minus |excess above 90| )
    latPlus = 180.0_deg - latPlus;
    // adjust lonPlus (= rotated 180 degrees)
    lonPlus += 180.0_deg;
    // ensure lonPlus is in range [0,360)
    if (lonPlus >= 360.0_deg)
      lonPlus -= 360.0_deg;
  }
  //=================================================================================================//
  // adjust latMinus and lonMinus if latMinus crosses S pole                                         //
  //-------------------------------------------------------------------------------------------------//
  else if (latMinus < -90.0_deg) {
    // adjust latMinus (= -90 plus |excess below -90| )
    latMinus = -180.0_deg - latMinus;
    // adjust lonPlus (= rotated 180 degrees)
    lonMinus += 180.0_deg;
    // ensure lonMinus is in range [0,360)
    if (lonMinus >= 360.0_deg)
      lonMinus -= 360.0_deg;
  }

  //=================================================================================================//
  //  get north-south slope at current position                                                      //
  //-------------------------------------------------------------------------------------------------//
  MOLATopography mola;
  mola.setAreoidRadiusCallback(getAreoidRadiusCallback);
  mola.setTopographicHeightCallback(getTopographicHeightCallback);
  mola.setCallbackData(callbackDataPointer);
  // get mola topo height at latPlus, lonPlus
  greal topoHeightPlus = mola.getTopographicHeight(latPlus, lonPlus);

  // get mola topo height at latMinus, lonMinus
  greal topoHeightMinus = mola.getTopographicHeight(latMinus, lonMinus);

  // N-S slope = change in height / change in distance
  nsSlope = (topoHeightPlus - topoHeightMinus) / (toRadians(0.5_deg) * latitudeRadius);

  //=================================================================================================//
  //  get east-west slope at current position                                                        //
  //-------------------------------------------------------------------------------------------------//
  // if lat is within 0.25 degrees of N or S pole, set E-W slope to 0
  if (abs(latitude) > 89.75_deg) {
    ewSlope = 0.0;
  }
  else { // else (lat not within 0.25 degrees of N or S pole) ...
    // current longitude +- 0.25 degrees to be used to compute eastward slope
    lonPlus = longitude + 0.25_deg;
    lonMinus = longitude - 0.25_deg;

    // ensure lonMinus is in range [0,360)
    if (lonMinus < 0.0_deg) {
      lonMinus += 360.0_deg;
    }
    // ensure lonPlus is in range [0,360)
    if (lonPlus >= 360.0_deg) {
      lonPlus -= 360.0_deg;
    }

    // get mola topo height at currrent latitude, lonPlus
    topoHeightPlus = mola.getTopographicHeight(latitude, lonPlus);

    // get mola topo height at currrent latitude, lonMinus
    topoHeightMinus = mola.getTopographicHeight(latitude, lonMinus);

    // E-W slope
    ewSlope = (topoHeightPlus - topoHeightMinus) / (toRadians(0.5_deg) * cos(toRadians(latitude)) * latitudeRadius);
  }
}

//! \brief Adjusts EW, NS winds with boundary layer winds.  Also, computes vertical winds.
//!
//! Boundary layer winds (EW, NS, and vertical) are computed based on solar time.
//! EW, NS winds are adjusted with boundary layer winds and user supplied scaling factors.
//! Vertical winds are computed from the boundary layer winds and user scaling.
//!
//! \b Inputs
//! \arg #solarTime
//! \arg #temperature
//! \arg #meanWindsScale
//! \arg #boundaryLayerWindsScale
//!
//! \retval #ewWind
//! \retval #nsWind
//! \retval #verticalWind
void SlopeWindsModel::updateWinds()
{
  // Get the boundary layer slope winds
  greal ewBoundaryLayerWinds = 0.0;
  greal nsBoundaryLayerWinds = 0.0;
  greal verticalBoundaryLayerWinds = 0.0;
  slopeWinds(solarTime, ewBoundaryLayerWinds, nsBoundaryLayerWinds, verticalBoundaryLayerWinds);

  // limit for wind speeds to approximately the speed of sound over sqrt(2).
  greal sosLimit = 0.7 * speedOfSound;

  // Add in boundary layer winds and apply user scaling.
  ewWind = meanWindsScale * ewWind + boundaryLayerWindsScale * ewBoundaryLayerWinds;
  nsWind = meanWindsScale * nsWind + boundaryLayerWindsScale * nsBoundaryLayerWinds;
  verticalWind = boundaryLayerWindsScale * verticalBoundaryLayerWinds;

  // limit winds to 0.7 * speed of sound
  ewWind = clamp(ewWind, -sosLimit, sosLimit);
  nsWind = clamp(nsWind, -sosLimit, sosLimit);

  //=================================================================================================//
  //  Compute daily average boundary layer wind effects and add them into the daily average winds.   //
  //-------------------------------------------------------------------------------------------------//
  greal ewSlopeWindSum = 0.0;
  greal nsSlopeWindSum = 0.0;

  // loop over 24 hours in a sol, with a 2-hour step (12 samples)
  for (int itime = 0; itime < 23; itime += 2) {
    // wind components
    greal u, v, w;
    // Get the slope wind components for itime.
    slopeWinds(greal(itime), u, v, w);

    // Sum up the components
    ewSlopeWindSum += u;
    nsSlopeWindSum += v;
  }

  // divide running of hourly values to get daily average
  greal ewBoundaryLayerWindsMean = ewSlopeWindSum / 12.0;
  greal nsBoundaryLayerWindsMean = nsSlopeWindSum / 12.0;

  // Add in boundary layer winds with user scaling.
  ewWindDailyAverage = meanWindsScale * ewWindDailyAverage + boundaryLayerWindsScale * ewBoundaryLayerWindsMean;
  nsWindDailyAverage = meanWindsScale * nsWindDailyAverage + boundaryLayerWindsScale * nsBoundaryLayerWindsMean;

  // limit winds to 0.7 * speed of sound
  ewWindDailyAverage = clamp(ewWindDailyAverage, -sosLimit, sosLimit);
  nsWindDailyAverage = clamp(nsWindDailyAverage, -sosLimit, sosLimit);
}

//! \brief Computes slope winds in boundary layer.                                                           //
//!
//! Analytical slope winds solution from Ye, Segal, and Pielke, J.Atmos.Sci., 47, 1990.  (YSP).
//!
//! \param time The local solar time in hours.
//! \param[out] ewSlopeWinds A greal.
//! \param[out] nsSlopeWinds A greal.
//! \param[out] verticalSlopeWinds A greal.
//!
//! \b Inputs
//! \arg #surfaceHeight
//! \arg #height
//! \arg #latitude
//! \arg #ewSlope
//! \arg #nsSlope
//!
//! \retval ewSlopeWinds The EW slope winds.
//! \retval nsSlopeWinds The NS slope winds.
//! \retval verticalSlopeWinds The vertical slope winds.
void SlopeWindsModel::slopeWinds(greal time, greal& ewSlopeWinds, greal& nsSlopeWinds, greal& verticalSlopeWinds)
{
  const greal omega = 7.0777e-5;  // rotational speed of Mars (radians/sec)
  const greal c0 = 0.07;          // constant coefficient from YSP (m/sec)

  // height above local topo surface
  greal hgt_sfc = height - surfaceHeight;

  // time of day cosine factor
  greal cos_fac = cos(toRadians(15.0_deg * (time - 15.0)));

  // boundary layer depth
  greal depth = 3.5 + cos_fac;

  //=================================================================================================//
  //  if height is at or below topographic surface, or above boundary layer(~4.5 km),
  //  set slope winds to zero
  //-------------------------------------------------------------------------------------------------//
  if (hgt_sfc <= 0.0 || hgt_sfc > depth) {
    ewSlopeWinds = 0.0;
    nsSlopeWinds = 0.0;
    verticalSlopeWinds = 0.0;
    return;
  }

  //=================================================================================================//
  //  height is between surface and ~ 4.5 km, define initial constants
  //-------------------------------------------------------------------------------------------------//
  // assumed diurnal variation in Qs
  greal qslam = 550.0 + 150.0 * cos_fac;

  // normalized BL height
  greal xi = hgt_sfc / depth;

  //=================================================================================================//
  // compute up-slope and cross-slope forces very near equator (where coriolis is negligible)
  //-------------------------------------------------------------------------------------------------//
  greal f_upslope = 0.0; // up-slope force
  greal f_xslope = 0.0;  // cross-slope force
  // if lat is within +/- 0.01 degree of equator...
  if (abs(latitude) < 0.01_deg) {
    // compute upslope force (YSP equation 12)
    f_upslope = (qslam / c0)*((2.0 + xi * xi) / 3.0 - xi)*xi;
    // no cross-slope forcing if coriolis is neglected
    f_xslope = 0.0;
  }
  //=================================================================================================//
  //  else compute up-slope and cross-slope forces, including coriolis, from YSP eq 13, 14           //
  //-------------------------------------------------------------------------------------------------//
  else {
    // coriolis parameter
    greal f = 2.0 * omega * sin(toRadians(latitude));

    // coefficients for YSP equations
    greal wb = 2.0 * qslam / (1000.0 * depth * abs(f));
    greal b = sqrt(1000.0 * abs(f) * depth / (2.0 * c0));
    greal p = exp(-4.0*b) - 2.0 * exp(-2.0*b) * cos(2.0 * b) + 1.0;
    greal b1 = exp(-b * xi);
    greal b2 = exp(-b * (2.0 + xi));
    greal b3 = exp(-b * (2.0 - xi));
    greal b4 = exp(-b * (4.0 - xi));

    // compute upslope force
    f_upslope = wb * ((b1 - b4) * sin(b * xi) + (b2 - b3) * sin(b * (2.0 - xi))) / p;
    // compute cross-slope force
    f_xslope = (f / abs(f)) * wb
      * (xi - 1.0 + ((b1 + b4) * cos(b * xi) - (b2 + b3) * cos(b * (2.0 - xi))) / p);
  }                                                                                      // end else (not near equator)

  //=================================================================================================//
  // conmpute U and V slope winds                                                                    //
  //-------------------------------------------------------------------------------------------------//
  // limit slope angles to ~3 degrees (ewSlope = tan(slopeAngle))
  greal ewSlopeCapped = clamp(ewSlope, -0.0525, 0.0525);
  greal nsSlopeCapped = clamp(nsSlope, -0.0525, 0.0525);

  // zonal slope wind
  ewSlopeWinds = f_upslope * ewSlopeCapped + f_xslope * nsSlopeCapped;
  // meridional slope wind
  nsSlopeWinds = f_upslope * nsSlopeCapped - f_xslope * ewSlopeCapped;

  //=================================================================================================//
  //  add assumed time of day dependence to U and V slope winds and compute vertical slope wind      //
  //-------------------------------------------------------------------------------------------------//
  // add time of day dependence to zonal slope wind
  ewSlopeWinds *= cos_fac;
  // add time of day dependence to meridional slope wind
  nsSlopeWinds *= sin(toRadians(15.0_deg * (time - 11.0)));

  // vetical factor, constant up to half BL height, then drops off as cosine
  greal z_fac = 1.0;
  if (xi > 0.5) {
    // adjust z_fac if over half way up in BL height
    z_fac = 0.5 * (1.0 + cos(TWO_PI *(xi - 0.5)));
  }
  // compute vertical wind component
  verticalSlopeWinds = z_fac * (ewSlope * (ewWind + ewSlopeWinds) + nsSlope * (nsWind + nsSlopeWinds));
}

//! \brief Compute specific heat at constant pressure.
//!
//!  This routine computes specific heat at constant pressure as a function of input temperature
//!  the source for this parameterization is unknown
//! \param t The current temperature (K).
//! \returns The specific heat capacity.
greal SlopeWindsModel::cp(greal t)
{
  return 639.5 + t * (0.123687 + t * 0.00200225);
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//


} // namespace