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
#include "MarsGCMBase.h"
#include "Interpolator.h"
#include "SlopeWindsModel.h"

using namespace std;

namespace GRAM {


//! \copydoc Atmosphere::Atmosphere()
MarsGCMBase::MarsGCMBase(MarsInterpolatorBase& sfc, MarsInterpolatorBase& low, MarsInterpolatorBase& high, MarsDustModelBase& dust)
: Atmosphere(), MarsCommon(this), marsSurface(sfc), marsLower(low), marsUpper(high), marsDustModel(dust)
{
  atmos.setPlanetSpecificMetrics(marsAtmos);
  lower.setPlanetSpecificMetrics(lowerAtmos);
  upper.setPlanetSpecificMetrics(upperAtmos);
}

//! \fn  MarsGCMBase::MarsGCMBase(const MarsGCMBase& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  MarsGCMBase::~MarsGCMBase()
//! \copydoc Atmosphere::~Atmosphere()

//! \fn  MarsGCMBase::setInputParameters(const MarsInputParameters& params)
//! \copydoc PerturbedAtmosphere::setInputParameters()

//! \fn  MarsGCMBase::updateDustModel()
//! \brief Computes a dust optical depth and the dust offset.
//! 
//! This pure virtual function must be implemented in the sub-class.  The implemention should
//! compute a dust optical depth and the dust offset.  It should set parameters coming from the dust model as appropriate.

//! \fn  MarsGCMBase::getGlobalMeanOffset()
//! \brief Computes and returns the global mean height offset.
//!
//! This pure virtual function must be implemented in the sub-class.  The implemention should
//! compute and return the global mean height offset needed for the MARS_GLOBAL_MEAN offset model.
//! \returns The global mean height offset.

//! \fn  MarsGCMBase::getSeasonalOffset()
//! \brief Computes and returns a seasonal height offset.
//!
//! This pure virtual function must be implemented in the sub-class.  The implemention should
//! compute and return the seasonal height offset needed for the MARS_SEASONAL offset model.
//! \returns The seasonal height offset.

//! \fn  MarsGCMBase::getMaxHeight()
//! \brief Returns the maximum height \units{km} in the upper atmosphere model.

//! \fn  MarsGCMBase::getSurfaceFairingHeight()
//! \brief Returns the lower bound \units{km} of the surface fairing region.

//! \copydoc Atmosphere::update()
void MarsGCMBase::update()
{
  // If height is too negative, use user input height.
  if (height < -8.7_km) {
    height = surfaceHeight + userHeightAboveSurface;
  }

  updateHeightOffsets();
  updateFairingHeights();
  updateUpperLower();
  updateScaleHeights();
  updateGasConstant();
  updateAtmos();
}

//! \brief Computes height offsets based on the offset model.
//!
//! This method computes the #heightOffset for the appropriate #offsetModel. Except for the 
//! MARS_GLOBAL_MEAN option, a #currentOffset is computed for the current lat, lon, and time.
//! A thermosphere base height is also computed.
//!
//! \b Inputs
//! \arg #offsetModel
//! \arg #position
//! \arg #ephem
//! \arg #constantHeightOffset
//! \arg #dustOffset
//!
//! \retval #currentOffset
//! \retval #heightOffset
//! \retval #thermosphereBaseHeight
void MarsGCMBase::updateHeightOffsets()
{
  thermosphereBaseHeight = 999.0;
  // Compute a local height offset ...
  if (offsetModel != MARS_GLOBAL_MEAN) {
    // Get densities and scale heights at the first layer of the upper atmos (80 to 85).
    // Use the current position and time, but override the height.
    Position pos = position;
    pos.height = 80.0;

    // Get the lower boundary of the layer (at 80 km).
    marsUpper.setPosition(pos);
    marsUpper.setEphemerisState(ephem);
    marsUpper.updateIndices(0);
    marsUpper.update();
    const AtmosphereState low = marsUpper.getAtmosphereState();
    const MarsAtmosphereState& lowMars = low.getPlanetSpecificMetrics<MarsAtmosphereState>();

    // Get the upper boundary of the layer (at 85 km).
    marsUpper.setPosition(pos);
    marsUpper.setEphemerisState(ephem);
    marsUpper.updateIndices(1);
    marsUpper.update();
    const AtmosphereState& high = marsUpper.getAtmosphereState();
    const MarsAtmosphereState& highMars = high.getPlanetSpecificMetrics<MarsAtmosphereState>();
    thermosphereBaseHeight = highMars.thermosphereBaseHeight;

    // Compute density scale height between the two MTGCM levels (80 to 85 km)
    greal densityScaleHeightAt80km = getScaleHeightFromDelta(5.0, low.density, high.density);
    greal densityScaleHeightDailyAt80km = getScaleHeightFromDelta(5.0, lowMars.densityDaily, highMars.densityDaily);

    // Get densities at upper MGCM height index (80 km)
    marsLower.setPosition(pos);
    marsLower.setEphemerisState(ephem);;
    marsLower.updateIndices(1);
    marsLower.update();
    const AtmosphereState& lowMGCM = marsLower.getAtmosphereState();
    const MarsAtmosphereState& lowMGCMMars = lowMGCM.getPlanetSpecificMetrics<MarsAtmosphereState>();

    // Compute height offset at current time and position
    if (offsetModel == MARS_CURRENT) {
      // current offset = local offset at current time
      currentOffset = densityScaleHeightAt80km * log(lowMGCM.density / low.density);
    }
    else { // offsetModel = MARS_CONSTANT, MARS_SEASONAL, and MARS_DAILY_AVERAGE
      // Use daily averages for offset computations
      currentOffset = densityScaleHeightDailyAt80km * log(lowMGCMMars.densityDaily / lowMars.densityDaily);
    }
  } 
  else {
    // When offsetModel == MARS_GLOBAL_MEAN
    currentOffset = 0.0;
  }

  // Compute the height offset.
  switch(offsetModel) {
  case MARS_CONSTANT:
    // Use user supplied constant offset.
    heightOffset = constantHeightOffset;
    break;
  case MARS_SEASONAL:
    // Use seasonally adjusted input z offset.
    heightOffset = getSeasonalOffset();
    break;
  case MARS_GLOBAL_MEAN:
    // The global mean offset is defined in each subclass.
    heightOffset = getGlobalMeanOffset();
    break;
  case MARS_DAILY_AVERAGE:
  case MARS_CURRENT:
    // Use the current offset.
    heightOffset = currentOffset;
    break;
  }

  // add dust offset to height offset term
  heightOffset += dustOffset;

  // Set the height offset in the upper interpolator.
  marsUpper.setHeightOffset(heightOffset);
}

//! \brief Establish regions for smoothing data.
//!
//! Smoothing (or fairing) occurs between the surface and lower atmosphere data, and between the
//! lower and upper atmosphere data.  This method establishes the heights that defines these regions.;
//!
//! \b Inputs
//! \arg #position
//! \arg #ephem
//!
//! \retval #upperFairingHeight
//! \retval #lowerFairingHeight
//! \retval #surfaceUpperFairingHeight
//! \retval #surfaceLowerFairingHeight
void MarsGCMBase::updateFairingHeights()
{
  // Update indices on an tesgcm interpolator so that we can get fairing heights.
  marsUpper.setPosition(position);
  marsUpper.setEphemerisState(ephem);
  marsUpper.updateIndices(0);
  upperFairingHeight = marsUpper.getUpperFairingHeight();
  lowerFairingHeight = marsUpper.getLowerFairingHeight();

  // Update indices on a surface interpolator so that we can get fairing heights.
  marsSurface.setPosition(position);
  marsSurface.setEphemerisState(ephem);
  marsSurface.updateIndices(0);
  surfaceUpperFairingHeight = marsSurface.getUpperFairingHeight();
  surfaceLowerFairingHeight = marsSurface.getLowerFairingHeight();
}

//! \brief Establish atmosphere states for the lower and upper height interpolation bounds.
//!
//! For the current position, compute the atmosphere states for the interpolations bounds with respect
//! to height.  The #lower and #upper states correspond to the #lowerHeight and #upperHeight levels.
//! A local height offset is also computed.
//!
//! \b Inputs
//! \arg #position
//! \arg #ephem
//! \arg #heightOffset
//! \arg #upperFairingHeight
//! \arg #lowerFairingHeight
//! \arg #surfaceUpperFairingHeight
//! \arg #surfaceLowerFairingHeight
//!
//! \retval #lower
//! \retval #upper
//! \retval #lowerHeight
//! \retval #upperHeight
//! \retval #thermosphereBaseHeight
//! \retval #localHeightOffset
void MarsGCMBase::updateUpperLower()
{
  //=================================================================================================//
  // Begin computing atmosphere states using methodologies that vary by height domain
  // Domains are defined by:
  //  First domain.
  //  --- Upper fairing height. ---
  //  Second domain.
  //  --- Lower fairing height. ---
  //  Third domain.
  //  --- Upper surface fairing height. ---
  //  Fourth domain.
  //  --- Lower surface fairing height. ---
  //  Fifth domain.
  //
  // First domain, upper and lower bounds are both in thermosphere:
  //-------------------------------------------------------------------------------------------------//
  if (height >= upperFairingHeight) {
    // Use the upper interpolator (MTGCM) model for the lower bound. 
    marsUpper.setPosition(position);
    marsUpper.setEphemerisState(ephem);
    marsUpper.updateIndices(0); // use bottom of data interval
    marsUpper.update();
    // Get the lower atmosphere state and the height level
    lower = marsUpper.getAtmosphereState();
    lowerHeight = marsUpper.getPosition().height + heightOffset;

    // Use the upper interpolator (MTGCM) model for the upper bound. 
    marsUpper.setPosition(position);
    marsUpper.setEphemerisState(ephem);
    marsUpper.updateIndices(1); // use top of data interval
    marsUpper.update();
    // Get the upper atmosphere state and the height level
    upper = marsUpper.getAtmosphereState();
    upperHeight = marsUpper.getPosition().height + heightOffset;

    // Adjust thermosphere base height with the height offset.
    thermosphereBaseHeight = upperAtmos.thermosphereBaseHeight + heightOffset;

    // Use the height offset as the local offset.
    localHeightOffset = heightOffset;
  }
  //=================================================================================================//
  // Second domain, lower level in MGCM region, upper in MTGCM region                                //
  //-------------------------------------------------------------------------------------------------//
  else if (height >= lowerFairingHeight) {
    // Use the lower interpolator (MGCM) model for the lower bound. 
    marsLower.setPosition(position);
    marsLower.setEphemerisState(ephem);
    marsLower.updateIndices(0); // use bottom of (uppermost) data interval
    marsLower.update();
    // Get the lower atmosphere state and the height level
    lower = marsLower.getAtmosphereState();
    lowerHeight = marsLower.getPosition().height; // second highest mgcm data level

    // Use the upper interpolator (MTGCM) model for the upper bound. 
    marsUpper.setPosition(position);
    marsUpper.setEphemerisState(ephem);
    marsUpper.updateIndices(0); // use bottom of data interval
    marsUpper.update();
    // Get the upper atmosphere state and the height level
    upper = marsUpper.getAtmosphereState();
    upperHeight = marsUpper.getPosition().height + heightOffset; // adjusted for offsets

    // Need to scale MGCM values based on offsets.
    if (offsetModel == MARS_CONSTANT || offsetModel == MARS_SEASONAL || dustOffset > 0.0) {
      greal offsetLowerBound = 60.0; // why is this 60???
      greal mgcmOffset = 0.0; // declare mgcm offset for interface with tesgcm data
      if (offsetModel == MARS_CONSTANT || offsetModel == MARS_SEASONAL) {
        mgcmOffset = heightOffset - currentOffset;
      }
      else {
        mgcmOffset = dustOffset;
      }
      greal mgcmScaledOffset = mgcmOffset * (lowerHeight - offsetLowerBound) / (upperHeight - offsetLowerBound);
      localHeightOffset = mgcmOffset * (height - offsetLowerBound) / (upperHeight - offsetLowerBound);

      // Get density value for scaling.
      marsLower.updateIndices(1); // use top of (uppermost) data interval
      marsLower.update();
      greal upperDensity = marsLower.getAtmosphereState().density;

      // Create a scaling factor for the lower boundary values
      greal offsetScalingFactor = exp(mgcmScaledOffset * log(lower.density / upperDensity) / 5.0);

      // adjust values at lower boundary
      lower.pressure *= offsetScalingFactor;
      lower.density *= offsetScalingFactor;
      lowerAtmos.pressureDaily *= offsetScalingFactor;
      lowerAtmos.densityDaily *= offsetScalingFactor;
      lowerAtmos.densityMax *= offsetScalingFactor;
      lowerAtmos.densityMin *= offsetScalingFactor;
    } // end if offsetModel == MARS_CONSTANT || offsetModel == MARS_SEASONAL || dustOffset > 0.0
  } 
  //=================================================================================================//
  // Third domain, height is fully above boundary layer, in MGCM region                              //
  //-------------------------------------------------------------------------------------------------//
  else if (height >= surfaceUpperFairingHeight) {
    // Use the lower interpolator (MGCM) model for the lower bound. 
    marsLower.setPosition(position);
    marsLower.setEphemerisState(ephem);
    marsLower.updateIndices(0); // use bottom of data interval
    marsLower.update();
    // Get the lower atmosphere state and the height level
    lower = marsLower.getAtmosphereState();
    lowerHeight = marsLower.getPosition().height;

    // Use the lower interpolator (MGCM) model for the upper bound. 
    marsLower.setPosition(position);
    marsLower.setEphemerisState(ephem);
    marsLower.updateIndices(1); // use top of data interval
    marsLower.update();
    // Get the upper atmosphere state and the height level
    upper = marsLower.getAtmosphereState();
    upperHeight = marsLower.getPosition().height;

    // Need to scale MGCM values based on offsets.
    if (offsetModel == MARS_CONSTANT || offsetModel == MARS_SEASONAL || dustOffset > 0.0) {
      greal offsetLowerBound = 60.0; // why 60 ???
      greal offsetUpperBound = upperFairingHeight; //  why this height???

      greal mgcmOffset = 0.0;
      if (offsetModel == MARS_CONSTANT || offsetModel == MARS_SEASONAL) {
        mgcmOffset = heightOffset - currentOffset;
      }
      else {
        mgcmOffset = dustOffset;
      }

      // If we are above the offeset boundary, then compute the local height offset 
      if (height > offsetLowerBound) {
        localHeightOffset = mgcmOffset * (height - offsetLowerBound) / (offsetUpperBound - offsetLowerBound);
      }

      // Compute a scaling factor for the lower boundary (MGCM) values
      greal lowerScalingFactor = 1.0;
      if (lowerHeight > offsetLowerBound) {
        // get local density scale hieght
        greal hden = getScaleHeightFromDelta(5.0, lower.density, upper.density);
        greal mgcmLowerScaledOffset = mgcmOffset * (lowerHeight - offsetLowerBound) / (offsetUpperBound - offsetLowerBound);
        lowerScalingFactor = exp(mgcmLowerScaledOffset / hden);
      }

      // Compute a scaling factor for the upper boundary (MGCM) values
      greal upperScalingFactor = 1.0;
      if (upperHeight > offsetLowerBound) {
        // get local density scale height
        greal hden = getScaleHeightFromDelta(5.0, lower.density, upper.density);
        greal mgcmUpperScaledOffset = mgcmOffset * (upperHeight - offsetLowerBound) / (offsetUpperBound - offsetLowerBound); // compute ofsz2
        upperScalingFactor = exp(mgcmUpperScaledOffset / hden);
      }

      // Adjust values at lower boundary
      lower.pressure *= lowerScalingFactor;
      lower.density *= lowerScalingFactor;
      lowerAtmos.pressureDaily *= lowerScalingFactor;
      lowerAtmos.densityDaily *= lowerScalingFactor;
      lowerAtmos.densityMax *= lowerScalingFactor;
      lowerAtmos.densityMin *= lowerScalingFactor;

      // Adjust values at upper boundary
      upper.pressure *= upperScalingFactor;
      upper.density *= upperScalingFactor;
      upperAtmos.pressureDaily *= upperScalingFactor;
      upperAtmos.densityDaily *= upperScalingFactor;
      upperAtmos.densityMax *= upperScalingFactor;
      upperAtmos.densityMin *= upperScalingFactor;
    } // end if offsetModel == MARS_CONSTANT || offsetModel == MARS_SEASONAL || dustOffset > 0.0
  } 
  //=================================================================================================//
  // Fourth domain,  lower level in boundary layer, upper level in MGCM region                       //
  //-------------------------------------------------------------------------------------------------//
  else if (height > surfaceLowerFairingHeight) {
    // Use the surface interpolator model for the lower bound. 
    marsSurface.setPosition(position);
    marsSurface.setEphemerisState(ephem);
    marsSurface.updateIndices(0); // use bottom of data interval
    marsSurface.update();
    // Get the lower atmosphere state and the height level
    lower = marsSurface.getAtmosphereState();
    pressureScaleHeightDaily = marsSurface.getPressureScaleHeightDaily();
    lowerHeight = marsSurface.getPosition().height;

    // Use the lower interpolator (MGCM) model for the upper bound. 
    marsLower.setPosition(position);
    marsLower.setEphemerisState(ephem);
    marsLower.updateIndices(0); // use bottom of data interval
    marsLower.update();
    // Get the upper atmosphere state and the height level
    upper = marsLower.getAtmosphereState();
    upperHeight = marsLower.getPosition().height;

    // No local height offset
    localHeightOffset = 0.0;
  }
  //=================================================================================================//
  // Fifth domain processing, height is in boundary layer                                            //
  //-------------------------------------------------------------------------------------------------//
  else {
    // Use the surface interpolator model for the lower bound. 
    marsSurface.setPosition(position);
    marsSurface.setEphemerisState(ephem);
    marsSurface.updateIndices(0); // use bottom of data interval
    marsSurface.update();
    // Get the lower atmosphere state and the height level
    lower = marsSurface.getAtmosphereState();
    pressureScaleHeightDaily = marsSurface.getPressureScaleHeightDaily();
    lowerHeight = marsSurface.getPosition().height;

    // Use the surface interpolator model for the upper bound. 
    marsSurface.setPosition(position);
    marsSurface.setEphemerisState(ephem);
    marsSurface.updateIndices(1); // use top of data interval
    marsSurface.update();
    // Get the upper atmosphere state and the height level
    upper = marsSurface.getAtmosphereState();
    upperHeight = marsSurface.getPosition().height;

    // No local height offset
    localHeightOffset = 0.0;
  } // end height domain processing

  // Don't know why.
  if (height > 80.0_km) {
    localHeightOffset = heightOffset;
  }
}

//! \brief Computes pressure and density scale heights.
//!
//! Given the #lower and #upper atmospheric states corresponding to the #lowerHeight and #upperHeight levels,
//! this method computes the pressure and density scale heights.
//!
//! \b Inputs
//! \arg #position
//! \arg #ephem
//! \arg #heightOffset
//! \arg #lower
//! \arg #upper
//! \arg #lowerHeight
//! \arg #upperHeight
//! \arg #surfaceUpperFairingHeight
//! \arg #surfaceLowerFairingHeight
//!
//! \retval #pressureScaleHeight
//! \retval #pressureScaleHeightDaily
//! \retval #densityScaleHeight
void MarsGCMBase::updateScaleHeights()
{
  greal deltaHeight = upperHeight - lowerHeight;

  //=================================================================================================//
  // Begin computing scale height values using methodologies that vary by height domain (see above). //
  // First, Second, and Third domains                                                                //
  //-------------------------------------------------------------------------------------------------//
  if (height >= surfaceUpperFairingHeight) {
    // Compute scale heights
    pressureScaleHeight = getScaleHeightFromDelta(deltaHeight, lower.pressure, upper.pressure);
    pressureScaleHeightDaily = getScaleHeightFromDelta(deltaHeight, lowerAtmos.pressureDaily, upperAtmos.pressureDaily);
    densityScaleHeight = getScaleHeightFromDelta(deltaHeight, lower.density, upper.density);
  } 
  //=================================================================================================//
  // Fourth domain,  lower level in boundary layer, upper level in MGCM region                       //
  //-------------------------------------------------------------------------------------------------//
  else if (height > surfaceLowerFairingHeight) {
    // Pressure scale heights come from the surface model
    pressureScaleHeight = lower.pressureScaleHeight;
//    pressureScaleHeightDaily = lower.getPressureScaleHeightDaily();

    // Get density scale height from pressure scale height
    greal dtdz = (upper.temperature - lower.temperature) / deltaHeight;
    // compute layer mean temperature
    greal tbar = (lower.temperature + upper.temperature) / 2.0;
    densityScaleHeight = pressureScaleHeight / (1.0 + (pressureScaleHeight / tbar) * dtdz);
  } 
  //=================================================================================================//
  // Fifth domain processing, height is in boundary layer                                            //
  //-------------------------------------------------------------------------------------------------//
  else {
    // MGCM values at Boundary Layer
    marsLower.setPosition(position);
    marsLower.setEphemerisState(ephem);
    marsLower.updateIndices(0);
    marsLower.update();
    AtmosphereState atBL = marsLower.getAtmosphereState();

    // Pressure scale heights come from the surface model
    pressureScaleHeight = upper.pressureScaleHeight;

    // Get density scale height from pressure scale height
    // Compute temperature rate of change
    greal dtdz = (atBL.temperature - lower.temperature) / (marsLower.getPosition().height - lowerHeight);
    // no temperature change below the surface
    if (height <= surfaceHeight) {
      dtdz = 0.0;
    }
    // get layer mean temperature
    greal tbar = (lower.temperature + upper.temperature) / 2.0;
    densityScaleHeight = upper.pressureScaleHeight / (1.0 + dtdz * (upper.pressureScaleHeight / tbar));
  } // end height domain processing

  //=================================================================================================//
  // Cap values if scale heights are out of range  (-9.99 to 99.99)                                  //
  // Specific bogus values are chosen to be obvious to the observer.                                 //
  //-------------------------------------------------------------------------------------------------//
  pressureScaleHeight = clamp(pressureScaleHeight, -9.99, 99.99);
  densityScaleHeight = clamp(densityScaleHeight, -9.99, 99.99);
}

//! \brief Interpolates specifc gas constant over height.
//!
//! Given the #lower and #upper atmospheric states corresponding to the #lowerHeight and #upperHeight levels,
//! this method computes the specifc gas constant and the daily average specifc gas constant.
//!
//! \b Inputs
//! \arg #position
//! \arg #lower
//! \arg #upper
//! \arg #lowerHeight
//! \arg #upperHeight
//!
//! \retval #specificGasConstant
//! \retval #specificGasConstantDaily
void MarsGCMBase::updateGasConstant()
{
  // If below the surface height, use the upper level values.
  if (height <= surfaceHeight) {
    // compute gas constants at upper boundary height (gas law)
    specificGasConstant = upper.pressure / (upper.density * upper.temperature);
    specificGasConstantDaily = upperAtmos.pressureDaily / (upperAtmos.densityDaily * upperAtmos.temperatureDaily);
  }
  else {
    // compute gas constants at lower boundary height (gas law)
    greal lowerR = lower.pressure / (lower.density * lower.temperature);
    greal lowerDailyR = lowerAtmos.pressureDaily / (lowerAtmos.densityDaily * lowerAtmos.temperatureDaily);

    // compute gas constants at tesgcm boundary height (gas law)
    greal upperR = upper.pressure / (upper.density * upper.temperature);
    greal upperDailyR = upperAtmos.pressureDaily / (upperAtmos.densityDaily * upperAtmos.temperatureDaily);

    // linearly interpolate gas constants over height
    Interpolator heightInterp;
    heightInterp.makeFraction(lowerHeight, upperHeight, height);
    specificGasConstant = heightInterp.linear(lowerR, upperR);
    specificGasConstantDaily = heightInterp.linear(lowerDailyR, upperDailyR);
  }
}

//! \brief Interpolate atmospheric states (TPD) over height.
//!
//! This method assumes the #lower and #upper atmosphere states have been computed.  The final temperature
//! is a linear interpolation of these two over height.  When the current height is near the boundary layer,
//! wind speeds factor into the temperature computations.  Pressure is computed using the scale height, and 
//! density uses the ideal gas law.
//!
//! \b Inputs
//! \arg #position
//! \arg #lower
//! \arg #upper
//! \arg #lowerHeight
//! \arg #upperHeight
//!
//! \retval #atmos
void MarsGCMBase::updateAtmos()
{
  // Default for polar ice is N/A (99).
  iceIsPresent = 99;
  // Initialize the ground temp to a bogus value.
  groundTemperature = 999.9;
  // The lowest boundary layer is at 5 m.
  greal boundaryLayer = marsSurface.getLowestBoundaryLayer();

  //=================================================================================================//
  // Begin computing final values using methodologies that vary by height                            //
  // Height is less than 5 m above topo height                                                       //
  //-------------------------------------------------------------------------------------------------//
  if (height < surfaceHeight + boundaryLayer) { 
    // set ground temperature equal to current temperature
    groundTemperature = lower.temperature; 

    // if surface temperature is near CO2 sublimation temperature
    // then set polar ice indicator on and reset surface roughness term
    greal surfaceRoughness;
    if (lower.temperature <= getCO2SublimationTemperature(upper.pressure) + 5.0) {
      // CO2 ice is present
      iceIsPresent = 1; 
      // surface roughness very low
      surfaceRoughness = icySurfaceRoughness; 
    }
    // if surface temperature is well above CO2 sublimation temperature
    else {
      // no CO2 ice
      iceIsPresent = 0; 
      // surface roughness set to 1 cm
      surfaceRoughness = minSurfaceRoughness; 
    }

    // compute evaluation height for winds (assuming no-slip below roughness length)
    greal windHeight = max(surfaceRoughness, height - surfaceHeight);
    // Boundary layer shape factor between surface and middle BL layer height (5 m)
    greal shapeFactor = log(windHeight / surfaceRoughness) / log(boundaryLayer / surfaceRoughness);

    // Adjust wind velocities by multiplying tesgcm boundary value by shape factor
    ewWind = upper.ewWind * shapeFactor;
    nsWind = upper.nsWind * shapeFactor;
    ewWindDaily = upperAtmos.ewWindDaily * shapeFactor;
    nsWindDaily = upperAtmos.nsWindDaily * shapeFactor;

    // convert heights to meters for call to bltp
    greal boundaryLayerMeters = boundaryLayer * 1000.0;
    greal windHeightMeters = windHeight * 1000.0;

    // Compute gravity for bltp routine.
    gravity = getGravity(latitude, latitudeRadius, height);

    // apply BL temperature profile to get current temperature
    temperature = bltp(lower.temperature, boundaryLayerMeters, upper.temperature, upper.ewWind, upper.nsWind, windHeightMeters, shapeFactor, gravity);
    temperatureDaily = bltp(lowerAtmos.temperatureDaily, boundaryLayerMeters, upperAtmos.temperatureDaily, upperAtmos.ewWindDaily, upperAtmos.nsWindDaily, windHeightMeters, shapeFactor, gravity);
    temperatureMax = bltp(lowerAtmos.temperatureMax, boundaryLayerMeters, upperAtmos.temperatureMax, upperAtmos.ewWindDaily, upperAtmos.nsWindDaily, windHeightMeters, shapeFactor, gravity);
    temperatureMin = bltp(lowerAtmos.temperatureMin, boundaryLayerMeters, upperAtmos.temperatureMin, upperAtmos.ewWindDaily, upperAtmos.nsWindDaily, windHeightMeters, shapeFactor, gravity);

    // compute current P by applying scale height to tesgcm bounding height current P
    pressure = upper.pressure * exp((upperHeight - height) / pressureScaleHeight);
    pressureDaily = upperAtmos.pressureDaily * exp((upperHeight - height) / pressureScaleHeightDaily);

    // compute current density from gas law using current P, T, and R
    density = pressure / (specificGasConstant * temperature);
    densityDaily = pressureDaily / (specificGasConstantDaily * temperatureDaily);

    // Min and Max density uses interpolated (using shapeFactor) ratios of min (or max) density to daily average
    Interpolator shapeInterp(shapeFactor);
    densityMax = densityDaily * shapeInterp.linear(lowerAtmos.densityMax / lowerAtmos.densityDaily, upperAtmos.densityMax / upperAtmos.densityDaily);
    densityMin = densityDaily * shapeInterp.linear(lowerAtmos.densityMin / lowerAtmos.densityDaily, upperAtmos.densityMin / upperAtmos.densityDaily);
  }
  //=================================================================================================//
  // Height is greater than 5 m above topo height                                                    //
  //-------------------------------------------------------------------------------------------------//
  else {
    Interpolator heightInterp;
    heightInterp.makeFraction(lowerHeight, upperHeight, height);

    // linearly interpolate winds
    ewWind = heightInterp.linear(lower.ewWind, upper.ewWind);
    nsWind = heightInterp.linear(lower.nsWind, upper.nsWind);
    ewWindDaily = heightInterp.linear(lowerAtmos.ewWindDaily, upperAtmos.ewWindDaily);
    nsWindDaily = heightInterp.linear(lowerAtmos.nsWindDaily, upperAtmos.nsWindDaily);

    // linearly interpolate temperature
    temperature = heightInterp.linear(lower.temperature, upper.temperature);
    temperatureDaily = heightInterp.linear(lowerAtmos.temperatureDaily, upperAtmos.temperatureDaily);
    temperatureMax = heightInterp.linear(lowerAtmos.temperatureMax, upperAtmos.temperatureMax);
    temperatureMin = heightInterp.linear(lowerAtmos.temperatureMin, upperAtmos.temperatureMin);

    // compute pressure by applying scale height
    pressure = upper.pressure * exp((upperHeight - height) / pressureScaleHeight);
    pressureDaily = upperAtmos.pressureDaily * exp((upperHeight - height) / pressureScaleHeightDaily);

    // comoute density from gas law using current T, P, and R
    density = pressure / (specificGasConstant * temperature);
    densityDaily = pressureDaily / (specificGasConstantDaily * temperatureDaily);

    // Min and Max density uses interpolated (using height) ratios of min (or max) density to daily average
    densityMax = densityDaily * heightInterp.linear(lowerAtmos.densityMax / lowerAtmos.densityDaily, upperAtmos.densityMax / upperAtmos.densityDaily);
    densityMin = densityDaily * heightInterp.linear(lowerAtmos.densityMin / lowerAtmos.densityDaily, upperAtmos.densityMin / upperAtmos.densityDaily);
  }

  //=================================================================================================//
  //  if height is < 500 m above local topography and latitude is near pole,                         //
  //  then interpolate to daily average at the pole                                                  //
  //-------------------------------------------------------------------------------------------------//
  if (height < surfaceHeight + 0.5_km && abs(latitude) >= 85.0_deg) {
    // pole interpolation factor (0 at 85, 1 at 90)
    Interpolator poleInterp((abs(latitude) - 85.0_deg) / 5.0);

    // interpolate between average and current pressure
    pressure = poleInterp.linear(pressure, pressureDaily);

    // interpolate between average and current density
    density = poleInterp.linear(density, densityDaily);
    densityMax = poleInterp.linear(densityMax, densityDaily);
    densityMin = poleInterp.linear(densityMin, densityDaily);
  }
}

//! \brief Computes Mars boundary layer temperature.
//!
//! Computes Mars boundary layer temperature from methods used in the NASA Ames Mars General
//! Circulation Model(MGCM), as described by Haberle et al., Jour.Geophys.Res. 104(E4),     
//! 8957 - 8974, 1999, (referred to as H99 below).                                          
//!
//! \param tg        surface temperature \units{K}                                
//! \param z5        height (wrt MOLA aeroid) of 5 m above topography level \units{km}
//! \param t5        temperature at height z5 \units{K}
//! \param u5        zonal wind component at z5 height \units{m/s}
//! \param v5        meridional wind component at z5 height \units{m/s}
//! \param zeval     height to evaluate temperature \units{km}
//! \param factor    asjustment factor    
//! \param grav      gravity \units{m/s^2}
//!
//! \returns The boundary layer temperature \units{K}.
greal MarsGCMBase::bltp(greal tg, greal z5, greal t5, greal u5, greal v5,	greal zeval, greal factor, greal grav) 
{
  // Constants
  const greal cp = this->cp(t5);                           // Specific heat for t5
  const greal theta5 = t5 + grav * z5 / cp;                // potential temperature at z5 height
  const greal udenom = max(0.1, pow(u5, 2) + pow(v5, 2));  // wind speed denominator for computing ri (at least 0.1)

  // Initial values
  greal sqrtfh = 1.0;          // square root of stability function Fh (section 4 H99)
  greal newTemperature = 0.0;  // computed return value temperature
  greal oldTemperature = t5;   // initialize comparison ('previous value') temperature to t5

  // Iterative temperature calculation (convergence limit = 10 iterations).
	for (int i = 0; i < 10; i++) {
    // compute new temperature
		newTemperature = tg + (theta5 - tg) * (1.0 + sqrtfh * factor) / (1.0 + sqrtfh) - grav * zeval / cp;

    // Convergence test: if delta temperature < 0.01 K
    if (abs(newTemperature - oldTemperature) < 0.01) {
      break;                                                          
    }

    // Save the new temperature.
		oldTemperature = newTemperature;  

    // compute Richardson number
    greal ri = (grav * sqrtfh / theta5) * ((theta5 - tg) / (1.0 + sqrtfh)) * z5 / udenom;

    // compute new stability function term based on value of ri (page 8959)
    if (ri < 0.0) {
      sqrtfh = pow(1.0 - 16.0 * ri, 0.25);                                         
    }
    else {
      sqrtfh = pow((1.0 + 15.0 * ri / sqrt(1.0 + 5.0 * ri)), -0.5);
    }
	} // End loop

  // Whether or not it converges, return the new temperature.
	return newTemperature; 
} 

//! \brief Computes CO2 sublimation temperature.
//!
//! Computes CO2 sublimation temperature as a funciton of pressure.
//! Taken from Kieffer and Jakosky, "Mars", 1992, U. of Ariz. Press, p.459.
//!
//! \param p Pressure \units{Pa}.
//! \returns CO2 sublimation temperature \units{K}.
greal MarsGCMBase::getCO2SublimationTemperature(greal p)
{
  return 3182.48 / (23.3494 - log(p / 100.0));
}

//! \brief Computes specific heat at constant pressure.
//!
//! This routine computes specific heat at constant pressure as a function of input temperature.
//! The source for this parameterization is unknown.
//!
//! \param t Temperature \units{K}.
//! \returns Specific heat at constant pressure.
greal MarsGCMBase::cp(greal t)
{ 
  return 639.5 + t * (0.123687 + t * 0.00200225); 
} 

} // namespace