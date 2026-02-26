//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//
// Original models developed by Dr. C. G. (Jere) Justus.
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>
#include "PerturbedAtmosphere.h"
#include "SpiceLoader.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
PerturbedAtmosphere::PerturbedAtmosphere()
{
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
PerturbedAtmosphere::PerturbedAtmosphere(const PerturbedAtmosphere& orig)
: Atmosphere(orig), AuxiliaryAdapter(orig)
{
  time = orig.time;
  perturbationAction = orig.perturbationAction;
  densityPerturbationScale = orig.densityPerturbationScale;
  ewWindPerturbationScale = orig.ewWindPerturbationScale;
  nsWindPerturbationScale = orig.nsWindPerturbationScale;
  verticalWindPerturbationScale = orig.verticalWindPerturbationScale;
  minRelativeStepSize = orig.minRelativeStepSize;
  previousPosition = orig.previousPosition;
  delta = orig.delta;
  pertState = orig.pertState;
  perturbationStep = orig.perturbationStep;
  RHOd = orig.RHOd;
  RHOu = orig.RHOu;
  RHOv = orig.RHOv;
  RHOw = orig.RHOw;
  randomNumberGenerator = orig.randomNumberGenerator;
  userPertState = orig.userPertState;
  hasVerticalWinds = orig.hasVerticalWinds;
  userSuppliedEphemeris = orig.userSuppliedEphemeris;
}

//! \fn  PerturbedAtmosphere::~PerturbedAtmosphere()
//! \copydoc Atmosphere::~Atmosphere()

//! \fn void PerturbedAtmosphere::setEphemerisFastModeOn(bool flag)
//! \copydoc Ephemeris::setFastModeOn()

//! \fn void PerturbedAtmosphere::setSubsolarUpdateTime(greal utime)
//! \copydoc Ephemeris::setSubsolarUpdateTime()

//! \brief Set the applicable input parameters.
//!
//! The routine copies the applicable members of the InputParameters
//! object into the appropriate members of this object.
//! \param params The input parameters.
void PerturbedAtmosphere::setInputParameters(const InputParameters& params)
{
  // First, initialize the Spice data paths (before time or ephemeris settings).
  SpiceLoader spiceLoader;
  spiceLoader.setInputParameters(params);

  // Set the start time.
  time.setStartTime(params.year, params.month, params.day,
    params.hour, params.minute, params.seconds, params.timeScale, params.timeFrame);

  // Notify ephemeris who we are.
  ephemeris.setBody(gramBody);

  // Set perturbation parameters.
  setMinRelativeStepSize(params.minRelativeStepSize);
  setDensityPerturbationScale(params.densityPerturbationScale);
  setEWWindPerturbationScale(params.ewWindPerturbationScale);
  setNSWindPerturbationScale(params.nsWindPerturbationScale);
  if (hasVerticalWinds) {
    setVerticalWindPerturbationScale(params.verticalWindPerturbationScale);
  }
}

//! \brief Sets the folder location of the NAIF SPICE data.
//!
//! This method initializes the SpiceLoader class with the folder location
//! of the NAIF SPICE data. Call this method before invoking the GramTime 
//! or Ephemeris objects.
//! \param path  The location of the NAIF SPICE folder.
void PerturbedAtmosphere::setSpiceDataPath(const std::string& path)
{
  SpiceLoader spiceLoader;
  spiceLoader.setSpiceDataPath(path);
}

//! \brief Initializes the perturbation state.
//!
//! This function seeds the random number generator and resets
//! the perturbation model to the initial state.
//! \param seed An integer betwee 0 and 30,000.
void PerturbedAtmosphere::setSeed(int seed)
{
  //this->seed = seed % 30000;
  //ix = this->seed;
  //iy = 172 * (ix % 176) - 35 * (ix / 176);
  //iz = 170 * (ix % 178) - 63 * (ix / 178);
  //if (iy < 0) iy = iy + 30307;
  //if (iz < 0) iz = iz + 30323;
  randomNumberGenerator.setSeed(seed);
  vector<greal> randomNumbers(3);
  if (hasVerticalWinds) {
    randomNumbers.resize(4);
  }
  randomNumberGenerator.getRandomNumbers(randomNumbers);
  densityRandomNumber = randomNumbers[0];
  ewWindRandomNumber = randomNumbers[1];
  nsWindRandomNumber = randomNumbers[2];
  RHOd = normalPercentagePoint(densityRandomNumber);
  RHOu = normalPercentagePoint(ewWindRandomNumber);
  RHOv = normalPercentagePoint(nsWindRandomNumber);
  if (hasVerticalWinds) {
    verticalWindRandomNumber = randomNumbers[3];
    RHOw = normalPercentagePoint(verticalWindRandomNumber);
  }
  perturbationStep = 0.0;
  previousPosition.height = 0.0;
  previousPosition.latitude = 0.0;
  previousPosition.longitude = 0.0;
  previousPosition.height = 0.0;
  previousPosition.elapsedTime = 0.0;
  delta = previousPosition;
  updateStatus = INITIAL_STATE;
}

//! \brief Set the minimum relative step size for perturbations.
//!
//! \param stepSize Between 0 and 1.
void PerturbedAtmosphere::setMinRelativeStepSize(greal stepSize) 
{ 
  minRelativeStepSize = clamp(stepSize, 0.0, 1.0); 
}

//! \brief Set the perturbation scale factor.
//!
//! \param scaleFactor Between 0 and 2.
void PerturbedAtmosphere::setPerturbationScales(greal scaleFactor)
{
  densityPerturbationScale = clamp(scaleFactor, 0.0, 2.0);
  ewWindPerturbationScale = densityPerturbationScale;
  nsWindPerturbationScale = densityPerturbationScale;
  verticalWindPerturbationScale = densityPerturbationScale;
}

//! \brief Set the perturbation scale factor.
//!
//! \param scaleFactor Between 0 and 2.
void PerturbedAtmosphere::setDensityPerturbationScale(greal scaleFactor)
{
  densityPerturbationScale = clamp(scaleFactor, 0.0, 2.0);
}

//! \brief Set the perturbation scale factor.
//!
//! \param scaleFactor Between 0 and 2.
void PerturbedAtmosphere::setEWWindPerturbationScale(greal scaleFactor)
{
  ewWindPerturbationScale = clamp(scaleFactor, 0.0, 2.0);
}

//! \brief Set the perturbation scale factor.
//!
//! \param scaleFactor Between 0 and 2.
void PerturbedAtmosphere::setNSWindPerturbationScale(greal scaleFactor)
{
  nsWindPerturbationScale = clamp(scaleFactor, 0.0, 2.0);
}

//! \brief Set the perturbation scale factor.
//!
//! \param scaleFactor Between 0 and 2.
void PerturbedAtmosphere::setVerticalWindPerturbationScale(greal scaleFactor)
{
  verticalWindPerturbationScale = clamp(scaleFactor, 0.0, 2.0);
}

//! \copydoc Atmosphere::setPosition()
void PerturbedAtmosphere::setPosition(const Position& pos)
{
  Atmosphere::setPosition(pos);
  if (!position.isPlanetoCentric) {
    planetodeticPosition = position;
    position.convertToPlanetocentric(polarRadius, equatorialRadius);
  }
  updatePosition();

  // Compute deltas based on input frame
  if (updateStatus != INITIAL_STATE) {
    delta.latitude = fabs(pos.latitude - previousPosition.latitude);
    delta.longitude = fabs(pos.longitude - previousPosition.longitude);
    delta.height = (pos.height) - (previousPosition.height);
    delta.elapsedTime = pos.elapsedTime - previousPosition.elapsedTime;

    // Correct for cases near 0/360 longitude
    if (delta.longitude > 180.0) {
      delta.longitude = 360.0 - delta.longitude;
    }
    if (delta.longitude < 0.0) {
      delta.longitude += 360.0;
    }

    // Correct for cases passing over poles
    if (delta.longitude > 90.0 && (fabs(pos.latitude) > 70.0 || fabs(previousPosition.latitude) > 70.0)) {
      delta.longitude = fabs(180.0 - delta.longitude);
      delta.latitude = 180.0 - fabs(pos.latitude) - fabs(previousPosition.latitude);
    }
    savedPosition = pos;
  }
  else {
    savedPosition = previousPosition = pos;
  }
}

//! \brief Set the change in position.
//!
//! Perturbations are correlated with previously computed perturbations
//! with a dependency on the change in position.  A call to setPostion()
//! will compute the deltas. Or the computed deltas can be overridden with
//! a call to this function.
//! \param delta A Position structure containing the change in position.
//! The fields height, latitude, longitude, elapsedTime, and isPlanetoCentric must be set.
void PerturbedAtmosphere::setDelta(const Position& delta)
{
  this->delta.latitude = delta.latitude;
  this->delta.longitude = delta.longitude;
  this->delta.height = delta.height;
  this->delta.elapsedTime = delta.elapsedTime;
}

//! \brief Overrides the computed perturbation state with a user provided state.
//!
//! The perturbation model will internally calculate new random numbers needed for perturbations.
//! However, if this method is called, then the supplied random numbers will be used instead.
//! \param pState A PerturbationState.
void PerturbedAtmosphere::setPerturbationState(PerturbationState& pState)
{ 
  pertState = pState; 
  userPertState = true;
}

void PerturbedAtmosphere::updateEphemeris()
{
  // Get ephemeris values from the user or from the Ephemeris object.
  if (userSuppliedEphemeris) {
    userSuppliedEphemeris = false;
  }
  else {
    // Make sure the ephemeris body is set.
    if (ephemeris.getBody() != gramBody) {
      ephemeris.setBody(gramBody);
    }
    // Time and longitude are needed for local solar time, and longitude of the sun
    ephemeris.setTime(time);
    ephemeris.setLongitude(longitude);
    // Latitutde is only needed for solar zenith angle
    ephemeris.setLatitude(latitude);
    // Compute the ephemeris state.
    ephemeris.update();
    ephem = ephemeris.getEphemerisState();
  }
}

//! \brief The primary interface for perturbation calculations.
//!
//! This method should be called after the atmospheric state has been updated
//! with Atmosphere::update().
void PerturbedAtmosphere::updatePerturbations(greal meanDensity, greal meanEWWind, greal meanNSWind)
{
  if (meanDensity < 0.0) {
    meanDensity = density;
  }
  if (meanEWWind < -9000.0) {
    meanEWWind = ewWind;
  }
  if (meanNSWind < -9000.0) {
    meanNSWind = nsWind;
  }

  // Compute speed of sound, if necessary.
  if (speedOfSound == 0.0) {
    updateSpecificHeatRatio();
    updateSpeedOfSound();
  }

  greal pertHigh, pertLow;
  getPerturbationFactors(pertLow, pertHigh);

  // Set vertical and horizontal scale parameters
  greal verticalScale, horizontalScale;
  getScaleParameters(verticalScale, horizontalScale);

  // Assumptions (checking occurs in debug mode only)
  assert(horizontalScale != 0.0);
  assert(verticalScale != 0.0);

  // Relative displacements between previous and current position
  greal relNS = totalRadius * toRadians(delta.latitude) / horizontalScale;
  greal relEW = totalRadius * cos(toRadians(latitude)) * toRadians(delta.longitude) / horizontalScale;
  greal relHeight = delta.height / verticalScale;
  greal relWinds = sqrt(meanEWWind * meanEWWind + meanNSWind * meanNSWind) * delta.elapsedTime / (1000.0 * horizontalScale);
  // Evaluate correlation and step size relative to accuracy limit
  if (perturbationAction == UPDATE_PERTS) {
    perturbationStep = perturbationStep + fabs(relNS) + fabs(relEW) + fabs(relHeight) + relWinds;
    updateStatus = STEP_UPDATED;
    previousPosition = savedPosition;
  }
  else {
    updateStatus = NO_UPDATES;
  }
  relativeStepSize = -perturbationStep / log(0.995);

  vector<greal> randomNumbers(3);
  if (hasVerticalWinds) {
    randomNumbers.resize(4);
  }

  greal corPrevious; // Weight of the previous Rho (CORREL)
  greal corRandom;   // Weight of the new random number (corbeta)
  // If the relative step size is too small or No Updates 
  if (relativeStepSize <= minRelativeStepSize || perturbationAction == DO_NOT_UPDATE_PERTS) {
    corPrevious = 1.0;
    corRandom = 0.0;
    densityRandomNumber = 0.5;
  }
  else {
    // Get uniform and Gaussian random numbers from PPND
    if (!userPertState) {
      randomNumberGenerator.getRandomNumbers(randomNumbers);
      densityRandomNumber = randomNumbers[0];
    }
    corPrevious = exp(-perturbationStep);
    corRandom = sqrt(1.0 - corPrevious * corPrevious);
    perturbationStep = 0.0;
    updateStatus = PERTS_UPDATED;
  }

  // If the user supplied a PerturbationState, then use their Rhos.
  // Otherwise, use the rhos from the previous update.
  if (!userPertState) {
    densityRho = RHOd;
    ewWindRho = RHOu;
    nsWindRho = RHOv;
    verticalWindRho = RHOw;
  }

  // Compute Density plus ~ 1 sigma
  greal densHigh = density * pertHigh;
  greal dPlus = densHigh - density;
  // Compute Density minus ~ 1 sigma
  greal densLow = density * pertLow;
  greal dMinus = density - densLow;

  // Current random density perturbation value, correlated with
  // previous random density perturbation
  RHOd = corPrevious * densityRho + corRandom * normalPercentagePoint(densityRandomNumber);
  greal densRand;
  if (RHOd < 0.0) {
    densRand = RHOd * dMinus * densityPerturbationScale;
  }
  else {
    densRand = RHOd * dPlus * densityPerturbationScale;
  }

  // Add random density perturbation
  // Save as perturbed density, for output
  // Check upper and lower bounds on density perturbations
  perturbedDensity = clamp(density + densRand, 0.1 * meanDensity, 10.0 * meanDensity);

  // Standard deviation in random density perturbation (% of unperturbed mean) for output
  densityStandardDeviation = densityPerturbationScale * 0.5 * fabs(densHigh - densLow) / meanDensity;
  if (getBody() == MARS) {
    densityStandardDeviation = densityPerturbationScale * fabs(densHigh - meanDensity) / meanDensity;
  }

  // Adjust random DENSHI, DENSLO for rpscale
  highDensity = density + densityPerturbationScale * (densHigh - density);
  lowDensity = density + densityPerturbationScale * (densLow - density);
  if (lowDensity < 0.1 * meanDensity) {
    lowDensity = 0.1 * meanDensity;
  }

  // Convert random density perturbation to % of (unperturbed) mean
  densityPerturbation = (perturbedDensity - density) / meanDensity;

  // Compute Speed of sound with perturbed density
  perturbedSpeedOfSound = sqrt(specificHeatRatio * pressure / perturbedDensity);

  // Standard deviation for wind perturbations
  getWindDeviations(ewStandardDeviation, nsStandardDeviation, verticalStandardDeviation);

  // Compute EW wind perturbations and total wind
  if (corRandom != 0.0 && !userPertState) {
    ewWindRandomNumber = randomNumbers[1];
  }
  RHOu = corPrevious * ewWindRho + corRandom * normalPercentagePoint(ewWindRandomNumber);

  // EW component of perturbation in wind and total wind
  greal ewPert = RHOu * ewStandardDeviation;
  greal ewTotal = ewWind + ewPert;

  // Compute NS wind perturbations
  if (corRandom != 0.0 && !userPertState) {
    nsWindRandomNumber = randomNumbers[2];
  }
  RHOv = corPrevious * nsWindRho + corRandom * normalPercentagePoint(nsWindRandomNumber);

  // NS component of perturbation in wind and total wind
  greal nsPert = RHOv * nsStandardDeviation;
  greal nsTotal = nsWind + nsPert;

  // Compute Vertical wind perturbations
  greal vertPert = 0.0;
  if (hasVerticalWinds) {
    if (corRandom != 0.0 && !userPertState) {
      verticalWindRandomNumber = randomNumbers[3];
    }
    RHOw = corPrevious * verticalWindRho + corRandom * normalPercentagePoint(verticalWindRandomNumber);

    // Vertical component of perturbation in wind and total wind
    vertPert = RHOw * verticalStandardDeviation;
  }
  else {
    RHOw = 0.0;
  }

  // Limit winds to sound speed
  greal ssp2 = speedOfSound * speedOfSound;
  greal fc = ssp2 - ewTotal * ewTotal - nsTotal * nsTotal;
  greal fact = 1.0;
  if (fc < 0.0) {
    // Find multiplier factor for wind perturbations that keep total
    // wind speed <= speed of sound
    fc = ssp2 - ewWind * ewWind - nsWind * nsWind;
    greal fa = ewPert * ewPert + nsPert * nsPert;
    if (fa <= 0.0) {
      fa = 1.0;
    }
    greal fb = ewPert * ewWind + nsPert * nsWind;
    fact = (-fb + sqrt(fabs(fb * fb + fa * fc))) / fa;
    fact = clamp(fact, 0.0, 1.0);
  }

  //   Recompute perturbations and total winds with required factor
  ewWindPerturbation = ewPert * fact;
  nsWindPerturbation = nsPert * fact;
  perturbedEWWind = ewWind + ewWindPerturbation;
  perturbedNSWind = nsWind + nsWindPerturbation;
  if (hasVerticalWinds) {
    verticalWindPerturbation = vertPert * fact;
    perturbedVerticalWind = verticalWind + verticalWindPerturbation;
  }
  userPertState = false;
}

//! \brief Percentage points of the normal distribution.
//! 
//! Algorithm AS 111 Appl. Statist. (1977) Vol. 26, p. 118
//! "Percentage Points of the Normal Distribution"
//! Produces normal deviate corresponding to lower tail area of P.
//! \param P A probability between 0 and 1.
//! \returns The lower tail area corresponding to P.
greal PerturbedAtmosphere::normalPercentagePoint(greal P)
{
  //      Single precision version with error epsilon = 2 ** (-31).
  //      For double precision version, change REAL to DOUBLE PRECISION
  //      in the FUNCTION statement and the declaration of variables;
  //      change E0 to D0 in the DATA statements and change ABS, ALOG
  //      and SQRT to abs, DLOG and DSQRT in the assignment statements.
  //      The hash sums are the sums of the moduli of the coefficients.
  //      They have no inherent meanings, but are included for use in
  //      checking transpositions.
  // 

  greal Q, R;
  // -----------------------------------------------
  // 
  const greal SPLIT = 0.42;
  // 
  const greal A0 = 2.50662823884;
  const greal A1 = -18.61500062529;
  const greal A2 = 41.39119773534;
  const greal A3 = -25.44106049637;
  const greal B1 = -8.47351093090;
  const greal B2 = 23.08336743743;
  const greal B3 = -21.06224101826;
  const greal B4 = 3.13082909833;
  // 
  //      Hash sum for a & b = 143.70383558076
  // 
  const greal C0 = -2.78718931138;
  const greal C1 = -2.29796479134;
  const greal C2 = 4.85014127135;
  const greal C3 = 2.32121276858;
  const greal D1 = 3.54388924762;
  const greal D2 = 1.63706781897;
  greal ppnd = 0.0;
  // 
  //      Hash sum for c & d = 17.43746520924
  // 
  // 
  Q = P - 0.5;
  if (fabs(Q) <= SPLIT) {
    R = Q * Q;
    ppnd = Q * (((A3 * R + A2) * R + A1) * R + A0) / ((((B4 * R + B3) * R + B2) * R + B1) * R + 1.0);
    return ppnd;
  }
  R = P;
  if (Q > 0.0) {
    R = 1.0 - P;
  }
  if (R >= 0.0) {
    R = sqrt((-log(R)));
    ppnd = (((C3 * R + C2) * R + C1) * R + C0) / ((D2 * R + D1) * R + 1.0);
    if (Q < 0.0) {
      ppnd = -ppnd;
    }
    return ppnd;
  }
  cout << " normalPercentagePoint ERROR";
  exit(1);
//  ppnd = 0.0;
//  return ppnd;
}

//! \copydoc Ephemeris::findDates
void PerturbedAtmosphere::findDates(greal targetLongitudeSun, greal targetSolarTime, GramTime gramTime[3], greal lonSun[3], greal tlst[3])
{
  // Make sure the ephemeris body is set.
  if (ephemeris.getBody() != gramBody) {
    ephemeris.setBody(gramBody);
  }
  ephemeris.setLongitude(longitude);
  ephemeris.setTime(time);
  ephemeris.findDates(targetLongitudeSun, targetSolarTime, gramTime, lonSun, tlst);
}

//! \fn void PerturbedAtmosphere::setStartTime(const GramTime& gramTime)
//! \brief Sets the start time (epoch) of the simulation.
//!
//! \param gramTime The start time as a GramTime object.

//! \fn const GramTime& PerturbedAtmosphere::getStartTime()
//! \brief Retieves the start time (epoch) of the simulation.
//!
//! \returns The start time as a GramTime object.

//! \fn void PerturbedAtmosphere::setPerturbationAction(PerturbationAction action)
//! \brief Sets the perturbation action.

//! \fn PerturbedAtmosphere::getInputParameters()
//! \brief Retieves the most recently set input parameters.
//!
//! \returns The most recently set InputParameters.

//! \fn PerturbedAtmosphere::getMinRelativeStepSize()
//! \brief Retieves the minimum relative step size for perturbations.
//!
//! \returns The minimum relative step size for perturbations.

//! \fn PerturbedAtmosphere::getPerturbationScaleFactor()
//! \brief Retieves the scale factor for perturbations.
//!
//! \returns The scale factor for perturbations.

//! \fn PerturbedAtmosphere::getSeed()
//! \brief Retieves the random number seed for perturbations.
//!
//! \returns The random number seed for perturbations.

//! \fn PerturbedAtmosphere::getPerturbationState()
//! \brief Retieves the perturbation state.
//!
//! \returns The perturbation state.

//! \fn PerturbedAtmosphere::getPerturbationFactors(greal& pertLow, greal& pertHigh)
//! \brief Retieves the low and high perturbation factors.
//!
//! This pure virtual function must be overridden by each subclass to provide the low and high
//! perturbation factors
//! \param pertLow The low perturbation factor.
//! \param pertHigh The high perturbation factor.

//! \fn PerturbedAtmosphere::getScaleParameters(greal& verticalScale, greal& horizontalScale)
//! \brief Retieves the scale parameters.
//!
//! This pure virtual function must be overridden by each subclass to provide the vertical and horizontal
//! scale factors
//! \param verticalScale The vertical scale factor.
//! \param horizontalScale The horizontal scale factor.

//! \fn PerturbedAtmosphere::getWindDeviations(greal& ewStdDev, greal& nsStdDev, greal& vertStdDev)
//! \brief Retieves the wind standard deviations.
//!
//! This pure virtual function must be overridden by each subclass to provide the wind deviations.
//! \param ewStdDev The East/West wind standard deviation.
//! \param nsStdDev The North/South wind standard deviation.
//! \param vertStdDev The vertical wind standard deviation.

} // namespace
