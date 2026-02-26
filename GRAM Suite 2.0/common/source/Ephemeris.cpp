//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>
#include "Ephemeris.h"
#include "SpiceUsr.h"
#include "SpiceLoader.h"
#include "error_strings.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
Ephemeris::Ephemeris()
{
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
Ephemeris::Ephemeris(const Ephemeris& orig)
{
  gramBody = orig.gramBody;
  naifId = orig.naifId;
  for (size_t i = 0; i < nameSize; ++i) {
    naifName[i] = orig.naifName[i];
    naifParentName[i] = orig.naifParentName[i];
    frameName[i] = orig.frameName[i];
  }
  longitude = orig.longitude;
  latitude = orig.latitude;
  startTime = orig.startTime;
  userInputs = orig.userInputs;
  state = orig.state;
}

//! \fn Ephemeris::~Ephemeris()
//! \copydoc Atmosphere::~Atmosphere()

//! \fn Ephemeris::setTime(const GramTime& tTime)
//! \brief Set the time for ephemeris computations.

//! \fn Ephemeris::setLongitude(greal longitude)
//! \brief Set the longitude in degrees.

//! \fn Ephemeris::setLatitude(greal latitude)
//! \brief Set the latitude in degrees.

//! \fn void Ephemeris::setFastModeOn(bool flag)
//! \brief Turn on/off the fast ephemeris mode.
//!
//! Ephemeris computations can be computationally expensive.  The fast mode will 
//! minimize ephemeris computations.  Normally, these values are computed with each
//! update cycle.  
//!
//! Some ephemeris values change very little over short time periods since they
//! are dependent on the location of the body relative to the sun or earth. 
//! In fast mode the one way light time, longitude of the sun,
//! and orbital radius are computed once. 
//!
//! Other ephemeris values are dependent on the rotation of the body. 
//! The subsolar latitude and longitude, and local solar time
//! are computed on a periodic basis in fast mode.  See setSubsolarUpdateTime() for 
//! more details.
//!
//! Once in fast mode, an update of the ephemeris values can be forced by
//! subsequent calls to setFastModeOn().  
//!
//! The fidelity of the atmospheric state may be affected by the ephemeris fast mode.
//! Fast mode should not be used if the the highest precision is required.  If computational
//! speed is the greater concern, then the fast mode may help.  But the user should 
//! analyze the effects of fast mode on their simulation.
//!
//! \param flag True turns fast mode on, false returns to normal mode.

//! \fn void Ephemeris::setSubsolarUpdateTime(greal utime) 
//! \brief Set the update cycle for subsolar lat/lon while in fast mode.
//!
//! In the ephemeris fast mode, subsolar latitude and longitude are computed
//! on a periodic basis depending on the elapsed time.  The cycle is controlled
//! by the subsolar update time.  If the elapsed time has increased more than the
//! subsolar update time since the last computations, then the subsolar values will
//! be updated.
//!
//! \param utime The update time in seconds.

//! \fn Ephemeris::getEphemerisState()
//! \brief Get the updated ephemeris values.

//! \brief Set user computed ephemeris values.
//!
//! If a user wishes to replace ephemeris computations with their own ephemeris
//! values, then this function should be called prior to each call to update().
//! \param eState The user supplied ephemeris state.
void Ephemeris::setEphemerisState(const EphemerisState& eState) {
  userInputs = true;
  state = eState;
}

//! \brief Computes the current ephemeris values.
//!
//!
void Ephemeris::update()
{
  // Use input ephemeris values, if supplied by user
  if (userInputs) {
    // Reset to force computations the next time update is called.
    // This forces the user to override ephemeris in each time step
    // in order to encourage correct ephemeris are supplied.
    userInputs = false;
  }
  else {
    try {
      // If the fast mode is on, then some values will not get computed 
      // on each update.  This is controlled by the computeFlag which is set
      // set below at the end of each update.
      if (computeFlag) {
        // If the start time was entered in Earth receive time, then the
        // one-way light time must be computed before any other ephemeris computations.
        updateOneWayLightTime();
        updateLongitudeOfTheSun();
        updateOrbitalRadius();
        updateSecondsPerSol(); // constant values
        updateSolarDeclination();
      }
      else {
        startTime.setOneWayLightTime(oneWayLightTime * 60.0);
      }

      // These values need to be computed more regularly even in fast mode.
      // Force computation based on time passage.
      if (computeFlag || startTime.getElapsedTime() - subsolarLastTime > subsolarUpdateTime) {
        updateSubsolarLatLon();
        subsolarLastTime = startTime.getElapsedTime();
        updateSolarZenithAngle();
      }
      // This a fast version base on subsolar longitude 
      // This call must appear after the call to updateSubsolarLatLon()
      updateLocalSolarTime();    

      // If the fast mode is on, turn off recomputes
      computeFlag = !fastMode;
    }
    catch (const string& msg) {
      throw(SPICE_ERROR_MESSAGE + msg
        + "\n       Output data may be corrupted at time: " +  to_string(startTime.getElapsedTime())
        + "\n       This is an unrecoverable error.");
    }
  }
}

//! \brief Sets the target body for ephemeris computations.
//!
//! \param body The target body.
void Ephemeris::setBody(GRAM_BODY body)
{
  // Save the body
  gramBody = body;

  // Use the SpiceLoader to load spice data files.
  SpiceLoader spiceLoader;
  spiceLoader.load(gramBody);

  // By default we will null out the parent name (used by moons).
  naifParentName[0] = 0;

  // Set the body name/id and the SPK file.
  // The SPK file defines body ephemeris data.
  switch (gramBody) {
  case VENUS:
    setBodyName("Venus");
    break;
  case EARTH:
    setBodyName("Earth");
    break;
  case MARS:
    setBodyName("Mars");
    break;
  case JUPITER:
    setBodyName("Jupiter");
    break;
  case URANUS:
    setBodyName("Uranus");
    break;
  case NEPTUNE:
    setBodyName("Neptune");
    break;
  case SATURN:
    setBodyName("Saturn");
    break;
  case TITAN:
    setBodyName("Titan");
    // Since Titan is a moon, set the parent body name/id.
    setBodyParentName("Saturn");
    break;
  default:
    break;
  }
}

//! \brief Sets the target body name and Spice id.
//!
//! \param name The target body name.
void Ephemeris::setBodyName(std::string name)
{
  // Save the NAIF (Spice) name.
  name.copy(naifName, nameSize);

  // Declere return values for Spice call.
  // Need to use Spice types to avoid compiler warnings.
  SpiceBoolean found;
  SpiceInt naifIdx;
  // SPICE CALL: Get the NAIF (Spice) id given the body name.
  bodn2c_c(naifName, &naifIdx, &found);
  // Copy SpiceInt(naifIdx) into long int(naifId)
  naifId = naifIdx;
  if (!found || failed_c()) {
    throw(SPICE_ERROR_MESSAGE + "Cannot locate a SPICE id for body: " + naifName
      + "\n       This is an unrecoverable error.");
  }

  if (name == "Earthx") {
    std::string frame = "ITRF93";
    frame.copy(frameName, nameSize);
  }
  else {
    // Fixed body reference frames are all of the form IAU_BODYNAME.
    // See the PCK file pck00010.tpc for more information.
    std::string frame = "IAU_";
    frame.append(name);
    frame.copy(frameName, nameSize);
  }
}

//! \brief Sets the parent body name (for moons).
//!
//! \param name The parent body name.
void Ephemeris::setBodyParentName(std::string name)
{
  name.copy(naifParentName, nameSize);
}

//! \brief Computes one way light time from earth.
//!
//! If the start time was entered as earth receive time, then this
//! function must be called before other ephemeris computations are 
//! made so that the start time can be adjusted to planet event time.
//!
//! \b Inputs
//! \arg #naifId    
//! \arg #startTime 
//!
//! \retval #oneWayLightTime 
void Ephemeris::updateOneWayLightTime()
{
  SpiceDouble  arrive;
  SpiceDouble  howlng;
  SpiceInt     earthBodyId = 399;

  startTime.setOneWayLightTime(0.0);

  // SPICE CALL: Gets the owlt in seconds
  // Args: time at observer,  observer id, signal direction, target id, time at arrival, transmit time
  ltime_c(startTime.getSpiceTime(), naifId, "->", earthBodyId, &arrive, &howlng);
 
  // Test for SPICE failure
  if (failed_c()) {
    throw string("Error in one way light time calculation.");
  }

  // Set the owlt in the time object
  // This is needed if the start time was set as Earth receive time.
  startTime.setOneWayLightTime(howlng);

  // convert to minutes
  oneWayLightTime = howlng / 60.0;
}

//! \brief Computes subsolar latitude and longitude.
//!
//! \b Inputs
//! \arg #naifId   
//! \arg #frameName
//! \arg #startTime
//!
//! \retval #subsolarLatitude 
//! \retval #subsolarLongitude
void Ephemeris::updateSubsolarLatLon()
{
  const size_t SIZE = 256;
  SpiceChar    method[SIZE] = "NEAR POINT/ELLIPSOID";
  SpiceChar    abcorr[SIZE] = "NONE";
  SpiceDouble  spoint[3];
  SpiceDouble  trgepc;
  SpiceDouble  srfvec[3];
  SpiceDouble  spcrad;
  SpiceDouble  spclon;
  SpiceDouble  spclat;

  // SPICE CALL: Compute the subsolar point at the specified time
  // Args: Computation method, target body name, time, body fixed frame, aberration correction, observing body name, 
  //       sub-solar point on target, sub-solar point epoch, vector from observer to target
  subslr_c(method, naifName, startTime.getSpiceTime(), frameName, abcorr, naifName, spoint, &trgepc, srfvec);

  // SPICE CALL: convert rectangular coordinates to lat/lon
  // Args: rectangular coordinates, distance from origin to point, longitude, latitude
  reclat_c(spoint, &spcrad, &spclon, &spclat);

  // Test for SPICE failure
  if (failed_c()) {
    throw string("Error in subsolar lat/lon calculation.");
  }

  // Save the results.  Make sure longitude is in range [0, 360).
  subsolarLatitude = toDegrees(spclat);
  subsolarLongitude = wrapDegrees(toDegrees(spclon));
}

//! \brief Computes subsolar latitude and longitude.
//!
//! \b Inputs
//! \arg #naifId   
//! \arg #frameName
//! \arg #startTime
//!
//! \retval #subsolarLatitude 
//! \retval #subsolarLongitude
void Ephemeris::updateSolarDeclination()
{
  const size_t SIZE = 256;

  SpiceChar    target[SIZE] = "SUN";
  SpiceDouble  et = startTime.getSpiceTime();
  SpiceChar    abcorr[SIZE] = "NONE";
  SpiceChar    frame[SIZE] = "J2000";
  SpiceDouble  ptarg[3];
  SpiceDouble  lt;

  // SPICE CALL: Return the position of a target body relative to an observing body.
  // Args: Target body name, time, body fixed frame, aberration correction, observing body name, 
  //       rectangular coordinates of the target, one way light time between observer and target
  spkpos_c(target, et, frame, abcorr, naifName, ptarg, &lt);

  SpiceDouble      range;
  SpiceDouble      ra;
  SpiceDouble      dec;

  // SPICE CALL: convert rectangular coordinates to ange, right ascension, and declination.
  // Args: rectangular coordinates, distance from origin to point,right ascension, declination
  recrad_c(ptarg, &range, &ra, &dec);
  
  solarDeclination = toDegrees(dec);
  solarRightAscension = toDegrees(ra);

  // Test for SPICE failure
  if (failed_c()) {
    throw string("Error in solar declination calculation.");
  }
}

//! \brief Computes target body longitude of the sun.
//!
//! \b Inputs
//! \arg #naifName
//! \arg #naifParentName 
//! \arg #startTime 
//!
//! \retval #longitudeSun 
greal Ephemeris::updateLongitudeOfTheSun()
{
  SpiceChar abcorr[10] = "NONE";
  SpiceChar* bodyName;
  if (naifParentName[0]) {
    bodyName = naifParentName;
  }
  else {
    bodyName = naifName;
  }

  // SPICE CALL: Get the planetocentric longitude of the sun.
  // Args: target body name, time, aberration correction
  SpiceDouble slon = lspcn_c(bodyName, startTime.getSpiceTime(), abcorr);

  // Test for SPICE failure
  if (failed_c()) {
    throw string("Error in longitude of the sun calculation.");
  }

  // Save the result.  Make sure longitude is in range [0, 360).
  longitudeSun = wrapDegrees(toDegrees(slon));
  return toDegrees(slon);
}

//! \brief Computes target body solar time at specified longitude.
//!
//! \b Inputs
//! \arg #naifId    
//! \arg #startTime 
//! \arg #longitude 
//!
//! \retval #solarTime 
void Ephemeris::updateLocalSolarTime()
{
  //const size_t SIZE = 256;
  //SpiceChar    type[SIZE] = "PLANETOCENTRIC";
  //SpiceInt     timlen = SIZE;
  //SpiceInt     ampmlen = SIZE;
  //SpiceInt     hr;
  //SpiceInt     mn;
  //SpiceInt     sc;
  //SpiceChar    time[SIZE];
  //SpiceChar    ampm[SIZE];
  //// The Spice convention for Planetocentric is East longitude positive
  //SpiceDouble  lon = toRadians(longitude);

  //// SPICE CALL: Compute the local solar time of the target at the specified longitude.
  //// Args: time, target body id, longitude, type of longitude, string length of time,  string length of ampm,
  ////       local hour (24-hour scale), minutes, seconds, 24-hour time string, am/pm time string
  //et2lst_c(startTime.getSpiceTime(), naifId, lon, type, timlen, ampmlen,
  //  &hr, &mn, &sc, time, ampm);

  //// Test for SPICE failure
  //if (failed_c()) {
  //  throw string("Error in local solar time calculation.");
  //}

  //// Save as local solar time in hours (24-hour scale)
  //double solarTimeSpice = hr + mn / 60.0 + sc / 3600.0;

  solarTime = 12.0 + (longitude - subsolarLongitude) / 15.0;
  if (gramBody == URANUS || gramBody == VENUS) {
    solarTime = 24.0 - solarTime;
  }
  if (solarTime > 24.0) {
    solarTime -= 24.0;
  }
  else if (solarTime < 0.0) {
    solarTime += 24.0;
  }
  solarHourAngle = 15.0_deg * (solarTime - 12.0);
}

//! \brief Computes orbital radius of the target body.
//!
//! \b Inputs
//! \arg #naifId    
//! \arg #startTime 
//!
//! \retval #orbitalRadius 
void Ephemeris::updateOrbitalRadius()
{
  ConstSpiceChar    ref[10] = "IAU_SUN";
  SpiceInt          obs = 10;
  SpiceDouble       state[6];
  SpiceDouble       lt;

  // SPICE CALL: Compute the position/velocity of the target body relative to the Sun.
  // Args: target body id, time, fixed body frame, observing body id,
  //       position/velocity state vector, light time
  spkgeo_c(naifId, startTime.getSpiceTime(), ref, obs, state, &lt);

  // Test for SPICE failure
  if (failed_c()) {
    throw string("Error in orbital radius calculation.");
  }

  // Compute the distance of the target body from the Sun.
  orbitalRadius = sqrt(state[0] * state[0] + state[1] * state[1] + state[2] * state[2]) * (1.0 / KILOMETERS_PER_AU);
}

//! \brief Computes the solar zenith angle at the specified lat/lon location on the target body.
//!
//! If the latitude has been set to a value in [-PI, PI], the the solar zenith angle will be computed.
//! If the latitude is outside of the accepted range, then this computation will be skipped.
//!
//! \b Inputs
//! \arg #naifName 
//! \arg #frameName
//! \arg #startTime
//! \arg #longitude
//! \arg #latitude 
//!
//! \retval #solarZenithAngle 
void Ephemeris::updateSolarZenithAngle()
{
  if (latitude > 90.0 || latitude < -90.0) {
    return;
  }

  const size_t SIZE = 256;

  SpiceChar       smethod[SIZE] = "ELLIPSOID";
  SpiceInt        npts = 1;
  SpiceDouble     lonlat[1][2] = { { (double)toRadians(longitude), (double)toRadians(latitude) } };
  SpiceDouble     srfpts[1][3];

  // SPICE CALL: Convert lat/lon to a cartesian surface point on the target body
  // Args: Computation method, time, target fixed frame, number of points, 
  //       lon/lat array, surface point array
  latsrf_c(smethod, naifName, startTime.getSpiceTime(), frameName, npts, lonlat, srfpts);

  SpiceChar    abcorr[SIZE] = "NONE";
  SpiceChar    imethod[SIZE] = "ELLIPSOID";
  SpiceChar    iobsrvr[SIZE] = "10";
  SpiceDouble  itrgepc;
  SpiceDouble  isrfvec[3];
  SpiceDouble  iphase;
  SpiceDouble  iincdnc;
  SpiceDouble  iemissn;

  // SPICE CALL: Find the ilumination angle at the surface point on the target body.
  // Args: Computation method, body name, time, target fixed frame, aberration correction, observing body name, surface point array, 
  //       surface point epoch, vector from observer to surface point, phase angle at surface point, solar incidence angle, emission angle
  ilumin_c(imethod, naifName, startTime.getSpiceTime(), frameName, abcorr, iobsrvr, srfpts[0],
    &itrgepc, isrfvec, &iphase, &iincdnc, &iemissn);

  // Test for SPICE failure
  if (failed_c()) {
    throw string("Error in solar zenith angle calculation.");
  }

  // Out of all that, we want the solar incidence angle at the surface point.
  solarZenithAngle = toDegrees(iincdnc);
}

//! \brief Sets the seconds per sol for the current body.
//!
//! Currently, this is a lookup of the seconds per sol for the current body.
//! It can be converted to a SPICE calculation if that is possible. 
//! \retval #secondsPerSol 
void Ephemeris::updateSecondsPerSol()
{
  switch (gramBody) {
  case VENUS:
    secondsPerSol = 1.00872e7;  // From NSSDCA Planetary Fact Sheet
    break;
  case EARTH:
    secondsPerSol = 86400.00;   // From NSSDCA Planetary Fact Sheet
    break;
  case MARS:
    secondsPerSol = 88774.92;   // From NSSDCA Planetary Fact Sheet
    break;
  case JUPITER:
    secondsPerSol = 35733.24;   // From NSSDCA Planetary Fact Sheet
    break;
  case URANUS:
    secondsPerSol = 62064.0;    // From NSSDCA Planetary Fact Sheet
    break;
  case NEPTUNE:
    secondsPerSol = 57996.0;    // From NSSDCA Planetary Fact Sheet
    break;
  case SATURN:
    secondsPerSol = 38361.6;    // From NSSDCA Planetary Fact Sheet
    break;
  case TITAN:
    secondsPerSol = 1377648.0;  // From NSSDCA Planetary Fact Sheet
    break;
  default:
    break;
  }
}

greal Ephemeris::getTropicalPeriod()
{
  greal tropicalPeriod;
  switch (gramBody) {
  case VENUS:
    tropicalPeriod = 224.695;  // From NSSDCA Planetary Fact Sheet
    break;
  case EARTH:
    tropicalPeriod = 365.242;  // From NSSDCA Planetary Fact Sheet
    break;
  case MARS:
    tropicalPeriod = 686.973;  // From NSSDCA Planetary Fact Sheet
    break;
  case JUPITER:
    tropicalPeriod = 4330.595;   // From NSSDCA Planetary Fact Sheet
    break;
  case URANUS:
    tropicalPeriod = 30588.740;    // From NSSDCA Planetary Fact Sheet
    break;
  case NEPTUNE:
    tropicalPeriod = 59799.9;    // From NSSDCA Planetary Fact Sheet
    break;
  case SATURN:
    tropicalPeriod = 10746.94;    // From NSSDCA Planetary Fact Sheet
    break;
  case TITAN:
    tropicalPeriod = 10746.94;  // From NSSDCA Planetary Fact Sheet
    break;
  default:
    tropicalPeriod = 0.0;
    break;
  }
  return tropicalPeriod;
}

//! \brief Find dates for tartget longitude of the sun and solar time.
//!
//! This routine will find the date closest to the startTime and the 
//! targetLongitudeSun at the current longitude.  Given that time, the 
//! occurences of the targetSolarTime immediately prior and following 
//! the LS time are found.  Values are returned in an arrays of size three containing
//! the prior TLST occurence, the LS occurence, and the latter TLST occrurence.
//!
//! \b Inputs
//! \arg startTime
//! \arg longitude
//! 
//! \param targetLongitudeSun  The desired longitude of the sun.
//! \param targetSolarTime     The desired solar time.
//! \param[out] gramTime       A GramTime array of size 3.
//! \param[out] lonSun         A greal array of size 3.
//! \param[out] tlst           A greal array of size 3.
//! \retval gramTime The times of the three desired states.
//! \retval lonSun   The longitude of the sun for each time.
//! \retval tlst     The true local solar time for each time.
void Ephemeris::findDates(greal targetLongitudeSun, greal targetSolarTime, GramTime gramTime[3], greal lonSun[3], greal tlst[3])
{
  // Set the accuracy desired.
  const greal epsilon = 1.0e-6;
  // Maximum loop iterations
  const int maxIter = 100;

  // Save the current time and state (reset on exit).
  GramTime saveTime = startTime;
  EphemerisState saveState = state;

  // Computations are performed in PET.
  // So zero out the OWLT and compute it later.
  startTime.setOneWayLightTime(0.0);

  // First, locate the target LS.
  // Start by finding the LS for the current time.
  greal ls = updateLongitudeOfTheSun();

  // Find out how far off the target LS is from the current LS.
  greal deltaLs = targetLongitudeSun - ls;
  // Shift to between -180 and 180 since the closest is desired.
  if (deltaLs < -180.0) {
    deltaLs += 360.0;
  }
  else if (deltaLs > 180.0) {
    deltaLs -= 360.0;
  }

  // Estimate the number of seconds per LS degree
  greal rateLs = getTropicalPeriod() * (86400.0 / 360.0); 

  // Convergence loop for LS
  int count = 0;
  // Stopping criteria: LS within tolerance or beyond max iterataions.
  while (abs(deltaLs) > epsilon && count < maxIter)
  {
    // If convergence is slow, then rateLs is probably too large.
    // This causes a yo-yo effect: -.05, 0.0499, -.0498, 0.0497, etc.
    // Reduce rateLs to avoid the yo-yo.
    if (count == 10) {
      rateLs *= 0.7;
    }

    // Make a new time estimate by adding the delta.
    startTime.setElapsedTime(startTime.getElapsedTime() + rateLs * deltaLs);

    // Find the new LS and update the delta
    ls = updateLongitudeOfTheSun();
    deltaLs = targetLongitudeSun - ls;
    if (deltaLs < -180.0) {
      deltaLs += 360.0;
    }
    else if (deltaLs > 180.0) {
      deltaLs -= 360.0;
    }
    // Increment the iteration counter.
    ++count;
  }

  // Compute TLST and OWLT and save the data.
  updateSubsolarLatLon();
  updateLocalSolarTime();
  updateOneWayLightTime();
  gramTime[1] = startTime;
  lonSun[1] = longitudeSun;
  tlst[1] = solarTime;

  // Next, find the target TLST occuring before the LS.
  // Computations are performed in PET. So zero out the OWLT.
  startTime.setOneWayLightTime(0.0);

  // The current solar time is at the desired LS.
  // See how far off this is from the target TLST.
  greal deltaSolarTime;
  if (targetSolarTime < solarTime) {
    deltaSolarTime = targetSolarTime - solarTime;
  }
  else {
    deltaSolarTime = targetSolarTime - solarTime - 24.0;
  }

  // Estimate the number of seconds per solar hour.
  updateSecondsPerSol();
  greal rateLTST = secondsPerSol / 24.0;
  size_t index = 0;
  // Retrograde planets reverse the solar clock.
  //if (gramBody == URANUS || gramBody == VENUS) {
  //  rateLTST *= -1.0;
  //  index = 2;
  //}

  // Convergence loop for prior LTST.
  count = 0;
  // Stopping criteria: TLST within tolerance or beyond max iterataions.
  while (abs(deltaSolarTime) > epsilon && count < maxIter)
  {
    // Make a new time estimate by adding the delta.
    startTime.setElapsedTime(startTime.getElapsedTime() + deltaSolarTime * rateLTST);
    
    // Update the solar time and recompute the delta.
    updateSubsolarLatLon();
    updateLocalSolarTime();
    deltaSolarTime = targetSolarTime - solarTime;

    // Increment the iteration counter.
    ++count;
  }

  // Compute LS and OWLT and save the data.
  updateLongitudeOfTheSun();
  updateOneWayLightTime();
  gramTime[index] = startTime;
  lonSun[index] = longitudeSun;
  tlst[index] = solarTime;

  // Finally, find the target TLST occuring after the LS.
  // Computations are performed in PET. So zero out the OWLT.
  startTime.setOneWayLightTime(0.0);

  // Estimate the new solar time by adding 24 solar hours to the time.
  startTime.setElapsedTime(startTime.getElapsedTime() + 24.0 * rateLTST);

  // Compute the TLST and see how far off this is.
  updateSubsolarLatLon();
  updateLocalSolarTime();
  deltaSolarTime = targetSolarTime - solarTime;

  // Convergence loop for the latter LTST.
  count = 0;
  // Stopping criteria: TLST within tolerance or beyond max iterataions.
  while (abs(deltaSolarTime) > epsilon && count < maxIter)
  {
    // Make a new time estimate by adding the delta.
    startTime.setElapsedTime(startTime.getElapsedTime() + deltaSolarTime * rateLTST);

    // Update the solar time and recompute the delta.
    updateSubsolarLatLon();
    updateLocalSolarTime();
    deltaSolarTime = targetSolarTime - solarTime;

    // Increment the iteration counter.
    ++count;
  }

  // Compute LS and OWLT and save the data.
  updateLongitudeOfTheSun();
  updateOneWayLightTime();
  gramTime[2-index] = startTime;
  lonSun[2-index] = longitudeSun;
  tlst[2-index] = solarTime;

  // Reset the time and state.
  startTime = saveTime;
  state = saveState;
}

} // namespace
