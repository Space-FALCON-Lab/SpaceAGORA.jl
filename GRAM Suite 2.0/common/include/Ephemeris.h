//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include "gram.h"
#include "EphemerisState.h"
#include "GramTime.h"

namespace GRAM {

//! \brief Computes ephemeris values for a target body.
//!
//! For a specified target body and time, the class updates the ephemeris state using Spice.
//! First, identify the target body with the setBody() method.  Then, prior to each
//! call to update(), set the time and longitude with setTime() and
//! setLongitude().  If the solar zenith angle is desired, then one must also
//! set the desired latitude with setLatitude().
//!
//! The updated ephemeris values can be accessed via getEphemerisState().
//! If a call to update() is preceeded by a call to setEphemerisState(), then the user
//! supplied ephemeris values will replace the Spice computations.
//! \ingroup CommonGRAM
class Ephemeris
{
public:
  Ephemeris();
  Ephemeris(const Ephemeris& orig);
  virtual ~Ephemeris() = default;

  void setBody(GRAM_BODY body);
  void setTime(const GramTime& tTime) { startTime = tTime; }
  void setLongitude(greal longitude) { this->longitude = wrapDegrees(longitude); }
  void setLatitude(greal latitude) { this->latitude = clamp(latitude, -90.0, 90.0); }
  void setEphemerisState(const EphemerisState& eState);
  void setFastModeOn(bool flag) { fastMode = flag; computeFlag = true; }
  void setComputeOn(bool flag) { computeFlag = flag; }
  void setSubsolarUpdateTime(greal utime) { subsolarUpdateTime = clamp(utime, 0.01, 3600.0); }

  virtual void update();

  GRAM_BODY getBody() const { return gramBody; }
  const EphemerisState& getEphemerisState() const { return state; }

  void findDates(greal targetLongitudeSun, greal targetSolarTime, GramTime gramTime[3], greal lonSun[3], greal tlst[3]);

protected:
  virtual void updateSubsolarLatLon();
  virtual void updateSolarDeclination();
  virtual void updateOneWayLightTime();
  virtual greal updateLongitudeOfTheSun();
  virtual void updateLocalSolarTime();
  virtual void updateOrbitalRadius();
  virtual void updateSolarZenithAngle();
  virtual void updateSecondsPerSol();
  greal getTropicalPeriod();
  void setBodyName(std::string name);
  void setBodyParentName(std::string name);

  GRAM_BODY gramBody = NO_BODY;         //!< The body this object will use for computations.
  long naifId = 0;                      //!< A body id used by NAIF Spice.
  static const int nameSize = 100;      //!< The max size of all name strings.
  char naifName[nameSize] = "";         //!< The Spice body name.
  char naifParentName[nameSize] = "";   //!< The Spice parent body name (for moons).
  char frameName[nameSize] = "";        //!< The Spice frame name  ("IAU_bodyname").

  greal longitude = 0.0;                //!< The longitude for solar time computations.
  greal latitude = 1000.0;              //!< The latitude for solar time computations.
  GramTime startTime;                   //!< The time for computations.
  bool userInputs = false;              //!< If true, user inputs override computations.
  bool fastMode = false;                //!< If true, ephemeris computations will be minimized.
  bool computeFlag = true;              //!< Controls computation of some ephemeris values. Values are computed when true.
  greal subsolarUpdateTime = 60.0;      //!< Interval in seconds in which to update subsolar lat and lon.
  greal subsolarLastTime = 0.0;         //!< Last time (in elapsed seconds) that subsolar lat and lon were updated.

  // Outputs
  EphemerisState state;                                    //!< The computed ephemeris state.
  greal& solarTime = state.solarTime;                      //!< \copydoc EphemerisState::solarTime
  greal& solarHourAngle = state.solarHourAngle;            //!< \copydoc EphemerisState::solarHourAngle
  greal& longitudeSun = state.longitudeSun;                //!< \copydoc EphemerisState::longitudeSun
  greal& subsolarLatitude = state.subsolarLatitude;        //!< \copydoc EphemerisState::subsolarLatitude
  greal& subsolarLongitude = state.subsolarLongitude;      //!< \copydoc EphemerisState::subsolarLongitude
  greal& solarDeclination = state.solarDeclination;        //!< \copydoc EphemerisState::solarDeclination
  greal& solarRightAscension = state.solarRightAscension;  //!< \copydoc EphemerisState::solarRightAscension
  greal& orbitalRadius = state.orbitalRadius;              //!< \copydoc EphemerisState::orbitalRadius
  greal& oneWayLightTime = state.oneWayLightTime;          //!< \copydoc EphemerisState::oneWayLightTime
  greal& solarZenithAngle = state.solarZenithAngle;        //!< \copydoc EphemerisState::solarZenithAngle
  greal& secondsPerSol = state.secondsPerSol;              //!< \copydoc EphemerisState::secondsPerSol
};

} // namespace
