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
#include <algorithm>
#include "gram.h"

namespace GRAM {

//! \brief This class stores time and manages time conversions.
//!
//! The GramTime class manages the initial and elaspsed time of a simulation.
//! It also performs time scale and time frame conversions as necessary.
//!
//! This class uses the Spice library routines for conversions of time scale and frame.
//! Thus a call must be made to loadSpiceData() during the simulation initialization and
//! before any other calls to GramTime methods.
//!
//! The start time is set using setStartTime().  The initial time may be provided in a variety of
//! time scales, GRAM_TIME_SCALE, and time frames, GRAM_TIME_FRAME.
//!
//! Once the start time is set, time is propogated using setElapsedTime().  The
//! elasped time is the number of seconds since the start time.
//! \ingroup CommonGRAM Cpp_Venus Cpp_Earth Cpp_Mars Cpp_Jupiter Cpp_Saturn Cpp_Uranus Cpp_Neptune Cpp_Titan
class GramTime
{
public:
  GramTime();
  GramTime(const GramTime& orig) = default;
  virtual ~GramTime() = default;

  void setStartTime(int year, int month, int day, int hour, int minute, double seconds, GRAM_TIME_SCALE scale, GRAM_TIME_FRAME frame);
  void setStartTime(double timeInDays, GRAM_TIME_SCALE scale, GRAM_TIME_FRAME frame);
  void setElapsedTime(double seconds) { elapsedTime = seconds; }
  double getElapsedTime() { return elapsedTime; }

  void getStartTime(GRAM_TIME_SCALE scale, GRAM_TIME_FRAME frame, int &year, int &month, int &day, int &hour, int &minute, double &seconds) const;
  void getStartTime(GRAM_TIME_SCALE scale, GRAM_TIME_FRAME frame, double &timeInDays) const;

  void getTime(GRAM_TIME_SCALE scale, GRAM_TIME_FRAME frame, int &year, int &month, int &day, int &hour, int &minute, double &seconds) const;
  void getTime(GRAM_TIME_SCALE scale, GRAM_TIME_FRAME frame, double &timeInDays) const;
  double getDayOfYear() const;

  double getSpiceTime() const; 
  void setSpiceTime(double time) { spiceEphemerisTime = time; }

  double getOneWayLightTime() const { return oneWayLightTime; }
  void setOneWayLightTime(double owlt) { oneWayLightTime = std::max(owlt, 0.0); }

  GRAM_TIME_SCALE getTimeScale() const { return timeScale; }
  GRAM_TIME_FRAME getTimeFrame() const { return timeFrame; }

  const std::string& getScaleString() const { return getScaleString(timeScale); }
  const std::string& getFrameString() const { return getFrameString(timeFrame); }

private:
  const std::string& getScaleString(GRAM_TIME_SCALE scale) const;
  const std::string& getFrameString(GRAM_TIME_FRAME frame) const;

  GRAM_TIME_SCALE timeScale = COORDINATED_UNIVERSAL_TIME;  //!< The start time scale.
  GRAM_TIME_FRAME timeFrame = EARTH_RECEIVE_TIME;          //!< The start time frame.
  double spiceEphemerisTime = 0.0;    //!< TDB seconds past J2000 in the PET frame.
  double oneWayLightTime = 0.0;       //!< Seconds for light to pass from the body to Earth.
  double elapsedTime = 0.0;           //!< Seconds of elapsed mission time since the start time.

  static constexpr double secondsToDays = 1.0 / 86400.0;    //!< Conversion factor.
  static constexpr double minutesToDays = 1.0 / 1440.0;     //!< Conversion factor.
  static constexpr double hoursToDays = 1.0 / 24.0;         //!< Conversion factor.
  static constexpr double J2000 = 2451545.0;                //!< Noon on 1/1/2000 in Julian days.
  static constexpr double daysToCenturies = 1.0 / 36525.0;  //!< Conversion factor.

};

} // namespace
