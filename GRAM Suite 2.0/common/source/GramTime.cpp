//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>
#include <iomanip>
#include "GramTime.h"
#include "SpiceUsr.h"
#include "error_strings.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
GramTime::GramTime()
{
}

//! \fn GramTime::GramTime(const GramTime& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn GramTime::~GramTime()
//! \copydoc Atmosphere::~Atmosphere()

//! \brief Sets the simulation start time (epoch).
//!
//! \param year   The simulation start year (4 digits).
//! \param month  The simulation start month (1 to 12).
//! \param day    The simulation start day (1 to 31).
//! \param hour   The simulation start hour (0 to 23).
//! \param minute  The simulation start minute (0 to 59).
//! \param seconds  The simulation start seconds (decimal).
//! \param scale  The input time scale.
//! \param frame  The input time frame.
void GramTime::setStartTime(int year, int month, int day, int hour, int minute, double seconds, GRAM_TIME_SCALE scale, GRAM_TIME_FRAME frame)
{
  // Save the time scale and frame
  timeScale = scale;
  timeFrame = frame;

  // Build a string containing the time
  // (yup, that's the way Spice wants it.  Don't blame me.)
  std::ostringstream timeString;
  timeString << year << "-" << month << "-" << day
    << " " << hour << ":" << minute << ":" << fixed << seconds
    << " " << getScaleString(scale);

  // SPICE CALL: convert the time string to Spice time (TBD/PET in seconds).
  str2et_c(timeString.str().c_str(), &spiceEphemerisTime);
  if (failed_c()) {
    throw(SPICE_ERROR_MESSAGE + "Cannot convert time string."
      + "\n       Time string: " + timeString.str()
      + "\n       This is an unrecoverable error.");
  }
}

//! \brief Sets the simulation start time (epoch).
//!
//! \param timeInDays The simulation start time in days since J2000.
//! \param scale  The input time scale.
//! \param frame  The input time frame.
void GramTime::setStartTime(double timeInDays, GRAM_TIME_SCALE scale, GRAM_TIME_FRAME frame)
{
  // Save the time scale and frame
  timeScale = scale;
  timeFrame = frame;

  // Build a J2000 Julian date string containing the time
  // (yup, that's the way Spice wants it.  Don't blame me.)
  std::ostringstream timeString;
  timeString << "JD " << fixed << setprecision(10) << timeInDays 
    << " " << getScaleString(scale);

  // SPICE CALL: convert the time string to Spice time (TBD/PET in seconds).
  str2et_c(timeString.str().c_str(), &spiceEphemerisTime);
  if (failed_c()) {
    throw(SPICE_ERROR_MESSAGE + "Cannot convert time string."
      + "\n       Time string: " + timeString.str()
      + "\n       This is an unrecoverable error.");
  }
}

//! \brief Gets the simulation start time (epoch) in month, day, ..., seconds format.
//!
//! \param scale The desired time scale.
//! \param frame The desired time frame.
//! \param[out] year    An integer.
//! \param[out] month   An integer.
//! \param[out] day     An integer.
//! \param[out] hour    An integer.
//! \param[out] minute  An integer.
//! \param[out] seconds A double.
//! \retval year   The simulation start year (4 digits).
//! \retval month  The simulation start month (1 to 12).
//! \retval day    The simulation start day (1 to 31).
//! \retval hour   The simulation start hour (0 to 23).
//! \retval minute  The simulation start minute (0 to 59).
//! \retval seconds  The simulation start seconds (decimal).
void GramTime::getStartTime(GRAM_TIME_SCALE scale, GRAM_TIME_FRAME frame, int &year, int &month, int &day, int &hour, int &minute, double &seconds) const
{
  SpiceChar spiceString[50];
  // Describe how Spice should return the time.
  // Spice only returns time as strings.
  std::string picture = "YYYY MM DD HR MN SC.#### ::";
  // Tack on the time scale to the format.
  picture += getScaleString(scale);

  // SPICE CALL: Get the start time as a string.
  timout_c(spiceEphemerisTime, picture.c_str(), 50, spiceString);
  if (failed_c()) {
    throw(SPICE_ERROR_MESSAGE + "Cannot convert to time string."
      + "\n       Ephemeris time: " + to_string(spiceEphemerisTime)
      + "\n       This is an unrecoverable error.");
  }

  // Now parse the string into the return values.
  std::istringstream timeString(spiceString);
  timeString >> year >> month >> day
    >> hour >> minute >> seconds;
}

//! \brief Gets the simulation start time (epoch) in the Julian day format.
//!
//! \param scale The desired time scale.
//! \param frame The desired time frame.
//! \param[out] timeInDays A double.
//! \retval timeInDays   The simulation start time in days since J2000.
void GramTime::getStartTime(GRAM_TIME_SCALE scale, GRAM_TIME_FRAME frame, double &timeInDays) const
{
  SpiceChar spiceString[50];
  // Describe how Spice should return the time.
  // Spice really only returns time as strings.
  std::string picture = "JULIAND.######## ::";
  // Tack on the time scale to the format.
  picture += getScaleString(scale);

  // SPICE CALL: Get the start time as a string.
  timout_c(spiceEphemerisTime, picture.c_str(), 50, spiceString);
  if (failed_c()) {
    throw(SPICE_ERROR_MESSAGE + "Cannot convert to time string."
      + "\n       Ephemeris time: " + to_string(spiceEphemerisTime)
      + "\n       This is an unrecoverable error.");
  }

  // Now parse the string into the return values.
  std::istringstream timeString(spiceString);
  timeString >> timeInDays;
}

//! \brief Gets the current simulation time in month, day, ..., seconds format.
//!
//! \param scale The desired time scale.
//! \param frame The desired time frame.
//! \param[out] year    An integer.
//! \param[out] month   An integer.
//! \param[out] day     An integer.
//! \param[out] hour    An integer.
//! \param[out] minute  An integer.
//! \param[out] seconds A double.
//! \retval year   The current year (4 digits).
//! \retval month  The current month (1 to 12).
//! \retval day    The current day (1 to 31).
//! \retval hour   The current hour (0 to 23).
//! \retval minute  The current minute (0 to 59).
//! \retval seconds  The current seconds (decimal).
void GramTime::getTime(GRAM_TIME_SCALE scale, GRAM_TIME_FRAME frame, int &year, int &month, int &day, int &hour, int &minute, double &seconds) const
{
  // Describe how Spice should return the time.
  // Spice really ONLY returns time as strings.
  SpiceChar spiceString[50];
  std::string picture = "YYYY MM DD HR MN SC.#### ::";
  // Tack on the time scale to the format.
  picture += getScaleString(scale);

  // SPICE CALL: Get the current time as a string.
  timout_c(getSpiceTime(), picture.c_str(), 50, spiceString);
  if (failed_c()) {
    throw(SPICE_ERROR_MESSAGE + "Cannot convert to time string."
      + "\n       Ephemeris time: " + to_string(getSpiceTime())
      + "\n       This is an unrecoverable error.");
  }

  // Now parse the string into the return values.
  std::istringstream timeString(spiceString);
  timeString >> year >> month >> day
    >> hour >> minute >> seconds;
}

//! \brief Gets the zero-based day of the year.
//!
//! The calendar day of the year is computed with possible fractional component.
//! Note that Jan 1 is returned as 0.0.
//!
//! \returns The day of the year (0.0 - 366.0).
double GramTime::getDayOfYear() const
{
  // Describe how Spice should return the time.
  // Spice really ONLY returns time as strings.
  SpiceChar spiceString[50];
  std::string picture = "DOY.############ ::";
  // Tack on the time scale to the format.
  picture += getScaleString(timeScale);

  // SPICE CALL: Get the current time as a string.
  timout_c(getSpiceTime(), picture.c_str(), 50, spiceString);
  if (failed_c()) {
    throw(SPICE_ERROR_MESSAGE + "Cannot convert to day of year."
      + "\n       Ephemeris time: " + to_string(getSpiceTime())
      + "\n       This is an unrecoverable error.");
  }

  // Now parse the string into the return values.
  std::istringstream timeString(spiceString);
  greal dayOfYear;
  timeString >> dayOfYear;
  return dayOfYear - 1.0;
}

//! \brief Gets the current simulation time (epoch) in the Julian day format.
//!
//! \param scale The desired time scale.
//! \param frame The desired time frame.
//! \param[out] timeInDays A double.
//! \retval timeInDays   The current time in days since J2000.
void GramTime::getTime(GRAM_TIME_SCALE scale, GRAM_TIME_FRAME frame, double &timeInDays) const
{
  // Describe how Spice should return the time.
  // Spice unbelievably ONLY returns time as strings.
  SpiceChar spiceString[50];
  std::string picture = "JULIAND.######## ::";
  // Tack on the time scale to the format.
  picture += getScaleString(scale);

  // SPICE CALL: Get the current time as a string.
  timout_c(getSpiceTime(), picture.c_str(), 50, spiceString);
  if (failed_c()) {
    throw(SPICE_ERROR_MESSAGE + "Cannot convert to time string."
      + "\n       Ephemeris time: " + to_string(getSpiceTime())
      + "\n       This is an unrecoverable error.");
  }

  // Now parse the string into the return values.
  std::istringstream timeString(spiceString);
  timeString >> timeInDays;
}

//! \brief Return the time scale as a string.
//!
//! \param scale The time scale enum value.
//! \returns A three character string representation of the time scale.
const std::string& GramTime::getScaleString(GRAM_TIME_SCALE scale) const
{
  const static string scaleString[3] = { "UTC", "TDT", "TDB" };
  switch (scale) {
  case UTC:
  case COORDINATED_UNIVERSAL_TIME:
    return scaleString[0];
  case TDT:
  case TERRESTRIAL_DYNAMICAL_TIME:
    return scaleString[1];
  default:
    return scaleString[2];
  }
}

//! \brief Return the time frame as a string.
//!
//! \param frame The time frame enum value.
//! \returns A three character string representation of the time frame.
const std::string& GramTime::getFrameString(GRAM_TIME_FRAME frame) const
{
  const static string frameString[3] = { "PET", "ERT" };
  switch (frame) {
  case PET:
  case PLANET_EVENT_TIME:
    return frameString[0];
  case ERT:
  case EARTH_RECEIVE_TIME:
    return frameString[1];
  default:
    return frameString[1];
  }
}

//! \brief Get the current simulation time for Spice ephemeris computations.
//!
//! Spice uses the Barycentric Dynamical Time scale.  And the Ephemeris
//! routines should use the Planet Event Time frame.  This routine returns
//! the current simulation time (start time + elapsed time) in TDB/PET.
//! \returns The current time in TDB seconds past J2000 in the PET frame. 
double GramTime::getSpiceTime() const
{
  if (timeFrame == ERT || timeFrame == EARTH_RECEIVE_TIME) {
    // We are looking at the planet from earth and we see where
    // the planet was oneWayLightTime seconds ago.  So adjust to
    // make sure computations are done in planet event time.
    return spiceEphemerisTime + elapsedTime - oneWayLightTime;
  }
  return spiceEphemerisTime + elapsedTime;
}

//! \fn GramTime::setSpiceTime(double time)
//! \brief Set the start time using the NAIF Spice ephemeris time.
//!
//! Spice internally uses the Barycentric Dynamical Time scale.  The start time
//! can be set by providing the start time in TDB seconds past J2000 in the PET frame.
//! \param time The current time in TDB seconds past J2000 in the PET frame.

//! \fn GramTime::setElapsedTime(double seconds)
//! \brief Set the seconds since the start time.
//!
//! The current time is computed by adding the elapse time (in seconds)
//! to the start time.
//! \param seconds The number of seconds since the start time.

//! \fn GramTime::getElapsedTime()
//! \brief Gets the number of seconds since the start time.
//!
//! The current time is computed by adding the elapse time (in seconds)
//! to the start time.
//! \returns The number of seconds since the start time.

//! \fn GramTime::setOneWayLightTime(double owlt)
//! \brief Set the one-way light time for ERT computations.
//!
//! If the start time is entered as an Earth receive time, then the one-way light time
//! must be set before any call to getSpiceTime().
//! \param owlt The one-way light time (in seconds).

//! \fn GramTime::getOneWayLightTime()
//! \brief Gets the one-way light time used in ERT computations.
//!
//! If the start time is entered as an Earth receive time, then the one-way light time
//! must be set before any call to getSpiceTime().
//! \returns The one-way light time (in seconds).

//! \fn GramTime::getTimeScale()
//! \brief Gets the time scale (UTC, TDB) for the current time.
//!
//! \returns The time scale (UTC, TDB) for the current time.

//! \fn GramTime::getTimeFrame()
//! \brief Gets the time frame (ERT, PET) for the current time.
//!
//! \returns The time frame (ERT, PET) for the current time.

//! \fn GramTime::getScaleString() const
//! \brief Gets the time scale (UTC, TDB) for the current time as a string.
//!
//! \returns A three character string representation of the time scale.

//! \fn GramTime::getFrameString() const
//! \brief Gets the time frame (ERT, PET) for the current time as a string.
//!
//! \returns A three character string representation of the time frame.

} // namespace
