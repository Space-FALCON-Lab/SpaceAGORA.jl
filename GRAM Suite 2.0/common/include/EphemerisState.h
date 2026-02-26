//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"

namespace GRAM {

//! \brief The current ephemeris values.
//!
//! This data class contains the current planetary ephemeris values computed by the Ephemeris class.
//! This class can also be used to supply ephemeris values to the model.
//! \ingroup CommonGRAM Cpp_Venus Cpp_Earth Cpp_Mars Cpp_Jupiter Cpp_Saturn Cpp_Uranus Cpp_Neptune Cpp_Titan
class EphemerisState
{
public:
  EphemerisState();
  EphemerisState(const EphemerisState& orig) = default;
  virtual ~EphemerisState() = default;

  greal solarTime = 0.0;            //!< Local solar time                                                  \note was TLOCAL, LST.
  greal solarHourAngle = 0.0;       //!< Angular measure of solar time. \units{\text{degrees}}.
  greal longitudeSun = 0.0;         //!< The planetocentric orbital longitude of the Sun \units{\text{degrees}}.
  greal subsolarLatitude = 0.0;     //!< Latitude of sub-solar point on the surface \units{\text{degrees}} \note was sunlat.
  greal subsolarLongitude = 0.0;    //!< Longitude of sub-solar point \units{\text{degrees}}               \note was sunlon.
  greal solarDeclination = 0.0;     //!< Declination of the sun \units{\text{degrees}}.
  greal solarRightAscension = 0.0;  //!< Right ascension of the sun \units{\text{degrees}}.
  greal solarZenithAngle = 0.0;     //!< Solar zenith angle \units{\text{degrees}}.
  greal orbitalRadius = 0.0;        //!< orbital radius \units{AU}                                         \note was rpl.
  greal oneWayLightTime = 0.0;      //!< Earth-Planet one-way light-time \units{\text{minutes}}            \note was owlt.
  greal secondsPerSol = 0.0;        //!< Seconds per sol \units{\text{seconds}}.

  greal getSubsolarLongitude(bool eastPositive = true) const { return ((eastPositive || subsolarLongitude == 0.0) ? subsolarLongitude : 360.0 - subsolarLongitude); }
  void reset();

  greal longitudeAscendingNode = 0;  //!< deprecated
  greal inclination = 0;  //!< deprecated
  greal argumentPerihelion = 0;  //!< deprecated
  greal semimajorAxis = 0;  //!< deprecated
  greal eccentricity = 0;  //!< deprecated
  greal meanAnomaly = 0;  //!< deprecated

};

} // namespace
