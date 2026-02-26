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

//! \brief An input state for an Atmosphere model
//!
//! The input state of an Atmosphere model is a point in space/time.
//! \ingroup CommonGRAM Cpp_Venus Cpp_Earth Cpp_Mars Cpp_Jupiter Cpp_Saturn Cpp_Uranus Cpp_Neptune Cpp_Titan
class Position
{
public:
  Position();
  Position(const Position& orig) = default;
  virtual ~Position() = default;

  greal height = 0.0;           //!< Height Above Reference geoid/ellipsoid \units{km}.
  greal latitude = 0.0;         //!< Geocentric latitude \units{\text{degrees}}.
  greal longitude = 0.0;        //!< East longitude positive \units{\text{degrees}}.
  greal elapsedTime = 0.0;      //!< Number of seconds from start time \units{\text{seconds}}.
  bool isPlanetoCentric = true; //!< True if planeto-centric, false if planeto-graphic.

  greal totalRadius = 0.0;      //!< Radius at current latitude and height \units{km}.
  greal latitudeRadius = 0.0;   //!< Radius at current latitude \units{km}.
  greal surfaceHeight = 0.0;    //!< Height of topographic surface above/below geoid at current lat/lon \units{km}.
  greal gravity = 0.0;          //!< The force of gravity at current position \units{m/s^2}.

  void setLongitude(greal lon, bool eastPositive = true) { longitude = (eastPositive ? lon : 360.0_deg - lon); }
  void setLongitudeDegrees(greal lon, bool eastPositive = true) { setLongitude(lon, eastPositive); }
  void setLongitudeRadians(greal lon, bool eastPositive = true) { setLongitude(toDegrees(lon), eastPositive); }

  void setLatitude(greal lat) { latitude = lat; }
  void setLatitudeDegrees(greal lat) { setLatitude(lat); }
  void setLatitudeRadians(greal lat) { setLatitude(toDegrees(lat)); }

  greal getLongitude(bool eastPositive = true) const { return ((eastPositive || longitude == 0.0) ? longitude : 360.0 - longitude); }
  greal getLongitudeDegrees(bool eastPositive = true) const { return getLongitude(eastPositive); }
  greal getLongitudeRadians(bool eastPositive = true) const { return toRadians(getLongitude(eastPositive)); }

  greal getLatitude() const { return latitude; }
  greal getLatitudeDegrees() const { return getLatitude(); }
  greal getLatitudeRadians() const { return toRadians(getLatitude()); }

  void reset();
  void convertToPlanetocentric(greal polarRadius, greal equatorialRadius);
  void convertToPlanetodetic(greal polarRadius, greal equatorialRadius);
};

} // namespace
