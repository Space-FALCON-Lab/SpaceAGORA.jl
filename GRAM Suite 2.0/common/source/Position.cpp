//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "Position.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
Position::Position()
{
}

//! \fn Position::Position(const Position& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn Position::~Position()
//! \copydoc Atmosphere::~Atmosphere()


//! \brief Resets values to zero.
void Position::reset()
{
  height = 0.0;          
  latitude = 0.0;        
  longitude = 0.0;       
  elapsedTime = 0.0;     
  isPlanetoCentric = true;    
  totalRadius = 0.0;     
  latitudeRadius = 0.0;  
  surfaceHeight = 0.0;   
  gravity = 0.0;
}

//! \brief Converts a planetodetic position to a planetocentric position.
//!
//! If the current position is flagged as planetodetic, then this method will convert
//! the position values to planetocentric.  The planetodetic values will be lost.
//! The latitude, height, totalRadius, and latitudeRadius are updated.
//!
//! \b Inputs
//! \arg #height      
//! \arg #latitude      
//! \arg #isPlanetoCentric      
//!
//! \param polarRadius       Planetary polar radius \units{km}.
//! \param equatorialRadius  Planetary equatorial radius \units{km}.
//!
//! \retval #height
//! \retval #latitude
//! \retval #totalRadius
//! \retval #latitudeRadius
void Position::convertToPlanetocentric(greal polarRadius, greal equatorialRadius)
{
  // Assumptions (checking occurs in debug mode only)
  assert(equatorialRadius != 0.0);
  assert(polarRadius != 0.0);

  // If the position is already planetocentric, then just return.
  if (!isPlanetoCentric) {
    // Special case for position at pole (do nothing)
    // Assume we have a non-polar latitude...
    if (abs(latitude) != 90.0_deg) {
      double omf = polarRadius / equatorialRadius; // one minus flattening
      double sphi = sin(toRadians(latitude));
      double cphi = cos(toRadians(latitude));
      double c = 1.0 / sqrt(square(cphi) + square(omf * sphi));

      // z-coordinate of position in rectangualr coordinates
      greal z_comp = (equatorialRadius * c * square(omf) + height) * sphi;

      // equatorial-conmponent of position in rectangualr coordfinates
      greal eq_comp = (equatorialRadius * c + height) * cphi;

      // Total radius and latitude from the rectangular coordinates
      totalRadius = sqrt(square(z_comp) + square(eq_comp));
      latitude = toDegrees(atan(z_comp / eq_comp));

      if (height == 0.0) {
        latitudeRadius = totalRadius;
      }
      else {
        // Planetocentric radius at latitude.
        greal s2lat = square(sin(toRadians(latitude)));
        latitudeRadius = equatorialRadius / sqrt(1.0 + (square(1.0 / omf) - 1.0) * s2lat);

        // Planetocentric height.
        height = totalRadius - latitudeRadius;
      }
    }
    else {
      latitudeRadius = polarRadius;
      totalRadius = latitudeRadius + height;
    }
    isPlanetoCentric = true;
  }
}

//! \brief Converts a planetocentric position to a planetodetic position.
//!
//! If the current position is flagged as planetocentric, then this method will convert
//! the position values to planetodetic.  The planetocentric values will be lost.
//! The latitude and height are updated.  The  totalRadius and latitudeRadius are set to zero.
//!
//! \b Inputs
//! \arg #latitude      
//! \arg #totalRadius      
//! \arg #isPlanetoCentric      
//!
//! \param polarRadius       Planetary polar radius \units{km}.
//! \param equatorialRadius  Planetary equatorial radius \units{km}.
//!
//! \retval #latitude
//! \retval #height
void Position::convertToPlanetodetic(greal polarRadius, greal equatorialRadius)
{
  if (isPlanetoCentric) {
    greal r = totalRadius * cos(toRadians(latitude));
    greal zin = totalRadius * sin(toRadians(latitude));
    if (abs(r) < 0.001) {
      latitude = copysign(90.0, zin);
      height = abs(zin) - polarRadius;
    }
    else {
      greal v;
      // Analytical solution for non - polar case
      // See also page K12 of Astronomical Almanac for iterative solution
      greal z = abs(zin);
      greal e = ((z + polarRadius) * polarRadius / equatorialRadius - equatorialRadius) / r;
      greal f = ((z - polarRadius) * polarRadius / equatorialRadius + equatorialRadius) / r;
      greal p = (e * f + 1.0) * 4.0 / 3.0;
      greal q = (e * e - f * f) * 2.0;
      greal d = p * p * p + q * q;
      if (d >= 0.0) {
        greal s = sqrt(d) + q;
        s = copysign(exp(log(abs(s)) / 3.0), s);
        v = p / s - s;
        v = -(q + q + v * v * v) / (3.0 * p);
      }
      else {
        v = 2.0 * sqrt((-p)) * cos(acos(q / p / sqrt((-p))) / 3.0);
      }
      greal g = 0.5 * (e + sqrt(e * e + v));
      greal t = sqrt(g * g + (f - v * g) / (g + g - e)) - g;
      greal lat = atan((1.0 - t * t) * equatorialRadius / (2.0 * polarRadius * t));
      height = (r - equatorialRadius * t) * cos(lat) + (z - polarRadius) * sin(lat);
      // Convert to degrees
      latitude = copysign(toDegrees(lat), zin);
    }
    latitudeRadius = totalRadius = 0.0;
    isPlanetoCentric = false;
  }
}

//! \fn Position::setLongitude(greal lon, bool eastPositive)
//! \brief Set the longitude in degrees.
//!
//! The longitude will be converted to east positive and stored in the Position object.
//! \param lon The longitude in degrees.
//! \param eastPositive Convention flag. True is east positive, false is west positive.

//! \fn Position::setLongitudeDegrees(greal lon, bool eastPositive)
//! \brief Set the longitude in degrees.
//!
//! The longitude will be converted to east positive and stored in the Position object.
//! \param lon The longitude in degrees.
//! \param eastPositive Convention flag. True is east positive, false is west positive.

//! \fn Position::setLongitudeRadians(greal lon, bool eastPositive)
//! \brief Set the longitude in radians.
//!
//! The longitude will be converted to degrees and stored in the Position object.
//! \param lon The longitude in radians.
//! \param eastPositive Convention flag. True is east positive, false is west positive.

//! \fn Position::setLatitude(greal lat)
//! \brief Set the latitude in degrees.
//!
//! \param lat The latitude in degrees.

//! \fn Position::setLatitudeDegrees(greal lat)
//! \brief Set the latitude in degrees.
//!
//! \param lat The latitude in degrees.

//! \fn Position::setLatitudeRadians(greal lat)
//! \brief Set the latitude in radians.
//!
//! The latitude will be converted to degrees and stored in the Position object.
//! \param lat The latitude in degrees.

//! \fn Position::getLongitude(bool eastPositive) const
//! \brief Gets the longitude in degrees.
//!
//! The longitude is retrieved using the desired east/west positive convention.
//! \param eastPositive Convention flag. True is east positive, false is west positive.
//! \returns The longitude in degrees.

//! \fn Position::getLongitudeDegrees(bool eastPositive) const
//! \brief Gets the longitude in degrees.
//!
//! The longitude is retrieved using the desired east/west positive convention.
//! \param eastPositive Convention flag. True is east positive, false is west positive.
//! \returns The longitude in degrees.

//! \fn Position::getLongitudeRadians(bool eastPositive) const
//! \brief Gets the longitude in radians.
//!
//! The longitude is retrieved using the desired east/west positive convention.
//! \param eastPositive Convention flag. True is east positive, false is west positive.
//! \returns The longitude in radians.

//! \fn Position::getLatitude()
//! \brief Gets the latitude in degrees.
//!
//! \returns The latitude in degrees.

//! \fn Position::getLatitudeDegrees()
//! \brief Gets the latitude in degrees.
//!
//! \returns The latitude in degrees.

//! \fn Position::getLatitudeRadians()
//! \brief Gets the latitude in radians.
//!
//! \returns The latitude in radians.

} // namespace