//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "EphemerisProfile.h"

using namespace std;

namespace GRAM {

//! \copydoc Profile::Profile()
EphemerisProfile::EphemerisProfile()
{
}

//! \copydoc Profile::Profile(const Profile& orig)
EphemerisProfile::EphemerisProfile(const EphemerisProfile& orig) :
  Profile(orig)
{
  initialPosition = orig.initialPosition;
  deltaPosition = orig.deltaPosition;
  initialEphemeris = orig.initialEphemeris;
  deltaEphemeris = orig.deltaEphemeris;
  numPoints = orig.numPoints;
}

//! \copydoc Profile::~Profile()
EphemerisProfile::~EphemerisProfile()
{
}

//! \copydoc Profile::setInputParameters()
void EphemerisProfile::setInputParameters(const InputParameters& params)
{
  // Make sure we have an atmosphere.
  if (atmosphere == nullptr) {
    throw(string("Error: No atmosphere model present in EphemerisProfile.\n"
                 "       Call setAtmosphere() prior to calling setInputParameters().\n"
                 "       This is an unrecoverable error."));
  }

  // Pass the input parameters on to the atmosphere.
  atmosphere->setInputParameters(params);
}

//!  \brief Generate inputs and run the model.
//!
//! This function will generate the positional inputs for a vertical height profile,
//! provide those inputs to the atmosphere model, and store the outputs.
//!
void EphemerisProfile::generate() 
{
  // Make sure we have an atmosphere to run.
  if (atmosphere == nullptr) {
    throw(string("Error: No atmosphere model present in EphemerisProfile.\n"
      "       Call setAtmosphere() prior to calling generate().\n"
      "       This is an unrecoverable error."));
  }

  // Make sure the profile data structure is cleared out.
  profile.clear();

  // Start with the initial position.
  EphemerisState ephem;
  Position position;

  // Start the generation loop.
  for (int step = 0; step < numPoints; ++step) {
    // Generate the new ephemeris by incrementing the values.
    ephem.longitudeSun = initialEphemeris.longitudeSun + step * deltaEphemeris.longitudeSun;
    ephem.solarTime = initialEphemeris.solarTime + step * deltaEphemeris.solarTime;
    ephem.subsolarLatitude = initialEphemeris.subsolarLatitude + step * deltaEphemeris.subsolarLatitude;
    ephem.subsolarLongitude = initialEphemeris.subsolarLongitude + step * deltaEphemeris.subsolarLongitude;
    ephem.solarZenithAngle = initialEphemeris.solarZenithAngle + step * deltaEphemeris.solarZenithAngle;

    // Give the ephemeris state to the atmosphere.
    atmosphere->setEphemerisState(ephem);

    // Compute the next position (p = p0 + step * delta)
    // Special case for the initial position.
    if (step == 0) {
      position.longitude = initialPosition.longitude;
      position.latitude = initialPosition.latitude;
      position.height = initialPosition.height;
      position.elapsedTime = initialPosition.elapsedTime;
    }
    // General case add the delta to the current position.
    else {
      position.longitude = wrapDegrees(position.longitude + deltaPosition.longitude);
      position.latitude += deltaPosition.latitude;
      position.height += deltaPosition.height;
      position.elapsedTime += deltaPosition.elapsedTime;
    }
    // Passage over the poles inverts the delta and flips the longitude.
    if (abs(position.latitude) > 90.0_deg) {
      position.latitude = copysign(180.0_deg, position.latitude) - position.latitude;
      position.longitude = wrapDegrees(position.longitude + 180.0_deg);
      deltaPosition.latitude *= -1.0;
    }

    // Decide if initial perturbations are to be updated.
    if (initialPerturbationsUpdated == false && step == 0) {
      atmosphere->setPerturbationAction(DO_NOT_UPDATE_PERTS);
    }
    else {
      atmosphere->setPerturbationAction(UPDATE_PERTS);
    }

    // Give the position to the atmosphere.
    atmosphere->setPosition(position);

    // Run the model.
    atmosphere->update();

    // Store the results.
    profile.emplace_back(atmosphere->getPosition(), atmosphere->getAtmosphereState(), atmosphere->getEphemerisState() );
  }
}

//! \fn void EphemerisProfile::setPosition(const Position& p)
//! \brief Sets the (unchanging) position of the profile.
//!
//! This method sets the initial position and zeroes out the delta position.
//! Call this method instead of calling setInitialPosition and setDeltaPosition.
//! \param p A Position object.

//! \fn void EphemerisProfile::setDeltaPosition(const Position& p)
//! \brief Sets the change in position of the profile.
//!
//! The position of a profile is computed as the initialPosition
//! plus step times the deltaPosition where step is the current step.
//! \param p A Position object.

//! \fn void EphemerisProfile::setEphemeris(const EphemerisState& eph)
//! \brief Sets the (unchanging) ephemeris of the profile.
//!
//! This method sets the initial ephemeris and zeroes out the delta ephemeris.
//! Call this method instead of calling setInitialEphemeris and setDeltaEphemeris.
//! \param eph An EphemerisState object.

//! \fn void EphemerisProfile::setInitialEphemeris(const EphemerisState& eph)
//! \brief Sets the starting ephemeris of the profile.
//!
//! \param eph An EphemerisState object.

//! \fn void EphemerisProfile::setDeltaEphemeris(const EphemerisState& eph)
//! \brief Sets the change in ephemeris of the profile.
//!
//! The ephemeris values of a profile is computed as the initialEphemeris
//! plus step times the deltaEphemeris where step is the current step.
//! \param eph An EphemerisState object.

//! \fn void EphemerisProfile::setNumberOfPoints(int numPts)
//! \brief Sets the number of data points to generate in the profile.
//!
//! \param numPts The number of data points to generate.

} // namespace
