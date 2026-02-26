//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include "VerticalProfile.h"

using namespace std;

namespace GRAM {

//! \copydoc Profile::Profile()
VerticalProfile::VerticalProfile()
{
}

//! \copydoc Profile::Profile(const Profile& orig)
VerticalProfile::VerticalProfile(const VerticalProfile& orig) :
  Profile(orig)
{
  initialPosition = orig.initialPosition;
  deltaHeight = orig.deltaHeight;
  numPoints = orig.numPoints;
}


//! \copydoc Profile::~Profile()
VerticalProfile::~VerticalProfile()
{
}

//! \brief Set the change in height and number of data points.
//!
//! Each point in the profile will be computed from the initialPosition
//! by adding step times the deltaHeight where step is the current step.
//! \param deltaHeight The change in height for each step in the profile.
//! \param numPts The number of data points to generate.
void VerticalProfile::setHeightConditions(greal deltaHeight, int numPts)
{
  this->deltaHeight = deltaHeight;
  numPoints = max(numPts, 1);
}

//! \copydoc Profile::setInputParameters()
void VerticalProfile::setInputParameters(const InputParameters& params)
{
  // Make sure we have an atmosphere.
  if (atmosphere == nullptr) {
    throw(string("Error: No atmosphere model present in EphemerisProfile.\n"
      "       Call setAtmosphere() prior to calling setInputParameters().\n"
      "       This is an unrecoverable error."));
  }

  // Set the inital position
  Position initialPosition;
  initialPosition.height = params.initialHeight;
  initialPosition.setLatitudeDegrees(params.initialLatitude);
  initialPosition.setLongitude(params.initialLongitude, params.isEastLongitudePositiveOnInput);
  initialPosition.elapsedTime = 0.0;
  setInitialPosition(initialPosition);

  // Set the change in height
  setHeightConditions(params.deltaHeight, params.numberOfPositions);

  // Pass the input parameters on to the atmosphere.
  atmosphere->setInputParameters(params);
}

//!  \brief Generate inputs and run the model.
//!
//! This function will generate the positional inputs for a vertical height profile,
//! provide those inputs to the atmosphere model, and store the outputs.
//!
void VerticalProfile::generate() 
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
  Position position = initialPosition;

  // Start the generation loop.
  for (int step = 0; step < numPoints; ++step) {
    // Generate the new position by incrementing the height.
    position.height = initialPosition.height + step * deltaHeight;

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

} // namespace
