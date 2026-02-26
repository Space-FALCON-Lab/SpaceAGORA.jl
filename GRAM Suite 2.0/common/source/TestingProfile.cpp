//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include "TestingProfile.h"

using namespace std;

namespace GRAM {

//! \copydoc Profile::Profile()
TestingProfile::TestingProfile()
{
  deltaPosition.height = 10.0;
}

//! \copydoc Profile::Profile(const Profile& orig)
TestingProfile::TestingProfile(const TestingProfile& orig) :
  Profile(orig)
{
  useTrajFile = orig.useTrajFile;
  initialPosition = orig.initialPosition;
  deltaPosition = orig.deltaPosition;
  numPoints = orig.numPoints;
  trajectoryFileName = orig.trajectoryFileName;
  trajectoryData = orig.trajectoryData;
}

//! \copydoc Profile::~Profile()
TestingProfile::~TestingProfile()
{
  trajectoryData.clear();
}

//! \copydoc Profile::setInputParameters()
void TestingProfile::setInputParameters(const InputParameters& params)
{
  // Make sure we have an atmosphere.
  if (atmosphere == nullptr) {
    throw(string("Error: No atmosphere model present in EphemerisProfile.\n"
      "       Call setAtmosphere() prior to calling setInputParameters().\n"
      "       This is an unrecoverable error."));
  }

  // If there is no traj file, then assume they will supply the position info.
  if (params.numberOfPositions > 0 || params.trajectoryFileName.empty()) {
    // Set the inital position
    Position initialPosition;
    initialPosition.height = params.initialHeight;
    initialPosition.setLatitudeDegrees(params.initialLatitude);
    initialPosition.setLongitude(params.initialLongitude, params.isEastLongitudePositiveOnInput);
    initialPosition.elapsedTime = 0.0;
    setInitialPosition(initialPosition);

    // Set the change in position
    Position deltaPosition;
    if (params.isEastLongitudePositiveOnInput) {
      deltaPosition.setLongitude(params.deltaLongitude);
    }
    else {
      deltaPosition.setLongitude(-params.deltaLongitude);
    }
    deltaPosition.setLatitudeDegrees(params.deltaLatitude);
    deltaPosition.height = params.deltaHeight;
    deltaPosition.elapsedTime = params.deltaTime;
    setDeltaPosition(deltaPosition);

    // The number of data points in each trajectory
    setNumberOfPoints(params.numberOfPositions);
  }
  else {
    // HACK!
        // Set the change in position
    Position deltaPosition;
    if (params.isEastLongitudePositiveOnInput) {
      deltaPosition.setLongitude(params.deltaLongitude);
    }
    else {
      deltaPosition.setLongitude(-params.deltaLongitude);
    }
    deltaPosition.setLatitudeDegrees(params.deltaLatitude);
    deltaPosition.height = params.deltaHeight;
    deltaPosition.elapsedTime = params.deltaTime;
    setDeltaPosition(deltaPosition);
    // If a traj file name is present, use it.
    setDataFile(params.trajectoryFileName);


  }
  
  // Pass the input parameters on to the atmosphere.
  atmosphere->setInputParameters(params);
}

//! \brief Set the trajectory data file and read the data.
//!
//! The specified trajectory file should contains the following fields separated by
//! whitespace (no commas): elapsed time in seconds, height in km, latitude in degrees,
//! and longitude in degrees.
//!
void TestingProfile::setDataFile(const std::string& fileName)
{
  try {
    // Make sure the traj vector is cleared out.
    trajectoryData.clear();

    // Save the file name.
    trajectoryFileName = fileName;
    // Open the file for reading.
    ifstream trajFile(trajectoryFileName);
    if (!trajFile.is_open()) {
      throw string("Unable to open file.");
    }

    // Read until the end of the file is reached.
    while (!trajFile.eof()) {
      Position p;
      EphemerisState e;
      // Read in a point.
      string lineBuffer = "";
      getline(trajFile, lineBuffer);

      // Replace CSV with space delimitted.
      replace(lineBuffer.begin(), lineBuffer.end(), ',', ' ');

      if (trajFile.bad()) {
        throw string("Unable to parse line after time = ") + to_string(trajectoryData.back().first.elapsedTime);
      }

      if (!lineBuffer.empty() && lineBuffer.find_first_not_of(" \t") != lineBuffer.npos) {
        istringstream lineInput(lineBuffer);
        lineInput >> p.elapsedTime >> p.height >> p.latitude >> p.longitude
          >> e.solarTime >> e.longitudeSun >> e.subsolarLatitude >> e.subsolarLongitude 
          >> e.orbitalRadius >> e.oneWayLightTime >> e.solarZenithAngle >> e.secondsPerSol;
        if (!lineInput.fail()) {
          // Store the point.
          trajectoryData.emplace_back(make_pair(p,e));
        }
        else {
          throw string("Unable to parse line: ") + lineBuffer;
        }
      }
    }

    // Close the trajectory file.
    trajFile.close();
    // Save the size of the data.
    numPoints = (int)trajectoryData.size();
    // Mark the source as coming from a traj file.
    useTrajFile = true;
  }
  catch (string msg) {
    throw(string("Error: StateProfile file loading error.\n       ") + msg
      + "\n       Trajectory file: " + trajectoryFileName
      + "\n       This is an unrecoverable error.");
  }
}

//!  \brief Get inputs and run the model.
//!
//! This function will get the positional inputs for a StateProfile,
//! provide those inputs to the atmosphere model, and store the outputs.
//!
void TestingProfile::generate()
{
  // Make sure we have an atmosphere to run.
  if (atmosphere == nullptr) {
    throw(string("Error: No atmosphere model present in EphemerisProfile.\n"
      "       Call setAtmosphere() prior to calling generate().\n"
      "       This is an unrecoverable error."));
  }

  // Make sure the profile data structure is cleared out.
  profile.clear();

  // Prime the deltas with the elapsed time.
  Position position;
  EphemerisState ephem;
  if (false) {
//  if (useTrajFile) {
    if (trajectoryData.size() > 1) {
      // Find the delta between the first two positions.
      getPosition(0, position, ephem);
      Position position2;
      EphemerisState ephem2;
      getPosition(1, position2, ephem2);
      position.elapsedTime = position2.elapsedTime - position.elapsedTime;
    }
    else {
      // Anomolous trajectory file.
      position.elapsedTime = 0.0;
    }
    // Only keep the elapsed time.
    position.height = 0.0;
    position.latitude = 0.0;
    position.longitude = 0.0;
  }
  else {
    // For generated trajectories, just use the delta provided.
    position.elapsedTime = deltaPosition.elapsedTime;
    position = deltaPosition;
  }
  atmosphere->setDelta(position);

  // Start the generation loop.
  for (int step = 0; step < numPoints; ++step) {
    // Get the next position
    getPosition(step, position, ephem);

    // Decide if initial perturbations are to be updated.
    if (initialPerturbationsUpdated == false && step == 0) {
      atmosphere->setPerturbationAction(DO_NOT_UPDATE_PERTS);
    }
    else {
      atmosphere->setPerturbationAction(UPDATE_PERTS);
    }

    // Give the position to the atmosphere.
    atmosphere->setPosition(position);
    atmosphere->setEphemerisState(ephem);

    // Run the model.
    atmosphere->update();

    // Store the results.
    profile.emplace_back(atmosphere->getPosition(), atmosphere->getAtmosphereState(), atmosphere->getEphemerisState() );
  }
}

//! \brief Read or generate a position.
//!
//! A new position for the specified step will either be read from a trajectory file
//! or generated at a specified delta from the previous position.  A step value of -1
//! may be used to generate a position prior to the initial position.  Positions that
//! are generated using a specified delta must be generated in order.
//!
//! /param step The index of the point in the trajectory.
void TestingProfile::getPosition(int step, Position& position, EphemerisState& ephem)
{
  // When data is loaded from a trajectory file:
  if (useTrajFile) {
    position = trajectoryData.at(step).first;
    ephem = trajectoryData.at(step).second;
    if (atmosphere->getInputParameters().isEastLongitudePositiveOnInput == false) {
      position.longitude = 360.0_deg - position.longitude;
      ephem.subsolarLongitude = 360.0_deg - ephem.subsolarLongitude;
//      ephem.longitudeSun = 360.0_deg - ephem.longitudeSun;
    }
  }
  // When data is generated at fixed deltas:
  else {
    // Compute the next position (p = p0 + step * delta)
    // Special case for the initial position.
    if (step == 0) {
      position.longitude = initialPosition.longitude;
      position.latitude = initialPosition.latitude;
      position.height = initialPosition.height;
      position.elapsedTime = initialPosition.elapsedTime;
    }
    // Special case for position prior to initial position.
    else if (step == -1) {
      position.longitude = wrapDegrees(initialPosition.longitude - deltaPosition.longitude);
      position.latitude = initialPosition.latitude - deltaPosition.latitude;
      position.height = initialPosition.height - deltaPosition.height;
      position.elapsedTime = initialPosition.elapsedTime - deltaPosition.elapsedTime;
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
  }
}

//! \fn void TestingProfile::setDeltaPosition(const Position& p)
//! \brief Sets the change in position of the trajectory.
//!
//! The position of a trajectory is computed as the initialPosition
//! plus i times the deltaPosition where i is the current step.
//! \param p A Position object.

//! \fn void TestingProfile::setNumberOfPoints(int numPts)
//! \brief Sets the number of data points to generate in the trajectory.
//!
//! \param numPts The number of data points to generate.


} // namespace
