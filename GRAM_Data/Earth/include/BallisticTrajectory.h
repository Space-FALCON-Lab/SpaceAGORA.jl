//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <algorithm>
#include <vector>
#include "gram.h"
#include "Profile.h"
#include "StateCorrelator.h"

namespace GRAM {

//! \brief Generate a ballistic trajectory profile.
//!
//! A ballistic trajectory profile is defined by specifying an initial position, 
//! a change in position, and the number of data points desired.  The trajectory will
//! increase in height until either a maximum specified height is reached, or half of 
//! the total number of desired data points has been reached.
//! Alternatively, an input file can be specified that contains the positions for the trajectory.
//! This ingested trajectory is assumed to be initially strictly increasing and then
//! strictly decreasing.
//!
//! \ingroup EarthGRAM
class BallisticTrajectory : public Profile
{
public:
  BallisticTrajectory();
  BallisticTrajectory(const BallisticTrajectory& orig);
  virtual ~BallisticTrajectory() override;

  void setCorrelator(StateCorrelator& corr) { correlator = &corr; }
  void setMaximumHeight(greal height) { maxHeight = height; }
  void setInitialPosition(const Position& p) override { initialPosition = p; useTrajFile = false; }
  void setDeltaPosition(const Position& p) { deltaPosition = p; useTrajFile = false; }
  void setNumberOfPoints(int numPts) { numPoints = std::max(numPts, 1); useTrajFile = false; }

  void setDataFile(const std::string& fileName);

  virtual void setInputParameters(const InputParameters& params) override;

  virtual void generate() override;

protected:
  virtual void getPosition(int step, Position& position);
  const ProfileData& getAscendingData(greal height);

  StateCorrelator* correlator = nullptr;  //!< Pointer to a state correlator object.

  bool isPlantoCentric = true;  //!< True for planeto-centric input data, false for planeto-graphic.
  bool useTrajFile = false;     //!< If true, read positions from an input file.
  Position deltaPosition;       //!< The delta position for each step of the trajectory.
  int numPoints = 21;           //!< The number of data points in the trajectory.
  greal maxHeight = 0.0;        //!< The peak height of the trajectory \units{km}.
  size_t index = 0;             //!< The last recorded index of an ascent point.
  bool reachedApex = false;     //!< If false, the trajectory is increasing. When true, it is decreasing.

  std::string trajectoryFileName;       //!< The full path to a trajectory input file.
  std::vector<Position> trajectoryData; //!< The input positions from a trajectory file.
};

} // namespace
