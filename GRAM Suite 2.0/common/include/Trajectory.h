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

namespace GRAM {

//! \brief Generate a trajectory profile.
//!
//! A trajectory profile is defined by specifying an initial position, a change in position, and
//! the number of data points desired.  Alternatively, an input file can be specified that contains
//! the positions for the trajectory.
//! \ingroup CommonGRAM
class Trajectory : public Profile
{
public:
  Trajectory();
  Trajectory(const Trajectory& orig);
  virtual ~Trajectory() override;

  void setDeltaPosition(const Position& p) { deltaPosition = p; useTrajFile = false; }
  void setNumberOfPoints(int numPts) { numPoints = std::max(numPts, 1); useTrajFile = false; }
  void setInitialPerturbationsUpdated(bool flag) { initialPerturbationsUpdated = flag; }

  void setDataFile(const std::string& fileName);

  void setInputParameters(const InputParameters& params) override;

  virtual void generate() override;

protected:
  virtual void getPosition(int step, Position& position);

  bool isPlantoCentric = true;     //!< True for planeto-centric input data, false for planeto-graphic.
  bool useTrajFile = false;        //!< If true, read positions from an input file.
  Position deltaPosition;          //!< The delta position for each step of the trajectory.
  int numPoints = 21;              //!< The number of data points in the trajectory.
  bool initialPerturbationsUpdated = true; //!< If false, initial perturbations will not be updated.

  std::string trajectoryFileName;       //!< The full path to a trajectory input file.
  std::vector<Position> trajectoryData; //!< The input positions from a trajectory file.
};

} // namespace
