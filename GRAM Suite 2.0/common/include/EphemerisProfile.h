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
#include "gram.h"
#include "Profile.h"
#include "Position.h"

namespace GRAM {

//! \brief Generate a profile using a sweep of ephemeris values.
//!
//! An ephemeris profile is defined by specifying an initial position and initial ephemeris values.
//! A change in ephemeris values and the number of data points must also be provided.  
//! Please note that these ephemeris values override the computed Spice values for that position and time.
//! The EphemerisProfile can be used to observe the affects of various ephemeris on the atmosphere model.
//! For example, one might sweep the longitude of the sun from 0 to 360 to see how the model
//! varies over the seasons.
//! \ingroup CommonGRAM
class EphemerisProfile : public Profile
{
public:
  EphemerisProfile();
  EphemerisProfile(const EphemerisProfile& orig);
  virtual ~EphemerisProfile() override;

  void setPosition(const Position& p) { initialPosition = p; deltaPosition.reset(); }
  void setDeltaPosition(const Position& p) { deltaPosition = p; }

  void setEphemeris(const EphemerisState& eph) { initialEphemeris = eph;  deltaEphemeris.reset();  }
  void setInitialEphemeris(const EphemerisState& eph) { initialEphemeris = eph; }
  void setDeltaEphemeris(const EphemerisState& eph) { deltaEphemeris = eph; }

  void setNumberOfPoints(int numPts) { numPoints = std::max(numPts, 1); }
  void setInitialPerturbationsUpdated(bool flag) { initialPerturbationsUpdated = flag; }

  void setInputParameters(const InputParameters& params) override;

  void generate() override;

protected:
  Position deltaPosition;           //!< The delta position for each step of the profile.
  EphemerisState initialEphemeris;  //!< The starting ephemeris of the profile.
  EphemerisState deltaEphemeris;    //!< The delta ephemeris for each step of the profile.
  int numPoints = 21;               //!< The number of data points in the profile.
  bool initialPerturbationsUpdated = true; //!< If false, initial perturbations will not be updated.
};

} // namespace
