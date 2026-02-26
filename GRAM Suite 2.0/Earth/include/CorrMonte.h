//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include "gram.h"
#include "Profile.h"
#include "ProfilePrinter.h"
#include "Position.h"
#include "InputParameters.h"
#include "StateCorrelator.h"

namespace GRAM {

//! \brief Generate a Monte Carlo of correlated profiles.
//!
//! The CorrMonte class will generate a number of profiles each with a 
//! different random pertubation seed and correlated to a nominal profile.  
//! Hourly dispersions can be created by providing a fixed time increment.
//! The class requires a Profile, such as a Trajectory, and a PrinterProfile 
//! which have been initialized for use.
//! \ingroup EarthGRAM
class CorrMonte
{
public:
  CorrMonte();
  CorrMonte(const CorrMonte& orig) = delete;
  virtual ~CorrMonte() = default;

  void setCorrelator(StateCorrelator& corr) { correlator = &corr; }
  void setProfile(Profile& prof) { profile = &prof; }
  void setProfilePrinter(ProfilePrinter& profPrinter) { printer = &profPrinter; }

  void setInitialSeed(int seed) { initialSeed = seed; }
  void setSampleSize(size_t size) { sampleSize = size; }
  void setCorrDeltaHours(greal deltaHours) { deltaSeconds = deltaHours * 3600.0; }

  void setInputParameters(const InputParameters& params);

  virtual void generate();

protected:
  StateCorrelator *correlator = nullptr;   //!< Pointer to a state correlator object.
  Profile* profile = nullptr;              //!< The Profile used for data generation.
  ProfilePrinter* printer = nullptr;       //!< The ProfilePrinter used for data output.
  size_t sampleSize = 1;                   //!< The number of profiles to generate.
  int initialSeed = 1001;                  //!< The initial perturbation seed.
  greal deltaSeconds = 0.0;                //!< A time increment for createing hourly dispersions.
  InputParameters inputParameters;         //!< User supplied parameters.
};

} // namespace
