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

namespace GRAM {

//! \brief Generate a Monte Carlo of profiles.
//!
//! The MonteCarlo class will generate a number of profiles each with a 
//! different random pertubation seed.  The class requires a Profile, such
//! as a Trajectory, and a PrinterProfile which have been initialized for use.
//! \ingroup CommonGRAM
class MonteCarlo
{
public:
  MonteCarlo();
  MonteCarlo(const MonteCarlo& orig) = delete;
  virtual ~MonteCarlo() = default;

  void setProfile(Profile& prof) { profile = &prof; }
  void setProfilePrinter(ProfilePrinter& profPrinter) { printer = &profPrinter; }

  void setInitialSeed(int seed) { initialSeed = seed;  useSeedList = false; }
  void setSampleSize(size_t size) { sampleSize = size;  }
  void setSeedList(const std::vector<int>& list) { seedList = list; useSeedList = true; }

  const std::vector<int>& getSeedList() const { return seedList; }

  void setInputParameters(const InputParameters& params);
  void setSeedParameters(int increment, int maximum) { seedIncrement = increment; maxSeed = maximum; }

  virtual void generate();

protected:
  virtual int getNextSeed(int seed, size_t index);

  Profile* profile = nullptr;         //!< The Profile used for data generation.
  ProfilePrinter* printer = nullptr;  //!< The ProfilePrinter used for data output.
  size_t sampleSize = 1;              //!< The number of profiles to generate.
  int initialSeed = 1001;             //!< The initial perturbation seed.
  std::vector<int> seedList;          //!< A list of pertubation seeds. Both input and output.
  bool useSeedList = false;           //!< Set to true if the user provides a seed list.  Otherwise, seeds will be auto-generated.
  int seedIncrement = 11;             //!< Increment used to arrive at the next random generator seed.
  int maxSeed = 16777216;             //!< Maximum allowed seed value.
};

} // namespace
