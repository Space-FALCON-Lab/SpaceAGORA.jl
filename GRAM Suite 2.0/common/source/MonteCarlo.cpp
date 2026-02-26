//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "MonteCarlo.h"
#include "PerturbationState.h"
#include "ProfilePrinter.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
MonteCarlo::MonteCarlo()
{
}

//! \fn  MonteCarlo::~MonteCarlo()
//! \copydoc Atmosphere::~Atmosphere()

//! \brief Set the applicable input parameters.
//!
//! The routine copies the applicable members of the InputParameters
//! object into the appropriate members of this object.
//! \param params The input parameters.
void MonteCarlo::setInputParameters(const InputParameters& params)
{
  // Make sure we have a profile.
  if (profile == nullptr) {
    throw(string("Error: No profile present in MonteCarlo.\n"
                 "       Call setProfile() prior to calling setInputParameters().\n"
                 "       This is an unrecoverable error."));
  }
  // Make sure we have a profile printer.
  if (printer == nullptr) {
    throw(string("Error: No profile printer present in MonteCarlo.\n"
                 "       Call setProfilePrinter() prior to calling setInputParameters().\n"
                 "       This is an unrecoverable error.\n"));
  }

  // Set the initial random number seed.
  setInitialSeed(params.initialRandomSeed);

  // Set the number of samples to generate.
  setSampleSize(params.numberOfMonteCarloRuns);

  // Pass the input parameters on to the profile and profile printer.
  profile->setInputParameters(params);
  printer->setInputParameters(params);
}

//! \brief Generate the Monte Carlo data.
//!
//! This method will generate the specified number of profiles.
//! Each profile will be given a different perturbation seed.
//! The profile data will be passed to the ProfilePrinter for output.
//! This routine does not store the generated data.
void MonteCarlo::generate()
{
  // Make sure we have a profile.
  if (profile == nullptr) {
    throw(string("Error: No profile present in MonteCarlo.\n"
                 "       Call setProfile() prior to calling generate().\n"
                 "       This is an unrecoverable error."));
  }
  // Make sure we have a profile printer.
  if (printer == nullptr) {
    throw(string("Error: No profile printer present in MonteCarlo.\n"
                 "       Call setProfilePrinter() prior to calling generate().\n"
                 "       This is an unrecoverable error.\n"));
  }

  // Initialize
  if (printer != nullptr) {
    printer->openOutput();
    printer->printFileHeader(profile->getAtmosphere());
  }

  int seed = 0;
  if (useSeedList) {
    sampleSize = seedList.size();
  }
  else {
    seedList.clear();
  }

  // Step through number of Monte Carlo runs
  for (size_t run = 0; run < sampleSize; ++run) {

    // generate the next seed
    seed = getNextSeed(seed, run);

    // Initialize random number, position and time for each run
    profile->getAtmosphere().setSeed(seed);

    // Run the profile
    profile->generate();

    // Save/print the profile
    if (printer != nullptr) {
      printer->printSectionHeader(profile->getAtmosphere());
      printer->printData(profile->getProfile());
    }
  }

  // Close output
  if (printer != nullptr) {
    printer->closeOutput();
  }
}

//! \brief Returns a new pertubation seed.
//!
//! Given the previous seed and the current index, this method will
//! either generate a new seed or look up the next seed in the seed list.
//! If seeds are being generated they will be saved in the seed list.
int MonteCarlo::getNextSeed(int seedin, size_t index)
{
  int seed;
  if (useSeedList) {
    seed = seedList.at(index);
  }
  else {
    if (index == 0) {
      seed = initialSeed;
    }
    else {
      // Note that it is not necessary to randomly select each seed value
      // NR1 in order to get a random sequence of output.  Any regular
      // progression of selected NR1 values will do for this process.
      seed = seedin + seedIncrement;
      if (seed > maxSeed) {
        seed = seed % maxSeed;
      }
    }
    // Save seed value to list
    seedList.push_back(seed);
  }
  return seed;
}

//! \fn void MonteCarlo::setProfile(Profile& prof)
//! \brief Sets the profile used to generate data.
//!
//! \param prof A configured subclass of Profile.

//! \fn void MonteCarlo::setProfilePrinter(ProfilePrinter& profPrinter)
//! \brief Sets the profile printer used to output data.
//!
//! \param profPrinter A configured ProfilePrinter.

//! \fn void MonteCarlo::setInitialSeed(int seed)
//! \brief Sets the first pertubation seed.
//!
//! If the first perturbation seed is set via this method, then
//! subsequent seeds will be auto-generated.
//! \param seed An integer between 0 and 30000.

//! \fn void MonteCarlo::setSampleSize(size_t size)
//! \brief Sets the number of profiles to generate.
//!
//! \param size The number of profiles.

//! \fn void MonteCarlo::setSeedList(const std::vector<int>& list)
//! \brief Sets a list of pertubation seeds to use.
//!
//! If a list of perturbations seeds (integers between 0 and 30000)
//! are provided, then the sample size is assumed to be equivalent to
//! the list size.
//! \param list A vector of integers between 0 and 30000.

//! \fn MonteCarlo::getSeedList()
//! \brief Gets the list of pertubation seeds used in the last run.
//!
//! \returns A vector of integers.

} // namespace
