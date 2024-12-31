//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: EarthGRAM
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "CorrMonte.h"
#include "PerturbationState.h"
#include "RandomNumberGenerator.h"

// Remove when this class is moved to common.
#include "EarthInputParameters.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
CorrMonte::CorrMonte()
{
}

//! \fn  CorrMonte::CorrMonte(const CorrMonte& orig)
//! \brief Copying this object is discouraged.

//! \fn  CorrMonte::~CorrMonte()
//! \copydoc Atmosphere::~Atmosphere()

//! \fn  CorrMonte::setCorrelator(StateCorrelator& corr)
//! \brief Set the state correlator.
//!
//! This object requires a StateCorrelator that has been subclassed for the 
//! specific planet.
//! \param corr A subclassed StateCorrelator object.

//! \brief Set the applicable input parameters.
//!
//! The routine copies the applicable members of the InputParameters
//! object into the appropriate members of this object.
//! \param params The input parameters.
void CorrMonte::setInputParameters(const InputParameters& params)
{
  // Make sure we have a profile.
  if (profile == nullptr) {
    throw(string("Error: No profile present in CorrMonte.\n"
                 "       Call setProfile() prior to calling setInputParameters().\n"
                 "       This is an unrecoverable error."));
  }
  // Make sure we have a profile printer.
  if (printer == nullptr) {
    throw(string("Error: No profile printer present in CorrMonte.\n"
                 "       Call setProfilePrinter() prior to calling setInputParameters().\n"
                 "       This is an unrecoverable error.\n"));
  }

  // For now, there are some Earth only parameters
  // Remove when this class is moved to common.
  const EarthInputParameters& eParams = static_cast<const EarthInputParameters&>(params);

  // Verify that a Corr Monte is desired.
  if (!eParams.corrMonte) {
    throw(string("Error: Inputs must specify use of CorrMonte.\n"
                 "       Add (CorrMonte = 1) to your namelist inputs.\n"
                 "       This is an unrecoverable error."));
  }

  // Set the time offset for hourly dispersions.
  setCorrDeltaHours(eParams.corrDeltaHours);

  // Set the initial random number seed.
  setInitialSeed(params.initialRandomSeed);

  // Set the number of samples to generate.
  setSampleSize(params.numberOfMonteCarloRuns);

  // Pass the input parameters on to the profile and profile printer.
  profile->setInputParameters(params);
  printer->setInputParameters(params);
  inputParameters = params;
}

//! \brief Generate the Corr Monte data.
//!
//! This method will generate the specified number of profiles.
//! Each profile will be given a different perturbation seed.
//! The profile data will be passed to the ProfilePrinter for output.
//! This routine does not store the generated data.
void CorrMonte::generate()
{
  // Make sure we have a profile.
  if (profile == nullptr) {
    throw(string("Error: No profile present in CorrMonte.\n"
                 "       Call setProfile() prior to calling generate().\n"
                 "       This is an unrecoverable error."));
  }
  // Make sure we have a profile printer.
  if (printer == nullptr) {
    throw(string("Error: No profile printer present in CorrMonte.\n"
                 "       Call setProfilePrinter() prior to calling generate().\n"
                 "       This is an unrecoverable error.\n"));
  }

  // Initialize
  if (printer != nullptr) {
    printer->setFileNamePrefix("Nominal_");
    printer->openOutput();
    printer->printFileHeader(profile->getAtmosphere());
  }

  // Generate the initial profile
  profile->getAtmosphere().setSeed(initialSeed);

  // Run the profile
  profile->generate();

  // Save/print the profile
  vector<ProfileData> nominalRun = profile->getProfile();
  if (printer != nullptr) {
    printer->printSectionHeader(profile->getAtmosphere());
    printer->printData(nominalRun);
  }

   // Close output
  if (printer != nullptr) {
    printer->closeOutput();
  }
  
  //GramTime time = profile->getAtmosphere().getStartTime();
  //double timeInDays;
  //time.getStartTime(UTC, PET, timeInDays);

  // Initialize
  if (printer != nullptr) {
    printer->setFileNamePrefix("");
    printer->openOutput();
    printer->printFileHeader(profile->getAtmosphere());
  }

  // Set the initial seed in the state correlator.
  correlator->setInitialSeed(initialSeed);

  // Step through number of Monte Carlo runs
  for (int run = 0; run < (int)sampleSize; ++run) {
    // Reset the delta time (in case of polar traversal)
    profile->setInputParameters(inputParameters);

    // Set the initial elapsed time
    Position initialPosition = profile->getInitialPosition();
    initialPosition.elapsedTime = deltaSeconds;
    profile->setInitialPosition(initialPosition);
    
    // Initialize random number for each run
    profile->getAtmosphere().setSeed(initialSeed);

    // Run the profile
    profile->generate();

    // Loop through profile
    vector<ProfileData>& runData = profile->getProfile();
    for (size_t step = 0; step < runData.size(); ++step) {
      // Reference the nominal state for this step.
      AtmosphereState& atmosNom = nominalRun[step].atmos;

      // Reference the atmosphere state for this step.
      ProfileData& stepData = runData[step];
      AtmosphereState& atmos = stepData.atmos;
            
      // Correlate the atmosphere data with the nominal data.
      correlator->randomize(step, atmos);
      correlator->updateCoefficients(nominalRun[step].position, stepData.position);
      correlator->correlate(atmosNom, atmos);
    }

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

//! \fn void CorrMonte::setProfile(Profile& prof)
//! \brief Sets the profile used to generate data.
//!
//! \param prof A configured subclass of Profile.

//! \fn void CorrMonte::setProfilePrinter(ProfilePrinter& profPrinter)
//! \brief Sets the profile printer used to output data.
//!
//! \param profPrinter A configured ProfilePrinter.

//! \fn void CorrMonte::setInitialSeed(int seed)
//! \brief Sets the first pertubation seed.
//!
//! If the first perturbation seed is set via this method, then
//! subsequent seeds will be auto-generated.
//! \param seed An integer between 0 and 30000.

//! \fn void CorrMonte::setSampleSize(size_t size)
//! \brief Sets the number of profiles to generate.
//!
//! \param size The number of profiles.

//! \fn void CorrMonte::setCorrDeltaHours(greal deltaHours)
//! \brief Sets the time offset for hourly dispersions.
//!
//! \param deltaHours The time offset \units{\text{hours}}.

} // namespace
