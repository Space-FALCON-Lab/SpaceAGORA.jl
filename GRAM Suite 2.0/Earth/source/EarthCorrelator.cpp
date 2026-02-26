//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include "EarthCorrelator.h"
#include "EarthAtmosphereState.h"
#include "EarthAtmosphere.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
EarthCorrelator::EarthCorrelator()
  : StateCorrelator()
{
}

//! \fn  EarthCorrelator::EarthCorrelator(const EarthCorrelator& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  EarthCorrelator::~EarthCorrelator()
//! \copydoc Atmosphere::~Atmosphere()

//! \fn  EarthCorrelator::setMean(bool flag)
//! \brief Turn on correlation using means.
//!
//! When the mean flag is set, the correlation of the target will be with the unperturbed means.
//! A base position and atmosphere state must still be provided.
//! \param flag  Set to true to enable correlation with the means.

//! \brief Randomize the target perturbations. 
//!
//! Only the first step of each trajectory is randomized in EarthCorrMonte.  The random number 
//! generator is initialized once on the first trajectory. Si the first
//! trajectory will randomize to the same sequence at the EarthAtmosphere.  But the same generator is
//! used by the successive trajectories without reinitialize the sequence.  So they will disperse the
//! small scale perturbations in the first step.
//!
//! \param step    The step number of the profile.
//! \param target  The target atmosphere state.
//! \returns The modified target atmosphere state.
void EarthCorrelator::randomize(size_t step, AtmosphereState& target)
{
  // Only randomize the first step in each trajectory.
  if (step > 0) {
    return;
  }

  // Create a vector for the random numbers.
  std::vector<double> randomNumbers(6);

  // On the first pass, initialize the random number generator.
  // This random number generator will be synchronized with the generator in EarthAtmosphere.
  if (!randomNumberGenerator.initialized()) {
    // Set the initial seed.
    randomNumberGenerator.setSeed(initialSeed);

    // Get 10 random numbers. This corresponds to the 10 random numbers generated 
    // by a call to EarthAtmosphere::initializePerturbations().
    randomNumbers.resize(10);
    randomNumberGenerator.getRandomNumbers(randomNumbers);

    // Resize back to 6.
    randomNumbers.resize(6);
  }

  // Get six random numbers for randomizing small scale perts
  randomNumberGenerator.getRandomNumbers(randomNumbers);

  // Randomize the small scale perturbations.
  EarthAtmosphereState& earth = target.getPlanetSpecificMetrics<EarthAtmosphereState>(); 
  const greal sr3 = sqrt(3);
  earth.presPertSmall *= sr3 * (2.0 * randomNumbers[0] - 1.0);
  earth.densPertSmall *= sr3 * (2.0 * randomNumbers[1] - 1.0);
  earth.tempPertSmall *= sr3 * (2.0 * randomNumbers[2] - 1.0);
  earth.ewWindPertSmall *= sr3 * (2.0 * randomNumbers[3] - 1.0);
  earth.nsWindPertSmall *= sr3 * (2.0 * randomNumbers[4] - 1.0);
  target.verticalWindPerturbation *= sr3 * (2.0 * randomNumbers[5] - 1.0);
}

//! \brief Updates correlation coefficients based on position. 
//!
//! The base coefficient uses horizontal distance, vertical distance, and time.  The incoming target position
//! must have a valid totalRadius in order to compute the horizontal distance.  Horizontal and vertical distance
//! scales are acquired from the EarthAtmosphere scales used for small scale perturbations.
//!
//! \param base    The position of the base atmosphere state.
//! \param target  The position of the target atmosphere state.
void EarthCorrelator::updateCoefficients(const Position& base, const Position& target)
{
  // Since this computation requires data in EarthAtmosphere, it was turned into a static EarthAtmosphere method.
  EarthAtmosphere::getCorrelationCoefficients(base, target, baseCoefficient, targetCoefficient);
}

//! \brief Correlates the target perturbations with the base perturbations.  
//!
//! This method correlates the small scale perturbations of the pressure, density, temperature, and winds.  It
//! then recomputes the total of the mean, large scale peturbations, and small scale peturbations.
//!
//! \param base    The base atmosphere state.
//! \param target  The target atmosphere state. 
//! \returns   The correlated target atmosphere state.
void EarthCorrelator::correlate(const AtmosphereState& base, AtmosphereState& target)
{
  // Get EarthAtmosphereStates from the base and target.  These contain the small scale perts.
  const EarthAtmosphereState& earthBase = base.getPlanetSpecificMetrics<EarthAtmosphereState>();
  EarthAtmosphereState& earthTarget = target.getPlanetSpecificMetrics<EarthAtmosphereState>();

  // If the mean flag is set, then correlate to the mean values.  That is, assume the base has no perturbations.
  if (corrMean) {
    // Correlate the small scale perts.
    earthTarget.presPertSmall *= targetCoefficient;
    earthTarget.densPertSmall *= targetCoefficient;
    earthTarget.tempPertSmall *= targetCoefficient;
    earthTarget.ewWindPertSmall *= targetCoefficient;
    earthTarget.nsWindPertSmall *= targetCoefficient;
    target.verticalWindPerturbation *= targetCoefficient;

    // Recompute the total perturbed values.
    earthTarget.perturbedPressure    = base.pressure     + earthTarget.presPertSmall * target.pressure;
    target.perturbedDensity          = base.density      + earthTarget.densPertSmall * target.density;
    earthTarget.perturbedTemperature = base.temperature  + earthTarget.tempPertSmall * target.temperature;
    target.perturbedEWWind           = base.ewWind       + earthTarget.ewWindPertSmall;
    target.perturbedNSWind           = base.nsWind       + earthTarget.nsWindPertSmall;
    target.perturbedVerticalWind     = base.verticalWind + target.verticalWindPerturbation;
  }
  // Otherwise perform correlation of the target with the base.
  else {
    // Compute large perts for winds in case they were limited by SOS.
    // Note that the existing LS perts won't contain that limitation.  It is in the perturbed winds.
    greal ewWindPertLarge = target.perturbedEWWind - target.ewWind - earthTarget.ewWindPertSmall;
    greal nsWindPertLarge = target.perturbedNSWind - target.nsWind - earthTarget.nsWindPertSmall;

    // Correlate the small scale perts.
    earthTarget.presPertSmall  = (baseCoefficient * earthBase.presPertSmall * base.pressure
                                + targetCoefficient * earthTarget.presPertSmall * target.pressure) / target.pressure;
    earthTarget.densPertSmall = (baseCoefficient * earthBase.densPertSmall * base.density
                               + targetCoefficient * earthTarget.densPertSmall * target.density) / target.density;
    earthTarget.tempPertSmall = (baseCoefficient * earthBase.tempPertSmall * base.temperature
                               + targetCoefficient * earthTarget.tempPertSmall * target.temperature) / target.temperature;
    earthTarget.ewWindPertSmall = baseCoefficient * earthBase.ewWindPertSmall
                                + targetCoefficient * earthTarget.ewWindPertSmall;
    earthTarget.nsWindPertSmall = baseCoefficient * earthBase.nsWindPertSmall
                                + targetCoefficient * earthTarget.nsWindPertSmall;
    target.verticalWindPerturbation = baseCoefficient * base.verticalWindPerturbation
                                    + targetCoefficient * target.verticalWindPerturbation;

    // Recompute the total perturbed values.
    earthTarget.perturbedPressure = target.pressure + earthTarget.presPertLarge * target.pressure
                                                    + earthTarget.presPertSmall * target.pressure;
    target.perturbedDensity = target.density + earthTarget.densPertLarge * target.density
                                             +  earthTarget.densPertSmall * target.density;
    earthTarget.perturbedTemperature = target.temperature + earthTarget.tempPertLarge * target.temperature
                                                          +  earthTarget.tempPertSmall * target.temperature;
    target.perturbedEWWind = target.ewWind + ewWindPertLarge + earthTarget.ewWindPertSmall;
    target.perturbedNSWind = target.nsWind + nsWindPertLarge + earthTarget.nsWindPertSmall;
    target.perturbedVerticalWind = target.verticalWind + target.verticalWindPerturbation;
  }
}

} // namespace
