//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: EarthGRAM
//////////////////////////////////////////////////////////////////////////

#include "StateCorrelator.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
StateCorrelator::StateCorrelator()
{
}

//! \fn  StateCorrelator::StateCorrelator(const StateCorrelator& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  StateCorrelator::~StateCorrelator()
//! \copydoc Atmosphere::~Atmosphere()

//! \fn void StateCorrelator::setInitialSeed(int seed)
//! \brief Sets the randomization seed.
//!
//! If the state correlator is also used to introduce greater randomization in the target
//! perturbations, then the random number seed can be set via this method.
//! \param seed The random number seed.

//! \fn void StateCorrelator::randomize(size_t step, AtmosphereState& target)
//! \brief Randomize the target perturbations. 
//!
//! Override this method if further randomization of the target perturbations is required.
//! Typically, a random number generator will be needed and must be provided in the sub-class.
//! Since randomization may need to match those of the planetary model, the step number of the
//! profile/trajectory must be passed as the first argument.
//! If no randomization is required, then implement an empty method.
//! \param step    The step number of the profile.
//! \param target  The target atmosphere state.
//! \returns The modified target atmosphere state.

//! \fn void StateCorrelator::updateCoefficients(const Position& base, const Position& target)
//! \brief Updates correlation coefficients based on position. 
//!
//! The correlation of perturbations is typically in the form of
//! \f$ newTargetPert = b * basePert + c * targetPert \f$
//! where the coefficents \p b and \p c range from 0 to 1 and are complementary.  This method is used
//! to compute the base and target coefficients based on the base and target positions.  The coefficients
//! should be computed so that the base coefficent is near 1 whenever the target position in near the base
//! position.  And the base coefficent should be near 0 whenever the target position in far from the base
//! position.
//! This method must be overriden by the sub-class.
//! \param base    The position of the base atmosphere state.
//! \param target  The position of the target atmosphere state.

//! \fn void StateCorrelator::correlate(const AtmosphereState& base, AtmosphereState& target)
//! \brief Correlates the target perturbations with the base perturbations.  
//!
//! Override this method to perform the correlations of the target perturbations using the coefficients
//! produced in updateCoefficients().  The correlation of perturbations is typically in the form of 
//! \f$ newTargetPert = b * basePert + c * targetPert \f$.  This method should perform the correlations
//! for all applicable perturbations in the target state.  It will also need to update the total
//! perturbed value for each of the modified metrics.
//! \param base    The base atmosphere state.
//! \param target  The target atmosphere state. 
//! \returns   The correlated target atmosphere state.

} // namespace
