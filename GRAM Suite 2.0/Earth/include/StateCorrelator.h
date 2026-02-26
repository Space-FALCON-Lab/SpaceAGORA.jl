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
#include "AtmosphereState.h"
#include "Position.h"

namespace GRAM {

//! \brief Correlate perturbations of atmosphere states.
//!
//! State correlators are used to correlate the perturbations of a target atmosphere state 
//! with those of a base atmosphere state.  The correlation is typically based on space-time 
//! distance. So if the target position is near the base position, then the correlated 
//! target perturbations will be similar to the base perturbations.  The state correlator
//! also allows the plantary models to optionally introduce randomization of the perturbations.
//! The StateCorrelator is an abstract class that defines the interface
//! for planet specific state correlators.  The subclasses must implement 
//! three member functions: \p randomize, \p updateCoefficients, and \p correlate.
//! \ingroup EarthGRAM
class StateCorrelator
{
public:
  StateCorrelator();
  StateCorrelator(const StateCorrelator& orig) = default;
  virtual ~StateCorrelator() = default;

  virtual void setInitialSeed(int seed) { initialSeed = seed; }
  virtual void randomize(size_t step, AtmosphereState& target) = 0;

  virtual void updateCoefficients(const Position& base, const Position& target) = 0;
  virtual void correlate(const AtmosphereState& base, AtmosphereState& target) = 0;

protected:
  int initialSeed = 0;   //!< The seed for random number generation.

};

} // namespace
