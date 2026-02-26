//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"
#include "StateCorrelator.h"
#include "RandomNumberGenerator.h"

namespace GRAM {

//! \brief Correlate perturbations of Earth atmosphere states.
//!
//! State correlators are used to correlate the perturbations of a target atmosphere state 
//! with those of a base atmosphere state.  The correlation is typically based on space-time 
//! distance. So if the target position is near the base position, then the correlated 
//! target perturbations will be similar to the base perturbations.  The state correlator
//! also allows the plantary models to optionally introduce randomization of the perturbations.
//! The EarthCorrelator implements \p randomize, \p updateCoefficients, and \p correlate.
//! \ingroup EarthGRAM
class EarthCorrelator : public StateCorrelator
{
//! \example Earth/examples/EarthCorrMulti.cpp
public:
  EarthCorrelator();
  EarthCorrelator(const EarthCorrelator& orig) = default;
  virtual ~EarthCorrelator() = default;

  void setMean(bool flag) { corrMean = flag; }

  void randomize(size_t step, AtmosphereState& target) override;
  void updateCoefficients(const Position& base, const Position& target) override;
  void correlate(const AtmosphereState& base, AtmosphereState& target) override;

protected:
  bool corrMean = false;
  greal baseCoefficient = 0.0;
  greal targetCoefficient = 0.0;

  RandomNumberGenerator randomNumberGenerator;
  int runSeed = 0;
};

} // namespace
