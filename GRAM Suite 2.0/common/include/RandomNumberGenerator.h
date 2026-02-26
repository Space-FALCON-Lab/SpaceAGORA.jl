//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Earth-GRAM
// Point of Contact: P. White
//////////////////////////////////////////////////////////////////////////

#pragma once

#include<vector>
#include "gram.h"

namespace GRAM {

//! \brief This class generates random numbers.
//!
//! A portable pseudorandom number generator with period of about (1/48)*(2^24)^24 =
//! 2^570 = 10^171.  Author F James, CERN, 1989.  Algorithm due to:  G. Marsgalia
//! and A. Zaman.  Code as in James, Comp. Phys. Comm., 60, 329-344 (1990).
//!
//! \ingroup EarthGRAM
class RandomNumberGenerator
{
public:
  RandomNumberGenerator();
  RandomNumberGenerator(const RandomNumberGenerator& orig) = default;
  virtual ~RandomNumberGenerator() = default;

  void setSeed(int seed) { rcarin(seed); inputSeed = seed; }
  void getRandomNumbers(std::vector<double>& rvec) { rcarry(rvec); }

  int getSeed() const { return inputSeed; }
  bool initialized() const { return (inputSeed != 0); }
  void reset() { inputSeed = 0; }

private:
  void rcarin(int ijkl);
  void rcarry(std::vector<double>& rvec);

  // Random numbers
  double mseeds[24] = { 0.0 };
  double mcarry = 0.0;
  int mi24 = 0;
  int mj24 = 0;
  int inputSeed = 0;

  static constexpr double twom24 = 1.0 / 16777216.0;
};

} // namespace
