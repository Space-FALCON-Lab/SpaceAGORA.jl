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
//! Returns a pseudo-random number rectangularly distributed between 0 and 1.
//! Algorithm AS 183 Appl. Statist. (1982) Vol. 31, p.188
//! "Algorithm AS 183: An Efficient and Portable Pseudo-Random Number Generator"
//! by B. A. Wichmann and I. D. Hill.
//!
//! \ingroup EarthGRAM
class WHRandomNumberGenerator
{
public:
  WHRandomNumberGenerator();
  WHRandomNumberGenerator(const WHRandomNumberGenerator& orig) = default;
  virtual ~WHRandomNumberGenerator() = default;

  void setSeed(int seed);
  void getRandomNumbers(std::vector<double>& rvec);

  int getSeed() { return inputSeed; }
  bool initialized() { return (inputSeed != 0); }

//private:
  greal uniformRandomNumber();

  // Random numbers
  int ix = 0;
  int iy = 0;
  int iz = 0;
  int inputSeed = 0;
};

} // namespace
