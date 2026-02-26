//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Earth-GRAM
// Point of Contact: P. White
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include "RandomNumberGenerator.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
RandomNumberGenerator::RandomNumberGenerator()
{
}

//! \fn RandomNumberGenerator::RandomNumberGenerator(const RandomNumberGenerator& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn RandomNumberGenerator::~RandomNumberGenerator()
//! \copydoc Atmosphere::~Atmosphere()

//! \fn RandomNumberGenerator::setSeed()
//! \brief Initializes the random number generator.
//!
//! \param seed This value determines the random number sequence.

//! \fn RandomNumberGenerator::getRandomNumbers()
//! \brief Fills a vector with random numbers.
//!
//! \param[out] rvec A vector of doubles. The size of the vector should be equal to the number of desired random numbers.
//!
//! \retval rvec Each entry in the vector is populated with a random double.

//! \brief Initializes the random number generator.
//!
//! Initializing routine for RCARRY.  Must be called before generating any pseudorandom
//! numbers with RCARRY.  The input parameter should be in the range: 0 <= seed <= 2^24 = 16,777,216.
//! Adapted from RMARIN subroutine of James, Comp. Phys. Comm. 60, 329-344 (1990).
//!
//! \param seed This value determines the random number sequence.
void RandomNumberGenerator::rcarin(int seed)
{
  int ij = seed / 30082;
  int kl = seed - 30082 * ij;
  int i = (ij / 177) % 177 + 2;
  int j = ij % 177 + 2;
  int k = (kl / 169) % 178 + 1;
  int l = kl % 169;
  for (int ii = 0; ii < 24; ii++) {
    double s = 0.0;
    double t = 0.5;
    for (int jj = 0; jj < 24; jj++) {
      int m = (((i * j) % 179) * k) % 179;
      i = j;
      j = k;
      k = m;
      l = (53 * l + 1) % 169;
      if ((l * m) % 64 >= 32) {
        s = s + t;
      }
      t = 0.5 * t;
    }
    mseeds[ii] = s;
  }
  mi24 = 23;
  mj24 = 9;
  mcarry = 0.0;
}

//! \brief Fills a vector with random numbers.
//!
//! Portable Pseudorandom number generator with period of about (1/48)*(2^24)^24 =
//! 2^570 = 10^171.  Author F James, CERN, 1989.  Algorithm due to:  G. Marsgalia
//! and A. Zaman.  Code as in James, Comp. Phys. Comm., 60, 329-344 (1990).
//!
//! \param[out] rvec A vector of doubles. The size of the vector should be equal to the number of desired random numbers.
//!
//! \retval rvec Each entry in the vector is populated with a random double.
void RandomNumberGenerator::rcarry(std::vector<double>& rvec)
{
  for (size_t ivec = 0; ivec < rvec.size(); ivec++) {
    double uni = mseeds[mi24] - mseeds[mj24] - mcarry;
    if (uni < 0.0) {
      uni = uni + 1.0;
      mcarry = twom24;
    }
    else {
      mcarry = 0.0;
    }
    mseeds[mi24] = uni;
    mi24 = mi24 - 1;
    if (mi24 == -1) {
      mi24 = 23;
    }
    mj24 = mj24 - 1;
    if (mj24 == -1) {
      mj24 = 23;
    }
    // Avoid random number of exactly zero (see James, p.344)
    if (uni == 0.0) {
      uni = mseeds[mi24] * twom24;

      if (uni == 0.0) {
        uni = pow(2.0, -48);
      }
    }
    rvec[ivec] = uni;
  }
}

} // namespace
