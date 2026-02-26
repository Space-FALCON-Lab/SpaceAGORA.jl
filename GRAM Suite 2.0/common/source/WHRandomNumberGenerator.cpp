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
#include "WHRandomNumberGenerator.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
WHRandomNumberGenerator::WHRandomNumberGenerator()
{
}

//! \fn WHRandomNumberGenerator::WHRandomNumberGenerator(const WHRandomNumberGenerator& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn WHRandomNumberGenerator::~WHRandomNumberGenerator()
//! \copydoc Atmosphere::~Atmosphere()

//! \brief Initializes the random number generator.
//!
//! Initializing routine which must be called before generating any pseudorandom numbers.
//!
//! \param seed This value determines the random number sequence.
void WHRandomNumberGenerator::setSeed(int seed)
{
  inputSeed = seed % 30000;
  ix = inputSeed;
  iy = 172 * (ix % 176) - 35 * (ix / 176);
  iz = 170 * (ix % 178) - 63 * (ix / 178);
  if (iy < 0) iy = iy + 30307;
  if (iz < 0) iz = iz + 30323;
}

//! \brief Fills a vector with random numbers.
//!
//! Wichmann-Hill
//! Algorithm AS 183 Appl. Statist. (1982) Vol. 31, p.188
//!
//! \param[out] rvec A vector of doubles. The size of the vector should be equal to the number of desired random numbers.
//!
//! \retval rvec Each entry in the vector is populated with a random double.
void WHRandomNumberGenerator::getRandomNumbers(std::vector<double>& rvec)
{
  for (size_t ivec = 0; ivec < rvec.size(); ivec++) {
    rvec[ivec] = uniformRandomNumber();
  }
}

//! \brief A uniform random number generator.
//!
//! Wichmann-Hill
//! Algorithm AS 183 Appl. Statist. (1982) Vol. 31, p.188
//! \returns A pseudo-random number rectangularly distributed between 0 and 1.
greal WHRandomNumberGenerator::uniformRandomNumber()
{
  // IX, IY and IZ should be set to integer values between
  // 1 and 30,000 before first entry.

  // Integer arithmetic up to 5,212,632 is required.
  // Returns L = 0 unless random = 0 or random = 1, in which
  // case L = 1
  greal retval = 0;
  int L = 1;
  while (L == 1) {
    ix = (171 * ix) % 30269;
    iy = (172 * iy) % 30307;
    iz = (170 * iz) % 30323;

    retval = fmod(double(ix) / 30269.0 + double(iy) / 30307.0 + double(iz) / 30323.0, 1.0);
    if (retval > 0.0 && retval < 1.0) {
      L = 0;
    }
  }
  return retval;
}

} // namespace
