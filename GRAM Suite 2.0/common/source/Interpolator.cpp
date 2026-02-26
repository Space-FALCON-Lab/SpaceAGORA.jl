//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <cmath>
#include "Interpolator.h"

using namespace std;

namespace GRAM {

//! \brief The default constructor
//!
//! This basic constructor is called when objects are created without parameters.
//! In the default constructor the interpolation fraction defaults to 0.
Interpolator::Interpolator()
{
}

//! \brief The copy constructor
//!
//! This constructor enables the copying of objects.  This can be done by using
//! the constructor: Object newobject(oldobject).  Or by assignment: newobject = oldobject.
//! \param orig The source object to copy.
Interpolator::Interpolator(const Interpolator& orig)
{
  fractions = orig.fractions;
  coefficents = orig.coefficents;
}

//! \brief The destructor.
//!
//! The destructor is called when an object is deleted or goes out of scope.
Interpolator::~Interpolator()
{
  fractions.clear();
  coefficents.clear();
}

//! \brief This constructor sets the interpolation fraction.
//!
//! \param fraction This should be a value between -1 and 1 (inclusive).
Interpolator::Interpolator(greal fraction)
{
  setFraction(fraction);
}

//! \brief This constructor sets the interpolation fraction.
//!
//! \param fraction1 This should be a value between -1 and 1 (inclusive).
//! \param fraction2 This should be a value between -1 and 1 (inclusive).
Interpolator::Interpolator(greal fraction1, greal fraction2)
{
  setFraction(vector<greal> {fraction1, fraction2});
}

//! \brief This constructor sets the interpolation fraction.
//!
//! \param fraction1 This should be a value between -1 and 1 (inclusive).
//! \param fraction2 This should be a value between -1 and 1 (inclusive).
//! \param fraction3 This should be a value between -1 and 1 (inclusive).
Interpolator::Interpolator(greal fraction1, greal fraction2, greal fraction3)
{
  setFraction(vector<greal> {fraction1, fraction2, fraction3});
}

//! \brief This constructor sets the interpolation fraction.
//!
//! \param fraction1 This should be a value between -1 and 1 (inclusive).
//! \param fraction2 This should be a value between -1 and 1 (inclusive).
//! \param fraction3 This should be a value between -1 and 1 (inclusive).
//! \param fraction4 This should be a value between -1 and 1 (inclusive).
Interpolator::Interpolator(greal fraction1, greal fraction2, greal fraction3, greal fraction4)
{
  setFraction(vector<greal> {fraction1, fraction2, fraction3, fraction4});
}

//! \brief This constructor sets the interpolation fraction for multiple dimensions.
//!
//! \param fractions This vector should be populated a value between -1 and 1 (inclusive) for each dimension.
Interpolator::Interpolator(const std::vector<greal>& fractions)
{
  setFraction(fractions);
}

//! \brief Computes and sets the interpolation fraction based on three data points.
//!
//! The fraction is set according to \f$ fraction = \frac{middle - end}{start - end} \f$.
//! \param start The beginning value of the interpolation interval.
//! \param end The ending value of the interpolation interval.
//! \param middle This should be a value between start and end (inclusive).
void Interpolator::makeFraction(greal start, greal end, greal middle)
{
  // Assumptions (checking occurs in debug mode only)
  assert(start != end);

  setFraction((middle - start) / (end - start));
}

//! \brief Computes and sets the interpolation fraction based on three data points.
//!
//! The fraction is set according to \f$ fraction = \frac{middle - end}{start - end} \f$.
//! \param start The beginning value of the interpolation interval.
//! \param end The ending value of the interpolation interval.
//! \param middle This should be a value between start and end (inclusive).
void Interpolator::makeCosineSquaredFraction(greal start, greal end, greal middle)
{
  // Assumptions (checking occurs in debug mode only)
  assert(start != end);

  fractions.resize(1);
  fractions[0] = (middle - start) / (end - start);
  coefficents.resize(2);
  coefficents[0] = pow(cos(HALF_PI * fractions[0]), 2);    // cos^2
  coefficents[1] = 1.0 - coefficents[0];                   // sin^2
}

//! \brief Sets the interpolation fraction.
//!
//! \param fraction This should be a value between -1 and 1 (inclusive).
void Interpolator::setFraction(greal fraction) 
{ 
  fractions.resize(1);
  fractions[0] = fraction;
  coefficents.resize(2);
  coefficents[0] = 1.0 - fraction;
  coefficents[1] = fraction;
}

//! \brief Sets the interpolation fraction for multiple dimensions.
//!
//! \param fractions This vector should be populated with a value between -1 and 1 (inclusive) for each dimension.
void Interpolator::setFraction(const std::vector<greal>& fractions)
{ 
  this->fractions = fractions;
  size_t dimensions = fractions.size();
  size_t power = (size_t(1) << dimensions);
  coefficents.resize(power);
  for (size_t i = 0; i < power; ++i) {
    coefficents[i] = 1.0;
    for (size_t dim = 0; dim < dimensions; ++dim) {
      if (i & (size_t(1) << (dimensions - dim - 1))) {
        coefficents[i] *= fractions[dim];
      }
      else {
        coefficents[i] *= 1.0 - fractions[dim];
      }
    }
  }
}

//! \brief Performs linear interpolation.
//!
//! Performs linear interpolation: \f[ result = (1-fraction) start + (fraction) end \f]
//! \param start The beginning value of the y interpolation interval.
//! \param end The ending value of the y interpolation interval.
//! \returns The interpolated y value.
greal Interpolator::linear(greal start, greal end)
{
  // Assumptions (checking occurs in debug mode only)
  assert(fractions.size() == 1);

  return coefficents[0] * start + coefficents[1] * end;
}

//! \brief Performs linear interpolation.
//!
//! Performs linear interpolation: \f[ result = (1-fraction) start + (fraction) end \f]
//! \param cube The start (cube[0]) and end (cube[1]) values of the y interpolation interval.
//! \returns The interpolated y value.
greal Interpolator::linear(const greal cube[2])
{
  // Assumptions (checking occurs in debug mode only)
  assert(fractions.size() == 1);

  return coefficents[0] * cube[0]
       + coefficents[1] * cube[1];
  }

//! \brief Performs 2-dimensional linear interpolation.
//!
//! \param cube The corner values of the y interpolation cube.
//! \returns The interpolated y value.
greal Interpolator::linear(const greal cube[2][2])
{
  // Assumptions (checking occurs in debug mode only)
  assert(fractions.size() == 2);

  return coefficents[0] * cube[0][0]
    + coefficents[1] * cube[0][1]
    + coefficents[2] * cube[1][0]
    + coefficents[3] * cube[1][1];
}

//! \brief Performs 2-dimensional linear interpolation.
//!
//! \param cube00 The corner values of the y interpolation cube.
//! \param cube01 The corner values of the y interpolation cube.
//! \param cube10 The corner values of the y interpolation cube.
//! \param cube11 The corner values of the y interpolation cube.
//! \returns The interpolated y value.
greal Interpolator::linear(greal cube00, greal cube01, greal cube10, greal cube11)
{
  // Assumptions (checking occurs in debug mode only)
  assert(fractions.size() == 2);

  return coefficents[0] * cube00
       + coefficents[1] * cube01
       + coefficents[2] * cube10
       + coefficents[3] * cube11;
}

//! \brief Performs 3-dimensional linear interpolation.
//!
//! \param cube The corner values of the y interpolation cube.
//! \returns The interpolated y value.
greal Interpolator::linear(const greal cube[2][2][2])
{
  // Assumptions (checking occurs in debug mode only)
  assert(fractions.size() == 3);

  return coefficents[0] * cube[0][0][0]
       + coefficents[1] * cube[0][0][1]
       + coefficents[2] * cube[0][1][0]
       + coefficents[3] * cube[0][1][1]
       + coefficents[4] * cube[1][0][0]
       + coefficents[5] * cube[1][0][1]
       + coefficents[6] * cube[1][1][0]
       + coefficents[7] * cube[1][1][1];
}

//! \brief Performs 4-dimensional linear interpolation.
//!
//! \param cube The corner values of the y interpolation cube.
//! \returns The interpolated y value.
greal Interpolator::linear(const greal cube[2][2][2][2])
{
  // Assumptions (checking occurs in debug mode only)
  assert(fractions.size() == 4);

  return coefficents[ 0] * cube[0][0][0][0]
       + coefficents[ 1] * cube[0][0][0][1]
       + coefficents[ 2] * cube[0][0][1][0]
       + coefficents[ 3] * cube[0][0][1][1]
       + coefficents[ 4] * cube[0][1][0][0]
       + coefficents[ 5] * cube[0][1][0][1]
       + coefficents[ 6] * cube[0][1][1][0]
       + coefficents[ 7] * cube[0][1][1][1]
       + coefficents[ 8] * cube[1][0][0][0]
       + coefficents[ 9] * cube[1][0][0][1]
       + coefficents[10] * cube[1][0][1][0]
       + coefficents[11] * cube[1][0][1][1]
       + coefficents[12] * cube[1][1][0][0]
       + coefficents[13] * cube[1][1][0][1]
       + coefficents[14] * cube[1][1][1][0]
       + coefficents[15] * cube[1][1][1][1];
}

//! \brief Performs interpolation of logs.
//!
//! Performs interpolation of logs: \f[ result = e^{[(1-fraction)\ln(start) + (fraction) \ln(end)]} \f]
//! \param start The beginning value of the y interpolation interval.
//! \param end The ending value of the y interpolation interval.
//! \returns The interpolated y value.
greal Interpolator::log(greal start, greal end)
{
  if (start <= 0.0 || end <= 0.0) {
    if (end > 0.0) {
      start = 0.1 * end;
    }
    else if (start > 0.0) {
      end = 0.1 * start;
    }
    else {
      return 0.0;
    }
  }
  return std::exp(coefficents[0] * std::log(start) + coefficents[1] * std::log(end));
}

//! \brief Performs interpolation of logs.
//!
//! Performs interpolation of logs: \f[ result = e^{(1-fraction)\ln(start) + (fraction) \ln(end)} \f]
//! \param cube The start (cube[0]) and end (cube[1]) values of the y interpolation interval.
//! \returns The interpolated y value.
greal Interpolator::log(const greal cube[2])
{
  // Assumptions (checking occurs in debug mode only)
  assert(fractions.size() == 1);
  assert(cube[0] > 0.0);
  assert(cube[1] > 0.0);

  return exp(coefficents[0] * std::log(cube[0])
           + coefficents[1] * std::log(cube[1]));
  }

//! \brief Performs 2-dimensional interpolation of logs.
//!
//! \param cube The corner values of the y interpolation cube.
//! \returns The interpolated y value.
greal Interpolator::log(const greal cube[2][2])
{
  // Assumptions (checking occurs in debug mode only)
  assert(fractions.size() == 2);
  assert(cube[0][0] > 0.0);
  assert(cube[0][1] > 0.0);
  assert(cube[1][0] > 0.0);
  assert(cube[1][1] > 0.0);

  return exp(coefficents[0] * std::log(cube[0][0])
           + coefficents[1] * std::log(cube[0][1])
           + coefficents[2] * std::log(cube[1][0])
           + coefficents[3] * std::log(cube[1][1]));
}

//! \brief Performs 2-dimensional interpolation of logs.
//!
//! \param cube00 The corner values of the y interpolation cube.
//! \param cube01 The corner values of the y interpolation cube.
//! \param cube10 The corner values of the y interpolation cube.
//! \param cube11 The corner values of the y interpolation cube.
//! \returns The interpolated y value.
greal Interpolator::log(greal cube00, greal cube01, greal cube10, greal cube11)
{
  // Assumptions (checking occurs in debug mode only)
  assert(fractions.size() == 2);
  assert(cube00 > 0.0);
  assert(cube01 > 0.0);
  assert(cube10 > 0.0);
  assert(cube11 > 0.0);

  return exp(coefficents[0] * std::log(cube00)
           + coefficents[1] * std::log(cube01)
           + coefficents[2] * std::log(cube10)
           + coefficents[3] * std::log(cube11));
}

//! \brief Performs 3-dimensional interpolation of logs.
//!
//! \param cube The corner values of the y interpolation cube.
//! \returns The interpolated y value.
greal Interpolator::log(const greal cube[2][2][2])
{
  // Assumptions (checking occurs in debug mode only)
  assert(fractions.size() == 3);
  assert(cube[0][0][0] > 0.0);
  assert(cube[0][0][1] > 0.0);
  assert(cube[0][1][0] > 0.0);
  assert(cube[0][1][1] > 0.0);
  assert(cube[1][0][0] > 0.0);
  assert(cube[1][0][1] > 0.0);
  assert(cube[1][1][0] > 0.0);
  assert(cube[1][1][1] > 0.0);

  return exp(coefficents[0] * std::log(cube[0][0][0])
           + coefficents[1] * std::log(cube[0][0][1])
           + coefficents[2] * std::log(cube[0][1][0])
           + coefficents[3] * std::log(cube[0][1][1])
           + coefficents[4] * std::log(cube[1][0][0])
           + coefficents[5] * std::log(cube[1][0][1])
           + coefficents[6] * std::log(cube[1][1][0])
           + coefficents[7] * std::log(cube[1][1][1]));
}

//! \brief Performs 4-dimensional interpolation of logs.
//!
//! \param cube The corner values of the y interpolation cube.
//! \returns The interpolated y value.
greal Interpolator::log(const greal cube[2][2][2][2])
{
  // Assumptions (checking occurs in debug mode only)
  assert(fractions.size() == 4);
  assert(cube[0][0][0][0] > 0.0);
  assert(cube[0][0][0][1] > 0.0);
  assert(cube[0][0][1][0] > 0.0);
  assert(cube[0][0][1][1] > 0.0);
  assert(cube[0][1][0][0] > 0.0);
  assert(cube[0][1][0][1] > 0.0);
  assert(cube[0][1][1][0] > 0.0);
  assert(cube[0][1][1][1] > 0.0);
  assert(cube[1][0][0][0] > 0.0);
  assert(cube[1][0][0][1] > 0.0);
  assert(cube[1][0][1][0] > 0.0);
  assert(cube[1][0][1][1] > 0.0);
  assert(cube[1][1][0][0] > 0.0);
  assert(cube[1][1][0][1] > 0.0);
  assert(cube[1][1][1][0] > 0.0);
  assert(cube[1][1][1][1] > 0.0);

  return exp(coefficents[ 0] * std::log(cube[0][0][0][0])
           + coefficents[ 1] * std::log(cube[0][0][0][1])
           + coefficents[ 2] * std::log(cube[0][0][1][0])
           + coefficents[ 3] * std::log(cube[0][0][1][1])
           + coefficents[ 4] * std::log(cube[0][1][0][0])
           + coefficents[ 5] * std::log(cube[0][1][0][1])
           + coefficents[ 6] * std::log(cube[0][1][1][0])
           + coefficents[ 7] * std::log(cube[0][1][1][1])
           + coefficents[ 8] * std::log(cube[1][0][0][0])
           + coefficents[ 9] * std::log(cube[1][0][0][1])
           + coefficents[10] * std::log(cube[1][0][1][0])
           + coefficents[11] * std::log(cube[1][0][1][1])
           + coefficents[12] * std::log(cube[1][1][0][0])
           + coefficents[13] * std::log(cube[1][1][0][1])
           + coefficents[14] * std::log(cube[1][1][1][0])
           + coefficents[15] * std::log(cube[1][1][1][1]));
}

//! \brief Performs interpolation of inverses.
//!
//! Performs interpolation of inverses: \f[ result = \frac{1}{\frac{1-fraction}{start} + \frac{fraction}{end}} \f]
//! \param start The beginning value of the y interpolation interval.
//! \param end The ending value of the y interpolation interval.
//! \returns The interpolated y value.
greal Interpolator::inverse(greal start, greal end)
{
  // Assumptions (checking occurs in debug mode only)
  assert(fractions.size() == 1);
  assert(start != 0.0);
  assert(end != 0.0);

  return 1.0 / (coefficents[0] / start + coefficents[1] / end);
}

} // namespace
