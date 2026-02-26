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
#include <cmath>
#include "gram.h"

namespace GRAM {

/*! \brief A utility class for the interpolation of data.
   
    This utility class performs interpolation of data for up to four dimensions.
    Each interpolator requires an interpolation fraction, f, between 0 and 1.
    For example, linear interpolation between y1 and y2 yields the point y that is
    a fraction f of the distance between y1 and y2.
    \verbatim
          |   f*(y1-y2)  |
 ---------+--------------*-------------------------+-------
          y1      y = y1 + f*(y2-y1)               y2
    \endverbatim
    Multidimensional interpolation will require a fraction for each dimension.
    \verbatim
      z2  +---------------------------+
          |                           |
          |                           |
          |                           |
       z  *--------------*            |
          |              |            |
      f2  |              |            |
          |              |            |
      z1  +--------------*------------+
          y1     f1      y            y2
    \endverbatim 
    \par Usage
       First set the interpolation fraction either using the constuctor, setFraction(), or makeFraction().
       Then invoke linear(), log(), or inverse() interpolation.  All calls to interpolate
       will use the same fraction until a subsequent call to setFraction() is made.
    \code
    Interpolator interp;
    interp.makeFraction(lowHeight, highHeight, height);
    temperature = interp.linear(lowTemp, highTemp);
    pressure = interp.log(lowPressure, highPressure);
    numberDensity = interp.log(lowND, highND);
    \endcode
*/
//! \ingroup CommonGRAM
class Interpolator
{
public:
  Interpolator();
  Interpolator(greal fraction);
  Interpolator(greal fraction1, greal fraction2);
  Interpolator(greal fraction1, greal fraction2, greal fraction3);
  Interpolator(greal fraction1, greal fraction2, greal fraction3, greal fraction4);
  Interpolator(const std::vector<greal>& fractions);
  Interpolator(const Interpolator& orig);
  virtual ~Interpolator();

  void makeFraction(greal start, greal end, greal middle);
  void makeCosineSquaredFraction(greal start, greal end, greal middle);
  void setFraction(greal fraction);
  void setFraction(const std::vector<greal>& fractions);
  greal linear(greal start, greal end);
  greal linear(const greal cube[2]);
  greal linear(const greal cube[2][2]);
  greal linear(const greal cube[2][2][2]);
  greal linear(const greal cube[2][2][2][2]);
  greal linear(greal cube00, greal cube01, greal cube10, greal cube11);
  greal log(greal start, greal end);
  greal log(const greal cube[2]);
  greal log(const greal cube[2][2]);
  greal log(const greal cube[2][2][2]);
  greal log(const greal cube[2][2][2][2]);
  greal log(greal cube00, greal cube01, greal cube10, greal cube11);
  greal inverse(greal start, greal end);

private:
  std::vector<greal> fractions;   //!< The interpolation fractions.
  std::vector<greal> coefficents; //!< The interpolation coefficients.

};

} // namespace
