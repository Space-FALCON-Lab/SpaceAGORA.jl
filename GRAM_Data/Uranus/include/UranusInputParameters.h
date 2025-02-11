//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Uranus-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"
#include "InputParameters.h"

namespace GRAM {

enum UranusModelType { Unew, Uold };

//! \brief Uranus input parameters.
//!
//! This class contains Uranus specific input parameters as well
//! as all common input parameters. At this point, there are no
//! Uranus specific parameters.
//! \ingroup Cpp_Uranus UranusGRAM
class UranusInputParameters : public InputParameters
{
public:
  UranusInputParameters();
  UranusInputParameters(const UranusInputParameters& orig) = default;
  virtual ~UranusInputParameters() = default;
  UranusInputParameters& operator=(const UranusInputParameters& orig) = default;

private:

};

} // namespace

