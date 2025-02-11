//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"
#include "InputParameters.h"

namespace GRAM {

//! \brief VenusInputParameters.
//!
//! VenusInputParameters.
//! \ingroup Cpp_Venus VenusGRAM
class VenusInputParameters : public InputParameters
{
public:
  VenusInputParameters();
  VenusInputParameters(const VenusInputParameters& orig) = default;
  virtual ~VenusInputParameters() = default;


private:

};

} // namespace

