//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Jupiter-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"
#include "InputParameters.h"

namespace GRAM {

//! \brief Jupiter input parameters.
//!
//! This class contains Jupiter specific input parameters as well
//! as all common input parameters. At this point, there are no
//! Jupiter specific parameters.
//! \ingroup Cpp_Jupiter JupiterGRAM
class JupiterInputParameters : public InputParameters
{
public:
  JupiterInputParameters();
  JupiterInputParameters(const JupiterInputParameters& orig) = default;
  virtual ~JupiterInputParameters() = default;
  JupiterInputParameters& operator=(const JupiterInputParameters& orig) = default;

private:

};

} // namespace

