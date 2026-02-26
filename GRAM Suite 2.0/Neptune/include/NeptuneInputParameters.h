//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Neptune-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"
#include "InputParameters.h"

namespace GRAM {

//! \brief Neptune input parameters.
//!
//! This class contains Neptune specific input parameters as well
//! as all common input parameters. See the User Guide for details
//! on each paramter.
//! \ingroup Cpp_Neptune NeptuneGRAM
class NeptuneInputParameters : public InputParameters
{
public:
  NeptuneInputParameters();
  NeptuneInputParameters(const NeptuneInputParameters& orig) = default;
  virtual ~NeptuneInputParameters() = default;

  greal minMaxFactor = 0.0;             //!< \brief For NeptuneAtmosphere (also FMINMAX).
  bool computeMinMaxFactor = true;      //!< \brief For NeptuneAtmosphere (also IFMM).
  greal dinitrogenMoleFraction = 0.0;   //!< \brief For NeptuneAtmosphere (also FMOLNITRO).
private:

};

} // namespace

