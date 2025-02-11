//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Titan-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"
#include "InputParameters.h"

namespace GRAM {

enum TitanModelType { Yelle97, GCM95 };

//! \brief Titan input parameters.
//!
//! This class contains Titan specific input parameters as well
//! as all common input parameters. See the User Guide for details
//! on each paramter.
//! \ingroup Cpp_Titan TitanGRAM
class TitanInputParameters : public InputParameters
{
public:
  TitanInputParameters();
  TitanInputParameters(const TitanInputParameters& orig) = default;
  virtual ~TitanInputParameters() = default;

  greal minMaxFactor = 0.0;             //!< For TitanAtmosphere (also FMINMAX).
  bool computeMinMaxFactor = true;      //!< For TitanAtmosphere (also IFMM).
  TitanModelType modelType = Yelle97;   //!< For TitanAtmosphere (also IFMM).
  greal userMethaneMoleFraction = 0.0;  //!< For TitanAtmosphere (also fmolmeth).

private:

};

} // namespace

