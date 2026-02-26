//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"

namespace GRAM {

//! \brief Contains the randomized values of a perturbation.
//!
//! Use this structure to persist the random state of a perturbation.
//! \ingroup CommonGRAM Cpp_Venus Cpp_Earth Cpp_Mars Cpp_Jupiter Cpp_Saturn Cpp_Uranus Cpp_Neptune Cpp_Titan
class PerturbationState
{
public:
  PerturbationState();
  PerturbationState(const PerturbationState& orig) = default;
  virtual ~PerturbationState() = default;

  greal densityRandomNumber = 0.0;      //!< \brief A uniform random number in [0,1]
  greal ewWindRandomNumber = 0.0;       //!< \brief A uniform random number in [0,1]
  greal nsWindRandomNumber = 0.0;       //!< \brief A uniform random number in [0,1]
  greal verticalWindRandomNumber = 0.0; //!< \brief A uniform random number in [0,1]
  greal densityRho = 0.0;               //!< \brief The previous density rho.
  greal ewWindRho = 0.0;                //!< \brief The previous EW rho.
  greal nsWindRho = 0.0;                //!< \brief The previous NS rho.
  greal verticalWindRho = 0.0;          //!< \brief The previous vertical rho.
};

} // namespace
