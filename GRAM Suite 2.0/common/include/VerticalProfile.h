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
#include "Profile.h"
#include "Position.h"

namespace GRAM {

//! \brief Generate a vertical height profile over a specified position.
//!
//! A vertical profile is defined by specifying an initial position, a height delta, and
//! the number of data points desired.  
//! \par Usage
//!     Instantiate an atmosphere for a planet and set its input parameters.  Next, provide 
//!     that object to the VerticalProfile with setAtmosphere().  Set the VerticalProfile's 
//!     initial position and height conditions with setInitialPostion() and setHeightConditions().
//!     Call the generate function to run the model.  Access the results with getProfile().
//! \ingroup CommonGRAM
class VerticalProfile : public Profile
{
public:
  VerticalProfile();
  VerticalProfile(const VerticalProfile& orig);
  virtual ~VerticalProfile() override;

  void setHeightConditions(greal deltaHeight, int numPoints);
  void setInitialPerturbationsUpdated(bool flag) { initialPerturbationsUpdated = flag; }

  void setInputParameters(const InputParameters& params) override;

  void generate() override;

protected:
  greal deltaHeight = 0;    //!< The change in height for each step of the profile.
  int numPoints = 21;       //!< The number of data points to generate in the profile.
  bool initialPerturbationsUpdated = true; //!< If false, initial perturbations will not be updated.
};

} // namespace
