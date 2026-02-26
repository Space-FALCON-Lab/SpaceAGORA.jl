//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include "Profile.h"

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
Profile::Profile()
{
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
Profile::Profile(const Profile& orig)
{
  atmosphere = orig.atmosphere;
  profile = orig.profile;
}

//! \copydoc Atmosphere::~Atmosphere()
Profile::~Profile()
{
  profile.clear();
}

//! \fn  Profile::~Profile()
//! \copydoc Atmosphere::~Atmosphere()

//! \fn void Profile::setAtmosphere(PerturbedAtmosphere& atmos)
//! \brief Sets the atmosphere model used to generate data.
//!
//! \param atmos A PerturbedAtmosphere subclass.

//! \fn Profile::getAtmosphere()
//! \brief Gets the atmosphere model used to generate data.
//!
//! \returns A PerturbedAtmosphere subclass.

//! \fn void Profile::setInputParameters(const InputParameters& params)
//! \brief Set the applicable input parameters.
//!
//! The routine copies the applicable members of the InputParameters
//! object into the appropriate members of this object.
//! \param params The input parameters.

//! \fn void Profile::setInitialPosition(const Position& p)
//! \brief Sets the starting position of the trajectory.
//!
//! \param p A Position object.

//! \fn void Profile::generate()
//! \brief Generates the profile data.
//!
//! This pure virtual method should control the generation of the profile data.
//! Typically, the override will loop over a number of steps, provide a position
//! to the atmosphere, and collect the atmospheric data.

//! \fn Profile::getProfile()
//! \brief Gets the generated list of ProfileData.
//!
//! \returns A std::vector of ProfileData.


} // namespace
