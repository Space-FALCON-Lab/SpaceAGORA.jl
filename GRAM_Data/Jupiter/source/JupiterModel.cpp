//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Jupiter-GRAM
//////////////////////////////////////////////////////////////////////////

#include "JupiterModel.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
JupiterModel::JupiterModel()
: HeightModel(), JupiterCommon(this)
{
  initializeGases();
  initializeData();
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
JupiterModel::JupiterModel(const JupiterModel& orig)
  : HeightModel(orig), JupiterCommon(this)
{
  initializeGases();
}

//! \fn  JupiterModel::~JupiterModel()
//! \copydoc Atmosphere::~Atmosphere()

//! \copydoc Atmosphere::update()
void JupiterModel::update()
{
  // Atmosphere state comes from the HeightModel.
  updateAtmosphereState();

  // Update wind components.
  updateWinds();
}

//! \brief Compute zonal and meridional winds at height and latitude.
//!
//! At this point, there is no Jupiter winds model.
//! \retval #ewWind       \copybrief ewWind
//! \retval #nsWind       \copybrief nsWind
void JupiterModel::updateWinds()
{
  ewWind = 0.0;
  nsWind = 0.0;
}

} // namespace