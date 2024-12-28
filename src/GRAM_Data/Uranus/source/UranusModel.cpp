//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Uranus-GRAM
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "UranusModel.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
UranusModel::UranusModel()
: HeightModel(), UranusCommon(this)
{
  initializeGases();
  initializeData();
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
UranusModel::UranusModel(const UranusModel& orig)
  : HeightModel(orig), UranusCommon(this)
{
  initializeGases();
}

//! \fn  UranusModel::~UranusModel()
//! \copydoc Atmosphere::~Atmosphere()

//! \copydoc Atmosphere::update()
void UranusModel::update()
{
  // Atmosphere state comes from the HeightModel.
  updateAtmosphereState();

  // Update wind components.
  updateWinds();
}

//! \brief Compute zonal and meridional winds at height and latitude.
//!
//! At this point, there is no Uranus winds model.
//! \retval #ewWind       \copybrief ewWind
//! \retval #nsWind       \copybrief nsWind
void UranusModel::updateWinds()
{
  ewWind = 0.0;
  nsWind = 0.0;
}

} // namespace