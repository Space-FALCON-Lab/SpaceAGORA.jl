//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include "HeightModel.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
HeightModel::HeightModel()
  : Atmosphere()
{
}

//! \fn HeightModel::HeightModel(const HeightModel& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  HeightModel::~HeightModel()
//! \copydoc Atmosphere::~Atmosphere()

//! \brief Adds all gases into a map for easy processing.
//!
//! This routine should be called by the constuctor of any subclass.
void HeightModel::initializeGases()
{
  allGasData.clear();
  allGasData.emplace(ARGON, mdArgonND());
  allGasData.emplace(CARBON_DIOXIDE, mdCarbonDioxideND());
  allGasData.emplace(CARBON_MONOXIDE, mdCarbonMonoxideND());
  allGasData.emplace(DIHYDROGEN, mdDihydrogenND());
  allGasData.emplace(DINITROGEN, mdDinitrogenND());
  allGasData.emplace(DIOXYGEN, mdDioxygenND());
  allGasData.emplace(HELIUM, mdHeliumND());
  allGasData.emplace(HYDROGEN, mdHydrogenND());
  allGasData.emplace(METHANE, mdMethaneND());
  allGasData.emplace(NITROGEN, mdNitrogenND());
  allGasData.emplace(OXYGEN, mdOxygenND());
  allGasData.emplace(OZONE, mdOzoneND());
  allGasData.emplace(NITROUS_OXIDE, mdNitrousOxideND());
}

//! \brief Compute the atmospheric state.
//!
//! Evaluate the atmospheric state for the current position.
//!
//! \b Inputs
//! \arg #position
//!
//! \retval #atmos
void HeightModel::updateAtmosphereState()
{
  // Convert height to z (km) with bounds check
  greal z = clamp(height, mdHeight()[0], mdHeight()[mdHeight().size() - 1]);

  // Find height index for vertical interpolation
  size_t i;
  for (i = 0; i < mdHeight().size() - 1; ++i) {
    if (z >= mdHeight()[i] && z <= mdHeight()[i + 1])
      break;
  }
  // index of the lower bound for height
  size_t lowZ = i;
  // index of the upper bound for height
  size_t highZ = i + 1;

  greal deltaZ = mdHeight()[highZ] - mdHeight()[lowZ];
  Interpolator zInterp;
  zInterp.makeFraction(mdHeight()[lowZ], mdHeight()[highZ], z);

  // Linear height interpolation for temperature
  temperature = zInterp.linear(mdTemperature()[lowZ], mdTemperature()[highZ]);

  // Assumptions (checking occurs in debug mode only)
  assert(mdPressure()[highZ] != 0.0);
  assert(mdPressure()[lowZ] / mdPressure()[highZ] > 0.0);

  // Interpolation on log pressure keeps hydrostatics consistent
  // Pressure scale height and vertical interpolation for pressure
  pressure = zInterp.log(mdPressure()[lowZ], mdPressure()[highZ]);
  pressureScaleHeight = deltaZ / log(mdPressure()[lowZ] / mdPressure()[highZ]);

  // Gas constant at lower height index
  greal rLow = mdPressure()[lowZ] / (mdDensity()[lowZ] * mdTemperature()[lowZ]);
  // Gas constant at upper height index
  greal rHigh = mdPressure()[highZ] / (mdDensity()[highZ] * mdTemperature()[highZ]);
  // Linear height interpolation for gas constant
  specificGasConstant = zInterp.linear(rLow, rHigh);

  // Assumptions (checking occurs in debug mode only)
  assert(specificGasConstant != 0.0);
  assert(temperature != 0.0);

  // Density from perfect gas law
  density = pressure / (specificGasConstant * temperature);

  // Assumptions (checking occurs in debug mode only)
  assert(mdDensity()[highZ] != 0.0);
  assert(mdDensity()[lowZ] / mdDensity()[highZ] > 0.0);

  // Density scale height between two index levels
  densityScaleHeight = deltaZ / log(mdDensity()[lowZ] / mdDensity()[highZ]);

  //============================================================================

  // Compute number densities for active gases
  // Iterate over the set of all activated gases.
  // This macro defines "gas" and "gasType" for each activated gas.
  FOR_ALL_GASES_PRESENT(
    // User the type to find the appropriate gas data
    auto& gasData = allGasData.at(gasType).get();
    if (gasData.empty()) {
      // No gas data found, zero out the number density.
      // This situation should only occur if the model has activated
      // a gas for which there is no data.  The assumption here is
      // that an analytical model will compute this value later.
      gas.numberDensity = 0;
    }
    else {
      // Use log interpolation on gas data
      gas.numberDensity = zInterp.log(gasData[lowZ], gasData[highZ]);
    }
  ELSE_NOT_PRESENT
    // If the gas is not active, then zero out the number density.
    gas.numberDensity = 0;
  );

  // Compute total number density
  updateTotalNumberDensity();

  // If no gas data is present, then there may be no numberDensity
  if (totalNumberDensity != 0.0) {

    // Get average molecular weight
    averageMolecularWeight = getTotalMass() / totalNumberDensity;

    // Compute mole and mass fractions
    updateMoleFractions();
    updateMassFractions();
  }
  else {
    averageMolecularWeight = 0.0;
  }
}


//! \brief Evaluates atmospheric reference values.
//!
//! From average HeightModel model in: "Engineering Models for Titan's Atmosphere", "Huygens
//! Science, Payload and Mission", ESA Report SP-1177, Aug. 1977, pages 243-256.
//! \param height The height (km) of the desired reference values.
//! \param[out] refTemperature A greal.
//! \param[out] refPressure A greal.
//! \param[out] refDensity A greal.
//! \retval refTemperature The reference temperature.
//! \retval refPressure The reference pressure.
//! \retval refDensity The reference density.
void HeightModel::getReferenceValues(greal height, greal& refTemperature, greal& refPressure, greal& refDensity)
{
  // Convert height to z (km) with bounds check
  greal z = clamp(height, mdHeight()[0], mdHeight()[mdHeight().size() - 1]);

  // Find height index for vertical interpolation
  size_t i;
  for (i = 0; i < mdHeight().size() - 1; ++i) {
    if (z >= mdHeight()[i] && z <= mdHeight()[i + 1])
      break;
  }
  // index of the lower bound for height
  size_t lowZ = i;
  // index of the upper bound for height
  size_t highZ = i + 1;


  Interpolator zInterp;
  zInterp.makeFraction(mdHeight()[lowZ], mdHeight()[highZ], z);

  // Linear height interpolation on temperature
  refTemperature = zInterp.linear(mdTemperature()[lowZ], mdTemperature()[highZ]);

  // Pressure scale height and vertical pressure interpolation
  refPressure = zInterp.log(mdPressure()[lowZ], mdPressure()[highZ]);

  // Assumptions (checking occurs in debug mode only)
  assert(mdDensity()[lowZ] != 0.0);
  assert(mdTemperature()[lowZ] != 0.0);
  assert(mdDensity()[highZ] != 0.0);
  assert(mdTemperature()[highZ] != 0.0);

  // Linear height interpolation for gas constant
  greal rLow = mdPressure()[lowZ] / (mdDensity()[lowZ] * mdTemperature()[lowZ]);
  greal rHigh = mdPressure()[highZ] / (mdDensity()[highZ] * mdTemperature()[highZ]);
  greal R = zInterp.linear(rLow, rHigh);

  // Assumptions (checking occurs in debug mode only)
  assert(R != 0.0);
  assert(refTemperature != 0.0);

  // Density from perfect gas law
  refDensity = refPressure / (R * refTemperature);
}

//! \brief Evaluate the pressure at the surface.
//!
//! \returns Pressure at surface.
greal HeightModel::getPressureAtSurface() const
{
  // find the index of the first non-negative height
  size_t i = 0;
  while (i < mdHeight().size()) {
    if (mdHeight()[i] < 0.0) {
      ++i;
    }
    else {
      // return pressure at that index
      return mdPressure()[i];
    }
  }
  // No non-negative heights!  Really?
  return mdPressure()[i-1];
}

} // namespace