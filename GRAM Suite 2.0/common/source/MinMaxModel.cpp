//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//
// Original models developed by Dr. C. G. (Jere) Justus.
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include "MinMaxModel.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
MinMaxModel::MinMaxModel()
: Atmosphere()
{
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
MinMaxModel::MinMaxModel(const MinMaxModel& orig)
  : Atmosphere(orig)
{
  userMinMaxFactor = orig.userMinMaxFactor;
  computeMinMaxFactor = orig.computeMinMaxFactor;
}

//! \fn  MinMaxModel::~MinMaxModel()
//! \copydoc Atmosphere::~Atmosphere()

//! \brief Adds all gases into a map for easy processing.
//!
//! This routine should be called by the constuctor of any subclass.
void MinMaxModel::initializeGases()
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

//! \fn  MinMaxModel::updateMinMaxFactor();
//! \brief Computes a min/max factor based on seasonal and diurnal data.
//!
//! This method is a pure virtual function and should be overridden in each subclass.
//! The minMaxFactor should be computed as a value between -1 and 1 with the following
//! correspondences: -1 corresponds to the minimum model, 0 corresponds to the average model,
//! and 1 corresponds to the maximum model.

//! \fn  MinMaxModel::getMinMaxFactor();
//! \brief Retrieves the computed minMaxFactor.
//!
//! \returns The minMaxFactor.

//! \brief Set the interpolation factor between models.
//!
//! This factor is used to interpolate between the minimum model (-1)
//! and the average model (0) or between the average model (0)
//! and the maximum model (1).
//! \param factor The user supplied value between -1 and 1.
//! \param computeFlag If true, the the minMaxFactor will be evaluated and
//!        scaled with the input value each update. If false, then the input value is used.
void MinMaxModel::setMinMaxFactor(greal factor, bool computeFlag)
{
  computeMinMaxFactor = computeFlag;
  // make sure massDensityFactor is in [-1, 1]
  userMinMaxFactor = clamp(factor, -1.0, 1.0);
  // If we are not computing the minMaxFactor, then we should set it
  // to the user supplied factor.  This can also be done in a model's
  // override of updateMinMaxFactor().
  if (!computeFlag) {
    minMaxFactor = userMinMaxFactor;
  }
}

//! \brief Compute the atmospheric state.
//!
//! Evaluate the atmospheric state for the current position.
//! In this model, values are the interpolated between two of three models (min, avg, max). 
//! A min-max factor determines whether the
//! the interpolation is between the min and average models or between the average
//! and max models.
//!
//! \b Inputs
//! \arg #position 
//! \arg #minMaxFactor 
//!
//! \retval #atmos 
void MinMaxModel::updateAtmosphereState()
{
  // Interpolate between avg and max if Fminmax > 0
  size_t minMax = MAX_IDX;
  // Interpolate between avg and min if Fminmax < 0
  if (minMaxFactor < 0.0) {
    minMax = MIN_IDX;
  }

  greal modelHeight = clamp(height, mdHeight()[0], mdHeight()[mdHeight().size() - 1]);

  // Find height index for vertical interpolation
  size_t i;
  for (i = 0; i < mdHeight().size() - 1; ++i) {
    if (modelHeight >= mdHeight()[i] && modelHeight <= mdHeight()[i + 1])
      break;
  }
  // index of the lower bound for height
  size_t lowZ = i;
  // index of the upper bound for height
  size_t highZ = i + 1;

  // Interpolate between avg-min or avg-max, based on minMaxFactor value
  Interpolator fInterp(abs(minMaxFactor));

  greal deltaZ = mdHeight()[highZ] - mdHeight()[lowZ];
  Interpolator zInterp;
  zInterp.makeFraction(mdHeight()[lowZ], mdHeight()[highZ], modelHeight);

  // Interpolation on inverse temperature preserves right scale height
  // Linear height interpolation for temperature
  greal t_low = fInterp.inverse(mdTemperature()[AVG_IDX][lowZ], mdTemperature()[minMax][lowZ]);
  greal t_high = fInterp.inverse(mdTemperature()[AVG_IDX][highZ], mdTemperature()[minMax][highZ]);
  temperature = zInterp.linear(t_low, t_high);

  // Interpolation on log pressure keeps hydrostatics consistent
  greal p_low = fInterp.log(mdPressure()[AVG_IDX][lowZ], mdPressure()[minMax][lowZ]);
  greal p_high = fInterp.log(mdPressure()[AVG_IDX][highZ], mdPressure()[minMax][highZ]);

  // Assumptions (checking occurs in debug mode only)
  assert(p_high != 0.0);
  assert(p_low / p_high > 0.0);

  // Pressure scale height and vertical interpolation for pressure
  pressure = zInterp.log(p_low, p_high);
  pressureScaleHeight = deltaZ / log(p_low / p_high);

  // Assumptions (checking occurs in debug mode only)
  assert(mdDensity()[AVG_IDX][lowZ] != 0.0);
  assert(mdTemperature()[AVG_IDX][lowZ] != 0.0);
  assert(mdDensity()[minMax][lowZ] != 0.0);
  assert(mdTemperature()[minMax][lowZ] != 0.0);
  assert(mdDensity()[AVG_IDX][highZ] != 0.0);
  assert(mdTemperature()[AVG_IDX][highZ] != 0.0);
  assert(mdDensity()[minMax][highZ] != 0.0);
  assert(mdTemperature()[minMax][highZ] != 0.0);

  // Interpolate on gas law constant to get density
  greal r_mean = mdPressure()[AVG_IDX][lowZ] / (mdDensity()[AVG_IDX][lowZ] * mdTemperature()[AVG_IDX][lowZ]);
  greal r_minmax = mdPressure()[minMax][lowZ] / (mdDensity()[minMax][lowZ] * mdTemperature()[minMax][lowZ]);
  // Gas constant at lower height index
  greal r_low = fInterp.linear(r_mean, r_minmax);
  r_mean = mdPressure()[AVG_IDX][highZ] / (mdDensity()[AVG_IDX][highZ] * mdTemperature()[AVG_IDX][highZ]);
  r_minmax = mdPressure()[minMax][highZ] / (mdDensity()[minMax][highZ] * mdTemperature()[minMax][highZ]);
  // Gas constant at upper height index
  greal r_high = fInterp.linear(r_mean, r_minmax);
  // Linear height interpolation for gas constant
  specificGasConstant = zInterp.linear(r_low, r_high);
  // Density from perfect gas law
  density = pressure / (specificGasConstant * temperature);

  // Assumptions (checking occurs in debug mode only)
  assert(r_low != 0.0);
  assert(t_low != 0.0);
  assert(r_high != 0.0);
  assert(t_high != 0.0);

  // Density scale height between two index levels
  greal d_low = p_low / (r_low * t_low);
  greal d_high = p_high / (r_high * t_high);

  // Assumptions (checking occurs in debug mode only)
  assert(d_high != 0.0);
  assert(d_low / d_high > 0.0);

  densityScaleHeight = deltaZ / log(d_low / d_high);

  //============================================================================

  // avg-min or avg-max interpolation of number densities
  // vertical interpolations for number densities
  greal nd_low;
  greal nd_high;

  // Iterate over the set of all activated gases.
  // This macro defines "gas" and "gasType" for each activated gas.
  FOR_ALL_GASES_PRESENT(
    // User the type to find the appropriate gas data
    auto& gasData = allGasData.at(gasType).get();
    // Be sure data is present before interpolating
    if (gasData[0].empty()) {
      // No gas data found, zero out the number density.
      // This situation should only occur if the model has activated
      // a gas for which there is no data.  The assumption here is
      // that an analytical model will compute this value later.
      gas.numberDensity = 0;
    }
    else {
      // Compute number density
      nd_low = fInterp.log(gasData[AVG_IDX][lowZ], gasData[minMax][lowZ]);
      nd_high = fInterp.log(gasData[AVG_IDX][highZ], gasData[minMax][highZ]);
      gas.numberDensity = zInterp.log(nd_low, nd_high);
    }
  ELSE_NOT_PRESENT
    // If the gas is not active, then zero out the number density.
    gas.numberDensity = 0;
  );

  // Compute total number density
  updateTotalNumberDensity();

  // Assumptions (checking occurs in debug mode only)
  assert(totalNumberDensity != 0.0);

  // Get average molecular weight
  averageMolecularWeight = getTotalMass() / totalNumberDensity;

  // Update mole and mass fractions
  updateMoleFractions();
  updateMassFractions();
}


//! \brief Evaluates atmospheric reference values.
//!
//! From average MinMaxModel model in: "Engineering Models for Titan's Atmosphere", "Huygens
//! Science, Payload and Mission", ESA Report SP-1177, Aug. 1977, pages 243-256.
//! \param height \copydoc Position::height
//! \param[out] refTemperature A greal.
//! \param[out] refPressure A greal.
//! \param[out] refDensity A greal.
//! \retval refTemperature The reference temperature.
//! \retval refPressure The reference pressure.
//! \retval refDensity The reference density.
void MinMaxModel::getReferenceValues(greal height, greal& refTemperature, greal& refPressure, greal& refDensity)
{
  // Convert height to z (km) within bounds
  greal z = clamp(height, mdHeight()[0], mdHeight()[mdHeight().size() - 1]);

  // Find height index for vertical interpolation
  size_t i = 0;
  for( size_t iz = 0; iz < mdHeight().size() - 1; ++iz) {
    if (mdHeight()[iz] <= z && z <= mdHeight()[iz+1]) {
      i = iz;
      break;
    }
  }

  // index of the lower bound for height
  size_t lowZ = i;
  // index of the upper bound for height
  size_t highZ = i + 1;

  Interpolator zInterp;
  zInterp.makeFraction(mdHeight()[i], mdHeight()[i + 1], z);

  // Linear height interpolation on temperature
  refTemperature = zInterp.linear(mdTemperature()[AVG_IDX][lowZ], mdTemperature()[AVG_IDX][highZ]);

  // Pressure scale height and vertical pressure interpolation
  refPressure = zInterp.log(mdPressure()[AVG_IDX][lowZ], mdPressure()[AVG_IDX][highZ]);

  // Assumptions (checking occurs in debug mode only)
  assert(mdDensity()[AVG_IDX][lowZ] != 0.0);
  assert(mdTemperature()[AVG_IDX][lowZ] != 0.0);
  assert(mdDensity()[AVG_IDX][highZ] != 0.0);
  assert(mdTemperature()[AVG_IDX][highZ] != 0.0);

  // Linear height interpolation for gas constant
  greal lowR = mdPressure()[AVG_IDX][lowZ] / (mdDensity()[AVG_IDX][lowZ] * mdTemperature()[AVG_IDX][lowZ]);
  greal highR = mdPressure()[AVG_IDX][highZ] / (mdDensity()[AVG_IDX][highZ] * mdTemperature()[AVG_IDX][highZ]);
  greal R = zInterp.linear(lowR, highR);

  // Assumptions (checking occurs in debug mode only)
  assert(R != 0.0);
  assert(refTemperature != 0.0);

  // Density from perfect gas law
  refDensity = refPressure / (R * refTemperature);
}

//! \brief Evaluate the pressure at the surface.
//!
//! \b Inputs
//! \arg #minMaxFactor 
//!
//! \returns Pressure at surface.
greal MinMaxModel::getPressureAtSurface() const
{
    int j = (minMaxFactor < 0.0) ? MIN_IDX : MAX_IDX;
    Interpolator fInterp( abs(minMaxFactor) );
    // find the index of the first non-negative height
    size_t i = 0;
    while (i < mdHeight().size()) {
      if (mdHeight()[i] < 0.0) {
        ++i;
      }
      else {
        // return pressure at that index
        // interpolate over the minmax factor
        return fInterp.log(mdPressure()[AVG_IDX][i], mdPressure()[j][i]);
      }
    }

    // No non-negative heights!  Really?
    return fInterp.log(mdPressure()[AVG_IDX][i - 1], mdPressure()[j][i - 1]);
}

} // namespace