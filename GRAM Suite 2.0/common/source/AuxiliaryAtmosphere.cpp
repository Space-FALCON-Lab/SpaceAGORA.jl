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

#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include "AuxiliaryAtmosphere.h"
#include "Interpolator.h"
#include "error_strings.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
  AuxiliaryAtmosphere::AuxiliaryAtmosphere()
{
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
AuxiliaryAtmosphere::AuxiliaryAtmosphere(const AuxiliaryAtmosphere& orig)
{
  auxDataTable = orig.auxDataTable;
  inputState = orig.inputState;
  auxFileName = orig.auxFileName;
  innerRadius = orig.innerRadius;
  outerRadius = orig.outerRadius;
  eastLongitudePositive = orig.eastLongitudePositive;
  active = orig.active;
  auxId = orig.auxId;
  profileWeight = orig.profileWeight;
  hasStandardDeviations = orig.hasStandardDeviations;
  auxValues = orig.auxValues;
}

//! \copydoc Atmosphere::~Atmosphere()
AuxiliaryAtmosphere::~AuxiliaryAtmosphere()
{
  auxDataTable.clear();
}

//! \brief Set the inner lat/lon radius for fairing with the auxiliary data.
//!
//! Within the inner radius, the auxiliary data will replace the incoming state.
//! Outside of the outer radius, the auxiliary data will have no effect on the incoming state.
//! The incoming and auxiliary data will be faired in between the inner and outer
//! radii with auxiliary data having weight 0 at the outer radius and weight 1 at the inner radius.
//! \param inner The inner lat/lon radius.
//! \param outer The outer lat/lon radius.
void AuxiliaryAtmosphere::setBounds(greal inner, greal outer)
{
  if (inner >= outer) {
    throw string("Error: Invalid auxiliary atmosphere bounds.\n"
                 "       The inner radius must be less than the outer radius.");
  }
  innerRadius = max(inner, (greal)0.0);
  outerRadius = max(outer, innerRadius);
}

//! \brief Sets the auxiliary data file name.
//!
//! Sets the path to the auxiliary data file and marks the auxiliary object as active.
//! \param fileName The complete path to the auxiliary data.
void AuxiliaryAtmosphere::setDataFile(const std::string& fileName)
{
  auxFileName = fileName;
}

//! \brief Sets the auxiliary data.
//!
//! Typically auxiliary data is loaded via loadData().  Use this function to
//! provide auxiliary data from an alternative data source.
//! \param dataTable A vector of auxiliary data.
void AuxiliaryAtmosphere::setData(const std::vector<AuxiliaryData>& dataTable)
{
  auxDataTable = dataTable;
  // Data should be increasing  by height
  if (auxDataTable.at(1).height < auxDataTable.at(0).height) {
    sort(auxDataTable.begin(), auxDataTable.end(), [](AuxiliaryData a, AuxiliaryData b) { return a.height < b.height; });
  }
  active = true;
}

void AuxiliaryAtmosphere::setValues(greal dens, greal pres, greal temp, greal ew, greal ns)
{
  active = true;
  useValues = true;
  auxValues.density = dens;
  auxValues.pressure = pres;
  auxValues.temperature = temp;
  auxValues.ewWind = ew;
  auxValues.nsWind = ns;
}

//! \brief Reads the auxiliary data.
//!
//! If a file name has been set, then this function reads in the auxiliary data.
//! Subclass and override this function to use an alternative data source.  Or set
//! the auxiliary data directly with setData().
void AuxiliaryAtmosphere::loadData()
{
  if (active || auxFileName.empty()) {
    return;
  }
 
  auxDataTable.clear();
  AuxiliaryData auxData;
  ifstream auxFile(auxFileName);
  if (!auxFile) {
    throw string(FILE_OPEN_ERROR_MESSAGE + auxFileName);
  }

  try {
    auxFile.ignore(256, '\n');
    greal lastHeight = -9999999.9;
    while (!auxFile.eof()) {
      string lineBuffer = "";
      getline(auxFile, lineBuffer);
      if (auxFile.bad()) {
        throw string("Unable to parse line after height = ") + to_string(auxDataTable.back().height);
      }

      if (!lineBuffer.empty() && lineBuffer.find_first_not_of(" \t") != lineBuffer.npos) {
        istringstream lineInput(lineBuffer);
        lineInput >> auxData.height >> auxData.latitude >> auxData.longitude
          >> auxData.temperature >> auxData.pressure >> auxData.density
          >> auxData.ewWind >> auxData.nsWind;
        if (hasStandardDeviations) {
          lineInput >> auxData.temperatureStandardDeviation >> auxData.pressureStandardDeviation
            >> auxData.densityStandardDeviation >> auxData.ewStandardDeviation >> auxData.nsStandardDeviation;
        }
        if (!lineInput.fail()) {
          if (auxDataTable.size() > 0 && auxData.height == lastHeight) {
            throw string("Auxiliary data is not monotonic at height = ") + to_string(auxData.height);
          }
          lastHeight = auxData.height;
          // Convert negative longitudes
          auxData.longitude = wrapDegrees(auxData.longitude);
          // Convert to East Longitude Positive if necessary
          if (!eastLongitudePositive) {
            auxData.longitude = 360.0_deg - auxData.longitude;
          }
          if (hasStandardDeviations) {
            if (auxData.temperature <= 0.0) {
              auxData.temperatureStandardDeviation = 0.0;
            }
            else if (auxData.temperatureStandardDeviation > 0.3 * auxData.temperature) {
              auxData.temperatureStandardDeviation = 0.3;
            }
            else {
              auxData.temperatureStandardDeviation /= auxData.temperature;
            }
            if (auxData.pressure <= 0.0) {
              auxData.pressureStandardDeviation = 0.0;
            }
            else if (auxData.pressureStandardDeviation > 0.3 * auxData.pressure) {
              auxData.pressureStandardDeviation = 0.3;
            }
            else {
              auxData.pressureStandardDeviation /= auxData.pressure;
            }
            if (auxData.density <= 0.0) {
              auxData.densityStandardDeviation = 0.0;
            }
            else if (auxData.densityStandardDeviation > 0.3 * auxData.density) {
              auxData.densityStandardDeviation = 0.3;
            }
            else {
              auxData.densityStandardDeviation /= auxData.density;
            }
          }
          auxDataTable.push_back(auxData);
        }
        else {
          throw string("Unable to parse line: ") + lineBuffer;
        }
      }
    }
    auxFile.close();
    // Data should be increasing  by height
    if (auxDataTable.at(1).height < auxDataTable.at(0).height) {
      sort(auxDataTable.begin(), auxDataTable.end(), [](AuxiliaryData a, AuxiliaryData b) { return a.height < b.height; });
    }
    active = true;
  }
  catch (const string& msg) {
    active = false;
    throw(string("Error: Auxiliary atmosphere parsing error.\n       ")
          + msg + "\n       Attempting to parse: " + auxFileName
          + "\n       Deactivating this auxiliary atmosphere.");
  }
}

//! \copydoc Atmosphere::update()
void AuxiliaryAtmosphere::update()
{
  atmos = inputState;
  if (!active) {
    return;
  }
  AuxiliaryData auxData;
  if (useValues) {
    useValues = false;
    active = false;
    density = auxValues.density;
    pressure = auxValues.pressure;
    temperature = auxValues.temperature;
    ewWind = auxValues.ewWind;
    nsWind = auxValues.nsWind;
    profileWeight = 1.0;
    return;
  }

  if (height < auxDataTable.front().height ) {
    auxData = auxDataTable.front();
    profileWeight = 0.0;
  }
  else if (height > auxDataTable.back().height) {
    auxData = auxDataTable.back();
    profileWeight = 0.0;
  }
  else {
    AuxiliaryData auxDataLow, auxDataHigh;
    // Find index values bracketing current height
    size_t low;
    for (low = 0; low < auxDataTable.size() - 1; ++low) {
      if (height >= auxDataTable[low].height && height <= auxDataTable[low+1].height) {
        break;
      }
    }
    size_t high = low + 1;
    Interpolator zInterp;
    zInterp.makeFraction(auxDataTable[low].height, auxDataTable[high].height, height);
    auxData.latitude = zInterp.linear(auxDataTable[low].latitude, auxDataTable[high].latitude);
    auxData.longitude = zInterp.linear(auxDataTable[low].longitude, auxDataTable[high].longitude);
    auxData.temperature = zInterp.linear(auxDataTable[low].temperature, auxDataTable[high].temperature);
    auxData.nsWind = zInterp.linear(auxDataTable[low].nsWind, auxDataTable[high].nsWind);
    auxData.ewWind = zInterp.linear(auxDataTable[low].ewWind, auxDataTable[high].ewWind);
    // Power-law interpolation for density (unless profile density
    // is zero, for which zero weight will be used)
    if (auxDataTable[low].density > 0.0) {
      auxData.density = zInterp.log(auxDataTable[low].density, auxDataTable[high].density);
    }
    else {
      auxData.density = 0.0;
    }
    // Power-law interpolation for pressure (unless profile pressure
    // is zero, for which zero weight will be used)
    if (auxDataTable[low].pressure > 0.0) {
      auxData.pressure = zInterp.log(auxDataTable[low].pressure, auxDataTable[high].pressure);
    }
    else {
      auxData.pressure = 0.0;
    }

    // Initialize weighting factor components for height and lat-lon
    greal facthgt = 1.0;
    greal factll = 1.0;
    size_t last = auxDataTable.size() - 1;
    if (height <= auxDataTable[1].height) {
      // Sin-squared variation of height weighting from 0 at 1st point to 1 at 2nd point
      facthgt = (height - auxDataTable[0].height) / (auxDataTable[1].height - auxDataTable[0].height);
      facthgt = pow(sin(HALF_PI * facthgt), 2);
    }
    else if (height >= auxDataTable[last - 1].height) {
      // Sin-squared variation of height weighting from 0 at next-to-
      // last point to 1 at last point
      facthgt=(height - auxDataTable[last].height) / (auxDataTable[last - 1].height - auxDataTable[last].height);
      facthgt = pow(sin(HALF_PI * facthgt), 2);
    }

    // Get the radius (in degrees) from site center to current location.
    greal radius = getArcAngle(latitude, longitude, auxData.latitude, auxData.longitude);

    // Use weight=0 if radius>outerRadius, weight=1 if radius<profnear,
    // with sin-squared variation between outerRadius and innerRadius
    if (radius >= outerRadius) {
      factll = 0.0;
    }
    else if (radius <= innerRadius) {
      factll = 1.0;
    }
    else {
      // Assumptions (checking occurs in debug mode only)
      assert(outerRadius != innerRadius);

      factll = (outerRadius - radius)/(outerRadius - innerRadius);
      factll = pow(sin(HALF_PI * factll), 2);
    }
    // Total weight = product of weights for lat-lon and height
    profileWeight = factll * facthgt;
  }
  greal tpdwgt = profileWeight;
  greal uvwgt = profileWeight;
  // Set profile weight to zero for p,d, & t if profile values are 0
  if (auxData.temperature * auxData.pressure * auxData.density == 0.0) {
    tpdwgt = 0.0;
  }
  // Set profile weight to zero for u & v if profile values are 0
  if (fabs(auxData.nsWind) + fabs(auxData.ewWind) == 0.0) {
    uvwgt = 0.0;
  }
  // Apply weighted averaging of profile values with input values
  Interpolator tpdInterp(tpdwgt);
  atmos.temperature = tpdInterp.linear(inputState.temperature, auxData.temperature);
  atmos.pressure = tpdInterp.linear(inputState.pressure, auxData.pressure);
  atmos.density = tpdInterp.linear(inputState.density, auxData.density);

  Interpolator uvInterp(uvwgt);
  atmos.nsWind = uvInterp.linear(inputState.nsWind, auxData.nsWind);
  atmos.ewWind = uvInterp.linear(inputState.ewWind, auxData.ewWind);
  atmos.profileWeight[auxId] = profileWeight;
}

void AuxiliaryAtmosphere::updateStandardDeviations()
{
  atmos = inputState;
  if (!active || useValues) {
    return;
  }
  AuxiliaryData auxData;

  if (height < auxDataTable.front().height) {
    auxData = auxDataTable.front();
  }
  else if (height > auxDataTable.back().height) {
    auxData = auxDataTable.back();
  }
  else {
    AuxiliaryData auxDataLow, auxDataHigh;
    // Find index values bracketing current height
    size_t low;
    for (low = 0; low < auxDataTable.size() - 1; ++low) {
      if (height >= auxDataTable[low].height && height <= auxDataTable[low + 1].height) {
        break;
      }
    }
    size_t high = low + 1;
    Interpolator zInterp;
    zInterp.makeFraction(auxDataTable[low].height, auxDataTable[high].height, height);
    auxData.latitude = zInterp.linear(auxDataTable[low].latitude, auxDataTable[high].latitude);
    auxData.longitude = zInterp.linear(auxDataTable[low].longitude, auxDataTable[high].longitude);
    auxData.temperatureStandardDeviation = zInterp.linear(auxDataTable[low].temperatureStandardDeviation, auxDataTable[high].temperatureStandardDeviation);
    auxData.pressureStandardDeviation = zInterp.linear(auxDataTable[low].pressureStandardDeviation, auxDataTable[high].pressureStandardDeviation);
    auxData.densityStandardDeviation = zInterp.linear(auxDataTable[low].densityStandardDeviation, auxDataTable[high].densityStandardDeviation);
    auxData.nsStandardDeviation = zInterp.linear(auxDataTable[low].nsStandardDeviation, auxDataTable[high].nsStandardDeviation);
    auxData.ewStandardDeviation = zInterp.linear(auxDataTable[low].ewStandardDeviation, auxDataTable[high].ewStandardDeviation);

    // Initialize weighting factor components for height and lat-lon
    greal facthgt = 1.0;
    greal factll = 1.0;
    size_t last = auxDataTable.size() - 1;
    if (height <= auxDataTable[1].height) {
      // Sin-squared variation of height weighting from 0 at 1st point to 1 at 2nd point
      facthgt = (height - auxDataTable[0].height) / (auxDataTable[1].height - auxDataTable[0].height);
      facthgt = pow(sin(HALF_PI * facthgt), 2);
    }
    else if (height >= auxDataTable[last - 1].height) {
      // Sin-squared variation of height weighting from 0 at next-to-
      // last point to 1 at last point
      facthgt = (height - auxDataTable[last].height) / (auxDataTable[last - 1].height - auxDataTable[last].height);
      facthgt = pow(sin(HALF_PI * facthgt), 2);
    }

    // Compute absolute lat-lon difference of current position from
    // profile lat-lon
    greal dlat = fabs(latitude - auxData.latitude);
    greal dlon = fabs(longitude - auxData.longitude);

    // Adjust lon difference for wrap at lon 360
    if (dlon > 180.0_deg) {
      dlon = 360.0_deg - dlon;
    }

    // Lat-lon radius of current position from profile lat-lon
    greal radius = sqrt(dlat * dlat + dlon * dlon);

    // Use weight=0 if radius>proffar, weight=1 if radius<profnear,
    // with sin-squared variation between proffar and profnear
    if (radius >= outerRadius) {
      factll = 0.0;
    }
    else if (radius <= innerRadius) {
      factll = 1.0;
    }
    else {
      // Assumptions (checking occurs in debug mode only)
      assert(outerRadius != innerRadius);

      factll = (outerRadius - radius) / (outerRadius - innerRadius);
      factll = pow(sin(HALF_PI * factll), 2);
    }
    // Total weight = product of weights for lat-lon and height
    profileWeight = factll * facthgt;
  }
  greal tpdwgt = profileWeight;
  greal uvwgt = profileWeight;
  // Set profile weight to zero for p,d, & t if profile values are 0
  if (auxData.temperatureStandardDeviation * auxData.pressureStandardDeviation * auxData.densityStandardDeviation == 0.0) {
    tpdwgt = 0.0;
  }
  // Set profile weight to zero for u & v if profile values are 0
  if (fabs(auxData.nsStandardDeviation) + fabs(auxData.ewStandardDeviation) == 0.0) {
    uvwgt = 0.0;
  }
  // Apply weighted averaging of profile values with input values
  Interpolator tpdInterp(tpdwgt);
  atmos.temperatureStandardDeviation = tpdInterp.linear(inputState.temperatureStandardDeviation, auxData.temperatureStandardDeviation);
  atmos.pressureStandardDeviation = tpdInterp.linear(inputState.pressureStandardDeviation, auxData.pressureStandardDeviation);
  atmos.densityStandardDeviation = tpdInterp.linear(inputState.densityStandardDeviation, auxData.densityStandardDeviation);

  Interpolator uvInterp(uvwgt);
  atmos.nsStandardDeviation = uvInterp.linear(inputState.nsStandardDeviation, auxData.nsStandardDeviation);
  atmos.ewStandardDeviation = uvInterp.linear(inputState.ewStandardDeviation, auxData.ewStandardDeviation);
  atmos.profileWeight[auxId] = profileWeight;
}

//! \fn void AuxiliaryAtmosphere::setInputState(const AtmosphereState& state)
//! \brief Set the incoming atmospheric state.
//!
//! The incoming atmospheric state will be faired appropriately with the auxiliary
//! data provided.
//! \param state Then incoming atmospheric state.

//! \fn void AuxiliaryAtmosphere::setEastLongitudePositive(bool flag)
//! \brief Set the positive longitude convention.
//!
//! \param flag True for east, false for west.

//! \fn void AuxiliaryAtmosphere::setWestLongitudePositive()
//! \brief Set the positive longitude convention to west.

//! \fn void AuxiliaryAtmosphere::setAuxId(size_t id)
//! \brief Set the index of the profileWeight in the AtmosphereState.
//!
//! \param id The index for the profileWeight.


} // namespace
