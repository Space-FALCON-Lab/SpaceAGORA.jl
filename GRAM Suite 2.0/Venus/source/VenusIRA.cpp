//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//
// Adapted from Venus-GRAM 2005 developed by Dr. C. G. (Jere) Justus.
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include "VenusIRA.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {

bool VenusIRA::isInitialized = false;
vector<greal> VenusIRA::zref;
vector<greal> VenusIRA::pref;
vector<greal> VenusIRA::tref;
vector<greal> VenusIRA::dref;


//! \copydoc Atmosphere::Atmosphere()
VenusIRA::VenusIRA()
: Atmosphere(), VenusCommon(this)
{
  initializeData();
}

//! \fn  VenusIRA::VenusIRA(const VenusIRA& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn  VenusIRA::~VenusIRA()
//! \copydoc Atmosphere::~Atmosphere()

//! \copydoc Atmosphere::update()
void VenusIRA::update()
{
  updateAtmosphereState();
  updateWinds();
  updateMoleFractions();
  updateMassFractions();
}

//! \brief Compute dynamics and species concentrations.
//!
//! This model uses four layer models: low, mid, high, and thermos.  Fairing is performed between the layers.
//!
//! \b Inputs
//! \arg #position   
//!
//! \retval #pressure
//! \retval #density
//! \retval #temperature
//! \retval #compressibilityFactor
//! \retval #specificGasConstant
//! \retval #averageMolecularWeight
//! \retval #totalNumberDensity
//! \retval #pressureScaleHeight
//! \retval #densityScaleHeight
//! \retval numberDensity  The number density of species.
void VenusIRA::updateAtmosphereState()
{
  // Get the heights for the following indices for each VIRA model:
  // A = first (lowest), B = second, Y = penultimate, Z = ultimate (highest)
  greal lowA, lowB, lowY, lowZ;
  greal midA, midB, midY, midZ;
  greal highA, highB, highY, highZ;
  venusLow.getInterpolationLimits(lowA, lowB, lowY, lowZ);
  venusMid.getInterpolationLimits(midA, midB, midY, midZ);
  venusHigh.getInterpolationLimits(highA, highB, highY, highZ);

  // Convert to height z (km)
  greal z = height;
  // Limit heights to > 0 km
  if (z < lowA) {
    z = lowA;
  }

  // For heights 0-98 km use interpolation of low altitude VIRA data
  if (z <= lowY) {
    // Set the position in the low model.
    venusLow.setPosition(position);

    // Get the height at lower bounding index.
    greal posLow = venusLow.setLowerHeightIndex();

    // Get the atmosphere state at the lower bounding height.
    venusLow.update();
    const AtmosphereState atmosLow = venusLow.getAtmosphereState();

    // Get the height at upper bounding index.
    greal posHigh = venusLow.setUpperHeightIndex();

    // Get the atmosphere state at the upper bounding height.
    venusLow.update();
    const AtmosphereState atmosHigh = venusLow.getAtmosphereState();

    // Use height interpolation
    interpolateHeight(z, posLow, posHigh, atmosLow, atmosHigh);
  }

  // For heights 98-100 km, interpolate between
  // Low: low altitude VIRA data
  // High: a mean of the highest low altitude and lowest mid altitude data
  else if (z <= lowZ) {
    // Set the position in the low and mid models.
    venusLow.setPosition(position);
    venusMid.setPosition(position);

    // Set the ephemeris in the mid model.
    venusMid.setEphemerisState(ephem);

    // Get the height at lower bounding index (should be 98 km).
    greal posLow = venusLow.setLowerHeightIndex();

    // Get the atmosphere state at the lower bounding height.
    venusLow.update();
    const AtmosphereState atmosLow = venusLow.getAtmosphereState();

    // Now we want to average the overlapping heights in the low and mid models.
    // First, get the height at upper bounding index from the low model (should be 100 km).
    greal posHigh = venusLow.setUpperHeightIndex();

    // Get the (low) atmosphere state at the upper bounding height.
    venusLow.update();
    AtmosphereState atmosHigh1 = venusLow.getAtmosphereState();

    // Next, get the height at upper bounding index from the mid model (should be 100 km).
    posHigh = venusMid.setLowerHeightIndex();

    // Get the (mid) atmosphere state at the upper bounding height.
    venusMid.update();
    const AtmosphereState atmosHigh2 = venusMid.getAtmosphereState();

    // Assume some numberDensities from the mid model.
    atmosHigh1.helium.numberDensity = atmosHigh2.helium.numberDensity;
    atmosHigh1.nitrogen.numberDensity = atmosHigh2.nitrogen.numberDensity;

    // Compute the atmosphere state for the upper bounding height by
    // averaging the low and mid model data.
    AtmosphereState atmosHigh = averageAtmospheres(atmosHigh1, atmosHigh2);

    // Use height interpolation
    interpolateHeight(z, posLow, posHigh, atmosLow, atmosHigh);
  }

  // For height 100-105 km, interpolate between
  // Low: a mean of the highest low altitude and lowest mid altitude data  
  // High: middle altitude VIRA data
  else if (z <= midB) {
    // Set the position in the mid and low models
    venusLow.setPosition(position);
    venusMid.setPosition(position);

    // Set the ephemeris in the mid model
    venusMid.setEphemerisState(ephem);

    // Now we want to average the overlapping heights in the low and mid models.
    // First, get the height at lower bounding index from the low model (should be 100 km).
    greal posLow = venusLow.setLowerHeightIndex();

    // Get the (low) atmosphere state at the lower bounding height.
    venusLow.update();
    AtmosphereState atmosLow1 = venusLow.getAtmosphereState();

    // Next, get the height at lower bounding index from the mid model (should be 100 km).
    posLow = venusMid.setLowerHeightIndex();

    // Get the (mid) atmosphere state at the lower bounding height.
    venusMid.update();
    const AtmosphereState atmosLow2 = venusMid.getAtmosphereState();

    // Assume some numberDensities from the mid model.
    atmosLow1.helium.numberDensity = atmosLow2.helium.numberDensity;
    atmosLow1.nitrogen.numberDensity = atmosLow2.nitrogen.numberDensity;

    // Compute the atmosphere state for the lower bounding height by
    // averaging the low and mid model data.
    AtmosphereState atmosLow = averageAtmospheres(atmosLow1, atmosLow2);

    // Get the height at upper bounding index from the mid model.
    greal posHigh = venusMid.setUpperHeightIndex();

    // Get the atmosphere state at the upper bounding height.
    venusMid.update();
    const AtmosphereState atmosHigh = venusMid.getAtmosphereState();

    // Use height interpolation
    interpolateHeight(z, posLow, posHigh, atmosLow, atmosHigh);
  }

  // For heights 105-145 km use interpolation of middle altitude VIRA data
  else if (z <= midY) {
    // Set the position and ephemeris in the mid model.
    venusMid.setPosition(position);
    venusMid.setEphemerisState(ephem);

    // Get the height at lower bounding index.
    greal posLow = venusMid.setLowerHeightIndex();

    // Get the atmosphere state at the lower bounding height.
    venusMid.update();
    AtmosphereState atmosLow = venusMid.getAtmosphereState();

    // Get the height at upper bounding index.
    greal posHigh = venusMid.setUpperHeightIndex();

    // Get the atmosphere state at the upper bounding height.
    venusMid.update();
    AtmosphereState atmosHigh = venusMid.getAtmosphereState();

    // Use height interpolation
    interpolateHeight(z, posLow, posHigh, atmosLow, atmosHigh);
  }

  // For heights 145-150 km, interpolate between
  // Low: mid altitude VIRA data
  // High: a mean of the highest mid altitude and lowest high altitude data
  else if (z <= midZ) {
    // Set the position in the high and mid models.
    venusMid.setPosition(position);
    venusHigh.setPosition(position);

    // Set the ephemeris in the high and mid models.
    venusMid.setEphemerisState(ephem);
    venusHigh.setEphemerisState(ephem);

    // Get the height at lower bounding index (should be 145 km).
    greal posLow = venusMid.setLowerHeightIndex();

    // Get the atmosphere state at the lower bounding height.
    venusMid.update();
    const AtmosphereState atmosLow = venusMid.getAtmosphereState();

    // Now we want to average the overlapping heights in the mid and high models.
    // First, get the height at upper bounding index from the mid model (should be 150 km).
    greal posHigh = venusMid.setUpperHeightIndex();

    // Get the (mid) atmosphere state at the upper bounding height.
    venusMid.update();
    AtmosphereState atmosHigh1 = venusMid.getAtmosphereState();

    // Next, get the height at upper bounding index from the high model (should be 150 km).
    posHigh = venusHigh.setLowerHeightIndex();

    // Get the (high) atmosphere state at the upper bounding height.
    venusHigh.update();
    const AtmosphereState atmosHigh2 = venusHigh.getAtmosphereState();

    // Assume some numberDensities from the high model.
    atmosHigh1.hydrogen.numberDensity = atmosHigh2.hydrogen.numberDensity;

    // Compute the atmosphere state for the upper bounding height by
    // averaging the mid and high model data.
    AtmosphereState atmosHigh = averageAtmospheres(atmosHigh1, atmosHigh2);

    // Use height interpolation
    interpolateHeight(z, posLow, posHigh, atmosLow, atmosHigh);
  }

  // For height 150-155 km, interpolate between
  // Low: a mean of the highest mid altitude and lowest high altitude data  
  // High: high altitude VIRA data
  else if (z <= highB) {
    // Set the position in the high and mid models.
    venusMid.setPosition(position);
    venusHigh.setPosition(position);

    // Set the ephemeris in the high and mid models.
    venusMid.setEphemerisState(ephem);
    venusHigh.setEphemerisState(ephem);

    // Now we want to average the overlapping heights in the mid and high models.
    // First, get the height at lower bounding index from the mid model (should be 150 km).
    greal posLow = venusMid.setUpperHeightIndex();

    // Get the (mid) atmosphere state at the lower bounding height.
    venusMid.update();
    AtmosphereState atmosLow1 = venusMid.getAtmosphereState();

    // Next, get the height at lower bounding index from the high model (should be 150 km).
    posLow = venusHigh.setLowerHeightIndex();

    // Get the (high) atmosphere state at the lower bounding height.
    venusHigh.update();
    const AtmosphereState atmosLow2 = venusHigh.getAtmosphereState();

    // Assume some numberDensities from the high model.
    atmosLow1.hydrogen.numberDensity = atmosLow2.hydrogen.numberDensity;
    
    // Compute the atmosphere state for the lower bounding height by
    // averaging the mid and high model data.
    AtmosphereState atmosLow = averageAtmospheres(atmosLow1, atmosLow2);

    // Get the height at upper bounding index from the high model.
    greal posHigh = venusHigh.setUpperHeightIndex();

    // Get the atmosphere state at the upper bounding height.
    venusHigh.update();
    AtmosphereState atmosHigh = venusHigh.getAtmosphereState();

    // Use height interpolation
    interpolateHeight(z, posLow, posHigh, atmosLow, atmosHigh);
  }

  // For height 155-250 km, interpolate using high altitude VIRA data
  else if (z <= highZ) {
    // Set the position and ephemeris in the high model.
    venusHigh.setPosition(position);
    venusHigh.setEphemerisState(ephem);

    // Get the height at lower bounding index.
    greal posLow = venusHigh.setLowerHeightIndex();

    // Get the atmosphere state at the lower bounding height.
    venusHigh.update();
    const AtmosphereState atmosLow = venusHigh.getAtmosphereState();

    // Get the height at upper bounding index.
    greal posHigh = venusHigh.setUpperHeightIndex();

    // Get the atmosphere state at the upper bounding height.
    venusHigh.update();
    const AtmosphereState atmosHigh = venusHigh.getAtmosphereState();

    // Use height interpolation
    interpolateHeight(z, posLow, posHigh, atmosLow, atmosHigh);
  }

  // Above 250 km, evaluate atmosphere with constant temperature thermosphere       
  else {
    // Get high altitude VIRA data at maximum height                         
    venusHigh.setPosition(position);
    venusHigh.setEphemerisState(ephem);
    greal maxHeight = venusHigh.setUpperHeightIndex();
    venusHigh.update();

    // Get the atmosphere state at 250 km
    AtmosphereState atmosHigh = venusHigh.getAtmosphereState();

    // The thermo model needs to know the conditions at 250 km
    venusThermos.setBase(maxHeight, atmosHigh);

    // Get thermospheric model values at current height
    venusThermos.setPosition(position);
    venusThermos.update();
    atmos = venusThermos.getAtmosphereState();
  }
}

//! \brief Use height interpolation between two atmosphere states.
//!
//! Given two atmosphere states at two heights, interpolation is performed for a specified height.
//!
//! \param z The height of the desired atmosphere state.
//! \param height1  The height of the first atmosphere state.
//! \param height2 The height of the second atmosphere state.
//! \param atmos1 The first atmosphere state.
//! \param atmos2 The second atmosphere state.
//!
//! \retval #pressure
//! \retval #density
//! \retval #temperature
//! \retval #compressibilityFactor
//! \retval #specificGasConstant
//! \retval #averageMolecularWeight
//! \retval #totalNumberDensity
//! \retval #pressureScaleHeight
//! \retval #densityScaleHeight
//! \retval numberDensity  The number density of species.
void VenusIRA::interpolateHeight(greal z, greal height1, greal height2,
  const AtmosphereState& atmos1, const AtmosphereState& atmos2)
{
  // Create a height interpolator
  Interpolator zInterp;
  zInterp.makeFraction(height1, height2, z);

  // Interpolation of temperature, pressure, and compressibility
  temperature = zInterp.linear(atmos1.temperature, atmos2.temperature);
  pressure = zInterp.log(atmos1.pressure, atmos2.pressure);
  compressibilityFactor = zInterp.linear(atmos1.compressibilityFactor, atmos2.compressibilityFactor);
  
  // Guard against missing gas constants.
  greal R1 = atmos1.specificGasConstant;
  greal R2 = atmos2.specificGasConstant;
  if (R1 == 0.0 || R2 == 0.0) {
    R1 = atmos1.pressure / (atmos1.density * atmos1.temperature);
    R2 = atmos2.pressure / (atmos2.density * atmos2.temperature);
  }
  specificGasConstant = zInterp.linear(R1, R2);

  // Density from gas law
  density = pressure / (compressibilityFactor * specificGasConstant * temperature);

  // Mean molecular weight
  averageMolecularWeight = UNIVERSAL_GAS * 1000.0 / specificGasConstant;

  // Logarithmic interpolation of number densities
  carbonDioxide.numberDensity = zInterp.log(atmos1.carbonDioxide.numberDensity, atmos2.carbonDioxide.numberDensity);
  dinitrogen.numberDensity = zInterp.log(atmos1.dinitrogen.numberDensity, atmos2.dinitrogen.numberDensity);
  oxygen.numberDensity = zInterp.log(atmos1.oxygen.numberDensity, atmos2.oxygen.numberDensity);
  carbonMonoxide.numberDensity = zInterp.log(atmos1.carbonMonoxide.numberDensity, atmos2.carbonMonoxide.numberDensity);
  helium.numberDensity = zInterp.log(atmos1.helium.numberDensity, atmos2.helium.numberDensity);
  nitrogen.numberDensity = zInterp.log(atmos1.nitrogen.numberDensity, atmos2.nitrogen.numberDensity);
  hydrogen.numberDensity = zInterp.log(atmos1.hydrogen.numberDensity, atmos2.hydrogen.numberDensity);

  // Sum the number densities
  totalNumberDensity = carbonDioxide.numberDensity + dinitrogen.numberDensity
    + oxygen.numberDensity + carbonMonoxide.numberDensity + helium.numberDensity
    + nitrogen.numberDensity + hydrogen.numberDensity;

  // Scale heights for pressure and density
  if (atmos1.pressure > atmos2.pressure) {
    pressureScaleHeight = (height2 - height1) / log(atmos1.pressure / atmos2.pressure);
  }
  else {
    pressureScaleHeight = 0.0;
  }
  if (atmos1.density > atmos2.density) {
    densityScaleHeight = (height2 - height1) / log(atmos1.density / atmos2.density);
  }
  else {
    densityScaleHeight = 0.0;
  }
}

//! \brief Compute the average of two atmosphere states.
//!
//! Given two atmosphere states, values are averaged.
//!
//! \param atmos1 The first atmosphere state.
//! \param atmos2 The second atmosphere state.
//!
//! \returns An atmosphere state containing the averaged values for pressure, density, temperature, and number densities.
const AtmosphereState VenusIRA::averageAtmospheres(const AtmosphereState& atmos1, const AtmosphereState& atmos2)
{
  AtmosphereState avg;
  // Logarithmic average pressure and density
  avg.pressure = sqrt(atmos1.pressure * atmos2.pressure);
  avg.density = sqrt(atmos1.density * atmos2.density);

  // Linear average temperature 
  avg.temperature = (atmos1.temperature + atmos2.temperature) / 2.0;

  // Logarithmic average number densities
  avg.carbonDioxide.numberDensity = sqrt(atmos1.carbonDioxide.numberDensity * atmos2.carbonDioxide.numberDensity);
  avg.dinitrogen.numberDensity = sqrt(atmos1.dinitrogen.numberDensity * atmos2.dinitrogen.numberDensity);
  avg.oxygen.numberDensity = sqrt(atmos1.oxygen.numberDensity * atmos2.oxygen.numberDensity);
  avg.carbonMonoxide.numberDensity = sqrt(atmos1.carbonMonoxide.numberDensity * atmos2.carbonMonoxide.numberDensity);
  avg.helium.numberDensity = sqrt(atmos1.helium.numberDensity * atmos2.helium.numberDensity);
  avg.nitrogen.numberDensity = sqrt(atmos1.nitrogen.numberDensity * atmos2.nitrogen.numberDensity);
  avg.hydrogen.numberDensity = sqrt(atmos1.hydrogen.numberDensity * atmos2.hydrogen.numberDensity);

  return avg;
}

//! \brief Compute zonal and meridional winds at height and latitude.
//!
//! Approximations to data in Fig. 3, page 466 and Fig. 5, page 469 of "Venus II", and Fig. 8,
//! page 696 of "Venus".
//!
//! \b Inputs
//! \arg #height   
//! \arg #latitude 
//! \arg #solarTime
//!
//! \retval #ewWind
//! \retval #nsWind
void VenusIRA::updateWinds()
{
  // Note: Many references adopt the convention for Venus that         
  // super-rotating (westward, or retrograde) zonal winds are          
  // positive.  We retain the traditional right-handed coordinate      
  // convention, whereby zonal winds are positive eastward and         
  // meridional winds are positive northward.                          

  // Retrograde zonal wind magnitude versus height (- for westward superrotation)                                                    
  greal ue = -1.0 - 99.0 * height / 80.0;
  // Decrease in ue above 80 km parameterized from page 333 of         
  // Lellouch et al., Icarus 110, 315-319 (1994) and Fig 2. of Hou and 
  // Farrell, J. Atmos. Sci. 44, 1049-1061 (1987).                     
  if (height >= 80.0) {
    ue = -100.0 * (1.0 - (height - 80.0) / 30.0);
  }
  ue = min(ue, (greal)0.0);

  // Meridional wind magnitude versus height for retrograde wind       
  greal ve = 0.1 * ue * cos(toRadians(2.25_deg * height));

  // Sub-solar to anti-solar diurnal wind (u peaks at terminator),     
  // parameterized from Fig. 2 of Zhang et al. J. Geophys. Res.        
  // 101(E10), 23,195-205 (1996) and Fig 4. of Bougher et al. Icarus,  
  // 73, 545-573 (1988)                                                
  greal ut = -3.0 * (height - 81.0);
  if (ut > 0.0) {
    ut = 0.0;
  }
  ut = max(ut, (greal)-237.0);

  // Time variation of sub-solar to anti-solar wind components         
  greal usas = ut * sin(toRadians(15.0_deg * (solarTime - 12.0)));
  greal vsas = -ut * cos(toRadians(15.0_deg * (solarTime - 12.0)));

  // Latitude variation of zonal and meridional wind                   
  ewWind = (ue + usas) * (1.0 - pow(latitude / 90.0_deg, 4.0));
  nsWind = ve * sin(toRadians(2.0 * latitude)) + vsas * sin(toRadians(latitude));
}

//! \brief Get reference data for a specific height.
//!
//! Evaluates Venus reference atmospheric parameters from approximate average VIRA data.
//!
//! \param height Height in km.
//! \param[out] refTemperature A greal.
//! \param[out] refPressure A greal.
//! \param[out] refDensity A greal.
//!
//! \retval refTemperature Reference temperature.
//! \retval refPressure Reference pressure.
//! \retval refDensity Reference density.
void VenusIRA::getReferenceValues(greal height, greal& refTemperature, greal& refPressure, greal& refDensity)
{
  // Get the height above which starts the thermosphere
  greal highZ = zref.back();

  // Convert to height z (km) with bounds check
  greal z = max(height, zref[0]);

  // If height is in the thermosphere
  if (z > highZ) {
    // Evaluate thermosphere at height z    
    // Data from high data table for SZA = 90 deg and 250 km boundary values                                                 
    AtmosphereState atmosBase;
    atmosBase.carbonDioxide.numberDensity = 3.07e6;
    atmosBase.dinitrogen.numberDensity = 8.61e8;
    atmosBase.oxygen.numberDensity = 1.05e12;
    atmosBase.carbonMonoxide.numberDensity = 1.13e9;
    atmosBase.helium.numberDensity = 3.07e12;
    atmosBase.nitrogen.numberDensity = 2.18e10;
    atmosBase.hydrogen.numberDensity = 1.05e12;
    atmosBase.temperature = 230.0;

    // The thermo model needs to know the conditions at 250 km
    venusThermos.setBase(highZ, atmosBase);

    // The thermo model only requires height and latitude
    Position pos;
    pos.height = z;
    pos.latitude = 0.0;
    venusThermos.setPosition(pos);
    venusThermos.update();                                   

    // Set the reference values
    refPressure = venusThermos.getAtmosphereState().pressure;
    refDensity = venusThermos.getAtmosphereState().density;
    refTemperature = 230.0;
  }

  // Non-thermosphere reference computations.
  else {
    // Find height index for vertical interpolation                      
    size_t i = 0;
    for (size_t iz = 0; iz < 110; ++iz) {
      if (zref[iz] <= z && zref[iz + 1] >= z) {
        i = iz;
        break;
      }
    }
    // Get Venus average temperature, pressure, and gas law constant     
    // at upper and lower height indexes                                 
    greal T1 = tref[i];
    greal T2 = tref[i + 1];
    greal p1 = pref[i];
    greal p2 = pref[i + 1];
    greal R1 = pref[i] / (dref[i] * tref[i]);
    greal R2 = pref[i + 1] / (dref[i + 1] * tref[i + 1]);

    // Linear height interpolation on temperature                        
    greal delz = z - zref[i];
    greal dztot = zref[i + 1] - zref[i];
    refTemperature = T1 + (T2 - T1) * delz / dztot;

    // Pressure scale height and vertical pressure interpolation         
    greal HSCALE = dztot / log(p1 / p2);
    refPressure = p1 * exp(-delz / HSCALE);

    // Linear height interpolation for gas constant                      
    greal R = R1 + (R2 - R1) * delz / dztot;

    // Density from perfect gas law
    refDensity = refPressure / (R * refTemperature);
  }
}

//! \copydoc Atmosphere::updatePressureAtSurface()
greal VenusIRA::getPressureAtSurface() const
{
  // Return the reference value of pressure at height=0 (index 0).
  return pref[0];
}

//! \brief Initializes reference values.
//!
//! This routine inializes static data tables of reference atmosphere values.
//!
//! \b Inputs
//! \arg isInitialized Initialization state flag.
//!
//! \retval zref Reference heights.
//! \retval pref Reference pressure data.
//! \retval dref Reference density data.
//! \retval tref Reference temperature data.
void VenusIRA::initializeData()
{
  // Static data only needs to be initialized once.
  if (isInitialized) {
    return;
  }
  isInitialized = true;

  // Resize tables for number of heights.  (subtract 2 for overlaps, low to mid, mid to high)
  size_t heightSize = venusLow.heightSize + venusMid.heightSize + venusHigh.heightSize - 2;
  zref.resize(heightSize);
  pref.resize(heightSize);
  dref.resize(heightSize);
  tref.resize(heightSize);

  // Create a height index.
  size_t zIndex = 0;

  // Save reference values of pressure, density, temperature          
  // Reference values for low altitudes (0-98 km) (0-30 lat data)     
  while (zIndex < venusLow.heightSize - 1) {
    zref[zIndex] = venusLow.zlo[zIndex];
    pref[zIndex] = venusLow.plo[0][zIndex];
    dref[zIndex] = venusLow.dlo[0][zIndex];
    tref[zIndex] = venusLow.tlo[0][zIndex];
    ++zIndex;
  }

  // Reference values for 100 km altitude        
  // Average from low and middle altitude data sets at 100 km.
  // Low data is at 30 degree latitude.
  // Middle altitude data is average of two LST values at 100 km.
  zref[zIndex] = venusMid.zmd[0];
  pref[zIndex] = sqrt(sqrt(venusMid.pmd[0][0] * venusMid.pmd[1][0]) * venusLow.plo[0][venusLow.heightSize - 1]);
  dref[zIndex] = sqrt(sqrt(venusMid.dmd[0][0] * venusMid.dmd[1][0]) * venusLow.dlo[0][venusLow.heightSize - 1]);
  tref[zIndex] = 0.5*(0.5*(venusMid.tmd[0][0] + venusMid.tmd[1][0]) + venusLow.tlo[0][venusLow.heightSize - 1]);
  ++zIndex;

  // Reference values for middle altitudes (105-145 km)       
  // Average of LST = 0 and LST = 12.                                         
  for (size_t i = 1; i < venusMid.heightSize - 1; ++i) {
    zref[zIndex] = venusMid.zmd[i];
    pref[zIndex] = sqrt(venusMid.pmd[0][i] * venusMid.pmd[1][i]);
    dref[zIndex] = sqrt(venusMid.dmd[0][i] * venusMid.dmd[1][i]);
    tref[zIndex] = 0.5 * (venusMid.tmd[0][i] + venusMid.tmd[1][i]);
    ++zIndex;
  }

  // Reference values for 150 km.
  // Average of middle and high altitude data sets at 150 km.
  // Middle altitude data is average of two LST values at 150 km.
  // High altitude data is at 90 degree solar zenith angle.
  zref[zIndex] = venusHigh.zhi[0];
  pref[zIndex] = sqrt(venusHigh.phi[3][0] * sqrt(venusMid.pmd[0][venusMid.heightSize - 1] * venusMid.pmd[1][venusMid.heightSize - 1]));
  dref[zIndex] = sqrt(venusHigh.dhi[3][0] * sqrt(venusMid.dmd[0][venusMid.heightSize - 1] * venusMid.dmd[1][venusMid.heightSize - 1]));
  tref[zIndex] = 0.5 * (venusHigh.thi[3][0] + 0.5 * (venusMid.tmd[0][venusMid.heightSize - 1] + venusMid.tmd[1][venusMid.heightSize - 1]));
  ++zIndex;

  // Reference values for high altitudes (155-250 km).   
  // Solar zenith angle 90 data.                          
  for (size_t i = 1; i < venusHigh.heightSize; ++i) {
    zref[zIndex] = venusHigh.zhi[i];
    pref[zIndex] = venusHigh.phi[3][i];
    dref[zIndex] = venusHigh.dhi[3][i];
    tref[zIndex] = venusHigh.thi[3][i];
    ++zIndex;
  }

}

} // namespace