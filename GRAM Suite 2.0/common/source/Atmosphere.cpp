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
#include <cmath>
#include <algorithm>
#include "Atmosphere.h"

using namespace std;

#if defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

namespace GRAM {

//! \brief Get the current version as a string.
//!
//! \returns The current version as a string.
const std::string& Atmosphere::getVersionString()
{
  static const std::string version = "GRAMLib 2023b :: GRAM Suite 2.0";
  return version;
}

//! \brief The default constructor.
//!
//! This basic constructor is called when objects are created without parameters.
Atmosphere::Atmosphere()
{
  mu = 0.0;
  equatorialRadius = 0.0;
  polarRadius = 0.0;
  J2 = 0.0;
  period = 0.0;
  buildGasMap();
}

//! \brief The copy constructor.
//!
//! This constructor enables the copying of objects.  This can be done by using
//! the constructor: Object newobject(oldobject).  Or by assignment: newobject = oldobject.
//! \param orig The source object to copy.
Atmosphere::Atmosphere(const Atmosphere& orig)
{
  position = orig.position;
  atmos = orig.atmos;
  ephem = orig.ephem;
  mu = orig.mu;
  equatorialRadius = orig.equatorialRadius;
  polarRadius = orig.polarRadius;
  J2 = orig.J2;
  period = orig.period;
  buildGasMap();
}

//! \brief A constructor helper function.
//!
//! This helper function builds a map of gases that helps locate the constituent
//! gas members of this object.
void Atmosphere::buildGasMap()
{
  // Search Key: All gases
  gases.clear();
  gases.emplace(ARGON, argon);
  gases.emplace(CARBON_DIOXIDE, carbonDioxide);
  gases.emplace(CARBON_MONOXIDE, carbonMonoxide);
  gases.emplace(DIHYDROGEN, dihydrogen);
  gases.emplace(DINITROGEN, dinitrogen);
  gases.emplace(DIOXYGEN, dioxygen);
  gases.emplace(HELIUM, helium);
  gases.emplace(HYDROGEN, hydrogen);
  gases.emplace(METHANE, methane);
  gases.emplace(NITROGEN, nitrogen);
  gases.emplace(OXYGEN, oxygen);
  gases.emplace(OZONE, ozone);
  gases.emplace(NITROUS_OXIDE, nitrousOxide);
  gases.emplace(WATER, water);
}


//! \fn  Atmosphere::~Atmosphere()
//! \brief The destructor.
//!
//! The destructor is called when an object is deleted or goes out of scope.

//! \brief The assignment operator.
//!
//! The assignment operator allows the copying of data from one object
//! into another with an assignment statment (LHS = RHS;).
//! \param rhs The right hand side of the assignment operator.
//! \returns A reference to the left hand side object.
const Atmosphere& Atmosphere::operator=(const Atmosphere& rhs)
{
  position = rhs.position;
  atmos = rhs.atmos;
  ephem = rhs.ephem;
  return *this;
}

//! \fn void Atmosphere::update()
//! \brief Interface for the primary atmosphere computations.
//!
//! This is a pure virtual function that defines the interface for the primary
//! atmosphere computations. All atmosphere models must implement this method to
//! populate the AtmosphereState variables.

//! \brief Get the planetocentric radius of the planet at the specified latitude.
//!
//! Relies on planetary constants which must be set in the model.
//! \param lat Latitude in degrees.
//! \returns The radius in km.
greal Atmosphere::getRadius(greal lat) const
{
  // Assumptions (checking occurs in debug mode only)
  assert(polarRadius != 0.0);
  assert(equatorialRadius >= polarRadius);

  // Flattening term
  greal flat = (equatorialRadius - polarRadius) / polarRadius;
  // Planetary radius at latitude
  greal s2lat = square(sin(toRadians(lat)));
  greal r = equatorialRadius / sqrt(1.0 + flat * (flat + 2.0) * s2lat);
  return r;
}

//! \brief Prepares position data including radius and gravity.
//!
//! This routine ensures position data is in geocentric frame and within bounds.
//! It also completes the Position data by computing latitude radius, 
//! total radius, and the force of gravity at the current position.
//!
//! \b Inputs
//! \arg #latitude      
//! \arg #height
//!
//! \retval #latitudeRadius
//! \retval #totalRadius   
//! \retval #gravity       
void Atmosphere::updatePosition()
{
  longitude = wrapDegrees(longitude);
  latitude = clamp(latitude, -90.0_deg, 90.0_deg);
  latitudeRadius = getRadius(latitude);
  totalRadius = latitudeRadius + height;
  gravity = getGravity(latitude, latitudeRadius, height);
}

//! \brief Get the gravity of the planet at the specified latitude and height.
//!
//! Relies on planetary constants which must be set in the model.
//! \param lat Latitude in degrees.
//! \param r   Radius of the planet at lat in km.
//! \param ht  Height in km.
//! \returns The gravity in m/s^2.
greal Atmosphere::getGravity(greal lat, greal r, greal ht) const
{
  // Assumptions (checking occurs in debug mode only)
  assert(period != 0.0);
  assert(r > 0.0);
  assert(r + ht != 0.0);

  // Rotation rate (rad/sec)
  greal omega = TWO_PI / abs(period);
  greal latr = toRadians(lat);
  greal P2 = 1.5 * pow(sin(latr), 2) - 0.5;
  // Gravity at height hgt
  greal gz = 1.0e-6 * mu / pow(r + ht, 2);
  // J2 correction term
  greal gz2 = gz - gz * 3.0 * J2 * P2 * (pow(equatorialRadius / (r + ht), 2));
  // Rotation correction term
  greal g = gz2 - 1.0e3 * (r + ht) * pow(omega * cos(latr), 2);
  return g;
}

greal Atmosphere::getScaleHeightFromDelta(greal deltaHeight, greal lowerHeightValue, greal upperHeightValue)
{
  return deltaHeight / log(lowerHeightValue / upperHeightValue);
}

//! \brief Sets a number of planetary constants.
//!
//! \param mu     The standard gravitational parameter mu = G M (m**3/s**2)
//! \param eqrad  The equatorial radius in km.
//! \param prad   The polar radius in km.
//! \param J2     The First gravity field non-spherical term (unitless).
//! \param period Rotation period (sec) (negative for retrograde).
//! \param specificHeatRatio The adiabatic index, Poisson constant, or the heat capacity ratio.
void Atmosphere::setPlanetaryConstants(greal mu, greal eqrad, greal prad, greal J2, greal period, greal specificHeatRatio)
{
  assert(mu != 0.0);
  assert(eqrad != 0.0);
  assert(prad != 0.0);
  assert(J2 != 0.0);
  assert(period != 0.0);
  assert(specificHeatRatio != 0.0);
  this->mu = mu;
  this->equatorialRadius = eqrad;
  this->polarRadius = prad;
  this->J2 = J2;
  this->period = period;
  this->specificHeatRatio = specificHeatRatio;
}

//! \brief Activates the gas computations and sets the molecular weight.
//!
//! All available gas types are listed in the GasType enumeration.  Only those
//! gases that are given a non-zero molecular weight will be used in computations.
//! Setting a molecular weight to zero will de-activate the computations for that gas.
//! \param gasType The name of the desired gas.
//! \param amw The average molecular weight of the gas.
void Atmosphere::setGasConstant(GasType gasType, greal amw)
{
  // Get the gas element.
  ConstituentGas& gas = gases.at(gasType);
  // Set the molecular weight.
  gas.averageMolecularWeight = amw;
  // Set the activation flag.
  gas.isPresent = (gas.averageMolecularWeight != 0.0);
}

//! \brief Tests for active gas types.
//!
//! Since not all gases are used in some atmosphere models, gases are flagged as active
//! witht the isPresent flag.  This method returns that flag for the specified gas type.
//! \param gasType A GasType to query.
bool Atmosphere::isGasPresent(GasType gasType) const
{
	// Get the gas element.
	ConstituentGas& gas = gases.at(gasType);
	// Return the boolean
	return gas.isPresent;
}


//! \brief Compute the total number density from the constituent gases.
//!
//! The routine adds up the total number densities of all active gases.
//!
//! \b Inputs
//! \arg #gases
//!
//! \retval #totalNumberDensity
void Atmosphere::updateTotalNumberDensity()
{
  totalNumberDensity = 0.0;
  // Iterate over the set of all activated gases.
  // This macro defines "gas" and "gasType" for each activated gas.
  FOR_ALL_GASES_PRESENT(
    totalNumberDensity += gas.numberDensity;
  );
}

//! \brief Compute the total mass from the constituent gases.
//!
//! The routine adds up the total mass using the number density
//! and molecular weight of each of the active gases.
//!
//! \b Inputs
//! \arg #gases
//!
//! \returns The total mass.
greal Atmosphere::getTotalMass()
{
  greal totalMass = 0.0;
  // Iterate over the set of all activated gases.
  // This macro defines "gas" and "gasType" for each activated gas.
  FOR_ALL_GASES_PRESENT(
    totalMass += gas.numberDensity * gas.averageMolecularWeight;
  );
  return totalMass;
}

//! \brief Compute the mole fraction of each constituent gases.
//!
//! The routine computes the mole fraction of each of the active gases.
//!
//! \b Inputs
//! \arg #gases
//!
//! \retval #gases
void Atmosphere::updateMoleFractions()
{
  // Assumptions (checking occurs in debug mode only)
  assert(totalNumberDensity != 0.0);

  // Iterate over the set of all activated gases.
  // This macro defines "gas" and "gasType" for each activated gas.
  FOR_ALL_GASES_PRESENT(
    gas.moleFraction = gas.numberDensity / totalNumberDensity;
  );
}

//! \brief Compute the mass fraction of each constituent gases.
//!
//! The routine computes the mass fraction for each of the active gases.
//!
//! \b Inputs
//! \arg #gases
//!
//! \retval #gases
void Atmosphere::updateMassFractions()
{
  // First get the total mass.
  greal totalMass = getTotalMass();

  // Assumptions (checking occurs in debug mode only)
  assert(totalMass != 0.0);

  // Iterate over the set of all activated gases.
  // This macro defines "gas" and "gasType" for each activated gas.
  FOR_ALL_GASES_PRESENT(
    gas.massFraction = gas.numberDensity * gas.averageMolecularWeight / totalMass;
  );
}

//! \brief Update a metrics that depend on the atmosphere state.
//!
//! This is a convenience method that calls updates for a number of metrics
//! that are derived from the current atmosphere state.
void Atmosphere::updateMetrics()
{
  updateCompressibilityFactor();
  updateReferenceValues();
  updateDeviations();
  updatePressureAtSurface();
  updatePressureAltitude();
}

//! \brief Computes compressibility factor (Z).
//!
//! The compressibility factor measures the deviation from the ideal gas behavior.
//!
//! \b Inputs
//! \arg #pressure   
//! \arg #density    
//! \arg #temperature
//! \arg #averageMolecularWeight 
//!
//! \retval #compressibilityFactor 
void Atmosphere::updateCompressibilityFactor()
{
  // Assumptions (checking occurs in debug mode only)
  assert(density != 0.0);
  assert(UNIVERSAL_GAS != 0.0);
  assert(temperature != 0.0);

  compressibilityFactor = pressure * averageMolecularWeight / (density * UNIVERSAL_GAS * 1.0e3 * temperature);
}

//! \brief Computes the reference atmosphere values.
//!
//! Computes pressure, density, and temperature for the reference atmosphere.
//! The routine will return null values if it has not been overridden.
//! \retval #referenceTemperature
//! \retval #referencePressure   
//! \retval #referenceDensity    
void Atmosphere::updateReferenceValues()
{
  // This base class method should be overridden by the subclass.
  // By default we make the values null.
  referenceTemperature = 0.0;
  referencePressure = 0.0;
  referenceDensity = 0.0;
}

//! \brief Computes the density deviations from the reference atmosphere values.
//!
//! Computes the density deviations of the atmosphere state from the reference atmosphere.
//! The routine will return -99.9 if reference values are undefined.
//!
//! \b Inputs
//! \arg #density          
//! \arg #lowDensity       
//! \arg #highDensity      
//! \arg #perturbedDensity 
//! \arg #referenceDensity
//!
//! \retval #densityDeviation          
//! \retval #lowDensityDeviation       
//! \retval #highDensityDeviation      
//! \retval #perturbedDensityDeviation 
void Atmosphere::updateDeviations()
{
  // If reference values are undefined, so are the deviations.
  if (referenceDensity <= 0.0) {
    densityDeviation = -0.999;
    lowDensityDeviation = -0.999;
    highDensityDeviation = -0.999;
    perturbedDensityDeviation = -0.999;
  }
  else {
    // Compute percent deviations from the average values
    densityDeviation = (density - referenceDensity) / referenceDensity;
    lowDensityDeviation = (lowDensity - referenceDensity) / referenceDensity;
    highDensityDeviation = (highDensity - referenceDensity) / referenceDensity;
    perturbedDensityDeviation = (perturbedDensity - referenceDensity) / referenceDensity;
  }
}

//! \brief Computes the pressure at the suface.
//!
//! This routine computes the pressure at the surface below (or above) the 
//! current position.  If not overridden, this routine will return null values.
//! \retval #pressureAtSurface
void Atmosphere::updatePressureAtSurface()
{
  pressureAtSurface = 0.0;
}

//! \brief Computes the pressure altitude.
//!
//! This routine computes the sigma level and the pressure altitude at the 
//! current position. 
//!
//! \b Inputs
//! \arg #pressure          
//! \arg #pressureAtSurface 
//!
//! \retval #sigmaLevel        
//! \retval #pressureAltitude  
void Atmosphere::updatePressureAltitude()
{
  if (pressureAtSurface != 0.0) {
    sigmaLevel = pressure / pressureAtSurface;
    pressureAltitude = -pressureScaleHeight * log(sigmaLevel);
  }
  else {
    // If pressureAtSurface is null, use null values here.
    sigmaLevel = 0.0;
    pressureAltitude = 0.0;
  }
}

//! \brief Computes the specific heat ratio.
//!
//! This routine computes the specific heat ratio of the mixture of 
//! constituent gases at the current temperature. 
//!
//! \b Inputs
//! \arg #temperature          
//! \arg #gases 
//!
//! \retval #specificHeatRatio        
void Atmosphere::updateSpecificHeatRatio()
{
  greal cpSum = 0.0; // Isobaric heat capacity of the gas mixture.
  greal cvSum = 0.0; // Isochoric heat capacity of the gas mixture.

  // Iterate over the set of all activated gases.
  // This macro defines "gas" and "gasType" for each activated gas.
#ifndef OLDCP
  FOR_ALL_GASES_PRESENT(
    gas.updateSpecificHeatCapacity(gasType, temperature, pressure);
    cpSum += gas.specificHeatPressure * gas.massFraction;
    cvSum += gas.specificHeatVolume * gas.massFraction;
  );
#else
// The Shomate method has been deprecated.  
  FOR_ALL_GASES_PRESENT(
    // Get Cp for the gas and convert from (J/mol K) to (J/g K).
    gas.specificHeatCapacity = gas.getSpecificHeatCapacity(gasType, temperature) / gas.averageMolecularWeight;
    // Weight the heat capacity by the mass fraction for the mixture.
    cpSum += gas.specificHeatCapacity * gas.massFraction;

    // The isochoric heat capacity is derived from the isobaric
    greal cv = gas.specificHeatCapacity - UNIVERSAL_GAS / gas.averageMolecularWeight;
    // Weight the heat capacity by the mass fraction for the mixture.
    cvSum += cv * gas.massFraction;
    );
#endif

  // Compute the specific heat ratio whenever summations above are non-trivial.
  // Otherwise, the default constant is used.
  if (cpSum > 0.0 && cvSum > 0.0) {
    specificHeatRatio = cpSum / cvSum;
  }
}

//! \brief Computes the speed of sound.
//!
//! The routine computes the speed of sound given the current atmosphere state.
//!
//! \b Inputs
//! \arg #pressure          
//! \arg #density           
//! \arg #specificHeatRatio 
//!
//! \retval #speedOfSound     
void Atmosphere::updateSpeedOfSound()
{
  // Assumptions (checking occurs in debug mode only)
  assert(density != 0.0);

  speedOfSound = sqrt(specificHeatRatio * pressure / density);
}

//! \brief Update the winds component.
//!
//! Winds are typically updated in the atmosphere model classes.
void Atmosphere::updateWinds()
{
  ewWind = 0.0;
  nsWind = 0.0;
  verticalWind = 0.0;
}

//! Returns great-circle distance (degrees) between two input lat-lon positions  (in degrees)
greal Atmosphere::getArcAngle(greal phi1, greal theta1, greal phi2, greal theta2)
{
  greal r1[3], r2[3], r1xr2[3];
  greal phi1r = toRadians(phi1);
  greal theta1r = toRadians(theta1);
  greal phi2r = toRadians(phi2);
  greal theta2r = toRadians(theta2);
  greal cphi1 = cos(phi1r);
  greal cphi2 = cos(phi2r);

  // Get components of unit-magnitude vector toward 1st lat-lon
  r1[0] = cphi1 * cos(theta1r);
  r1[1] = cphi1 * sin(theta1r);
  r1[2] = sin(phi1r);

  // Get components of unit-magnitude vector toward 2nd lat-lon
  r2[0] = cphi2 * cos(theta2r);
  r2[1] = cphi2 * sin(theta2r);
  r2[2] = sin(phi2r);

  // Get cross product vector components from these two unit vectors
  r1xr2[0] = (r1[1] * r2[2]) - (r1[2] * r2[1]);
  r1xr2[1] = (r1[2] * r2[0]) - (r1[0] * r2[2]);
  r1xr2[2] = (r1[0] * r2[1]) - (r1[1] * r2[0]);

  // Get magnitude of cross product vector (Sine of great-circle distance)
  greal r1xr2mag = sqrt(square(r1xr2[0]) + square(r1xr2[1]) + square(r1xr2[2]));

  // Get great-circle distance from cross-product magnitude
  greal arcLength;
  if (r1xr2mag >= 1.0) {
    arcLength = 90.0;
  }
  else {
    arcLength = toDegrees(asin(r1xr2mag));
  }

  if ((r1[0] * r2[0] + r1[1] * r2[1] + r1[2] * r2[2]) < 0.0) {
    arcLength = 180.0 - arcLength;
  }

  return arcLength;
}

//! \fn void Atmosphere::setPosition(const Position& pos)
//! \brief Sets the position prior to update.
//!
//! Before update() is called, the Position information must be provided.
//! The fields height, latitude, longitude, elapsedTime, and isPlanetoCentric must be set.
//! \param pos The new position. The fields height, latitude, longitude, elapsedTime, and isPlanetoCentric must be set.

//! \fn const Position& Atmosphere::getPosition()
//! \brief Get the atmosphere's current position.
//! \returns A Position structure.

//! \fn void Atmosphere::setEphemerisState(const EphemerisState& state)
//! \brief Optional input prior to an update.
//!
//! In order to bypass internal ephemeris calculations, this method must be called
//! prior to each call to update().  See model documentation for details about
//! which ephemeris values must be provided.
//! \param state The desired ephemeris state.

//! \fn const EphemerisState& Atmosphere::getEphemerisState()
//! \brief Get the atmosphere's current ephemeris values.
//!
//! The ephemeris values, after an update, will either be the values provided by
//! a call to setEphemerisState() or those computed by the internal ephemeris model.

//! \fn const AtmosphereState& Atmosphere::getAtmosphereState()
//! \brief Get the current atmosphere values after on update.
//!
//! After an update, this call returns the atmospheric values corresponding to the
//! position provided in setPosition().

} // namespace
