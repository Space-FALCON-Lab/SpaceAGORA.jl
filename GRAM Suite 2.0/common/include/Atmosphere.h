//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <vector>
#include <map>
#include <functional>
#include "unittest_friend.h"
#include "gram.h"
#include "Position.h"
#include "AtmosphereState.h"
#include "EphemerisState.h"
#include "GramTime.h"
#include "InputParameters.h"

//! \defgroup CommonGRAM The Full List of Common GRAM Framework Classes.
//! \brief Class references for all of the GRAM common framework classes.
//!
//! These are the common framework classes needed for building the GRAM common library.

namespace GRAM {

//! \brief An abstract base class for all GRAM atmosphere models.
//!
//! This defines an interface class for all GRAM atmosphere models. Typically, each planet will
//! implement an interface class with planet specific functionality.  Each model
//! for the planet would then inherit the planet atmosphere interface.
//! \ingroup CommonGRAM
class Atmosphere
{
public:
  Atmosphere();
  Atmosphere(const Atmosphere& orig);
  virtual ~Atmosphere() = default;

  virtual void update() = 0;

  virtual void setPosition(const Position& pos) { position = pos; }
  const Position& getPosition() const { return position; }

  virtual void setEphemerisState(const EphemerisState& state) { ephem = state; }
  const EphemerisState& getEphemerisState() const { return ephem; }

  const AtmosphereState& getAtmosphereState() const { return atmos; }

  void setPlanetaryConstants(greal mu, greal eqrad, greal prad, greal j2, greal per, greal specificHeatRatio);
  greal getGravitationalParameter() const { return mu; }
  greal getEquatorialRadius() const { return equatorialRadius; }
  greal getPolarRadius() const { return polarRadius; }
  greal getJ2() const { return J2; }
  greal getPeriod() const { return period; }
  greal getSpecificHeatRatio() const { return specificHeatRatio; }

  void setGasConstant(GasType gasType, greal amw);
  bool isGasPresent(GasType gasType) const;
  const ConstituentGas& getConstituentGas(GasType gasType) const { return gases.at(gasType).get(); }

  const Atmosphere& operator= (const Atmosphere& rhs);

  static const std::string& getVersionString();
  static greal getArcAngle(greal phi1, greal theta1, greal phi2, greal theta2);

protected:
  virtual greal getRadius(greal lat) const;
  virtual greal getGravity(greal lat, greal r, greal ht) const;
  virtual greal getTotalMass();
  virtual greal getScaleHeightFromDelta(greal deltaHeight, greal lowerHeightValue, greal upperHeightValue);
  virtual void updatePosition();
  virtual void updateTotalNumberDensity();
  virtual void updateMoleFractions();
  virtual void updateMassFractions();
  virtual void updateCompressibilityFactor();
  virtual void updateDeviations();
  virtual void updatePressureAltitude();
  virtual void updateSpecificHeatRatio();
  virtual void updateSpeedOfSound();
  virtual void updateReferenceValues();
  virtual void updatePressureAtSurface();
  virtual void updateMetrics();
  virtual void updateWinds();


  Position position;               //!< Set the position prior to update.
  EphemerisState ephem;            //!< Can be set prior to update or computed internally.
  AtmosphereState atmos;           //!< Values to be computed during an update.

  // Planetary constants
  greal mu = 0;                    //!< The standard gravitational parameter \f$ \mu = G M\f$ (\f$m^3/s^2\f$).
  greal equatorialRadius = 0;      //!< Equatorial radius (km).
  greal polarRadius = 0;           //!< Polar radius (km).
  greal J2 = 0;                    //!< First gravity field non-spherical term (unitless).
  greal period = 0;                //!< Rotation period (sec) (negative for retrograde).

  // The variables below are references into the Position, EphemerisState, and
  // AtmosphereState structures.  They are provided as a convenience to the model
  // developers.  As an example, position.height may be referenced more simply as height.

  // Convenience variables (Position)
  greal& height         = position.height;          //!< \copydoc Position::height
  greal& latitude       = position.latitude;        //!< \copydoc Position::latitude
  greal& longitude      = position.longitude;       //!< \copydoc Position::longitude
  greal& elapsedTime    = position.elapsedTime;     //!< \copydoc Position::elapsedTime
  bool& isPlanetoCentric    = position.isPlanetoCentric;    //!< \copydoc Position::isPlanetoCentric

  greal& totalRadius    = position.totalRadius;     //!< \copydoc Position::totalRadius
  greal& latitudeRadius = position.latitudeRadius;  //!< \copydoc Position::latitudeRadius
  greal& surfaceHeight  = position.surfaceHeight;   //!< \copydoc Position::surfaceHeight
  greal& gravity        = position.gravity;         //!< \copydoc Position::gravity

  // Convenience variables (EphemerisState)
  greal& solarTime              = ephem.solarTime;               //!< \copydoc EphemerisState::solarTime
  greal& solarHourAngle         = ephem.solarHourAngle;          //!< \copydoc EphemerisState::solarHourAngle
  greal& longitudeSun           = ephem.longitudeSun;            //!< \copydoc EphemerisState::longitudeSun
  greal& subsolarLatitude       = ephem.subsolarLatitude;        //!< \copydoc EphemerisState::subsolarLatitude
  greal& subsolarLongitude      = ephem.subsolarLongitude;       //!< \copydoc EphemerisState::subsolarLongitude
  greal& solarDeclination       = ephem.solarDeclination;        //!< \copydoc EphemerisState::solarDeclination
  greal& solarRightAscension    = ephem.solarRightAscension;     //!< \copydoc EphemerisState::solarRightAscension
  greal& orbitalRadius          = ephem.orbitalRadius;           //!< \copydoc EphemerisState::orbitalRadius
  greal& oneWayLightTime        = ephem.oneWayLightTime;         //!< \copydoc EphemerisState::oneWayLightTime

  greal& longitudeAscendingNode = ephem.longitudeAscendingNode;  //!< deprecated
  greal& inclination            = ephem.inclination;             //!< deprecated
  greal& argumentPerihelion     = ephem.argumentPerihelion;      //!< deprecated
  greal& semimajorAxis          = ephem.semimajorAxis;           //!< deprecated
  greal& eccentricity           = ephem.eccentricity;            //!< deprecated
  greal& meanAnomaly            = ephem.meanAnomaly;             //!< deprecated

  // Convenience variables (AtmosphereState)
  greal& temperature                = atmos.temperature;                 //!< \copydoc AtmosphereState::temperature
  greal& pressure                   = atmos.pressure;                    //!< \copydoc AtmosphereState::pressure
  greal& density                    = atmos.density;                     //!< \copydoc AtmosphereState::density
  greal& pressureScaleHeight        = atmos.pressureScaleHeight;         //!< \copydoc AtmosphereState::pressureScaleHeight
  greal& densityScaleHeight         = atmos.densityScaleHeight;          //!< \copydoc AtmosphereState::densityScaleHeight
  greal& totalNumberDensity         = atmos.totalNumberDensity;          //!< \copydoc AtmosphereState::totalNumberDensity
  greal& averageMolecularWeight     = atmos.averageMolecularWeight;      //!< \copydoc AtmosphereState::averageMolecularWeight
  greal& compressibilityFactor      = atmos.compressibilityFactor;       //!< \copydoc AtmosphereState::compressibilityFactor
  greal& specificGasConstant        = atmos.specificGasConstant;         //!< \copydoc AtmosphereState::specificGasConstant
  greal& speedOfSound               = atmos.speedOfSound;                //!< \copydoc AtmosphereState::speedOfSound
  greal& specificHeatRatio          = atmos.specificHeatRatio;           //!< \copydoc AtmosphereState::specificHeatRatio 

  // Search Key: All gases
  ConstituentGas& argon             = atmos.argon;                       //!< \copydoc AtmosphereState::argon
  ConstituentGas& helium            = atmos.helium;                      //!< \copydoc AtmosphereState::helium          
  ConstituentGas& dihydrogen        = atmos.dihydrogen;                  //!< \copydoc AtmosphereState::dihydrogen      
  ConstituentGas& hydrogen          = atmos.hydrogen;                    //!< \copydoc AtmosphereState::hydrogen        
  ConstituentGas& methane           = atmos.methane;                     //!< \copydoc AtmosphereState::methane         
  ConstituentGas& dinitrogen        = atmos.dinitrogen;                  //!< \copydoc AtmosphereState::dinitrogen      
  ConstituentGas& nitrogen          = atmos.nitrogen;                    //!< \copydoc AtmosphereState::nitrogen        
  ConstituentGas& oxygen            = atmos.oxygen;                      //!< \copydoc AtmosphereState::oxygen          
  ConstituentGas& dioxygen          = atmos.dioxygen;                    //!< \copydoc AtmosphereState::dioxygen        
  ConstituentGas& carbonMonoxide    = atmos.carbonMonoxide;              //!< \copydoc AtmosphereState::carbonMonoxide  
  ConstituentGas& carbonDioxide     = atmos.carbonDioxide;               //!< \copydoc AtmosphereState::carbonDioxide   
  ConstituentGas& ozone             = atmos.ozone;                       //!< \copydoc AtmosphereState::ozone        
  ConstituentGas& nitrousOxide      = atmos.nitrousOxide;                //!< \copydoc AtmosphereState::nitrousOxide        
  ConstituentGas& water             = atmos.water;                       //!< \copydoc AtmosphereState::water
  std::map<GasType, std::reference_wrapper<ConstituentGas>> gases;       //!< \brief A map of gases by type (for ease of processing).

  greal& nsWind                     = atmos.nsWind;                      //!< \copydoc AtmosphereState::nsWind
  greal& ewWind                     = atmos.ewWind;                      //!< \copydoc AtmosphereState::ewWind
  greal& verticalWind               = atmos.verticalWind;                //!< \copydoc AtmosphereState::verticalWind
  greal& lowDensity                 = atmos.lowDensity;                  //!< \copydoc AtmosphereState::lowDensity
  greal& highDensity                = atmos.highDensity;                 //!< \copydoc AtmosphereState::highDensity
  greal& densityPerturbation        = atmos.densityPerturbation;         //!< \copydoc AtmosphereState::densityPerturbation
  greal& perturbedDensity           = atmos.perturbedDensity;            //!< \copydoc AtmosphereState::perturbedDensity
  greal& densityStandardDeviation   = atmos.densityStandardDeviation;    //!< \copydoc AtmosphereState::densityStandardDeviation
  greal& perturbedSpeedOfSound      = atmos.perturbedSpeedOfSound;       //!< \copydoc AtmosphereState::perturbedSpeedOfSound

  greal& referenceDensity           = atmos.referenceDensity;            //!< \copydoc AtmosphereState::referenceDensity
  greal& referenceTemperature       = atmos.referenceTemperature;        //!< \copydoc AtmosphereState::referenceTemperature
  greal& referencePressure          = atmos.referencePressure;           //!< \copydoc AtmosphereState::referencePressure
  greal& densityDeviation           = atmos.densityDeviation;            //!< \copydoc AtmosphereState::densityDeviation         
  greal& lowDensityDeviation        = atmos.lowDensityDeviation;         //!< \copydoc AtmosphereState::lowDensityDeviation      
  greal& highDensityDeviation       = atmos.highDensityDeviation;        //!< \copydoc AtmosphereState::highDensityDeviation     
  greal& perturbedDensityDeviation  = atmos.perturbedDensityDeviation;   //!< \copydoc AtmosphereState::perturbedDensityDeviation

  greal& sigmaLevel                 = atmos.sigmaLevel;                  //!< \copydoc AtmosphereState::sigmaLevel
  greal& pressureAtSurface          = atmos.pressureAtSurface;           //!< \copydoc AtmosphereState::pressureAtSurface
  greal& pressureAltitude           = atmos.pressureAltitude;            //!< \copydoc AtmosphereState::pressureAltitude

  greal& ewWindPerturbation         = atmos.ewWindPerturbation;          //!< \copydoc AtmosphereState::ewWindPerturbation
  greal& nsWindPerturbation         = atmos.nsWindPerturbation;          //!< \copydoc AtmosphereState::nsWindPerturbation
  greal& verticalWindPerturbation   = atmos.verticalWindPerturbation;    //!< \copydoc AtmosphereState::verticalWindPerturbation
  greal& perturbedEWWind            = atmos.perturbedEWWind;             //!< \copydoc AtmosphereState::perturbedEWWind
  greal& perturbedNSWind            = atmos.perturbedNSWind;             //!< \copydoc AtmosphereState::perturbedNSWind
  greal& perturbedVerticalWind      = atmos.perturbedVerticalWind;       //!< \copydoc AtmosphereState::perturbedVerticalWind
  greal& ewStandardDeviation        = atmos.ewStandardDeviation;         //!< \copydoc AtmosphereState::ewStandardDeviation
  greal& nsStandardDeviation        = atmos.nsStandardDeviation;         //!< \copydoc AtmosphereState::nsStandardDeviation
  greal& verticalStandardDeviation  = atmos.verticalStandardDeviation;   //!< \copydoc AtmosphereState::verticalStandardDeviation
  UpdateStatus& updateStatus        = atmos.updateStatus;                //!< \copydoc AtmosphereState::updateStatus

private:
  void buildGasMap();

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(Atmosphere, getRadius);
  FRIEND_TEST(Atmosphere, getGravity);
  FRIEND_TEST(Atmosphere, getTotalMass);
  FRIEND_TEST(Atmosphere, updateTotalNumberDensity);
  FRIEND_TEST(Atmosphere, updateMoleFractions);
  FRIEND_TEST(Atmosphere, updateMassFractions);
#endif // GRAM_UNIT_TEST
};

} // namespace

