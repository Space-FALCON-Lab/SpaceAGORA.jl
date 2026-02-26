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
#include "Many.h"
#include "ConstituentGas.h"

namespace GRAM {

constexpr greal MIN_MAX_NOT_USED = -999.9;

//! \brief The output state of an Atmosphere model.
//!
//! For a given Position, an Atmosphere model will produce an AtmosphereState as output.
//! \ingroup CommonGRAM
class AtmosphereStateBase
{
public:
  AtmosphereStateBase() = default;
  AtmosphereStateBase(const AtmosphereStateBase& orig) = default;
  virtual ~AtmosphereStateBase() = default;

  // Dynamics
  greal temperature = 0.0;              //!< \units{K}.
  greal pressure = 0.0;                 //!< \units{Pa}.
  greal density = 0.0;                  //!< \units{kg/m^3}.
  greal pressureScaleHeight = 0.0;      //!< \units{km}.
  greal densityScaleHeight = 0.0;       //!< \units{km}.
  greal speedOfSound = 0.0;             //!< \units{m/s}.
  greal referenceDensity = 0.0;         //!< \units{kg/m^3}.
  greal referenceTemperature = 0.0;     //!< \units{K}.
  greal referencePressure = 0.0;        //!< \units{Pa}.
  greal sigmaLevel = 0.0;               //!< Ratio of pressure to surface pressure (unitless).
  greal pressureAtSurface = 0.0;        //!< \units{Pa}.
  greal pressureAltitude = 0.0;         //!< \units{km}.

  greal profileWeight[10] = { 0.0 };    //!< Auxiliary atmosphere weighting factor.

  // Gases
  ConstituentGas argon;                 //!< Ar.
  ConstituentGas carbonDioxide;         //!< CO<sub>2</sub>.
  ConstituentGas carbonMonoxide;        //!< CO.
  ConstituentGas dihydrogen;            //!< H<sub>2</sub> (Molecular hydrogen).
  ConstituentGas dinitrogen;            //!< N<sub>2</sub> (Molecular nitrogen).
  ConstituentGas dioxygen;              //!< O<sub>2</sub> (Molecular oxygen).
  ConstituentGas helium;                //!< He.
  ConstituentGas hydrogen;              //!< H (Atomic hydrogen).
  ConstituentGas methane;               //!< CH<sub>4</sub>.
  ConstituentGas nitrogen;              //!< N (Atomic nitrogen).
  ConstituentGas oxygen;                //!< O (Atomic oxygen).
  ConstituentGas ozone;                 //!< O<sub>3</sub> (Triatomic oxygen).
  ConstituentGas nitrousOxide;          //!< N<sub>2</sub>O (Dinitrogen monoxide).
  ConstituentGas water;                 //!< H<sub>2</sub>O (Water vapor).
  greal totalNumberDensity = 0.0;       //!< Particles per \units{m^3}
  greal averageMolecularWeight = 0.0;   //!< Or average molecular mass (unitless).
  greal compressibilityFactor = 1.0;    //!< Z (unitless).
  greal specificGasConstant = 0.0;      //!< \units{J/(kg K)}.
  greal specificHeatRatio = 0.0;        //!< The adiabatic index, Poisson constant, or the heat capacity ratio (unitless).

  // Density
  greal lowDensity = 0.0;                    //!< \units{kg/m^3}.
  greal highDensity = 0.0;                   //!< \units{kg/m^3}.
  greal densityPerturbation = 0.0;           //!< \units{\%}.   (was densp, drh)
  greal perturbedDensity = 0.0;              //!< \units{kg/m^3}.  (was denstot, dmpert)
  greal densityStandardDeviation = 0.0;      //!< \units{\%}.   (was sigd, sdh)
  greal perturbedSpeedOfSound = 0.0;         //!< \units{m/s}.
  greal relativeStepSize = 0.0;              //!< (unitless).
  greal densityDeviation = 0.0;              //!< \units{\%}.   (was devav, dghp)
  greal lowDensityDeviation = 0.0;           //!< \units{\%}.   (was devlo)
  greal highDensityDeviation = 0.0;          //!< \units{\%}.   (was devhi)
  greal perturbedDensityDeviation = 0.0;     //!< \units{\%}.   (was devtot, dhp)
  UpdateStatus updateStatus = INITIAL_STATE; //!< Status of the most recent update.

  // Winds
  greal ewWind = 0.0;                    //!< East/west winds \units{m/s}. 
  greal nsWind = 0.0;                    //!< North/south winds \units{m/s}.
  greal verticalWind = 0.0;              //!< Vertical winds \units{m/s}.
  greal ewWindPerturbation = 0.0;        //!< \units{m/s}. (was ewpert, urh)
  greal nsWindPerturbation = 0.0;        //!< \units{m/s}. (was nspert, vrh)
  greal verticalWindPerturbation = 0.0;  //!< \units{m/s}. (was wrh)
  greal perturbedEWWind = 0.0;           //!< \units{m/s}. (was ewtot, umpert)
  greal perturbedNSWind = 0.0;           //!< \units{m/s}. (was nstot, vmpert)
  greal perturbedVerticalWind = 0.0;     //!< \units{m/s}. (was wmpert)
  greal ewStandardDeviation = 0.0;       //!< \units{m/s}. (was sigu, suh)
  greal nsStandardDeviation = 0.0;       //!< \units{m/s}. (was sigv, svh)
  greal verticalStandardDeviation = 0.0; //!< \units{m/s}. (was sigw, swh)

  // Special cases
  greal minMaxFactor = MIN_MAX_NOT_USED;      //!< Interpolation factor for the min-max model.
  greal pressureStandardDeviation = 0.0;      //!< \units{\%}. (was sph)
  greal temperatureStandardDeviation = 0.0;   //!< \units{\%}. (was sth)
};


//! \brief The output state of an Atmosphere model.
//!
//! For a given Position, an Atmosphere model will produce an AtmosphereState as output.
//! \note Most of the public attributes are inherited from AtmosphereStateBase.  The class makes it possible
//! to extend the common AtmosphereStateBase with planetary specific metrics (e.g. dust on Mars).
//! \ingroup CommonGRAM Cpp_Venus Cpp_Earth Cpp_Mars Cpp_Jupiter Cpp_Saturn Cpp_Uranus Cpp_Neptune Cpp_Titan
class AtmosphereState : public AtmosphereStateBase
{
public:
  AtmosphereState();
  AtmosphereState(const AtmosphereState& orig);
  virtual ~AtmosphereState();
  AtmosphereState& operator=(const AtmosphereState &value);

  Many planetSpecificMetrics;  //!< A variant memory container to hold planet specific information.

  template<typename T>
  void setPlanetSpecificMetrics(T &value);

  template<typename T>
  const T& getPlanetSpecificMetrics() const;

  template<typename T>
  T& getPlanetSpecificMetrics();

};

//! \brief Set the value of the planet specific metrics.
//!
//! This method sets the planet specific metrics to point at the \p value supplied as an arguement. The value will typically
//! be a planet's atmosphere state, such as MarsAtmosphereState or EarthAtmosphereState.  No memory is copied
//! by this method.  It only saves the reference to \p value.  If \p value is modified external to this object, then
//! it will also appear modified within this object.  If the \p value is deleted externally, then there could be 
//! some trouble.  So don't do that.
//! \param value The source value to set.
//! \tparam T The type of the source value.  This type does not need to be specified since it can be implied from the
//!           type of the value arguement.
template<typename T>
void AtmosphereState::setPlanetSpecificMetrics(T &value)
{
  planetSpecificMetrics.set<T>(value);
}

//! \brief Get a constant reference to the value of the planet specific metrics.
//!
//! This method returns a type appropriate reference to the planet specific metrics.  Since this class does not know the type
//! of the planet specific metrics, the get template must specify this type (e.g. \p get<EarthAtmosphereState>()).
//! \code
//! Many manyObject;
//! EarthAtmosphereState eState;
//! manyObject.set(eState);
//! // Now get the value from the variant memory.
//! const EarthAtmosphereState& eAtmos = manyObject.get<EarthAtmosphereState>();
//! \endcode
//! The result of the code above is that eState and eAtmos are equivalent.
//! Note that the reference returned by this method is constant (non-modifiable).
//! \tparam T The type of the planet specific metrics.
template<typename T>
const T& AtmosphereState::getPlanetSpecificMetrics() const
{
  return planetSpecificMetrics.get<T>();
}

//! \brief Get a reference to the value of the planet specific metrics.
//!
//! This method returns a type appropriate reference to the planet specific metrics.  Since this class does not know the type
//! of the planet specific metrics, the get template must specify this type (e.g. \p get<EarthAtmosphereState>()).
//! \code
//! Many manyObject;
//! EarthAtmosphereState eState;
//! manyObject.set(eState);
//! // Now get the value from the variant memory.
//! EarthAtmosphereState& eAtmos = manyObject.get<EarthAtmosphereState>();
//! \endcode
//! The result of the code above is that eState and eAtmos are equivalent.
//! Note that the reference returned by this method is non-constant (modifiable).
//! \tparam T The type of the planet specific metrics.
template<typename T>
T& AtmosphereState::getPlanetSpecificMetrics()
{
  return planetSpecificMetrics.get<T>();
}

} // namespace