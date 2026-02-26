//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <array>
#include <map>
#include "gram.h"

namespace GRAM {

//! \brief The output state of an Atmosphere model.
//!
//! For a given Position, an Atmosphere model will produce an ConstituentGas as output.
//! \ingroup CommonGRAM Cpp_Venus Cpp_Earth Cpp_Mars Cpp_Jupiter Cpp_Saturn Cpp_Uranus Cpp_Neptune Cpp_Titan
class ConstituentGas
{
public:
  ConstituentGas();
  ConstituentGas(const ConstituentGas& orig) = default;
  virtual ~ConstituentGas() = default;

  void updateSpecificHeatCapacity(GasType type, greal temperature, greal pressure);

  bool isPresent = false;             //!< Set to true when the gas is used in a model.
  greal numberDensity = 0.0;          //!< Number density of the gas.
  greal moleFraction = 0.0;           //!< Molecular mole fraction.
  greal massFraction = 0.0;           //!< Molecular mass fraction.
  greal averageMolecularWeight = 0.0; //!< Average molecular weight.
  greal specificHeatCapacity = 0.0;   //!< Specific heat capacity, Cp, at constant pressure /f$[J/gK]/f$.
  greal specificHeatPressure = 0.0;   //!< Specific heat capacity, Cp, at constant pressure /f$[J/gK]/f$.
  greal specificHeatVolume = 0.0;     //!< Specific heat capacity, Cv, at constant volume /f$[J/gK]/f$.

private:
  bool update1D(GasType type, greal temperature);
  void update2D(GasType type, greal temperature, greal pressure);
  void updateO2(GasType type, greal temperature, greal pressure);
  void updateH2O(GasType type, greal temperature, greal pressure);
  size_t getTemperatureIndex(greal temperature);
  size_t getPressureIndex(greal pressure);

  static const size_t tSize = 20;
  static const size_t pSize = 10;

  static greal tLevel[tSize];
  static greal pLevel[pSize];

  typedef const std::array<greal, tSize> T_Array;
  typedef const std::array<std::array<greal, pSize>, tSize> TP_Array;
  typedef const std::map<GasType, T_Array> T_Map;
  typedef const std::map<GasType, TP_Array > TP_Map;

  static T_Map cpDataByT;
  static T_Map cvDataByT;

  static TP_Map cpDataByTP;
  static TP_Map cvDataByTP;

// The Shomate method has been deprecated.  
#ifdef OLDCP
public:
  static greal getSpecificHeatCapacity(GasType type, greal temperature);
private:
  struct Shomate {
    // These are the traditional names for the coefficeints.  
    // cp = A + Bt + Ct^2 + Dt^3 + Et^-2
    greal A, B, C, D, E;
  };

  const static Shomate coefficents[static_cast<size_t>(WATER) + 1][4];  //!< The Shomate spline coefficients.
  const static greal bounds[static_cast<size_t>(WATER) + 1][5];       //!< The temperature bounds for the splines.
#endif

};

} // namespace