//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Earth-GRAM
// Point of Contact: P. White
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"
#include "Atmosphere.h"

namespace GRAM {

//! \brief Common values for the Earth atmosphere models.
//!
//! Common values for the Earth atmosphere models. All Earth
//! atmosphere models should inherit this class.
//!
//! \anchor EC1 (1) Data from NGA.STND.0036_1.0.0_WGS84 updated 7/8/2014.
//!
//! \anchor EC2 (2) Data from NAIF Spice files (gm_de431.tpc).
//! \ingroup EarthGRAM
class EarthCommon
{
public:
  EarthCommon(Atmosphere* atmos);
  virtual ~EarthCommon() = default;

  greal dedt(greal t, greal p);
  greal d2edt2(greal t, greal p);
  greal wexler(greal t, greal p);
  void caltojul(int iy, int im, int id, int ihour, int imin, double sec, double &xjd);
  int getNextInt(const std::string&, size_t& pos);
  void fair(int i, greal hg, greal pg, greal dg, greal tg, greal ug, greal vg, greal wg,
    greal hj, greal pj, greal dj, greal tj, greal uj, greal vj, greal wj, greal h,
    greal &p, greal &d, greal &t, greal &u, greal &v, greal &w);
  void fair(int i, greal hg, const AtmosphereState& atmosG,
    greal hj, const AtmosphereState& atmosJ, greal h, AtmosphereState& atmosX);
  void getEffectiveRadius(greal lat, greal &g0, greal &re);

protected:
	const greal referenceGravity = 9.80665;  //!< reference value of surface gravity (m / s^2)
  const greal dryAirMolWt = 28.9649;       //!< Average molecular weight of dry air on Earth (360 ppm CO2).
  const greal dryAirGasConstant = 287.055; //!< Gas constant for dry air.

  static std::string dataPath; //!< Path to the folder containing the data files.
  static std::string atmPath;  //!< Path to the folder containing the data files.
  static std::string NCEPPath; //!< Path to the folder containing the data files.
  static std::string M2Path;   //!< Path to the folder containing the data files.
  static std::string rraPath;  //!< Path to the folder containing the data files.

private:
  Atmosphere& atmosphere;                            //!< A reference to the parent Atmosphere object.

  const greal muDefault = 3.9860043543609598E+14;    //!< Gravitational constant \f$m^3/s^2\f$.  \ref EC2 "(2)"
  const greal equatorialRadiusDefault = 6378.137;    //!< \f$km\f$.  \ref EC1 "(1)"
  const greal polarRadiusDefault = 6356.752314;      //!< \f$km\f$.  \ref EC1 "(1)"
  const greal J2Default = 0.00108263;                //!< Gravitational coefficient.   \ref EC1 "(1)"
  const greal periodDefault = 86164.0905;            //!< Sidereal period in seconds.

  // Using hydrogen
  const greal specificHeatRatioDefault = 1.41;       //!< \f$C_P/C_V\f$ (unitless).

  const greal amwNitrogen = 14.0067;        //!< Average molecular weight of Nitrogen atom (N) on Earth.
  const greal amwDinitrogen = 28.0134;      //!< Average molecular weight of Nitrogen diatom (N2) on Earth.
  const greal amwArgon = 39.948;            //!< Average molecular weight of Argon (AR) on Earth.
  const greal amwMethane = 16.043;          //!< Average molecular weight of Methane (CH4) on Earth.
  const greal amwDioxygen = 31.9988;        //!< Average molecular weight of Oxygen diatom (O2) on Earth.
  const greal amwCarbonDioxide = 44.01;     //!< Average molecular weight of Carbon Dioxide (CO2) on Earth.
  const greal amwCarbonMonoxide = 28.01;    //!< Average molecular weight of Carbon Monoxide (CO) on Earth.
  const greal amwOxygen = 15.9994;          //!< Average molecular weight of Oxygen atom (O) on Earth.
  const greal amwHelium = 4.0026;           //!< Average molecular weight of helium (HE) on Earth.
  const greal amwHydrogen = 1.00797;        //!< Average molecular weight of hydrogen atom (H) on Earth.
  const greal amwWater = 18.01528;          //!< Average molecular weight of water vapor (H2O) on Earth.
  const greal amwOzone = 47.9982;           //!< Average molecular weight of Trioxygen (O3) on Earth.
  const greal amwNitrousOxide = 44.0129;    //!< Average molecular weight of Nitrous oxide (O) on Earth.

};

} // namespace

