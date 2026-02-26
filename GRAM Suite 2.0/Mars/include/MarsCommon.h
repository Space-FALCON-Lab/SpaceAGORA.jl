//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"
#include "Atmosphere.h"

namespace GRAM {

enum MarsOffsetModel {
  MARS_CONSTANT,        //!< Use constant input height offset.
  MARS_SEASONAL,        //!< Use constant input with LS seasonal effect.
  MARS_GLOBAL_MEAN,     //!< Use global mean height offset.
  MARS_DAILY_AVERAGE,   //!< Use daily average height offset.
  MARS_CURRENT          //!< Use height offset at current time and position.
};

typedef double(*TopoCallback)(double, double, void*);

//! \brief Common values for the Mars atmosphere models.
//!
//! Common values for the Mars atmosphere models. All Mars
//! atmosphere models should inherit this class.
//!
//! \anchor MC1 (1) Data from NSSDCA Planetary Fact Sheet for Mars (https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html).
//! Page last updated on Sept 27, 2018.
//!
//! \anchor MC2 (2) Data from NAIF Spice files (mar098.cmt).
class MarsCommon
{
public:
  MarsCommon(Atmosphere* atmos);
  virtual ~MarsCommon() = default;

protected:
  // from atmos2.f
  const greal amw = 43.49;  //!<  Constant average molecular weight used in the homosphere.

  static std::string dataPath; //!< Path to the folder containing the binary data files.
  const std::string molaFileName = "MOLA_data.bin";                 //!< Default data file name.
  const std::string mgcmLowerFileName = "MGCM_lower_data.bin";      //!< Default data file name.
  const std::string mgcmUpperFileName = "MGCM_upper_data.bin";      //!< Default data file name.
  const std::string mgcmSurfaceFileName = "MGCM_surface_data.bin";  //!< Default data file name.
  const std::string tesLowerFileName = "TES_lower_data.bin";        //!< Default data file name.
  const std::string tesUpperFileName = "TES_upper_data.bin";        //!< Default data file name.
  const std::string tesSurfaceFileName = "TES_surface_data.bin";    //!< Default data file name.
  const std::string tesDustFileName = "TES_dust_data.bin";          //!< Default data file name.

private:
 
  Atmosphere* atmosphere;                             //!< Pointer to the parent atmosphere.
  const greal muDefault = 4.282837362069909e13;       //!< Gravitational constant \units{m^3/s^2}.  \ref MC2 "(2)"
  const greal J2Default = 0.00196045;                 //!< Gravitational coefficient.  \ref MC1 "(1)"
  const greal periodDefault = 88642.44;               //!< Sidereal rotation period in seconds.  \ref MC1 "(1)"
  const greal equatorialRadiusDefault = 3396.2;       //!< \units{km}.  \ref MC1 "(1)"
  const greal polarRadiusDefault = 3376.20;           //!< \units{km}.  \ref MC1 "(1)"

  // Assume specific heat ratio = 4 / 3  (datastep.f DSTP417)
  const greal specificHeatRatioDefault = 4.0 / 3.0;   //!< The default values for the specific heat ratio.

  // from species.f
  // Molecular weights computed using isotope ratios applicable for
  // Mars atmosphere. See Kieffer et al., editors, "Mars" (1992) and
  // Table A - 6 of NASA / TM - 2001 - 210935 (2001)
  const greal amwDinitrogen = 28.01781;     //!< Average molecular weight of Nitrogen diatom (N2) on Mars.
  const greal amwArgon = 39.96093;          //!< Average molecular weight of Argon (AR) on Mars.
  const greal amwDioxygen = 31.99799;       //!< Average molecular weight of Oxygen diatom (O2) on Mars.
  const greal amwCarbonDioxide = 44.00903;  //!< Average molecular weight of Carbon Dioxide (CO2) on Mars.
  const greal amwCarbonMonoxide = 28.01003; //!< Average molecular weight of Carbon Monoxide (CO) on Mars.
  const greal amwOxygen = 15.99900;         //!< Average molecular weight of Oxygen atom (O) on Mars.
  const greal amwHelium = 4.00260;          //!< Average molecular weight of helium (HE) on Mars.
  const greal amwDihydrogen = 2.01746;      //!< Average molecular weight of hydrogen diatom (H2) on Mars.
  const greal amwHydrogen = 1.00873;        //!< Average molecular weight of hydrogen atom (H) on Mars.
  const greal amwWater = 18.01646;          //!< Average molecular weight of water vapor (H2O) on Mars.

};

} // namespace

