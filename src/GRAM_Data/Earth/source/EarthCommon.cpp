//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Earth-GRAM
// Point of Contact: P. White
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cctype>
#include "EarthCommon.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {

std::string EarthCommon::dataPath;
std::string EarthCommon::atmPath;
std::string EarthCommon::NCEPPath;
std::string EarthCommon::M2Path;
std::string EarthCommon::rraPath;

//! \brief The default constructor
//!
//! This basic constructor is called when objects are created without parameters.
EarthCommon::EarthCommon(Atmosphere* atmos)
  : atmosphere(*atmos)
{
  if (atmos != nullptr) {
    atmosphere.setPlanetaryConstants(muDefault, equatorialRadiusDefault, polarRadiusDefault,
      J2Default, periodDefault, specificHeatRatioDefault);

    atmosphere.setGasConstant(NITROGEN, amwNitrogen);
    atmosphere.setGasConstant(DINITROGEN, amwDinitrogen);
    atmosphere.setGasConstant(ARGON, amwArgon);
    atmosphere.setGasConstant(METHANE, amwMethane);
    atmosphere.setGasConstant(DIOXYGEN, amwDioxygen);
    atmosphere.setGasConstant(CARBON_DIOXIDE, amwCarbonDioxide);
    atmosphere.setGasConstant(CARBON_MONOXIDE, amwCarbonMonoxide);
    atmosphere.setGasConstant(OXYGEN, amwOxygen);
    atmosphere.setGasConstant(HELIUM, amwHelium);
    atmosphere.setGasConstant(HYDROGEN, amwHydrogen);
    atmosphere.setGasConstant(WATER, amwWater);
    atmosphere.setGasConstant(OZONE, amwOzone);
    atmosphere.setGasConstant(NITROUS_OXIDE, amwNitrousOxide);
  }
}

//! \fn  EarthCommon::EarthCommon(const EarthCommon& orig)
//! \brief The copy constructor
//!
//! This constructor enables the copying of objects.  This can be done by using
//! the constructor: Object newobject(oldobject).  Or by assignment: newobject = oldobject.

//! \fn  EarthCommon::~EarthCommon()
//! \brief Destructor.
//!
//! The destructor is called when an object is deleted or goes out of scope.

//! \brief First derivative of saturation vapor pressure
//!
//! Wexler formulation for the first derivative of saturation vapor 
//! pressure (Pa/K) as function of temperature T (kelvin), and pressure
//! (Pa), as given by Flatau et al., J. Appl. Meteorol., 31(12), 1507, 
//! Dec., 1992.
//! \param t Temperature (K)
//! \param p Pressure (Pa)
//! \returns The first derivative of saturation vapor pressure
greal EarthCommon::dedt(greal t, greal p)
{

  constexpr greal g0 = -0.29912729E4;
  constexpr greal g1 = -0.60170128E4;
  constexpr greal g3 = -0.028354721;
  constexpr greal g4 = 0.17838301E-4;
  constexpr greal g5 = -0.84150417E-9;
  constexpr greal g6 = 0.44412543E-12;
  constexpr greal g7 = 2.858487;

  greal esat = wexler(t, p);

  greal dedt_out;
  if (esat <= 0.0) {
    dedt_out = 0.0;
  }
  else {
    dedt_out = esat
               * (((g7 + (g3 + (2.0 * g4 + (3.0 * g5 + 4.0 * g6 * t) * t) * t) * t) * t - g1) * t
                  - 2.0 * g0)
               / pow(t, 3);
  }

  return dedt_out;
}

//! \brief Second derivative of saturation vapor pressure
//!
//! Wexler formulation for the 2nd derivative of saturation vapor pressure
//! (Pa/K^2) as a function of temperature T (kelvin), and pressure (Pa)
//! as given by the method of Flatau et al., J. Appl. Meteorol., 31(12),
//! 1507, Dec., 1992.
//! \param t Temperature (K)
//! \param p Pressure (Pa)
//! \returns The second derivative of saturation vapor pressure
greal EarthCommon::d2edt2(greal t, greal p)
{
  constexpr greal g0 = -0.29912729e4;
  constexpr greal g1 = -0.60170128e4;
  constexpr greal g4 = 0.17838301e-4;
  constexpr greal g5 = -0.84150417e-9;
  constexpr greal g6 = 0.44412543e-12;
  constexpr greal g7 = 2.858487;

  greal esat = wexler(t, p);
  greal fde = dedt(t, p);

  greal d2edt2_out;
  if (esat <= 0.0) {
    d2edt2_out = 0.0;
  }
  else {
    d2edt2_out= (esat
                 * (((-g7 + ((2.0 * g4 + (6.0 * g5 + 12.0 * g6 * t) * t) * t) * t) * t + 2.0 * g1) * t
                    + 6.0 * g0)
                 / pow(t, 4))
                + (pow(fde, 2) / esat);
  }

  return d2edt2_out;
}


//! \brief Saturation vapor pressure
//!
//! Wexler formulation for saturation vapor pressure (Pa) as a function
//! of temperature tk (kelvin) and pressure pres (Pa), as given by Flatau
//! et al., J. Appl. Meteorol., 31(12), 1507, Dec., 1992, with 
//! correction factor (f3) as given by Buck, J. Appl. Meteorol., 
//! 20(12), 1527, Dec. 1981.
//! wexler(T dry bulb) gives saturation vapor pressure.
//! wexler(T dew point) gives actual vapor pressure.
//! Relative humidity (0-1) is wexler(Td)/wexler(T)
//! \param tk Temperature (K)
//! \param pres Pressure (Pa)
//! \returns The saturation vapor pressure
greal EarthCommon::wexler(greal tk, greal pres)
{
  constexpr greal g0 = -0.29912729e4;
  constexpr greal g1 = -0.60170128e4;
  constexpr greal g2 = 18.87643854e0;
  constexpr greal g3 = -0.028354721e0;
  constexpr greal g4 = 0.17838301e-4;
  constexpr greal g5 = -0.84150417e-9;
  constexpr greal g6 = 0.44412543e-12;
  constexpr greal g7 = 2.858487e0;
  constexpr greal a = 1.0007e0;
  constexpr greal b = 3.46e-08;
  greal wexler_out;

  if (tk <= 75.0) {
    wexler_out = 1.0e-23;
  }
  else {
    wexler_out = exp((g0
               + (g1 + (g2 + g7 * log(tk) + (g3 + (g4 + (g5 + g6 * tk) * tk) * tk) * tk) * tk) * tk)
              / pow(tk, 2))
          * (a + b * pres);
  }

  return wexler_out;
}

//! \brief Effective radius and surface gravity from latitude.
//!
//! Computes surface gravity g0 (m/s^2) and effective Earth radius re (km) from 
//! input geocentric latitude phi (degrees).  
//! \param lat Latitude (degrees)
//! \param[out] g0 Surface gravity (m/s^2)
//! \param[out] re Effective radius (km)
void EarthCommon::getEffectiveRadius(greal lat, greal &g0, greal &re)
{
  const greal g0a = 9.80616;
  const greal g0b = 0.0026373;
  const greal g0c = 0.0000059;
  const greal rea = 3.085462e-3;
  const greal reb = 2.27e-6;
  const greal rec = 2.0e-9;

  //Parameters for computing effective Earth radius
  //Set up cosines:  cphi2 = [cos(phir)]**2, ...
  greal cphi2 = square(cos(toRadians(lat)));
  greal c2phi = 2.0 * cphi2 - 1.0;
  greal c4phi = 8.0 * cphi2 * (cphi2 - 1.0) + 1.0;

  //Compute surface gravity
  g0 = g0a * (1.0 - g0b * c2phi + g0c * square(c2phi));

  //Compute effective Earth radius, re
  re = 2.0 * g0 / (rea + reb * c2phi - rec * c4phi);
}

//! \brief Gregorian to Julian date conversion.
//!
//! Compute Julian day by method of Meeus, Astronomical Algorithms, 
//! 2nd Edition, 1998, page 61. 
//! \param iy Year.
//! \param im Month.
//! \param id Day.
//! \param ihour Hour.
//! \param imin Minute.
//! \param sec Seconds.
//! \param[out] xjd Julian date.
void EarthCommon::caltojul(int iy, int im, int id, int ihour, int imin, double sec, double &xjd)
{

  int y, m, a, b;
  double d;

  y = iy;
  m = im;

  //Consider Jan or Feb as if months 13 and 14 of previous year
  if (im <= 2) {
    y = iy - 1;
    m = im + 12;
  }

  //Compute day of month plus fractional part
  d = id + (ihour / 24) + (imin / 1440) + (sec / 86400);
  a = int(0.01*y);
  b = 2 - a + int(0.25*a);

  //Compute julian day with fractional part
  xjd = int(365.25*(y + 4716)) + int(30.6001*(m + 1)) + d + b - 1524.5;
}

//! \brief Parse an integer from a string.
//!
//! Parses the first integer found after a specified position from a string.
//! \param data The string to be parsed.
//! \param[in,out] pos <tt>[in]</tt> The starting position within the string. <br><tt>[out]</tt> The location just after the parsed integer.
//! \returns The parsed integer.
int EarthCommon::getNextInt(const std::string& data, size_t& pos) {
  int value = 0;
  bool negate = false;
  pos = data.find_first_not_of(" ", pos);
  if (data[pos] == '-') {
    negate = true;
    ++pos;
  }
  while (pos < data.size() && std::isdigit(data[pos])) {
    value = 10 * value + (int)(data[pos] - '0');
    ++pos;
  }
  if (negate) {
    value *= -1;
  }
  return value;
}

//! \brief Fair between two atmosphere states.
//!
//! Fairs at height h, between pg,dg,tg,ug,vg,wg and pj,dj,tj,uj,vj,wj.  Fairing is cosine**2,
//! such that all "g" values obtain at height hg and all "j" values obtain at height hj.  The
//! faired results are p,d,t,u,v,w.   
//! \param i 1 means density (d) is faired and pressure (p) is computed from a faired 
//!          value of the gas-law constant.  If i not equal to 1, both p and d are faired directly.
//! \param hg The height of atmosG.
//! \param atmosG The first bounding atmosphere state.
//! \param hj The height of atmosJ.
//! \param atmosJ The second bounding atmosphere state.
//! \param h The fairing height (between hg and hj).
//! \param[out] atmosX The faired atmosphere state.
void EarthCommon::fair(int i, greal hg, const AtmosphereState& atmosG,
  greal hj, const AtmosphereState& atmosJ, greal h, AtmosphereState& atmosX)
{
  fair(i,
    hg, atmosG.pressure, atmosG.density, atmosG.temperature, atmosG.ewWind, atmosG.nsWind, atmosG.verticalWind,
    hj, atmosJ.pressure, atmosJ.density, atmosJ.temperature, atmosJ.ewWind, atmosJ.nsWind, atmosJ.verticalWind,
    h,  atmosX.pressure, atmosX.density, atmosX.temperature, atmosX.ewWind, atmosX.nsWind, atmosX.verticalWind);
}

//! \brief Fair pressure, density, temperature, and winds between two atmospheres.
//!
//! Fairs at height h, between pg,dg,tg,ug,vg,wg and pj,dj,tj,uj,vj,wj.  Fairing is cosine**2,
//! such that all "g" values obtain at height hg and all "j" values obtain at height hj.  The
//! faired results are p,d,t,u,v,w. 
//! \param i 1 means density (d) is faired and pressure (p) is computed from a faired 
//!          value of the gas-law constant.  If i not equal to 1, both p and d are faired directly.
//! \param hg The height of the G parameters.
//! \param pg The first bounding pressure.
//! \param dg The first bounding density.
//! \param tg The first bounding temperature.
//! \param ug The first bounding E/W wind.
//! \param vg The first bounding N/S wind.
//! \param wg The first bounding vertical wind.
//! \param hj The height of the J parameters.
//! \param pj The second bounding pressure.
//! \param dj The second bounding density.
//! \param tj The second bounding temperature.
//! \param uj The second bounding E/W wind.
//! \param vj The second bounding N/S wind.
//! \param wj The second bounding vertical wind.
//! \param h The fairing height (between hg and hj).
//! \param[out] p The faired pressure.
//! \param[out] d The faired density.
//! \param[out] t The faired temperature.
//! \param[out] u The faired E/W wind.
//! \param[out] v The faired N/S wind.
//! \param[out] w The faired vertical wind.
void EarthCommon::fair(int i, greal hg, greal pg, greal dg, greal tg, greal ug, greal vg, greal wg,
  greal hj, greal pj, greal dj, greal tj, greal uj, greal vj, greal wj, greal h,
  greal &p, greal &d, greal &t, greal &u, greal &v, greal &w)
{
  Interpolator cosInterp;
  cosInterp.makeCosineSquaredFraction(hg, hj, h);

  //Faired temperature
  t = cosInterp.linear(tg, tj);

  if (i == 1) {
    //Faired density
    d = cosInterp.log(dg, dj);

    //Faired gas constant and pressure
    greal rg = pg / (dg * tg);
    greal rj = pj / (dj * tj);
    greal r = cosInterp.linear(rg, rj);
    p = r * d * t;
  }
  else {
    d = cosInterp.linear(dg, dj);
    p = cosInterp.linear(pg, pj);
  }

  //Faired wind components
  u = cosInterp.linear(ug, uj);
  v = cosInterp.linear(vg, vj);
  w = cosInterp.linear(wg, wj);
}

} // namespace
