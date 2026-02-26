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
#include "MAP.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
MAP::MAP()
  : Atmosphere(), EarthCommon(this)
{
}

//! \fn MAP::MAP(const MAP& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn MAP::~MAP()
//! \copydoc Atmosphere::~Atmosphere()

//! \copydoc PerturbedAtmosphere::setInputParameters()
void MAP::setInputParameters(const EarthInputParameters& params)
{
  ewWindPerturbationScale = params.horizontalWindPerturbationScale;
  nsWindPerturbationScale = params.horizontalWindPerturbationScale;
  densityPerturbationScale = params.randomPerturbationScale;
  pressurePerturbationScale = params.randomPerturbationScale;
  temperaturePerturbationScale = params.randomPerturbationScale;
  month = params.month;
}

//! \fn MAP::getTemperatureGradient()
//! \brief Gets the temperature gradient with respect to height.
//!
//! \returns Temperature gradient with respect to height \units{K}.

//! \brief Compute mean pressure, density, temperature, wind components from the MAP data
//!
//! This routine will read the MAP data, if necessary, and perform interpolations based on the
//! current position.
//!
//! \b Inputs
//! \arg #position          
//! \arg MAP data          
//!
//! \retval #pressure
//! \retval #temperature
//! \retval #density
//! \retval #ewWind
//! \retval #nsWind
//! \retval #verticalWind
//! \retval #temperatureGradient
void MAP::update()
{
  initializeData();
  // Distances for 5 degrees of lat, lon
  double dy5 = toRadians(5000.0 * totalRadius);
  double dx5 = dy5 * cos(toRadians(latitude));
  dx5 = max(2000.0, dx5);

  // The following section is for zonal mean or mixed zonal mean
  // Upper height level
  int ihgb = 5 * (int(height) / 5) + 5;
  // Upper height
  double hgb = double(ihgb);

  // Zonal mean at upper height
  double pgb, dgb, tgb, ugb, dpygb, dtygb;
  gterp(ihgb, latitude, pgb, dgb, tgb, ugb, dpygb, dtygb);

  // Upper stationary perturbation height
  int ihsb = min(90, ihgb);
  double hsb = double(ihsb);

  // Stationary perturbations at upper height
  double psb, dsb, tsb, usb, vsb, dtxsb, dtysb;
  pdtuv(ihsb, psb, dsb, tsb, usb, vsb, dtxsb, dtysb);

  // Lower height level
  int ihga = ihgb - 5;
  double hga = double(ihga);
  // Zonal mean at lower height
  double pga, dga, tga, uga, dpyga, dtyga;
  gterp(ihga, latitude, pga, dga, tga, uga, dpyga, dtyga);

  // Lower stationary perturbation height
  int ihsa = ihsb - 5;
  double hsa = double(ihsa);
  double hs = min(90.0, height);

  // Stationary perturbations at lower height
  double psa, dsa, usa, vsa, tsa, dtxsa, dtysa;
  pdtuv(ihsa, psa, dsa, tsa, usa, vsa, dtxsa, dtysa);

  // Interpolate for winds
  Interpolator interpHeight;
  interpHeight.makeFraction(hga, hgb, height);
  double umb = interpHeight.linear(uga, ugb);
  double dtyg = interpHeight.linear(dtyga, dtygb);

  interpHeight.makeFraction(hsa, hsb, hs);
  double spu = interpHeight.linear(usa, usb);
  double spv = interpHeight.linear(vsa, vsb);
  double ush = umb + spu;
  double vsh = spv;

  // Stationary perturbations height interpolation
  double psh = interpHeight.linear(psa, psb);
  double dsh = interpHeight.linear(dsa, dsb);
  double tsh = interpHeight.linear(tsa, tsb);

  // Zonal mean values height interpolation
  double pmb = 0, dmb = 0, tmb = 0;
  heightInterpolation(pga, dga, tga, ihga, pgb, dgb, tgb, ihgb, pmb, dmb, tmb, height);

  // Height interpolation of stationary perturbation, dt/dx and dt/dy
  interpHeight.makeFraction(hsa, hsb, height);
  double dtxs = interpHeight.linear(dtxsa, dtxsb);
  double dtys = interpHeight.linear(dtysa, dtysb);

  // Unperturbed (monthly mean) values for output
  temperature = tmb * (1.0 + tsh);
  pressure = pmb * (1.0 + psh);
  density = dmb * (1.0 + dsh);
  ewWind = ush;
  nsWind = vsh;

  // Total dt/dx
  double dtx = dtxs * temperature;
  // Total dt/dy
  double dty = temperature * dtys + dtyg * (1.0 + tsh + dtys);

  double cp = 7.0 * pressure / (2.0 * density * temperature);
  temperatureGradient = (tgb * (1.0 + tsb) - tga * (1.0 + tsa)) / 5000.0;

  // Vertical mean wind
  verticalWind = -cp * (ewWind * dtx / dx5 + nsWind * dty / dy5) / (gravity + cp * temperatureGradient);
}

//! \brief Computes random perturbation standard deviations for PDT.
//!
//! Computes random perturbation standard deviations pressure, density, temperture 
//! at height h (km), geocentric latitude phi (degrees) from sigma arrays pr, dr, tr.
//!
//! \b Inputs
//! \arg #height          
//! \arg #latitude          
//! \arg #pr          
//! \arg #dr          
//! \arg #tr          
//! 
//! \param[out] pressureSD    A greal.
//! \param[out] densitySD     A greal.
//! \param[out] temperatureSD A greal.
//! 
//! \retval pressureSD     Pressure standard deviation \units{\%}.
//! \retval densitySD      Density standard deviation \units{\%}.
//! \retval temperatureSD  Temperature standard deviation \units{\%}.
void MAP::getStandardDeviations(greal& pressureSD, greal& densitySD, greal& temperatureSD)
{
  initializeData();
  // i = lower height index
  int i;
  if (height < 120.0) {
    i = int(0.2 * (height + 5.0)) - 1;
  }
  else {
    i = 24 + int(0.05 * (height - 120.0));
  }
  i = max(0, i);
  i = min(28, i);

  // ip = upper height index
  int ip = i + 1;
  ip = min(28, ip);

  // Lower latitude index
  int j = int((latitude + 100.0) / 10.0) - 1;
  int jp = j + 1;
  jp = min(18, jp);

  // Lower height for pr,tr,dr arrays
  greal z1;
  if (i > 24) {
    z1 = 120.0 + 20.0 * (i - 24);
  }
  else {
    z1 = 5.0 * (i + 1) - 5.0;
  }

  // Upper height for pr,tr,dr arrays
  greal z2;
  if (ip > 24) {
    z2 = 120.0 + 20.0 * (ip - 24);
  }
  else {
    z2 = 5.0 * (ip + 1) - 5.0;
  }

  greal phi1 = (-100.0) + 10.0 * (j + 1);
  greal phi2 = (-100.0) + 10.0 * (jp + 1);

  greal heightFraction;
  if (i == ip) {
    heightFraction = 1.0;
  }
  else {
    heightFraction = (height - z1) / (z2 - z1);
  }

  greal latitudeFraction;
  if (j == jp) {
    latitudeFraction = 1.0;
  }
  else {
    latitudeFraction = (latitude - phi1) / (phi2 - phi1);
  }

  //Interpolate on latitude and height
  Interpolator interp2d(heightFraction, latitudeFraction);

  int m = month - 1;
  pressureSD = square(pressurePerturbationScale)
      * interp2d.linear(pr[m][i][j],  pr[m][i][jp], 
                        pr[m][ip][j], pr[m][ip][jp]);

  densitySD = square(densityPerturbationScale)
      * interp2d.linear(dr[m][i][j],  dr[m][i][jp], 
                        dr[m][ip][j], dr[m][ip][jp]);

  temperatureSD = square(temperaturePerturbationScale)
      * interp2d.linear(tr[m][i][j],  tr[m][i][jp],
                        tr[m][ip][j], tr[m][ip][jp]);
}

//! \brief Computes random perturbation standard deviations for winds.
//!
//! Computes random perturbation standard deviations for winds 
//! at height h (km), geocentric latitude phi (degrees) from sigma arrays ur, vr.
//!
//! \b Inputs
//! \arg #height          
//! \arg #latitude          
//! \arg #ur          
//! \arg #vr          
//! 
//! \param[out] ewSD    A greal.
//! \param[out] nsSD    A greal.
//! 
//! \retval ewSD      East/west wind standard deviation \units{m/s}.
//! \retval nsSD      North/south wind standard deviation \units{m/s}.
void MAP::getWindStandardDeviations(greal& ewSD, greal& nsSD)
{
  assert(height >= 0.0);

  // Lower height index
  int i = int(height) / 5;
  if (height >= 125.0) {
    i = 24 + (int(height) - 120) / 20;
  }
  i = min(28, i);

  // Upper height index
  int ip = i + 1;
  ip = min(28, ip);

  assert(latitude <= 90.0);
  assert(latitude >= -90.0);

  // Lower latitude index
  int j = int(latitude + 100.0) / 10 - 1;

  // Upper latitude index
  int jp = j + 1;
  jp = min(18, jp);

  // phi1 - lower latitude for ur and vr array values
  greal phi1 = (-100.0) + 10.0 * (j + 1);

  // phi2 - upper latitude for ur and vr array values
  greal phi2 = (-100.0) + 10.0 * (jp + 1);

  // Lower height for ur and vr arrays
  greal z1;
  if (i > 24) {
    z1 = 20.0 * (i - 18);
  }
  else {
    z1 = 5.0 * i;
  }

  // Upper height for ur and vr array values
  greal z2;
  if (ip > 24) {
    z2 = 20.0 * (ip - 18);
  }
  else {
    z2 = 5.0 * ip;
  }

  greal heightFraction;
  if (i == ip) {
    heightFraction = 1.0;
  }
  else {
    heightFraction = (height - z1) / (z2 - z1);
  }

  greal latitudeFraction;
  if (j == jp) {
    latitudeFraction = 1.0;
  }
  else {
    latitudeFraction = (latitude - phi1) / (phi2 - phi1);
  }

  // Interpolate on geocentric latitude and height
  int m = month - 1;
  Interpolator interp2d(heightFraction, latitudeFraction);
  ewSD = interp2d.linear(ur[m][i][j], ur[m][i][jp],
                        ur[m][ip][j], ur[m][ip][jp]);
  nsSD = interp2d.linear(vr[m][i][j], vr[m][i][jp],
                        vr[m][ip][j], vr[m][ip][jp]);
  ewSD = sqrt(abs(ewSD)) * ewWindPerturbationScale;
  nsSD = sqrt(abs(nsSD)) * nsWindPerturbationScale;
}

//! \brief The main driver for the species concentration data.
//!
//! Computes species concentrations.
//!
//! \b Inputs
//! \arg #height          
//! 
//! \param iyr             The current (4-digit) year.
//! \param pres            Pressure \units{Pa}.
//! \param ppmWaterLow     Concentration of water from the lower atmosphere model \units{ppm}.
//! \param[out] waterSD    A greal.
//! \param[out] oxygenND   A greal.
//! \param[out] ppm        An EarthPPM struct.
//! 
//! \retval waterSD      Water standard deviation \units{m/s}.
//! \retval oxygenND     Oxygen (atomic) number density.
//! \retval ppm          Species concentrations \units{ppm}.
void MAP::getSpeciesConcentrations(int iyr, greal pres, greal ppmWaterLow,
                                    greal& waterSD, greal& oxygenND, EarthPPM& ppm)
{
  // Fractional rate of change (per year) for constituents
  constexpr greal rco2 = 0.5 / 100.0;
  constexpr greal rch4 = 0.9 / 100.0;
  constexpr greal rn2o = 0.3 / 100.0;
  constexpr greal rco = 0.7 / 100.0;
  constexpr greal ro30 = 0.23 / 100.0;
  constexpr greal ro340 = -0.5 / 100.0;
  constexpr greal pg1 = 3.0e04;
  constexpr greal pg2 = 2.5e04;

  static const greal al0p3 = log(0.3);
  static const greal al0p4 = log(0.4);
  static const greal al10 = log(10.0);
  static const greal al12p5 = log(12.5);
  static const greal al1500 = log(1500.0);
  static const greal al2000 = log(2000.0);
  static const greal al1 = log(1.0);
  static const greal al1p5 = log(1.5);
  static const greal alpg1 = log(pg1);
  static const greal alpg2 = log(pg2);

  initializeData();

  // Evaluate concentrations from MAP Handbook vol 31
  greal o3m, h2om, n2om, ch4m, oxynd;
  mapconc(pres, o3m, h2om, waterSD, n2om, ch4m, oxynd);

  // Evaluate LaRC and AFGL concentrations unless z > 120 km
  if (height > hgta[aHeightSize - 1]) {
    return;
  }

  greal ppmWaterLARC = larcwat();
  greal ppmWaterAFGL = afglconc(ppm);

  // Update concentrations from 1976 to year of output:
  // CO2, CH4, CO, and N2O annual rates based on Table 14.5, page 306
  // of Graedel and Crutzen, "Atmospheric Change", 1993.
  int iypt0 = iyr - 1976;
  ppm.nitrousOxide = ppm.nitrousOxide * pow(1.0 + rn2o, iypt0);
  ppm.carbonMonoxide = ppm.carbonMonoxide * pow(1.0 + rco, iypt0);
  ppm.methane = ppm.methane * pow(1.0 + rch4, iypt0);
  ppm.carbonDioxide = ppm.carbonDioxide * pow(1.0 + rco2, iypt0);

  // Height-dependent annual rate for ozone change based on
  // Figure 17.1, page 373 of Graedel and Crutzen.
  greal ro3;
  if (height < 15.0) {
    ro3 = ro30 * ((15.0 - height) / 15.0);
  }
  else if (height < 30.0) {
    ro3 = 0.0;
  }
  else if (height < 40.0) {
    ro3 = ro340 * ((height - 30.0) / 10.0);
  }
  else {
    ro3 = ro340 * ((120.0 - height) / 80.0);
  }

  ppm.ozone = ppm.ozone * pow(1.0 + ro3, iypt0);

  // Rates of change for water vapor, oxygen and nitrogen are assumed
  // to be zero.  Update concetrations to year of output for MAP data.
  iypt0 = iyr - 1981;
  n2om = n2om * pow(1.0 + rn2o, iypt0);
  ch4m = ch4m * pow(1.0 + rch4, iypt0);
  o3m = o3m * pow(1.0 + ro3, iypt0);

  // Get log of p for interpolations
  greal alp = log(pres);

  // For water vapor use the following:
  // (a) AFGL value (ppmWaterAFGL) above the 0.01 mb height level
  // (b) Fair between AFGL and MAP vol 31 (h2om) for 0.015-0.01 bm
  // (c) MAP vol 31 value between 40.5 km and 0.015 mb
  // (d) Fair between MAP vol 31 and LaRC (ppmWaterLARC) for 39.5-40.5 km
  // (e) LaRC value between 39.5 km and 250 mb level
  // (f) Fair between LaRC and NCEP (ppmWaterNCEP) values 300-250 mb
  // (g) NCEP (or RH-extended) (ppmWaterNCEP) for surface to 300 mb
  if (pres <= 1.5) {
    if (pres <= 1.0) {
      ppm.water = ppmWaterAFGL;
    }
    else {
      ppm.water = ((alp - al1) * h2om + (al1p5 - alp) * ppmWaterAFGL) / (al1p5 - al1);
    }
  }
  else if (height >= 39.5) {
    if (height >= 40.5) {
      ppm.water = h2om;
    }
    else {
      ppm.water = (height - 39.5) * h2om + (40.5 - height) * ppmWaterLARC;
    }
  }
  else if (pres <= pg1) {
    if (pres <= pg2) {
      ppm.water = ppmWaterLARC;
    }
    else {
      ppm.water = ((alp - alpg2) * ppmWaterLow + (alpg1 - alp) * ppmWaterLARC) / (alpg1 - alpg2);
    }
  }
  else {
    ppm.water = ppmWaterLow;
  }

  // Use AFGL O3, N2O, and CH4 if pressure > 2000 Pa; fair 1500 - 2000 Pa
  // Use AFGL O3 if p < 0.3 Pa; fair 0.3 - 0.4 Pa
  // Use AFGL N2O and CH4 if p < 10 Pa; fair 10 - 12.5 Pa
  // Otherwise use MAP concentrations
  if (pres <= 2000.0) {
    if (pres >= 1500.0) {
      ppm.nitrousOxide = ((alp - al1500) * ppm.nitrousOxide + (al2000 - alp) * n2om) / (al2000 - al1500);
      ppm.methane = ((alp - al1500) * ppm.methane + (al2000 - alp) * ch4m) / (al2000 - al1500);
      ppm.ozone = ((alp - al1500) * ppm.ozone + (al2000 - alp) * o3m) / (al2000 - al1500);
    }
    else if (pres <= 12.5) {
      if (pres >= 10.0) {
        ppm.nitrousOxide = ((alp - al10) * n2om + (al12p5 - alp) * ppm.nitrousOxide) / (al12p5 - al10);
        ppm.methane = ((alp - al10) * ch4m + (al12p5 - alp) * ppm.methane) / (al12p5 - al10);
      }

      if (pres <= 0.4) {
        if (pres >= 0.3) {
          ppm.ozone = ((alp - al0p3) * o3m + (al0p4 - alp) * ppm.ozone) / (al0p4 - al0p3);
        }
      }
      else {
        ppm.ozone = o3m;
      }
    }
    else {
      ppm.nitrousOxide = n2om;
      ppm.methane = ch4m;
      ppm.ozone = o3m;
    }
  }

  // Use MET oxygen z > 100; MAP oxygen z < 90; fair otherwise
  if (height < 100.0) {
    if (height < 90.0) {
      oxygenND = oxynd;
    }
    else {
      oxygenND = ((height - 90.0) * (oxygenND) + (100.0 - height) * oxynd) / 10.0;
    }
  }
}


//! \brief Interpolates arrays of middle atmosphere values.
//!
//! Interpolates pg, dg, tg, ug arrays of middle atmosphere values to p,
//! d, t, u, at height ih (integer km) and geocentric latitude phi
//! (degrees).  Also computes N-S pressure and temperature gradients dp/dy, dt/dy.
//! 
//! \b Inputs
//! \arg #dg          
//! \arg #tg          
//! \arg #pg          
//! \arg #ug          
//! 
//! \param ih         A height level (multiple of 5).
//! \param phi        Latitude \units{\text{degrees}}.
//! \param[out] p     A greal.
//! \param[out] d     A greal.
//! \param[out] t     A greal.
//! \param[out] u     A greal.
//! \param[out] dpy   A greal.
//! \param[out] dty   A greal.
//! 
//! \retval p     Pressure \units{Pa}.
//! \retval d     Density \units{kg/m^3}.
//! \retval t     Temperature \units{K}.
//! \retval u     E/W winds \units{m/s}.
//! \retval dpy   dp/dy for geostrophic winds.
//! \retval dty   dt/dy for vertical winds.
void MAP::gterp(int ih, greal phi, greal& p, greal& d, greal& t, greal& u, greal& dpy, greal& dty)
{
  // Interpolates zonal mean data to height ih and geocentric latitude phi
  // Height index
  int i = (ih - 20) / 5;
  i = min(20, i);

  // Lower latitude index
  int j = int((phi + 100.0) / 10.0) - 1;
  j = max(0, j);
  j = min(17, j);
  // Upper latitude index
  int jp = j + 1;

  int m = month - 1;

  // Check for density or temperature leq 0
  greal chk = dg[m][i][j] * tg[m][i][j] * dg[m][i][jp] * tg[m][i][jp];
  if (chk <= 0.0) {
    p = pg[m][i][j];
    d = dg[m][i][j];
    t = tg[m][i][j];
  }

  // Geocentric latitude deviation from zonal mean position
  greal phif = (phi + 100.0 - 10.0 * (j + 1)) / 10.0;

  // Latitude interpolation
  Interpolator interpLat(phif);
  greal tl = interpLat.linear(tg[m][i][j], tg[m][i][jp]);
  greal dl = interpLat.linear(dg[m][i][j], dg[m][i][jp]);
  u = interpLat.linear(ug[m][i][j], ug[m][i][jp]);

  greal r1 = pg[m][i][j] / (dg[m][i][j] * tg[m][i][j]);
  greal r2 = pg[m][i][jp] / (dg[m][i][jp] * tg[m][i][jp]);
  // Interpolated gas constant
  greal r = interpLat.linear(r1, r2);

  // Pressure computed from interpolated gas constant
  p = dl * r * tl;
  d = dl;
  t = tl;
  // dt/dy for vertical wind
  dty = (tg[m][i][jp] - tg[m][i][j]) * 0.5;
  // dp/dy for geostrophic winds
  dpy = (pg[m][i][jp] - pg[m][i][j]) * 0.5;
}

//! \brief Interpolates stationary perturbations on geocentric latitude and longitude at height ih.
//! 
//! \b Inputs
//! \arg #longitude
//! \arg #latitude          
//! \arg #psp          
//! \arg #dsp          
//! \arg #tsp          
//! \arg #usp          
//! \arg #vsp          
//! 
//! \param ih          A height level (multiple of 5).
//! \param[out] ps     A greal.
//! \param[out] ds     A greal.
//! \param[out] ts     A greal.
//! \param[out] us     A greal.
//! \param[out] vs     A greal.
//! \param[out] dtx    A greal.
//! \param[out] dty    A greal.
//! 
//! \retval ps     Pressure \units{Pa}.
//! \retval ds     Density \units{kg/m^3}.
//! \retval ts     Temperature \units{K}.
//! \retval us     E/W winds \units{m/s}.
//! \retval vs     N/S winds \units{m/s}.
//! \retval dtx    dt/dx for vertical wind.
//! \retval dty    dt/dy for vertical wind.
void MAP::pdtuv(int ih, greal& ps, greal& ds, greal& ts, greal& us, greal& vs, greal& dtx, greal& dty)
{
  //Height index k
  int k = (ih - 20) / 5;

  greal xlon = longitude;
  // dlon - relative longitude deviation from corner reference point
  if (xlon >= 180.0) {
    xlon = xlon - 360.0;
  }
  if (xlon < (-180.0)) {
    xlon = xlon + 360.0;
  }

  // Lower longitude index j
  int j = int((xlon + 200.0) / 20.0) - 1;
  greal dlon = (xlon + 200.0 - 20.0 * (j + 1)) / 20.0;
  // Upper longitude index jp
  int jp = j + 1;
  if (jp > 17) {
    jp = 0;
  }

  // Lower latitude index i
  int i = int((latitude + 100.0) / 10.0) - 1;
  // Upper latitude index ip
  int ip = i + 1;
  ip = min(18, ip);
  // dlat - relative latitude deviation from corner reference location
  greal dlat = (latitude - 10.0 * (i + 1) + 100.0) / 10.0;

  Interpolator interp2D(dlat, dlon);
  int m = month - 1;

  // Pressure lat-lon interpolation
  ps = interp2D.linear(psp[m][k][i][j], psp[m][k][i][jp], psp[m][k][ip][j], psp[m][k][ip][jp]);
  // Density lat-lon interpolation
  ds = interp2D.linear(dsp[m][k][i][j], dsp[m][k][i][jp], dsp[m][k][ip][j], dsp[m][k][ip][jp]);
  // Temperature lat-lon interpolation
  ts = interp2D.linear(tsp[m][k][i][j], tsp[m][k][i][jp], tsp[m][k][ip][j], tsp[m][k][ip][jp]);
  // Zonal wind lat-lon interpolation
  us = interp2D.linear(usp[m][k][i][j], usp[m][k][i][jp], usp[m][k][ip][j], usp[m][k][ip][jp]);
  // Meridional wind lat-lon interpolation
  vs = interp2D.linear(vsp[m][k][i][j], vsp[m][k][i][jp], vsp[m][k][ip][j], vsp[m][k][ip][jp]);
  // dtx - dt/dx for vertical wind
  greal dtxa = (tsp[m][k][i][jp] - tsp[m][k][i][j]) / 4.0;
  dtx = dtxa + ((tsp[m][k][ip][jp] - tsp[m][k][ip][j]) / 4.0 - dtxa) * dlat;
  // dty - dt/dy for vertical wind
  greal dtya = (tsp[m][k][ip][j] - tsp[m][k][i][j]) / 2.0;
  dty = dtya + ((tsp[m][k][ip][jp] - tsp[m][k][i][jp]) / 2.0 - dtya) * dlon;
}

//! \brief Performs height interpolation for PDT values.
//!
//! Interpolates between p1, d1, t1 at height z1 and p2, d2, t2 at 
//! height z2 to output values of p, d, t at height z.
//! Checks for t1, d1, t2, d2 product = 0, for gas interpolation.
//! 
//! \param p1          Pressure for lower level \units{Pa}.
//! \param d1          Density for lower level \units{kg/m^3}.
//! \param t1          Temperature for lower level \units{K}.
//! \param z1          Height of lower level \units{km}.
//! \param p2          Pressure for upper level \units{Pa}.
//! \param d2          Density for upper level \units{kg/m^3}.
//! \param t2          Temperature for upper level \units{K}.
//! \param z2          Height of upper level \units{km}.
//! \param[out] p      A greal.
//! \param[out] d      A greal.
//! \param[out] t      A greal.
//! \param z           Interpolation height \units{km}.
//! 
//! \retval p     Pressure \units{Pa}.
//! \retval d     Density \units{kg/m^3}.
//! \retval t     Temperature \units{K}.
void MAP::heightInterpolation(greal p1, greal d1, greal t1, greal z1,
                 greal p2, greal d2, greal t2, greal z2, 
                 greal& p, greal& d, greal& t, greal z)
{
	// Set p=d=t=0 if some input values are negative or zero
  if (p1 * d1 * t1 * p2 * d2 * t2 <= 0.0 || p2 * p1 <= 0.0) {
		p = 0.0;
		d = 0.0;
		t = 0.0;
		return;
	}

	// Sets p,d,t = p1,d1,t1 if z1=z2
	if (abs(z1 - z2) <= 0.001){
		p = p1;
		d = d1;
		t = t1;
		return;
	}

  // Height based interpolation
  Interpolator interpHeight;
  interpHeight.makeFraction(z1, z2, z);

	// Linear interpolation on t
	t = interpHeight.linear(t1, t2);

  // Gas constants from perfect gas law
  greal r1 = p1 / (d1 * t1);
  greal r2 = p2 / (d2 * t2);
  // Linear interpolation on gas constant r
  greal r = interpHeight.linear(r1, r2);

	// For pressure, use logarithmic interpolation if temperature gradient = 0
	if (abs(t2 - t1) <= 0.001){
    p = interpHeight.log(p1, p2);
	}
	else {
		// Otherwise use a power law interpolation on pressure
    greal exponent = log(p2 / p1) / log(t1 / t2);
    greal base = (t1 / t);
		p = p1 * pow(base, exponent);
	}

  // Density from perfect gas law
  d = p / (r * t);
}

//! \brief Evaluate the AFGL concentration data for a given location.
//! 
//! \b Inputs
//! \arg #height
//! \arg #latitude          
//! \arg #hgta          
//! \arg #h2oa          
//! \arg #o3          
//! \arg #n2o          
//! \arg #co          
//! \arg #ch4          
//! \arg #co2          
//! \arg #o2          
//! \arg #n2          
//! 
//! \param[out] ppm    Species concentrations \units{ppm}.
//! 
//! \returns  Water concentration \units{ppm}.
double MAP::afglconc(EarthPPM& ppm)
{
	const greal xlat[8] = { -90.0, -60.0, -45.0, -15.0, 15.0, 45.0, 60.0, 90.0 };

	if (height < hgta[0] || height > hgta[aHeightSize - 1]){
    greal ppmWater = 0.0;
		ppm.ozone = 0.0;
		ppm.nitrousOxide = 0.0;
		ppm.carbonMonoxide = 0.0;
		ppm.methane = 0.0;
		ppm.carbonDioxide = 0.0;
		ppm.dioxygen = 0.0;
		ppm.dinitrogen = 0.0;
		return ppmWater;
	}

  size_t upperHeight = 1;
  greal heightFraction = 0.0;
  // Find height index and height interpolation coefficient
  while (upperHeight < aHeightSize - 1 && height > hgta[upperHeight]) {
    ++upperHeight;
  }
  size_t lowerHeight = upperHeight - 1;
  heightFraction = (height - hgta[lowerHeight]) / (hgta[upperHeight] - hgta[lowerHeight]);

  size_t upperLat = 1;
  greal latFraction = 0.0;
  // Find latitude index and latitude interpolation coefficient
  while (upperLat < aLatSize && latitude > xlat[upperLat]) {
    ++upperLat;
  }
  size_t lowerLat = upperLat - 1;
  latFraction = (latitude - xlat[lowerLat]) / (xlat[upperLat] - xlat[lowerLat]);

  int m = month - 1;
  // Latitude-height interpolation for H2O, O3, N2O, CO, and CH4
  Interpolator interp2D(latFraction, heightFraction);
  greal ppmWater = interp2D.log(h2oa[m][lowerLat][lowerHeight], h2oa[m][lowerLat][upperHeight],
                                h2oa[m][upperLat][lowerHeight], h2oa[m][upperLat][upperHeight]);
  ppm.ozone = interp2D.log(o3[m][lowerLat][lowerHeight], o3[m][lowerLat][upperHeight],
                           o3[m][upperLat][lowerHeight], o3[m][upperLat][upperHeight]);
  ppm.nitrousOxide = interp2D.log(n2o[m][lowerLat][lowerHeight], n2o[m][lowerLat][upperHeight],
                                  n2o[m][upperLat][lowerHeight], n2o[m][upperLat][upperHeight]);
  ppm.carbonMonoxide = interp2D.log(co[m][lowerLat][lowerHeight], co[m][lowerLat][upperHeight],
                                    co[m][upperLat][lowerHeight], co[m][upperLat][upperHeight]);
  ppm.methane = interp2D.log(ch4[m][lowerLat][lowerHeight], ch4[m][lowerLat][upperHeight],
                             ch4[m][upperLat][lowerHeight], ch4[m][upperLat][upperHeight]);

  // Height interpolation only for CO2, O2, N2
  Interpolator interpHeight(heightFraction);
  ppm.carbonDioxide = interpHeight.log(co2[lowerHeight], co2[upperHeight]);
  ppm.dioxygen = interpHeight.log(o2[lowerHeight], o2[upperHeight]);
  ppm.dinitrogen = interpHeight.log(n2[lowerHeight], n2[upperHeight]);

  return ppmWater;
}

//! \brief Evaluates the MAP concentration data at a given point
//! 
//! \b Inputs
//! \arg #month
//! \arg #latitude          
//! \arg #ph2omap          
//! \arg #h2omap          
//! \arg #sh2omap          
//! \arg #po3map          
//! \arg #o3map          
//! \arg #pmap31          
//! \arg #n2omap          
//! \arg #ch4map          
//! \arg #mapzox          
//! \arg #oxymap          
//! 
//! \param pres        Pressure \units{Pa}.
//! \param[out] o3m    A greal.
//! \param[out] h2om   A greal.
//! \param[out] sh2om  A greal.
//! \param[out] n2om   A greal.
//! \param[out] ch4m   A greal.
//! \param[out] oxynd  A greal.
//! 
//! \retval o3m    Ozone concentration \units{ppm}.
//! \retval h2om   Water vapor concentration \units{ppm}.
//! \retval sh2om  Water vapor concentration standard deviation \units{ppm}.
//! \retval n2om   Nitrous concentration \units{ppm}.
//! \retval ch4m   Methane concentration \units{ppm}.
//! \retval oxynd  Atomic oxygen number density.
void MAP::mapconc(greal pres, greal& o3m, greal& h2om, greal& sh2om, greal& n2om, greal& ch4m,
                  greal& oxynd)
{
  static const greal xlat[8] = {-90.0, -60.0, -45.0, -15.0, 15.0, 45.0, 60.0, 90.0};

  o3m = 0.0;
  h2om = 0.0;
  sh2om = 0.0;
  n2om = 0.0;
  ch4m = 0.0;
  oxynd = 0.0;

  // Convert pressure to log scale
  greal logPressure = log(pres);

  int m = month - 1;

  // Find H2O pressure level index and interpolation coefficient
  size_t upperHeight = getUpperIndex(logPressure, 19, ph2omap);
  size_t lowerHeight = upperHeight - 1;
  greal h2oHeightFraction = (logPressure - ph2omap[lowerHeight]) / (ph2omap[upperHeight] - ph2omap[lowerHeight]);

  // Find latitude index and latitude interpolation coefficient Index and coefficient for H2O
  size_t upperLat = getUpperIndex(latitude, 8, xlat);
  size_t lowerLat = upperLat - 1;
  greal h2oLatFraction = (latitude - xlat[lowerLat]) / (xlat[upperLat] - xlat[lowerLat]);

  // Find MAP water vapor value and standard deviation
  if (upperHeight > 0) {
    Interpolator interp2D(h2oLatFraction, h2oHeightFraction);
    h2om = interp2D.log(h2omap[m][lowerLat][lowerHeight], h2omap[m][lowerLat][upperHeight],
                        h2omap[m][upperLat][lowerHeight], h2omap[m][upperLat][upperHeight]);
    sh2om = interp2D.log(sh2omap[m][lowerLat][lowerHeight], sh2omap[m][lowerLat][upperHeight],
                         sh2omap[m][upperLat][lowerHeight], sh2omap[m][upperLat][upperHeight]);
  }

  // Index and coefficient for all other MAP species
  upperLat = int((latitude + 90.0) / 10.0) + 1;
  lowerLat = upperLat - 1;

  if (upperLat > 18) {
    upperLat = 18;
  }
  greal latFraction = (latitude + 100.0 - 10.0 * upperLat) / 10.0;

  // Find ozone pressure level index and interpolation coefficient
  upperHeight = getUpperIndex(logPressure, 24, po3map);
  lowerHeight = upperHeight > 0 ? upperHeight - 1 : 0;
  greal ozoneHeightFraction = (logPressure - po3map[lowerHeight]) / (po3map[upperHeight] - po3map[lowerHeight]);

  // Find MAP ozone value
  if (upperHeight > 0) {
    Interpolator interp2D(latFraction, ozoneHeightFraction);
    o3m = interp2D.log(o3map[m][lowerLat][lowerHeight], o3map[m][lowerLat][upperHeight],
                       o3map[m][upperLat][lowerHeight], o3map[m][upperLat][upperHeight]);
  }

  // Find N2O and CH4 level index and interpolation coefficient
  upperHeight = getUpperIndex(logPressure, 17, pmap31);
  lowerHeight = upperHeight > 0 ? upperHeight - 1 : 0;
  greal beta2 = (logPressure - pmap31[lowerHeight]) / (pmap31[upperHeight] - pmap31[lowerHeight]);

  // Find MAP N2O and CH4 values (convert N2O to ppm)
  if (upperHeight > 0) {
    Interpolator interp2D(latFraction, beta2);
    n2om = interp2D.log(n2omap[m][lowerLat][lowerHeight], n2omap[m][lowerLat][upperHeight],
                        n2omap[m][upperLat][lowerHeight], n2omap[m][upperLat][upperHeight]) / 1000.0;
    ch4m = interp2D.log(ch4map[m][lowerLat][lowerHeight], ch4map[m][lowerLat][upperHeight],
                        ch4map[m][upperLat][lowerHeight], ch4map[m][upperLat][upperHeight]);
  }

  // Find atmoic oxygen height level and interpolation coefficient
  upperHeight = getUpperIndex(height, 19, mapzox);
  lowerHeight = upperHeight > 0 ? upperHeight - 1 : 0;
  greal oxyHeightFraction = (height - mapzox[lowerHeight]) / (mapzox[upperHeight] - mapzox[lowerHeight]);

  // Find MAP atomic oxygen values (convert to #/m^3)
  if (upperHeight > 0) {
    Interpolator interp2D(latFraction, oxyHeightFraction);
    oxynd = interp2D.log(oxymap[m][lowerLat][lowerHeight], oxymap[m][lowerLat][upperHeight],
                         oxymap[m][upperLat][lowerHeight], oxymap[m][upperLat][upperHeight]) * 1.0e06;
  }
}

//! \brief Water vapor concentration from NASA Langley climatology.
//!
//! Water vapor concentration from NASA Langley climatology, as a function
//! of height (km) and latitude (degrees)
//! 
//! \b Inputs
//! \arg #month
//! \arg #height
//! \arg #latitude          
//! \arg #hgtl          
//! \arg #h2ol          
//! 
//! \returns Water vapor concentration \units{ppm}.
greal MAP::larcwat()
{
  greal ppmWater = 0.0;
  if (height >= hgtl[0] && height <= hgtl[34]) {
    // Height index
    int upperHeight = 1 + int(height - hgtl[0]);
    if (upperHeight >= 35) {
      upperHeight = 34;
    }
    int lowerHeight = upperHeight - 1;

    // Latitude index
    greal xlat = (latitude + 90.0) / 20.0;
    int upperLat = 1 + int(xlat);
    if (upperLat >= 10) {
      upperLat = 9;
    }
    int lowerLat = upperLat - 1;

    // Interpolation fractions.
    greal latFraction = xlat + 1.0 - upperLat;
    greal heightFraction = height - hgtl[lowerHeight];

    int m = month - 1;
    Interpolator interp(latFraction, heightFraction);
    ppmWater = interp.log(h2ol[m][lowerLat][lowerHeight], h2ol[m][lowerLat][upperHeight],
                          h2ol[m][upperLat][lowerHeight], h2ol[m][upperLat][upperHeight]);
  }

  return ppmWater;
}

} // namespace
