//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United
// States without explicit approval by NASA Marshall Space Flight Center.
//
// Module: Earth-GRAM
// Point of Contact: P. White
//////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cmath>
#include "MSIS.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
MSIS::MSIS()
  : Atmosphere(), EarthCommon(this)
{
}

//! \fn MSIS::MSIS(const MSIS& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn MSIS::~MSIS()
//! \copydoc Atmosphere::~Atmosphere()

//! \copydoc PerturbedAtmosphere::setInputParameters()
void MSIS::setInputParameters(const EarthInputParameters &params)
{
  year = params.year;
  f10 = params.dailyF10;
  f10b = params.meanF10;
  ap = params.ap;

  month = params.month;
  day = params.day;
  hour = params.hour;
  minute = params.minute;
  seconds = params.seconds;
}

//! \fn MSIS::setDayOfYear()
//! \brief Sets the day of the year.
//!
//! \param doy   Day of the year (0-367).

//! \fn MSIS::getTemperatureGradient()
//! \brief Gets the temperature gradient with respect to height.
//!
//! \returns Temperature gradient with respect to height \units{K}.

//! \fn MSIS::getMolecularWeightGradient()
//! \brief Gets the molecular weight gradient with respect to height.
//!
//! \returns Molecular weight gradient with respect to height.

//! \fn MSIS::getCswSw(int index)
//! \brief Get the sw value at spedified index.

//! \fn MSIS::getCswSwc(int index)
//! \brief Get the swc value at spedified index.

//! \fn MSIS::zeta(greal zz, greal zl, greal re)
//! \brief Unknown.
//!
//! \param zz
//! \param zl
//! \param re
//!
//! \returns Unknown.


//! \brief MSIS driver routine to evaluate mean p, d, t, u, v, w
//!
//! \b Inputs
//! \arg #position          
//!
//! \returns  pressure, density, temperature, ewWind, nsWind, verticalWind, 
//! gas number densities, averageMolecularWeight, temperature gradient, molecular weight gradient
void MSIS::update()
{
  int imin = int(minute + int(elapsedTime) / 60);
  greal sec = seconds + elapsedTime;
  sec = int(sec) % 60;
  imin = imin % 60;

  // Distances for 5 degrees of latitude, longitude
  greal dy5 = toRadians(5000.0 * totalRadius);
  greal dx5 = dy5 * cos(toRadians(latitude));

  if (dx5 < 2000.0) {
    dx5 = 2000.0;
  }

  // Following is the pure thermosphere height range section
  // MSIS and HWM values at current position
  greal apArray[7] = {ap, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  greal tex;
  msishwm(year, height, latitude, longitude, f10b, f10, apArray, pressure, density,
          temperature, dinitrogen.numberDensity, dioxygen.numberDensity, oxygen.numberDensity,
          argon.numberDensity, helium.numberDensity, hydrogen.numberDensity, nitrogen.numberDensity,
          averageMolecularWeight, tex, ewWind, nsWind);

  // Set latitude increment for temperature gradients
  greal dphi = 5.0;

  if (latitude >= 85.0) {
    dphi = -5.0;
  }

  greal phn, dhn, thn, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11;

  // MSIS temperature at current position latitude increment
  msishwm(year, height, latitude + dphi, longitude, f10b, f10, apArray, phn, dhn, thn,
          d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11);

  // MSIS temperature at current position+5 degrees lon
  greal phe, dhe, the;
  msishwm(year, height, latitude, longitude + 5.0, f10b, f10, apArray, phe, dhe, the,
          d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11);

  // dt/dx, dt/dy, and dt/dz for vertical wind
  greal dtx = the - temperature;
  greal dty = thn - temperature;

  if (dphi < 0.0) {
    dty = -dty;
  }

  // MSIS temperature and molecular weight 1 km higher
  greal pb, db, tb;
  msishwm(year, height + 1.0, latitude, longitude, f10b, f10, apArray, pb, db, tb, d1,
          d2, d3, d4, d5, d6, d7, d8, d9, d10, d11);

  // Gradients for temperature and molecular weight
  dtz = (tb - temperature) / 1000.0;
  dmdz = (d8 - averageMolecularWeight) / 1000.0;

  // Compute vertical mean wind Specific heat
  greal cp = 7.0 * (pressure) / (2.0 * (density) * (temperature));

  // Mean vertical wind from Montgomery stream function
  verticalWind = -cp * (ewWind * dtx / dx5 + nsWind * dty / dy5) / (gravity + cp * dtz);
}

//! \brief Calculate temperature and density profiles for MSIS models
//!
//! New lower thermo polynomial 10/30/89
//!
//! \param alt
//! \param dlb
//! \param tinf
//! \param tlb
//! \param xm
//! \param alpha
//! \param[out] tz
//! \param zlb
//! \param s2
//! \param mn1
//! \param zn1
//! \param[out] tn1
//! \param[in,out] tgn1
//!
//! \retval tz Temperature at altitude.
//! \retval tn1 Temperature at nodes.
//! \retval tgn1 Temperature gradient at end nodes.
//! \returns Density at altitude
greal MSIS::densu(greal alt, greal dlb, greal tinf, greal tlb, greal xm, greal alpha,
                   greal &tz, greal zlb, greal s2, int mn1, const greal zn1[5], greal tn1[5],
                   greal tgn1[2])
{
  constexpr greal rgas = UNIVERSAL_GAS * 100.0;

  // Joining altitudes of Bates and spline
  greal za = zn1[0];
  greal z = fmax(alt, za);

  // Geopotantial altitude difference from zlb
  greal zg2 = zeta(z, zlb, re);

  // Bates temperature
  greal tt = tinf - (tinf - tlb) * exp(-s2 * zg2);
  greal ta = tt;
  tz = tt;
  greal densu_out = tz;

  greal z1, t1, zgdiff, x;
  greal xs[5], ys[5], y2out[5];
  int mn;
  if (alt < za) {
    // Calculate temperature below za
    // Temperature gradient at za from Bates profile
    greal dta = (tinf - ta) * s2 * pow(((re + zlb) / (re + za)), 2.0);
    tgn1[0] = dta;
    tn1[0] = ta;
    z = fmax(alt, zn1[mn1 - 1]);
    mn = mn1;
    z1 = zn1[0];
    greal z2 = zn1[mn - 1];
    t1 = tn1[0];
    greal t2 = tn1[mn - 1];
    // Geopotential difference from z1
    greal zg = zeta(z, z1, re);
    zgdiff = zeta(z2, z1, re);

    // Set up spline nodes
    for (int k = 0; k < mn; k++) {
      xs[k] = zeta(zn1[k], z1, re) / zgdiff;
      ys[k] = 1.0 / tn1[k];
    }

    // End node derivatives
    greal yd1 = -tgn1[0] / (t1 * t1) * zgdiff;
    greal yd2 = -tgn1[1] / (t2 * t2) * zgdiff * pow((re + z2) / (re + z1), 2.0);

    // Calculate spline coefficients
    hwm.spline(xs, ys, mn, yd1, yd2, y2out);
    x = zg / zgdiff;

    greal y;
    hwm.splint(xs, ys, y2out, mn, x, y);

    // Temperature at altitude
    tz = 1.0 / y;

    densu_out = tz;
  }

  if (xm != 0.0) {
    // Calculate density above za
    greal glb = gsurf / pow(1.0 + zlb / re, 2.0);

    greal gamma1 = xm * glb / (s2 * rgas * tinf);
    greal expll = exp(-s2 * gamma1 * zg2);

    if ((expll > 50.0) || (0.0 >= tt)) {
      expll = 50.0;
    }

    // Density at altitude
    densu_out = dlb * pow(tlb / tt, 1.0 + alpha + gamma1) * expll;

    if (alt < za) {
      // Calculate density below za
      glb = gsurf / pow(1.0 + z1 / re, 2.0);

      greal gamm = xm * glb * zgdiff / rgas;

      // Integrate spline temperatures
      greal yi;
      splini(xs, ys, y2out, mn, x, yi);
      expll = gamm * yi;

      if ((expll > 50.0) || (tz <= 0.0)) {
        expll = 50.0;
      }

      densu_out = densu_out * pow(t1 / (tz), 1.0 + alpha) * exp(-expll);
    }
  }

  return densu_out;
}

//! \brief Calculate temperature and density profiles for lower atmos
//!
//! \param alt
//! \param d0
//! \param xm
//! \param[out] tz
//! \param mn3
//! \param zn3
//! \param tn3
//! \param tgn3
//! \param mn2
//! \param zn2
//! \param tn2
//! \param tgn2
//!
//! \retval tz Temperature at altitude.
//! \returns Density at altitude
greal MSIS::densm(greal alt, greal d0, greal xm, greal &tz, int mn3, greal zn3[], greal tn3[],
                  greal tgn3[2], int mn2, greal zn2[], greal tn2[], greal tgn2[2])
{
  constexpr greal rgas = UNIVERSAL_GAS * 100.0;

  greal densm_out = d0;

  if (alt <= zn2[0]) {
    // Stratosphere/mesosphere temperature
    greal z = fmax(alt, zn2[mn2 - 1]);
    int mn = mn2;
    greal z1 = zn2[0];
    greal z2 = zn2[mn - 1];
    greal t1 = tn2[0];
    greal t2 = tn2[mn - 1];
    greal zg = zeta(z, z1, re);
    greal zgdif = zeta(z2, z1, re);

    // Set up slpine nodes
    greal xs[10], ys[10];
    for (int k = 0; k < mn; k++) {
      xs[k] = zeta(zn2[k], z1, re) / zgdif;
      ys[k] = 1.0 / tn2[k];
    }

    greal yd1 = -tgn2[0] / (t1 * t1) * zgdif;
    greal yd2 = -tgn2[1] / (t2 * t2) * zgdif * pow((re + z2) / (re + z1), 2);

    // Calculate spline coefficients
    greal y2out[10];
    hwm.spline(xs, ys, mn, yd1, yd2, y2out);
    greal x = zg / zgdif;
    greal y;
    hwm.splint(xs, ys, y2out, mn, x, y);

    // Temperature at altitude
    tz = 1.0 / y;

    greal glb, gamm, yi, expll;
    if (xm != 0.0) {
      // Calculate stratosphere/mesosphere density
      glb = gsurf / pow(1.0 + z1 / re, 2);
      gamm = xm * glb * zgdif / rgas;

      // Integrate temperature profile
      splini(xs, ys, y2out, mn, x, yi);
      expll = gamm * yi;

      if (expll > 50.0) {
        expll = 50.0;
      }

      // Density at altitude
      densm_out = densm_out * (t1 / (tz)) * exp(-expll);
    }

    if (alt <= zn3[0]) {
      // Troposphere/stratosphere temperature
      z = alt;
      mn = mn3;
      z1 = zn3[0];
      z2 = zn3[mn - 1];
      t1 = tn3[0];
      t2 = tn3[mn - 1];
      zg = zeta(z, z1, re);
      zgdif = zeta(z2, z1, re);

      // Set up spline nodes
      for (int k = 0; k < mn; k++) {
        xs[k] = zeta(zn3[k], z1, re) / zgdif;
        ys[k] = 1.0 / tn3[k];
      }

      yd1 = -tgn3[0] / (t1 * t1) * zgdif;
      yd2 = -tgn3[1] / (t2 * t2) * zgdif * pow((re + z2) / (re + z1), 2);

      // Calculate spline coefficients
      hwm.spline(xs, ys, mn, yd1, yd2, y2out);
      x = zg / zgdif;
      hwm.splint(xs, ys, y2out, mn, x, y);

      // Temperature at altitude
      tz = 1.0 / y;

      if (xm != 0.0) {
        // Calculate tropospheric/stratosphere density
        glb = gsurf / pow(1.0 + z1 / re, 2);
        gamm = xm * glb * zgdif / rgas;

        // Integrate temperature profile
        splini(xs, ys, y2out, mn, x, yi);
        expll = gamm * yi;

        if (expll > 50) {
          expll = 50.0;
        }

        densm_out = densm_out * (t1 / (tz)) * exp(-expll);
      }
    }
  }

  if (xm == 0.0) {
    densm_out = tz;
  }

  return densm_out;
}

//! \brief Turbopause correction for MSIS model
//!
//! Root mean density
//! 8/20/80
//!
//! \param dd diffusive density
//! \param dm full mixed density
//! \param zhm transition scale length
//! \param xmm full mixed molecular weight
//! \param xm species molecular weight
//!
//! \returns  combined density
greal MSIS::dnet(greal dd, greal dm, greal zhm, greal xmm, greal xm)
{
  greal dnet_out = 0.0;

  greal a = zhm / (xmm - xm);

  if ((dm <= 0.0) && (dd <= 0.0)) {

    if ((dd != 0.0) && (dm != 0.0)) {
      dd = 1.0;
    }

    if (dm == 0.0) {
      dnet_out = dd;
    }
    else if (dd == 0.0) {
      dnet_out = dm;
    }
  }
  else {
    greal ylog = a * log(dm / dd);

    if (ylog < -10.0) {
      dnet_out = dd;
    }
    else if (ylog > 10.0) {
      dnet_out = dm;
    }
    else {
      dnet_out = dd * pow(1.0 + exp(ylog), 1 / a);
    }
  }
  return dnet_out;
}

//! \brief Chemistry/dissociation correction for MSIS models
//!
//! \param alt  altitude
//! \param r    target ratio
//! \param h1   transition scale length
//! \param zh   altitude of 1/2 r
//!
//! \returns Chemistry/dissociation correction for MSIS models.
greal MSIS::ccor(greal alt, greal r, greal h1, greal zh)
{
  greal ccor_exp;

  greal e = (alt - zh) / h1;

  if (e > 70.0) {
    ccor_exp = 0.0;
  }
  else if (e < -70.0) {
    ccor_exp = r;
  }
  else {
    greal ex = exp(e);
    ccor_exp = r / (1.0 + ex);
  }

  greal ccor_out = exp(ccor_exp);

  return ccor_out;
}

//! \brief O and O2 chemistry dissociation correction for MSIS models
//!
//! \param alt
//! \param r
//! \param h1
//! \param zh
//! \param h2
//!
//! \returns O and O2 chemistry dissociation correction for MSIS models.
greal MSIS::ccor2(greal alt, greal r, greal h1, greal zh, greal h2)
{
  greal e1 = (alt - zh) / h1;
  greal e2 = (alt - zh) / h2;

  greal ccor2_exp;
  if ((e1 > 70.0) || (e2 > 70.0)) {
    ccor2_exp = 0.0;
  }
  else if ((e1 < -70.0) && (e2 < -70.0)) {
    ccor2_exp = r;
  }
  else {
    greal ex1 = exp(e1);
    greal ex2 = exp(e2);
    ccor2_exp = r / (1.0 + 0.5 * (ex1 + ex2));
  }

  greal ccor2_out = exp(ccor2_exp);

  return ccor2_out;
}

//! \brief Convert outputs to kg and meters if meter true
//!
//! Not used.
void MSIS::meters(bool meter)
{
  imr = 0;

  if (meter) {
    imr = 1;
  }
}

//! \brief Integrate cubic spline function from xa[0] to x
//!
//! \param xa  array of tabulated function in ascending order by x
//! \param ya  array of tabulated function in ascending order by x
//! \param y2a array of second derivatives
//! \param n   size of arrays xa, ya, y2a
//! \param x   abscissa endpoint for integration
//! \param[out] yi A greal.
//!
//! \retval yi  output value
void MSIS::splini(greal xa[], greal ya[], greal y2a[], int n, greal x, greal &yi)
{

  // int klo, khi;
  // greal xx, h, a, b, a2, b2;

  yi = 0;
  int klo = 1;
  int khi = 2;

  while ((x > xa[klo - 1]) && (khi <= n)) {
    greal xx = x;
    if (khi < n) {
      xx = fmin(x, xa[khi - 1]);
    }

    greal h = xa[khi - 1] - xa[klo - 1];
    greal a = (xa[khi - 1] - xx) / h;
    greal b = (xx - xa[klo - 1]) / h;
    greal a2 = a * a;
    greal b2 = b * b;
    yi = yi
         + h
               * ((1.0 - a2) * ya[klo - 1] / 2.0 + b2 * ya[khi - 1] / 2.0
                  + ((-(1.0 + a2 * a2) / 4.0 + a2 / 2.0) * y2a[klo - 1]
                     + (b2 * b2 / 4.0 - b2 / 2.0) * y2a[khi - 1])
                        * h * h / 6.0);

    klo++;
    khi++;

    if (klo > 1000) {
      break;
    }
  }
}

//! \brief Calculates latitude variable gravity and effective radius
//!
//! \param lat Latitude.
//! \param[out] gv     A greal.
//! \param[out] refff  A greal.
//!
//! \retval gv  variable gravity
//! \retval refff  effective radius
void MSIS::glatf(greal lat, greal &gv, greal &refff)
{
  greal c2 = cos(2.0 * toRadians(lat));
  gv = 980.616 * (1.0 - 0.0026373 * c2);
  refff = 2.0 * gv / (3.085462e-6 + 2.27e-9 * c2) * 1e-5;
}

//! \brief Find altitude of pressure surface (press) from gtd7.
//!
//!	\param iyd   year and day as yyddd
//!	\param sec   ut(sec)
//! \param[out]	alt   altitude(km)
//!	\param glat   geodetic latitude(deg)
//!	\param glong   geodetic longitude(deg)
//!	\param stl   local apparent solar time(hrs)
//!	\param f107a   3 month average of f10.7 flux
//!	\param f107   daily f10.7 flux for previous day
//!	\param ap   magnetic index(daily) or when sw(9) = -1. :
//!		array containing:
//!		(1) daily ap
//!		(2) 3 hr ap index for current time
//!		(3) 3 hr ap index for 3 hrs before current time
//!		(4) 3 hr ap index for 6 hrs before current time
//!		(5) 3 hr ap index for 9 hrs before current time
//!		(6) average of eight 3 hr ap indicies from 12 to 33 hrs
//!		(7) average of eight 3 hr ap indicies from 36 to 59 hrs
//! \param[out] d A greal array.
//! \param[out] t A greal array.
//!	\param press   pressure level(mb)
//! 
//!	\retval alt - altitude(km)
//! \retval d 
//!	  \li d(1) - he number density(cm-3)
//!	  \li d(2) - o number density(cm-3)
//!	  \li d(3) - n2 number density(cm-3)
//!	  \li d(4) - o2 number density(cm-3)
//!	  \li d(5) - ar number density(cm-3)
//!	  \li d(6) - total mass density(gm/cm3)
//!	  \li d(7) - h number density(cm-3)
//!	  \li d(8) - n numner density(cm-3)
//!	  \li d(9) - hot o number density(cm-3)
//! \retval t
//!	  \li t(1) - exospheric temperature
//!	  \li t(2) - temperature at alt
void MSIS::ghp7(int iyd, greal sec, greal &alt, greal glat, greal glong, greal stl,
                greal f107a, greal f107, greal ap[], greal d[], greal t[], greal press)
{
  constexpr greal bm = BOLTZMANN * 1.0e4;
  constexpr greal rgas = UNIVERSAL_GAS * 100.0;
  constexpr greal test = 0.00043;
  constexpr int ltest = 12;

//  greal p, cl, cl2, cd, z, ca, xn, diff, xm, g, sh;
  greal z; 
  greal pl = log10(press);

  // Initial altitude estimate
  if (pl >= -5.0) {
    greal zi;
    if (pl > 2.5) {
      zi = 18.06 * (3.0 - pl);
    }
    else if (pl > 0.75 && pl <= 2.5) {
      zi = 14.98 * (3.08 - pl);
    }
    else if (pl > -1.0 && pl <= 0.75) {
      zi = 17.8 * (2.72 - pl);
    }
    else if (pl > -2.0 && pl <= -1.0) {
      zi = 14.28 * (3.64 - pl);
    }
    else if (pl > -4.0 && pl <= -2.0) {
      zi = 12.72 * (4.32 - pl);
    }
    else { // if (pl <= -4.0)
      zi = 25.3 * (0.11 - pl);
    }

    int iday = iyd % 1000;
    greal cl = glat / 90.0;
    greal cl2 = cl * cl;

    greal cd; 
    if (iday < 182) {
      cd = 1.0 - iday / 91.25;
    }
    else { //if (iday >= 182) {
      cd = iday / 91.25 - 3.0;
    }

    greal ca = 0.0;

    if (pl > -1.11 && pl <= -0.23) {
      ca = 1.0;
    }
    else if (pl > -0.23) {
      ca = (2.79 - pl) / (2.79 + 0.23);
    }
    else if (pl <= -1.11 && pl > -3.0) {
      ca = (-2.93 - pl) / (-2.93 + 1.11);
    }

    z = zi - 4.87 * cl * cd * ca - 1.64 * cl2 * ca + 0.31 * ca * cl;
  }
  else { //if (pl < -5.0) {
    z = 22.0 * pow(pl + 4.0, 2) + 110.0;
  }

  // Iteration Loop
  int l = 0;
  while (l <= ltest) {
    l++;

    gtd7(iyd, sec, z, glat, glong, stl, f107a, f107, ap, 48, d, t);

    greal xn = d[0] + d[1] + d[2] + d[3] + d[4] + d[5] + d[6] + d[7];

    greal p = bm * xn * t[1];

    if (imr == 1) {
      p = p * 1.0e-06;
    }

    greal diff = pl - log10(p);

    if (fabs(diff) < test || l == ltest) {
      break;
    }

    greal xm = d[5] / (xn * 1.66e-24);

    if (imr == 1) {
      xm = xm * 1.0e03;
    }

    greal g = gsurf / pow(1.0 + z / re, 2);
    greal sh = rgas * t[1] / (xm * g);

    // New altitude estimate using scale height
    if (l < 6) {
      z = z - sh * diff * 2.302;
    }
    else {
      z = z - sh * diff;
    }
  }

  alt = z;
}

//! \brief Evaluate MSIS and HWM model for specified position.
//!
//! \param iyr
//! \param z
//! \param xlat
//! \param xlon
//! \param f10b
//! \param f10
//! \param ap
//! \param[out] ph Pressure
//! \param[out] dh Density
//! \param[out] th Temperature
//! \param[out] n2nd    H2 number density.
//! \param[out] o2nd    O2 number density.
//! \param[out] ond     O number density.
//! \param[out] arnd    AR number density.
//! \param[out] hend    HE number density.
//! \param[out] hnd     H number density.
//! \param[out] nnd     N number density.
//! \param[out] wtmol Molecular weight
//! \param[out] tex Exospheric temperature
//! \param[out] uh E/W winds.
//! \param[out] vh N/S winds.
void MSIS::msishwm(int iyr, greal z, greal xlat,
                   greal xlon, greal f10b, greal f10, greal ap[], greal &ph, greal &dh,
                   greal &th, greal &n2nd, greal &o2nd, greal &ond, greal &arnd, greal &hend,
                   greal &hnd, greal &nnd, greal &wtmol, greal &tex, greal &uh, greal &vh)
{
  // Molecular weights
  const greal ei[9] = {4.0026, 15.9994, 28.0134, 31.9988, 39.948, 0.0, 1.00797, 14.0067, 15.9994};
  const greal r0 = UNIVERSAL_GAS * 1000.0;
  greal ddd = dayOfYear;
  int iday = 1000 * iyr + int(ddd + 1.0);
  greal dut = 86400 * (ddd - int(ddd));
  greal utsec = dut;
  greal dlon = xlon - longitude;
  greal xlst = solarTime + dlon / 15.0;

  // Call MSIS thermosphere subroutine
  greal d[9], t[2];
  gtd7(iday, utsec, z, xlat, xlon, xlst, f10b, f10, ap, 48, d, t);

  // MSIS output variables:
  //  d[0] - He number density                (cm-3)
  //  d[1] - O number density                 (cm-3)
  //  d[2] - N2 number density                (cm-3)
  //  d[3] - O2 number density                (cm-3)
  //  d[4] - Ar number density                (cm-3)
  //  d[5] - Total mass density               (GM/cm3)
  //  d[6] - H number density                 (cm-3)
  //  d[7] - N number density                 (cm-3)
  //  d[8] - Anomalous oxygen number density  (cm-3)
  //  t[0] - exospheric temperature
  //  t[1] - temperature at alt

  greal sumn = 0.0;
  greal summn = 0.0;

  // Convert MSIS outputs to SI units
  for (int i = 0; i < 9; i++) {
    d[i] = 1.0e06 * d[i];

    if (i != 5) {
      sumn = sumn + d[i];
    }

    summn = summn + d[i] * ei[i];
  }

  // Change notation for outputs
  dh = d[5] / 1000.0;
  th = t[1];
  tex = t[0];
  wtmol = summn / sumn;
  ph = dh * r0 * (th) / (wtmol);
  hend = d[0];
  ond = d[1] + d[8];
  n2nd = d[2];
  o2nd = d[3];
  arnd = d[4];
  hnd = d[6];
  nnd = d[7];

  // Call MSIS Harmonic Wind Model (HWM) subroutine
  greal w[2] = {0.0};
  hwm.gws5(iday, utsec, z, xlat, xlon, xlst, f10b, f10, ap, w);

  // Change notation for outputs
  uh = w[1];
  vh = w[0];
}

//! \brief Neutral Atmosphere Empirical Model fro the surface to lower exosphere (NRLMSISE-00)
//!
//! \param iyd
//! \param sec
//! \param alt
//! \param glat
//! \param glong
//! \param stl
//! \param f107a
//! \param f107
//! \param ap
//! \param mass
//! \param[out] d
//! \param[out] t
void MSIS::gtd7(int iyd, greal sec, greal alt, greal glat, greal glong, greal stl,
                greal f107a, greal f107, greal ap[], greal mass, greal d[], greal t[])
{
  const int mn2 = 4, mn3 = 5;
  const greal zmix = 62.5;

  // alast should probably be a member variable, but we are forcing thermosphere
  // calculations on every pass.  Variable is retained to match original code.
  greal alast = 99999.0; 
  greal zn2[4] = { 72.5, 55.0, 45.0, 32.5 }, zn3[5] = { 32.5, 20.0, 15.0, 10.0, 0.0 };
  //greal sv[25] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
  //                1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  greal  ds[9], ts[2], dm28m, tz;
  int mssl = -999;

  // Test for changed input
  //int v1 = vtst7(iyd, sec, glat, glong, stl, f107a, f107, ap, 1);
  // Force thermosphere computations
  int v1 = 1;

  // Latitude variation of gravity (none for getCswSw(1) = 0.0
  greal xlat = glat;

  if (getCswSw(1) == 0.0) {
    xlat = 45.0;
  }

  glatf(xlat, gsurf, re);

  greal xmm = pdm[4][2];

  // Thermosphere/mesosphere (above zn2[0])
  greal altt = fmax(alt, zn2[0]);
  int mss = int(mass);

  // Only calculate thermosphere if input parameters changed or altitude
  // above zn2[0] in mesosphere
  if (v1 == 1 || alt > zn2[0] || alast > zn2[0] || mss != mssl) {
    gts7(iyd, sec, altt, glat, glong, stl, f107a, f107, ap, mss, dm28m, ds, ts);

    // Metric adjustment
    if (imr == 1) {
      dm28m = dm28m * 1.0e06;
    }

    mssl = mss;
  }

  t[0] = ts[0];
  t[1] = ts[1];

  if (alt >= zn2[0]) {

    for (int j = 0; j < 9; j++) {
      d[j] = ds[j];
    }

    alast = alt;
    return;
  }

  // The rest of the code in this method never gets called.
  // Only calculate nodes if input changed
  if (v1 == 1.0 || alast >= zn2[0]) {
    tgn2[0] = tgn1[1];
    tn2[0] = tn1[4];
    tn2[1] = pma[0][0] * pavgm[0] / (1.0 - getCswSw(19) * glob7s(&pma[0][0]));
    tn2[2] = pma[1][0] * pavgm[1] / (1.0 - getCswSw(19) * glob7s(&pma[1][0]));
    tn2[3] = pma[2][0] * pavgm[2] / (1.0 - getCswSw(19) * getCswSw(21) * glob7s(&pma[2][0]));
    tgn2[1] = pavgm[8] * pma[9][0] * (1.0 + getCswSw(19) * getCswSw(21) * glob7s(&pma[9][0]))
              * tn2[3] * tn2[3] / pow(pma[2][0] * pavgm[2], 2);
    tn3[0] = tn2[3];
  }

  if (alt < zn3[0]) {
    // LOWER STRATOSPHERE AND TROPOSPHERE [below ZN3(1)]
    // Temperature at nodes and gradients at end nodes
    // Inverse temperature a linear function of spherical harmonics
    // Only calculate nodes if input changed
    if (v1 == 1.0 || alast >= zn3[0]) {
      tgn3[0] = tgn2[1];
      tn3[1] = pma[3][0] * pavgm[3] / (1.0 - getCswSw(21) * glob7s(&pma[3][0]));
      tn3[2] = pma[4][0] * pavgm[4] / (1.0 - getCswSw(21) * glob7s(&pma[4][0]));
      tn3[3] = pma[5][0] * pavgm[5] / (1.0 - getCswSw(21) * glob7s(&pma[5][0]));
      tn3[4] = pma[6][0] * pavgm[6] / (1.0 - getCswSw(21) * glob7s(&pma[6][0]));
      tgn3[1] = pma[7][0] * pavgm[7] * (1.0 + getCswSw(21) * glob7s(&pma[7][0])) * tn3[4] * tn3[4]
                / pow(pma[6][0] * pavgm[6], 2);
    }
  }

  if (mass == 0) {
    dd = densm(alt, 1.0, 0.0, tz, mn3, zn3, tn3, tgn3, mn2, zn2, tn2, tgn2);

    t[1] = tz;
    alast = alt;
    return;
  }

  // Linear transition to full mixing below zn2[0]
  greal dmc = 0.0;

  if (alt > zmix) {
    dmc = 1.0 - (zn2[0] - alt) / (zn2[0] - zmix);
  }

  greal dz28 = ds[2];

  // N2 Density
  greal dmr = ds[2] / dm28m - 1.0;
  d[2] = densm(alt, dm28m, xmm, tz, mn3, zn3, tn3, tgn3, mn2, zn2, tn2, tgn2);
  d[2] = d[2] * (1.0 + dmr * dmc);

  // He Density
  d[0] = 0.0;

  if (mass == 4 || mass == 48) {
    dmr = ds[0] / (dz28 * pdm[1][0]) - 1.0;
    d[0] = d[2] * pdm[1][0] * (1.0 + dmr * dmc);
  }

  // O Density
  d[1] = 0.0;
  d[8] = 0.0;

  // O2 Density
  d[3] = 0.0;

  if (mass == 32 || mass == 48) {
    dmr = ds[3] / (dz28 * pdm[1][3]) - 1.0;
    d[3] = d[2] * pdm[1][3] * (1.0 + dmr * dmc);
  }

  // Ar Density
  d[4] = 0.0;

  if (mass == 40 || mass == 48) {
    dmr = ds[4] / (dz28 * pdm[1][4]) - 1.0;
    d[4] = d[2] * pdm[1][4] * (1.0 + dmr * dmc);
  }

  // H Density
  d[6] = 0.0;

  // N Density
  d[7] = 0.0;

  // Total Mass Density
  if (mass == 48) {
    d[5] = 1.66e-24 
           * (4.0 * d[0] + 16.0 * d[1] + 28.0 * d[2] + 32.0 * d[3] + 40.0 * d[4] + d[6]
              + 14.0 * d[7]);
  }

  if (imr == 1) {
    d[5] = d[5] / 1000.0;
  }

  t[1] = tz;

  alast = alt;
}

//! \brief Thermospheric portion of NRLMSISE-00
//!
//! \param iyd
//! \param sec
//! \param alt
//! \param glat
//! \param glong
//! \param stl
//! \param f107a
//! \param f107
//! \param ap
//! \param mass
//! \param[out] dm28
//! \param[out] d
//! \param[out] t
void MSIS::gts7(int iyd, greal sec, greal alt, greal glat, greal glong, greal stl, greal f107a,
                greal f107, greal ap[], int mass, greal &dm28, greal d[], greal t[])
{
  const greal altl[8] = { 200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0 };
  const greal dr = 1.72142e-02;
  const greal alpha[9] = {-0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0};
  const int mn1 = 5;
  // Mass is always 48 in EG, so this is not used.
  //const int mt[11] = { 48, 0, 4, 16, 28, 32, 40, 1, 49, 14, 17 };

  greal zn1[5] = {120.0, 110.0, 100.0, 90.0, 72.5};
  greal alast = -999.0;
  greal tlb = 0.0;
  
  greal g28, yday, zhf, z, zh28, zhm28, b28, xmd, tz, xmm, g4, zh04, b04, zhm04, zc04, hc04,
      g16, zh16, b16, zhm16, hc16, zc16, hc216, hcc16, zcc16, rc16, g32, zh32, b32, zhm32, hc32,
      zc32, hcc32, hcc232, zcc32, rc32, g40, zh40, b40, zhm40, hc40, zc40, g1, zh01, b01, zhm01,
      hc01, zc01, hcc01, zcc01, rc01, g14, zh14, b14, zhm14, hc14, zc14, hcc14, zcc14, rc14, g16h,
      db16h, tho, t2, zsht, zhmo, zsho, dm04, dm16, dm32, dm40, dm01, dm14, db04, db16, db28, db32,
      db40, db01, rl, db14;
  //greal db48;  // set but not used
  int lenyr, iyr;

  b28 = 0.0;
  zhm28 = 0.0;

  glatf(glat, gsurf, re);

  // Test for changed input
  //int v2 = vtst7(iyd, sec, glat, glong, stl, f107a, f107, ap, 2);
  // Force computations.
  int v2 = 1;

  int yrd = iyd;

  greal za = pdl[15][1];
  zn1[0] = za;

  for (int i = 0; i < 9; i++) {
    d[i] = 0.0;
  }

  // Tinf variations not important below za or zn1[0]
  greal tinf;
  if (alt > zn1[0]) {

    if ((v2 == 1) || (alast <= zn1[0])) {
      tinf = ptm[0] * pt[0]
             * (1.0 + getCswSw(15) * globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, pt));
    }
  }
  else {
    tinf = ptm[0] * pt[0];
  }

  t[0] = tinf;

  // Gradient variations not impontant below zn1[4]
  greal g0;
  if (alt > zn1[4]) {
    if ((v2 == 1) || (alast <= zn1[4])) {
      g0 = ptm[3] * ps[0]
           * (1.0 + getCswSw(18) * globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, ps));
    }
  }
  else {
    g0 = ptm[3] * ps[0];
  }

  // Calculate these temperatures only if input changed
  if ((v2 == 1) || (alt < 300.0)) {
    tlb = ptm[1]
          * (1.0 + getCswSw(16) * globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, &pd[3][0]))
          * pd[3][0];
  }
  greal s = g0 / (tinf - tlb);

  // Lower themosphere temp variations not significant for density above 300 km
  if (alt < 300.0) {
    if (v2 == 1 || alast >= 300.0) {
      tn1[1] = ptm[6] * ptl[0][0] / (1.0 - getCswSw(17) * glob7s(&ptl[0][0]));
      tn1[2] = ptm[2] * ptl[1][0] / (1.0 - getCswSw(17) * glob7s(&ptl[1][0]));
      tn1[3] = ptm[7] * ptl[2][0] / (1.0 - getCswSw(17) * glob7s(&ptl[2][0]));
      tn1[4] = ptm[4] * ptl[3][0] / (1.0 - getCswSw(17) * getCswSw(19) * glob7s(&ptl[3][0]));
      tgn1[1] = ptm[8] * pma[8][0] * (1.0 + getCswSw(17) * getCswSw(19) * glob7s(&pma[8][0]))
                * tn1[4] * tn1[4] / pow(ptm[4] * ptl[3][0], 2);
    }
  }
  else {
    tn1[1] = ptm[6] * ptl[0][0];
    tn1[2] = ptm[2] * ptl[1][0];
    tn1[3] = ptm[7] * ptl[2][0];
    tn1[4] = ptm[4] * ptl[3][0];
    tgn1[1] = ptm[8] * pma[8][0] * tn1[4] * tn1[4] / pow(ptm[4] * ptl[3][0], 2);
  }

  //greal z0 = zn1[3];  // set but not used
  //greal t0 = tn1[3];  // set but not used
  //greal tr12 = 1.0;   // set but not used

  // N2 variation factor at Zlb
  g28 = getCswSw(20) * globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, &pd[2][0]);
  yday = fmod(yrd, 1000.0);
  iyr = int(yrd) / 1000;

  // Adjust day-of-year dependency for leap year vs non-leap year
  lenyr = 365;

  if (iyr % 4 == 0) {
    lenyr = 366;
  }

  yday = (yday - 1.0 + sec / 86400.0) * (365.0 / lenyr);

  // Variation of turbopause height
  zhf = pdl[24][1] * (1.0 + getCswSw(4) * pdl[24][0] * sin(toRadians(glat)) * cos(dr * (yday - pt[13])));
  yrd = iyd;
  t[0] = tinf;
  xmm = pdm[4][2];
  z = alt;

  if ((z <= altl[5]) || (mass == 28) || (mass == 48)) {
    // N2 Density
    // Diffusive density at zlb
    db28 = pdm[0][2] * exp(g28) * pd[2][0];

    // Diffusive density at alt
    d[2] = densu(z, db28, tinf, tlb, 28.0, alpha[2], t[1], ptm[5], s, mn1, zn1, tn1, tgn1);

    dd = d[2];

    // Turbopause
    zh28 = pdm[2][2] * zhf;
    zhm28 = pdm[3][2] * pdl[5][1];
    xmd = 28.0 - xmm;

    // Mixed density at zlb
    b28 = densu(zh28, db28, tinf, tlb, xmd, alpha[2] - 1.0, tz, ptm[5], s, mn1, zn1, tn1, tgn1);

    if ((z <= altl[2]) && (getCswSw(14) != 0.0)) {
      // Mixed density at alt
      dm28 = densu(z, b28, tinf, tlb, xmm, alpha[2], tz, ptm[5], s, mn1, zn1, tn1, tgn1);

      // Net density at alt
      d[2] = dnet(d[2], dm28, zhm28, xmm, 28.0);
    }
  }


  // He Density
  // Density variation factor at zlb
  g4 = getCswSw(20) * globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, &pd[0][0]);

  // Diffusive density at zlb
  db04 = pdm[0][0] * exp(g4) * pd[0][0];

  // Diffusive density at alt
  d[0] = densu(z, db04, tinf, tlb, 4.0, alpha[0], t[1], ptm[5], s, mn1, zn1, tn1, tgn1);
  dd = d[0];

  if ((z <= altl[0]) && (getCswSw(14) != 0.0)) {
    // Turbopause
    zh04 = pdm[2][0];

    // Mixed density at zlb
    b04 = densu(zh04, db04, tinf, tlb, 4.0 - xmm, alpha[0] - 1.0, t[1], ptm[5], s, mn1, zn1, tn1, tgn1);

    // Mixed density at alt
    dm04 = densu(z, b04, tinf, tlb, xmm, 0.0, t[1], ptm[5], s, mn1, zn1, tn1, tgn1);
    zhm04 = zhm28;
    // Net density at alt
    d[0] = dnet(d[0], dm04, zhm04, xmm, 4.0);

    // Correction to specified mixing ratio at ground
    rl = log(b28 * pdm[1][0] / b04);
    zc04 = pdm[4][0] * pdl[0][1];
    hc04 = pdm[5][0] * pdl[1][1];

    // Net density corrected at alt
    d[0] = d[0] * ccor(z, rl, hc04, zc04);
  }


  // O Density
  g16 = getCswSw(20) * globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, &pd[1][0]);

  // Diffusive density at zlb
  db16 = pdm[0][1] * exp(g16) * pd[1][0];

  // Diffusive density at Alt
  d[1] = densu(z, db16, tinf, tlb, 16.0, alpha[1], t[1], ptm[5], s, mn1, zn1, tn1, tgn1);

  dd = d[1];

  if (z <= altl[1] && getCswSw(14) != 0.0) {
    // Corrected from low.pdm[2][0] to low.pdm[2][1] 12/2/85
    // Turbopause
    zh16 = pdm[2][1];

    // Mixed density at zlb
    b16 = densu(zh16, db16, tinf, tlb, 16.0 - xmm, alpha[1] - 1.0, t[1], ptm[5], s, mn1, zn1, tn1, tgn1);

    // Mixed density at alt
    dm16 = densu(z, b16, tinf, tlb, xmm, 0.0, t[1], ptm[5], s, mn1, zn1, tn1, tgn1);
    zhm16 = zhm28;

    // Net density at alt
    d[1] = dnet(d[1], dm16, zhm16, xmm, 16.0);

    // 3/16/99 change form to match O2 departure from diff equil near 150 km
    // and add dependence on F10.7
    rl = pdm[1][1] * pdl[16][1] * (1.0 + getCswSw(0) * pdl[23][0] * (f107a - 150.0));
    hc16 = pdm[5][1] * pdl[3][1];
    zc16 = pdm[4][1] * pdl[2][1];
    hc216 = pdm[5][1] * pdl[4][1];
    d[1] = d[1] * ccor2(z, rl, hc16, zc16, hc216);

    // Chemistry correction
    hcc16 = pdm[7][1] * pdl[13][1];
    zcc16 = pdm[6][1] * pdl[12][1];
    rc16 = pdm[3][1] * pdl[14][1];

    // Net density corrected at alt
  d[1] = d[1] * ccor(z, rc16, hcc16, zcc16);
  }


  // O2 Density
  // Density variation factor at zlb
  g32 = getCswSw(20) * globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, &pd[4][0]);

  // Diffusive density at zlb
  db32 = pdm[0][3] * exp(g32) * pd[4][0];

  // Diffusive density at alt
  d[3] = densu(z, db32, tinf, tlb, 32.0, alpha[3], t[1], ptm[5], s, mn1, zn1, tn1, tgn1);

  if (mass == 49) {
    dd = dd + 2.0 * d[3];
  }
  else {
    dd = d[3];
  }

  if (getCswSw(14) != 0.0) {

    if (z <= altl[3]) {
      // Turbopause
      zh32 = pdm[2][3];

      // Mixed density at zlb
      b32 = densu(zh32, db32, tinf, tlb, 32.0 - xmm, alpha[3] - 1.0, t[1], ptm[5], s, mn1, zn1, tn1, tgn1);

      // Mixed density at alt
      dm32 = densu(z, b32, tinf, tlb, xmm, 0.0, t[1], ptm[5], s, mn1, zn1, tn1, tgn1);

      zhm32 = zhm28;

      // Net density at Alt
      d[3] = dnet(d[3], dm32, zhm32, xmm, 32.0);

      // Correction to specified mixing ratio at ground
      rl = log(b28 * pdm[1][3] / b32);
      hc32 = pdm[5][3] * pdl[7][1];
      zc32 = pdm[4][3] * pdl[6][1];
      d[3] = d[3] * ccor(z, rl, hc32, zc32);
    }
  
    // Correction for general deprature from diffusive equilibrium above zlb
    hcc32 = pdm[7][3] * pdl[22][1];
    hcc232 = pdm[7][3] * pdl[22][0];
    zcc32 = pdm[6][3] * pdl[21][1];
    rc32 = pdm[3][3] * pdl[23][1] * (1.0 + getCswSw(0) * pdl[23][0] * (f107a - 150.0));

    // Net density corrected at alt
    d[3] = d[3] * ccor2(z, rc32, hcc32, zcc32, hcc232);
  }


  // AR Density
  // Density variation factor at zlb
  g40 = getCswSw(20) * globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, &pd[5][0]);

  // Diffusive density at zlb
  db40 = pdm[0][4] * exp(g40) * pd[5][0];

  // Diffusive density at alt
  d[4] = densu(z, db40, tinf, tlb, 40.0, alpha[4], t[1], ptm[5], s, mn1, zn1, tn1, tgn1);
  dd = d[4];

  if ((z <= altl[4]) && (getCswSw(14) != 0.0)) {
    // Turbopause
    zh40 = pdm[2][4];

    // Mixed density at zlb
    b40 = densu(zh40, db40, tinf, tlb, 40.0 - xmm, alpha[4] - 1.0, t[1], ptm[5], s, mn1, zn1, tn1, tgn1);

    // Mixed density at alt
    dm40 = densu(z, b40, tinf, tlb, xmm, 0.0, t[1], ptm[5], s, mn1, zn1, tn1, tgn1);
    zhm40 = zhm28;

    // Net density at alt
    d[4] = dnet(d[4], dm40, zhm40, xmm, 40.0);

    // Correction to specified mixing ratio at ground
    rl = log(b28 * pdm[1][4] / b40);
    hc40 = pdm[5][4] * pdl[9][1];
    zc40 = pdm[4][4] * pdl[8][1];

    // Net density corrected at alt
    d[4] = d[4] * ccor(z, rl, hc40, zc40);
  }


  // Hydrogen Density
  // Density variation factor at zlb

  g1 = getCswSw(20) * globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, &pd[6][0]);

  // Diffusive density at zlb
  db01 = pdm[0][5] * exp(g1) * pd[6][0];

  // Diffusive density at alt
  d[6] = densu(z, db01, tinf, tlb, 1.0, alpha[6], t[1], ptm[5], s, mn1, zn1, tn1, tgn1);
  dd = d[6];

  if ((z <= altl[6]) && (getCswSw(14) != 0.0)) {
    // Turbopause
    zh01 = pdm[2][5];

    // Mixed density at zlb
    b01 = densu(zh01, db01, tinf, tlb, 1.0 - xmm, alpha[6] - 1.0, t[1], ptm[5], s, mn1, zn1, tn1,
                tgn1);

    // Mixed density at alt
    dm01 = densu(z, b01, tinf, tlb, xmm, 0.0, t[1], ptm[5], s, mn1, zn1, tn1, tgn1);
    zhm01 = zhm28;
    // Net density at alt
    d[6] = dnet(d[6], dm01, zhm01, xmm, 1.0);

    // Correction to specified mixing ratio at ground
    rl = log(b28 * pdm[1][5] * fabs(pdl[17][1]) / b01);
    hc01 = pdm[5][5] * pdl[11][1];
    zc01 = pdm[4][5] * pdl[10][1];
    d[6] = d[6] * ccor(z, rl, hc01, zc01);

    // Chemistry correction
    hcc01 = pdm[7][5] * pdl[19][1];
    zcc01 = pdm[6][5] * pdl[18][1];
    rc01 = pdm[3][5] * pdl[20][1];

    // Net density corrected at alt
    d[6] = d[6] * ccor(z, rc01, hcc01, zcc01);
  }


  // Atomic Nitrogen Density
  // Density variation at zlb
  g14 = getCswSw(20) * globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, &pd[7][0]);

  // Diffusive density at zlb
  db14 = pdm[0][6] * exp(g14) * pd[7][0];

  // Diffusive density at alt
  d[7] = densu(z, db14, tinf, tlb, 14.0, alpha[7], t[1], ptm[5], s, mn1, zn1, tn1, tgn1);
  dd = d[7];

  if ((z <= altl[7]) && (getCswSw(14) != 0.0)) {
    // Turbopause
    zh14 = pdm[2][6];

    // Mixed density at zlb
    b14 = densu(zh14, db14, tinf, tlb, 14.0 - xmm, alpha[7] - 1.0, t[1], ptm[5], s, mn1, zn1, tn1, tgn1);

    // Mixed density at alt
    dm14 = densu(z, b14, tinf, tlb, xmm, 0.0, t[1], ptm[5], s, mn1, zn1, tn1, tgn1);
    zhm14 = zhm28;
    // Net density at alt
    d[7] = dnet(d[7], dm14, zhm14, xmm, 14.0);

    // Correction to specified mixing ratio at ground
    rl = log(b28 * pdm[1][6] * fabs(pdl[2][0]) / b14);
    hc14 = pdm[5][6] * pdl[1][0];
    zc14 = pdm[4][6] * pdl[0][0];
    d[7] = d[7] * ccor(z, rl, hc14, zc14);

    // Chemistry correction
    hcc14 = pdm[7][6] * pdl[4][0];
    zcc14 = pdm[6][6] * pdl[3][0];
    rc14 = pdm[3][6] * pdl[5][0];

    // Net density corrected at alt
    d[7] = d[7] * ccor(z, rc14, hcc14, zcc14);
  }

  // Anomalous Oxygen Density
  g16h = getCswSw(20) * globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, &pd[8][0]);
  db16h = pdm[0][7] * exp(g16h) * pd[8][0];
  tho = pdm[9][7] * pdl[6][0];
  dd = densu(z, db16h, tho, tho, 16.0, alpha[8], t2, ptm[5], s, mn1, zn1, tn1, tgn1);
  zsht = pdm[5][7];
  zhmo = pdm[4][7];
  zsho = scalh(zhmo, 16.0, tho);
  d[8] = dd * exp(-zsht / zsho * (exp(-(z - zhmo) / zsht) - 1.0));


  // Total mass Density
  d[5] = 1.66e-24 * (4.0 * d[0] + 16.0 * d[1] + 28.0 * d[2] + 32.0 * d[3] + 40.0 * d[4] + d[6] + 14.0 * d[7]);
  //db48 = 1.66e-24 * (4.0 * db04 + 16.0 * db16 + 28.0 * db28 + 32.0 * db32 + 40.0 * db40 + db01 + 14.0 * db14);

//  // Temperature at altitude
//  z = fabs(alt);
//  ddum = densu(z, 1.0, tinf, tlb, 0.0, 0.0, t[1], ptm[5], s, mn1, zn1, tn1, tgn1);

  // Adjust densities from cgs to kgm
  if (imr == 1) {
    for (int i = 0; i < 9; i++) {
      d[i] = d[i] * 1.0e+06;
    }
    d[5] = d[5] / 1000.0;
  }

  alast = alt;
}

//! \brief Calculate scale height.
//!
//! \param alt
//! \param xm
//! \param temp
//!
//! \returns  scale height \units{km}
greal MSIS::scalh(greal alt, greal xm, greal temp)
{
  constexpr greal rgas = UNIVERSAL_GAS * 100.0;

  greal g = gsurf / pow(1.0 + alt / re, 2);

  return rgas * temp / (g * xm);
}

//! \brief Calculate g(l) function Upper Thermosphere Parameters
//!
//! \param yrd
//! \param sec
//! \param lat
//! \param lon
//! \param tloc
//! \param f107a
//! \param f107
//! \param ap
//! \param[out] p
//!
//! \returns tinf
greal MSIS::globe7(greal yrd, greal sec, greal lat, greal lon, greal tloc, greal f107a,
                    greal f107, greal ap[], const greal p[])
{
  const greal dr = 1.72142e-02, hr = 0.2618, sr = 7.2722e-05;
  const int nsw = 14;

  // These values, used to speed up code, should be member variables to enable optimizations.
  // They are left as local variables to force calculations.
  greal xl = 1000.0, tll = 1000.0, dayl = -1.0, p14 = -1000.0, p18 = -1000.0, p32 = -1000.0, p39 = -1000.0;

  greal sv[25] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  hwm.tselec(sv);

  greal sw9 = 1.0;
  //if (getCswSw(8) > 0) {
  //  sw9 = 1.0;
  //}
  if (getCswSw(8) < 0) {
    sw9 = -1.0;
  }

  int iyr = int(yrd / 1000.0);
  day_yr = yrd - iyr * 1000.0;

  // Adjust day-of-year dependency of coefficients for leap year vs non-leap year
  int lenyr = 365;
  if (iyr % 4 == 0) {
    lenyr = 366;
  }

  day_yr = (day_yr - 1.0 + sec / 86400.0) * (365.0 / lenyr);

  xlong = lon;

  // Eq. A22 (remainder of code)
  if (xl != lat) {
    // Calculate Legendre polynomials
    greal c = sin(toRadians(lat));
    greal s = cos(toRadians(lat));
    greal c2 = c * c;
    greal c4 = c2 * c2;
    greal s2 = s * s;
    plg[1][0] = c;
    plg[2][0] = 0.5 * (3.0 * c2 - 1.0);
    plg[3][0] = 0.5 * (5.0 * c * c2 - 3.0 * c);
    plg[4][0] = (35.0 * c4 - 30.0 * c2 + 3.0) / 8.0;
    plg[5][0] = (63.0 * c2 * c2 * c - 70.0 * c2 * c + 15.0 * c) / 8.0;
    plg[6][0] = (11.0 * c * plg[5][0] - 5.0 * plg[4][0]) / 6.0;
    plg[1][1] = s;
    plg[2][1] = 3.0 * c * s;
    plg[3][1] = 1.5 * (5.0 * c2 - 1.0) * s;
    plg[4][1] = 2.5 * (7.0 * c2 * c - 3.0 * c) * s;
    plg[5][1] = 1.875 * (21.0 * c4 - 14.0 * c2 + 1.0) * s;
    plg[6][1] = (11.0 * c * plg[5][1] - 6.0 * plg[4][1]) / 5.0;
    plg[2][2] = 3.0 * s2;
    plg[3][2] = 15.0 * s2 * c;
    plg[4][2] = 7.5 * (7.0 * c2 - 1.0) * s2;
    plg[5][2] = 3.0 * c * plg[4][2] - 2.0 * plg[3][2];
    plg[6][2] = (11.0 * c * plg[5][2] - 7.0 * plg[4][2]) / 4.0;
    plg[7][2] = (13.0 * c * plg[6][2] - 8.0 * plg[5][2]) / 5.0;
    plg[3][3] = 15.0 * s2 * s;
    plg[4][3] = 105.0 * s2 * s * c;
    plg[5][3] = (9.0 * c * plg[4][3] - 7.0 * plg[3][3]) / 2.0;
    plg[6][3] = (11.0 * c * plg[5][3] - 8.0 * plg[4][3]) / 3.0;
    xl = lat;
  }

  if ((tll != tloc) || ((getCswSw(6) != 0.0) && (getCswSw(7) != 0.0) && (getCswSw(13) != 0.0))) {
    stloc = sin(hr * tloc);
    ctloc = cos(hr * tloc);
    s2tloc = sin(2.0 * hr * tloc);
    c2tloc = cos(2.0 * hr * tloc);
    s3tloc = sin(3.0 * hr * tloc);
    c3tloc = cos(3.0 * hr * tloc);
    tll = tloc;
  }

  greal cd14 = 0.0, cd18 = 0.0, cd32 = 0.0, cd39 = 0.0;
  if ((day_yr != dayl) || (p[13] != p14)) {
    cd14 = cos(dr * (day_yr - p[13]));
  }

  if ((day_yr != dayl) || (p[17] != p18)) {
    cd18 = cos(2.0 * dr * (day_yr - p[17]));
  }

  if ((day_yr != dayl) || (p[31] != p32)) {
    cd32 = cos(dr * (day_yr - p[31]));
  }

  if ((day_yr != dayl) || (p[38] != p39)) {
    cd39 = cos(2.0 * dr * (day_yr - p[38]));
  }

  dayl = day_yr;
  p14 = p[13];
  p18 = p[17];
  p32 = p[31];
  p39 = p[38];

  // f10.7 effect
  greal t[15] = {0.0};
  greal df = f107 - f107a;
  dfa = f107a - 150.0;
  t[0] = p[19] * df * (1.0 + p[59] * dfa) + p[20] * df * df + p[21] * dfa + p[29] * pow(dfa, 2.0);
  greal f1 = 1.0 + (p[47] * dfa + p[19] * df + p[20] * df * df) * getCswSwc(0);
  greal f2 = 1.0 + (p[49] * dfa + p[19] * df + p[20] * df * df) * getCswSwc(0);

  // Time independent
  t[1] = (p[1] * plg[2][0] + p[2] * plg[4][0] + p[22] * plg[6][0])
         + (p[14] * plg[2][0]) * dfa * getCswSwc(0) + p[26] * plg[1][0];

  // Symmetrical annual
  t[2] = (p[18]) * cd32;

  // Symmetrical semiannual
  t[3] = (p[15] + p[16] * plg[2][0]) * cd18;

  // Asymmetrical annual
  t[4] = f1 * (p[9] * plg[1][0] + p[10] * plg[3][0]) * cd14;

  // Asymmetrical semiannual
  t[5] = p[37] * plg[1][0] * cd39;

  // Diurnal
  if (getCswSw(6) != 0.0) {
    greal t71 = (p[11] * plg[2][1]) * cd14 * getCswSwc(4);
    greal t72 = (p[12] * plg[2][1]) * cd14 * getCswSwc(4);
    t[6] = f2 * ( (p[3] * plg[1][1] + p[4] * plg[3][1] + p[27] * plg[5][1] + t71) * (ctloc)
                + (p[6] * plg[1][1] + p[7] * plg[3][1] + p[28] * plg[5][1] + t72) * (stloc));
  }

  // Semidiurnal
  if (getCswSw(7) != 0.0) {
    greal t81 = (p[23] * plg[3][2] + p[35] * plg[5][2]) * cd14 * getCswSwc(4);
    greal t82 = (p[33] * plg[3][2] + p[36] * plg[5][2]) * cd14 * getCswSwc(4);

    t[7] = f2 * ( (p[5] * plg[2][2] + p[41] * plg[4][2] + t81) * (c2tloc)
                + (p[8] * plg[2][2] + p[42] * plg[4][2] + t82) * (s2tloc));
  }

  // Terdiurnal
  if (getCswSw(13) != 0.0) {
    t[13] = f2 * ( (p[39] * plg[3][3] + (p[93] * plg[4][3] + p[46] * plg[6][3]) * cd14 * getCswSwc(4)) * s3tloc
                 + (p[40] * plg[3][3] + (p[94] * plg[4][3] + p[48] * plg[6][3]) * cd14 * getCswSwc(4)) * c3tloc);
  }

  // Magnetic activity based on daily ap
  if (sw9 != -1.0) {
    greal apd = ap[0] - 4.0;
    greal p44 = p[43];
    greal p45 = p[44];

    if (p44 < 0.0) {
      p44 = 1.0e-05;
    }

    apdf = apd + (p45 - 1.0) * (apd + (exp(-p44 * apd) - 1.0) / p44);

    if (0.0 != getCswSw(8)) {
      t[8] = apdf * (p[32] + p[45] * plg[2][0] + p[34] * plg[4][0]
                    + (p[100] * plg[1][0] + p[101] * plg[3][0] + p[102] * plg[5][0]) * cd14 * getCswSwc(4)
                    + (p[121] * plg[1][1] + p[122] * plg[3][1] + p[123] * plg[5][1]) * getCswSwc(6)
                       * cos(hr * (tloc - p[124])));
    }
  }
  else {
    if (p[51] != 0.0) {

      greal exp1 = exp(-10800.0 * fabs(p[51]) / (1.0 + p[138] * (45.0 - fabs(lat))));

      if (0.99999 < exp1) {
        exp1 = 0.99999;
      }

      //if (1.0e-04 > p[24]) {
      //  p[24] = 1.0e-04;
      //}

      apt[0] = sg0(exp1, ap, p);

      if (getCswSw(8) != 0.0) {
        t[8] = apt[0] * (p[50] + p[96] * plg[2][0] + p[54] * plg[4][0]
                        + (p[125] * plg[1][0] + p[126] * plg[3][0] + p[127] * plg[5][0]) * cd14 * getCswSwc(4)
                        + (p[128] * plg[1][1] + p[129] * plg[3][1] + p[130] * plg[5][1]) * getCswSwc(6)
                              * cos(hr * (tloc - p[131])));
      }
    }
  }

  if ((getCswSw(9) != 0.0) || (lon > -1000.0)) {

    // Longitudinal
    if (0.0 != getCswSw(10)) {
      t[10] = (1.0 + p[80] * dfa * getCswSwc(0))
              * ((p[64] * plg[2][1] + p[65] * plg[4][1] + p[66] * plg[6][1] 
                  + p[103] * plg[1][1] + p[104] * plg[3][1] + p[105] * plg[5][1]
                  + getCswSwc(4) * (p[109] * plg[1][1] + p[110] * plg[3][1] + p[111] * plg[5][1]) * cd14)
                     * cos(toRadians(lon))
                 + (p[90] * plg[2][1] + p[91] * plg[4][1] + p[92] * plg[6][1] 
                    + p[106] * plg[1][1] + p[107] * plg[3][1] + p[108] * plg[5][1]
                    + getCswSwc(4) * (p[112] * plg[1][1] + p[113] * plg[3][1] + p[114] * plg[5][1]) * cd14)
                       * sin(toRadians(lon)));
    }

    // Ut and mixed ut, longitude
    if (getCswSw(11) != 0.0) {
      t[11] = (1.0 + p[95] * plg[1][0]) * (1.0 + p[81] * dfa * getCswSwc(0))
              * (1.0 + p[119] * plg[1][0] * getCswSwc(4) * cd14)
              * ((p[68] * plg[1][0] + p[69] * plg[3][0] + p[70] * plg[5][0]) * cos(sr * (sec - p[71])));

      t[11] = t[11]
              + getCswSwc(10) * (p[76] * plg[3][2] + p[77] * plg[5][2] + p[78] * plg[7][2])
                    * cos(sr * (sec - p[79]) + 2.0 * toRadians(lon))
                    * (1.0 + p[137] * dfa * getCswSwc(0));
    }

    // Ut, longitude magnetic activity
    if (getCswSw(12) != 0.0) {

      if (sw9 != -1.0) {
        t[12] = apdf * getCswSwc(10) * (1.0 + p[120] * plg[1][0])
                    * ((p[60] * plg[2][1] + p[61] * plg[4][1] + p[62] * plg[6][1])
                       * cos(toRadians(lon - p[63])))
                + apdf * getCswSwc(10) * getCswSwc(4)
                      * (p[115] * plg[1][1] + p[116] * plg[3][1] + p[117] * plg[5][1]) * cd14
                      * cos(toRadians(lon - p[118]))
                + apdf * getCswSwc(11) * (p[83] * plg[1][0] + p[84] * plg[3][0] + p[85] * plg[5][0])
                      * cos(sr * (sec - p[75]));
      }

      else if (p[51] != 0.0) {
        t[12] = apt[0] * getCswSwc(10) * (1.0 + p[132] * plg[1][0])
                    * ((p[52] * plg[2][1] + p[98] * plg[4][1] + p[67] * plg[6][1])
                       * cos(toRadians(lon - p[97])))
                + apt[0] * getCswSwc(10) * getCswSwc(4)
                      * (p[133] * plg[1][1] + p[134] * plg[3][1] + p[135] * plg[5][1]) * cd14
                      * cos(toRadians(lon - p[136]))
                + apt[0] * getCswSwc(11)
                      * (p[55] * plg[1][0] + p[56] * plg[3][0] + p[57] * plg[5][0])
                      * cos(sr * (sec - p[58]));
      }
    }
  }

  // Parms not used: 83, 90, 100, 140-150
  greal tinf = p[30];

  for (int i = 0; i < nsw; ++i) {
    tinf = tinf + fabs(getCswSw(i)) * t[i];
  }

  return tinf;
}

//! \brief Version of globe for lower atmosphere
//!
//! \param p
//!
//! \returns tt
greal MSIS::glob7s(const greal p[])
{
  const greal dr = 1.72142e-02, pset = 2.0;

  // These values, used to speed up code, should be member variables to enable optimizations.
  // They are left as local variables to force calculations.
  greal dayl = -1.0, p14 = -1000.0, p18 = -1000.0, p32 = -1000.0, p39 = -1000.0;

  //if (p[99] == 0.0) {
  //  p[99] = pset;
  //}

  assert(pset == p[99]);
  //  cout << "Wrong parameter set for glob7s" << '\n';

  greal cd14 = 0.0, cd18 = 0.0, cd32 = 0.0, cd39 = 0.0;
  if ((day_yr != dayl) || (p[13] != p14)) {
    cd14 = cos(dr * (day_yr - p[13]));
  }
  if ((day_yr != dayl) || (p[17] != p18)) {
    cd18 = cos(2.0 * dr * (day_yr - p[17]));
  }
  if ((day_yr != dayl) || (p[31] != p32)) {
    cd32 = cos(dr * (day_yr - p[31]));
  }
  if ((day_yr != dayl) || (p[38] != p39)) {
    cd39 = cos(2.0 * dr * (day_yr - p[38]));
  }

  dayl = day_yr;
  p14 = p[13];
  p18 = p[17];
  p32 = p[31];
  p39 = p[38];

  // F10.7
  greal t[14] = {0.0};
  t[0] = p[21] * (dfa);

  // Time independent
  t[1] = p[1] * plg[2][0] + p[2] * plg[4][0] + p[22] * plg[6][0] + p[26] * plg[1][0]
         + p[14] * plg[3][0] + p[59] * plg[5][0];

  // Symmetrical annual
  t[2] = (p[18] + p[47] * plg[2][0] + p[29] * plg[4][0]) * cd32;

  // Symmetrical semiannual
  t[3] = (p[15] + p[16] * plg[2][0] + p[30] * plg[4][0]) * cd18;

  // Asymmetrical annual
  t[4] = (p[9] * plg[1][0] + p[10] * plg[3][0] + p[20] * plg[5][0]) * cd14;

  // Asymmetrical semiannual
  t[5] = (p[37] * plg[1][0]) * cd39;

  // Diurnal
  if (getCswSw(6) != 0.0) {
    greal t71 = p[11] * plg[2][1] * cd14 * getCswSwc(4);
    greal t72 = p[12] * plg[2][1] * cd14 * getCswSwc(4);
    t[6] = ((p[3] * plg[1][1] + p[4] * plg[3][1] + t71) * ctloc
            + (p[6] * plg[1][1] + p[7] * plg[3][1] + t72) * stloc);
  }

  // Semidiurnal
  if (getCswSw(7) != 0.0) {
    greal t81 = (p[23] * plg[3][2] + p[35] * plg[5][2]) * cd14 * getCswSwc(4);
    greal t82 = (p[33] * plg[3][2] + p[36] * plg[5][2]) * cd14 * getCswSwc(4);
    t[7] = ((p[5] * plg[2][2] + p[41] * plg[4][2] + t81) * c2tloc
            + (p[8] * plg[2][2] + p[42] * plg[4][2] + t82) * s2tloc);
  }

  // Terdiurnal
  if (getCswSw(13) != 0.0) {
    t[13] = p[39] * plg[3][3] * s3tloc + p[40] * plg[3][3] * c3tloc;
  }

  // Magnetic activity
  if (getCswSw(8) != 0.0) {
    if (1.0 == getCswSw(8)) {
      t[8] = apdf * (p[32] + p[45] * plg[2][0] * getCswSwc(1));
    }
    else if (getCswSw(8) == -1.0) {
      t[8] = (p[50] * apt[0] + p[96] * plg[2][0] * apt[0] * getCswSwc(1));
    }
  }

  if ((getCswSw(9) != 0.0) || (getCswSw(10) != 0.0) || (xlong > -1000.0)) {
    // Longitudinal
    t[10] = (1.0 + plg[1][0] * (p[80] * getCswSwc(4) * cos(dr * (day_yr - p[81])) 
                                + p[85] * getCswSwc(5) * cos(2.0 * dr * (day_yr - p[86])))
                 + p[83] * getCswSwc(2) * cos(dr * (day_yr - p[84])) 
                 + p[87] * getCswSwc(3) * cos(2.0 * dr * (day_yr - p[88])))
            * ((p[64] * plg[2][1] + p[65] * plg[4][1] + p[66] * plg[6][1] 
                 + p[74] * plg[1][1] + p[75] * plg[3][1] + p[76] * plg[5][1]) * cos(toRadians(xlong))
               + (p[90] * plg[2][1] + p[91] * plg[4][1] + p[92] * plg[6][1]  
                 + p[77] * plg[1][1] + p[78] * plg[3][1] + p[79] * plg[5][1]) * sin(toRadians(xlong)));
  }

  greal tt = 0.0;

  for (int i = 0; i < 14; ++i) {
    tt = tt + fabs(getCswSw(i)) * t[i];
  }

  return tt;
}

//! \brief Test if geophysical variables or switches changed and save.
//!
//! \param iyd
//! \param sec
//! \param glat
//! \param glong
//! \param stl
//! \param f107a
//! \param f107
//! \param ap
//! \param ic
//!
//! \returns 0 if unchanged and 1 if changed
int MSIS::vtst7(int iyd, greal sec, greal glat, greal glong, greal stl, greal f107a,
                   greal f107, greal ap[], int ic)
{
  int iydl[2] = {-999, -999};
  greal secl[2] = {-999.0, -999.0}, glatl[2] = {-999.0, -999.0}, gll[2] = {-999.0, -999.0},
         stll[2] = {-999.0, -999.0}, fal[2] = {-999.0, -999.0}, fl[2] = {-999.0, -999.0}, apl[7][2],
         swl[25][2], swcl[25][2];

  for (int i = 0; i < 25; i++) {
    for (int j = 0; j < 2; j++) {
      swl[i][j] = -999.0;
      swcl[i][j] = -999.0;

      if (i < 7) {
        apl[i][j] = -999.0;
      }
    }
  }

  ic = ic - 1;

  int vtst7_out = 0;

  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < 25; ++j) {
      if ((iyd != iydl[ic]) || (sec != secl[ic]) || (glat != glatl[ic]) || (glong != gll[ic])
          || (stl != stll[ic]) || (f107a != fal[ic]) || (f107 != fl[ic]) || (ap[i] != apl[i][ic])
          || (getCswSw(j) != swl[j][ic]) || (swcl[j][ic] != getCswSwc(j))) {

        vtst7_out = 1;
        iydl[ic] = iyd;
        secl[ic] = sec;
        glatl[ic] = glat;
        gll[ic] = glong;
        stll[ic] = stl;
        fal[ic] = f107a;
        fl[ic] = f107;

        for (i = 0; i < 7; ++i) {
          apl[i][ic] = ap[i];
        }

        for (i = 0; i < 25; ++i) {
          swl[i][ic] = getCswSw(i);
          swcl[i][ic] = getCswSwc(i);
        }

        return vtst7_out;
      }
    }
  }

  return vtst7_out;
}

//! \brief Unknown
//!
//! \param ex
//!
//! \returns unknown
greal MSIS::sumex(greal ex) 
{
  return (1.0 + (1.0 - pow(ex, 19)) / (1.0 - ex)*pow(ex, 0.5)); 
}

//! \brief Unknown
//!
//! \param a
//! \param p
//!
//! \returns unknown
greal MSIS::g0(greal a, const greal p[])
{
  return (a - 4.0 + (p[25] - 1.0)*(a - 4.0 + (exp(-fabs(p[24])*(a - 4.0)) - 1.0) / fabs(p[24])));
}

//! \brief Unknown
//!
//! \param ex
//! \param ap
//! \param p
//!
//! \returns unknown
greal MSIS::sg0(greal ex, const greal ap[], const greal p[])
{
  return (g0(ap[1], p)
    + (g0(ap[2], p) * ex + g0(ap[3], p) * ex * ex + g0(ap[4], p) * pow(ex, 3)
      + (g0(ap[5], p) * pow(ex, 4) + g0(ap[6], p) * pow(ex, 12)) * (1.0 - pow(ex, 8))
      / (1.0 - ex)))
    / sumex(ex);
}


// void MSIS::gtd7bk()
//{
// MSISE-00 01-Feb-02
// Tabulated MSIS data

// isdate[0] = "01-F";
// isdate[1] = "EB-0";
// isdate[2] = "2   ";
// istime[0] = "15:4";
// istime[1] = "9:27";
// name[0] = "MSIS";
// name[1] = "E-00";
//}


} // namespace GRAM
