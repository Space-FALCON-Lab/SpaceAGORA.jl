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
#include "gram.h"
#include "HWM.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
HWM::HWM()
{

}

//! \fn HWM::HWM(const HWM& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn HWM::~HWM()
//! \copydoc Atmosphere::~Atmosphere()

//! \fn HWM::getCswSw(int index)
//! \brief Get the sw value at spedified index.

//! \fn HWM::getCswSwc(int index)
//! \brief Get the swc value at spedified index.

//! \brief Horizontal wind model HWM93
//!
//! Horizontal wind model HWM93 covering all altitude regions, A.E. Hedin (1/25/93) (4/9/93)
//! Calling argument list made similar to GTs5 subroutine for MSIS-86 density model and GWs4 
//! for thermospheric winds.
void HWM::gws5(int iyd, greal sec, greal alt, greal glat, greal glon, greal stl, greal f107a, greal f107, 
  const greal ap[2], greal w[2])
{
  const greal s = 0.016, zl = 200.0;
  const int nnn = 3, mn1 = 5, mn2 = 14;
  const greal zn1[5] = { 200.0, 150.0, 130.0, 115.0, 100.0 };
  const greal zn2[14] = { 100.0, 90.0, 82.5, 75.0, 67.5, 60.0, 52.5, 45.0, 37.5, 30.0, 22.5, 15.0, 7.5, 0.0 };
  const greal sv[25] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

  int mn2s = 1, mn2m = 1;
  greal windf[2] = { 0.0 }, ww[2] = { 0.0 }, wzl[2] = { 0.0 }, wdzl[2] = { 0.0 };
  greal un1[2][mn1] = { { 0.0 } }, un2[2][mn2] = { { 0.0 } }, ugn1[2][2] = { { 0.0 } }, ugn2[2][2] = { { 0.0 } };

  tselec(sv);

  greal yrd = iyd;
  ww[0] = w[0];
  ww[1] = w[1];

  if (alt > zn1[mn1 - 1]) {

    //Exospheric wind
    glbw5e(yrd, sec, glat, glon, stl, f107a, f107, ap, pwb, pwc, windf);
    windf[0] = sw[15] * windf[0];
    windf[1] = sw[15] * windf[1];

    //Wind at z1
    glbw5m(yrd, sec, glat, glon, stl, f107a, f107, ap, pwbl, pwcl, ww);
    wzl[0] = (pwbl[0] * windf[0] + ww[0])*sw[16] * sw[17];
    wzl[1] = (pwbl[0] * windf[1] + ww[1])*sw[16] * sw[17];
    un1[0][0] = wzl[0];
    un1[1][0] = wzl[1];

    //wind derivative at z1
    ww[0] = 0.0;
    ww[1] = 0.0;
    glbw5m(yrd, sec, glat, glon, stl, f107a, f107, ap, pwbld, pwcld, ww);
    wdzl[0] = (pwbld[0] * windf[0] + ww[0])*sw[18] * sw[17];
    wdzl[1] = (pwbld[0] * windf[1] + ww[1])*sw[18] * sw[17];
    ugn1[0][0] = wdzl[0] * s;
    ugn1[1][0] = wdzl[1] * s;

    if (alt < zl) {
      //wind at zn1[1] (150)
      glbw5m(yrd, sec, glat, glon, stl, f107a, f107, ap, pb12, pc12, ww);
      un1[0][1] = (pb12[0] * windf[0] + ww[0])*sw[17];
      un1[1][1] = (pb12[0] * windf[1] + ww[1])*sw[17];

      //wind at zn1[2] (130)
      glbw5m(yrd, sec, glat, glon, stl, f107a, f107, ap, pb13, pc13, ww);
      un1[0][2] = ww[0] * sw[17];
      un1[1][2] = ww[1] * sw[17];

      //wind at zn1[3] (115)
      glbw5m(yrd, sec, glat, glon, stl, f107a, f107, ap, pb14, pc14, ww);
      un1[0][3] = ww[0] * sw[17];
      un1[1][3] = ww[1] * sw[17];
    }
  }

  // Unwinding "goto" statements led to the condition (alt <= zn1[mn1 - 1] || alt < zl)
  // This simplifies to the condition below since zn1[mn1 - 1] < zl
  if (alt < zl) {
    int mnn = max(1, min(mn2, nnn + 1));

    if (alt >= zn2[mnn - 1]) {
      //wind at zn1[4] (100)
      glbw5m(yrd, sec, glat, glon, stl, f107a, f107, ap, pb15, pc15, ww);
      un1[0][4] = ww[0] * sw[17];
      un1[1][4] = ww[1] * sw[17];

      //wind derivative at zn1[4] (100)
      glbw5m(yrd, sec, glat, glon, stl, f107a, f107, ap, pb15d, pc15d, ww);
      ugn1[0][1] = ww[0] * sw[17];
      ugn1[1][1] = ww[1] * sw[17];

      if (alt < zn1[mn1 - 1]) {
        ugn2[0][0] = ugn1[0][1];
        ugn2[1][0] = ugn1[1][1];
        un2[0][0] = un1[0][4];
        un2[1][0] = un1[1][4];
      }
    }
    else {
      ugn2[0][0] = 1.0e30;
      ugn2[1][0] = 1.0e30;
      un2[0][0] = 0.0;
      un2[1][0] = 0.0;
    }

    // Unwinding "goto" statements led to the condition (alt < zn2[mnn - 1] || alt < zn1[mn1 - 1])
    // This simplifies to the condition below since zn2[mnn - 1] < zn1[mn1 - 1]
    if (alt < zn1[mn1 - 1]) {
      int iz = 1;
      for (iz = 1; iz < mn2 + 1; iz++) {
        if (alt > zn2[iz - 1]) {
          break;
        }
      }

      mn2s = max(1, min(iz - 1, iz - nnn));
      int mn2e = min(mn2, max(mn2s + 1, iz - 1 + nnn));

      for (int i = mn2s; i < mn2e + 1; i++) {
        int ii = 2 * (i - 2) + 1;

        if (i > 1) {
          glbw5s(iyd, glat, glon, stl, &pwp[ii - 1][0], &pwp[ii][0], ww);
          un2[0][i - 1] = ww[0] * sw[19];
          un2[1][i - 1] = ww[1] * sw[19];
        }
      }

      mn2m = mn2e - mn2s + 1;
      ugn2[0][1] = 1.0e30;
      ugn2[1][1] = 1.0e30;
    }
  }

  //wind at altitude
  if (w[0] != 9898.0) {
    w[0] = wprof(alt, zl, s, windf[0], wzl[0], wdzl[0], mn1, zn1, 
                 &un1[0][0], &ugn1[0][0], mn2m, &zn2[mn2s - 1], &un2[0][mn2s - 1], &ugn2[0][0]);
  }

  if (w[1] != 9898.0) {
    w[1] = wprof(alt, zl, s, windf[1], wzl[1], wdzl[1], mn1, zn1, 
                 &un1[1][0], &ugn1[1][0], mn2m, &zn2[mn2s - 1], &un2[1][mn2s - 1], &ugn2[1][0]);
  }
}

//! \brief Calculate 2nd derivative of cubic spline interp function
void HWM::spline(const greal x[], const greal y[], int n, greal yp1, greal ypn, greal y2[])
{
  const int nmax = 100;
  greal u[nmax];

  if (yp1 > 0.99e30) {
    y2[0] = 0.0;
    u[0] = 0.0;
  }
  else {
    y2[0] = -0.5;
    u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
  }

  for (int i = 1; i < n - 1; i++) {
    greal sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
    greal p = sig * y2[i - 1] + 2.0;
    y2[i] = (sig - 1.0) / p;
    u[i] = (6.0 * ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]))
                / (x[i + 1] - x[i - 1])
            - sig * u[i - 1])
           / p;
  }

  greal qn, un;
  if (ypn > 0.99e30) {
    qn = 0.0;
    un = 0.0;
  }
  else {
    qn = 0.5;
    un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
  }

  y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);

  for (int k = n - 2; k > -1; k--) {
    y2[k] = y2[k] * y2[k + 1] + u[k];
  }
}

//! \brief Calculate cubic spline interp value adapted from numerical recipes
void HWM::splint(const greal xa[], const greal ya[], const greal y2[], int n, greal x, greal &y)
{
  int klo = 1;
  int khi = n;

  int count = 0;
  while (khi - klo > 1) {
    count++;
    int k = (khi + klo) / 2;
    if (x < xa[k - 1]) {
      khi = k;
    }
    else {
      klo = k;
    }
    if (count > 20) {
      break;
    }
  }

  greal h = xa[khi - 1] - xa[klo - 1];
  if (h == 0.0) {
    h = 0.1;
  }

  greal a = (xa[khi - 1] - x) / h;
  greal b = (x - xa[klo - 1]) / h;
  y = a * ya[klo - 1] + b * ya[khi - 1]
      + ((a * a * a - a) * y2[klo - 1] + (b * b * b - b) * y2[khi - 1]) * h * h / 6.0;
}

//! \brief Turns on and off particular variations.
//!
//! To turn on and off particular variations call tselec(sv), where
//! sv is a 25 element array containing:  0.0 for off, 1.0 for on, or
//!  2.0 for main effects off but cross terms on to get current values
//! of sw: call tretrv(sw)
void HWM::tselec(const greal sv[25])
{
  if (callTselec) {
    for (int i = 0; i < 25; i++) {
      sw[i] = fmod(sv[i], 2.0);

      if ((fabs(sv[i]) == 1.0) || (fabs(sv[i]) == 2.0)) {
        swc[i] = 1.0;
      }
      else {
        swc[i] = 0.0;
      }
    }
  }
  callTselec = false;
}

//! \brief Compute winds at altitude
greal HWM::wprof(greal z, greal zl, greal s, greal uinf, greal ulb, greal ulbd, 
  int mn1, const greal zn1[], const greal un1[], const greal ugn1[], 
  int mn2, const greal zn2[], const greal un2[], const greal ugn2[])
{
  greal xs[15], ys[15], y2out[15];

  greal wprof_out = 0.0;

  if (z >= zl) {
    greal x = s * (z - zl);
    greal f = exp(-x);

    wprof_out = uinf + (ulb - uinf) * f + (ulb - uinf + ulbd) * x * f;
    return wprof_out;
  }

  if ((z >= zn1[mn1 - 1]) && (z < zn1[0])) {
    int mn = mn1;
    greal z1 = zn1[0];
    greal z2 = zn1[mn - 1];
    greal zdif = z2 - z1;

    for (int k = 0; k < mn; ++k) {
      xs[k] = (zn1[k] - z1) / zdif;
      ys[k] = un1[k];
    }

    greal yd1 = ugn1[0] * zdif;
    greal yd2 = ugn1[1] * zdif;
    spline(xs, ys, mn, yd1, yd2, y2out);

    greal x = (z - z1) / zdif;
    splint(xs, ys, y2out, mn, x, wprof_out);
    return wprof_out;
  }

  if (z < zn2[0]) {
    int mn = mn2;
    greal z1 = zn2[0];
    greal z2 = zn2[mn - 1];
    greal zdif = z2 - z1;

    for (int k = 0; k < mn; ++k) {
      xs[k] = (zn2[k] - z1) / zdif;
      ys[k] = un2[k];
    }

    greal yd1 = ugn2[0];
    if (1.0e30 > ugn2[0]) {
      yd1 = ugn2[0] * zdif;
    }

    greal yd2 = ugn2[1];
    if (1.0e30 > ugn2[1]) {
      yd2 = ugn2[1] * zdif;
    }

    spline(xs, ys, mn, yd1, yd2, y2out);

    greal x = (z - z1) / zdif;
    splint(xs, ys, y2out, mn, x, wprof_out);
    return wprof_out;
  }

  return wprof_out;
}

//! \brief Calculate legendre polynomials plg[l+1][m+1] through order l, m
//! for cosine c and sine s of colatitude.
void HWM::legpl1(greal c, greal s, int l, int m, greal plg[][20], int lmax)
{
  assert(m <= l);
  assert(l <= lmax - 1);
  //if ((m > l) || (l > lmax - 1)) {
  //  cout << "Illegal indices to legpl1" << '\n';
  //  return;
  //}
  plg[0][0] = 1.0;

  if ((l == 0) && (m == 0)) {
    return;
  }

  //Calculate L=M case and L=M+1
  for (int mm = 0; mm < m + 1; mm++) {
    if (mm > 0) {
      plg[mm][mm] = plg[mm - 1][mm - 1] * (2.0 * mm - 1.0) * s;
    }

    if (l > mm) {
      plg[mm + 1][mm] = plg[mm][mm] * (2.0 * mm + 1) * c;
    }
  }

  if (l == 1) {
    return;
  }

  int mmx = min(m, l - 2);

  for (int mm = 0; mm < mmx + 1; mm++) {
    for (int ll = mm + 2; ll < l + 1; ll++) {
      plg[ll][mm] = ((2.0 * ll - 1.0) * c * plg[ll - 1][mm] - (ll + mm - 1.0) * plg[ll - 2][mm])
                    / (ll - mm);
    }
  }
}

//! \brief Calculates vector spherical harmonic.
//!
//! Calculates vector spherical harmonic b field theta and phi functions bt, 
//! bp through order l, m for colatitude (theta) with cosine c and sine s
//! of colatitude
void HWM::vsphr1(greal c, greal s, int l, int m, greal bt[][20], greal bp[][20], int lmax)
{
  greal plg[20][20] = { { 0.0 } };

  if ((m > l) || (l > lmax - 1)) {
    return;
  }

  bt[0][0] = 0.0;
  bp[0][0] = 0.0;

  if ((l == 0) && (m == 0)) {
    return;
  }

  legpl1(c, s, l + 1, m, plg, 20);

  greal ic = 0;
  if (fabs(s) < 1.0e-5) {
    ic = copysign(1.0, s);
    s = 0;
  }

  for (int ll = 1; ll < l + 1; ll++) {
    greal sqt = sqrt(greal(ll) * (greal(ll) + 1.0));
    int lmx = min(ll, m);

    for (int mm = 0; mm < lmx + 1; mm++) {
      if (s == 0.0) {
        if (mm != 1) {
          bt[ll][mm] = 0.0;
          bp[ll][mm] = 0.0;
        }
        else {
          bt[ll][mm] = (ll * (ll + 1.0) * (ll + 2.0) * 0.5 * pow(ic, ll + 2.0)
                        - (ll + 1.0) * c * ll * (ll + 1.0) * 0.5 * pow(ic, ll + 1.0))
                       / sqt;
          bp[ll][mm] = mm * ll * (ll + 1.0) * 0.5 * pow(ic, ll + 1.0) / sqt;

        }
      }
      else {
        bt[ll][mm] = ((ll - mm + 1.0) * plg[ll + 1][mm] - (ll + 1.0) * c * plg[ll][mm]) / (s * sqt);
        bp[ll][mm] = mm * plg[ll][mm] / (s * sqt);
      }
    }
  }
}


//! \brief glbw5s member function from HWM class
void HWM::glbw5s(int iyd, greal lat, greal lon, greal stl, const greal *pb, const greal *pc, greal ww[])
{
  const greal hr = 0.2618;
  const greal dr = 1.72142e-02;
  const greal pset = 5.0;
  const int nsw = 14;
  const int nsv = 2;
  const int ngv = 1;

  //Confirm parameter set
  assert(pb[99] == pset);

  int iyr = iyd / 1000;
  greal day = iyd - iyr * 1000.0;

  int lv = 7;
  int mv = 2;

  if ((xvl != lat) || (lv > lvl) || (mv > mvl)) {
    slat = sin(toRadians(lat));
    clat = cos(toRadians(lat));

    vsphr1(slat, clat, lv, mv, bt, bp, 20);

    xvl = lat;
    lvl = lv;
    mvl = mv;
  }

  if ((tll != stl) || (nsv > nsvl)) {
    sstl = sin(hr * stl);
    cstl = cos(hr * stl);
    s2stl = sin(2.0 * hr * stl);
    c2stl = cos(2.0 * hr * stl);
    tll = stl;
    nsvl = nsv;
  }

  greal cd14b = cos(dr * (day - pb[13]));
  greal cd14c = cos(dr * (day - pc[13]));
  greal cd18b = cos(2.0 * dr * (day - pb[17]));
  greal cd18c = cos(2.0 * dr * (day - pc[17]));
  greal cd19b = cos(2.0 * dr * (day - pb[18]));
  greal cd25b = cos(dr * (day - pb[24]));
  greal cd26b = cos(dr * (day - pb[25]));
  greal cd32c = cos(dr * (day - pc[31]));
  greal cd39c = cos(2.0 * dr * (day - pc[38]));
  greal cd64c = cos(dr * (day - pc[63]));
  greal cd87c = cos(2.0 * dr * (day - pc[86]));

  if (xll != lon || ngv > ngvl) {
    slong = sin(toRadians(lon));
    clong = cos(toRadians(lon));
    xll = lon;
    ngvl = ngv;
  }

  greal wb[2][15] = { { 0.0 } };
  greal wc[2][15] = { { 0.0 } };

  // Time independent
  constexpr greal f1b = 1.0;

  if (ww[0] != 9898.0) {
    wb[0][1] = (pb[1] * bt[2][0] + pb[2] * bt[4][0] + pb[22] * bt[6][0]) * f1b;
  }

  wb[1][1] = 0.0;
  wc[0][1] = 0.0;

  constexpr greal f1c = 1.0;

  if (ww[1] != 9898.0) {
    wc[1][1] = -(pc[1] * bt[1][0] + pc[2] * bt[3][0] + pc[22] * bt[5][0]) * f1c
               - (pc[26] * bt[2][0] + pc[14] * bt[4][0] + pc[59] * bt[6][0]) * f1c;
  }

  // symmetrical annual
  if (ww[1] != 9898.0) {
    wc[1][2] = -(pc[47] * bt[1][0] + pc[29] * bt[3][0]) * cd32c;
  }

  // symmetrical semiannual
  if (ww[0] != 9898.0) {
    wb[0][3] = (pb[16] * bt[2][0] + pb[30] * bt[4][0]) * cd18b;
  }

  wb[1][3] = 0.0;
  wc[0][3] = 0.0;

  if (ww[1] != 9898.0) {
    wc[1][3] = -(pc[16] * bt[1][0] + pc[30] * bt[3][0] + pc[49] * bt[5][0]) * cd18c;
  }


  // Asymmetrical annual
  constexpr greal f5b = 1.0;

  if (ww[0] != 9898.0) {
    wb[0][4] = (pb[9] * bt[1][0] + pb[10] * bt[3][0]) * cd14b * f5b;
  }

  wb[1][4] = 0.0;
  wc[0][4] = 0.0;

  constexpr greal f5c = 1.0;

  if (ww[1] != 9898.0) {
    wc[1][4] = -(pc[9] * bt[2][0] + pc[10] * bt[4][0] + pc[20] * bt[6][0]) * cd14c * f5c;
  }

  // Asymmetrical semiannual
  if (ww[1] != 9898.0) {
    wc[1][5] = -(pc[37] * bt[2][0] + pc[98] * bt[4][0]) * cd39c;
  }

  // Diurnal
  if (sw[6] != 0.0) {
    constexpr greal f7b = 1.0;
    constexpr greal f75b = 1.0;

    if (ww[0] != 9898.0) {
      wb[0][6] = (pb[6] * bt[1][1] + pb[7] * bt[3][1] + pb[28] * bt[5][1] + pb[88] * bt[2][1]) * sstl * f7b
                 + (pb[12] * bt[2][1]) * cd25b * sstl * f75b * swc[4]
                 + (pb[3] * bt[1][1] + pb[4] * bt[3][1] + pb[27] * bt[5][1] + pb[87] * bt[2][1]) * cstl * f7b
                 + (pb[11] * bt[2][1]) * cd25b * cstl * f75b * swc[4];
    }

    if (ww[1] != 9898.0) {
      wb[1][6] = -(pb[3] * bp[1][1] + pb[4] * bp[3][1] + pb[27] * bp[5][1] + pb[87] * bp[2][1]) * sstl * f7b
                 - (pb[11] * bp[2][1]) * cd25b * sstl * f75b * swc[4]
                 + (pb[6] * bp[1][1] + pb[7] * bp[3][1] + pb[28] * bp[5][1] + pb[88] * bp[2][1]) * cstl * f7b
                 + (pb[12] * bp[2][1]) * cd25b * cstl * f75b * swc[4];
    }

    constexpr greal f7c = 1.0;
    constexpr greal f75c = 1.0;

    if (ww[0] != 9898.0) {
      wc[0][6] = -(pc[3] * bp[2][1] + pc[4] * bp[4][1] + pc[27] * bp[6][1] + pc[87] * bp[1][1]) * sstl * f7c
                 - (pc[11] * bp[1][1]) * cd25b * sstl * f75c * swc[4]
                 + (pc[6] * bp[2][1] + pc[7] * bp[4][1] + pc[28] * bp[6][1] + pc[88] * bp[1][1]) * cstl * f7c
                 + (pc[12] * bp[1][1]) * cd25b * cstl * f75c * swc[4];
    }

    if (ww[1] != 9898.0) {
      wc[1][6] = -(pc[6] * bt[2][1] + pc[7] * bt[4][1] + pc[28] * bt[6][1] + pc[88] * bt[1][1]) * sstl * f7c
                 - (pc[12] * bt[1][1]) * cd25b * sstl * f75c * swc[4]
                 - (pc[3] * bt[2][1] + pc[4] * bt[4][1] + pc[27] * bt[6][1] + pc[87] * bt[1][1]) * cstl * f7c
                 - (pc[11] * bt[1][1]) * cd25b * cstl * f75c * swc[4];
    }
  }

  // semidiurnal
  if (sw[7] != 0.0) {
    constexpr greal f8b = 1.0;

    if (ww[0] != 9898.0) {
      wb[0][7] = (pb[8] * bt[2][2] + pb[42] * bt[4][2] + pb[34] * bt[6][2] + pb[97] * bt[3][2]
                  + (pb[33] * bt[3][2]) * cd26b * swc[4] + (pb[36] * bt[3][2]) * cd19b * swc[5])
                     * s2stl * f8b
                 + (pb[5] * bt[2][2] + pb[41] * bt[4][2] + pb[32] * bt[6][2] + pb[95] * bt[3][2]
                    + (pb[23] * bt[3][2]) * cd26b * swc[4] + (pb[35] * bt[3][2]) * cd19b * swc[5])
                       * c2stl * f8b;
    }

    if (ww[1] != 9898.0) {
      wb[1][7] = -(pb[5] * bp[2][2] + pb[41] * bp[4][2] + pb[32] * bp[6][2] + pb[95] * bp[3][2]
                   + (pb[23] * bp[3][2]) * cd26b * swc[4] + (pb[35] * bp[3][2]) * cd19b * swc[5])
                     * s2stl * f8b
                 + (pb[8] * bp[2][2] + pb[42] * bp[4][2] + pb[34] * bp[6][2] + pb[97] * bp[3][2]
                    + (pb[33] * bp[3][2]) * cd26b * swc[4] + (pb[36] * bp[3][2]) * cd19b * swc[5])
                       * c2stl * f8b;
    }

    constexpr greal f8c = 1.0;

    if (ww[0] != 9898.0) {
      wc[0][7] = -(pc[5] * bp[3][2] + pc[41] * bp[5][2] + pc[32] * bp[7][2] + pc[95] * bp[2][2]
                   + (pc[23] * bp[2][2]) * cd26b * swc[4] + (pc[35] * bp[2][2]) * cd19b * swc[5])
                     * s2stl * f8c
                 + (pc[8] * bp[3][2] + pc[42] * bp[5][2] + pc[34] * bp[7][2] + pc[97] * bp[2][2]
                    + (pc[33] * bp[2][2]) * cd26b * swc[4] + (pc[36] * bp[2][2]) * cd19b * swc[5])
                       * c2stl * f8c;
    }

    if (ww[1] != 9898.0) {
      wc[1][7] = -(pc[8] * bt[3][2] + pc[42] * bt[5][2] + pc[34] * bt[7][2] + pc[97] * bt[2][2]
                   + (pc[33] * bt[2][2]) * cd26b * swc[4] + (pc[36] * bt[2][2]) * cd19b * swc[5])
                     * s2stl * f8c
                 - (pc[5] * bt[3][2] + pc[41] * bt[5][2] + pc[32] * bt[7][2] + pc[95] * bt[2][2]
                    + (pc[23] * bt[2][2]) * cd26b * swc[4] + (pc[35] * bt[2][2]) * cd19b * swc[5])
                       * c2stl * f8c;
    }
  }

  // Longitudinal
  if ((sw[9] != 0.0) || (sw[10] != 0.0)) {

    if (ww[0] != 9898.0) {
      wc[0][10] = -(pc[64] * bp[1][1] + pc[65] * bp[3][1] + pc[66] * bp[5][1] + pc[74] * bp[2][1]
                    + pc[75] * bp[4][1] + pc[76] * bp[6][1]
                    + (pc[56] * bp[1][1] + pc[58] * bp[3][1] + pc[61] * bp[5][1] + pc[50] * bp[2][1]
                       + pc[52] * bp[4][1] + pc[54] * bp[6][1])
                          * cd64c * swc[2]
                    + (pc[73] * bp[1][1] + pc[81] * bp[3][1] + pc[84] * bp[5][1] + pc[67] * bp[2][1]
                       + pc[69] * bp[4][1] + pc[71] * bp[6][1])
                          * cd87c * swc[3])
                      * slong
                  + (pc[90] * bp[1][1] + pc[91] * bp[3][1] + pc[92] * bp[5][1] + pc[77] * bp[2][1]
                     + pc[78] * bp[4][1] + pc[79] * bp[6][1]
                     + (pc[57] * bp[1][1] + pc[60] * bp[3][1] + pc[62] * bp[5][1]
                        + pc[51] * bp[2][1] + pc[53] * bp[4][1] + pc[55] * bp[6][1])
                           * cd64c * swc[2]
                     + (pc[80] * bp[1][1] + pc[83] * bp[3][1] + pc[85] * bp[5][1]
                        + pc[68] * bp[2][1] + pc[70] * bp[4][1] + pc[72] * bp[6][1])
                           * cd87c * swc[3])
                        * clong;
    }

    if (ww[1] != 9898.0) {
      wc[1][10] = -(pc[90] * bt[1][1] + pc[91] * bt[3][1] + pc[92] * bt[5][1] + pc[77] * bt[2][1]
                    + pc[78] * bt[4][1] + pc[79] * bt[6][1]
                    + (pc[57] * bt[1][1] + pc[60] * bt[3][1] + pc[62] * bt[5][1] + pc[51] * bt[2][1]
                       + pc[53] * bt[4][1] + pc[55] * bt[6][1])
                          * cd64c * swc[2]
                    + (pc[80] * bt[1][1] + pc[83] * bt[3][1] + pc[85] * bt[5][1] + pc[68] * bt[2][1]
                       + pc[70] * bt[4][1] + pc[72] * bt[6][1])
                          * cd87c * swc[3])
                      * slong
                  - (pc[64] * bt[1][1] + pc[65] * bt[3][1] + pc[66] * bt[5][1] + pc[74] * bt[2][1]
                     + pc[75] * bt[4][1] + pc[76] * bt[6][1]
                     + (pc[56] * bt[1][1] + pc[58] * bt[3][1] + pc[61] * bt[5][1]
                        + pc[50] * bt[2][1] + pc[52] * bt[4][1] + pc[54] * bt[6][1])
                           * cd64c * swc[2]
                     + (pc[73] * bt[1][1] + pc[81] * bt[3][1] + pc[84] * bt[5][1]
                        + pc[67] * bt[2][1] + pc[69] * bt[4][1] + pc[71] * bt[6][1])
                           * cd87c * swc[3])
                        * clong;
    }
  }

  greal wbt[2] = { 0.0, 0.0 };
  greal wct[2] = { 0.0, 0.0 };

  //sum winds and change meridional sign to + North
  for (int k = 0; k < nsw; ++k) {
    wbt[0] = wbt[0] - fabs(sw[k]) * wb[0][k];
    wct[0] = wct[0] - fabs(sw[k]) * wc[0][k];
    wbt[1] = wbt[1] + fabs(sw[k]) * wb[1][k];
    wct[1] = wct[1] + fabs(sw[k]) * wc[1][k];
  }

  if (ww[0] != 9898.0) {
    ww[0] = wbt[0] * sw[23] + wct[0] * sw[24];
  }

  if (ww[1] != 9898.0) {
    ww[1] = wbt[1] * sw[23] + wct[1] * sw[24];
  }
}

//! \brief glbw5e member function from HWM class
void HWM::glbw5e(greal yrd, greal sec, greal lat, greal lon, greal stl, greal f107a, greal f107, 
  const greal ap[], const greal pb[], const greal pc[], greal ww[])
{
  const greal sr = 7.2722e-05, hr = 0.2618, dr = 1.72142e-02, pset = 3.0;
  const int nsw = 14, lv = 12, mv = 3, nsv = 3;

  greal wb[2][15] = { { 0.0 } }, wc[2][15] = { { 0.0 } };

  //Confirm parameter set
  assert(pb[99] == pset);

  for (int j = 0; j < nsw; j++) {
    wb[0][j] = 0.0;
    wb[1][j] = 0.0;
    wc[0][j] = 0.0;
    wc[1][j] = 0.0;
  }

  greal sw9 = 1.0;
  if (sw[8] > 0.0) {
    sw9 = 1.0;
  }
  if (sw[8] < 0.0) {
    sw9 = -1.0;
  }

  int iyr = int(yrd) / 1000;
  greal day = yrd - iyr * 1000.0;

  if ((xvl != lat) || (lv > lvl) || (mv > mvl)) {
    slat = sin(toRadians(lat));
    clat = cos(toRadians(lat));

    vsphr1(slat, clat, lv, mv, bt, bp, 20);

    xvl = lat;
    lvl = lv;
    mvl = mv;
  }

  if ((tll != stl) || (nsv > nsvl)) {
    sstl = sin(hr * stl);
    cstl = cos(hr * stl);
    s2stl = sin(2.0 * hr * stl);
    c2stl = cos(2.0 * hr * stl);
    s3stl = sin(3.0 * hr * stl);
    c3stl = cos(3.0 * hr * stl);
    tll = stl;
    nsvl = nsv;
  }

  greal cd14 = cos(dr * (day - pb[13]));
  greal cd18 = cos(2.0 * dr * (day - pb[17]));

  if (xll != lon) {
    slong = sin(toRadians(lon));
    clong = cos(toRadians(lon));
    s2long = sin(2.0*toRadians(lon));
    c2long = cos(2.0*toRadians(lon));
    xll = lon;
    ngvl = 2;
  }

  // f10.7 effect
  greal df = f107 - f107a;
  greal dfa = f107a - 150.0;
  greal dfc = dfa + pb[19] * df;

  // Time independent
  greal f1b = 1.0 + pb[21] * dfc * swc[0];
  if (ww[0] != 9898.0) {
    wb[0][1] = (pb[1] * bt[2][0] + pb[2] * bt[4][0] + pb[22] * bt[6][0]) * f1b;
  }

  wb[1][1] = 0.0;
  wc[0][1] = 0.0;

  greal f1c = 1.0 + pc[21] * dfc * swc[0];

  if (ww[1] != 9898.0) {
    wc[1][1] = -(pc[1] * bt[1][0] + pc[2] * bt[3][0] + pc[22] * bt[5][0]) * f1c
               - (pc[26] * bt[2][0] + pc[14] * bt[4][0] + pc[59] * bt[6][0] + pc[160] * bt[8][0]
                  + pc[161] * bt[10][0] + pc[162] * bt[12][0]) * f1c;
  }

  // symmetrical annual
  // symmetrical semiannual
  if (ww[0] != 9898.0) {
    wb[0][3] = (pb[16] * bt[2][0] + pb[30] * bt[4][0]) * cd18;
  }

  wb[1][3] = 0.0;
  wc[0][3] = 0.0;

  if (ww[1] != 9898.0) {
    wc[1][3] = -(pc[16] * bt[1][0] + pc[30] * bt[3][0]) * cd18;
  }

  // Asymmetrical annual
  greal f5b = 1.0 + pb[47] * dfc*swc[0];

  if (ww[0] != 9898.0) {
    wb[0][4] = (pb[9] * bt[1][0] + pb[10] * bt[3][0]) * cd14 * f5b;
  }

  wb[1][4] = 0.0;
  wc[0][4] = 0.0;

  greal f5c = 1.0 + pc[47] * dfc*swc[0];

  if (ww[1] != 9898.0) {
    wc[1][4] = -(pc[9] * bt[2][0] + pc[10] * bt[4][0]) * cd14 * f5c;
  }

  // Asymmetrical semiannual
  //     none
  // Diurnal
  if (sw[6] != 0.0) {
    greal f7b = 1.0 + pb[49] * dfc*swc[0];
    greal f75b = 1.0 + pb[82] * dfc*swc[0];

    if (ww[0] != 9898.0) {
      wb[0][6] = (pb[6] * bt[1][1] + pb[7] * bt[3][1] + pb[28] * bt[5][1] + pb[141] * bt[7][1]
                  + pb[143] * bt[9][1] + pb[181] * bt[2][1] + pb[183] * bt[4][1])
                     * sstl * f7b
                 + (pb[12] * bt[2][1] + pb[145] * bt[4][1]) * cd14 * sstl * f75b * swc[4]
                 + (pb[170] * bt[1][1] + pb[172] * bt[3][1]) * cd18 * sstl * f75b * swc[3]
                 + (pb[3] * bt[1][1] + pb[4] * bt[3][1] + pb[27] * bt[5][1] + pb[140] * bt[7][1]
                    + pb[142] * bt[9][1] + pb[180] * bt[2][1] + pb[182] * bt[4][1])
                       * cstl * f7b
                 + (pb[11] * bt[2][1] + pb[144] * bt[4][1]) * cd14 * cstl * f75b * swc[4]
                 + (pb[169] * bt[1][1] + pb[171] * bt[3][1]) * cd18 * cstl * f75b * swc[3];
    }

    if (ww[1] != 9898.0) {
      wb[1][6] = -(pb[3] * bp[1][1] + pb[4] * bp[3][1] + pb[27] * bp[5][1] + pb[140] * bp[7][1]
                   + pb[142] * bp[9][1] + pb[180] * bp[2][1] + pb[182] * bp[4][1])
                     * sstl * f7b
                 - (pb[11] * bp[2][1] + pb[144] * bp[4][1]) * cd14 * sstl * f75b * swc[4]
                 - (pb[169] * bp[1][1] + pb[171] * bp[3][1]) * cd18 * sstl * f75b * swc[3]
                 + (pb[6] * bp[1][1] + pb[7] * bp[3][1] + pb[28] * bp[5][1] + pb[141] * bp[7][1]
                    + pb[143] * bp[9][1] + pb[181] * bp[2][1] + pb[183] * bp[4][1])
                       * cstl * f7b
                 + (pb[12] * bp[2][1] + pb[145] * bp[4][1]) * cd14 * cstl * f75b * swc[4]
                 + (pb[170] * bp[1][1] + pb[172] * bp[3][1]) * cd18 * cstl * f75b * swc[3];
    }

    greal f7c = 1.0 + pc[49] * dfc * swc[0];
    greal f75c = 1.0 + pc[82] * dfc * swc[0];

    if (ww[0] != 9898.0) {
      wc[0][6] = -(pc[3] * bp[2][1] + pc[4] * bp[4][1] + pc[27] * bp[6][1] + pc[140] * bp[8][1]
                   + pc[142] * bp[10][1] + pc[180] * bp[1][1] + pc[182] * bp[3][1]
                   + pc[184] * bp[5][1] + pc[186] * bp[7][1] + pc[188] * bp[9][1])
                     * sstl * f7c
                 - (pc[11] * bp[1][1] + pc[144] * bp[3][1]) * cd14 * sstl * f75c * swc[4]
                 - (pc[169] * bp[2][1] + pc[171] * bp[4][1]) * cd18 * sstl * f75c * swc[3]
                 + (pc[6] * bp[2][1] + pc[7] * bp[4][1] + pc[28] * bp[6][1] + pc[141] * bp[8][1]
                    + pc[143] * bp[10][1] + pc[181] * bp[1][1] + pc[183] * bp[3][1]
                    + pc[185] * bp[5][1] + pc[187] * bp[7][1] + pc[189] * bp[9][1])
                       * cstl * f7c
                 + (pc[12] * bp[1][1] + pc[145] * bp[3][1]) * cd14 * cstl * f75c * swc[4]
                 + (pc[170] * bp[2][1] + pc[172] * bp[4][1]) * cd18 * cstl * f75c * swc[3];
    }

    if (ww[1] != 9898.0) {
      wc[1][6] = -(pc[6] * bt[2][1] + pc[7] * bt[4][1] + pc[28] * bt[6][1] + pc[141] * bt[8][1]
                   + pc[143] * bt[10][1] + pc[181] * bt[1][1] + pc[183] * bt[3][1]
                   + pc[185] * bt[5][1] + pc[187] * bt[7][1] + pc[189] * bt[9][1])
                     * sstl * f7c
                 - (pc[12] * bt[1][1] + pc[145] * bt[3][1]) * cd14 * sstl * f75c * swc[4]
                 - (pc[170] * bt[2][1] + pc[172] * bt[4][1]) * cd18 * sstl * f75c * swc[3]
                 - (pc[3] * bt[2][1] + pc[4] * bt[4][1] + pc[27] * bt[6][1] + pc[140] * bt[8][1]
                    + pc[142] * bt[10][1] + pc[180] * bt[1][1] + pc[182] * bt[3][1]
                    + pc[184] * bt[5][1] + pc[186] * bt[7][1] + pc[188] * bt[9][1])
                       * cstl * f7c
                 - (pc[11] * bt[1][1] + pc[144] * bt[3][1]) * cd14 * cstl * f75c * swc[4]
                 - (pc[169] * bt[2][1] + pc[171] * bt[4][1]) * cd18 * cstl * f75c * swc[3];
    }
  }

  // semidiurnal
  if (sw[7] != 0.0) {
    greal f8b = 1.0 + pb[89] * dfc * swc[0];

    if (ww[0] != 9898.0) {
      wb[0][7] = (pb[8] * bt[2][2] + pb[42] * bt[4][2] + pb[110] * bt[6][2]
                  + (pb[33] * bt[3][2] + pb[147] * bt[5][2]) * cd14 * swc[4]
                  + (pb[133] * bt[2][2]) * cd18 * swc[3] + pb[151] * bt[3][2] + pb[153] * bt[5][2]
                  + pb[155] * bt[7][2] + pb[157] * bt[9][2])
                     * s2stl * f8b
                 + (pb[5] * bt[2][2] + pb[41] * bt[4][2] + pb[109] * bt[6][2]
                    + (pb[23] * bt[3][2] + pb[146] * bt[5][2]) * cd14 * swc[4]
                    + (pb[134] * bt[2][2]) * cd18 * swc[3] + pb[150] * bt[3][2] + pb[152] * bt[5][2]
                    + pb[154] * bt[7][2] + pb[156] * bt[9][2])
                       * c2stl * f8b;
    }

    if (ww[1] != 9898.0) {
      wb[1][7] = -(pb[5] * bp[2][2] + pb[41] * bp[4][2] + pb[109] * bp[6][2]
                   + (pb[23] * bp[3][2] + pb[146] * bp[5][2]) * cd14 * swc[4]
                   + (pb[134] * bp[2][2]) * cd18 * swc[3] + pb[150] * bp[3][2] + pb[152] * bp[5][2]
                   + pb[154] * bp[7][2] + pb[156] * bp[9][2])
                     * s2stl * f8b
                 + (pb[8] * bp[2][2] + pb[42] * bp[4][2] + pb[110] * bp[6][2]
                    + (pb[33] * bp[3][2] + pb[147] * bp[5][2]) * cd14 * swc[4]
                    + (pb[133] * bp[2][2]) * cd18 * swc[3] + pb[151] * bp[3][2] + pb[153] * bp[5][2]
                    + pb[155] * bp[7][2] + pb[157] * bp[9][2])
                       * c2stl * f8b;
    }

    greal f8c = 1.0 + pc[89] * dfc * swc[0];

    if (ww[0] != 9898.0) {
      wc[0][7] = -(pc[5] * bp[3][2] + pc[41] * bp[5][2] + pc[109] * bp[7][2]
                   + (pc[23] * bp[2][2] + pc[146] * bp[4][2]) * cd14 * swc[4]
                   + (pc[134] * bp[3][2]) * cd18 * swc[3] + pc[150] * bp[2][2] + pc[152] * bp[4][2]
                   + pc[154] * bp[6][2] + pc[156] * bp[8][2])
                     * s2stl * f8c
                 + (pc[8] * bp[3][2] + pc[42] * bp[5][2] + pc[110] * bp[7][2]
                    + (pc[33] * bp[2][2] + pc[147] * bp[4][2]) * cd14 * swc[4]
                    + (pc[133] * bp[3][2]) * cd18 * swc[3] + pc[151] * bp[2][2] + pc[153] * bp[4][2]
                    + pc[155] * bp[6][2] + pc[157] * bp[8][2])
                       * c2stl * f8c;
    }

    if (ww[1] != 9898.0) {
      wc[1][7] = -(pc[8] * bt[3][2] + pc[42] * bt[5][2] + pc[110] * bt[7][2]
                   + (pc[33] * bt[2][2] + pc[147] * bt[4][2]) * cd14 * swc[4]
                   + (pc[133] * bt[3][2]) * cd18 * swc[3] + pc[151] * bt[2][2] + pc[153] * bt[4][2]
                   + pc[155] * bt[6][2] + pc[157] * bt[8][2])
                     * s2stl * f8c
                 - (pc[5] * bt[3][2] + pc[41] * bt[5][2] + pc[109] * bt[7][2]
                    + (pc[23] * bt[2][2] + pc[146] * bt[4][2]) * cd14 * swc[4]
                    + (pc[134] * bt[3][2]) * cd18 * swc[3] + pc[150] * bt[2][2] + pc[152] * bt[4][2]
                    + pc[154] * bt[6][2] + pc[156] * bt[8][2])
                       * c2stl * f8c;
    }
  }

  // Terdiurnal
  if (sw[13] != 0.0) {
    greal f14b = 1.0;

    if (ww[0] != 9898.0) {
      wb[0][13] = (pb[39] * bt[3][3] + pb[148] * bt[5][3] + pb[113] * bt[7][3]
                   + (pb[93] * bt[4][3] + pb[46] * bt[6][3]) * cd14 * swc[4])
                      * s3stl * f14b
                  + (pb[40] * bt[3][3] + pb[149] * bt[5][3] + pb[114] * bt[7][3]
                     + (pb[94] * bt[4][3] + pb[48] * bt[6][3]) * cd14 * swc[4])
                        * c3stl * f14b;
    }

    if (ww[1] != 9898.0) {
      wb[1][13] = -(pb[40] * bp[3][3] + pb[149] * bp[5][3] + pb[114] * bp[7][3]
                    + (pb[94] * bp[4][3] + pb[48] * bp[6][3]) * cd14 * swc[4])
                      * s3stl * f14b
                  + (pb[39] * bp[3][3] + pb[148] * bp[5][3] + pb[113] * bp[7][3]
                     + (pb[93] * bp[4][3] + pb[46] * bp[6][3]) * cd14 * swc[4])
                        * c3stl * f14b;
    }

    greal f14c = 1.0;

    if (ww[0] != 9898.0) {
      wc[0][13] = -(pc[40] * bp[4][3] + pc[149] * bp[6][3] + pc[114] * bp[8][3]
                    + (pc[94] * bp[3][3] + pc[48] * bp[5][3]) * cd14 * swc[4])
                      * s3stl * f14c
                  + (pc[39] * bp[4][3] + pc[148] * bp[6][3] + pc[113] * bp[8][3]
                     + (pc[93] * bp[3][3] + pc[46] * bp[5][3]) * cd14 * swc[4])
                        * c3stl * f14c;
    }

    if (ww[1] != 9898.0) {
      wc[1][13] = -(pc[39] * bt[4][3] + pc[148] * bt[6][3] + pc[113] * bt[8][3]
                    + (pc[93] * bt[3][3] + pc[46] * bt[5][3]) * cd14 * swc[4])
                      * s3stl * f14c
                  - (pc[40] * bt[4][3] + pc[149] * bt[6][3] + pc[114] * bt[8][3]
                     + (pc[94] * bt[3][3] + pc[48] * bt[5][3]) * cd14 * swc[4])
                        * c3stl * f14c;
    }
  }

  // Magnetic activity
  greal apd = 0.0, apdf = 0.0, apdfc = 0.0, apt = 0.0;
  if (sw[8] != 0.0) {
    if (sw9 != -1.0) {

      // Daily ap
      apd = ap[0] - 4.0;
      apdf = (apd + (pb[44] - 1.0) * (apd + (exp(-pb[43] * apd) - 1.0) / pb[43]));
      // apdfc=(apd + (pc[44] - 1.0)*(apd + (exp( - pc[43]*apd) - 1.0)/pc[43]));
      apdfc = apdf;

      if (apd != 0) {

        if (ww[0] != 9898.0) {
          wb[0][8] = (pb[45] * bt[2][0] + pb[34] * bt[4][0] + pb[32] * bt[6][0]) * apdf
                     + (pb[174] * bt[2][2] + pb[176] * bt[4][2]) * s2stl * apdf
                     + (pb[173] * bt[2][2] + pb[175] * bt[4][2]) * c2stl * apdf;
        }

        if (ww[1] != 9898.0) {
          wb[1][8] = -(pb[173] * bp[2][2] + pb[175] * bp[4][2]) * s2stl * apdf
                     + (pb[174] * bp[2][2] + pb[176] * bp[4][2]) * c2stl * apdf;
        }

        if (ww[0] != 9898.0) {
          wc[0][8] = swc[6] * wc[0][6] * pc[121] * apdfc
                     - (pc[173] * bp[3][2] + pc[175] * bp[5][2]) * s2stl * apdfc
                     + (pc[174] * bp[3][2] + pc[176] * bp[5][2]) * c2stl * apdfc;
        }

        wc[1][8] = -(pc[45] * bt[1][0] + pc[34] * bt[3][0] + pc[32] * bt[5][0]) * apdfc
                   + swc[6] * wc[1][6] * pc[121] * apdfc
                   - (pc[174] * bt[3][2] + pc[176] * bt[5][2]) * s2stl * apdfc
                   - (pc[173] * bt[3][2] + pc[175] * bt[5][2]) * c2stl * apdfc;
      }
    }				
    if (ww[1] != 9898.0) {

    }
    else {
      //if (pb[24] < 1.0e-04) {
      //	pb[24] = 1.0e-04;
      //}

      // apt = g0(ap[1], pb);
      apt = 0.0;
      if (apt != 0.0) {

        if (ww[0] != 9898.0) {
          wb[0][8] = (pb[96] * bt[2][0] + pb[54] * bt[4][0] + pb[50] * bt[6][0]) * apt
                     + (pb[159] * bt[2][2] + pb[178] * bt[4][2]) * s2stl * apt
                     + (pb[158] * bt[2][2] + pb[177] * bt[4][2]) * c2stl * apt;
        }

        if (ww[1] != 9898.0) {
          wb[1][8] = -(pb[158] * bp[2][2] + pb[177] * bp[4][2]) * s2stl * apt
                     + (pb[159] * bp[2][2] + pb[178] * bp[4][2]) * c2stl * apt;
        }

        if (ww[0] != 9898.0) {
          wc[0][8] = swc[6] * wc[0][6] * pc[128] * apt
                     - (pc[158] * bp[3][2] + pc[177] * bp[5][2]) * s2stl * apt
                     + (pc[159] * bp[3][2] + pc[178] * bp[5][2]) * c2stl * apt;
        }

        if (ww[1] != 9898.0) {
          wc[1][8] = -(pc[96] * bt[1][0] + pc[54] * bt[3][0] + pc[50] * bt[5][0]) * apt
                     + swc[6] * wc[1][6] * pc[128] * apt
                     - (pc[159] * bt[3][2] + pc[178] * bt[5][2]) * s2stl * apt
                     - (pc[158] * bt[3][2] + pc[177] * bt[5][2]) * c2stl * apt;
        }
      }
    }
  }

  if (sw[9] != 0.0) {

    // Longitudinal
    greal dbasy1 = 1.0 + pb[198] * slat;
    greal dbasy2 = 1.0 + pb[199] * slat;
    greal f11b = 1.0 + pb[80] * dfc * swc[0];

    if (sw[10] != 0.0) {

      if (ww[0] != 9898.0) {
        wb[0][10]
            = (pb[90] * bt[2][1] + pb[91] * bt[4][1] + pb[92] * bt[6][1]) * slong * dbasy1 * f11b
              + (pb[64] * bt[2][1] + pb[65] * bt[4][1] + pb[66] * bt[6][1]) * clong * dbasy1 * f11b
              + (pb[190] * bt[2][2] + pb[192] * bt[4][2] + pb[194] * bt[6][2] + pb[196] * bt[8][2])
                    * s2long * dbasy2 * f11b
              + (pb[191] * bt[2][2] + pb[193] * bt[4][2] + pb[195] * bt[6][2] + pb[197] * bt[8][2])
                    * c2long * dbasy2 * f11b;
      }

      if (ww[1] != 9898.0) {
        wb[1][10]
            = -(pb[64] * bp[2][1] + pb[65] * bp[4][1] + pb[66] * bp[6][1]) * slong * dbasy1 * f11b
              + (pb[90] * bp[2][1] + pb[91] * bp[4][1] + pb[92] * bp[6][1]) * clong * dbasy1 * f11b
              - (pb[191] * bp[2][2] + pb[193] * bp[4][2] + pb[195] * bp[6][2] + pb[197] * bp[8][2])
                    * s2long * dbasy2 * f11b
              + (pb[190] * bp[2][2] + pb[192] * bp[4][2] + pb[194] * bp[6][2] + pb[196] * bp[8][2])
                    * c2long * dbasy2 * f11b;
      }

      greal dcasy1 = 1.0 + pc[198] * slat;
      greal dcasy2 = 1.0 + pc[199] * slat;
      greal f11c = 1.0 + pc[80] * dfc * swc[0];

      if (ww[0] != 9898.0) {
        wc[0][10]
            = -(pc[64] * bp[1][1] + pc[65] * bp[3][1] + pc[66] * bp[5][1] + pc[72] * bp[7][1]
                + pc[73] * bp[9][1])
                  * slong * dcasy1 * f11c
              + (pc[90] * bp[1][1] + pc[91] * bp[3][1] + pc[92] * bp[5][1] + pc[86] * bp[7][1]
                 + pc[87] * bp[9][1])
                    * clong * dcasy1 * f11c
              - (pc[191] * bp[3][2] + pc[193] * bp[5][2] + pc[195] * bp[7][2] + pc[197] * bp[9][2])
                    * s2long * dcasy2 * f11c
              + (pc[190] * bp[3][2] + pc[192] * bp[5][2] + pc[194] * bp[7][2] + pc[196] * bp[9][2])
                    * c2long * dcasy2 * f11c;
      }

      if (ww[1] != 9898.0) {
        wc[1][10]
            = -(pc[90] * bt[1][1] + pc[91] * bt[3][1] + pc[92] * bt[5][1] + pc[86] * bt[7][1]
                + pc[87] * bt[9][1])
                  * slong * dcasy1 * f11c
              - (pc[64] * bt[1][1] + pc[65] * bt[3][1] + pc[66] * bt[5][1] + pc[72] * bt[7][1]
                 + pc[73] * bt[9][1])
                    * clong * dcasy1 * f11c
              - (pc[190] * bt[3][2] + pc[192] * bt[5][2] + pc[194] * bt[7][2] + pc[196] * bt[9][2])
                    * s2long * dcasy2 * f11c
              - (pc[191] * bt[3][2] + pc[193] * bt[5][2] + pc[195] * bt[7][2] + pc[197] * bt[9][2])
                    * c2long * dcasy2 * f11c;
      }
    }

    // Ut & mixed ut/lon
    greal utbasy = 1.0;
    greal f12b = 1.0 + pb[81] * dfc * swc[0];

    if (sw[11] != 0.0) {

      if (ww[0] != 9898.0) {
        wb[0][11] = (pb[68] * bt[1][0] + pb[69] * bt[3][0] + pb[70] * bt[5][0] + pb[115] * bt[7][0]
                     + pb[116] * bt[9][0] + pb[117] * bt[11][0])
                        * cos(sr * (sec - pb[71])) * utbasy * f12b
                    + (pb[76] * bt[3][2] + pb[77] * bt[5][2] + pb[78] * bt[7][2])
                          * cos(sr * (sec - pb[79]) + 2.0 * toRadians(lon)) * utbasy * f12b
                          * swc[10];
      }

      if (ww[1] != 9898.0) {
        wb[1][11] = (pb[76] * bp[3][2] + pb[77] * bp[5][2] + pb[78] * bp[7][2])
                    * cos(sr * (sec - pb[79] + 21600.0) + 2.0 * toRadians(lon)) * utbasy * f12b
                    * swc[10];
      }

      greal utcasy = 1.0;
      greal f12c = 1.0 + pc[81] * dfc * swc[0];

      if (ww[0] != 9898.0) {
        wc[0][11] = (pc[76] * bp[2][2] + pc[77] * bp[4][2] + pc[78] * bp[6][2] + pc[164] * bp[8][2]
                     + pc[165] * bp[10][2] + pc[166] * bp[12][2])
                    * cos(sr * (sec - pc[79]) + 2.0 * toRadians(lon)) * utcasy * f12c * swc[10];
      }

      if (ww[1] != 9898.0) {
        wc[1][11] = -(pc[68] * bt[2][0] + pc[69] * bt[4][0] + pc[70] * bt[6][0] + pc[115] * bt[8][0]
                      + pc[116] * bt[10][0] + pc[117] * bt[12][0])
                        * cos(sr * (sec - pc[71])) * utcasy * f12c
                    + (pc[76] * bt[2][2] + pc[77] * bt[4][2] + pc[78] * bt[6][2]
                       + pc[164] * bt[8][2] + pc[165] * bt[10][2] + pc[166] * bt[12][2])
                          * cos(sr * (sec - pc[79] + 21600.0) + 2.0 * toRadians(lon)) * utcasy
                          * f12c * swc[10];
      }
    }

    // Mixed lon, ut, ap
    if ((sw[12] != 0.0) || (apd != 0.0)) {

      if (sw9 != -1.0) {

        if (ww[0] != 9898.0) {
          wb[0][12] = (pb[60] * bt[2][1] + pb[61] * bt[4][1] + pb[62] * bt[6][1])
                          * cos(toRadians(lon - pb[63])) * apdf * swc[10]
                      + (pb[83] * bt[1][0] + pb[84] * bt[3][0] + pb[85] * bt[5][0])
                            * cos(sr * (sec - pb[75])) * apdf * swc[11];
        }

        if (ww[1] != 9898.0) {
          wb[1][12] = (pb[60] * bp[2][1] + pb[61] * bp[4][1] + pb[62] * bp[6][1])
                      * cos(toRadians(lon - pb[63] + 90.)) * apdf * swc[10];
        }

        if (ww[0] != 9898.0) {
          wc[0][12] = swc[10] * wc[0][10] * pc[60] * apdfc + swc[11] * wc[0][11] * pc[83] * apdfc;
        }

        if (ww[1] != 9898.0) {
          wc[1][12] = swc[10] * wc[1][10] * pc[60] * apdfc + swc[11] * wc[1][11] * pc[83] * apdfc;
        }
      }
      else {

        if (apt != 0.0) {

          if (ww[0] != 9898.0) {
            wb[0][12] = (pb[52] * bt[2][1] + pb[98] * bt[4][1] + pb[67] * bt[6][1])
                            * cos(toRadians(lon - pb[97])) * apt * swc[10]
                        + (pb[55] * bt[1][0] + pb[56] * bt[3][0] + pb[57] * bt[5][0])
                              * cos(sr * (sec - pb[58])) * apt * swc[11];
          }

          if (ww[1] != 9898.0) {
            wb[1][12] = (pb[52] * bp[2][1] + pb[98] * bp[4][1] + pb[67] * bp[6][1])
                        * cos(toRadians(lon - pb[97] + 90.0)) * apt * swc[10];
          }

          if (ww[0] != 9898.0) {
            wc[0][12] = swc[10] * wc[0][10] * pc[52] * apt + swc[11] * wc[0][11] * pc[55] * apt;
          }

          if (ww[1] != 9898.0) {
            wc[1][12] = swc[10] * wc[1][10] * pc[52] * apt + swc[11] * wc[1][11] * pc[55] * apt;
          }
        }
      }
    }
  }

  greal wbt[2] = {0.0, 0.0};
  greal wct[2] = {0.0, 0.0};

  // sum winds and change meridional sign to + North
  for (int k = 0; k < nsw; ++k) {
    wbt[0] = wbt[0] - fabs(sw[k]) * wb[0][k];
    wct[0] = wct[0] - fabs(sw[k]) * wc[0][k];
    wbt[1] = wbt[1] + fabs(sw[k]) * wb[1][k];
    wct[1] = wct[1] + fabs(sw[k]) * wc[1][k];
  }

  if (ww[0] != 9898.0) {
    ww[0] = wbt[0] * sw[23] + wct[0] * sw[24];
  }

  if (ww[1] != 9898.0) {
    ww[1] = wbt[1] * sw[23] + wct[1] * sw[24];
  }
}

//! \brief glbw5m member function from HWM class
void HWM::glbw5m(greal yrd, greal sec, greal lat, greal lon, greal stl, greal f107a, greal f107,
                 const greal ap[], const greal pb[], const greal pc[], greal ww[])
{
  const greal hr = 0.2618, dr = 1.72142e-02, pset = 4.0;
  greal wb[2][15] = {{0.0}}, wc[2][15] = {{0.0}};
  const int nsw = 14, lv = 10, mv = 2, nsv = 2;

  // Confirm parameter set
  assert(pb[99] == pset);

  for (int j = 0; j < nsw; ++j) {
    wb[0][j] = 0.0;
    wb[1][j] = 0.0;
    wc[0][j] = 0.0;
    wc[1][j] = 0.0;
  }

  greal sw9 = 1.0;
  if (sw[8] > 0.0) {
    sw9 = 1.0;
  }

  if (sw[8] < 0.0) {
    sw9 = -1.0;
  }

  int iyr = int(yrd) / 1000;
  greal day = yrd - iyr * 1000.0;

  if ((xvl != lat) || (lv > lvl) || (mv > mvl)) {
    slat = sin(toRadians(lat));
    clat = cos(toRadians(lat));

    vsphr1(slat, clat, lv, mv, bt, bp, 20);

    xvl = lat;
    lvl = lv;
    mvl = mv;
  }

  if ((tll != stl) || (nsv > nsvl)) {
    sstl = sin(hr * stl);
    cstl = cos(hr * stl);
    s2stl = sin(2.0 * hr * stl);
    c2stl = cos(2.0 * hr * stl);
    tll = stl;
    nsvl = nsv;
  }

  greal cd14 = cos(dr * (day - pb[13]));
  greal cd18 = cos(2.0 * dr * (day - pb[17]));
  greal cd19b = cos(2.0 * dr * (day - pb[18]));

  // f10.7 effect
  greal df = f107 - f107a;
  greal dfa = f107a - 150.0;
  greal dfc = dfa + pb[19] * df;

  // Time independent
  greal f1b = 1.0;

  if (ww[0] != 9898.0) {
    wb[0][1] = (pb[1] * bt[2][0] + pb[2] * bt[4][0] + pb[22] * bt[6][0]) * f1b;
  }

  wb[1][1] = 0.0;
  wc[0][1] = 0.0;

  greal f1c = 1.0;

  if (ww[1] != 9898.0) {
    wc[1][1] = -(pc[1] * bt[1][0] + pc[2] * bt[3][0] + pc[22] * bt[5][0]) * f1c
               - (pc[26] * bt[2][0] + pc[14] * bt[4][0] + pc[59] * bt[6][0]) * f1c;
  }

  // symmetrical annual
  // symmetrical semiannual
  if (ww[0] != 9898.0) {
    wb[0][3] = (pb[16] * bt[2][0] + pb[30] * bt[4][0]) * cd18;
  }

  wb[1][3] = 0.0;
  wc[0][3] = 0.0;

  if (ww[1] != 9898.0) {
    wc[1][3] = -(pc[16] * bt[1][0] + pc[30] * bt[3][0]) * cd18;
  }

  // Asymmetrical annual
  greal f5b = 1.0;

  if (ww[0] != 9898.0) {
    wb[0][4] = (pb[9] * bt[1][0] + pb[10] * bt[3][0]) * cd14 * f5b;
  }

  wb[1][4] = 0.0;
  wc[0][4] = 0.0;

  greal f5c = 1.0;

  if (ww[1] != 9898.0) {
    wc[1][4] = -(pc[9] * bt[2][0] + pc[10] * bt[4][0]) * cd14 * f5c;
  }

  // Asymmetrical semiannual
  // Diurnal
  if (sw[6] != 0.0) {
    greal f7b = 1.0;
    greal f75b = 1.0;

    if (ww[0] != 9898.0) {
      wb[0][6] = (pb[6] * bt[1][1] + pb[7] * bt[3][1] + pb[28] * bt[5][1] + pb[88] * bt[2][1])
                     * sstl * f7b
                 + (pb[12] * bt[2][1] + pb[145] * bt[4][1]) * cd14 * sstl * f75b * swc[4]
                 + (pb[3] * bt[1][1] + pb[4] * bt[3][1] + pb[27] * bt[5][1] + pb[87] * bt[2][1])
                       * cstl * f7b
                 + (pb[11] * bt[2][1] + pb[144] * bt[4][1]) * cd14 * cstl * f75b * swc[4];
    }

    if (ww[1] != 9898.0) {
      wb[1][6] = -(pb[3] * bp[1][1] + pb[4] * bp[3][1] + pb[27] * bp[5][1] + pb[87] * bp[2][1])
                     * sstl * f7b
                 - (pb[11] * bp[2][1] + pb[144] * bp[4][1]) * cd14 * sstl * f75b * swc[4]
                 + (pb[6] * bp[1][1] + pb[7] * bp[3][1] + pb[28] * bp[5][1] + pb[88] * bp[2][1])
                       * cstl * f7b
                 + (pb[12] * bp[2][1] + pb[145] * bp[4][1]) * cd14 * cstl * f75b * swc[4];
    }

    greal f7c = 1.0;
    greal f75c = 1.0;

    if (ww[0] != 9898.0) {
      wc[0][6] = -(pc[3] * bp[2][1] + pc[4] * bp[4][1] + pc[27] * bp[6][1] + pc[87] * bp[1][1]
                   + pc[140] * bp[8][1] + pc[142] * bp[10][1])
                     * sstl * f7c
                 - (pc[11] * bp[1][1] + pc[144] * bp[3][1]) * cd14 * sstl * f75c * swc[4]
                 + (pc[6] * bp[2][1] + pc[7] * bp[4][1] + pc[28] * bp[6][1] + pc[88] * bp[1][1]
                    + pc[141] * bp[8][1] + pc[143] * bp[10][1])
                       * cstl * f7c
                 + (pc[12] * bp[1][1] + pc[145] * bp[3][1]) * cd14 * cstl * f75c * swc[4];
    }

    if (ww[1] != 9898.0) {
      wc[1][6] = -(pc[6] * bt[2][1] + pc[7] * bt[4][1] + pc[28] * bt[6][1] + pc[88] * bt[1][1]
                   + pc[141] * bt[8][1] + pc[143] * bt[10][1])
                     * sstl * f7c
                 - (pc[12] * bt[1][1] + pc[145] * bt[3][1]) * cd14 * sstl * f75c * swc[4]
                 - (pc[3] * bt[2][1] + pc[4] * bt[4][1] + pc[27] * bt[6][1] + pc[87] * bt[1][1]
                    + pc[140] * bt[8][1] + pc[142] * bt[10][1])
                       * cstl * f7c
                 - (pc[11] * bt[1][1] + pc[144] * bt[3][1]) * cd14 * cstl * f75c * swc[4];
    }
  }

  // semidiurnal
  if (sw[7] != 0.0) {
    greal f8b = 1.0 + pb[89] * dfc * swc[0];

    if (ww[0] != 9898.0) {
      wb[0][7] = (pb[8] * bt[2][2] + pb[42] * bt[4][2] + pb[110] * bt[6][2] + pb[97] * bt[3][2]
                  + (pb[33] * bt[3][2] + pb[147] * bt[5][2]) * cd14 * swc[4]
                  + (pb[36] * bt[3][2]) * cd19b * swc[5])
                     * s2stl * f8b
                 + (pb[5] * bt[2][2] + pb[41] * bt[4][2] + pb[109] * bt[6][2] + pb[95] * bt[3][2]
                    + (pb[23] * bt[3][2] + pb[146] * bt[5][2]) * cd14 * swc[4]
                    + (pb[35] * bt[3][2]) * cd19b * swc[5])
                       * c2stl * f8b;
    }

    if (ww[1] != 9898.0) {
      wb[1][7] = -(pb[5] * bp[2][2] + pb[41] * bp[4][2] + pb[109] * bp[6][2] + pb[95] * bp[3][2]
                   + (pb[23] * bp[3][2] + pb[146] * bp[5][2]) * cd14 * swc[4]
                   + (pb[35] * bp[3][2]) * cd19b * swc[5])
                     * s2stl * f8b
                 + (pb[8] * bp[2][2] + pb[42] * bp[4][2] + pb[110] * bp[6][2] + pb[97] * bp[3][2]
                    + (pb[33] * bp[3][2] + pb[147] * bp[5][2]) * cd14 * swc[4]
                    + (pb[36] * bp[3][2]) * cd19b * swc[5])
                       * c2stl * f8b;
    }

    greal f8c = 1.0 + pc[89] * dfc * swc[0];

    if (ww[0] != 9898.0) {
      wc[0][7] = -(pc[5] * bp[3][2] + pc[41] * bp[5][2] + pc[109] * bp[7][2] + pc[95] * bp[2][2]
                   + (pc[23] * bp[2][2] + pc[146] * bp[4][2]) * cd14 * swc[4]
                   + (pc[35] * bp[2][2]) * cd19b * swc[5])
                     * s2stl * f8c
                 + (pc[8] * bp[3][2] + pc[42] * bp[5][2] + pc[110] * bp[7][2] + pc[97] * bp[2][2]
                    + (pc[33] * bp[2][2] + pc[147] * bp[4][2]) * cd14 * swc[4]
                    + (pc[36] * bp[2][2]) * cd19b * swc[5])
                       * c2stl * f8c;
    }

    if (ww[1] != 9898.0) {
      wc[1][7] = -(pc[8] * bt[3][2] + pc[42] * bt[5][2] + pc[110] * bt[7][2] + pc[97] * bt[2][2]
                   + (pc[33] * bt[2][2] + pc[147] * bt[4][2]) * cd14 * swc[4]
                   + (pc[36] * bt[2][2]) * cd19b * swc[5])
                     * s2stl * f8c
                 - (pc[5] * bt[3][2] + pc[41] * bt[5][2] + pc[95] * bt[2][2] + pc[109] * bt[7][2]
                    + (pc[23] * bt[2][2] + pc[146] * bt[4][2]) * cd14 * swc[4]
                    + (pc[35] * bt[2][2]) * cd19b * swc[5])
                       * c2stl * f8c;
    }
  }

  // Terdiurnal
  // Magnetic activity
  greal apd, apdf, apdfc, apt;
  if (sw[8] != 0.0) {
    if (sw9 != -1.0) {

      // Daily ap
      apd = ap[0] - 4.0;
      apdf = (apd + (pb[44] - 1.0) * (apd + (exp(-pb[43] * apd) - 1.0) / pb[43]));
      // apdfc=(apd  +  (pc[44]  -  1.0)*(apd  +  (exp(  -  pc[43]*apd)  -  1.0)/pc[43]));
      apdfc = apdf;

      if (apd != 0) {

        if (ww[0] != 9898.0) {
          wb[0][8] = (pb[45] * bt[2][0] + pb[34] * bt[4][0]) * apdf
                     + (pb[121] * bt[1][1] + pb[122] * bt[3][1] + pb[123] * bt[5][1])
                           * cos(hr * (stl - pb[124])) * apdf * swc[6];
        }

        if (ww[1] != 9898.0) {
          wb[1][8] = (pb[121] * bp[1][1] + pb[122] * bp[3][1] + pb[123] * bp[5][1])
                     * cos(hr * (stl - pb[124] + 6.0)) * apdf * swc[6];
        }

        if (ww[0] != 9898.0) {
          wc[0][8] = (pc[121] * bp[2][1] + pc[122] * bp[4][1] + pc[123] * bp[6][1])
                     * cos(hr * (stl - pc[124])) * apdfc * swc[6];
        }

        if (ww[1] != 9898.0) {
          wc[1][8] = -(pc[45] * bt[1][0] + pc[34] * bt[3][0]) * apdfc
                     + (pc[121] * bt[2][1] + pc[122] * bt[4][1] + pc[123] * bt[6][1])
                           * cos(hr * (stl - pc[124] + 6.0)) * apdfc * swc[6];
        }
      }
    }
    else {
      //if (pb[24] < 1.0e-04) {
      //	pb[24] = 1.0e-04;
      //}

      //apt = g0(ap[1], pb);
      apt = 0.0;
      if (apt != 0.0) {

        if (ww[0] != 9898.0) {
          wb[0][8] = (pb[96] * bt[2][0] + pb[54] * bt[4][0]) * apt
                     + (pb[128] * bt[1][1] + pb[129] * bt[3][1] + pb[130] * bt[5][1])
                           * cos(hr * (stl - pb[131])) * apt * swc[6];
        }

        if (ww[1] != 9898.0) {
          wb[1][8] = (pb[128] * bp[1][1] + pb[129] * bp[3][1] + pb[130] * bp[5][1])
                     * cos(hr * (stl - pb[131] + 6.0)) * apt * swc[6];
        }

        if (ww[0] != 9898.0) {
          wc[0][8] = (pc[128] * bp[2][1] + pc[129] * bp[4][1] + pc[130] * bp[6][1])
                     * cos(hr * (stl - pc[131])) * apt * swc[6];
        }

        if (ww[1] != 9898.0) {
          wc[1][8] = -(pc[96] * bt[1][0] + pc[54] * bt[3][0]) * apt
                     + (pc[128] * bt[2][1] + pc[129] * bt[4][1] + pc[130] * bt[6][1])
                           * cos(hr * (stl - pc[131] + 6.0)) * apt * swc[6];
        }
      }
    }
  }

  greal wbt[2] = {0.0, 0.0};
  greal wct[2] = {0.0, 0.0};

  // sum winds and change meridional sign to + North
  for (int k = 0; k < nsw; ++k) {
    wbt[0] = wbt[0] - fabs(sw[k]) * wb[0][k];
    wct[0] = wct[0] - fabs(sw[k]) * wc[0][k];
    wbt[1] = wbt[1] + fabs(sw[k]) * wb[1][k];
    wct[1] = wct[1] + fabs(sw[k]) * wc[1][k];
  }

  if (ww[0] != 9898.0) {
    ww[0] = wbt[0] * sw[23] + wct[0] * sw[24];
  }

  if (ww[1] != 9898.0) {
    ww[1] = wbt[1] * sw[23] + wct[1] * sw[24];
  }
}

} // namespace
