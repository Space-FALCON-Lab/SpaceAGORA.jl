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
#include <algorithm>
#include "JB2008.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
JB2008::JB2008()
  : Atmosphere(), EarthCommon(this)
{
}

//! \fn JB2008::JB2008(const JB2008& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn JB2008::~JB2008()
//! \copydoc Atmosphere::~Atmosphere()

//! \copydoc PerturbedAtmosphere::setInputParameters()
void JB2008::setInputParameters(const EarthInputParameters& params)
{
  caltojul(params.year, params.month, params.day, 0, 0, 0.0, xmjd);

  year = params.year;
  month = params.month;
  day = params.day;
  hour = params.hour;       
  minute = params.minute;   
  seconds = params.seconds; 

  f10 = params.dailyF10;
  f10b = params.meanF10;
  ap = params.ap;
  s10 = params.dailyS10;
  s10b = params.meanS10;
  xm10 = params.dailyXM10;
  xm10b = params.meanXM10;
  y10 = params.dailyY10;
  y10b = params.meanY10;
  dstdtc = params.dstdtc;

  if (s10 <= 0.0 || s10b == 0.0) {
    s10 = f10;
    s10b = f10b;
  }

  if (xm10 <= 0.0 || xm10b == 0.0) {
    xm10 = f10;
    xm10b = f10b;
  }

  if (y10 <= 0.0 || y10b == 0.0) {
    y10 = f10;
    y10b = f10b;
  }

  if (dstdtc <= 0.0) {
    dstdtc = 0.0;
  }
}

//! \fn JB2008::setDayOfYear()
//! \brief Sets the day of the year.
//!
//! \param doy   Day of the year (0-367).

//! \fn JB2008::setJulianDay()
//! \brief Sets the Julian day.
//!
//! \param jDay   Julian day.

//! \fn JB2008::getTemperatureGradient()
//! \brief Gets the temperature gradient with respect to height.
//!
//! \returns Temperature gradient with respect to height \units{K}.

//! \fn JB2008::getMolecularWeightGradient()
//! \brief Gets the molecular weight gradient with respect to height.
//!
//! \returns Molecular weight gradient with respect to height.

//! \brief JB2008 driver routine to evaluate mean p, d, t, u, v, w
//!
//! \b Inputs
//! \arg #position          
//!
//! \returns  pressure, density, temperature, ewWind, nsWind, verticalWind, 
//! gas number densities, averageMolecularWeight, temperature gradient, molecular weight gradient
void JB2008::update()
{
  int imin = (int)minute + int(elapsedTime) / 60;
  greal sec = seconds + elapsedTime;
  sec = int(sec) % 60;
  imin = imin % 60;
  greal xmin = imin + sec / 60.0;

  // Distances for 5 degrees of geocentric latitude, longitude
  greal dy5 = 5000.0 * toRadians(totalRadius);
  greal dx5 = dy5 * cos(toRadians(latitude));
  dx5 = max(2000.0, dx5);

  // Following is the pure thermospheric height range section
  // JB2008 and HWM values at current position

  greal tex;
  JB08HWM(xmin, height, latitude, longitude, pressure, density, temperature,
          dinitrogen.numberDensity, dioxygen.numberDensity, oxygen.numberDensity,
          argon.numberDensity, helium.numberDensity, hydrogen.numberDensity, nitrogen.numberDensity,
          averageMolecularWeight, tex, ewWind, nsWind);
 
  // Set geocentric latitude increment for temperature gradients
  greal dphij = 5.0;
  if (latitude > 85.0) {
    dphij = -5.0;
  }

  // JB2008 temperature at current position plus geocentric latitude increment
  greal thn, dhn, phn, d1, d2, d4, d6, d3, d5, d7, d8, d10, d9, d11;
  JB08HWM(xmin, height, latitude + dphij, longitude, phn, dhn, thn, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11);

  greal the, phe, dhe;
  // JB2008 temperature at current position plus 5 degrees lon
  JB08HWM(xmin, height, latitude, longitude + 5.0, phe, dhe, the, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11);

  // dt/dx, dt/dy and dt/dz for vertical wind
  greal dtx = the - temperature;
  greal dty = thn - temperature;
  if (dphij < 0.0) {
    dty = -dty;
  }

  // JB2008 temperature and molecular weight at current lat-lon, 1 km higher
  greal tb, pb, db;
  JB08HWM(xmin, height + 1.0, latitude, longitude, pb, db, tb, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11);

  // Gradients for temperature and molecular weight
  dtz = (tb - temperature) / 1000.0;
  dmdz = (d8 - averageMolecularWeight) / 1000.0;

  // Compute vertical mean wind
  // Specific heat
  greal cp = 7.0 * pressure / (2.0 * density * temperature);

  // Mean vertical wind from Montgomery stream function
  verticalWind = -cp * (ewWind * dtx / dx5 + nsWind * dty / dy5) / (gravity + cp * dtz);
}

//! \brief Unknown
//!
//! \param z Height.
//!
//! \returns unknown
greal JB2008::xambar(greal z)
{
	//Evalutes equation (1) 
	const greal c[7] = { 28.15204, -8.5586e-2, 1.284e-4, -1.0056e-5, -1.021e-5, 1.5044e-6, 9.9826e-8 };

	greal dz = z - 100.0;
	greal amb = c[6];

	for (int i = 0; i < 6; i++){
		int j = 5 - i;
    amb = dz * amb + c[j];
	}
	return amb;
}

//! \brief Compute semiannual variation (delta log rho).
//!
//! \param iyr
//! \param dayOfYear
//! \param ht
//! \param f10b
//! \param s10b
//! \param xm10b
//! \param[out] fzz
//! \param[out] gtz
//! \param[out] drlog
void JB2008::semian08(int iyr, greal dayOfYear, greal ht, greal f10b, greal s10b, greal xm10b,
                      greal &fzz, greal &gtz, greal &drlog)
{
  greal const gtm[10] = {-0.3633,     0.8506e-01, 0.2401,      -0.1897,     -0.2554,
                          -0.1790e-01, 0.565e-03,  -0.6407e-03, -0.3418e-02, -0.1252e-02};
  greal const fzm[5] = {0.2689, -0.1176e-01, 0.2782e-01, -0.2782e-01, 0.3470e-03};

	//Compute new 81-day centered solar index for fz
  greal fsmb = 1.00 * f10b - 0.7 * s10b - 0.04 * xm10b;

	greal htz = ht / 1000.0;

	fzz = fzm[0] + fzm[1] * fsmb + fzm[2] * fsmb * htz + fzm[3] * fsmb * pow(htz, 2)
        + fzm[4] * pow(fsmb, 2) * htz;

	//Compute daily 81-day centered solar index for gt
	fsmb = 1.00*f10b - 0.75*s10b - 0.37*xm10b;

	greal yrln = 365.0;
  if (iyr % 4 == 0) {
    yrln = 366.0;
  }

	greal tau = (dayOfYear - 1.0) / yrln;
  greal sin1p = sin(TWO_PI * tau);
  greal cos1p = cos(TWO_PI * tau);
  greal sin2p = sin(2.0 * TWO_PI * tau);
  greal cos2p = cos(2.0 * TWO_PI * tau);

	gtz = gtm[0] + gtm[1] * sin1p + gtm[2] * cos1p + gtm[3] * sin2p +
		gtm[4] * cos2p + gtm[5] * fsmb + gtm[6] * fsmb*sin1p + gtm[7] * fsmb*cos1p
		+ gtm[8] * fsmb*sin2p + gtm[9] * fsmb*cos2p;

	if (fzz < 1.0e-6) {
    fzz = 1.0e-06;
  }

	drlog = fzz * gtz;
}

//! \brief Compute dTc correction for Jacchia-Bowman model.
//!
//! \param f10
//! \param xlst
//! \param xlat
//! \param zht
//! \param[out] dtc
void JB2008::dtsub(greal f10, greal xlst, greal xlat, greal zht, greal &dtc)
{

	greal const b[19] = { -0.457512297e+01, -0.512114909e+01, -0.693003609e+02, 0.203716701e+03,
       0.703316291e+03, -0.194349234e+04,  0.110651308e+04, -0.174378996e+03, 0.188594601e+04,
      -0.709371517e+04,  0.922454523e+04, -0.384508073e+04, -0.645841789e+01, 0.409703319e+02,
      -0.482006560e+03,  0.181870931e+04, -0.237389204e+04,  0.996703815e+03, 0.361416936e+02 };

	greal const c[23] = { -0.155986211e+02, -0.512114909e+01, -0.693003609e+02,  0.203716701e+03,
       0.703316291e+03, -0.194349234e+04,  0.110651308e+04, -0.220835117e+03,  0.143256989e+04,
      -0.318481844e+04,  0.328981513e+04, -0.135332119e+04,  0.199956489e+02, -0.127093998e+02,
       0.212825156e+02, -0.275555432e+01,  0.110234982e+02,  0.148881951e+03, -0.751640284e+03,
       0.637876542e+03,  0.127093998e+02, -0.212825156e+02,  0.275555432e+01 };

	dtc = 0.0;
  greal ycs = cos(xlat);
  greal f = (f10 - 100.0) / 100.0;
  greal tx = xlst / 24.0;
  greal tx2 = tx * tx;
  greal tx3 = tx * tx2;
  greal tx4 = tx * tx3;
  greal tx5 = tx * tx4;

	//Calculate dTc
	if ((zht >= 120.0) && (zht <= 200.0)) {
    //greal h = (zht - 200.0) / 50.0;
    greal dtc200 = c[16] 
                    + (c[17] * tx + c[18] * tx2 + c[19] * tx3) * ycs
                    + (c[20] + c[21] * tx + c[22] * tx2) * f * ycs;
    greal sum = c[0] 
                 + (b[1] + c[2] * tx + c[3] * tx2 + c[4] * tx3 + c[5] * tx4 + c[6] * tx5) * f
                 + (c[7] * tx + c[8] * tx2 + c[9] * tx3 + c[10] * tx4 + c[11] * tx5) * ycs
                 + c[12] * ycs 
                 + (c[13] + c[14] * tx + c[15] * tx2) * f * ycs;
    greal dtc200dz = sum;
    greal cc = 3.0 * dtc200 - dtc200dz;
    greal dd = dtc200 - cc;
    greal zp = (zht - 120.0) / 80.0;
    dtc = cc * zp * zp + dd * zp * zp * zp;

	}

	if ((zht > 200.0) & (zht <= 240.0)) {
    greal h = (zht - 200.0) / 50.0;
    greal sum = c[0] * h
                 + (b[1] + c[2] * tx + c[3] * tx2 + c[4] * tx3 + c[5] * tx4 + c[6] * tx5) * f * h
                 + (c[7] * tx + c[8] * tx2 + c[9] * tx3 + c[10] * tx4 + c[11] * tx5) * ycs * h
                 + c[12] * ycs * h 
                 + (c[13] + c[14] * tx + c[15] * tx2) * f * ycs * h 
                 + c[16]
                 + (c[17] * tx + c[18] * tx2 + c[19] * tx3) * ycs
                 + (c[20] + c[21] * tx + c[22] * tx2) * f * ycs;
    dtc = sum;
  }

	if ((zht > 240.0) & (zht <= 300.0)) {
    greal h = (40.0) / 50.0;
    greal sum = c[0] * h
                 + (b[1] + c[2] * tx + c[3] * tx2 + c[4] * tx3 + c[5] * tx4 + c[6] * tx5) * f * h
                 + (c[7] * tx + c[8] * tx2 + c[9] * tx3 + c[10] * tx4 + c[11] * tx5) * ycs * h
                 + c[12] * ycs * h 
                 + (c[13] + c[14] * tx + c[15] * tx2) * f * ycs * h
                 + c[16]
                 + (c[17] * tx + c[18] * tx2 + c[19] * tx3) * ycs
                 + (c[20] + c[21] * tx + c[22] * tx2) * f * ycs;
    greal aa = sum;
    greal bb = c[0] 
                + (b[1] + c[2] * tx + c[3] * tx2 + c[4] * tx3 + c[5] * tx4 + c[6] * tx5) * f
                + (c[7] * tx + c[8] * tx2 + c[9] * tx3 + c[10] * tx4 + c[11] * tx5) * ycs
                + c[12] * ycs 
                + (c[13] + c[14] * tx + c[15] * tx2) * f * ycs;
    h = 300.0 / 100.0;
    sum = b[0] 
          + (b[1] + b[2] * tx + b[3] * tx2 + b[4] * tx3 + b[5] * tx4 + b[6] * tx5) * f
          + (b[7] * tx + b[8] * tx2 + b[9] * tx3 + b[10] * tx4 + b[11] * tx5) * ycs
          + (b[12] + b[13] * tx + b[14] * tx2 + b[15] * tx3 + b[16] * tx4 + b[17] * tx5) * h * ycs
          + b[18] * ycs;
    greal dtc300 = sum;
    sum = (b[12] + b[13] * tx + b[14] * tx2 + b[15] * tx3 + b[16] * tx4 + b[17] * tx5) * ycs;
    greal dtc300dz = sum;
    greal cc = 3.0 * dtc300 - dtc300dz - 3.0 * aa - 2.0 * bb;
    greal dd = dtc300 - aa - bb - cc;
    greal zp = (zht - 240.0) / 60.0;
    dtc = aa + bb * zp + cc * zp * zp + dd * zp * zp * zp;
  }

	if ((zht > 300) & (zht <= 600.0)) {
    greal h = zht / 100.0;
    greal sum = b[0] 
                 + (b[1] + b[2] * tx + b[3] * tx2 + b[4] * tx3 + b[5] * tx4 + b[6] * tx5) * f
                 + (b[7] * tx + b[8] * tx2 + b[9] * tx3 + b[10] * tx4 + b[11] * tx5) * ycs
                 + b[12] * h * ycs
                 + (b[13] * tx + b[14] * tx2 + b[15] * tx3 + b[16] * tx4 + b[17] * tx5) * h * ycs
                 + b[18] * ycs;
    dtc = sum;
  }

	if ((zht > 600.0) & (zht <= 800.0)) {
    greal zp = (zht - 600.0) / 100.0;
    greal hp = 600.0 / 100.0;
    greal aa = b[0] 
                + (b[1] + b[2] * tx + b[3] * tx2 + b[4] * tx3 + b[5] * tx4 + b[6] * tx5) * f
                + (b[7] * tx + b[8] * tx2 + b[9] * tx3 + b[10] * tx4 + b[11] * tx5) * ycs
                + b[12] * hp * ycs
                + (b[13] * tx + b[14] * tx2 + b[15] * tx3 + b[16] * tx4 + b[17] * tx5) * hp * ycs
                + b[18] * ycs;
    greal bb = (b[12] + b[13] * tx + b[14] * tx2 + b[15] * tx3 + b[16] * tx4 + b[17] * tx5) * ycs;
    greal cc = -(3.0 * aa + 4.0 * bb) / 4.0;
    greal dd = (aa + bb) / 4.0;
    dtc = aa + bb * zp + cc * zp * zp + dd * zp * zp * zp;
  }
}

//! \brief Compute day and year from the d1950 (days since 1950).
//!
//! \param d1950
//! \param[out] iyr
//! \param[out] days
void JB2008::tmoutd(greal d1950, int &iyr, greal &days)
{
	int iyday = (int)d1950;
  greal frac = d1950 - iyday;
	iyday = iyday + 364;
  int itemp = iyday / 1461;
	iyday = iyday - itemp * 1461;
	iyr = 1949 + 4 * itemp;
  itemp = iyday / 365;
	if (itemp >= 3) {
    itemp = 3;
  }
	iyr = iyr + itemp;
	iyday = iyday - 365 * itemp + 1;
	iyr = iyr - 1900;
	days = iyday + frac;
	if (iyr >= 100) {
    iyr = iyr - 100;
  }
}

//! \brief Evaluates equation (10) or equation (13), depending on z
//!
//! \param z
//! \param tc
//! 
//! \returns unknown
greal JB2008::xlocal(greal z, const greal tc[])
{
	greal xlocal_out;
	greal dz = z - 125.0;

  if (dz > 0.0) {
    xlocal_out = tc[0] + tc[2] * atan(tc[3] * dz * (1.0 + 4.5e-06 * pow(dz, 2.5)));
  }
  else {
    xlocal_out = ((-9.8204695e-06 * dz - 7.3039742e-04) * pow(dz, 2) + 1) * dz * tc[1] + tc[0];
  }
  return xlocal_out;
}

//! \brief Evaluates equation (8)
//!
//! \param z
//! 
//! \returns unknown
greal JB2008::xgrav(greal z)
{
	return referenceGravity / pow((1.0 + z / 6356.766), 2);
}


//! \brief Jacchia-Bowman 2008 Model Atmosphere.
//!
//! This is the CIRA "Integration Form" of a Jacchia Model.  There are no tabular
//! values of density.  Instead, the barometric equation and diffusion equation are
//! integrated numerically using the Newton-Coates method to produce the density
//! profile up to the input position.
//!
//! \param amjd
//! \param sun
//! \param sat
//! \param f10
//! \param f10b
//! \param s10
//! \param s10b
//! \param xm10
//! \param xm10b
//! \param y10
//! \param y10b
//! \param dstdtc
//! \param[out] temp
//! \param[out] rho
//! \param[out] pres
//! \param[out] avgmw
//! \param[out] and1
//! \param[out] sumn
void JB2008::JB08(greal amjd, greal sun[], greal sat[], greal f10,
	greal f10b, greal s10, greal s10b, greal xm10, greal xm10b,
	greal y10, greal y10b, greal dstdtc, greal temp[], greal &rho,
	greal &pres, greal &avgmw, greal and1[], greal &sumn)
{

	//The alpha are the thermal diffusion coefficients in Eq. (6)
	const greal alpha[5] = { 0.0, 0.0, 0.0, 0.0, -0.38 };

	//al10 is log(10.0)
	const greal al10 = 2.3025851;

	//The AMW are the molecular weights in order:  N2, O2, O, Ar, He, and H
	const greal amw[6] = { 28.0134, 31.9988, 15.9994, 39.948, 4.0026, 1.00797 };

  //avogad is Avogadro's number in mks units (molecules/kmol)
  const double avogad = AVOGADRO * 1000.0; 
  
  //the frac are assumed sea-level volume fractions in order:  N2, O2, Ar, and He
	const greal frac[4] = { 0.7811, 0.20955, 9.34e-03, 1.289e-05 };

	//rstar is the universal gas-constant in mks units (joules/K/kmol)
  //const greal rstar = 8314.32;
  const greal rstar = UNIVERSAL_GAS * 1000.0;

	//The R# are values used to establish height step sizes in the regimes 90km
	//to 105km, 105km to 500km and 500km upward.
	const greal r1 = 0.01, r2 = 0.025, r3 = 0.075;

	//The wt are weights for the Newton-Coates five-point quad. formula
	const greal wt[5] = { 0.311111111111111, 1.422222222222222, 0.533333333333333,
		                    1.422222222222222, 0.311111111111111 };

	//The cht are coefficients for high altitude density correction
	const greal cht[4] = { 0.22, -0.20e-02, 0.115e-02, -0.211e-05 };

	//Equation (14)
  greal fn = pow((f10b / 240.0), (1.0 / 4.0));
	fn = min(fn, 1.0);
  greal fsb = f10b * fn + s10b * (1.0 - fn);
  greal tsubc = 392.4 + 3.227 * fsb + 0.298 * (f10 - f10b) + 2.259 * (s10 - s10b)
                 + 0.312 * (xm10 - xm10b) + 0.178 * (y10 - y10b);

	//Equation (15)
  greal eta = 0.5 * abs(sat[1] - sun[1]);
  greal theta = 0.5 * abs(sat[1] + sun[1]);

	//Equation (16)
  greal h = sat[0] - sun[0];
  greal tau = h - 0.64577182 + 0.10471976 * sin(h + 0.75049158);
  greal glat = sat[1];
  greal zht = sat[2];
  greal glst = h + PI;
  greal glsthr = toDegrees(glst) * (24.0 / 360.0);
  if (glsthr >= 24.0) {
    glsthr = glsthr - 24.0;
  }
  else if (glsthr < 0.0) {
    glsthr = glsthr + 24.0;
  }

	//Equation (17)
  greal c = pow(cos(eta), 2.5);
  greal s = pow(sin(theta), 2.5);
  greal df = s + (c - s) * pow(abs(cos(0.5 * tau)), 3);
  greal tsubl = tsubc * (1.0 + 0.31 * df);

	//Compute correction to dTc for local solar time and lat correction
  greal dtclst; 
	dtsub(f10, glsthr, glat, zht, dtclst);

	//Compute the local exospheric temperature.
	//Add geomagnetic storm effect from input dTc value

	temp[0] = tsubl + dstdtc;
  greal tinf = tsubl + dstdtc + dtclst;

	//Equation (9)
  greal tsubx = 444.3807 + 0.02385 * tinf - 392.8292 * exp(-0.0021357 * tinf);

	//Equation (11)
  greal gsubx = 0.054285714 * (tsubx - 183.0);

	//The tc array will be an argument in the call to xlocal, which evaluates 
	//Equation (10) and Equation (13)
  greal tc[4];
	tc[0] = tsubx;
	tc[1] = gsubx;

	//a and gsubx/a of equation (13)
  tc[2] = (tinf - tsubx) / HALF_PI;
	tc[3] = gsubx / tc[2];

	//Equation (5)
  greal z1 = 90.0;

	greal z2 = min(sat[2], 105.0);

	greal al = log(z2 / z1);

	greal oor = 1.0 / r1;

	int n = int(al * oor) + 1;

	greal zr = exp(al / greal(n));
	greal ambar1 = xambar(z1);
	greal tloc1 = xlocal(z1, tc);
	greal zend = z1;
	greal sum2 = 0.0;
	greal ain = ambar1 * xgrav(z1) / tloc1;

  greal z = 0.0, ambar2 = 0.0, tloc2 = 0.0, gravl = 0.0; 
	for (int i = 0; i < n; i++){
		z = zend;
		zend = zr * z;
    greal dz = 0.25 * (zend - z);
    greal sum1 = wt[0] * ain;
		for (int j = 1; j < 5; j++){
			z = z + dz;
			ambar2 = xambar(z);
			tloc2 = xlocal(z, tc);
			gravl = xgrav(z);
			ain = ambar2*gravl / tloc2;
			sum1 = sum1 + wt[j] * ain;
		}
		sum2 = sum2 + dz*sum1;
	}

	greal fact1 = 1000.0 / rstar;
	rho = 3.46e-6*ambar2*tloc1*exp(-fact1*sum2) / ambar1 / tloc2;

	//Equation (2)
  greal anm = avogad * rho;
  greal an = anm / ambar2;

	//Equation (3)
  greal fact2 = anm / dryAirMolWt;
  greal aln[6];
	aln[0] = log(frac[0] * fact2);

	aln[3] = log(frac[2] * fact2);
	aln[4] = log(frac[3] * fact2);

	//Equation (4)
	aln[1] = log(fact2 * (1.0 + frac[1]) - an);
	aln[2] = log(2.0 * (an - fact2));

	if (sat[2] <= 105.0) {
    temp[1] = tloc2;

    // Put in negligible hydrogen for use in loop 13
    aln[5] = aln[4] - 25.0;
  }
  else {
    // Equation(6) 
    greal z3 = min(sat[2], 500.0);
    al = log(z3 / z);
    oor = 1.0 / r2;
    n = int(al * oor) + 1;
    zr = exp(al / greal(n));
    sum2 = 0.0;
    ain = gravl / tloc2;

    greal tloc3 = 0.0; 
    for (int i = 0; i < n; i++) {
      z = zend;
      zend = zr * z;
      greal dz = 0.25 * (zend - z);
      greal sum1 = wt[0] * ain;
      for (int j = 1; j < 5; j++) {
        z = z + dz;
        tloc3 = xlocal(z, tc);
        gravl = xgrav(z);
        ain = gravl / tloc3;
        sum1 = sum1 + wt[j] * ain;
      }
      sum2 = sum2 + dz * sum1;
    }

    greal z4 = max(sat[2], 500.0);
    al = log(z4 / z);
    greal r = r2;
    if (sat[2] > 500.0)
      r = r3;
    oor = 1.0 / r;
    n = int(al * oor) + 1;
    zr = exp(al / greal(n));
    greal sum3 = 0.0;
    greal tloc4 = 0.0; 
    for (int i = 0; i < n; i++) {
      z = zend;
      zend = zr * z;
      greal dz = 0.25 * (zend - z);
      greal sum1 = wt[0] * ain;
      for (int j = 1; j < 5; j++) {
        z = z + dz;
        tloc4 = xlocal(z, tc);
        gravl = xgrav(z);
        ain = gravl / tloc4;
        sum1 = sum1 + wt[j] * ain;
      }
      sum3 = sum3 + dz * sum1;
    }
    greal altr, hsign; 
    if (sat[2] <= 500.0) {
      temp[1] = tloc3;
      altr = log(tloc3 / tloc2);
      fact2 = fact1 * sum2;
      hsign = 1.0;
    }
    else {
      temp[1] = tloc4;
      altr = log(tloc4 / tloc2);
      fact2 = fact1 * (sum2 + sum3);
      hsign = -1.0;
    }

    for (int i = 0; i < 5; i++) {
      aln[i] = aln[i] - (1.0 + alpha[i]) * altr - fact2 * amw[i];
    }

    // Equation (7) - Note that in CIRA72, al10t5 = log10(t500)
    greal al10t5 = log10(tinf);
    greal alnh5 = (5.5 * al10t5 - 39.40) * al10t5 + 73.13;
    aln[5] = al10 * (alnh5 + 6.0) + hsign * (log(tloc4 / tloc3) + fact1 * sum3 * amw[5]);
  }

	//Equation (24) - J70 Seasonal-Latitudinal Variation
  greal trash = (amjd - 36204.0) / 365.2422;
  greal capphi = fmod(trash, 1.0);

	greal dlrsl = 0.02 * (sat[2] - 90.0) * exp(-0.045 * (sat[2] - 90.0)) * copysign(1.0, sat[1])
                     * sin(TWO_PI * capphi + 1.72) * pow(sin(sat[1]), 2);

	//Equation (23) - Computes the semiannual variational
  greal dlrsa = 0.0;
	if (z < 2000.0){
    greal d1950 = amjd - 33281.0;
    greal yrday; 
    int iyr;
		tmoutd(d1950, iyr, yrday);
		//Use new semiannual model
    greal fzz, gtz;
		semian08(iyr, yrday, zht, f10b, s10b, xm10b, fzz, gtz, dlrsa);
		if (fzz < 0.0) dlrsa = 0.0;
	}

	//Sum the delta-log-rhos and apply to the number densities.
  greal dlr = al10 * (dlrsl + dlrsa);

	for (int i = 0; i < 6; i++){
		aln[i] = aln[i] + dlr;
	}

	//Compute mass-density and mean-molecular-weight and convert number density
	//logs from natural to common.
  sumn = 0.0;
  greal sumnm = 0.0;
  //greal al10n[6];
	for (int i = 0; i < 6; i++){
		an = exp(aln[i]);
		and1[i] = an;
		sumn = sumn + an;
		sumnm = sumnm + an*amw[i];
		//al10n[i] = aln[i] / al10;
	}

	rho = sumnm / avogad;
	avgmw = sumnm / sumn;

	//Compute the high altitude exospheric density correction factor
  greal fex = 1.0;
	if ((zht >= 1000.0) && (zht < 1500.0)){
    greal zeta = (zht - 1000.0) * 0.002;
    greal zeta2 = zeta * zeta;
    greal zeta3 = zeta * zeta2;
    greal f15c = cht[0] + cht[1] * f10b + cht[2] * 1500.0 + cht[3] * f10b * 1500;
    greal f15c_zeta = (cht[2] + cht[3] * f10b) * 500.0;
    greal fex2 = 3.0 * f15c - f15c_zeta - 3.0;
    greal fex3 = f15c_zeta - 2.0 * f15c + 2.0;
    fex = 1.0 + fex2 * zeta2 + fex3 * zeta3;
	}

	if (zht >= 1500.0){
    fex = cht[0] + cht[1] * f10b + cht[2] * zht + cht[3] * f10b * zht;
	}

	//Apply the exospheric density correction factor.
	rho = fex* rho;

	pres = rho * rstar * temp[1] / avgmw;
	for (int i = 0; i < 6; i++){
		and1[i] = and1[i] * fex;
	}

	sumn = sumn * fex;
}

//! \brief Evaluate JB2008 and HWM model for specified position.
//!
//! \param xmin
//! \param z
//! \param xlat
//! \param xlon
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
void JB2008::JB08HWM(greal xmin, greal z, greal xlat, greal xlon,
	greal &ph, greal &dh, greal &th, greal &n2nd, greal &o2nd,
	greal &ond, greal &arnd, greal &hend, greal &hnd, greal &nnd,
	greal &wtmol, greal &tex, greal &uh, greal &vh)
{
  greal dlon = xlon - longitude;
  greal sra = toRadians(solarRightAscension);
  greal sda = toRadians(solarDeclination);
  greal dd = dayOfYear;
  //Compute local RA
  greal raloc = toRadians(solarRightAscension + solarHourAngle + dlon);
  if (raloc < 0.0) raloc = raloc + (2 * PI);
  if (raloc > (2 * PI)) raloc = raloc - (2 * PI);

	//JB2008 thermosphere subroutine
  greal sun[2], sat[3];
	sun[0] = sra;
	sun[1] = sda;
	sat[0] = raloc;
	sat[1] = toRadians(xlat);
	sat[2] = z;


	greal  temp[2], rho, pres, avgmw, and1[6], sumn;
  JB08(xmjd, sun, sat, f10, f10b, s10, s10b, xm10, xm10b, y10, y10b, dstdtc, temp, rho, pres, avgmw,
       and1, sumn);

	//Change of notation for outputs
	dh = rho;
	th = temp[1];
	tex = temp[0];
	wtmol = avgmw;
	ph = pres;
	n2nd = and1[0];
	o2nd = and1[1];
	ond = and1[2];
	arnd = and1[3];
	hend = and1[4];
	hnd = and1[5];
	nnd = 0.0;

	//Compute inputs needed for HWM
	int iday = 1000 * year + int(dd + 1.0);
  greal dut = 86400.0 * (dd - int(dd));
  greal utsec = dut;
  greal xlst = solarTime;

	//Call harmonic wind model (HWM) subroutine
  greal aph[2];
  aph[0] = ap;
	aph[1] = 0.0;

  greal w[2]; 
	hwm.gws5(iday, utsec, z, xlat, xlon, xlst, f10b, f10, aph, w);

	//change notation for outputs
	uh = w[1];
	vh = w[0];
}

} // namespace
