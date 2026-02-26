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
#include "MET.h"
#include "Interpolator.h"

using namespace std;

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
MET::MET()
  : Atmosphere(), EarthCommon(this)
{
}

//! \fn MET::MET(const MET& orig)
//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)

//! \fn MET::~MET()
//! \copydoc Atmosphere::~Atmosphere()

//! \copydoc PerturbedAtmosphere::setInputParameters()
void MET::setInputParameters(const EarthInputParameters& params)
{
  year = params.year;
  daysInYear = 365.0;
  if (params.year % 4 == 0 && (params.year % 100 != 0 || params.year % 400 == 0))
  {
    daysInYear = 366.0;
  }

  month = params.month;
  day = params.day;
  hour = params.hour;
  minute = params.minute;
  seconds = params.seconds;


  f10 = params.dailyF10;
  f10b = params.meanF10;
  ap = params.ap;
}

//! \fn MET::setDayOfYear()
//! \brief Sets the day of the year.
//!
//! \param doy   Day of the year (0-367).

//! \fn MET::getTemperatureGradient()
//! \brief Gets the temperature gradient with respect to height.
//!
//! \returns Temperature gradient with respect to height \units{K}.

//! \fn MET::getMolecularWeightGradient()
//! \brief Gets the molecular weight gradient with respect to height.
//!
//! \returns Molecular weight gradient with respect to height.

//! \brief MET (Jacchia) model driver routine
//!
//! Evaluates pressure, temperature, density, winds, and constituent
//! number densities. Temperature gradient and molecular weight gradient
//! are also evaluated.  This routine was known as jacmod() in the previous
//! versions.
//!
//! \b Inputs
//! \arg #position     
//!
//! \returns  pressure, density, temperature, ewWind, nsWind, verticalWind, 
//! gas number densities, averageMolecularWeight, temperature gradient, molecular weight gradient
void MET::update()
{
  int imin = int(minute + int(elapsedTime) / 60);
  greal sec = seconds + elapsedTime;
  sec = fmod(sec, 60.0);
  imin = imin % 60;

  // Distances for 5 degrees of geocentric latitude, longitude
  greal dy5 = toRadians(5000.0 * totalRadius);
  greal dx5 = dy5 * cos(toRadians(latitude));
  dx5 = max(2000.0, dx5);

  // Following is the pure Jacchia height range section
  // Jacchia values at current position
  InData inData = { height, latitude, longitude, 2, f10, f10b, ap };
  OutData outData;
  met07(inData, outData);

  temperature = outData.th;
  pressure = outData.ph;
  density = outData.dh;

  // Atmospheric Constituents
  argon.numberDensity = outData.arnd;
  helium.numberDensity = outData.hend;
  hydrogen.numberDensity = outData.hnd;
  dioxygen.numberDensity = outData.o2nd;
  dinitrogen.numberDensity = outData.n2nd;
  oxygen.numberDensity = outData.ond;
  averageMolecularWeight = outData.wtmol;

  greal flat = 15;
  greal dphij = 5.0;
  if (latitude >= 85.0) {
    dphij = -5.0;
  }
  if (abs(latitude) > flat) {

    // Jacchia values at current position+5 degrees lat
    inData.lat = latitude + dphij;
    met07(inData, outData);
    greal phn = outData.ph;
    //greal dhn = outData.dh;
    greal thn = outData.th;

    // Jacchia values at current position+5 degrees lon
    inData.lat = latitude;
    inData.lng = longitude + 5.0;
    met07(inData, outData);
    greal phe = outData.ph;
    //greal dhe = outData.dh;
    greal the = outData.th;

    // dp/dy and dp/dx for geostrophic wind
    greal dpy = phn - pressure;
    greal dpx = phe - pressure;

    // dt/dx, dt/dy and dt/dz for vertical wind
    greal dtx = the - temperature;
    greal dty = thn - temperature;
    if (dphij < 0.0) {
      dpy = -dpy;
      dty = -dty;
    }

    inData.z = height + 1.0;
    inData.lat = latitude;
    inData.lng = longitude;
    met07(inData, outData);
    //greal pb = outData.ph;
    //greal db = outData.dh;
    greal tb = outData.th;
    greal wtm = outData.wtmol;

    dtz = (tb - temperature) / 1000.0;
    dmdz = (wtm - averageMolecularWeight) / 1000.0;

    greal ugh, vgh, wgh;
    wind(pressure, density, temperature, height, gravity, latitude, totalRadius, 
         dpx, dpy, dtx, dty, dtz, ugh, vgh, wgh);

    // change notation for output
    ewWind = ugh;
    nsWind = vgh;
    verticalWind = wgh;
  }
  // Latitude is between -15 and 15 degrees
  else {
    // Jacchia values at +flat and +flat+5
    inData.lat = flat;
    met07(inData, outData);
    greal pjx = outData.ph;
    //greal djx = outData.dh;
    greal tjx = outData.th;

    inData.lat = flat + dphij;
    met07(inData, outData);
    greal phn = outData.ph;
    //greal dhn = outData.dh;
    greal thn = outData.th;

    // Jacchia values at +flat and lon+5
    inData.lat = flat;
    inData.lng = longitude + 5.0;
    met07(inData, outData);
    greal phe = outData.ph;
    //greal dhe = outData.dh;
    greal the = outData.th;

    // dp/py and dp/dx for geostrophic wind
    greal dpy = phn - pjx;
    greal dpx = phe - pjx;
    // dt/dx, dt/dy and dt/dz for vertical wind
    greal dtx = the - tjx;
    greal dty = thn - tjx;
    if (dphij < 0.0) {
      dpy = -dpy;
      dty = -dty;
    }

    inData.z = height + 1.0;
    inData.lat = flat;
    inData.lng = longitude;
    met07(inData, outData);
    //greal pb = outData.ph;
    //greal db = outData.dh;
    greal tb = outData.th;
    greal wtm = outData.wtmol;

    dtz = (tb - tjx) / 1000.0;
    dmdz = (wtm - averageMolecularWeight) / 1000.0;
    greal ugh, vgh, wgh;
    wind(pressure, density, temperature, height, gravity, flat, totalRadius, 
         dpx, dpy, dtx, dty, dtz, ugh, vgh, wgh);

    greal uplus = ugh;
    greal vplus = vgh;
    greal wplus = wgh;

    // Jacchia values at -flat and -flat+5
    inData.z = height;
    inData.lat = -flat;
    inData.lng = longitude;
    met07(inData, outData);
    pjx = outData.ph;
    //djx = outData.dh;
    tjx = outData.th;

    inData.lat = -flat + dphij;
    met07(inData, outData);
    phn = outData.ph;
    //dhn = outData.dh;
    thn = outData.th;

    // Jacchia values at -fla and lon+5
    inData.lat = -flat;
    inData.lng = longitude + 5.0;
    met07(inData, outData);
    phe = outData.ph;
    //dhe = outData.dh;
    the = outData.th;

    // dp/dy and dp/dx for geostrophic wind
    dpy = phn - pjx;
    dpx = phe - pjx;
    // dt/dx, dt/dy and dt/dz for vertical wind
    dtx = the - tjx;
    dty = thn - tjx;
    if (dphij < 0.0) {
      dpy = -dpy;
      dty = -dty;
    }

    inData.z = height + 1.0;
    inData.lat = -flat;
    inData.lng = longitude;
    met07(inData, outData);
    //pb = outData.ph;
    //db = outData.dh;
    tb = outData.th;

    dtz = (tb - tjx) / 1000.0;
    wind(pressure, density, temperature, height, gravity, (-flat), totalRadius, 
         dpx, dpy, dtx, dty, dtz, ugh, vgh, wgh);

    // Change notation for output
    greal uminus = ugh;
    greal vminus = vgh;
    greal wminus = wgh;

    // Interpolate wind to latitude
    greal factint = 0.5 * (latitude / flat + 1.0);
    ewWind = uminus + factint * (uplus - uminus);
    nsWind = vminus + factint * (vplus - vminus);
    verticalWind = wminus + factint * (wplus - wminus);
  }
}


//! \brief met07 is the driving program for the MET model at a specified height.
//!
//! The atmospheric model is a modified Jacchia 1970 model and is given
//! in the subroutine J70.  All of the other subroutines were designed to
//! allow flexible use of the model so that various input parameters could
//! be varied within a driving program with very little software development.
//!
//! \param inData The input data structure.
//! \param[out] outData The output data structure.
void MET::met07(const InData &inData, OutData &outData)
{
  // Molecular weights
  const greal ei[6] = {28.0134, 31.9988, 15.9994, 39.948, 4.0026, 1.00797};

  // Set parameters to INDATA values
  greal z = inData.z;
  greal xlat = inData.lat;
  greal xlng = inData.lng;

  int i1 = inData.i1;
  greal f10 = inData.f10;
  greal f10b = inData.f10b;
  greal gi = inData.gi;

  greal gphi, rphi;
  getEffectiveRadius(xlat, gphi, rphi);

  greal dlon = xlng - longitude;
  greal sha = toRadians(solarHourAngle + dlon);
  greal sda = toRadians(solarDeclination);
  greal dd = dayOfYear;
  greal dy = dayOfYear / daysInYear;

  greal te;
  tinf(i1, f10, f10b, gi, toRadians(xlat), sda, sha, dy, te);

  greal a[6];
  greal tz, em, dens, dl;
  met_07_jac(z, te, tz, a[0], a[1], a[2], a[3], a[4], a[5], em, dens, dl, gphi, rphi);

  greal denlg = 0.0;
  slv(z, xlat, dd, denlg);

  // Fair helium number density between base fairing height (bfh) and 500 km

  if (z >= 500.0) {
    slvh(xlat, sda, dl, a[4]);
  }
  else if (z > bfh) {
    greal dhel1 = a[4];
    greal dhel2 = a[4];
    greal dlg1 = dl;
    greal dlg2 = dl;
    slvh(xlat, sda, dlg2, dhel2);

    fair5(dhel1, dhel2, dlg1, dlg2, z, a[4], dl);
  }

  dl = dl + denlg;
  dens = pow(10, dl);

  // Fill out outdata array
  outData.tex = te;
  outData.th = tz;
  for (int i = 0; i < 6; i++) {
    outData.gasnd[i] = 1.0e+6 * pow(10, a[i]);
  }
  outData.wtmol = em;
  outData.dh = dens * 1000.0;
  outData.unk = dl + 3.0;
  outData.ph = outData.dh *  UNIVERSAL_GAS * 1.0e3 * tz / em;

  // Compute sumn = sum of number densities and sumn = sum of product of
  // number density and molecular weight
  greal sumn = 0.0;
  greal summn = 0.0;
  for (int i = 0; i < 6; i++) {
    sumn += outData.gasnd[i];
    summn += ei[i] * outData.gasnd[i];
  }

  if (z < 170.0) {
    // Adjust number densities to be consistent with mass density
    for (int i = 0; i < 6; i++) {
      greal avx = AVOGADRO * 1.0e3 * outData.dh / summn;
      outData.gasnd[i] *= avx;
    }
  }
  else if (z > bfh) {
    // Adjust molecular weight to be consistent with number densities
    outData.wtmol = summn / sumn;
    // Re-compute pressure with revised molecular weight
    outData.ph = outData.dh *  UNIVERSAL_GAS * 1.0e3 * outData.th / outData.wtmol;
  }
}


//! \brief Calculates the molecular weight for altitudes between 90 and 105 km
//!
//! This function calculates the molecular weight for altitudes between 90 and 105
//! km according to equation (1) of SAO report 313.  Otherwise, molwt is 
//! set to unity.
//!
//! \param z Height in km.
//! \returns Mean molecular weight.
greal MET::molwt(greal z)
{
  // Coefficients for equation (1)
  const greal c[7] = { 28.15204, -0.085586, 1.284e-04, -1.0056e-05, -1.021e-05, 1.5044e-06, 9.9826e-08 };

  // Return 1 when height is not within bounds.
  greal molwt_out = 1.0;

	if (z <= 105.0) {
    greal u = z - 100.0;
    // molwt_out = c0 + c1*u + c2*u^2 + ... + c6*u^6
    molwt_out = c[6];
    for (int i = 5; i >= 0; --i) {
      molwt_out = molwt_out * u + c[i];
    }
	}

	return molwt_out;
}

//! \brief Calculates the temperature at altitude 
//!
//! Calculates the temperature at altitude z using equation (10) for 
//! altitudes between 90 and 125 km and equation (13) for altitudes greater
//! than 125 km, from SAO Report 313.
//!
//! \param z   Height \units{km}.
//! \param tx  See SAO Report 313.
//! \param t1  See SAO Report 313.
//! \param t3  See SAO Report 313.
//! \param t4  See SAO Report 313.
//! \param A   See SAO Report 313.
//!
//! returns The temperature at altitude z \units{K}.
greal MET::temp(greal z, greal tx, greal t1, greal t3, greal t4, greal A)
{
  const greal B = 4.5e-6;
	greal u = z - 125.0;
  greal temp_e10;
  if (u > 0.0){
		temp_e10 = tx + A * atan(t1 * u * (1.0 + B * pow(u, 2.5)) / A);
	}
	else {
    // Equation 10 (with c_i = ti)
		temp_e10 = tx + t1 * u + t3 * pow(u, 3.0) + t4 * pow(u, 4.0);
	}

	return temp_e10;
}

//! \brief Performs Gaussian Quadrature (integration).
//!
//! Subdivide total integration-altitude range into intervals suitable for 
//! applying Gaussian Quadrature, set the number of points for integration
//! for each sub-interval, and then perform Gaussian Quadrature.
//!
//! \param nmin   Undocumented.
//! \param z2  Height \units{km}.
//! \param tx  See SAO Report 313.
//! \param t1  See SAO Report 313.
//! \param t3  See SAO Report 313.
//! \param t4  See SAO Report 313.
//! \param a2  See SAO Report 313.
//! \param gphi  Undocumented.
//! \param rphi  Undocumented.
//! \param[out] r  The result of the integration.
void MET::gauss(int nmin, greal z2, greal tx, greal t1, greal t3, greal t4, greal a2, greal gphi, greal rphi, greal &r)
{
	const int ng[8] = { 4, 5, 6, 6, 6, 6, 6, 6 };
	const greal altmin[9] = { 90.0, 105.0, 125.0, 160.0, 200.0, 300.0, 500.0, 1500.0, 2500.0 };

	//coefficients for gaussian quadrature
	const greal c[8][6] = {
		{ 0.5555556, 0.3478548, 0.2369269, 0.1713245, 0.1294850, 0.1012285 },
		{ 0.8888889, 0.6521452, 0.4786287, 0.3607616, 0.2797054, 0.2223810 },
		{ 0.5555556, 0.6521452, 0.5688889, 0.4679139, 0.3818301, 0.3137067 },
		{ 0.0000000, 0.3478548, 0.4786287, 0.4679139, 0.4179592, 0.3626838 },
		{ 0.0000000, 0.0000000, 0.2369269, 0.3607616, 0.3818301, 0.3626838 },
		{ 0.0000000, 0.0000000, 0.0000000, 0.1713245, 0.2797054, 0.3137067 },
		{ 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.1294850, 0.2223810 },
		{ 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.1012285 }
	};

	//abscissas for gaussian quadrature
	const greal x[8][6] =
	{
		{ -0.7745967, -0.8611363, -0.9061798, -0.9324695, -0.9491079, -0.9602899 },
		{ 0.0000000, -0.3399810, -0.5384693, -0.6612094, -0.7415312, -0.7966665 },
		{ 0.7745967, 0.3399810, 0.0000000, -0.2386192, -0.4058452, -0.5255324 },
		{ 0.0000000, 0.8611363, 0.5384693, 0.2386192, 0.0000000, -0.1834346 },
		{ 0.0000000, 0.0000000, 0.9061798, 0.6612094, 0.4058452, 0.1834346 },
		{ 0.0000000, 0.0000000, 0.0000000, 0.9324695, 0.7415312, 0.5255324 },
		{ 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.9491079, 0.7966665 },
		{ 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.9602899 }
	};

	r = 0.0;

	for (int k = nmin - 1; k < 8; k++){
		int ngauss = ng[k];
		greal a = altmin[k];
		greal d = fmin(z2, altmin[k + 1]);
		greal rr = 0.0;
		greal del = 0.5*(d - a);
		int j = ngauss - 3;

		for (int i = 0; i < ngauss; i++){
			greal z = del * (x[i][j] + 1.0) + a;
			//calculation includes using molwt function and temp function
			rr += c[i][j] * molwt(z) *(gphi / (pow((1.0 + z / rphi), 2.0)) / temp(z, tx, t1, t3, t4, a2));
		}

		rr = del * rr;
		r += rr;

		if (d == z2) {
			break;
		}
	}
}

//! \brief Computes a number of metrics given height and exospheric temperature.
//!
//! jac calculates the temperature tz, the total density dens, and its 
//! logarithm dl, the mean molecular weight em, the individual species number
//! densities for n, o2, o, a, he, h at altitude z given the exospheric 
//! temperature t.  
void MET::met_07_jac(greal z, greal t, greal &tz, greal &an, greal &ao2, greal &ao, greal &aa, greal &ahe, greal &ah, greal &em, greal &dens, greal &dl, greal gphi, greal rphi)
{
  constexpr greal qn = 0.7811;
  constexpr greal qo2 = 0.20955;
  constexpr greal qa = 0.009343;
  constexpr greal qhe = 1.289e-05;
  constexpr greal t0 = 183.0;
  const greal alpha[6] = { 0.0, 0.0, 0.0, 0.0, -0.380, 0.0 };
  const greal ei[6] = { 28.0134, 31.9988, 15.9994, 39.948, 4.0026, 1.00797 };
  
  greal tx = 444.3807 + 0.02385 * t - 392.8292 * exp(-0.0021357 * t);
  greal a2 = 2.0 * (t - tx) / PI;
  greal tx_t0 = tx - t0;
  // Coefficients for equation (10) based on conditions (11)
  greal t1 = 1.9 * tx_t0 / 35.0;
  // t2 = 0
  greal t3 = -1.7 * tx_t0 / pow(35.0, 3);
  greal t4 = -0.8 * tx_t0 / pow(35.0, 4);
  tz = temp(z, tx, t1, t3, t4, a2);
	
	// Section 1

  greal d = min(z, 105.0);

  // Integrate gM/T from 90 to minimum of Z or 105 km

  greal r;
  gauss(1, d, tx, t1, t3, t4, a2, gphi, rphi, r);

  // The number 2.1926e-08 = density x temperature / mean molecular weight
  // at 90 km.

  em = molwt(d);
  greal td = temp(d, tx, t1, t3, t4, a2);
  dens = 2.1926e-08 * em * exp(-r / UNIVERSAL_GAS) / td;
  greal factor = AVOGADRO * dens;
  greal par = factor / (em);
  factor = factor / dryAirMolWt;

  // for altitudes below and at 105 calculate the individual specie number
  // densities from the mean molecular weight and total density

  if (z <= 105.0) {

    dl = log10(dens);
    an = log10(qn * factor);
    aa = log10(qa * factor);
    ahe = log10(qhe * factor);
    ao = log10(2.0 * par * (1.0 - em / dryAirMolWt));
    ao2 = log10(par * (em * ((1.0 + qo2) / dryAirMolWt) - 1.0));
    ah = 0.0;

    return;
  }

  // Section 2: This section is only for altitudes above 105 km

  // Note that having reached this section means that D in section 1 is
  // 105 km.

  // Calculate individual specie number densities from the total and mean
  // molecular weight at 105 km.

  greal di[6];
  di[0] = qn * factor;
  di[1] = par * (em * (1.0 + qo2) / 28.96 - 1.0);
  di[2] = 2.0 * par * (1 - em / 28.96);
  di[3] = qa * factor;
  di[4] = qhe * factor;

  // Integrate g/T from 105 km z km:

  gauss(2, z, tx, t1, t3, t4, a2, gphi, rphi, r);

  greal dit[6];
  for (int i = 0; i < 5; i++) {
    dit[i] = di[i] * pow((td / (tz)), (1.0 + alpha[i])) * exp((-ei[i] * r / UNIVERSAL_GAS));

    if (dit[i] <= 0.0) {
      dit[i] = 1.0e-06;
    }
  }

  // This section calculates atomic hydrogen densities above 500 km altitude.
  // Below this altitude, H densities are set to 10**6.

  // Section 3

  if (z > 500.0) {
    greal a1 = 500.0;
    greal s = temp(a1, tx, t1, t3, t4, a2);

    di[5] = pow(10.0, (73.13 - 39.4 * log10(s) + 5.5 * log10(s) * log10(s)));

    gauss(7, z, tx, t1, t3, t4, a2, gphi, rphi, r);

    dit[5] = di[5] * (s / (tz)) * exp((-ei[5] * r / UNIVERSAL_GAS));
  }
  else {
    dit[5] = 1.0;
  }

  // For altitudes greater than 105 km, calculate total density and mean
  // molecular weight from individual specie number densities.
  dens = 0;
  for (int i = 0; i < 6; i++) {
    dens = dens + ei[i] * dit[i] / AVOGADRO;
  }

  em = dens * AVOGADRO / (dit[0] + dit[1] + dit[2] + dit[3] + dit[4] + dit[5]);
  dl = log10(dens);
  an = log10(dit[0]);
  ao2 = log10(dit[1]);
  ao = log10(dit[2]);
  aa = log10(dit[3]);
  ahe = log10(dit[4]);
  ah = log10(dit[5]);
}

//! \brief Computes the seasonal-latitudinal variation of density in the lower thermosphere.
//!
//! This function computes the seasonal-latitudinal variation of density in the lower
//! thermosphere in accordance with L. Jacchia, SAO 332, 1971.  This affects
//! the densities between 90 and 170 km.  This member function need not be 
//! called for densities above 170 km, because no effect is observed.
//! 
//! The variation should be computed after the calculation of density due 
//! to temperature variations and the density must be in the form of a 
//! base 10 log.  No adjustments are made to the temperature or consituent
//! number densities in the region affected by this variation.
//!
//! \param alt Height in km.
//! \param lat Latitude in degrees.
//! \param dayOfTheYear Day of the year. ???
//! \param[out] deltaLog10Density A greal reference.
//! \retval deltaLog10Density The seasonal-latitudinal variation of density (log10).
void MET::slv(greal alt, greal lat, greal dayOfTheYear, greal &deltaLog10Density)
{
	deltaLog10Density = 0.0;

	// The lower thermosphere is defined as below 170 km
  if (alt <= 170.0) {

    // Compute density change in lower thermosphere  (equation 24)
    greal z = alt - 90.0;
    greal x = -0.0013 * z * z;
    greal S = 0.014 * z * exp(x);                  // height factor
    greal P = sin(0.0172 * dayOfTheYear + 1.72);      // day factor  (0.0172 used for 2*PI/365.2422)
    greal sinPhi = pow(sin(toRadians(lat)), 2.0);  // latitude factor
    deltaLog10Density = S * P * sinPhi;

    // The sign of the variation agrees with the sign of the latitude.
    if (lat < 0.0) {
      deltaLog10Density = -deltaLog10Density;
    }
  }
}

//! \brief Computes the seasonal-latitudinal variation of the helium number density.
//! 
//! This function computes the seasonal-latitudinal variation of the helium number
//! density according to L. Jacchia, SAO 332, 1971.  This correction is not
//! important below about 500 km.
//!
//! \param lat Latitude in degrees.
//! \param sda The solar declination angle in radians.
//! \param log10Density The density (log10).
//! \param log10HeND The helium number density (log10).
//! \retval log10Density The seasonal-latitudinal variation of density (log10).
//! \retval log10HeND The helium number density (log10).
void MET::slvh(greal lat, greal sda, greal &log10Density, greal &log10HeND)
{
  // Save the Helium number density before adjusting.
	greal initialHeND = pow(10.0, log10HeND);

  greal a = abs(0.65 * (sda / 0.40909079));  // (equation 25) (0.40909079 is used for 23.44 * PI / 180) 
  greal b = 0.5 * toRadians(lat);

	// The sign of b agrees with the sign of sda.  (this is faster than sda/abs(sda))
	if (sda < 0.0) b = -b;

	// Compute variation deltaLog10HeND (equation 25)
  greal y = pow(sin(0.7854 - b), 3);               // .7854 used for PI/4
  greal deltaLog10HeND = a * (y - 0.35356);        // .35356 used for (sin(PI/4))^3

	//Compute helium number density change, deltaHeND
  log10HeND = log10HeND + deltaLog10HeND;
  greal adjustedHeND = pow(10.0, log10HeND);
  greal deltaHeND = adjustedHeND - initialHeND;

  // Add the variation into the density (log10)
  greal rho = pow(10.0, log10Density);
  greal deltaRho = 6.646e-24 * deltaHeND;  // 6.646e-24 is atmoic mass of helium in grams
	rho = rho + deltaRho;
	log10Density = log10(rho);
}


//! \brief calculates the exospheric temperature.
//!
//! tinf calculates the exospheric temperature according to L. Jacchia SAO 313, 1970.
void MET::tinf(int i1, greal f10, greal f10b, greal gi, greal xlat, greal sda, greal sha, greal dy, greal &te)
{
  constexpr greal beta = -0.6457718;    // (Equation 17)
  constexpr greal gamma = 0.7504916;    // (Equation 17)
  constexpr greal p = 0.1047198;        // (Equation 17)
  constexpr greal xm = 2.5, xnn = 3.0;  // (Equation 17)
  constexpr greal R = 0.31; // (Equation 18, value in comments above 18)
  constexpr greal c1 = 383.0, c2 = 3.32, c3 = 1.8;  // (equation 14)
  constexpr greal d1 = 28.0, d2 = 0.03;             // (equation 21)
  constexpr greal d3 = 1.0, d4 = 100.0, d5 = -0.08; // (equation 22)
  // Ei are semiannual variation variables
  constexpr greal e1 = 2.41, e2 = 0.349, e3 = 0.206, e4 = TWO_PI;                 // (equation 23)
  constexpr greal e5 = 3.9531708, e6 = 2.0 * TWO_PI, e7 = 4.3214352, e8 = 0.1145; // (equation 23)
  constexpr greal e9 = 0.5, e10 = TWO_PI, e11 = 5.974262, e12 = 2.16;             // (equation 23)

  // solar activity variation  (equation 14)
  greal Tc = c1 + c2 * f10b + c3 * (f10 - f10b);

  // diurnal variation
  greal eta = 0.5 * abs(xlat - sda);             // (equation 15)
  greal theta = 0.5 * abs(xlat + sda);           // (equation 15)
  greal tau = sha + beta + p * sin(sha + gamma); // (equation 16)

  if (tau > PI)
    tau = tau - TWO_PI;
  if (tau < (-PI))
    tau = tau + TWO_PI;

  // (Equation 17)
  greal sinThetaM = pow(sin(theta), xm);
  greal cosEtaM = pow(cos(eta), xm);
  greal cosTauN = pow(cos(tau / 2.0), xnn);
  greal denom = 1.0 + R * sinThetaM;             // denominator in eq 17
  greal numer = (cosEtaM - sinThetaM);           // numerator in eq 17
  greal Tl = Tc * (denom + R * numer * cosTauN); // denom has been distributed in eq 17

	// geomagnetic variation
  greal deltaTg;
  if (i1 == 1) {
    // (Equation 21)
    deltaTg = d1 * gi + d2 * exp(gi);
  }
  else {
    // (Equation 22)
    deltaTg = d3 * gi + d4 * (1.0 - exp(d5 * gi));
  }

  // semiannual variation  (Equation 23)
  greal g3 = 0.5 * (1.0 + sin(e10 * dy + e11));
  g3 = pow(g3, e12);
  greal tau1 = dy + e8 * (g3 - e9);
  greal g1 = e2 + e3 * sin(e4 * tau1 + e5);
  greal g2 = sin(e6 * tau1 + e7);
  greal deltaTs = e1 + f10b * g1 * g2;

  // exospheric temperature
  te = Tl + deltaTg + deltaTs;
}

//! \brief Computes u, v, and w winds.
//!
//! Computes winds for MET (Jacchia) section (geostrophic winds with 
//! viscosity modifications at high altitude)
void MET::wind(greal ph, greal dh, greal th, greal h, greal g, greal phid, greal ri, greal dpx, greal dpy, greal dtx, greal dty, greal dtz, greal &ugh, greal &vgh, greal &wgh)
{
	constexpr greal beta = 1.458e-6, sval = 110.4;

	// Avoid cases with 0 temperature or density
	if (dh <= 0.0 || th <= 0.0) {
		ugh = 0.0;
		vgh = 0.0;
		wgh = 0.0;
	}
	else {
    // distances for 5 degrees separation points
    greal dy = 5000.0 * ri * TO_RADIANS;
    greal dx = dy * cos(toRadians(phid));

    // Coriolis factor
    greal coriol = TO_RADIANS * sin(toRadians(phid)) / 120.0;

    // Viscosity coefficients
    greal visc = beta * pow(th, 1.5) / (th + sval);
    greal vls = 5.3 + 0.0622 * pow(h, 1.5);
    vls = min(h, vls);
    greal viscfac = visc / (1.0e6 * dh * pow(vls, 2));

    // Horizontal pressure gradients
    greal dpdx = dpx / (dx * dh);
    greal dpdy = dpy / (dy * dh);
    greal denom = pow(coriol, 2) + pow(viscfac, 2);

    // Viscosity-modified geostrophic winds (horizontal components)
    ugh = -((coriol * dpdy) + (viscfac * dpdx)) / denom;
    vgh = ((coriol * dpdx) - (viscfac * dpdy)) / denom;

    // Insure wind component magnitudes < 0.7 times sound speed
    greal splim = 0.7 * sqrt(1.4 * ph / dh);
    if (abs(ugh) > splim) {
      ugh = copysign(splim, ugh);
    }

    if (abs(vgh) > splim) {
      vgh = copysign(splim, vgh);
    }

    greal cp = 7.0 * ph / (2.0 * dh * th);

    // Mean vertical wind from Montgomery stream function
    wgh = -cp * (ugh * dtx / dx + vgh * dty / dy) / (g + cp * dtz);
  }
}

//! \brief Fairs between the region above 500 km and the region below 440 km.
//!
//! This subroutine fairs between the region above 500 km, which invokes
//! the seasonal-latitudinal variation of the helium number density (slvh),
//! and the region below 440 km, which does not invoke any seasonal-latitudinal 
//! variation at all.
void MET::fair5(greal dhel1, greal dhel2, greal dlg1, greal dlg2, greal h, greal &fdhel, greal &fdlg)
{
  // height fairing factor
  greal hi = TO_RADIANS * 1.50 * (h - bfh);

  // Non-slvh fairing coefficient
  // greal czi = pow(cos(hi), 2);

  // slvh fairing factor
  Interpolator interp(square(sin(hi)));

  // Faired density
  fdlg = interp.linear(dlg1, dlg2);

  // Faired helium number density
  fdhel = interp.linear(dhel1, dhel2);
}


} // namespace
