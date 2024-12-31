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

#include "unittest_friend.h"
#include "gram.h"
#include "EarthCommon.h"
#include "Atmosphere.h"
#include "HWM.h"
#include "EarthInputParameters.h"

namespace GRAM {

//! \brief The model of a Earth upper atmosphere.
//!
//! The NRLMSIS-00 empirical atmosphere model was developed by Mike Picone, Alan Hedin, 
//! and Doug Drob based on the MSISE90 model. The main differences to MSISE90 involve
//! (1) the extensive use of drag and accelerometer data on total mass density, 
//! (2) the addition of a component to the total mass density that accounts for possibly 
//!     significant contributions of O + and hot oxygen at altitudes above 500 km, and 
//! (3) the inclusion of the SMM UV occultation data.
//!
//! The MSISE90 model describes the neutral temperature and densities in Earth's atmosphere 
//! from ground to thermospheric heights. Below 72.5 km the model is primarily based on the 
//! MAP Handbook (Labitzke et al., 1985) tabulation of zonal average temperature and pressure 
//! by Barnett and Corney, which was also used for the CIRA-86. Below 20 km these data were 
//! supplemented with averages from the National Meteorological Center (NMC). In addition, 
//! pitot tube, falling sphere, and grenade sounder rocket measurements from 1947 to 1972 
//! were taken into consideration. Above 72.5 km MSISE-90 is essentially a revised MSIS-86 
//! model taking into account data derived from space shuttle flights and newer incoherent scatter results.
//!
//! \ingroup EarthGRAM
class MSIS : public Atmosphere, public EarthCommon
{
public:
	MSIS();
  MSIS(const MSIS& orig) = default;
  virtual ~MSIS() = default;

  void setInputParameters(const EarthInputParameters& params);
  void setDayOfYear(greal doy) { dayOfYear = doy; }

  void update() override;

  greal getTemperatureGradient() const { return dtz; }
  greal getMolecularWeightGradient() const { return dmdz; }

private:
  inline greal getCswSw(int index) { return hwm.getCswSw(index); }
  inline greal getCswSwc(int index) { return hwm.getCswSwc(index); }

	void ghp7(int iyd, greal sec, greal &alt, greal glat, greal glong, greal stl, greal f107a,
            greal f107, greal ap[], greal d[], greal t[], greal press);
  void gtd7(int iyd, greal sec, greal alt, greal glat, greal glong, greal stl, greal f107a,
            greal f107, greal ap[], greal mass, greal d[], greal t[]);
  void gts7(int iyd, greal sec, greal alt, greal glat, greal glong, greal stl, greal f107a,
            greal f107, greal ap[], int mass, greal &dm28, greal d[], greal t[]);
  void msishwm(int iyr, greal z, greal xlat, greal xlon,
               greal f10b, greal f10, greal ap[], greal &ph, greal &dh, greal &th, greal &n2nd,
               greal &o2nd, greal &ond, greal &arnd, greal &hend, greal &hnd, greal &nnd,
               greal &wtmol, greal &tex, greal &uh, greal &vh);
  greal globe7(greal yrd, greal sec, greal lat, greal lon, greal tloc, greal f107a, greal f107,
               greal ap[], const greal p[]);
  greal densm(greal alt, greal d0, greal xm, greal &tz, int mn3, greal zn3[], greal tn3[],
              greal tgn3[2], int mn2, greal zn2[], greal tn2[], greal tgn2[2]);
  greal densu(greal alt, greal dlb, greal tinf, greal tlb, greal xm, greal alpha, greal &tz,
              greal zlb, greal s2, int mn1, const greal zn1[5], greal tn1[5], greal tgn1[2]);
	greal glob7s(const greal p[]);
	greal dnet(greal dd, greal dm, greal zhm, greal xmm, greal xm);
	greal ccor(greal alt, greal r, greal h1, greal zh);
	greal ccor2(greal alt, greal r, greal h1, greal zh, greal h2);
	greal zeta(greal zz, greal zl, greal re) { return (zz - zl)*(re + zl) / (re + zz); }
	greal scalh(greal alt, greal xm, greal temp);
  greal sumex(greal ex);
	int vtst7(int iyd, greal sec, greal glat, greal glong, greal stl, greal f107a, greal f107, greal ap[], int ic);
	void meters(bool meter);
	void splini(greal xa[], greal ya[], greal y2a[], int n, greal x, greal &yi);
	void glatf(greal lat, greal &gv, greal &refff);

  greal g0(greal a, const greal p[]);

  greal sg0(greal ex, const greal ap[], const greal p[]);


  greal plg[9][4] = { {0} };  //!< Legendre polynomials
  greal dd = 0.0;             //!< unknown
  greal ctloc = 0.0;          //!< Cosine of tloc
  greal stloc = 0.0;          //!< Sine of tloc
  greal c2tloc = 0.0;         //!< Cosine of 2 tloc
  greal s2tloc = 0.0;         //!< Sine of 2 tloc
  greal c3tloc = 0.0;         //!< Cosine of 3 tloc
  greal s3tloc = 0.0;         //!< Sine of 3 tloc
  greal day_yr = 0.0;         //!< unknown
  greal dfa = 0.0;            //!< unknown
  greal apdf = 0.0;           //!< unknown
  greal apt[4] = { 0.0 };     //!< unknown
  greal xlong = 0.0;          //!< longitude
  greal tn1[5] = { 0.0 };     //!< K Lower thermosphere temperature variations
  greal tn2[4] = { 0.0 };     //!< K Lower mesosphere/upper stratosphere temperature variations
  greal tn3[5] = { 0.0 };     //!< K Lower stratosphere and troposphere temperature variations
	greal tgn1[2] = { 0.0 };    //!< K Lower thermosphere temperature gradient
  greal tgn2[2] = { 0.0 };    //!< K Lower mesospherer/upper stratosphere temperature gradient
  greal tgn3[2] = { 0.0 };    //!< K Lower stratosphere and troposphere temperature gradient
  greal gsurf = 0.0;          //!< gravity.
  greal re = 0.0;             //!< effective radius.

	const static greal pavgm[10];     //!< Middle Atmosphere Averages
  const static greal ptm[10];       //!< Lower Boundary
  const static greal pdm[10][8];    //!< Model data.
  const static greal pdl[25][2];    //!< TURBO
  const static greal sam[100];      //!< SEMIANNUAL MULT SAM
  const static greal pt[150];       //!< TEMPERATURE
  const static greal ps[150];       //!< S PARAM
  const static greal pd[9][150];    //!< Gas DENSITY
  const static greal ptl[4][100];   //!< Model data.
  const static greal pma[10][100];  //!< Model data.

  int imr = 0;              //!< Conversion flag.
  greal dayOfYear = 0.0;    //!< Day of year.

  HWM hwm;                  //!< The harmonic winds model.
  greal dtz = 0;            //!< Temperature gradient.
  greal dmdz;               //!< Molecular weight gradient.

  // Input Parameters
  int year = 0;             //!< Input start date.
  int month = 0;            //!< Input start date.
  int day = 0;              //!< Input start date.
  int hour = 0;             //!< Input start date.
  int minute = 0;           //!< Input start date.
  greal seconds = 0.0;      //!< Input start date.
  greal f10 = 0.0;          //!< solar radio noise flux
  greal f10b = 0.0;         //!< 162-day average F10
  greal ap = 0.0;           //!< geomagnetic index

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(MSIS, update);
#endif // GRAM_UNIT_TEST
};

} // namespace
