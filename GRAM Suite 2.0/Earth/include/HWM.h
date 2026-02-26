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

namespace GRAM {

//! \brief The HWM93 Horizontal Wind Model.
//!
//! The HWM93 horizontal wind model extends the HWM90 wind model from 100 km to the ground
//! and is thus a companion in a sense for MSISE90. The thermosphere
//! portion is unchanged from HWM90. The extension is primarily based on 
//! MF / Meteor radar, rocket data, and gradient winds from  CIRA86.
//!
//! \ingroup EarthGRAM
class HWM
{
public:
	HWM();
  HWM(const HWM& orig) = default;
  virtual ~HWM() = default;


	void gws5(int iyd, greal sec, greal alt, greal glat, greal glon, greal stl,
		greal f107a, greal f107, const greal ap[2], greal w[2]);
	void spline(const greal x[], const greal y[], int n, greal yp1, greal ypn, greal y2[]);
	void splint(const greal xa[], const greal ya[], const greal y2[], int n, greal x, greal &y);
	void tselec(const greal sv[25]);

  inline greal getCswSw(int index) { return sw[index]; }
  inline greal getCswSwc(int index) { return swc[index]; }

private:

	void legpl1(greal c, greal s, int l, int m, greal plg[][20], int lmax);
  void vsphr1(greal c, greal s, int l, int m, greal bt[][20], greal bp[][20], int lmax);
  void glbw5s(int iyd, greal lat, greal lon, greal stl, const greal *pb, const greal *pc, greal ww[]);
  void glbw5e(greal yrd, greal sec, greal lat, greal lon, greal stl, greal f107a, greal f107,
              const greal ap[], const greal pb[], const greal pc[], greal ww[]);
  void glbw5m(greal yrd, greal sec, greal lat, greal lon, greal stl, greal f107a, greal f107,
              const greal ap[], const greal pb[], const greal pc[], greal ww[]);
  greal wprof(greal z, greal zl, greal s, greal uinf, greal ulb, greal ulbd, int mn1,
              const greal zn1[], const greal un1[], const greal ugn1[], int mn2, const greal zn2[],
              const greal un2[], const greal ugn2[]);

  greal clat = 0.0;    //!< Cosine of the latitude
  greal slat = 0.0;    //!< Sine of the latitude
  greal clong = 0.0;   //!< Cosine of longitude
  greal slong = 0.0;   //!< Sine of longitude
  greal c2long = 0.0;  //!< Cosine of 2 * longitude
  greal s2long = 0.0;  //!< Sine of 2 * longitude
  greal cstl = 0.0;    //!< Cosine of apparent local solar time
  greal sstl = 0.0;    //!< Sine of apparent local solar time
  greal c2stl = 0.0;   //!< Cosine of 2 * apparent local solar time
  greal s2stl = 0.0;   //!< Sine of 2 * apparent local solar time
  greal c3stl = 0.0;   //!< Cosine of 3 * apparent local solar time
  greal s3stl = 0.0;   //!< Sine of 3 * apparent local solar time
  greal sw[25] = {     //!< Main terms. Hard coded to 1.0.
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, };
  greal swc[25] = {    //!< Cross terms. Hard coded to 1.0.
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, };
  greal bt[20][20] = { { 0.0 } };
  greal bp[20][20] = { { 0.0 } };
  greal xvl = -999.0;
  greal lvl = -1.0;
  greal mvl = -1.0;
  greal tll = -999.0;
  greal nsvl = -1.0;
  greal xll = -999.0;
  greal ngvl = -1.0;
	bool callTselec = true;

  static const greal pwb[200];
  static const greal pwc[200];
  static const greal pwbl[150];
  static const greal pwcl[150];
  static const greal pwbld[150];
  static const greal pwcld[150];
  static const greal pb12[150];
  static const greal pc12[150];
  static const greal pb13[150];
  static const greal pc13[150];
  static const greal pb14[150];
  static const greal pc14[150];
  static const greal pb15[150];
  static const greal pc15[150];
  static const greal pb15d[150];
  static const greal pc15d[150];
  static const greal pwp[26][100];

};

} // namespace
