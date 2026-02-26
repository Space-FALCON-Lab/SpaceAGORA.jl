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
#include "HWM.h"
#include "Atmosphere.h"
#include "EarthCommon.h"
#include "EarthInputParameters.h"

namespace GRAM {

//! \brief The JB2008 empirical model of the Earth thermosphere.
//!
//! The Jacchia-Bowman 2008 Model Atmosphere.
//! This is the CIRA "Integration Form" of a Jacchia Model.
//! There are no tabular values of density. Instead, the barometric
//! equation and diffusion equation are integrated numerically using
//! the Newton-Coates method to produce the density profile up to the
//! input position.
//!
//! \ingroup EarthGRAM
class JB2008 : public Atmosphere, public EarthCommon
{
public:
  JB2008();
  JB2008(const JB2008& orig) = default;
  virtual ~JB2008() override = default;

  void setInputParameters(const EarthInputParameters& params);
  void setDayOfYear(greal doy) { dayOfYear = doy; }
  void setJulianDay(greal jDay) { xmjd = jDay - 2400000.5; }

  void update() override;

  greal getTemperatureGradient() const { return dtz; }
  greal getMolecularWeightGradient() const { return dmdz; }

private:
  greal xambar(greal z);
	void semian08(int iyr, greal day, greal ht, greal f10b, greal s10b, greal xm10b, 
                greal &fzz, greal &gtz, greal &drlog);
	void dtsub(greal f10, greal xlst, greal xlat, greal zht, greal &dtc);
	void tmoutd(greal d1950, int &iyr, greal &day);
	greal xlocal(greal z, const greal tc[]);
	greal xgrav(greal z);
  void JB08(greal amjd, greal sun[], greal sat[], greal f10, greal f10b, greal s10, greal s10b,
            greal xm10, greal xm10b, greal y10, greal y10b, greal dstdtc, greal temp[], greal &rho,
            greal &pres, greal &avgmw, greal and1[], greal &sumn);
  void JB08HWM(greal xmin, greal z, greal xlat, greal xlon, greal &ph, greal &dh, greal &th,
               greal &n2nd, greal &o2nd, greal &ond, greal &arnd, greal &hend, greal &hnd,
               greal &nnd, greal &wtmol, greal &tex, greal &uh, greal &vh);

  HWM hwm;               //!< The harmonic winds model.
  greal dtz = 0.0;       //!< Temperature gradient.
  greal dmdz = 0.0;      //!< Molecular weight gradient.
  greal dayOfYear = 0.0; //!< Day of year.
  greal xmjd = 0.0;      //!< Julian day.

  // Input Parameters
  int year = 0;          //!< Input start date.
  int month = 0;         //!< Input start date.
  int day = 0;           //!< Input start date.
  int hour = 0;          //!< Input start date.
  int minute = 0;        //!< Input start date.
  greal seconds = 0.0;   //!< Input start date.
  greal f10 = 0.0;       //!< Daily 10.7-cm flux \units{sfu}.
  greal f10b = 0.0;      //!< Mean 10.7-cm flux \units{sfu}.
  greal ap = 0.0;        //!< Geomagnetic index.
  greal s10 = 0.0;       //!< EUV index (26-34 nm) scaled to F10 units (0.0 -> dailyS10 = dailyF10).
  greal s10b = 0.0;      //!< EUV 81-day center-averaged index (0.0 -> meanS10 = meanF10).
  greal xm10 = 0.0;      //!< MG2 index scaled to F10 units (0.0 -> dailyXM10 = dailyF10).
  greal xm10b = 0.0;     //!< MG2 81-day center-averaged index (0.0 -> meanXM10 = meanF10).
  greal y10 = 0.0;       //!< Solar X-Ray & Lya index scaled to F10 (0.0 -> dailyY10 = dailyF10).
  greal y10b = 0.0;      //!< Solar X-Ray & Lya 81-day avg. centered index (0.0 -> meanY10 = meanF10).
  greal dstdtc = 0.0;    //!< Temperature change computed from Dst index.

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(JB2008, update);
#endif // GRAM_UNIT_TEST

};

} // namespace
