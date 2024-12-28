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
#include "Atmosphere.h"
#include "EarthCommon.h"
#include "EarthInputParameters.h"

namespace GRAM {

//! \brief The Marshall Engineering Thermosphere model (MET) for the Earth themosphere.
//!
//! The Marshall Engineering Thermosphere Model(MET) is essentially a modified 
//! Jacchia 1970 model that includes some spatial and temporal variation patterns 
//! of the Jacchia 1971 model.In addition to thermospheric densities and temperatures
//! the well-documented code provides also several often used parameters like 
//! gravitational acceleration and specific heat. MET was developed at NASA's Marshall 
//! Space Flight Center in Huntsville primarily for engineering applications. 
//! The MSIS model is generally considered superior to MET because of its larger data 
//! base and its more elaborate mathematical formalism.
//!
//! \ingroup EarthGRAM
class MET : public Atmosphere, public EarthCommon
{
public:
	MET();
  MET(const MET& orig) = default;
  virtual ~MET() override = default;

  void setInputParameters(const EarthInputParameters& params);
  void setDayOfYear(greal doy) { dayOfYear = doy; }

  void update() override;
	
  greal getTemperatureGradient() const { return dtz; }
  greal getMolecularWeightGradient() const { return dmdz; }

private:
  //! \brief Input structure for met07().
  struct InData {      
    greal z;           //!< Height
    greal lat;         //!< Latitude
    greal lng;         //!< Longitude
    int i1;            //!< geomagnetic index
    greal f10;         //!< solar radio noise flux
    greal f10b;        //!< 162-day average F10
    greal gi;          //!< geomagnetic activity index
  };

  //! \brief Output structure for met07().
  struct OutData {     
    greal tex;         //!< Exospheric temperature
    greal th;          //!< Temperature
    union {            //!< Allows dual use (array or struct)
      greal gasnd[6];  //!< Treats the gases below as an array.
      struct {         //!< Or use the struct.
        greal n2nd;    //!< N2 number density  
        greal o2nd;    //!< O2 number density  
        greal ond;     //!< O number density   
        greal arnd;    //!< A number density   
        greal hend;    //!< He number density  
        greal hnd;     //!< H number density   
      };               
    };                 
    greal wtmol;       //!< Average molecular wt.
    greal dh;          //!< Total mass density   
    greal unk;         //!< Log10 mass density   
    greal ph;          //!< Total pressure       
  };

  greal molwt(greal z);
  void gauss(int nmin, greal z2, greal tx, greal t1, greal t3, greal t4, greal a2, greal gphi, greal rphi, greal &r);
	greal temp(greal alt, greal tx, greal t1, greal t3, greal t4, greal A);
	void met_07_jac(greal z, greal t, greal &tz, greal &an, greal &ao2, greal &ao, greal &aa, greal &ahe, greal &ah, greal &em, greal &dens, greal &dl, greal gphi, greal rphi);
	void slv(greal alt, greal lat, greal day, greal &deltaLog10Density);
	void slvh(greal lat, greal sda, greal &log10Density, greal &log10HeND);
	void met07(const InData &inData, OutData& outData);
	void tinf(int i1, greal f10, greal f10b, greal gi, greal xlat, greal sda, greal sha, greal dy, greal &te);
	void wind(greal ph, greal dh, greal th, greal h, greal g, greal phid, greal ri, greal dpx, greal dpy, greal dtx, greal dty, greal dtz, greal &ugh, greal &vgh, greal &wgh);
	void fair5(greal dhel1, greal dhel2, greal dlg1, greal dlg2, greal h, greal &fdhel, greal &fdlg);

  // outputs
  greal dtz = 0.0;   //!< Temperature gradient.
  greal dmdz = 0.0;  //!< Molecular Weight gradient.

  // Input Parameters
  int year = 0;              //!< Input start date.
  int month = 0;             //!< Input start date.
  int day = 0;               //!< Input start date.
  int hour = 0;              //!< Input start date.
  int minute = 0;            //!< Input start date.
  greal seconds = 0.0;       //!< Input start date.
  greal f10 = 0.0;           //!< Daily solar flux.
  greal f10b = 0.0;          //!< Mean solar flux.
  greal ap = 0.0;            //!< Geomagnetic index.
  greal dayOfYear = 0.0;     //!< Day of year.
  greal daysInYear = 365.0;  //!< Days in a year.

  const greal bfh = 440.0;    //!< base fairing height \units{km}.

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(MET, molwt);
  FRIEND_TEST(MET, temp);
  FRIEND_TEST(MET, gauss);
  FRIEND_TEST(MET, met_07_jac);
  FRIEND_TEST(MET, slv);
  FRIEND_TEST(MET, slvh);
  FRIEND_TEST(MET, met07);
  FRIEND_TEST(MET, tinf);
  FRIEND_TEST(MET, wind);
  FRIEND_TEST(MET, fair5);
  FRIEND_TEST(MET, update);
#endif // GRAM_UNIT_TEST

};

} // namespace
