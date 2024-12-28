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

#include <fstream>
#include <string>
#include "unittest_friend.h"
#include "gram.h"
#include "EarthCommon.h"
#include "Atmosphere.h"
#include "EarthInputParameters.h"

namespace GRAM {

// Predeclaration.  See definition below.
struct EarthPPM;

//! \brief The Middle Atmosphere Program (MAP) model of Earth's middle atmosphere.
//!
//! The MAP class computes atmospheric values for the middle atmosphere region.
//!
//! \ingroup EarthGRAM
class MAP : public Atmosphere, public EarthCommon
{
public:
  MAP();
  MAP(const MAP& orig) = default;
  virtual ~MAP() override = default;

  void setInputParameters(const EarthInputParameters& params);

  void update() override;  // was mapmod

  void getStandardDeviations(greal& pressureSD, greal& densitySD, greal& temperatureSD);  // was rterp
  void getWindStandardDeviations(greal& ewSD, greal& nsSD);
  double getTemperatureGradient() const { return temperatureGradient; }

  void getSpeciesConcentrations(int iyr, greal pres, greal ppmWaterLow, 
                                greal& waterSD, greal &oxygenND, EarthPPM& ppm);  // was concvals

private:
  void gterp(int ih, greal phi, greal& p, greal& d, greal& t, greal& u, greal& dpy, greal& dty);
  void pdtuv(int ih, greal& ps, greal& ds, greal& ts, greal& us,greal& vs, greal& dtx, greal& dty);
  void heightInterpolation(greal p1, greal d1, greal t1, greal z1, 
              greal p2, greal d2, greal t2, greal z2, 
              greal& p, greal& d, greal& t, greal z);    // was inter2
  double afglconc(EarthPPM& ppm);
  void mapconc(greal pres, greal& o3m, greal& h2om, greal& sh2om, greal& n2om, greal& ch4m, greal& oxynd);
  greal larcwat();

  template<typename T>
  size_t getUpperIndex(greal value, size_t size, T &array);

  // MapData
  void initializeData();
  void concinit();
  void mapinit();

  // Input Parameters
  greal ewWindPerturbationScale = 1.0;      //!< Scales East/West wind perturbations.
  greal nsWindPerturbationScale = 1.0;      //!< Scales North/South wind perturbations.
  greal densityPerturbationScale = 1.0;     //!< Scales density standard deviations.
  greal pressurePerturbationScale = 1.0;    //!< Scales pressure standard deviations.
  greal temperaturePerturbationScale = 1.0; //!< Scales temperature standard deviations.
  int month = 0;                            //!< The start time month.

  double temperatureGradient = 0.0;         //!< Temperature gradient with respect to height. \note was dtz
  static bool initialized;                  //!< Set to true when NCEP data has been initialized.
  static const size_t monthSize = 12;       //!< The number of model months.
  static const size_t zHeightSize = 21;     //!< The number of zdata height levels.
  static const size_t zLatSize = 19;        //!< The number of zdata latitude levels.

  void readZDataBlock(std::ifstream& zdataFile, greal zdata[monthSize][zHeightSize][zLatSize]);

  static greal pg[monthSize][zHeightSize][zLatSize];  //!< Pressure \units{Pa}.
  static greal dg[monthSize][zHeightSize][zLatSize];  //!< Density \units{kg/m^3}.
  static greal tg[monthSize][zHeightSize][zLatSize];  //!< Temperature \units{K}.
  static greal ug[monthSize][zHeightSize][zLatSize];  //!< E/W winds \units{m/s}.

  static const size_t sHeightSize = 15;     //!< The number of sdata height levels.
  static const size_t sLatSize = 19;        //!< The number of sdata latitude levels.
  static const size_t sLonSize = 18;        //!< The number of sdata longitude levels.

  void readSDataBlock(std::ifstream& sdataFile, greal sdata[monthSize][sHeightSize][sLatSize][sLonSize], greal scale);

  static greal psp[monthSize][sHeightSize][sLatSize][sLonSize];  //!< Pa    Stationary perturbation data for pressure
  static greal dsp[monthSize][sHeightSize][sLatSize][sLonSize];  //!< kg/M3 Stationary perturbation data for density
  static greal tsp[monthSize][sHeightSize][sLatSize][sLonSize];  //!< K     Stationary perturbation data for temperature
  static greal usp[monthSize][sHeightSize][sLatSize][sLonSize];  //!< M/s   Stationary perturbation data for Eastward wind
  static greal vsp[monthSize][sHeightSize][sLatSize][sLonSize];  //!< M/s   Stationary perturbation data for Northward wind

  static const size_t rHeightSize = 29;     //!< The number of rdata height levels.
  static const size_t rLatSize = 19;        //!< The number of rdata latitude levels.

  static greal pr[monthSize][rHeightSize][rLatSize];  //!< Pressure \units{Pa}.
  static greal dr[monthSize][rHeightSize][rLatSize];  //!< Density \units{kg/m^3}.
  static greal tr[monthSize][rHeightSize][rLatSize];  //!< Temperature \units{K}.
  static greal ur[monthSize][rHeightSize][rLatSize];  //!< E/W winds \units{m/s}.
  static greal vr[monthSize][rHeightSize][rLatSize];  //!< N/S winds \units{m/s}.

  static int rpData[monthSize][rHeightSize][rLatSize];  //!< Pressure \units{Pa}.
  static int rdData[monthSize][rHeightSize][rLatSize];  //!< Density \units{kg/m^3}.
  static int rtData[monthSize][rHeightSize][rLatSize];  //!< Temperature \units{K}.
  static int ruData[monthSize][rHeightSize][rLatSize];  //!< E/W winds \units{m/s}.
  static int rvData[monthSize][rHeightSize][rLatSize];  //!< N/S winds \units{m/s}.

  static const size_t lHeightSize = 35;      //!< The number of ldata height levels.
  static const size_t lLatSize = 10;         //!< The number of ldata latitude levels.

  static greal lData[4][lHeightSize][lLatSize];         //!< LaRC water data
  static greal hgtl[lHeightSize];                       //!< LaRC data height levels
  static greal h2ol[monthSize][lLatSize][lHeightSize];  //!< LaRC water data

  static const size_t aHeightSize = 50;     //!< The number of adata height levels.
  static const size_t aLatSize = 8;         //!< The number of adata latitude levels.

  static greal hgta[aHeightSize];                        //!< AFGL data height levels
  static greal aData[5][aHeightSize][5];                 //!< AFGL data
  static greal h2oa[monthSize][aLatSize][aHeightSize];   //!< AFGL data water
  static greal o3[monthSize][aLatSize][aHeightSize];     //!< AFGL data ozone
  static greal n2o[monthSize][aLatSize][aHeightSize];    //!< AFGL data nitrous-oxide
  static greal co[monthSize][aLatSize][aHeightSize];     //!< AFGL data carbon-monoxide
  static greal ch4[monthSize][aLatSize][aHeightSize];    //!< AFGL data methane
  static greal co2[aHeightSize];                         //!< AFGL data carbon-dioxide
  static greal o2[aHeightSize];                          //!< AFGL data dioxygen
  static greal n2[aHeightSize];                          //!< AFGL data dinitrogen

  static const size_t o3PressureSize = 24;     //!< The number of pressure levels.
  static const size_t o3LatSize = 19;          //!< The number of latitude levels.

  static greal po3map[o3PressureSize];                          //!< MAP data pressure levels.
  static greal o3map[monthSize][o3LatSize][o3PressureSize];     //!< MAP data ozone.
  static int o3Data[monthSize][o3PressureSize][o3LatSize - 2];  //!< MAP data ozone.

  static const size_t h2oPressureSize = 19;     //!< The number of pressure levels.
  static const size_t h2oLatSize = 8;           //!< The number of latitude levels.

  static greal ph2omap[h2oPressureSize];                          //!< MAP data pressure levels.
  static greal h2omap[monthSize][h2oLatSize][h2oPressureSize];    //!< MAP data water
  static greal sh2omap[monthSize][h2oLatSize][h2oPressureSize];   //!< MAP data water standard deviation
  static greal h2oData[monthSize][11][10];                        //!< MAP data water
  static greal h2oLowData[8][10];                                 //!< MAP data water

  static const size_t n2oPressureSize = 17;     //!< The number of pressure levels.
  static const size_t n2oLatSize = 19;          //!< The number of latitude levels.

  static greal pmap31[n2oPressureSize];                             //!< MAP data pressure levels.
  static greal n2omap[monthSize][n2oLatSize][n2oPressureSize];      //!< MAP data nitrous-oxide
  static int n2oData[monthSize][n2oPressureSize][n2oLatSize - 4];   //!< MAP data nitrous-oxide
  static greal ch4map[monthSize][n2oLatSize][n2oPressureSize];      //!< MAP data methane
  static int ch4Data[monthSize][n2oPressureSize][n2oLatSize - 4];   //!< MAP data methane

  static const size_t oxyPressureSize = 19;     //!< The number of pressure levels.
  static const size_t oxyLatSize = 19;          //!< The number of latitude levels.

  static greal mapzox[oxyPressureSize];                             //!< MAP data pressure levels.
  static greal oxymap[monthSize][oxyLatSize][oxyPressureSize];      //!< MAP data oxygen
  static int oxyData[monthSize][oxyPressureSize][oxyLatSize - 1];   //!< MAP data oxygen

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(MAP, initializeData);
  FRIEND_TEST(MAP, concinit);
  FRIEND_TEST(MAP, mapinit);
  FRIEND_TEST(MAP, mapmod);
  FRIEND_TEST(MAP, heightInterpolation);
  FRIEND_TEST(MAP, afglconc);
  FRIEND_TEST(MAP, mapconc);
  FRIEND_TEST(MAP, larcwat);
  FRIEND_TEST(MAP, initializeData2);
  FRIEND_TEST(MAP, getStandardDeviations);
#endif // GRAM_UNIT_TEST
};

//! \brief Part per million structure.
struct EarthPPM {
  double water = 0.0;
  double ozone = 0.0;
  double nitrousOxide = 0.0;
  double carbonMonoxide = 0.0;
  double methane = 0.0;
  double carbonDioxide = 0.0;
  double dioxygen = 0.0;
  double dinitrogen = 0.0;
  double argon = 0.0;
  double nitrogen = 0.0;
  double oxygen = 0.0;
  double helium = 0.0;
  double hydrogen = 0.0;
};

//! \brief Gets the upper index of bracket in an array containing a specified value.
//!
//! \param value     The search value to be bracketed.     
//! \param size      The size of the array.
//! \param array     The data array.    
//! 
//! \returns Upper index of bracket in array containing value.
template<typename T>
size_t MAP::getUpperIndex(greal value, size_t size, T &array)
{
  size_t index = 0;
  if (value < array[0] || value > array[size - 1]) {
    return index;
  }
  while (index < size - 1 && value >= array[index]) {
    ++index;
  }
  return index;
}


} // namespace
