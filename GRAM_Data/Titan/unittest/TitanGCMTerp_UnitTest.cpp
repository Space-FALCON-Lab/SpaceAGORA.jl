#include "unittest.h"
#include "gram.h"
#include "TitanGCMTerp.h"
#include "Position.h"
#include "AtmosphereState.h"

namespace GRAM {

TEST( TitanGCMTerp, update)
{
  // SETUP
  Position p;
  AtmosphereState s;
  EphemerisState ephem;
  TitanGCMTerp gcmTerp;

  // GIVEN (INPUTS)
  p.longitude = 56.0_deg;
  p.latitude = 35.0_deg;
  p.height = 100.0;
  gcmTerp.setPosition(p);
  ephem.solarTime = 3.6;
  ephem.longitudeSun = 2.34_deg;
  gcmTerp.setEphemerisState(ephem);
  gcmTerp.updateGCMTables();

  // RUN
  gcmTerp.update();
  s = gcmTerp.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(141.05419044560941, s.temperature);
  EXPECT_DOUBLE_EQ(1090.8445012408342, s.pressure);
  EXPECT_DOUBLE_EQ(33.201889556795507, s.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(0.025866993573267771, s.density);
  EXPECT_DOUBLE_EQ(30.388910330146881, s.densityScaleHeight);
  EXPECT_DOUBLE_EQ(1.1203587193342481e+22, s.argon.numberDensity);
  EXPECT_DOUBLE_EQ(27.808000000000000, s.averageMolecularWeight);
  EXPECT_DOUBLE_EQ(1.6805380790013726e+22, s.methane.numberDensity);
  EXPECT_DOUBLE_EQ(5.3217039168376799e+23, s.dinitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(5.601793596671242e+23, s.totalNumberDensity);
  EXPECT_DOUBLE_EQ(62.911223184556043, s.ewWind);
  EXPECT_DOUBLE_EQ(0.0, s.nsWind);

  // TEAR-DOWN
}

TEST( TitanGCMTerp, updateGCMTables)
{
  // SETUP
  greal T03mb, z03mb;
  Position p;
  EphemerisState ephem;
  TitanGCMTerp gcmTerp;

  // GIVEN (INPUTS)
  p.longitude = 0; // not used
  p.height = 0; // not used
  p.latitude = 35.0_deg;
  gcmTerp.setPosition(p);
  ephem.solarTime = 13.6;
  ephem.longitudeSun = 2.34_deg;
  gcmTerp.setEphemerisState(ephem);

  // RUN
  gcmTerp.updateGCMTables();
  T03mb = gcmTerp.getTemperatureAt03mb();
  z03mb = gcmTerp.getHeightAt03mb();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(178.78353035414324, T03mb);
  EXPECT_DOUBLE_EQ(258.18287643769372, z03mb);

  // TEAR-DOWN
}

} // namespace
