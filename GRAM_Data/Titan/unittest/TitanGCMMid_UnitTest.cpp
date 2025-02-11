#include "unittest.h"
#include "gram.h"
#include "TitanGCMMid.h"
#include "Position.h"
#include "AtmosphereState.h"

namespace GRAM {

TEST( TitanGCMMid, update )
{
  // SETUP
  Position p;
  AtmosphereState s;
  EphemerisState ephem;
  TitanGCMTerp gcmTerp;
  TitanGCMMid gcmMid(gcmTerp);

  // GIVEN (INPUTS)
  p.longitude = 56.0_deg;
  p.latitude = 35.0_deg;
  p.height = 500.0;
  gcmMid.setPosition(p);
  ephem.solarTime = 13.6;
  ephem.longitudeSun = 2.34_deg;
  gcmMid.setEphemerisState(ephem);
  gcmMid.setExosphericTemperature(123.4);
  gcmTerp.setPosition(p);
  gcmTerp.setEphemerisState(ephem);
  gcmTerp.updateGCMTables();

  // RUN
  gcmMid.update();
  s = gcmMid.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(93.578639097033857, s.temperature);
  EXPECT_DOUBLE_EQ(0.073258840037147777, s.pressure);
  EXPECT_DOUBLE_EQ(29.567455422272808, s.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(2.6166006372860525e-06, s.density);
  EXPECT_DOUBLE_EQ(32.696105736244739, s.densityScaleHeight);
  EXPECT_DOUBLE_EQ(1.1050739864191535e+18, s.argon.numberDensity);
  EXPECT_DOUBLE_EQ(27.789985170613466, s.averageMolecularWeight);
  EXPECT_DOUBLE_EQ(1.7659106378407217e+18, s.methane.numberDensity);
  EXPECT_DOUBLE_EQ(5.3831231133325804e+19, s.dinitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(5.6702215757585678e+19, s.totalNumberDensity);
  EXPECT_DOUBLE_EQ(104.74813301070095, s.ewWind);
  EXPECT_DOUBLE_EQ(0.0000000000000000, s.nsWind);

  // TEAR-DOWN
}

} // namespace
