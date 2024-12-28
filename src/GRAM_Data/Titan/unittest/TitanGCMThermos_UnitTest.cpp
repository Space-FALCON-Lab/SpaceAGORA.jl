#include "unittest.h"
#include "gram.h"
#include "TitanGCMThermos.h"
#include "Position.h"
#include "AtmosphereState.h"

namespace GRAM {

TEST( TitanGCMThermos, update )
{
  // SETUP
  Position p;
  EphemerisState ephem;
  TitanGCMTerp gcmTerp;
  TitanGCMThermos gcmThermos(gcmTerp);

  // GIVEN (INPUTS)
  p.longitude = 56.0_deg;
  p.latitude = 35.0_deg;
  p.height = 500.0;
  gcmThermos.setPosition(p);
  ephem.solarTime = 13.6;
  ephem.longitudeSun = 2.34_deg;
  gcmTerp.setEphemerisState(ephem);
  gcmThermos.setExosphericTemperature(123.4);
  gcmTerp.setPosition(p);
  gcmTerp.setEphemerisState(ephem);
  gcmTerp.updateGCMTables();

  // RUN
  gcmThermos.update();
  AtmosphereState s = gcmThermos.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(123.40000000000000, s.temperature);
  EXPECT_DOUBLE_EQ(0.0080081992472871288, s.pressure);
  EXPECT_DOUBLE_EQ(36.250792089856631, s.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(2.3274868335460794e-07, s.density);
  EXPECT_DOUBLE_EQ(35.764068543352899, s.densityScaleHeight);
  EXPECT_DOUBLE_EQ(1.0225088684828785e+18, s.argon.numberDensity);
  EXPECT_DOUBLE_EQ(29.819592726797794, s.averageMolecularWeight);
  EXPECT_DOUBLE_EQ(3019094526034035.5, s.methane.numberDensity);
  EXPECT_DOUBLE_EQ(3.6748894628155412e+18, s.dinitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(4.7004174258244536e+18, s.totalNumberDensity);
  EXPECT_DOUBLE_EQ(104.74813301070095, s.ewWind);
  EXPECT_DOUBLE_EQ(0.0000000000000000, s.nsWind);

  // TEAR-DOWN
}

} // namespace
