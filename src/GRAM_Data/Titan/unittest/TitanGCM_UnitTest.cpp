#include "unittest.h"
#define GRAM_UNIT_TEST
#include "TitanGCM.h"

namespace GRAM {

TEST( TitanGCM, Jacchia )
{
  // SETUP
  TitanGCM gcmTitan;

  // GIVEN (INPUTS)
  greal latitude = 12.3_deg;
  greal solarTime = 23.4;

  // RUN
  greal Tex270 = gcmTitan.Jacchia(latitude, solarTime, -26.43_deg, 1.25, 3.7, 0.077,
          0.064, -5.0_deg, 0.0_deg, -25.0_deg, 147.5);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(148.19757989554054, Tex270);

  // TEAR-DOWN
}

TEST( TitanGCM, updateExosphericTemperature )
{
  // SETUP
  Position p;
  TitanGCM gcmTitan;

  // GIVEN (INPUTS)
  p.latitude = 12.3_deg;
  gcmTitan.setPosition(p);
  EphemerisState ephem;
  ephem.solarTime = 23.4;
  ephem.longitudeSun = 34.5_deg;
  gcmTitan.setEphemerisState(ephem);

  // RUN
  gcmTitan.updateExosphericTemperature();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(155.28047704575522, gcmTitan.getExosphericTemperature());

  // TEAR-DOWN
}

TEST( TitanGCM, update_terp )
{
  // SETUP
  Position p;
  TitanGCM gcmTitan;

  // GIVEN (INPUTS)
  p.height = 234.5;
  p.longitude = 56.7_deg;  // not used
  p.latitude = 12.3_deg;
  gcmTitan.setPosition(p);
  EphemerisState ephem;
  ephem.solarTime = 23.4;
  ephem.longitudeSun = 34.5_deg;
  gcmTitan.setEphemerisState(ephem);

  // RUN
  gcmTitan.update();
  AtmosphereState s = gcmTitan.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(160.84136571946357, s.temperature);
  EXPECT_DOUBLE_EQ(32.831392911518137, s.pressure);
  EXPECT_DOUBLE_EQ(42.184445242669661, s.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(0.00068272919725217334, s.density);
  EXPECT_DOUBLE_EQ(40.321971380716079, s.densityScaleHeight);
  EXPECT_DOUBLE_EQ(2.9570564778584533e+20, s.argon.numberDensity);
  EXPECT_DOUBLE_EQ(27.808000000000000, s.averageMolecularWeight);
  EXPECT_DOUBLE_EQ(4.4355847167876792e+20, s.methane.numberDensity);
  EXPECT_DOUBLE_EQ(1.4046018269827652e+22, s.dinitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(1.4785282389292266e+22, s.totalNumberDensity);
  EXPECT_DOUBLE_EQ(109.10174888309838, s.ewWind);
  EXPECT_DOUBLE_EQ(0.0000000000000000, s.nsWind);

  // TEAR-DOWN
}

TEST( TitanGCM, update_mid )
{
  // SETUP
  Position p;
  TitanGCM gcmTitan;

  // GIVEN (INPUTS)
  p.height = 834.5;
  p.longitude = 56.7_deg;  // not used
  p.latitude = 12.3_deg;
  gcmTitan.setPosition(p);
  EphemerisState ephem;
  ephem.solarTime = 23.4;
  ephem.longitudeSun = 34.5_deg;
  gcmTitan.setEphemerisState(ephem);

  // RUN
  gcmTitan.update();
  AtmosphereState s = gcmTitan.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(149.36077153487457, s.temperature);
  EXPECT_DOUBLE_EQ(2.5787019679449214e-05, s.pressure);
  EXPECT_DOUBLE_EQ(58.138509451939399, s.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(5.7396838442206511e-10, s.density);
  EXPECT_DOUBLE_EQ(55.022678272488662, s.densityScaleHeight);
  EXPECT_DOUBLE_EQ(200919189970555.06, s.argon.numberDensity);
  EXPECT_DOUBLE_EQ(27.64124196921145, s.averageMolecularWeight);
  EXPECT_DOUBLE_EQ(514497152766050.94, s.methane.numberDensity);
  EXPECT_DOUBLE_EQ(11789515261384344.0, s.dinitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(12504931604120950.0, s.totalNumberDensity);
  EXPECT_DOUBLE_EQ(109.32620929908805, s.ewWind);
  EXPECT_DOUBLE_EQ(0.0000000000000000, s.nsWind);

  // TEAR-DOWN
}

TEST( TitanGCM, update_thermos )
{
  // SETUP
  Position p;
  TitanGCM gcmTitan;

  // GIVEN (INPUTS)
  p.height = 1034.5;
  p.longitude = 56.7_deg;  // not used
  p.latitude = 12.3_deg;
  gcmTitan.setPosition(p);
  EphemerisState ephem;
  ephem.solarTime = 23.4;
  ephem.longitudeSun = 34.5_deg;
  gcmTitan.setEphemerisState(ephem);

  // RUN
  gcmTitan.update();
  AtmosphereState s = gcmTitan.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(155.28047704575522, s.temperature);
  EXPECT_DOUBLE_EQ(1.1423293040889387e-06, s.pressure);
  EXPECT_DOUBLE_EQ(69.624167828949268, s.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(2.3817754693987277e-11, s.density);
  EXPECT_DOUBLE_EQ(68.374514342564439, s.densityScaleHeight);
  EXPECT_DOUBLE_EQ(4393141993332.0396, s.argon.numberDensity);
  EXPECT_DOUBLE_EQ(26.919100372619514, s.averageMolecularWeight);
  EXPECT_DOUBLE_EQ(51070122439976.555, s.methane.numberDensity);
  EXPECT_DOUBLE_EQ(477369814556176.63, s.dinitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(532833078989485.25, s.totalNumberDensity);
  EXPECT_DOUBLE_EQ(109.32620929908805, s.ewWind);
  EXPECT_DOUBLE_EQ(0.0000000000000000, s.nsWind);

  // TEAR-DOWN
}

};