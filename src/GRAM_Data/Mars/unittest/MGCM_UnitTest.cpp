#include "gtest/gtest.h"
#define GRAM_UNIT_TEST
#include "MGCM.h"

namespace GRAM {

TEST(MGCM, updateHeightOffsets)
{
  // SETUP
  MGCM mgcm;
  Position position;
  EphemerisState ephem;
  MarsInputParameters params;

  // GIVEN (INPUTS)
  params.offsetModel = MARS_CONSTANT;
  params.constantHeightOffset = 23.4;
  mgcm.setInputParameters(params);
  mgcm.setDustParameters(0.0, 0.1);
  position.height = 100.0;
  position.latitude = 10.0_deg;
  position.setLongitude(10.0_deg, WEST_POSITIVE);
  mgcm.setPosition(position);
  ephem.longitudeSun = 10.0_deg;
  ephem.solarTime = 10.0;
  mgcm.setEphemerisState(ephem);

  // RUN
  mgcm.updateHeightOffsets();
  greal heightOffset = mgcm.heightOffset;
  greal currentOffset = mgcm.currentOffset;

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(23.4, heightOffset);
  EXPECT_NEAR(0.21490659239671187, currentOffset, 1.0e-7);

  // GIVEN (INPUTS)
  params.offsetModel = MARS_SEASONAL;
  params.constantHeightOffset = 23.4;
  mgcm.setInputParameters(params);
  mgcm.setDustParameters(0.0, 0.1);
  position.height = 123.4;
  position.latitude = 44.4_deg;
  position.setLongitude(177.7_deg, WEST_POSITIVE);
  mgcm.setPosition(position);
  ephem.longitudeSun = 188.8_deg;
  ephem.solarTime = 15.6;
  mgcm.setEphemerisState(ephem);

  // RUN
  mgcm.updateHeightOffsets();
  heightOffset = mgcm.heightOffset;
  currentOffset = mgcm.currentOffset;

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(24.353190667792944, heightOffset);
  EXPECT_NEAR(0.12230031308462858, currentOffset, 1.0e-7);

  // GIVEN (INPUTS)
  params.offsetModel = MARS_GLOBAL_MEAN;
  params.constantHeightOffset = 23.4;
  mgcm.setInputParameters(params);
  mgcm.setDustParameters(0.0, 0.3);
  position.height = 123.4;
  position.latitude = 44.4_deg;
  position.setLongitude(177.7_deg, WEST_POSITIVE);
  mgcm.setPosition(position);
  ephem.longitudeSun = 188.8_deg;
  ephem.solarTime = 15.6;
  mgcm.setEphemerisState(ephem);

  // RUN
  mgcm.updateHeightOffsets();
  heightOffset = mgcm.heightOffset;
  currentOffset = mgcm.currentOffset;

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-0.0027066666666666662, heightOffset);
  EXPECT_DOUBLE_EQ(0.0, currentOffset);

  // GIVEN (INPUTS)
  params.offsetModel = MARS_DAILY_AVERAGE;
  params.constantHeightOffset = 23.4;
  mgcm.setInputParameters(params);
  mgcm.setDustParameters(0.0, 1.3);
  position.height = 123.4;
  position.latitude = 88.8_deg;
  position.setLongitude(177.7_deg, WEST_POSITIVE);
  mgcm.setPosition(position);
  ephem.longitudeSun = 188.8_deg;
  ephem.solarTime = 15.6;
  mgcm.setEphemerisState(ephem);

  // RUN
  mgcm.updateHeightOffsets();
  heightOffset = mgcm.heightOffset;
  currentOffset = mgcm.currentOffset;

  // EXPECT (OUTPUTS)
  EXPECT_NEAR(0.89158479960381976, heightOffset, 1.0e-6);
  EXPECT_NEAR(0.89158479960381976, currentOffset, 1.0e-6);

  // GIVEN (INPUTS)
  params.offsetModel = MARS_CURRENT;
  params.constantHeightOffset = 23.4;
  mgcm.setInputParameters(params);
  mgcm.setDustParameters(0.0, 1.3);
  position.height = 123.4;
  position.latitude = 18.8_deg;
  position.setLongitude(177.7_deg, WEST_POSITIVE);
  mgcm.setPosition(position);
  ephem.longitudeSun = 188.8_deg;
  ephem.solarTime = 15.6;
  mgcm.setEphemerisState(ephem);

  // RUN
  mgcm.updateHeightOffsets();
  heightOffset = mgcm.heightOffset;
  currentOffset = mgcm.currentOffset;

  // EXPECT (OUTPUTS)
  EXPECT_NEAR(-0.36937734826354168, heightOffset, 1.0e-6);
  EXPECT_NEAR(-0.36937734826354168, currentOffset, 1.0e-6);

  // TEAR-DOWN
}

TEST(MGCM, updateFairingHeights)
{
  // SETUP
  MGCM mgcm;
  Position position;

  // GIVEN (INPUTS)
  position.height = 100.0_km;
  position.surfaceHeight = 3.45_km;
  mgcm.setPosition(position);
  mgcm.mgcmUpper.setHeightOffset(2.34_km);

  // RUN
  mgcm.updateFairingHeights();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(75.0, mgcm.lowerFairingHeight);
  EXPECT_DOUBLE_EQ(82.34, mgcm.upperFairingHeight);
  EXPECT_DOUBLE_EQ(3.48, mgcm.surfaceLowerFairingHeight);
  EXPECT_DOUBLE_EQ(5.0, mgcm.surfaceUpperFairingHeight);

  // GIVEN (INPUTS)
  position.height = 5.0_km;
  position.surfaceHeight = 7.45_km;
  mgcm.setPosition(position);
  mgcm.mgcmUpper.setHeightOffset(-5.34_km);

  // RUN
  mgcm.updateFairingHeights();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(75.0, mgcm.lowerFairingHeight);
  EXPECT_DOUBLE_EQ(79.66, mgcm.upperFairingHeight);
  EXPECT_DOUBLE_EQ(7.48, mgcm.surfaceLowerFairingHeight);
  EXPECT_DOUBLE_EQ(10.0, mgcm.surfaceUpperFairingHeight);

  // TEAR-DOWN
}

TEST(MGCM, updateUpperLower_1)
{
  // SETUP
  MGCM mgcm;
  Position position;
  MarsInputParameters params;
  mgcm.setInputParameters(params);  // set defaults

  // GIVEN (INPUTS)
  mgcm.upperFairingHeight = 80.0_km;
  mgcm.lowerFairingHeight = 75.0_km;
  mgcm.surfaceUpperFairingHeight = 10.0_km;
  mgcm.surfaceLowerFairingHeight = 5.0_km;
  position.height = 100.0_km;
  mgcm.setPosition(position);
  mgcm.heightOffset = 2.34_km;
  mgcm.mgcmUpper.setHeightOffset(mgcm.heightOffset);
  mgcm.setDustParameters(0.0, 0.1);
  
  // RUN
  mgcm.updateUpperLower();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(124.16673618398966, mgcm.thermosphereBaseHeight);
  EXPECT_DOUBLE_EQ(2.34, mgcm.localHeightOffset);
  EXPECT_DOUBLE_EQ(97.34, mgcm.lowerHeight);
  EXPECT_DOUBLE_EQ(102.34, mgcm.upperHeight);
  EXPECT_DOUBLE_EQ(120.83375842783347, mgcm.lower.temperature);
  EXPECT_DOUBLE_EQ(121.28080147945178, mgcm.upper.temperature);

  // TEAR-DOWN
}

TEST(MGCM, updateUpperLower_2)
{
  // SETUP
  MGCM mgcm;
  Position position;
  MarsInputParameters params;
  mgcm.setInputParameters(params);  // set defaults
  greal defaultThermosphereBase = mgcm.thermosphereBaseHeight;
  greal defaultLocal = mgcm.localHeightOffset;

  // GIVEN (between upperFairingHeight and lowerFairingHeight)
  mgcm.upperFairingHeight = 80.0_km;
  mgcm.lowerFairingHeight = 75.0_km;
  mgcm.surfaceUpperFairingHeight = 10.0_km;
  mgcm.surfaceLowerFairingHeight = 5.0_km;
  position.height = 76.0_km;
  mgcm.setPosition(position);
  mgcm.heightOffset = 2.34_km;
  mgcm.mgcmUpper.setHeightOffset(mgcm.heightOffset);
  mgcm.setDustParameters(0.0, 0.1);

  // RUN
  mgcm.updateUpperLower();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(defaultThermosphereBase, mgcm.thermosphereBaseHeight);
  EXPECT_DOUBLE_EQ(defaultLocal, mgcm.localHeightOffset);
  EXPECT_DOUBLE_EQ(75.0, mgcm.lowerHeight);
  EXPECT_DOUBLE_EQ(82.34, mgcm.upperHeight);
  EXPECT_DOUBLE_EQ(159.90222755257525, mgcm.lower.temperature);
  EXPECT_DOUBLE_EQ(128.82374290274095, mgcm.upper.temperature);

  // GIVEN (between lowerFairingHeight and surfaceUpperFairingHeight)
  position.height = 50.0_km;
  mgcm.setPosition(position);
  mgcm.heightOffset = 2.34_km;
  mgcm.mgcmUpper.setHeightOffset(mgcm.heightOffset);
  mgcm.setDustParameters(0.0, 0.1);

  // RUN
  mgcm.updateUpperLower();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(defaultThermosphereBase, mgcm.thermosphereBaseHeight);
  EXPECT_DOUBLE_EQ(defaultLocal, mgcm.localHeightOffset);
  EXPECT_DOUBLE_EQ(50.0, mgcm.lowerHeight);
  EXPECT_DOUBLE_EQ(55.0, mgcm.upperHeight);
  EXPECT_DOUBLE_EQ(137.66330833846732, mgcm.lower.temperature);
  EXPECT_DOUBLE_EQ(132.47998427901391, mgcm.upper.temperature);

  // GIVEN (between surfaceUpperFairingHeight and surfaceLowerFairingHeight)
  position.height = 7.0_km;
  position.setLongitude(100.0_deg, WEST_POSITIVE);
  position.surfaceHeight = 4.711244_km;
  mgcm.setPosition(position);
  mgcm.heightOffset = 2.34_km;
  mgcm.setDustParameters(0.0, 0.1);

  // RUN
  mgcm.updateUpperLower();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, mgcm.thermosphereBaseHeight);
  EXPECT_DOUBLE_EQ(0.0, mgcm.localHeightOffset);
  EXPECT_DOUBLE_EQ(4.741244, mgcm.lowerHeight);
  EXPECT_DOUBLE_EQ(10.0, mgcm.upperHeight);
  EXPECT_DOUBLE_EQ(205.42765536164569, mgcm.lower.temperature);
  EXPECT_DOUBLE_EQ(200.55036302835120, mgcm.upper.temperature);

  // GIVEN (below surfaceLowerFairingHeight)
  position.height = 1.0_km;
  position.setLongitude(100.0_deg, WEST_POSITIVE);
  position.surfaceHeight = 4.711244_km;
  mgcm.setPosition(position);
  mgcm.heightOffset = 2.34_km;
  mgcm.setDustParameters(0.0, 0.1);

  // RUN
  mgcm.updateUpperLower();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, mgcm.thermosphereBaseHeight);
  EXPECT_DOUBLE_EQ(0.0, mgcm.localHeightOffset);
  EXPECT_DOUBLE_EQ(4.711244, mgcm.lowerHeight);
  EXPECT_DOUBLE_EQ(4.716244, mgcm.upperHeight);
  EXPECT_DOUBLE_EQ(197.87714759008048, mgcm.lower.temperature);
  EXPECT_DOUBLE_EQ(201.95224349013708, mgcm.upper.temperature);

  // TEAR-DOWN
}

TEST(MGCM, updateScaleHeights)
{
  // SETUP
  MGCM mgcm;
  Position position;
  MarsAtmosphereState& lowerMars = mgcm.lower.getPlanetSpecificMetrics<MarsAtmosphereState>();
  MarsAtmosphereState& upperMars = mgcm.upper.getPlanetSpecificMetrics<MarsAtmosphereState>();

  // GIVEN (INPUTS)
  mgcm.upperFairingHeight = 80.0_km;
  mgcm.lowerFairingHeight = 75.0_km;
  mgcm.surfaceUpperFairingHeight = 10.0_km;
  mgcm.surfaceLowerFairingHeight = 5.0_km;

  position.height = 52.0_km;
  position.latitude = 0.0_deg;
  position.setLongitude(100.0_deg, WEST_POSITIVE);
  mgcm.setPosition(position);

  mgcm.lowerHeight = 50.0_km;
  mgcm.upperHeight = 55.0_km;

  mgcm.lower.pressure = 300.0;
  mgcm.upper.pressure = 250.0;
  lowerMars.pressureDaily = 300.0;
  upperMars.pressureDaily = 250.0;
  mgcm.lower.density = 0.008;
  mgcm.upper.density = 0.007;

  // RUN
  mgcm.updateScaleHeights();
  AtmosphereState atmos = mgcm.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(27.4240747387354, atmos.pressureScaleHeight);
//  EXPECT_DOUBLE_EQ(27.4240747387354, atmos.pressureScaleHeightDaily);
  EXPECT_DOUBLE_EQ(37.4443784470931, atmos.densityScaleHeight);

  // GIVEN (INPUTS)
  mgcm.upperFairingHeight = 80.0_km;
  mgcm.lowerFairingHeight = 75.0_km;
  mgcm.surfaceUpperFairingHeight = 10.0_km;
  mgcm.surfaceLowerFairingHeight = 5.0_km;

  position.height = 8.0_km;
  position.latitude = 0.0_deg;
  position.setLongitude(100.0_deg, WEST_POSITIVE);
  mgcm.setPosition(position);

  mgcm.lowerHeight = 5.0_km;
  mgcm.upperHeight = 10.0_km;

  mgcm.lower.pressureScaleHeight = 27.0;
//  mgcm.lower.pressureScaleHeightDaily = 28.0;
  mgcm.lower.temperature = 200.0;
  mgcm.upper.temperature = 211.0;

  // RUN
  mgcm.updateScaleHeights();
  atmos = mgcm.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(27.0, atmos.pressureScaleHeight);
//  EXPECT_DOUBLE_EQ(28.0, atmos.pressureScaleHeightDaily);
  EXPECT_DOUBLE_EQ(20.94563986409966, atmos.densityScaleHeight);

  // TEAR-DOWN
}

}