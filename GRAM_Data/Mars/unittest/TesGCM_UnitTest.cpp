#include "gtest/gtest.h"
#define GRAM_UNIT_TEST
#include "TesGCM.h"

namespace GRAM {

TEST(TesGCM, updateHeightOffsets)
{
  // SETUP
  TesGCM tesgcm;
  Position position;
  EphemerisState ephem;
  MarsInputParameters params;

  // GIVEN (INPUTS)
  params.offsetModel = MARS_CONSTANT;
  params.constantHeightOffset = 23.4;
  params.mapYear = 1;
  tesgcm.setInputParameters(params);
  position.height = 100.0;
  position.latitude = 10.0_deg;
  position.setLongitude(10.0_deg, WEST_POSITIVE);
  tesgcm.setPosition(position);
  ephem.longitudeSun = 10.0_deg;
  ephem.solarTime = 10.0;
  tesgcm.setEphemerisState(ephem);

  // RUN
  tesgcm.updateHeightOffsets();
  greal heightOffset = tesgcm.heightOffset;
  greal currentOffset = tesgcm.currentOffset;

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(23.4, heightOffset);
  EXPECT_NEAR(-1.3155907623503034, currentOffset, 1.0e-10);

  // GIVEN (INPUTS)
  params.offsetModel = MARS_SEASONAL;
  params.constantHeightOffset = 23.4;
  params.mapYear = 1;
  tesgcm.setInputParameters(params);
  position.height = 123.4;
  position.latitude = 44.4_deg;
  position.setLongitude(177.7_deg, WEST_POSITIVE);
  tesgcm.setPosition(position);
  ephem.longitudeSun = 188.8_deg;
  ephem.solarTime = 15.6;
  tesgcm.setEphemerisState(ephem);

  // RUN
  tesgcm.updateHeightOffsets();
  heightOffset = tesgcm.heightOffset;
  currentOffset = tesgcm.currentOffset;

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(23.762843868900532, heightOffset);
  EXPECT_DOUBLE_EQ(-0.79570503624574329, currentOffset);

  // GIVEN (INPUTS)
  params.offsetModel = MARS_GLOBAL_MEAN;
  params.constantHeightOffset = 23.4;
  params.mapYear = 1;
  tesgcm.setInputParameters(params);
  position.height = 123.4;
  position.latitude = 44.4_deg;
  position.setLongitude(177.7_deg, WEST_POSITIVE);
  tesgcm.setPosition(position);
  ephem.longitudeSun = 188.8_deg;
  ephem.solarTime = 15.6;
  tesgcm.setEphemerisState(ephem);

  // RUN
  tesgcm.updateHeightOffsets();
  heightOffset = tesgcm.heightOffset;
  currentOffset = tesgcm.currentOffset;

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-1.3302666666666665, heightOffset);
  EXPECT_DOUBLE_EQ(0.0, currentOffset);

  // GIVEN (INPUTS)
  params.offsetModel = MARS_DAILY_AVERAGE;
  params.constantHeightOffset = 23.4;
  params.mapYear = 2;
  tesgcm.setInputParameters(params);
  position.height = 123.4;
  position.latitude = 88.8_deg;
  position.setLongitude(177.7_deg, WEST_POSITIVE);
  tesgcm.setPosition(position);
  ephem.longitudeSun = 188.8_deg;
  ephem.solarTime = 15.6;
  tesgcm.setEphemerisState(ephem);

  // RUN
  tesgcm.updateHeightOffsets();
  heightOffset = tesgcm.heightOffset;
  currentOffset = tesgcm.currentOffset;

  // EXPECT (OUTPUTS)
  EXPECT_NEAR(-0.35094956007737005, heightOffset, 1.0e-14);
  EXPECT_NEAR(-0.35094956007737005, currentOffset, 1.0e-14);

  // GIVEN (INPUTS)
  params.offsetModel = MARS_CURRENT;
  params.constantHeightOffset = 23.4;
  params.mapYear = 2;
  tesgcm.setInputParameters(params);
  position.height = 123.4;
  position.latitude = 18.8_deg;
  position.setLongitude(177.7_deg, WEST_POSITIVE);
  tesgcm.setPosition(position);
  ephem.longitudeSun = 188.8_deg;
  ephem.solarTime = 15.6;
  tesgcm.setEphemerisState(ephem);

  // RUN
  tesgcm.updateHeightOffsets();
  heightOffset = tesgcm.heightOffset;
  currentOffset = tesgcm.currentOffset;

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.4457611532608123, heightOffset);
  EXPECT_DOUBLE_EQ(1.4457611532608123, currentOffset);

  // TEAR-DOWN
}

TEST(TesGCM, updateFairingHeights)
{
  // SETUP
  TesGCM tesgcm;
  Position position;

  // GIVEN (INPUTS)
  position.height = 100.0_km;
  position.surfaceHeight = 3.45_km;
  tesgcm.setPosition(position);
  tesgcm.tesUpper.setHeightOffset(2.34_km);

  // RUN
  tesgcm.updateFairingHeights();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(75.0, tesgcm.lowerFairingHeight);
  EXPECT_DOUBLE_EQ(82.34, tesgcm.upperFairingHeight);
  EXPECT_DOUBLE_EQ(3.48, tesgcm.surfaceLowerFairingHeight);
  EXPECT_DOUBLE_EQ(4.0, tesgcm.surfaceUpperFairingHeight);

  // GIVEN (INPUTS)
  position.height = 5.0_km;
  position.surfaceHeight = 7.45_km;
  tesgcm.setPosition(position);
  tesgcm.tesUpper.setHeightOffset(-5.34_km);

  // RUN
  tesgcm.updateFairingHeights();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(75.0, tesgcm.lowerFairingHeight);
  EXPECT_DOUBLE_EQ(79.66, tesgcm.upperFairingHeight);
  EXPECT_DOUBLE_EQ(7.48, tesgcm.surfaceLowerFairingHeight);
  EXPECT_DOUBLE_EQ(8.0, tesgcm.surfaceUpperFairingHeight);

  // TEAR-DOWN
}


}