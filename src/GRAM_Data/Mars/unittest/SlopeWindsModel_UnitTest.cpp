#include "gtest/gtest.h"
#define GRAM_UNIT_TEST
#include "SlopeWindsModel.h"

namespace GRAM {

TEST(SlopeWindsModel, updateSlopes)
{
  // SETUP
  SlopeWindsModel swm;

  // GIVEN (INPUTS)
  Position pos;
  pos.latitude = 0.0_deg;
  pos.setLongitude(0.0_deg, WEST_POSITIVE);
  pos.latitudeRadius = 3395.62774136_km;
  swm.setPosition(pos);

  // RUN
  swm.updateSlopes();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(4.9089980287415728E-003, swm.nsSlope);
  EXPECT_DOUBLE_EQ(1.2526466714391604E-003, swm.ewSlope);

  // GIVEN (INPUTS)
  pos.latitude = 48.9_deg;
  pos.setLongitude(123.4_deg, WEST_POSITIVE);
  pos.latitudeRadius = 3386.0684724960006_km;
  swm.setPosition(pos);

  // RUN
  swm.updateSlopes();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-2.4595873447851684E-003, swm.nsSlope);
  EXPECT_DOUBLE_EQ(6.4711785677948699E-004, swm.ewSlope);

  // TEAR-DOWN
}

TEST(SlopeWindsModel, slopeWinds)
{
  // SETUP
  SlopeWindsModel swm;
  greal u = 23.4;
  greal v = 23.4;
  greal w = 23.4;

  // GIVEN (INPUTS)
  greal time = 13.7;
  Position pos;
  pos.height = -2.0;
  pos.surfaceHeight = -1.0;
  swm.setPosition(pos);

  // RUN
  swm.slopeWinds(time, u, v, w);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, u);
  EXPECT_DOUBLE_EQ(0.0, v);
  EXPECT_DOUBLE_EQ(0.0, w);

  // GIVEN (INPUTS)
  swm.ewWind = 23.4;
  swm.nsWind = 23.4;
  time = 13.7;
  pos.height = 0.2;
  pos.surfaceHeight = -1.3321410400000002;
  pos.latitude = 0.0;
  swm.setPosition(pos);
  swm.ewSlope = 1.2526466714391604E-003;
  swm.nsSlope = 4.9089980287415728E-003;
  //swm.ewSlopeW = 1.2526466714391604E-003;
  //swm.nsSlopeW = 4.9089980287415728E-003;

  // RUN
  swm.slopeWinds(time, u, v, w);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.4537790635342068, u);
  EXPECT_DOUBLE_EQ(3.9251887322216822, v);
  EXPECT_DOUBLE_EQ(0.16527230123808809, w);

  // TEAR-DOWN
}



}