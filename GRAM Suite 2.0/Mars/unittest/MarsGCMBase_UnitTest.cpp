#include "gtest/gtest.h"
#define GRAM_UNIT_TEST
#include "MGCM.h"

namespace GRAM {

TEST(MarsGCMBase, bltp)
{
  // SETUP
  MGCM mgcm;  // Using MGCM to subclass MarsGCMBase

  // GIVEN (INPUTS)
  greal gz = 3.7;
  greal tg = 23.4;
  greal z5 = 23.4;
  greal t5 = 23.4;
  greal u5 = 23.4;
  greal v5 = 23.4;
  greal zeval = 23.4;
  greal factor = 23.4;

  // RUN
  greal temperature = mgcm.bltp(tg, z5, t5, u5, v5, zeval, factor, gz);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(24.905656831304885, temperature);

  // TEAR-DOWN
}

TEST(MarsGCMBase, updateGasConstant)
{
  // SETUP
  MGCM mgcm;  // Using MGCM to subclass MarsGCMBase
  Position position;

  // GIVEN (INPUTS)
  position.height = 2.0_km;
  position.surfaceHeight = 5.0_km;
  mgcm.setPosition(position);

  MarsAtmosphereState& lowerMars = mgcm.lower.getPlanetSpecificMetrics<MarsAtmosphereState>();
  MarsAtmosphereState& upperMars = mgcm.upper.getPlanetSpecificMetrics<MarsAtmosphereState>();

  mgcm.lower.pressure = 300.0;
  mgcm.upper.pressure = 250.0;
  lowerMars.pressureDaily = 300.0;
  upperMars.pressureDaily = 250.0;
  mgcm.lower.density = 0.008;
  mgcm.upper.density = 0.007;
  lowerMars.densityDaily = 0.008;
  upperMars.densityDaily = 0.007;
  mgcm.lower.temperature = 500.0;
  mgcm.upper.temperature = 450.0;
  lowerMars.temperatureDaily = 500.0;
  upperMars.temperatureDaily = 450.0;

  // RUN
  mgcm.updateGasConstant();
  AtmosphereState atmos = mgcm.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(79.3650793650794, atmos.specificGasConstant);
//  EXPECT_DOUBLE_EQ(79.3650793650794, atmos.specificGasConstantDaily);

  // GIVEN (INPUTS)
  position.height = 12.0_km;
  position.surfaceHeight = 5.0_km;
  mgcm.setPosition(position);

  mgcm.lowerHeight = 10.0;
  mgcm.upperHeight = 15.0;

  mgcm.lower.pressure = 300.0;
  mgcm.upper.pressure = 250.0;
  lowerMars.pressureDaily = 300.0;
  upperMars.pressureDaily = 250.0;
  mgcm.lower.density = 0.008;
  mgcm.upper.density = 0.007;
  lowerMars.densityDaily = 0.008;
  upperMars.densityDaily = 0.007;
  mgcm.lower.temperature = 500.0;
  mgcm.upper.temperature = 450.0;
  lowerMars.temperatureDaily = 500.0;
  upperMars.temperatureDaily = 450.0;

  // RUN
  mgcm.updateGasConstant();
  atmos = mgcm.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(76.746031746031747, atmos.specificGasConstant);

  // TEAR-DOWN
}


}