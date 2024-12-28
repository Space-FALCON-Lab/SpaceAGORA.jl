#include "gtest/gtest.h"
#include "gram.h"
#include "JupiterModel.h"
 
using namespace GRAM;

// The HeightModel has its own unit tests.  So this test is basically ensuring that
// we are getting DPT in the right fields.  No need to get fancy here.
TEST(JupiterModel, update )
{
  // SETUP
  JupiterModel jupiter;

  // GIVEN (INPUTS)
  Position pos;
  pos.height = 1002.100;
  jupiter.setPosition(pos);

  // RUN
  jupiter.update();

  // EXPECT (OUTPUTS)
  const AtmosphereState& atmos = jupiter.getAtmosphereState();
  EXPECT_NEAR(3.0380e-11, atmos.density, 1.0e-15);
  EXPECT_NEAR(1.1330e-04, atmos.pressure, 1.0e-7);
  EXPECT_NEAR(898.10, atmos.temperature, 1.0e-3);
  EXPECT_DOUBLE_EQ(0.0, atmos.totalNumberDensity);
  // Winds model is 0
  EXPECT_DOUBLE_EQ(0.0, atmos.nsWind);
  EXPECT_DOUBLE_EQ(0.0, atmos.ewWind);

  // GIVEN (INPUTS)
  pos.height = -123.11;
  jupiter.setPosition(pos);

  // RUN
  jupiter.update();

  // EXPECT (OUTPUTS)
  const AtmosphereState& atmos2 = jupiter.getAtmosphereState();
  EXPECT_NEAR(1.3379, atmos2.density, 1.0e-5);
  EXPECT_NEAR(1955800.0, atmos2.pressure, 0.1);
  EXPECT_NEAR(413.0, atmos2.temperature, 0.1);
  EXPECT_DOUBLE_EQ(0.0, atmos.totalNumberDensity);
  // Winds model is 0
  EXPECT_DOUBLE_EQ(0.0, atmos2.nsWind);
  EXPECT_DOUBLE_EQ(0.0, atmos2.ewWind);


  // TEAR-DOWN
}
