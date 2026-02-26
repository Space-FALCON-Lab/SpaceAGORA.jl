#include "gtest/gtest.h"
#include "gram.h"
#include "UranusModel.h"
 
using namespace GRAM;

// The HeightModel has its own unit tests.  So this test is basically ensuring that
// we are getting DPT in the right fields.  No need to get fancy here.
TEST(UranusModel, update )
{
  // SETUP
  UranusModel uranus;

  // GIVEN (INPUTS)
  Position pos;
  pos.height = 1000.000;
  uranus.setPosition(pos);

  // RUN
  uranus.update();

  // EXPECT (OUTPUTS)
  const AtmosphereState& atmos = uranus.getAtmosphereState();
  EXPECT_NEAR(8.8644e-009, atmos.density, 1.0e-15);
  EXPECT_NEAR(1.8871e-002, atmos.pressure, 1.0e-8);
  EXPECT_NEAR(521.40, atmos.temperature, 1.0e-3);
  EXPECT_NEAR(2.5909e+018, atmos.dihydrogen.numberDensity, 1.0e12);
  EXPECT_NEAR(2.5116e+016, atmos.helium.numberDensity, 1.0e10);
  EXPECT_NEAR(2.6163e+014, atmos.methane.numberDensity, 1.0e8);
  EXPECT_NEAR(2.6163e+018, atmos.totalNumberDensity, 1.0e14);
  // Winds model is 0
  EXPECT_DOUBLE_EQ(0.0, atmos.nsWind);
  EXPECT_DOUBLE_EQ(0.0, atmos.ewWind);

  // GIVEN (INPUTS)
  pos.height = -123.185;
  uranus.setPosition(pos);

  // RUN
  uranus.update();

  // EXPECT (OUTPUTS)
  const AtmosphereState& atmos2 = uranus.getAtmosphereState();
  EXPECT_NEAR(2.4367e+000, atmos2.density, 1.0e-6);
  EXPECT_NEAR(1.5078e+006, atmos2.pressure, 1.0e0);
  EXPECT_NEAR(176.06, atmos2.temperature, 1.0e-4);
  EXPECT_NEAR(4.6434e+026, atmos2.dihydrogen.numberDensity, 1.0e20);
  EXPECT_NEAR(8.6158e+025, atmos2.helium.numberDensity, 1.0e19);
  EXPECT_NEAR(1.1521e+025, atmos2.methane.numberDensity, 1.0e19);
  EXPECT_NEAR(5.6202e+026, atmos2.totalNumberDensity, 1.0e22);
  // Winds model is 0
  EXPECT_DOUBLE_EQ(0.0, atmos2.nsWind);
  EXPECT_DOUBLE_EQ(0.0, atmos2.ewWind);


  // TEAR-DOWN
}
