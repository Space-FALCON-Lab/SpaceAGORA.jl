#include "gtest/gtest.h"
#define GRAM_UNIT_TEST
#include "TesDustModel.h"

namespace GRAM {

TEST(TesDustModel, update)
{
  // SETUP
  TesDustModel dustModel;
  MarsInputParameters params;
  Position position;

  // GIVEN (INPUTS)
  params.mapYear = 1;
  dustModel.setInputParameters(params);
  position.latitude = 0.0_deg;
  position.setLongitude(0.0_deg, WEST_POSITIVE);
  dustModel.setPosition(position);
  dustModel.setLongitudeSun(0.0_deg);

  // RUN
  dustModel.update();
  greal tesod = dustModel.getDustOpticalDepth();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.3715, tesod);

  // GIVEN (INPUTS)
  params.mapYear = 1;
  dustModel.setInputParameters(params);
  position.latitude = -33.3_deg;
  position.setLongitude(123.4_deg, WEST_POSITIVE);
  dustModel.setPosition(position);
  dustModel.setLongitudeSun(234.5_deg);

  // RUN
  dustModel.update();
  tesod = dustModel.getDustOpticalDepth();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.52843324444444451, tesod);

  // GIVEN (INPUTS)
  params.mapYear = 2;
  dustModel.setInputParameters(params);
  position.latitude = 66.6_deg;
  position.setLongitude(195.7_deg, WEST_POSITIVE);
  dustModel.setPosition(position);
  dustModel.setLongitudeSun(230.4_deg);

  // RUN
  dustModel.update();
  tesod = dustModel.getDustOpticalDepth();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.0942056000000000, tesod);
}

}