#include "gtest/gtest.h"
#define GRAM_UNIT_TEST
#include "MGCMDustModel.h"

namespace GRAM {

TEST(MGCMDustModel, update)
{
  // SETUP
  MGCMDustModel dustModel;
  MarsInputParameters params;

  // GIVEN (INPUTS)
  params.mgcmConstantDustLevel = 0.0;
  params.mgcmMaxDustLevel = 100.0;
  params.mgcmMinDustLevel = 40.0;
  dustModel.setInputParameters(params);

  // RUN
  dustModel.setLongitudeSun(0.0_deg);
  dustModel.update();
  greal sod0 = dustModel.getDustOpticalDepth();

  dustModel.setLongitudeSun(90.0_deg);
  dustModel.update();
  greal sod90 = dustModel.getDustOpticalDepth();

  dustModel.setLongitudeSun(180.0_deg);
  dustModel.update();
  greal sod180 = dustModel.getDustOpticalDepth();

  dustModel.setLongitudeSun(270.0_deg);
  dustModel.update();
  greal sod270 = dustModel.getDustOpticalDepth();

  dustModel.setLongitudeSun(30.0_deg);
  dustModel.update();
  greal sod30 = dustModel.getDustOpticalDepth();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(70.0, sod0);
  EXPECT_DOUBLE_EQ(40.0, sod90);
  EXPECT_DOUBLE_EQ(70.0, sod180);
  EXPECT_DOUBLE_EQ(100.0, sod270);
  EXPECT_DOUBLE_EQ(55.0, sod30);

}

}