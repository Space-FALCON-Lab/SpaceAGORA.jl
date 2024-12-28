#include "gtest/gtest.h"
#define GRAM_UNIT_TEST
#include "MarsDustModelBase.h"

namespace GRAM {

class MarsDustModelTest : public MarsDustModelBase
{
public:
  MarsDustModelTest() = default;
  virtual ~MarsDustModelTest() = default;
  void updateDustOpticalDepth() {}
#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(MarsDustModelBase, getDustIntensityFactor_zero);
  FRIEND_TEST(MarsDustModelBase, getDustIntensityFactor_nonzero);
  FRIEND_TEST(MarsDustModelBase, getDustOffset);
#endif // GRAM_UNIT_TEST
};

TEST(MarsDustModelBase, getDustIntensityFactor_zero)
{
  // SETUP
  MarsDustModelTest dustModel;
  MarsInputParameters params;

  // GIVEN (INPUTS)
  dustModel.setLongitudeSun(45.0_deg);
  params.stormLongitudeSun = 100.0_deg;
  params.stormDuration = 50.0_deg;
  params.stormIntensity = 5.0;
  dustModel.setInputParameters(params);

  // RUN
  dustModel.updateDustIntensity();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, dustModel.intensity);

  // GIVEN (INPUTS)
  dustModel.setLongitudeSun(152.0_deg);
  params.stormLongitudeSun = 100.0_deg;
  params.stormDuration = 50.0_deg;
  params.stormIntensity = 5.0;
  dustModel.setInputParameters(params);

  // RUN
  dustModel.updateDustIntensity();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, dustModel.intensity);

  // GIVEN (INPUTS)
  dustModel.setLongitudeSun(130.0_deg);
  params.stormLongitudeSun = 100.0_deg;
  params.stormDuration = 50.0_deg;
  params.stormIntensity = 0.0;
  dustModel.setInputParameters(params);

  // RUN
  dustModel.updateDustIntensity();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, dustModel.intensity);

  // GIVEN (INPUTS)
  dustModel.setLongitudeSun(30.0_deg);
  params.stormLongitudeSun = 350.0_deg;
  params.stormDuration = 30.0_deg;
  params.stormIntensity = 0.0;
  dustModel.setInputParameters(params);

  // RUN
  dustModel.updateDustIntensity();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, dustModel.intensity);

  // TEAR-DOWN
}

TEST(MarsDustModelBase, getDustIntensityFactor_nonzero)
{
  // SETUP
  MarsDustModelTest dustModel;
  MarsInputParameters params;
  Position position;

  // GIVEN (INPUTS)
  dustModel.setLongitudeSun(105.0_deg);
  params.stormLongitudeSun = 100.0_deg;
  params.stormDuration = 50.0_deg;
  params.stormIntensity = 5.0;
  dustModel.setInputParameters(params);

  // RUN
  dustModel.updateDustIntensity();

  // EXPECT (OUTPUTS)
  EXPECT_NEAR(4.16666666666667, dustModel.intensity, 1.0e-14);

  // GIVEN (INPUTS)
  dustModel.setLongitudeSun(120.0_deg);
  position.latitude = 20.0_deg;
  position.setLongitude(50.0_deg, WEST_POSITIVE);
  position.latitudeRadius = 5.0_km;
  dustModel.setPosition(position);
  params.isEastLongitudePositiveOnInput = false;
  params.stormLongitudeSun = 100.0_deg;
  params.stormDuration = 50.0_deg;
  params.stormIntensity = 5.0;
  params.stormMaxRadius = 20.0;
  params.stormLatitude = 30.0_deg;
  params.stormLongitude = 40.0_deg;
  dustModel.setInputParameters(params);

  // RUN
  dustModel.updateDustIntensity();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(4.9889510805803026, dustModel.intensity);

  // GIVEN (INPUTS)
  dustModel.setLongitudeSun(140.0_deg);
  params.stormLongitudeSun = 100.0_deg;
  params.stormDuration = 50.0_deg;
  params.stormIntensity = 5.0;
  params.stormMaxRadius = 20.0;
  dustModel.setInputParameters(params);

  // RUN
  dustModel.updateDustIntensity();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.6337149464482652, dustModel.intensity);

  // GIVEN (INPUTS near 360 degrees)
  dustModel.setLongitudeSun(30.0_deg);
  params.stormLongitudeSun = 350.0_deg;
  params.stormDuration = 50.0_deg;
  params.stormIntensity = 5.0;
  params.stormMaxRadius = 20.0;
  dustModel.setInputParameters(params);

  // RUN
  dustModel.updateDustIntensity();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.6337149464482652, dustModel.intensity);
}

TEST(MarsDustModelBase, getDustOffset)
{
  // SETUP
  MarsDustModelTest dustModel;
  MarsInputParameters params;

  // GIVEN (INPUTS)
  params.stormIntensity = 9.0;
  dustModel.setInputParameters(params);
  dustModel.intensity = 0.0;

  // RUN
  greal offset = dustModel.getDustOffset();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, offset);

  // GIVEN (INPUTS)
  params.stormIntensity = 0.0;
  dustModel.setInputParameters(params);
  dustModel.intensity = 9.0;

  // RUN
  offset = dustModel.getDustOffset();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, offset);

  // GIVEN (INPUTS)
  params.stormIntensity = 9.0;
  dustModel.setInputParameters(params);
  dustModel.intensity = 0.97314947886294711;

  // RUN
  offset = dustModel.getDustOffset();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.6219157981049119, offset);
}


}