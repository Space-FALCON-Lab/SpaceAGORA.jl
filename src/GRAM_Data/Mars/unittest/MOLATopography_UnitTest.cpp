#include "gtest/gtest.h"
#define GRAM_UNIT_TEST
#include "MOLATopography.h"

namespace GRAM {

TEST(MOLATopography, updateIndices)
{
  // SETUP
  MOLATopography mola;

  // GIVEN (INPUTS)
  greal lat = 0.0_deg;
  greal lon = 0.0_deg;
  MOLATopography::MOLAParams topoParams;
  topoParams.latSize = MOLATopography::LAT_SIZE;
  topoParams.lonSize = MOLATopography::LON_SIZE;
  topoParams.stepSize = 0.5;

  // RUN
  mola.updateIndices(lat, lon, topoParams);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(180U, topoParams.latIndex);
  EXPECT_EQ(0U, topoParams.lonIndex);
  EXPECT_DOUBLE_EQ(0.5, topoParams.latDisp);
  EXPECT_DOUBLE_EQ(0.976, topoParams.lonDisp);

  // GIVEN (INPUTS)
  lat = -89.95_deg;
  lon = 360.0_deg - 0.010_deg;

  // RUN
  mola.updateIndices(lat, lon, topoParams);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(0U, topoParams.latIndex);
  EXPECT_EQ(0U, topoParams.lonIndex);
  EXPECT_NEAR(0.2, topoParams.latDisp, 1.0e-4);
  EXPECT_NEAR(0.996, topoParams.lonDisp, 1.0e-4);

  // GIVEN (INPUTS)
  lat = -89.75_deg;
  lon = 360.0_deg - 0.012_deg;

  // RUN
  mola.updateIndices(lat, lon, topoParams);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(1U, topoParams.latIndex);
  EXPECT_EQ(1U, topoParams.lonIndex);
  EXPECT_NEAR(0.0, topoParams.latDisp, 1.0e-4);
  EXPECT_NEAR(0.0, topoParams.lonDisp, 1.0e-4);

  // GIVEN (INPUTS)
  lat = -32.7_deg;
  lon = 360.0_deg - 45.6_deg;

  // RUN
  mola.updateIndices(lat, lon, topoParams);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(115U, topoParams.latIndex);
  EXPECT_EQ(92U, topoParams.lonIndex);
  EXPECT_NEAR(0.1, topoParams.latDisp, 1.0e-4);
  EXPECT_NEAR(0.176, topoParams.lonDisp, 1.0e-4);

  // GIVEN (INPUTS)
  lat = 89.75_deg;
  lon = 360.0_deg - 359.761_deg;

  // RUN
  mola.updateIndices(lat, lon, topoParams);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(360U, topoParams.latIndex);
  EXPECT_EQ(720U, topoParams.lonIndex);
  EXPECT_DOUBLE_EQ(0.0, topoParams.latDisp);
  EXPECT_NEAR(0.498, topoParams.lonDisp, 1.0e-4);

  // GIVEN (INPUTS)
  lat = 90.0_deg;
  lon = 0.0_deg;

  // RUN
  mola.updateIndices(lat, lon, topoParams);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(360U, topoParams.latIndex);
  EXPECT_EQ(0U, topoParams.lonIndex);
  EXPECT_DOUBLE_EQ(1.0, topoParams.latDisp);
  EXPECT_NEAR(0.976, topoParams.lonDisp, 1.0e-4);

}

TEST(MOLATopography, updateAlbedoIndices)
{
  // SETUP
  MOLATopography mola;
  MOLATopography::MOLAParams albedoParams;

  // GIVEN (INPUTS)
  greal lat = 0.0_deg;
  greal lon = 0.0_deg;
  albedoParams.latSize = MOLATopography::ALB_LAT_SIZE;
  albedoParams.lonSize = MOLATopography::ALB_LON_SIZE;
  albedoParams.stepSize = 1.0;

  // RUN
  mola.updateIndices(lat, lon, albedoParams);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(90U, albedoParams.latIndex);
  EXPECT_EQ(0U, albedoParams.lonIndex);
  EXPECT_DOUBLE_EQ(0.5, albedoParams.latDisp);
  EXPECT_DOUBLE_EQ(0.738, albedoParams.lonDisp);

  // GIVEN (INPUTS)
  lat = -89.95_deg;
  lon = 360.0_deg - 0.010_deg;

  // RUN
  mola.updateIndices(lat, lon, albedoParams);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(0U, albedoParams.latIndex);
  EXPECT_EQ(0U, albedoParams.lonIndex);
  EXPECT_NEAR(0.1, albedoParams.latDisp, 1.0e-4);
  EXPECT_NEAR(0.748, albedoParams.lonDisp, 1.0e-4);

  // GIVEN (INPUTS)
  lat = -89.5_deg;
  lon = 360.0_deg - 0.262_deg;

  // RUN
  mola.updateIndices(lat, lon, albedoParams);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(1U, albedoParams.latIndex);
  EXPECT_EQ(1U, albedoParams.lonIndex);
  EXPECT_NEAR(0.0, albedoParams.latDisp, 1.0e-4);
  EXPECT_NEAR(0.0, albedoParams.lonDisp, 1.0e-4);

  // GIVEN (INPUTS)
  lat = -32.7_deg;
  lon = 360.0_deg - 45.6_deg;

  // RUN
  mola.updateIndices(lat, lon, albedoParams);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(57U, albedoParams.latIndex);
  EXPECT_EQ(46U, albedoParams.lonIndex);
  EXPECT_NEAR(0.8, albedoParams.latDisp, 1.0e-4);
  EXPECT_NEAR(0.338, albedoParams.lonDisp, 1.0e-4);

  // GIVEN (INPUTS)
  lat = -35.7_deg;
  lon = 360.0_deg - 123.45_deg;

  // RUN
  mola.updateIndices(lat, lon, albedoParams);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(54U, albedoParams.latIndex);
  EXPECT_EQ(124U, albedoParams.lonIndex);
  EXPECT_NEAR(0.8, albedoParams.latDisp, 1.0e-4);
  EXPECT_NEAR(0.188, albedoParams.lonDisp, 1.0e-4);

  // GIVEN (INPUTS)
  lat = 89.5_deg;
  lon = 360.0_deg - 359.761_deg;

  // RUN
  mola.updateIndices(lat, lon, albedoParams);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(180U, albedoParams.latIndex);
  EXPECT_EQ(360U, albedoParams.lonIndex);
  EXPECT_DOUBLE_EQ(0.0, albedoParams.latDisp);
  EXPECT_NEAR(0.499, albedoParams.lonDisp, 1.0e-6);

  // GIVEN (INPUTS)
  lat = 90.0_deg;
  lon = 0.0_deg;

  // RUN
  mola.updateIndices(lat, lon, albedoParams);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(180U, albedoParams.latIndex);
  EXPECT_EQ(0U, albedoParams.lonIndex);
  EXPECT_DOUBLE_EQ(1.0, albedoParams.latDisp);
  EXPECT_NEAR(0.738, albedoParams.lonDisp, 1.0e-4);

}

TEST(MOLATopography, getAreoidRadius)
{
  // SETUP
  MOLATopography mola;

  // GIVEN (INPUTS)
  greal lat = 0.25_deg;
  greal lon = 360.0_deg - (0.25_deg - 0.238_deg);

  // RUN
  greal radius = mola.getAreoidRadius(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(3395.62726, radius);

  // GIVEN (INPUTS)
  lat = 0.25_deg;
  lon = 0.25_deg + 0.238_deg;

  // RUN
  radius = mola.getAreoidRadius(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(3395.63311, radius);

  // GIVEN (INPUTS)
  lat = 0.0_deg;
  lon = 360.0_deg - (-0.238_deg);

  // RUN
  radius = mola.getAreoidRadius(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(3395.630545, radius);

  // GIVEN (INPUTS)
  lat = -35.7_deg;
  lon = 360.0_deg - 123.45_deg;

  // RUN
  radius = mola.getAreoidRadius(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(3390.0488369160003, radius);

  // TEAR-DOWN
}

TEST(MOLATopography, getTopographicHeight)
{
  // SETUP
  MOLATopography mola;

  // GIVEN (INPUTS)
  greal lat = 0.25_deg;
  greal lon = 360.0_deg - (0.25_deg - 0.238_deg);

  // RUN
  greal radius = mola.getTopographicHeight(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-1.2621500000000001, radius);

  // GIVEN (INPUTS)
  lat = 0.25_deg;
  lon = 0.25_deg + 0.238_deg;

  // RUN
  radius = mola.getTopographicHeight(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-1.1479100000000000, radius);

  // GIVEN (INPUTS)
  lat = 0.0_deg;
  lon = 360.0_deg - (-0.238_deg);

  // RUN
  radius = mola.getTopographicHeight(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-1.2445974999999998, radius);

  // GIVEN (INPUTS)
  lat = -35.7_deg;
  lon = 360.0_deg - 123.45_deg;

  // RUN
  radius = mola.getTopographicHeight(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(3.1261365279999991, radius);

  // TEAR-DOWN
}

TEST(MOLATopography, getAlbedo)
{
  // SETUP
  MOLATopography mola;

  // GIVEN (INPUTS)
  greal lat = 0.25_deg;
  greal lon = 360.0_deg - (0.25_deg - 0.238_deg);

  // RUN
  greal albedo = mola.getAlbedo(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.19, albedo);

  // GIVEN (INPUTS)
  lat = 0.25_deg;
  lon = 0.25_deg + 0.238_deg;

  // RUN
  albedo = mola.getAlbedo(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.19, albedo);

  // GIVEN (INPUTS)
  lat = 0.0_deg;
  lon = 360.0_deg - (-0.238_deg);

  // RUN
  albedo = mola.getAlbedo(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.19, albedo);

  // GIVEN (INPUTS)
  lat = -35.7_deg;
  lon = 360.0_deg - 123.45_deg;

  // RUN
  albedo = mola.getAlbedo(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.16049600000000003, albedo);

  // TEAR-DOWN
}

}