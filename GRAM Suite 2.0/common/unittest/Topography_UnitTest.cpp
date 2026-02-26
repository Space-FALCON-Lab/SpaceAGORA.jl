#include "gtest/gtest.h"
#define GRAM_UNIT_TEST
#include "Topography.h"

namespace GRAM {

TEST(Topography, updateIndices)
{
  // SETUP
  Topography topo;
  topo.latStep = 1.0;
  topo.latSize = 182;
  topo.lonStep = 1.0;
  topo.lonSize = 360;

  // GIVEN (INPUTS)
  greal lat = 0.0_deg;
  greal lon = 0.0_deg;

  // RUN
  topo.updateIndices(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(90U, topo.latIndex);
  EXPECT_EQ(0U, topo.lonIndex);
  EXPECT_DOUBLE_EQ(0.5, topo.latDisp);
  EXPECT_DOUBLE_EQ(0.0, topo.lonDisp);

  // GIVEN (INPUTS)
  lat = -89.95_deg;
  lon = 359.9_deg;

  // RUN
  topo.updateIndices(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(0U, topo.latIndex);
  EXPECT_EQ(359U, topo.lonIndex);
  EXPECT_NEAR(0.1, topo.latDisp, 1.0e-4);
  EXPECT_NEAR(0.9, topo.lonDisp, 1.0e-4);

  // GIVEN (INPUTS)
  lat = -89.25_deg;
  lon = 180.25_deg;

  // RUN
  topo.updateIndices(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(1U, topo.latIndex);
  EXPECT_EQ(180U, topo.lonIndex);
  EXPECT_NEAR(0.25, topo.latDisp, 1.0e-4);
  EXPECT_NEAR(0.25, topo.lonDisp, 1.0e-4);

  // GIVEN (INPUTS)
  lat = 32.7_deg;
  lon = -45.6_deg; // 314.4

  // RUN
  topo.updateIndices(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(123U, topo.latIndex);
  EXPECT_EQ(314U, topo.lonIndex);
  EXPECT_NEAR(0.2, topo.latDisp, 1.0e-4);
  EXPECT_NEAR(0.4, topo.lonDisp, 1.0e-4);

  // GIVEN (INPUTS)
  lat = 89.75_deg;
  lon = 725.5_deg;

  // RUN
  topo.updateIndices(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(180U, topo.latIndex);
  EXPECT_EQ(5U, topo.lonIndex);
  EXPECT_DOUBLE_EQ(0.5, topo.latDisp);
  EXPECT_NEAR(0.5, topo.lonDisp, 1.0e-4);

  // GIVEN (INPUTS)
  lat = 90.0_deg;
  lon = 0.0_deg;

  // RUN
  topo.updateIndices(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(180U, topo.latIndex);
  EXPECT_EQ(0U, topo.lonIndex);
  EXPECT_DOUBLE_EQ(1.0, topo.latDisp);
  EXPECT_NEAR(0.0, topo.lonDisp, 1.0e-4);
}

TEST(Topography, updateIndices_oddSizes)
{
  // SETUP
  Topography topo;
  topo.latStep = 0.25;
  topo.latSize = 722;
  topo.lonStep = 0.3;
  topo.lonSize = 1200;

  // GIVEN (INPUTS)
  greal lat = 0.0_deg;
  greal lon = 0.0_deg;

  // RUN
  topo.updateIndices(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(360U, topo.latIndex);
  EXPECT_EQ(0U, topo.lonIndex);
  EXPECT_DOUBLE_EQ(0.5, topo.latDisp);
  EXPECT_DOUBLE_EQ(0.0, topo.lonDisp);

  // GIVEN (INPUTS)
  lat = -89.95_deg;
  lon = 359.9_deg;

  // RUN
  topo.updateIndices(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(0U, topo.latIndex);
  EXPECT_EQ(1199U, topo.lonIndex);
  EXPECT_NEAR(0.4, topo.latDisp, 1.0e-4);
  EXPECT_NEAR(0.666666, topo.lonDisp, 1.0e-4);

  // GIVEN (INPUTS)
  lat = -89.25_deg;
  lon = 180.1_deg;

  // RUN
  topo.updateIndices(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(3U, topo.latIndex);
  EXPECT_EQ(600U, topo.lonIndex);
  EXPECT_NEAR(0.5, topo.latDisp, 1.0e-4);
  EXPECT_NEAR(0.333333, topo.lonDisp, 1.0e-4);

  // GIVEN (INPUTS)
  lat = 32.7_deg;
  lon = -45.6_deg; // 314.4

  // RUN
  topo.updateIndices(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(491U, topo.latIndex);
  EXPECT_EQ(1048U, topo.lonIndex);
  EXPECT_NEAR(0.3, topo.latDisp, 1.0e-4);
  EXPECT_NEAR(0.0, topo.lonDisp, 1.0e-4);

  // GIVEN (INPUTS)
  lat = 89.95_deg;
  lon = 725.5_deg;

  // RUN
  topo.updateIndices(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(720U, topo.latIndex);
  EXPECT_EQ(18U, topo.lonIndex);
  EXPECT_NEAR(0.6, topo.latDisp, 1.0e-4);
  EXPECT_NEAR(0.333333, topo.lonDisp, 1.0e-4);

  // GIVEN (INPUTS)
  lat = 90.0_deg;
  lon = 0.0_deg;

  // RUN
  topo.updateIndices(lat, lon);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(720U, topo.latIndex);
  EXPECT_EQ(0U, topo.lonIndex);
  EXPECT_DOUBLE_EQ(1.0, topo.latDisp);
  EXPECT_NEAR(0.0, topo.lonDisp, 1.0e-4);
}

}