#include "unittest.h"
#include "EarthAtmosphere.h"

namespace GRAM {

TEST(EarthAtmosphere, initializeData)
{
    // SETUP
  EarthAtmosphere earth;

  // GIVEN (INPUTS)
  earth.month = 3;

  // RUN
  earth.initializeData();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.999, earth.plp[0][0]);
  EXPECT_DOUBLE_EQ(0.908, earth.dlp[0][0]);
  EXPECT_DOUBLE_EQ(0.908, earth.tlp[0][0]);
  EXPECT_DOUBLE_EQ(0.733, earth.ulp[0][0]);
  EXPECT_DOUBLE_EQ(0.733, earth.vlp[0][0]);
  EXPECT_DOUBLE_EQ(0.909, earth.plp[24][9]);
  EXPECT_DOUBLE_EQ(0.131, earth.dlp[24][9]);
  EXPECT_DOUBLE_EQ(0.130, earth.tlp[24][9]);
  EXPECT_DOUBLE_EQ(0.457, earth.ulp[24][9]);
  EXPECT_DOUBLE_EQ(0.457, earth.vlp[24][9]);

  EXPECT_DOUBLE_EQ(-0.070, earth.uds[0][0]);
  EXPECT_DOUBLE_EQ(-0.201, earth.vds[0][0]);
  EXPECT_DOUBLE_EQ(0.172, earth.uds[24][9]);
  EXPECT_DOUBLE_EQ(0.366, earth.vds[24][9]);

  EXPECT_DOUBLE_EQ(-0.869, earth.udl[0][0]);
  EXPECT_DOUBLE_EQ(0.313, earth.uvt[0][0]);
  EXPECT_DOUBLE_EQ(-0.132, earth.udl[24][9]);
  EXPECT_DOUBLE_EQ(-0.012, earth.uvt[24][9]);

  EXPECT_DOUBLE_EQ(7.8, earth.xlbar[0]);
  EXPECT_DOUBLE_EQ(0.63, earth.wr[0]);
  EXPECT_DOUBLE_EQ(840.0, earth.xlbar[28]);
  EXPECT_DOUBLE_EQ(12.01, earth.wr[28]);

  // TEAR-DOWN
}

TEST(EarthAtmosphere, updateSurfaceHeight)
{
  // SETUP
  EarthAtmosphere earth;
  Position pos;

  // GIVEN (INPUTS)
  earth.initializeTopographyData();
  pos.latitude = 34.7;
  pos.longitude = 254.8;
  earth.setPosition(pos);

  // RUN
  earth.updateSurfaceHeight();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.7499999999999969, earth.surfaceHeight);

  // GIVEN (INPUTS)
  pos.latitude = -34.7;
  pos.longitude = 144.1;
  earth.setPosition(pos);

  // RUN
  earth.updateSurfaceHeight();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.1, earth.surfaceHeight);

  // GIVEN (INPUTS)
  pos.latitude = -89.9;
  pos.longitude = 359.8;
  earth.setPosition(pos);

  // RUN
  earth.updateSurfaceHeight();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(2.7, earth.surfaceHeight);

  // GIVEN (INPUTS)
  pos.latitude = 51.7;
  pos.longitude = 0.2;
  earth.setPosition(pos);

  // RUN
  earth.updateSurfaceHeight();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.16500000000000054, earth.surfaceHeight);

  // GIVEN (INPUTS)
  pos.latitude = 89.8;
  pos.longitude = 0.2;
  earth.setPosition(pos);

  // RUN
  earth.updateSurfaceHeight();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, earth.surfaceHeight);

  // TEAR-DOWN
}

TEST(EarthAtmosphere, interpolateRSData)
{
  // SETUP
  EarthAtmosphere earth;
  Position pos;

  // GIVEN (INPUTS)
  greal height = 12.3;

  // RUN
  greal xbar = earth.interpolateRSData(earth.xlbar, height);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(80.89, xbar);

  // GIVEN (INPUTS)
  height = 47.8;

  // RUN
  xbar = earth.interpolateRSData(earth.xlbar, height);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(249.396, xbar);

  // TEAR-DOWN
}

TEST(EarthAtmosphere, interpolatePCData)
{
  // SETUP
  EarthAtmosphere earth;

  // GIVEN (INPUTS)
  greal height = 12.3;
  greal latitude = 42.7;

  // RUN
  greal u, v;
  earth.interpolatePCData(earth.udl, earth.uvt, height, latitude, u, v);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-0.4007293, u);
  EXPECT_DOUBLE_EQ(0.137439, v);

  // GIVEN (INPUTS)
  height = 47.8;
  latitude = -85.6;

  // RUN
  earth.interpolatePCData(earth.udl, earth.uvt, height, latitude, u, v);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-0.5801832, u);
  EXPECT_DOUBLE_EQ(-0.0168904, v);

  // GIVEN (INPUTS)
  height = 247.8;
  latitude = 80.0;

  // RUN
  earth.interpolatePCData(earth.udl, earth.uvt, height, latitude, u, v);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-0.1295, u);
  EXPECT_DOUBLE_EQ(-0.017, v);

  // GIVEN (INPUTS)
  height = 27.8;
  latitude = 90.0;

  // RUN
  earth.interpolatePCData(earth.udl, earth.uvt, height, latitude, u, v);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-0.6742, u);
  EXPECT_DOUBLE_EQ(-0.02632, v);

  // TEAR-DOWN
}

TEST(EarthAtmosphere, getUS76StandardAtmosphere)
{
  // SETUP
  EarthAtmosphere earth;
  earth.initializeData();
  double t, p, d;
  Position pos;

  // GIVEN (INPUTS)
  pos.height = -10.0;
  earth.setPosition(pos);

  // RUN
  earth.getUS76StandardAtmosphere(t, p, d);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, t);
  EXPECT_DOUBLE_EQ(0.0, p);
  EXPECT_DOUBLE_EQ(0.0, d);

  // GIVEN (INPUTS)
  pos.height = 1001.0;
  earth.setPosition(pos);

  // RUN
  earth.getUS76StandardAtmosphere(t, p, d);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, t);
  EXPECT_DOUBLE_EQ(0.0, p);
  EXPECT_DOUBLE_EQ(0.0, d);

  // GIVEN (INPUTS)
  pos.height = 0.0;
  earth.setPosition(pos);

  // RUN
  earth.getUS76StandardAtmosphere(t, p, d);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(288.14999999999998, t);
  EXPECT_DOUBLE_EQ(101325.0, p);
  EXPECT_DOUBLE_EQ(1.2249767612159046, d);

  // GIVEN (INPUTS)
  pos.height = 34.5;
  earth.setPosition(pos);

  // RUN
  earth.getUS76StandardAtmosphere(t, p, d);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(235.12836638033789, t);
  EXPECT_DOUBLE_EQ(617.28371675137953, p);
  EXPECT_DOUBLE_EQ(0.0091455463647337198, d);

  // GIVEN (INPUTS)
  pos.height = 97.6;
  earth.setPosition(pos);

  // RUN
  earth.getUS76StandardAtmosphere(t, p, d);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(191.16810482287642, t);
  EXPECT_DOUBLE_EQ(0.048248043640846672, p);
  EXPECT_DOUBLE_EQ(8.6754826661565395e-07, d);

  // TEAR-DOWN
}

}
