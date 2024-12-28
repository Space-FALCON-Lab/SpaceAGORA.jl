#include "unittest.h"
#include "NCEP.h"
 
namespace GRAM {

TEST(NCEP, getIndex)
{
  // SETUP
  NCEP ncep;

  // GIVEN (INPUTS)
  int i = 1;
  int j = 1;
  int k = 1;

  // RUN
  size_t index = ncep.getIndex(i, j, k);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(size_t(10731), index);

  // GIVEN (INPUTS)
  size_t ii = 18;
  size_t jj = 72;
  size_t kk = 144;

  // RUN
  index = ncep.getIndex(ii, jj, kk);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(size_t(201114), index);

  // TEAR-DOWN
}

TEST(NCEP, readNCEPFile)
{
  // SETUP
  NCEP ncep;

  // GIVEN (INPUTS)
  ncep.NCEPYear = 9715;
  ncep.NCEPHour = 1;
  ncep.month = 6;
  ncep.initialized = false;

  // RUN
  ncep.readNCEPFile();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(237.80416743499907, ncep.temp[0]);
  EXPECT_DOUBLE_EQ(2.7405161523045822, ncep.geop[10731]);
  EXPECT_DOUBLE_EQ(32.069499999999998, ncep.geop[201114]);

  // GIVEN (INPUTS)
  ncep.NCEPYear = 9715;
  ncep.NCEPHour = 4;
  ncep.month = 12;
  ncep.initialized = false;

  // RUN
  ncep.readNCEPFile();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(270.0787234042553, ncep.temp[0]);
  EXPECT_DOUBLE_EQ(2.738226996501095, ncep.geop[10731]);
  EXPECT_DOUBLE_EQ(28.7806, ncep.geop[201114]);

  // TEAR-DOWN
}

TEST(NCEP, wexler)
{
  // SETUP
  NCEP ncep;

  // GIVEN (INPUTS)

  // RUN
  greal wex = ncep.wexler(0.0, 1.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.0e-23, wex);

  // RUN
  wex = ncep.wexler(75.0, 1.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.0e-23, wex);

  // RUN
  wex = ncep.wexler(300.0, 101325.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(3550.2933550897073, wex);

  // RUN
  wex = ncep.wexler(200.0, 10000.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.33266340110751758, wex);

  // TEAR-DOWN
}

TEST(NCEP, getGasConstant)
{
  // SETUP
  NCEP ncep;

  // GIVEN (INPUTS)

  // RUN
  greal con = ncep.getGasConstant(0.0, 1.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(287.05500000000001, con);

  // RUN
  con = ncep.getGasConstant(75.0, 1.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(287.05500000000001, con);

  // RUN
  con = ncep.getGasConstant(300.0, 101325.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(287.09302739172239, con);

  // RUN
  con = ncep.getGasConstant(200.0, 10000.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(287.05503609910716, con);

  // TEAR-DOWN
}

TEST(NCEP, getGeopotentialHeight)
{
  // SETUP
  NCEP ncep;

  // GIVEN (INPUTS)

  // RUN
  greal h = ncep.getGeopotentialHeight(9.8, 6378.137, 0.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, h);

  // RUN
  h = ncep.getGeopotentialHeight(9.8, 6378.137, 10.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(9.9775754861621113, h);

  // RUN
  h = ncep.getGeopotentialHeight(9.9, 6377.0, 100.0);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(99.393283731027424, h);

  // TEAR-DOWN
}

TEST(NCEP, getHeights)
{
  // SETUP
  NCEP ncep;
  ncep.NCEPYear = 9715;
  ncep.NCEPHour = 1;
  ncep.month = 6;
  ncep.readNCEPFile();

  // GIVEN (INPUTS)

  // RUN
  greal high20mb, low10mb;
  ncep.getHeights(0.0, 0.0, high20mb, low10mb);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(26.659301555381152, high20mb);
  EXPECT_DOUBLE_EQ(31.328489672891852, low10mb);
  //EXPECT_DOUBLE_EQ(26.659097395099391, high20mb);
  //EXPECT_DOUBLE_EQ(31.328249579292237, low10mb);

  // RUN
  ncep.getHeights(45.0, 123.0, high20mb, low10mb);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(26.931608483457559, high20mb);
  EXPECT_DOUBLE_EQ(31.625084237678823, low10mb);
  //EXPECT_DOUBLE_EQ(26.932171327557953, high20mb);
  //EXPECT_DOUBLE_EQ(31.625745656779149, low10mb);

  // RUN
  ncep.getHeights(90.0, 213.0, high20mb, low10mb);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(27.273949137955285, high20mb);
  EXPECT_DOUBLE_EQ(32.142689062777769, low10mb);
  //EXPECT_DOUBLE_EQ(27.273913041495113, high20mb);
  //EXPECT_DOUBLE_EQ(32.14264649031346, low10mb);

  // RUN
  ncep.getHeights(-90.0, 321.0, high20mb, low10mb);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(23.707681883903589, high20mb);
  EXPECT_DOUBLE_EQ(27.537723693802935, low10mb);
  //EXPECT_DOUBLE_EQ(23.707650524788079, high20mb);
  //EXPECT_DOUBLE_EQ(27.537687246742163, low10mb);

  // TEAR-DOWN
}

TEST(NCEP, sortLevel)
{
  // SETUP
  NCEP ncep;

  // GIVEN (INPUTS)
  greal zsort[2][2][19] = {
    {{ 0.0, -2.0, -1.0, 1.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0},
     { 0.0, -5.5, -6.0, -3.0, 2.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0}},
    {{ 0.0, 15.5, -3.0, -2.0, -1.0, 1.0, 4.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0},
     { 0.0, 18.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0}}
  };

  size_t ans[2][2][19] = {
    {{ 1, 2, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18},
     { 2, 1, 3, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}},
    {{ 2, 3, 4, 0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 1, 15, 16, 17, 18},
     { 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 1, 18}}
  };

  // RUN
  size_t lsort[2][2][19] = { { { 0 } } };
  ncep.sortLevel(lsort, zsort);

  // EXPECT (OUTPUTS)
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 19; ++k) {
        EXPECT_EQ(ans[i][j][k], lsort[i][j][k]);
      }
    }
  }

  // TEAR-DOWN
}

TEST(NCEP, ncepmd)
{
  // SETUP
  NCEP ncep;
  ncep.NCEPYear = 9715;
  ncep.NCEPHour = 1;
  ncep.month = 6;
  ncep.readNCEPFile();

  // GIVEN (INPUTS)
  Position pos;
  pos.height = 0.0;
  pos.latitude = 0.0;
  pos.longitude = 0.0;
  ncep.setPosition(pos);
  ncep.updatePosition();   // for gravity and radii

  // RUN
  ncep.update();
  const AtmosphereState& atmos = ncep.getAtmosphereState();
  const EarthAtmosphereState& earth = atmos.getPlanetSpecificMetrics<EarthAtmosphereState>();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(101406.9162739202, atmos.pressure);
  EXPECT_DOUBLE_EQ(1.1686509658592317, atmos.density);
  EXPECT_DOUBLE_EQ(298.79890013198417, atmos.temperature);
  EXPECT_DOUBLE_EQ(-0.49535122451972435, atmos.ewWind);
  EXPECT_DOUBLE_EQ(5.5792051620472209, atmos.nsWind);
  EXPECT_DOUBLE_EQ(0.0, atmos.verticalWind);
  EXPECT_DOUBLE_EQ(296.96195923155892, earth.dewPoint);
  EXPECT_DOUBLE_EQ(1.5268979612024316e-06, atmos.pressureStandardDeviation);
  EXPECT_DOUBLE_EQ(1.9640682952710046e-05, atmos.densityStandardDeviation);
  EXPECT_DOUBLE_EQ(1.3552748566213166e-05, atmos.temperatureStandardDeviation);
  EXPECT_DOUBLE_EQ(2.0745563865669454, atmos.ewStandardDeviation);
  EXPECT_DOUBLE_EQ(1.7745563865669451, atmos.nsStandardDeviation);
  EXPECT_DOUBLE_EQ(0.91039741897638948, earth.dewPointSD);
  EXPECT_DOUBLE_EQ(89.261533949259416, earth.relativeHumidity);
  EXPECT_DOUBLE_EQ(3.7907024490394492, earth.relativeHumiditySD);
  EXPECT_DOUBLE_EQ(2979.2172103906028, earth.vaporPressure);
  EXPECT_DOUBLE_EQ(174.28693004113217, earth.vaporPressureSD);
  EXPECT_DOUBLE_EQ(5.9000000000000004, earth.windSpeed);
  EXPECT_DOUBLE_EQ(1.6924769027716675, earth.windSpeedStandardDeviation);
  EXPECT_DOUBLE_EQ(-0.013900278633230678, earth.windCorrelation);

  // GIVEN (INPUTS)
  pos.height = 45.6;
  pos.latitude = -23.4;
  pos.longitude = 167.8;
  ncep.setPosition(pos);
  ncep.updatePosition();   // for gravity and radii

  // RUN
  ncep.update();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(117.38322959127825, atmos.pressure);
  EXPECT_DOUBLE_EQ(0.0016544124916079343, atmos.density);
  EXPECT_DOUBLE_EQ(247.11722586542538, atmos.temperature);
  EXPECT_DOUBLE_EQ(37.810671732097639, atmos.ewWind);
  EXPECT_DOUBLE_EQ(-1.8812018640092028, atmos.nsWind);
  EXPECT_DOUBLE_EQ(-0.0006158309910499288, atmos.verticalWind);
  EXPECT_DOUBLE_EQ(208.62956638788665, earth.dewPoint);
  EXPECT_DOUBLE_EQ(0.0073664009703712126, atmos.pressureStandardDeviation);
  EXPECT_DOUBLE_EQ(0.0083033647040430605, atmos.densityStandardDeviation);
  EXPECT_DOUBLE_EQ(0.00036787002369342273, atmos.temperatureStandardDeviation);
  EXPECT_DOUBLE_EQ(22.542816776129218, atmos.ewStandardDeviation);
  EXPECT_DOUBLE_EQ(10.428816905472985, atmos.nsStandardDeviation);
  EXPECT_DOUBLE_EQ(-1.5781262129570379, earth.dewPointSD);
  EXPECT_DOUBLE_EQ(-4.9268612195520616, earth.relativeHumidity);
  EXPECT_DOUBLE_EQ(-6.7175990654221902, earth.relativeHumiditySD);
  EXPECT_DOUBLE_EQ(0.86591337958627679, earth.vaporPressure);
  EXPECT_DOUBLE_EQ(0.12211755641052817, earth.vaporPressureSD);
  EXPECT_DOUBLE_EQ(39.348471379423799, earth.windSpeed);
  EXPECT_DOUBLE_EQ(18.407519454288131, earth.windSpeedStandardDeviation);
  EXPECT_DOUBLE_EQ(-0.93003231846752055, earth.windCorrelation);

  // TEAR-DOWN
}

} // namespace