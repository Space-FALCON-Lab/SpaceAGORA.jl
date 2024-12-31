#include "unittest.h"
#include "TesLowerInterpolator.h"

namespace GRAM {

TEST(TesLowerInterpolator, initializeLowerData_getTideParameters)
{
  // SETUP
  TesLowerInterpolator tesInterp;

  // GIVEN (INPUTS)
  greal poleFactor = 1.0;
  tesInterp.index.hgt = 0;
  tesInterp.index.lat = 0;
  tesInterp.index.ls = 0;
  tesInterp.index.tyr = 0;

  // RUN
  MarsTideParameters t = tesInterp.getTideParameters(tesInterp.mgcmTemperature, tesInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(130.0, t.diurnalMean);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude1D);
  EXPECT_DOUBLE_EQ(12.6, t.phase1D);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude2D);
  EXPECT_DOUBLE_EQ(8.4, t.phase2D);

  // RUN
  t = tesInterp.getTideParameters(tesInterp.mgcmDensity, tesInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(4.139E-02, t.diurnalMean);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude1D);
  EXPECT_DOUBLE_EQ(11.3, t.phase1D);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude2D);
  EXPECT_DOUBLE_EQ(3.7, t.phase2D);

  // RUN
  t = tesInterp.getTideParameters(tesInterp.mgcmEWWind, tesInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.6, t.diurnalMean);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude1D);
  EXPECT_DOUBLE_EQ(5.9, t.phase1D);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude2D);
  EXPECT_DOUBLE_EQ(7.4, t.phase2D);

  // RUN
  t = tesInterp.getTideParameters(tesInterp.mgcmNSWind, tesInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, t.diurnalMean);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude1D);
  EXPECT_DOUBLE_EQ(0.0, t.phase1D);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude2D);
  EXPECT_DOUBLE_EQ(0.0, t.phase2D);

  // GIVEN (INPUTS)
  tesInterp.index.hgt = 7;  // 2
  tesInterp.index.lat = 7;  // -37.5
  tesInterp.index.ls = 7;   // 210
  tesInterp.index.tyr = 1;

  // RUN
  t = tesInterp.getTideParameters(tesInterp.mgcmTemperature, tesInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(231.75, t.diurnalMean);
  EXPECT_DOUBLE_EQ(17.16, t.amplitude1D);
  EXPECT_DOUBLE_EQ(17.3, t.phase1D);
  EXPECT_DOUBLE_EQ(5.84, t.amplitude2D);
  EXPECT_DOUBLE_EQ(3.4, t.phase2D);

  // RUN
  t = tesInterp.getTideParameters(tesInterp.mgcmDensity, tesInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.079E-02, t.diurnalMean);
  EXPECT_DOUBLE_EQ(10.02, t.amplitude1D);
  EXPECT_DOUBLE_EQ(5.5, t.phase1D);
  EXPECT_DOUBLE_EQ(4.95, t.amplitude2D);
  EXPECT_DOUBLE_EQ(8.1, t.phase2D);

  // RUN
  t = tesInterp.getTideParameters(tesInterp.mgcmEWWind, tesInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-3.0, t.diurnalMean);
  EXPECT_DOUBLE_EQ(17.0, t.amplitude1D);
  EXPECT_DOUBLE_EQ(7.4, t.phase1D);
  EXPECT_DOUBLE_EQ(13.8, t.amplitude2D);
  EXPECT_DOUBLE_EQ(1.8, t.phase2D);

  // RUN
  t = tesInterp.getTideParameters(tesInterp.mgcmNSWind, tesInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.6, t.diurnalMean);
  EXPECT_DOUBLE_EQ(16.9, t.amplitude1D);
  EXPECT_DOUBLE_EQ(13.0, t.phase1D);
  EXPECT_DOUBLE_EQ(12.9, t.amplitude2D);
  EXPECT_DOUBLE_EQ(5.2, t.phase2D);

  // TEAR-DOWN
}


TEST(TesLowerInterpolator, updateIndices)
{
  // SETUP
  TesLowerInterpolator lower;
  Position pos;

  // GIVEN (INPUTS)
  pos.height = 0.0;
  pos.surfaceHeight = 0.0_deg;
  lower.setPosition(pos);

  // RUN
  lower.updateIndices(0);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(6U, lower.getBaseIndex().hgt);

  // GIVEN (INPUTS)
  pos.height = 17.2;
  pos.surfaceHeight = 0.0_km;
  lower.setPosition(pos);

  // RUN
  lower.updateIndices(0);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(16U, lower.getBaseIndex().hgt);

  // GIVEN (INPUTS)
  pos.height = 78.0;
  pos.surfaceHeight = -12.4_km;
  lower.setPosition(pos);

  // RUN
  lower.updateIndices(0);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(28U, lower.getBaseIndex().hgt);

  // TEAR-DOWN
}

TEST(TesLowerInterpolator, update_0)
{
  // SETUP
  TesLowerInterpolator lower;
  Position pos;
  EphemerisState ephem;

  // GIVEN (INPUTS)
  pos.height = 0.0;
  pos.latitude = 0.0_deg;
  pos.setLongitude(0.0, WEST_POSITIVE);
  pos.surfaceHeight = 0.0;
  ephem.longitudeSun = 0.0_deg;
  ephem.solarTime = 0.0;
  lower.setPosition(pos);
  lower.setEphemerisState(ephem);
  lower.setMapYear(1);

  // RUN
  lower.updateIndices(0);
  lower.update();
  const AtmosphereState& atmos = lower.getAtmosphereState();
  const MarsAtmosphereState& atmosMars = atmos.getPlanetSpecificMetrics<MarsAtmosphereState>();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(221.32324545094451, atmos.temperature);
  EXPECT_DOUBLE_EQ(225.77, atmosMars.temperatureDaily);
  EXPECT_DOUBLE_EQ(215.17147683165450, atmosMars.temperatureMin);
  EXPECT_DOUBLE_EQ(239.84554397455946, atmosMars.temperatureMax);
  EXPECT_DOUBLE_EQ(536.81456136233157, atmos.pressure);
  EXPECT_DOUBLE_EQ(538.5, atmosMars.pressureDaily);
  EXPECT_DOUBLE_EQ(1.2792588325664331E-002, atmos.density);
  EXPECT_DOUBLE_EQ(1.2579999999999999E-002, atmosMars.densityDaily);
  EXPECT_DOUBLE_EQ(1.1663694818155912E-002, atmosMars.densityMin);
  EXPECT_DOUBLE_EQ(1.3327026689420681E-002, atmosMars.densityMax);
  EXPECT_DOUBLE_EQ(-7.1234474135167574, atmos.ewWind);
  EXPECT_DOUBLE_EQ(-3.3999999999999999, atmosMars.ewWindDaily);
  EXPECT_DOUBLE_EQ(-7.8286144411851444, atmos.nsWind);
  EXPECT_DOUBLE_EQ(-0.85, atmosMars.nsWindDaily);

  // TEAR-DOWN
}

TEST(TesLowerInterpolator, update_1)
{
  // SETUP
  TesLowerInterpolator lower;
  Position pos;
  EphemerisState ephem;

  // GIVEN (INPUTS)
  pos.height = 0.0;
  pos.latitude = 0.0_deg;
  pos.setLongitude(0.0, WEST_POSITIVE);
  pos.surfaceHeight = 0.0;
  ephem.longitudeSun = 0.0_deg;
  ephem.solarTime = 0.0;
  lower.setPosition(pos);
  lower.setEphemerisState(ephem);
  lower.setMapYear(1);

  // RUN
  lower.updateIndices(1);
  lower.update();
  const AtmosphereState& atmos = lower.getAtmosphereState();
  const MarsAtmosphereState& atmosMars = atmos.getPlanetSpecificMetrics<MarsAtmosphereState>();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(219.88901373042219, atmos.temperature);
  EXPECT_DOUBLE_EQ(224.24, atmosMars.temperatureDaily);
  EXPECT_DOUBLE_EQ(214.26923577892066, atmosMars.temperatureMin);
  EXPECT_DOUBLE_EQ(237.18928936030582, atmosMars.temperatureMax);
  EXPECT_DOUBLE_EQ(493.30981085197459, atmos.pressure);
  EXPECT_DOUBLE_EQ(493.5, atmosMars.pressureDaily);
  EXPECT_DOUBLE_EQ(1.1835166412634409E-002, atmos.density);
  EXPECT_DOUBLE_EQ(1.1610000000000001E-002, atmosMars.densityDaily);
  EXPECT_DOUBLE_EQ(1.0805990482265781E-002, atmosMars.densityMin);
  EXPECT_DOUBLE_EQ(1.2247176903399284E-002, atmosMars.densityMax);
  EXPECT_DOUBLE_EQ(-7.2682751489577804, atmos.ewWind);
  EXPECT_DOUBLE_EQ(-4.0, atmosMars.ewWindDaily);
  EXPECT_DOUBLE_EQ(-5.2480287909782941, atmos.nsWind);
  EXPECT_DOUBLE_EQ(0.59999999999999998, atmosMars.nsWindDaily);

  // TEAR-DOWN
}

} // namespace
