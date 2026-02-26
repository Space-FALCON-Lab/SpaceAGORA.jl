#include "unittest.h"
#include "MGCMLowerInterpolator.h"

namespace GRAM {

TEST(MGCMLowerInterpolator, initializeLowerData_getTideParameters)
{
  // SETUP
  MGCMLowerInterpolator marsInterp;

  // GIVEN (INPUTS)
  greal poleFactor = 1.0;
  marsInterp.index.hgt = 0;
  marsInterp.index.lat = 0;
  marsInterp.index.lon = 0;
  marsInterp.index.ls = 0;
  marsInterp.index.od = 1;

  // RUN
  MarsTideParameters t = marsInterp.getTideParameters(marsInterp.mgcmTemperature, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(165.49, t.diurnalMean);
  EXPECT_DOUBLE_EQ(1.28, t.amplitude1D);
  EXPECT_DOUBLE_EQ(1.0, t.phase1D);
  EXPECT_DOUBLE_EQ(0.55, t.amplitude2D);
  EXPECT_DOUBLE_EQ(0.4, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mgcmPressure, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(505.12958724133, t.diurnalMean);
  EXPECT_DOUBLE_EQ(0.28, t.amplitude1D);
  EXPECT_DOUBLE_EQ(12.9, t.phase1D);
  EXPECT_DOUBLE_EQ(0.14, t.amplitude2D);
  EXPECT_DOUBLE_EQ(2.1, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mgcmEWWind, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, t.diurnalMean);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude1D);
  EXPECT_DOUBLE_EQ(0.0, t.phase1D);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude2D);
  EXPECT_DOUBLE_EQ(0.0, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mgcmNSWind, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, t.diurnalMean);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude1D);
  EXPECT_DOUBLE_EQ(0.0, t.phase1D);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude2D);
  EXPECT_DOUBLE_EQ(0.0, t.phase2D);

  // GIVEN (INPUTS)
  marsInterp.index.hgt = 7;  // 35
  marsInterp.index.lat = 7;  // -37.5
  marsInterp.index.lon = 7;  // 63
  marsInterp.index.ls = 7;   // 210
  marsInterp.index.od = 2;

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mgcmTemperature, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(199.57, t.diurnalMean);
  EXPECT_DOUBLE_EQ(4.20, t.amplitude1D);
  EXPECT_DOUBLE_EQ(17.0, t.phase1D);
  EXPECT_DOUBLE_EQ(1.65, t.amplitude2D);
  EXPECT_DOUBLE_EQ(8.2, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mgcmPressure, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(24.626526935150132, t.diurnalMean);
  EXPECT_DOUBLE_EQ(12.25, t.amplitude1D);
  EXPECT_DOUBLE_EQ(17.9, t.phase1D);
  EXPECT_DOUBLE_EQ(1.68, t.amplitude2D);
  EXPECT_DOUBLE_EQ(11.4, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mgcmEWWind, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(5.2, t.diurnalMean);
  EXPECT_DOUBLE_EQ(14.7, t.amplitude1D);
  EXPECT_DOUBLE_EQ(19.7, t.phase1D);
  EXPECT_DOUBLE_EQ(4.4, t.amplitude2D);
  EXPECT_DOUBLE_EQ(3.6, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mgcmNSWind, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.2, t.diurnalMean);
  EXPECT_DOUBLE_EQ(57.4, t.amplitude1D);
  EXPECT_DOUBLE_EQ(0.2, t.phase1D);
  EXPECT_DOUBLE_EQ(7.1, t.amplitude2D);
  EXPECT_DOUBLE_EQ(2.4, t.phase2D);

  // TEAR-DOWN
}

TEST(MGCMLowerInterpolator, updateIndices)
{
  // SETUP
  MGCMLowerInterpolator lower;
  Position pos;

  // GIVEN (INPUTS)
  pos.height = 0.0_km;
  pos.surfaceHeight = 0.0_km;
  lower.setPosition(pos);

  // RUN
  lower.updateIndices(0);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(1U, lower.getBaseIndex().hgt);

  // GIVEN (INPUTS)
  pos.height = 17.2;
  pos.surfaceHeight = 0.0_km;
  lower.setPosition(pos);

  // RUN
  lower.updateIndices(0);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(3U, lower.getBaseIndex().hgt);

  // GIVEN (INPUTS)
  pos.height = 78.0;
  pos.surfaceHeight = -12.4_km;
  lower.setPosition(pos);

  // RUN
  lower.updateIndices(0);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(15U, lower.getBaseIndex().hgt);

  // TEAR-DOWN
}

TEST(MGCMLowerInterpolator, update_0)
{
  // SETUP
  MGCMLowerInterpolator lower;
  Position pos;
  EphemerisState ephem;

  // GIVEN (INPUTS)
  pos.height = 0.0;
  pos.latitude = 0.0_deg;
  pos.setLongitude(0.0, WEST_POSITIVE);
  pos.surfaceHeight = 0.0_km;
  ephem.longitudeSun = 0.0_deg;
  ephem.solarTime = 0.0;
  lower.setPosition(pos);
  lower.setEphemerisState(ephem);
  lower.setDustOpticalDepth(0.5);

  // RUN
  lower.updateIndices(0);
  lower.update();
  const AtmosphereState& atmos = lower.getAtmosphereState();
  const MarsAtmosphereState& atmosMars = atmos.getPlanetSpecificMetrics<MarsAtmosphereState>();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(215.97899785531393, atmos.temperature);
  EXPECT_DOUBLE_EQ(218.01355846628604, atmosMars.temperatureDaily);
  EXPECT_DOUBLE_EQ(213.77971072053151, atmosMars.temperatureMin);
  EXPECT_DOUBLE_EQ(223.11823626655752, atmosMars.temperatureMax);
  EXPECT_DOUBLE_EQ(362.73207886626744, atmos.pressure);
  EXPECT_DOUBLE_EQ(359.59986420757133, atmosMars.pressureDaily);
  EXPECT_DOUBLE_EQ(8.8892780233347304E-003, atmos.density);
  EXPECT_DOUBLE_EQ(8.7302777692353142E-003, atmosMars.densityDaily);
  EXPECT_DOUBLE_EQ(8.4773779357293806E-003, atmosMars.densityMin);
  EXPECT_DOUBLE_EQ(8.8894407330828446E-003, atmosMars.densityMax);
  EXPECT_DOUBLE_EQ(-1.5171479248005570, atmos.ewWind);
  EXPECT_DOUBLE_EQ(3.0594101885511291, atmosMars.ewWindDaily);
  EXPECT_DOUBLE_EQ(-0.20171139925298331, atmos.nsWind);
  EXPECT_DOUBLE_EQ(8.7858321246722493E-002, atmosMars.nsWindDaily);

  // TEAR-DOWN
}

TEST(MGCMLowerInterpolator, update_1)
{
  // SETUP
  MGCMLowerInterpolator lower;
  Position pos;
  EphemerisState ephem;

  // GIVEN (INPUTS)
  pos.height = 0.0;
  pos.latitude = 0.0_deg;
  pos.setLongitude(0.0, WEST_POSITIVE);
  pos.surfaceHeight = 0.0_km;
  ephem.longitudeSun = 0.0_deg;
  ephem.solarTime = 0.0;
  lower.setPosition(pos);
  lower.setEphemerisState(ephem);
  lower.setDustOpticalDepth(0.5);

  // RUN
  lower.updateIndices(1);
  lower.update();
  const AtmosphereState& atmos = lower.getAtmosphereState();
  const MarsAtmosphereState& atmosMars = atmos.getPlanetSpecificMetrics<MarsAtmosphereState>();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(205.99739548049979, atmos.temperature);
  EXPECT_DOUBLE_EQ(209.83137327412436, atmosMars.temperatureDaily);
  EXPECT_DOUBLE_EQ(205.63795642680051, atmosMars.temperatureMin);
  EXPECT_DOUBLE_EQ(214.00678038172458, atmosMars.temperatureMax);
  EXPECT_DOUBLE_EQ(228.63348013738204, atmos.pressure);
  EXPECT_DOUBLE_EQ(228.47610087406252, atmosMars.pressureDaily);
  EXPECT_DOUBLE_EQ(5.8734206580971947E-003, atmos.density);
  EXPECT_DOUBLE_EQ(5.7621341427331252E-003, atmosMars.densityDaily);
  EXPECT_DOUBLE_EQ(5.6667014485776972E-003, atmosMars.densityMin);
  EXPECT_DOUBLE_EQ(5.8862484932708584E-003, atmosMars.densityMax);
  EXPECT_DOUBLE_EQ(16.889257651973274, atmos.ewWind);
  EXPECT_DOUBLE_EQ(15.431029595741352, atmosMars.ewWindDaily);
  EXPECT_DOUBLE_EQ(-2.6194998195367525E-002, atmos.nsWind);
  EXPECT_DOUBLE_EQ(-0.17421581762073929, atmosMars.nsWindDaily);

  // TEAR-DOWN
}

} // namespace