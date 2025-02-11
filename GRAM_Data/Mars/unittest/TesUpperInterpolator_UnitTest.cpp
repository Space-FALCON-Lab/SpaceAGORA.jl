#include "unittest.h"
#include "TesUpperInterpolator.h"

namespace GRAM {

TEST(TesUpperInterpolator, initializeUpperData_getTideParameters)
{
  // SETUP
  TesUpperInterpolator marsInterp;

  // GIVEN (INPUTS)
  greal poleFactor = 1.0;
  marsInterp.index.hgt = 0;
  marsInterp.index.lat = 0;
  marsInterp.index.ls = 0;
  marsInterp.index.tyr = 0;
  marsInterp.index.f107 = 0;

  // RUN
  MarsTideParameters t = marsInterp.getTideParameters(marsInterp.mtgcmTemperature, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(149.64, t.diurnalMean);
  EXPECT_DOUBLE_EQ(2.16, t.amplitude1D);
  EXPECT_DOUBLE_EQ(21.3, t.phase1D);
  EXPECT_DOUBLE_EQ(0.2, t.amplitude2D);
  EXPECT_DOUBLE_EQ(13.2, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mtgcmPressure, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(2.082E-02, t.diurnalMean);
  EXPECT_DOUBLE_EQ(3.37, t.amplitude1D);
  EXPECT_DOUBLE_EQ(14.5, t.phase1D);
  EXPECT_DOUBLE_EQ(1.01, t.amplitude2D);
  EXPECT_DOUBLE_EQ(13.4, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mtgcmDensity, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(7.268E-07, t.diurnalMean);
  EXPECT_DOUBLE_EQ(3.85, t.amplitude1D);
  EXPECT_DOUBLE_EQ(12.7, t.phase1D);
  EXPECT_DOUBLE_EQ(0.93, t.amplitude2D);
  EXPECT_DOUBLE_EQ(13.0, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mtgcmEWWind, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(15.5, t.diurnalMean);
  EXPECT_DOUBLE_EQ(34.4, t.amplitude1D);
  EXPECT_DOUBLE_EQ(17.0, t.phase1D);
  EXPECT_DOUBLE_EQ(0.4, t.amplitude2D);
  EXPECT_DOUBLE_EQ(20.8, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mtgcmNSWind, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(2.5, t.diurnalMean);
  EXPECT_DOUBLE_EQ(36.0, t.amplitude1D);
  EXPECT_DOUBLE_EQ(22.2, t.phase1D);
  EXPECT_DOUBLE_EQ(0.2, t.amplitude2D);
  EXPECT_DOUBLE_EQ(12.3, t.phase2D);


  // GIVEN (INPUTS)
  poleFactor = 1.0;
  marsInterp.index.hgt = 7;  // 115
  marsInterp.index.lat = 7;  // -52.5
  marsInterp.index.ls = 7;   // 210
  marsInterp.index.tyr = 1;
  marsInterp.index.f107 = 1;

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mtgcmTemperature, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(123.11, t.diurnalMean);
  EXPECT_DOUBLE_EQ(8.59, t.amplitude1D);
  EXPECT_DOUBLE_EQ(15.3, t.phase1D);
  EXPECT_DOUBLE_EQ(9.57, t.amplitude2D);
  EXPECT_DOUBLE_EQ(17.4, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mtgcmPressure, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(2.018E-03, t.diurnalMean);
  EXPECT_DOUBLE_EQ(13.36, t.amplitude1D);
  EXPECT_DOUBLE_EQ(10.3, t.phase1D);
  EXPECT_DOUBLE_EQ(6.18, t.amplitude2D);
  EXPECT_DOUBLE_EQ(16.7, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mtgcmDensity, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(8.622E-08, t.diurnalMean);
  EXPECT_DOUBLE_EQ(13.32, t.amplitude1D);
  EXPECT_DOUBLE_EQ(32.8, t.phase1D);
  EXPECT_DOUBLE_EQ(2.7, t.amplitude2D);
  EXPECT_DOUBLE_EQ(13.4, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mtgcmEWWind, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-114.0, t.diurnalMean);
  EXPECT_DOUBLE_EQ(57.9, t.amplitude1D);
  EXPECT_DOUBLE_EQ(17.1, t.phase1D);
  EXPECT_DOUBLE_EQ(43.7, t.amplitude2D);
  EXPECT_DOUBLE_EQ(14.7, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mtgcmNSWind, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(6.2, t.diurnalMean);
  EXPECT_DOUBLE_EQ(63.6, t.amplitude1D);
  EXPECT_DOUBLE_EQ(24.6, t.phase1D);
  EXPECT_DOUBLE_EQ(17.3, t.amplitude2D);
  EXPECT_DOUBLE_EQ(17.4, t.phase2D);

  // TEAR-DOWN
}

TEST(TesUpperInterpolator, updateIndices)
{
  // SETUP
  TesUpperInterpolator upper;
  Position pos;

  // GIVEN (INPUTS)
  pos.height = 80.0;
  pos.latitude = 0.0_deg;
  upper.setPosition(pos);
  upper.setF107(100.0);

  // RUN
  upper.updateIndices(0);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(0U, upper.getBaseIndex().hgt);
  EXPECT_EQ(0U, upper.getBaseIndex().f107);
  EXPECT_DOUBLE_EQ(0.57617504528836971, upper.getDisplacements().f107);

  // GIVEN (INPUTS)
  pos.height = 107.2;
  pos.latitude = 40.0_deg;
  upper.setPosition(pos);
  upper.setF107(100.0);

  // RUN
  upper.updateIndices(0);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(5U, upper.getBaseIndex().hgt);
  EXPECT_EQ(0U, upper.getBaseIndex().f107);
  EXPECT_DOUBLE_EQ(0.57617504528836971, upper.getDisplacements().f107);

  // GIVEN (INPUTS)
  pos.height = 166.0;
  pos.latitude = -55.0_deg;
  upper.setPosition(pos);
  upper.setF107(100.0);

  // RUN
  upper.updateIndices(0);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(17U, upper.getBaseIndex().hgt);
  EXPECT_EQ(0U, upper.getBaseIndex().f107);
  EXPECT_DOUBLE_EQ(0.57617504528836971, upper.getDisplacements().f107);

  // TEAR-DOWN
}

TEST(TesUpperInterpolator, update_0)
{
  // SETUP
  TesUpperInterpolator upper;
  Position pos;
  EphemerisState ephem;

  // GIVEN (INPUTS)
  pos.height = 80.0;
  pos.latitude = -2.5_deg;
  pos.setLongitude(0.0, WEST_POSITIVE);
  ephem.longitudeSun = 0.0_deg;
  ephem.solarTime = 0.0;
  upper.setPosition(pos);
  upper.setEphemerisState(ephem);
  upper.setF107(100.0);
  upper.setMapYear(1);

  // RUN
  upper.updateIndices(0);
  upper.update();
  const AtmosphereState& atmos = upper.getAtmosphereState();
  const MarsAtmosphereState& atmosMars = atmos.getPlanetSpecificMetrics<MarsAtmosphereState>();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(140.30870264119750, atmos.temperature);
  EXPECT_DOUBLE_EQ(138.13728525135866, atmosMars.temperatureDaily);
  EXPECT_DOUBLE_EQ(118.45797373861039, atmosMars.temperatureMin);
  EXPECT_DOUBLE_EQ(163.80438053991955, atmosMars.temperatureMax);
  EXPECT_DOUBLE_EQ(3.3927523018682626E-002, atmos.pressure);
  EXPECT_DOUBLE_EQ(2.9685761339134954E-002, atmosMars.pressureDaily);
  EXPECT_DOUBLE_EQ(1.1889050573952355E-006, atmos.density);
  EXPECT_DOUBLE_EQ(1.1395759814047137E-006, atmosMars.densityDaily);
  EXPECT_DOUBLE_EQ(9.8897470879355332E-007, atmosMars.densityMin);
  EXPECT_DOUBLE_EQ(1.3426948944760350E-006, atmosMars.densityMax);
  EXPECT_DOUBLE_EQ(-72.057124791177202, atmos.ewWind);
  EXPECT_DOUBLE_EQ(-69.445705027173020, atmosMars.ewWindDaily);
  EXPECT_DOUBLE_EQ(17.913337311037893, atmos.nsWind);
  EXPECT_DOUBLE_EQ(3.4423824954711630, atmosMars.nsWindDaily);
  EXPECT_DOUBLE_EQ(117.44226674438839, atmosMars.thermosphereBaseHeight);

  // TEAR-DOWN
}

TEST(TesUpperInterpolator, update_1)
{
  // SETUP
  TesUpperInterpolator upper;
  Position pos;
  EphemerisState ephem;
  MarsInputParameters params;

  // GIVEN (INPUTS)
  pos.height = 80.0;
  pos.latitude = -2.5_deg;
  pos.setLongitude(0.0, WEST_POSITIVE);
  ephem.longitudeSun = 0.0_deg;
  ephem.solarTime = 0.0;
  upper.setPosition(pos);
  upper.setEphemerisState(ephem);
  upper.setF107(100.0);
  upper.setMapYear(1);

  // RUN
  upper.updateIndices(1);
  upper.update();
  const AtmosphereState& atmos = upper.getAtmosphereState();
  const MarsAtmosphereState& atmosMars = atmos.getPlanetSpecificMetrics<MarsAtmosphereState>();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(141.62379019880618, atmos.temperature);
  EXPECT_DOUBLE_EQ(134.50238249547115, atmosMars.temperatureDaily);
  EXPECT_DOUBLE_EQ(118.02074545785513, atmosMars.temperatureMin);
  EXPECT_DOUBLE_EQ(154.35546186247961, atmosMars.temperatureMax);
  EXPECT_DOUBLE_EQ(1.7701322500574365E-002, atmos.pressure);
  EXPECT_DOUBLE_EQ(1.5205760947422807E-002, atmosMars.pressureDaily);
  EXPECT_DOUBLE_EQ(6.7081616119051387E-007, atmos.density);
  EXPECT_DOUBLE_EQ(5.9580317817485022E-007, atmosMars.densityDaily);
  EXPECT_DOUBLE_EQ(4.9468347075469203E-007, atmosMars.densityMin);
  EXPECT_DOUBLE_EQ(7.2852652802446409E-007, atmosMars.densityMax);
  EXPECT_DOUBLE_EQ(-90.953109159516060, atmos.ewWind);
  EXPECT_DOUBLE_EQ(-78.733792549817196, atmosMars.ewWindDaily);
  EXPECT_DOUBLE_EQ(29.843470105117433, atmos.nsWind);
  EXPECT_DOUBLE_EQ(3.8847649909423261, atmosMars.nsWindDaily);
  EXPECT_DOUBLE_EQ(117.44226674438839, atmosMars.thermosphereBaseHeight);

  // TEAR-DOWN
}

} // namespace
