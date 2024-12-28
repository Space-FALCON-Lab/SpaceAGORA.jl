#include "unittest.h"
#include "MGCMUpperInterpolator.h"

namespace GRAM {

TEST(MGCMUpperInterpolator, initializeUpperData_getTideParameters)
{
  // SETUP
  MGCMUpperInterpolator marsInterp;  // This calls initializeUpperData();

  // GIVEN (INPUTS)
  greal poleFactor = 1.0;
  marsInterp.index.hgt = 0;
  marsInterp.index.lat = 0;
  marsInterp.index.lon = 0;
  marsInterp.index.ls = 0;
  marsInterp.index.od = 2;
  marsInterp.index.f107 = 0;

  // RUN
  MarsTideParameters t = marsInterp.getTideParameters(marsInterp.mtgcmTemperature, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(128.36, t.diurnalMean);
  EXPECT_DOUBLE_EQ(0.25, t.amplitude1D);
  EXPECT_DOUBLE_EQ(18.1, t.phase1D);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude2D);
  EXPECT_DOUBLE_EQ(17.8, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mtgcmPressure, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.078714847727580278, t.diurnalMean);
  EXPECT_DOUBLE_EQ(0.1, t.amplitude1D);
  EXPECT_DOUBLE_EQ(15.6, t.phase1D);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude2D);
  EXPECT_DOUBLE_EQ(19.1, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mtgcmDensity, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(3.0550881504635369e-06, t.diurnalMean);
  EXPECT_DOUBLE_EQ(0.1, t.amplitude1D);
  EXPECT_DOUBLE_EQ(31.2, t.phase1D);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude2D);
  EXPECT_DOUBLE_EQ(20.5, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mtgcmEWWind, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.0, t.diurnalMean);
  EXPECT_DOUBLE_EQ(2.5, t.amplitude1D);
  EXPECT_DOUBLE_EQ(17.0, t.phase1D);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude2D);
  EXPECT_DOUBLE_EQ(22.3, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mtgcmNSWind, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.1, t.diurnalMean);
  EXPECT_DOUBLE_EQ(2.3, t.amplitude1D);
  EXPECT_DOUBLE_EQ(23.1, t.phase1D);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude2D);
  EXPECT_DOUBLE_EQ(13.1, t.phase2D);


  // GIVEN (INPUTS)
  poleFactor = 1.0;
  marsInterp.index.hgt = 7;  // 115
  marsInterp.index.lat = 7;  // -52.5
  marsInterp.index.lon = 7;  // 63
  marsInterp.index.ls = 7;   // 210
  marsInterp.index.od = 0;
  marsInterp.index.f107 = 1;

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mtgcmTemperature, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(126.58, t.diurnalMean);
  EXPECT_DOUBLE_EQ(7.82, t.amplitude1D);
  EXPECT_DOUBLE_EQ(16.1, t.phase1D);
  EXPECT_DOUBLE_EQ(3.37, t.amplitude2D);
  EXPECT_DOUBLE_EQ(11.2, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mtgcmPressure, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.00033283134561237841, t.diurnalMean);
  EXPECT_DOUBLE_EQ(30.2, t.amplitude1D);
  EXPECT_DOUBLE_EQ(16.2, t.phase1D);
  EXPECT_DOUBLE_EQ(4.4, t.amplitude2D);
  EXPECT_DOUBLE_EQ(11.7, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mtgcmDensity, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.2882780429650712e-08, t.diurnalMean);
  EXPECT_DOUBLE_EQ(26.9, t.amplitude1D);
  EXPECT_DOUBLE_EQ(16.2, t.phase1D);
  EXPECT_DOUBLE_EQ(2.0, t.amplitude2D);
  EXPECT_DOUBLE_EQ(12.7, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mtgcmEWWind, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-9.8, t.diurnalMean);
  EXPECT_DOUBLE_EQ(64.9, t.amplitude1D);
  EXPECT_DOUBLE_EQ(16.8, t.phase1D);
  EXPECT_DOUBLE_EQ(12.1, t.amplitude2D);
  EXPECT_DOUBLE_EQ(17.4, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.mtgcmNSWind, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-2.3, t.diurnalMean);
  EXPECT_DOUBLE_EQ(80.6, t.amplitude1D);
  EXPECT_DOUBLE_EQ(22.6, t.phase1D);
  EXPECT_DOUBLE_EQ(10.5, t.amplitude2D);
  EXPECT_DOUBLE_EQ(20.6, t.phase2D);

  // TEAR-DOWN
}

TEST(MGCMUpperInterpolator, updateIndices)
{
  // SETUP
  MGCMUpperInterpolator upper;
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

TEST(MGCMUpperInterpolator, update_0)
{
  // SETUP
  MGCMUpperInterpolator upper;
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
  upper.setDustOpticalDepth(0.5);
  upper.setF107(100.0);

  // RUN
  upper.updateIndices(0);
  upper.update();
  const AtmosphereState& atmos = upper.getAtmosphereState();
  const MarsAtmosphereState& atmosMars = atmos.getPlanetSpecificMetrics<MarsAtmosphereState>();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(131.13670477945630, atmos.temperature);
  EXPECT_DOUBLE_EQ(133.32761668705120, atmosMars.temperatureDaily);
  EXPECT_DOUBLE_EQ(127.75389073159300, atmosMars.temperatureMin);
  EXPECT_DOUBLE_EQ(138.04048394837946, atmosMars.temperatureMax);
  EXPECT_DOUBLE_EQ(6.5266278642724065E-002, atmos.pressure);
  EXPECT_DOUBLE_EQ(6.4213082267221208E-002, atmosMars.pressureDaily);
  EXPECT_DOUBLE_EQ(2.4423321115260274E-006, atmos.density);
  EXPECT_DOUBLE_EQ(2.4014693931561993E-006, atmosMars.densityDaily);
  EXPECT_DOUBLE_EQ(2.3017768317715962E-006, atmosMars.densityMin);
  EXPECT_DOUBLE_EQ(2.4987892681721880E-006, atmosMars.densityMax);
  EXPECT_DOUBLE_EQ(-2.2063707576562441, atmos.ewWind);
  EXPECT_DOUBLE_EQ(-1.0909492104842728, atmosMars.ewWindDaily);
  EXPECT_DOUBLE_EQ(-2.9003041177737967, atmos.nsWind);
  EXPECT_DOUBLE_EQ(1.5143328498689004E-002, atmosMars.nsWindDaily);
  EXPECT_DOUBLE_EQ(122.86916265647130, atmosMars.thermosphereBaseHeight);

  // TEAR-DOWN
}

TEST(MGCMUpperInterpolator, update_1)
{
  // SETUP
  MGCMUpperInterpolator upper;
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
  upper.setDustOpticalDepth(0.5);
  upper.setF107(100.0);

  // RUN
  upper.updateIndices(1);
  upper.update();
  const AtmosphereState& atmos = upper.getAtmosphereState();
  const MarsAtmosphereState& atmosMars = atmos.getPlanetSpecificMetrics<MarsAtmosphereState>();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(126.52964086536295, atmos.temperature);
  EXPECT_DOUBLE_EQ(128.83761030835052, atmosMars.temperatureDaily);
  EXPECT_DOUBLE_EQ(122.97932881518649, atmosMars.temperatureMin);
  EXPECT_DOUBLE_EQ(133.71979506247214, atmosMars.temperatureMax);
  EXPECT_DOUBLE_EQ(4.0594431264620617E-002, atmos.pressure);
  EXPECT_DOUBLE_EQ(3.9928049093169285E-002, atmosMars.pressureDaily);
  EXPECT_DOUBLE_EQ(1.5865823719471072E-006, atmos.density);
  EXPECT_DOUBLE_EQ(1.5445898613410370E-006, atmosMars.densityDaily);
  EXPECT_DOUBLE_EQ(1.4424015974793494E-006, atmosMars.densityMin);
  EXPECT_DOUBLE_EQ(1.6424110430664369E-006, atmosMars.densityMax);
  EXPECT_DOUBLE_EQ(-3.1900113029344803, atmos.ewWind);
  EXPECT_DOUBLE_EQ(-2.1336295764931563, atmosMars.ewWindDaily);
  EXPECT_DOUBLE_EQ(-3.5353946311363091, atmos.nsWind);
  EXPECT_DOUBLE_EQ(-0.36971334300262204, atmosMars.nsWindDaily);
  EXPECT_DOUBLE_EQ(122.86916265647130, atmosMars.thermosphereBaseHeight);

  // TEAR-DOWN
}

} // namespace
