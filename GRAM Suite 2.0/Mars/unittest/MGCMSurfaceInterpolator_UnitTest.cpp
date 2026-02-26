#include "unittest.h"
#include "MGCMSurfaceInterpolator.h"

namespace GRAM {

TEST(MGCMSurfaceInterpolator, initializeSurfaceData_getTideParameters)
{
  // SETUP
  MGCMSurfaceInterpolator marsInterp;

  // GIVEN (INPUTS)
  greal poleFactor = 1.0;
  marsInterp.index.hgt = 0;
  marsInterp.index.lat = 0;
  marsInterp.index.lon = 0;
  marsInterp.index.ls = 0;
  marsInterp.index.od = 0;

  // RUN
  MarsTideParameters t = marsInterp.getTideParameters(marsInterp.surfaceTemperatures, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(144.9, t.diurnalMean);
  EXPECT_DOUBLE_EQ(0.02, t.amplitude1D);
  EXPECT_DOUBLE_EQ(15.8, t.phase1D);
  EXPECT_DOUBLE_EQ(0.02, t.amplitude2D);
  EXPECT_DOUBLE_EQ(2.3, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.surfaceEWWinds, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, t.diurnalMean);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude1D);
  EXPECT_DOUBLE_EQ(0.0, t.phase1D);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude2D);
  EXPECT_DOUBLE_EQ(0.0, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.surfaceNSWinds, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, t.diurnalMean);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude1D);
  EXPECT_DOUBLE_EQ(0.0, t.phase1D);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude2D);
  EXPECT_DOUBLE_EQ(0.0, t.phase2D);

  // GIVEN (INPUTS)
  poleFactor = 2.2;
  marsInterp.index.hgt = 0;
  marsInterp.index.lat = 0;
  marsInterp.index.lon = 0;
  marsInterp.index.ls = 0;
  marsInterp.index.od = 0;

  // RUN
  t = marsInterp.getTideParameters(marsInterp.surfaceTemperatures, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(144.9, t.diurnalMean);
  EXPECT_DOUBLE_EQ(0.02 * poleFactor, t.amplitude1D);
  EXPECT_DOUBLE_EQ(15.8, t.phase1D);
  EXPECT_DOUBLE_EQ(0.02 * poleFactor, t.amplitude2D);
  EXPECT_DOUBLE_EQ(2.3, t.phase2D);

  // GIVEN (INPUTS)
  poleFactor = 1.0;
  marsInterp.index.hgt = 0;
  marsInterp.index.lat = 7;  // -37.5
  marsInterp.index.lon = 7;  // 63
  marsInterp.index.ls = 7;   // 210
  marsInterp.index.od = 0;

  // RUN
  t = marsInterp.getTideParameters(marsInterp.surfaceTemperatures, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(229.67, t.diurnalMean);
  EXPECT_DOUBLE_EQ(48.92, t.amplitude1D);
  EXPECT_DOUBLE_EQ(13.3, t.phase1D);
  EXPECT_DOUBLE_EQ(13.46, t.amplitude2D);
  EXPECT_DOUBLE_EQ(0.4, t.phase2D);

  // GIVEN (INPUTS)
  poleFactor = 1.0;
  marsInterp.index.hgt = 1;
  marsInterp.index.lat = 7;  // -37.5
  marsInterp.index.lon = 7;  // 63
  marsInterp.index.ls = 7;   // 210
  marsInterp.index.od = 0;

  // RUN
  t = marsInterp.getTideParameters(marsInterp.surfaceTemperatures, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(219.95, t.diurnalMean);
  EXPECT_DOUBLE_EQ(27.51, t.amplitude1D);
  EXPECT_DOUBLE_EQ(15.0, t.phase1D);
  EXPECT_DOUBLE_EQ(3.56, t.amplitude2D);
  EXPECT_DOUBLE_EQ(2.0, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.surfaceEWWinds, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-1.8, t.diurnalMean);
  EXPECT_DOUBLE_EQ(3.1, t.amplitude1D);
  EXPECT_DOUBLE_EQ(11.2, t.phase1D);
  EXPECT_DOUBLE_EQ(1.1, t.amplitude2D);
  EXPECT_DOUBLE_EQ(9.6, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.surfaceNSWinds, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(3.0, t.diurnalMean);
  EXPECT_DOUBLE_EQ(4.3, t.amplitude1D);
  EXPECT_DOUBLE_EQ(16.1, t.phase1D);
  EXPECT_DOUBLE_EQ(1.4, t.amplitude2D);
  EXPECT_DOUBLE_EQ(1.4, t.phase2D);

  // GIVEN (INPUTS)
  poleFactor = 1.0;
  marsInterp.index.hgt = 2;
  marsInterp.index.lat = 7;  // -37.5
  marsInterp.index.lon = 7;  // 63
  marsInterp.index.ls = 7;   // 210
  marsInterp.index.od = 2;

  // RUN
  t = marsInterp.getTideParameters(marsInterp.surfaceTemperatures, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(224.51, t.diurnalMean);
  EXPECT_DOUBLE_EQ(21.22, t.amplitude1D);
  EXPECT_DOUBLE_EQ(15.7, t.phase1D);
  EXPECT_DOUBLE_EQ(6.35, t.amplitude2D);
  EXPECT_DOUBLE_EQ(2.5, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.surfaceEWWinds, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-7.5, t.diurnalMean);
  EXPECT_DOUBLE_EQ(5.2, t.amplitude1D);
  EXPECT_DOUBLE_EQ(8.7, t.phase1D);
  EXPECT_DOUBLE_EQ(4.4, t.amplitude2D);
  EXPECT_DOUBLE_EQ(2.2, t.phase2D);

  // RUN
  t = marsInterp.getTideParameters(marsInterp.surfaceNSWinds, marsInterp.index.lat, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(3.9, t.diurnalMean);
  EXPECT_DOUBLE_EQ(8.0, t.amplitude1D);
  EXPECT_DOUBLE_EQ(14.2, t.phase1D);
  EXPECT_DOUBLE_EQ(2.4, t.amplitude2D);
  EXPECT_DOUBLE_EQ(4.5, t.phase2D);

  // TEAR-DOWN
}

TEST(MGCMSurfaceInterpolator, updateIndices)
{
  // SETUP
  MGCMSurfaceInterpolator surface;
  Position pos;

  // GIVEN (INPUTS)
  pos.height = 0.0_km;
  pos.longitude = 0.0_deg;
  pos.surfaceHeight = 0.0_km;
  surface.setPosition(pos);

  // RUN
  surface.updateIndices(0);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(0U, surface.getBaseIndex().hgt);
  EXPECT_EQ(0U, surface.getBaseIndex().lon);
  EXPECT_EQ(1U, surface.oneKilometerIndex);
  EXPECT_DOUBLE_EQ(0.0, surface.getDisplacements().lon);

  // GIVEN (INPUTS)
  pos.height = 17.22_km;
  pos.longitude = 0.0_deg;
  pos.surfaceHeight = 17.2_km;
  surface.setPosition(pos);

  // RUN
  surface.updateIndices(0);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(1U, surface.getBaseIndex().hgt);
  EXPECT_EQ(0U, surface.getBaseIndex().lon);
  EXPECT_EQ(4U, surface.oneKilometerIndex);
  EXPECT_DOUBLE_EQ(0.0, surface.getDisplacements().lon);

  // GIVEN (INPUTS)
  pos.height = -2.0_km;
  pos.setLongitude(100.0, WEST_POSITIVE);
  pos.surfaceHeight = -12.4_km;
  surface.setPosition(pos);

  // RUN
  surface.updateIndices(0);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(2U, surface.getBaseIndex().hgt);
  EXPECT_EQ(11U, surface.getBaseIndex().lon);
  EXPECT_EQ(0U, surface.oneKilometerIndex);
  EXPECT_NEAR(1.0 / 9.0, surface.getDisplacements().lon, 1.0e-12);

  // TEAR-DOWN
}

TEST(MGCMSurfaceInterpolator, update_0)
{
  // SETUP
  MGCMSurfaceInterpolator surface;
  Position pos;
  EphemerisState ephem;

  // GIVEN (INPUTS)
  pos.height = 0.0_km;
  pos.latitude = 0.0_deg;
  pos.setLongitude(0.0, WEST_POSITIVE);
  pos.surfaceHeight = -0.0_km;
  ephem.longitudeSun = 0.0_deg;
  ephem.solarTime = 0.0;
  surface.setPosition(pos);
  surface.setEphemerisState(ephem);
  surface.setDustOpticalDepth(0.5);

  // RUN
  surface.updateIndices(0);
  surface.update();
  const AtmosphereState& atmos = surface.getAtmosphereState();
  const MarsAtmosphereState& atmosMars = atmos.getPlanetSpecificMetrics<MarsAtmosphereState>();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(176.57493189426651, atmos.temperature);
  EXPECT_DOUBLE_EQ(210.28040822879183, atmosMars.temperatureDaily);
  EXPECT_DOUBLE_EQ(164.29513178020713, atmosMars.temperatureMin);
  EXPECT_DOUBLE_EQ(286.92924879060388, atmosMars.temperatureMax);
  EXPECT_DOUBLE_EQ(849.70516311633094, atmos.pressure);
  EXPECT_DOUBLE_EQ(552.38511859284733, atmosMars.pressureDaily);
  EXPECT_DOUBLE_EQ(2.5470135541849903E-002, atmos.density);
  EXPECT_DOUBLE_EQ(1.3903854197821050E-002, atmosMars.densityDaily);
  EXPECT_DOUBLE_EQ(9.4343517992287428E-003, atmosMars.densityMin);
  EXPECT_DOUBLE_EQ(1.8859977732929582E-002, atmosMars.densityMax);
  EXPECT_DOUBLE_EQ(0.0, atmos.ewWind);
  EXPECT_DOUBLE_EQ(0.0, atmosMars.ewWindDaily);
  EXPECT_DOUBLE_EQ(0.0, atmos.nsWind);
  EXPECT_DOUBLE_EQ(0.0, atmosMars.nsWindDaily);
  EXPECT_DOUBLE_EQ(5.8738880784051792, atmos.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(11.648127132752066, surface.getPressureScaleHeightDaily());
  EXPECT_DOUBLE_EQ(6.2424084552087606, atmos.densityScaleHeight);
  EXPECT_DOUBLE_EQ(0.0, surface.temperatureAtFiveMeters);

  // TEAR-DOWN
}

TEST(MGCMSurfaceInterpolator, update_1)
{
  // SETUP
  MGCMSurfaceInterpolator surface;
  Position pos;
  EphemerisState ephem;

  // GIVEN (INPUTS)
  pos.height = 0.0_km;
  pos.latitude = 0.0_deg;
  pos.setLongitude(0.0, WEST_POSITIVE);
  pos.surfaceHeight = -0.0_km;
  ephem.longitudeSun = 0.0_deg;
  ephem.solarTime = 0.0;
  surface.setPosition(pos);
  surface.setEphemerisState(ephem);
  surface.setDustOpticalDepth(0.5);

  // RUN
  surface.updateIndices(1);
  surface.update();
  const AtmosphereState& atmos = surface.getAtmosphereState();
  const MarsAtmosphereState& atmosMars = atmos.getPlanetSpecificMetrics<MarsAtmosphereState>();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(179.78523768970709, atmos.temperature);
  EXPECT_DOUBLE_EQ(202.40870179478475, atmosMars.temperatureDaily);
  EXPECT_DOUBLE_EQ(171.03140318242191, atmosMars.temperatureMin);
  EXPECT_DOUBLE_EQ(243.56032813583516, atmosMars.temperatureMax);
  EXPECT_DOUBLE_EQ(576.93862072597926, atmos.pressure);
  EXPECT_DOUBLE_EQ(552.14805619963079, atmosMars.pressureDaily);
  EXPECT_DOUBLE_EQ(1.6985081288350235E-002, atmos.density);
  EXPECT_DOUBLE_EQ(1.4438378231578182E-002, atmosMars.densityDaily);
  EXPECT_DOUBLE_EQ(1.1488327666650384E-002, atmosMars.densityMin);
  EXPECT_DOUBLE_EQ(1.7962432908810805E-002, atmosMars.densityMax);
  EXPECT_DOUBLE_EQ(-1.3699401456718585, atmos.ewWind);
  EXPECT_DOUBLE_EQ(-0.83628995648820115, atmosMars.ewWindDaily);
  EXPECT_DOUBLE_EQ(6.1954649534615518, atmos.nsWind);
  EXPECT_DOUBLE_EQ(5.2287182922388569, atmosMars.nsWindDaily);
  EXPECT_DOUBLE_EQ(10.763430000653845, atmos.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(11.648127132752066, surface.getPressureScaleHeightDaily());
  EXPECT_DOUBLE_EQ(11.438714110019593, atmos.densityScaleHeight);
  EXPECT_DOUBLE_EQ(179.78523768970709, surface.temperatureAtFiveMeters);

  // TEAR-DOWN
}

} // namespace