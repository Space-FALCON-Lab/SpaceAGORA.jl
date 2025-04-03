#include "unittest.h"
#include "TesSurfaceInterpolator.h"

namespace GRAM {

TEST(TesSurfaceInterpolator, initializeSurfaceData_getTideParameters)
{
  // SETUP
  TesSurfaceInterpolator tesInterp;

  // GIVEN (INPUTS)
  greal poleFactor = 1.0;
  tesInterp.index.hgt = 0;
  tesInterp.index.lat = 0;
  tesInterp.index.lon = 0;
  tesInterp.index.ls = 0;
  tesInterp.index.tyr = 0;

  // RUN
  MarsTideParameters t = tesInterp.getTideParameters(tesInterp.surfaceTemperatures, tesInterp.index.lat, tesInterp.index.lon, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(143.93, t.diurnalMean);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude1D);
  EXPECT_DOUBLE_EQ(13.8, t.phase1D);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude2D);
  EXPECT_DOUBLE_EQ(1.9, t.phase2D);

  // RUN
  t = tesInterp.getTideParameters(tesInterp.surfaceEWWinds, tesInterp.index.lat, tesInterp.index.lon, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, t.diurnalMean);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude1D);
  EXPECT_DOUBLE_EQ(0.0, t.phase1D);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude2D);
  EXPECT_DOUBLE_EQ(0.0, t.phase2D);

  // RUN
  t = tesInterp.getTideParameters(tesInterp.surfaceNSWinds, tesInterp.index.lat, tesInterp.index.lon, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(0.0, t.diurnalMean);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude1D);
  EXPECT_DOUBLE_EQ(0.0, t.phase1D);
  EXPECT_DOUBLE_EQ(0.0, t.amplitude2D);
  EXPECT_DOUBLE_EQ(0.0, t.phase2D);

  // GIVEN (INPUTS)
  poleFactor = 1.0;
  tesInterp.index.hgt = 0;
  tesInterp.index.lat = 7;  // -37.5
  tesInterp.index.lon = 7;  // 63
  tesInterp.index.ls = 7;   // 210
  tesInterp.index.tyr = 0;

  // RUN
  t = tesInterp.getTideParameters(tesInterp.surfaceTemperatures, tesInterp.index.lat, tesInterp.index.lon, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(231.99, t.diurnalMean);
  EXPECT_DOUBLE_EQ(48.56, t.amplitude1D);
  EXPECT_DOUBLE_EQ(13.5, t.phase1D);
  EXPECT_DOUBLE_EQ(13.37, t.amplitude2D);
  EXPECT_DOUBLE_EQ(0.5, t.phase2D);

  // GIVEN (INPUTS)
  poleFactor = 2.2;
  tesInterp.index.hgt = 0;
  tesInterp.index.lat = 7;  // -37.5
  tesInterp.index.lon = 7;  // 63
  tesInterp.index.ls = 7;   // 210
  tesInterp.index.tyr = 0;

  // RUN
  t = tesInterp.getTideParameters(tesInterp.surfaceTemperatures, tesInterp.index.lat, tesInterp.index.lon, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(231.99, t.diurnalMean);
  EXPECT_DOUBLE_EQ(48.56 * poleFactor, t.amplitude1D);
  EXPECT_DOUBLE_EQ(13.5, t.phase1D);
  EXPECT_DOUBLE_EQ(13.37 * poleFactor, t.amplitude2D);
  EXPECT_DOUBLE_EQ(0.5, t.phase2D);

  // GIVEN (INPUTS)
  poleFactor = 1.0;
  tesInterp.index.hgt = 1;
  tesInterp.index.lat = 7;  // -37.5
  tesInterp.index.lon = 7;  // 63
  tesInterp.index.ls = 7;   // 210
  tesInterp.index.tyr = 0;

  // RUN
  t = tesInterp.getTideParameters(tesInterp.surfaceTemperatures, tesInterp.index.lat, tesInterp.index.lon, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(222.15, t.diurnalMean);
  EXPECT_DOUBLE_EQ(27.49, t.amplitude1D);
  EXPECT_DOUBLE_EQ(15.3, t.phase1D);
  EXPECT_DOUBLE_EQ(3.65, t.amplitude2D);
  EXPECT_DOUBLE_EQ(2.0, t.phase2D);

  // RUN
  t = tesInterp.getTideParameters(tesInterp.surfaceEWWinds, tesInterp.index.lat, tesInterp.index.lon, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.9, t.diurnalMean);
  EXPECT_DOUBLE_EQ(8.3, t.amplitude1D);
  EXPECT_DOUBLE_EQ(4.3, t.phase1D);
  EXPECT_DOUBLE_EQ(1.4, t.amplitude2D);
  EXPECT_DOUBLE_EQ(24.0, t.phase2D);

  // RUN
  t = tesInterp.getTideParameters(tesInterp.surfaceNSWinds, tesInterp.index.lat, tesInterp.index.lon, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-1.3, t.diurnalMean);
  EXPECT_DOUBLE_EQ(4.1, t.amplitude1D);
  EXPECT_DOUBLE_EQ(11.2, t.phase1D);
  EXPECT_DOUBLE_EQ(1.3, t.amplitude2D);
  EXPECT_DOUBLE_EQ(18.2, t.phase2D);

  // GIVEN (INPUTS)
  poleFactor = 1.0;
  tesInterp.index.hgt = 2;
  tesInterp.index.lat = 7;  // -37.5
  tesInterp.index.lon = 7;  // 63
  tesInterp.index.ls = 7;   // 210
  tesInterp.index.tyr = 0;

  // RUN
  t = tesInterp.getTideParameters(tesInterp.surfaceTemperatures, tesInterp.index.lat, tesInterp.index.lon, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(223.76, t.diurnalMean);
  EXPECT_DOUBLE_EQ(25.52, t.amplitude1D);
  EXPECT_DOUBLE_EQ(15.5, t.phase1D);
  EXPECT_DOUBLE_EQ(3.92, t.amplitude2D);
  EXPECT_DOUBLE_EQ(2.0, t.phase2D);

  // RUN
  t = tesInterp.getTideParameters(tesInterp.surfaceEWWinds, tesInterp.index.lat, tesInterp.index.lon, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(2.3, t.diurnalMean);
  EXPECT_DOUBLE_EQ(10.9, t.amplitude1D);
  EXPECT_DOUBLE_EQ(16.5, t.phase1D);
  EXPECT_DOUBLE_EQ(1.5, t.amplitude2D);
  EXPECT_DOUBLE_EQ(0.0, t.phase2D);

  // RUN
  t = tesInterp.getTideParameters(tesInterp.surfaceNSWinds, tesInterp.index.lat, tesInterp.index.lon, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-0.9, t.diurnalMean);
  EXPECT_DOUBLE_EQ(5.5, t.amplitude1D);
  EXPECT_DOUBLE_EQ(23.1, t.phase1D);
  EXPECT_DOUBLE_EQ(1.9, t.amplitude2D);
  EXPECT_DOUBLE_EQ(18.2, t.phase2D);

  // GIVEN (INPUTS)
  poleFactor = 1.0;
  tesInterp.index.hgt = 2;
  tesInterp.index.lat = 7;  // -37.5
  tesInterp.index.lon = 7;  // 63
  tesInterp.index.ls = 7;   // 210
  tesInterp.index.tyr = 1;

  // RUN
  t = tesInterp.getTideParameters(tesInterp.surfaceTemperatures, tesInterp.index.lat, tesInterp.index.lon, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(226.56, t.diurnalMean);
  EXPECT_DOUBLE_EQ(24.45, t.amplitude1D);
  EXPECT_DOUBLE_EQ(15.7, t.phase1D);
  EXPECT_DOUBLE_EQ(8.95, t.amplitude2D);
  EXPECT_DOUBLE_EQ(2.3, t.phase2D);

  // RUN
  t = tesInterp.getTideParameters(tesInterp.surfaceEWWinds, tesInterp.index.lat, tesInterp.index.lon, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(-3.2, t.diurnalMean);
  EXPECT_DOUBLE_EQ(10.1, t.amplitude1D);
  EXPECT_DOUBLE_EQ(17.3, t.phase1D);
  EXPECT_DOUBLE_EQ(10.0, t.amplitude2D);
  EXPECT_DOUBLE_EQ(6.6, t.phase2D);

  // RUN
  t = tesInterp.getTideParameters(tesInterp.surfaceNSWinds, tesInterp.index.lat, tesInterp.index.lon, poleFactor);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(1.0, t.diurnalMean);
  EXPECT_DOUBLE_EQ(6.6, t.amplitude1D);
  EXPECT_DOUBLE_EQ(1.6, t.phase1D);
  EXPECT_DOUBLE_EQ(9.4, t.amplitude2D);
  EXPECT_DOUBLE_EQ(22.6, t.phase2D);

  // TEAR-DOWN
}



TEST(TesSurfaceInterpolator, updateIndices)
{
  // SETUP
  TesSurfaceInterpolator surface;
  Position pos;

  // GIVEN (INPUTS)
  pos.height = 0.0;
  pos.longitude = 0.0_deg;
  pos.surfaceHeight = 3.45_km;
  surface.setPosition(pos);

  // RUN
  surface.updateIndices(0);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(0U, surface.getBaseIndex().hgt);
  EXPECT_EQ(0U, surface.getBaseIndex().lon);
  EXPECT_EQ(0U, surface.getBaseIndex().wlon);
  EXPECT_EQ(9U, surface.oneKilometerIndex);
  EXPECT_DOUBLE_EQ(0.0, surface.getDisplacements().lon);
  EXPECT_DOUBLE_EQ(0.5, surface.getDisplacements().wlon);

  // GIVEN (INPUTS)
  pos.height = 17.22;
  pos.surfaceHeight = 17.2_km;
  pos.longitude = 0.0_deg;
  surface.setPosition(pos);

  // RUN
  surface.updateIndices(0);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(1U, surface.getBaseIndex().hgt);
  EXPECT_EQ(0U, surface.getBaseIndex().lon);
  EXPECT_EQ(17U, surface.oneKilometerIndex);
  EXPECT_DOUBLE_EQ(0.0, surface.getDisplacements().lon);

  // GIVEN (INPUTS)
  pos.height = -2.0;
  pos.surfaceHeight = -12.4_km;
  pos.setLongitude(100.0, WEST_POSITIVE);
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

TEST(TesSurfaceInterpolator, update_0)
{
  // SETUP
  TesSurfaceInterpolator surface;
  Position pos;
  EphemerisState ephem;

  // GIVEN (INPUTS)
  pos.height = 0.0;
  pos.latitude = 0.0_deg;
  pos.setLongitude(0.0, WEST_POSITIVE);
  pos.surfaceHeight = 0.0;
  ephem.longitudeSun = 0.0_deg;
  ephem.solarTime = 0.0;
  surface.setPosition(pos);
  surface.setEphemerisState(ephem);
  surface.setMapYear(1);

  // RUN
  surface.updateIndices(0);
  surface.update();
  const AtmosphereState& atmos = surface.getAtmosphereState();
  const MarsAtmosphereState& atmosMars = atmos.getPlanetSpecificMetrics<MarsAtmosphereState>();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(198.45439004703826, atmos.temperature);
  EXPECT_DOUBLE_EQ(222.99, atmosMars.temperatureDaily);
  EXPECT_DOUBLE_EQ(185.46184263025577, atmosMars.temperatureMin);
  EXPECT_DOUBLE_EQ(281.60930301432143, atmosMars.temperatureMax);
  EXPECT_DOUBLE_EQ(640.27045036128891, atmos.pressure);
  EXPECT_DOUBLE_EQ(587.39999999999998, atmosMars.pressureDaily);
  EXPECT_DOUBLE_EQ(1.7016252244453706E-002, atmos.density);
  EXPECT_DOUBLE_EQ(1.3893437800588380E-002, atmosMars.densityDaily);
  EXPECT_DOUBLE_EQ(1.2457310016901004E-002, atmosMars.densityMin);
  EXPECT_DOUBLE_EQ(1.4447216799115967E-002, atmosMars.densityMax);
  EXPECT_DOUBLE_EQ(0.0, atmos.ewWind);
  EXPECT_DOUBLE_EQ(0.0, atmosMars.ewWindDaily);
  EXPECT_DOUBLE_EQ(0.0, atmos.nsWind);
  EXPECT_DOUBLE_EQ(0.0, atmosMars.nsWindDaily);
  EXPECT_DOUBLE_EQ(5.6741466627144783, atmos.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(11.505027640325896, surface.getPressureScaleHeightDaily());
  EXPECT_DOUBLE_EQ(6.4477614742218847, atmos.densityScaleHeight);
  EXPECT_DOUBLE_EQ(0.0, surface.temperatureAtFiveMeters);

  // TEAR-DOWN
}

TEST(TesSurfaceInterpolator, update_1)
{
  // SETUP
  TesSurfaceInterpolator surface;
  Position pos;
  EphemerisState ephem;

  // GIVEN (INPUTS)
  pos.height = 0.0;
  pos.latitude = 0.0_deg;
  pos.setLongitude(0.0, WEST_POSITIVE);
  pos.surfaceHeight = 0.0;
  ephem.longitudeSun = 0.0_deg;
  ephem.solarTime = 0.0;
  surface.setPosition(pos);
  surface.setEphemerisState(ephem);
  surface.setMapYear(1);

  // RUN
  surface.updateIndices(1);
  surface.update();
  const AtmosphereState& atmos = surface.getAtmosphereState();
  const MarsAtmosphereState& atmosMars = atmos.getPlanetSpecificMetrics<MarsAtmosphereState>();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(199.26106585382351, atmos.temperature);
  EXPECT_DOUBLE_EQ(217.55000000000001, atmosMars.temperatureDaily);
  EXPECT_DOUBLE_EQ(195.96813510924110, atmosMars.temperatureMin);
  EXPECT_DOUBLE_EQ(250.92466101354424, atmosMars.temperatureMax);
  EXPECT_DOUBLE_EQ(588.70802805089863, atmos.pressure);
  EXPECT_DOUBLE_EQ(587.14477576369291, atmosMars.pressureDaily);
  EXPECT_DOUBLE_EQ(1.5582555257457055E-002, atmos.density);
  EXPECT_DOUBLE_EQ(1.4234665949505157E-002, atmosMars.densityDaily);
  EXPECT_DOUBLE_EQ(1.2504312615025305E-002, atmosMars.densityMin);
  EXPECT_DOUBLE_EQ(0.014413457672824584, atmosMars.densityMax);
  EXPECT_DOUBLE_EQ(-2.3676411421720367, atmos.ewWind);
  EXPECT_DOUBLE_EQ(-2.25, atmosMars.ewWindDaily);
  EXPECT_DOUBLE_EQ(-3.8511743529662412E-002, atmos.nsWind);
  EXPECT_DOUBLE_EQ(0.89999999999999991, atmosMars.nsWindDaily);
  EXPECT_DOUBLE_EQ(10.782676991373533, atmos.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(11.505027640325896, surface.getPressureScaleHeightDaily());
  EXPECT_DOUBLE_EQ(12.252790318376629, atmos.densityScaleHeight);
  EXPECT_DOUBLE_EQ(199.26106585382351, surface.temperatureAtFiveMeters);

  // TEAR-DOWN
}

} // namespace
