#include "gtest/gtest.h"
#include "gram.h"
#include "VenusIRA.h"
#include "Position.h"
#include "AtmosphereState.h"

using namespace GRAM;


TEST(VenusIRA, update_low)
{
  // SETUP
  Position p;
  AtmosphereState s;
  VenusIRA venusIRA;

  // GIVEN (INPUTS)
  p.longitude = 56.0_deg;
  p.latitude = 67.3_deg;
  p.height = 89.3;
  p.elapsedTime = 0;
  venusIRA.setPosition(p);

  // RUN
  venusIRA.update();
  s = venusIRA.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(177.48493333333334, s.temperature);
  EXPECT_DOUBLE_EQ(45.998147537793351, s.pressure);
  EXPECT_DOUBLE_EQ(3.9607117120445077, s.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(0.001352954185219342, s.density);
  EXPECT_DOUBLE_EQ(4.1627345148682728, s.densityScaleHeight);
  EXPECT_DOUBLE_EQ(0.0000000000000000, s.helium.numberDensity);
  EXPECT_DOUBLE_EQ(0.0000000000000000, s.nitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(0.0000000000000000, s.hydrogen.numberDensity);
  EXPECT_DOUBLE_EQ(1.8074790794205518e+22, s.carbonDioxide.numberDensity);
  EXPECT_DOUBLE_EQ(6.9468734448383715e+20, s.dinitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(5.2340053753127552e+17, s.oxygen.numberDensity);
  EXPECT_DOUBLE_EQ(2.6170026876565709e+18, s.carbonMonoxide.numberDensity);
  EXPECT_DOUBLE_EQ(1.8772618541914542e+22, s.totalNumberDensity);
  EXPECT_DOUBLE_EQ(43.404866555994069, s.averageMolecularWeight);
  // TEAR-DOWN
}

TEST(VenusIRA, update_low_top)
{
  // SETUP
  Position p;
  AtmosphereState s;
  VenusIRA venusIRA;

  // GIVEN (INPUTS)
  p.longitude = 56.0_deg;
  p.latitude = 67.3_deg;
  p.height = 99.3;
  p.elapsedTime = 0;
  venusIRA.setPosition(p);

  // RUN
  venusIRA.update();
  s = venusIRA.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(168.56041666666667, s.temperature);
  EXPECT_DOUBLE_EQ(3.2951155109242856, s.pressure);
  EXPECT_DOUBLE_EQ(3.5349866912440029, s.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(1.0184057665520507E-004, s.density);
  EXPECT_DOUBLE_EQ(3.5210691536952012, s.densityScaleHeight);
  EXPECT_DOUBLE_EQ(1.3542642805830462e+21, s.carbonDioxide.numberDensity);
  EXPECT_DOUBLE_EQ(5.8667592008477917e+19, s.dinitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(2.1443499817072058e+17, s.oxygen.numberDensity);
  EXPECT_DOUBLE_EQ(1.0689681892903542e+18, s.carbonMonoxide.numberDensity);
  EXPECT_DOUBLE_EQ(4466835921509598.0, s.helium.numberDensity);
  EXPECT_DOUBLE_EQ(768295778499.65369, s.nitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(0.0000000000000000, s.hydrogen.numberDensity);
  EXPECT_DOUBLE_EQ(1.4142197433832026e+21, s.totalNumberDensity);
  EXPECT_DOUBLE_EQ(43.315166436209076, s.averageMolecularWeight);

  // TEAR-DOWN
}

TEST(VenusIRA, update_mid)
{
  // SETUP
  Position p;
  AtmosphereState s;
  VenusIRA venusIRA;

  // GIVEN (INPUTS)
  p.longitude = 56.0_deg;
  p.latitude = 67.3_deg;
  p.height = 127.3;
  p.elapsedTime = 0;
  venusIRA.setPosition(p);

  // RUN
  venusIRA.update();
  s = venusIRA.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(127.81600000000000, s.temperature);
  EXPECT_DOUBLE_EQ(0.0004615176779649479, s.pressure);
  EXPECT_DOUBLE_EQ(3.0676419144600753, s.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(1.699944692743425e-08, s.density);
  EXPECT_DOUBLE_EQ(2.9691865560964139, s.densityScaleHeight);
  EXPECT_DOUBLE_EQ(2.0430179348723581E+017, s.carbonDioxide.numberDensity);
  EXPECT_DOUBLE_EQ(15642168870768338., s.dinitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(28117952741338068., s.oxygen.numberDensity);
  EXPECT_DOUBLE_EQ(11392833553566016., s.carbonMonoxide.numberDensity);
  EXPECT_DOUBLE_EQ(55088980568581.617, s.helium.numberDensity);
  EXPECT_DOUBLE_EQ(154094284667402.28, s.nitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(0.0000000000000000, s.hydrogen.numberDensity);
  EXPECT_DOUBLE_EQ(2.5966393191814413E+017, s.totalNumberDensity);
  EXPECT_DOUBLE_EQ(39.144059086388999, s.averageMolecularWeight);

  // TEAR-DOWN
}

TEST(VenusIRA, update_mid_bottom)
{
  // SETUP
  Position p;
  AtmosphereState s;
  VenusIRA venusIRA;

  // GIVEN (INPUTS)
  p.longitude = 56.0_deg;
  p.latitude = 67.3_deg;
  p.height = 102.3;
  p.elapsedTime = 0;
  venusIRA.setPosition(p);

  // RUN
  venusIRA.update();
  s = venusIRA.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(162.88122000000001, s.temperature);
  EXPECT_DOUBLE_EQ(1.4145688501224025, s.pressure);
  EXPECT_DOUBLE_EQ(3.5515975801970061, s.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(4.5262970013581138E-005, s.density);
  EXPECT_DOUBLE_EQ(3.7621453435317895, s.densityScaleHeight);
  EXPECT_DOUBLE_EQ(6.025008616134852e+20, s.carbonDioxide.numberDensity);
  EXPECT_DOUBLE_EQ(2.4448191592643764e+19, s.dinitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(2.2357358679264016e+17, s.oxygen.numberDensity);
  EXPECT_DOUBLE_EQ(9.0300299837786291e+17, s.carbonMonoxide.numberDensity);
  EXPECT_DOUBLE_EQ(6731573142092238.0, s.helium.numberDensity);
  EXPECT_DOUBLE_EQ(4742667483283.3730, s.nitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(0.0000000000000000, s.hydrogen.numberDensity);
  EXPECT_DOUBLE_EQ(6.2808236610710903e+20, s.totalNumberDensity);
  EXPECT_DOUBLE_EQ(43.333538707001615, s.averageMolecularWeight);

  // TEAR-DOWN
}

TEST(VenusIRA, update_mid_top)
{
  // SETUP
  Position p;
  AtmosphereState s;
  VenusIRA venusIRA;
  EphemerisState ephem;

  // GIVEN (INPUTS)
  p.longitude = 56.0_deg;
  p.latitude = 67.3_deg;
  p.height = 147.3;
  p.elapsedTime = 0;
  venusIRA.setPosition(p);
  ephem.solarZenithAngle = 65.5_deg;
  venusIRA.setEphemerisState(ephem);

  // RUN
  venusIRA.update();
  s = venusIRA.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(154.00306896551737, s.temperature);
  EXPECT_DOUBLE_EQ(5.4909257809287227e-06, s.pressure);
  EXPECT_DOUBLE_EQ(0.0, s.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(1.0296406925814949e-10, s.density);
  EXPECT_DOUBLE_EQ(0.0, s.densityScaleHeight);
  EXPECT_DOUBLE_EQ(465033488545312.25, s.carbonDioxide.numberDensity);
  EXPECT_DOUBLE_EQ(189671729534985.62, s.dinitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(1539941319348112.0, s.oxygen.numberDensity);
  EXPECT_DOUBLE_EQ(238669199723162.75, s.carbonMonoxide.numberDensity);
  EXPECT_DOUBLE_EQ(15221360327962.426, s.helium.numberDensity);
  EXPECT_DOUBLE_EQ(11083577081264.332, s.nitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(119802452116.76515, s.hydrogen.numberDensity);
  EXPECT_DOUBLE_EQ(2459740477012917.0, s.totalNumberDensity);
  EXPECT_DOUBLE_EQ(24.010637171123939, s.averageMolecularWeight);

  // TEAR-DOWN
}

TEST(VenusIRA, update_high)
{
  // SETUP
  Position p;
  AtmosphereState s;
  VenusIRA venusIRA;
  EphemerisState ephem;

  // GIVEN (INPUTS)
  p.longitude = 56.0_deg;
  p.latitude = 67.3_deg;
  p.height = 211.3;
  p.elapsedTime = 0;
  venusIRA.setPosition(p);
  ephem.solarZenithAngle = 65.5_deg;
  venusIRA.setEphemerisState(ephem);

  // RUN
  venusIRA.update();
  s = venusIRA.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(281.85831034482754, s.temperature);
  EXPECT_DOUBLE_EQ(2.529359122361732e-07, s.pressure);
  EXPECT_DOUBLE_EQ(18.111842951857614, s.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(1.6750882615394788e-12, s.density);
  EXPECT_DOUBLE_EQ(16.611276324744708, s.densityScaleHeight);
  EXPECT_NEAR(155656405407.56903, s.carbonDioxide.numberDensity, 1.0e-3);
  EXPECT_DOUBLE_EQ(1007904972013.4441, s.dinitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(56531232564685.398, s.oxygen.numberDensity);
  EXPECT_DOUBLE_EQ(1378094735260.1990, s.carbonMonoxide.numberDensity);
  EXPECT_DOUBLE_EQ(1881245333175.7634, s.helium.numberDensity);
  EXPECT_DOUBLE_EQ(939145458948.02637, s.nitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(296650778017.00903, s.hydrogen.numberDensity);
  EXPECT_DOUBLE_EQ(62189930247507.391, s.totalNumberDensity);
  EXPECT_DOUBLE_EQ(15.520018296381624, s.averageMolecularWeight);

  // TEAR-DOWN
}

TEST(VenusIRA, update_high_bottom)
{
  // SETUP
  Position p;
  AtmosphereState s;
  VenusIRA venusIRA;
  EphemerisState ephem;

  // GIVEN (INPUTS)
  p.longitude = 56.0_deg;
  p.latitude = 67.3_deg;
  p.height = 153.3;
  p.elapsedTime = 0;
  venusIRA.setPosition(p);
  ephem.solarZenithAngle = 65.5_deg;
  venusIRA.setEphemerisState(ephem);

  // RUN
  venusIRA.update();
  s = venusIRA.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(230.99920689655190, s.temperature);
  EXPECT_DOUBLE_EQ(1.1706028120610132E-005, s.pressure);
  EXPECT_DOUBLE_EQ(0.0, s.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(1.6937723349835547E-010, s.density);
  EXPECT_DOUBLE_EQ(0.0, s.densityScaleHeight);
  EXPECT_NEAR(1092664654679021.9, s.carbonDioxide.numberDensity, 1.0e2);
  EXPECT_DOUBLE_EQ(300714741164634.06, s.dinitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(1629784273250475.0, s.oxygen.numberDensity);
  EXPECT_DOUBLE_EQ(404016942188736.12, s.carbonMonoxide.numberDensity);
  EXPECT_DOUBLE_EQ(5929662916792.1338, s.helium.numberDensity);
  EXPECT_DOUBLE_EQ(16342505205077.727, s.nitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(400728034124.14874, s.hydrogen.numberDensity);
  EXPECT_NEAR(3449853507438860.5, s.totalNumberDensity, 1.0e2);
  EXPECT_DOUBLE_EQ(27.790102326132399, s.averageMolecularWeight);

  // TEAR-DOWN
}

TEST(VenusIRA, update_thermos)
{
  // SETUP
  Position p;
  AtmosphereState s;
  VenusIRA venusIRA;
  EphemerisState ephem;

  // GIVEN (INPUTS)
  p.longitude = 56.0_deg;
  p.latitude = 67.3_deg;
  p.height = 353.3;
  p.elapsedTime = 0;
  venusIRA.setPosition(p);
  ephem.solarZenithAngle = 65.5_deg;
  venusIRA.setEphemerisState(ephem);

  // RUN
  venusIRA.update();
  s = venusIRA.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(282.37931034482756, s.temperature);
  EXPECT_DOUBLE_EQ(1.8172508407303739e-09, s.pressure);
  EXPECT_DOUBLE_EQ(86.473060831962798, s.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(2.6539671918346713e-15, s.density);
  EXPECT_DOUBLE_EQ(46.837249002393136, s.densityScaleHeight);
  EXPECT_DOUBLE_EQ(64.518599900320865, s.carbonDioxide.numberDensity);
  EXPECT_DOUBLE_EQ(1078328.8552546897, s.dinitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(21888568784.363403, s.oxygen.numberDensity);
  EXPECT_DOUBLE_EQ(1474028.2720919757, s.carbonMonoxide.numberDensity);
  EXPECT_DOUBLE_EQ(263120272664.33408, s.helium.numberDensity);
  EXPECT_DOUBLE_EQ(967357997.48212624, s.nitrogen.numberDensity);
  EXPECT_DOUBLE_EQ(180142305074.54443, s.hydrogen.numberDensity);
  EXPECT_DOUBLE_EQ(466121056942.37, s.totalNumberDensity);
  EXPECT_DOUBLE_EQ(3.4288440230950106, s.averageMolecularWeight);

  // TEAR-DOWN
}
TEST(VenusIRA, getReferenceValues)
{
  // SETUP
  VenusIRA venusIRA;
  greal refT, refP, refD;

  // GIVEN (INPUTS)
  greal height = 0.0;

  // RUN
  venusIRA.getReferenceValues(height, refT, refP, refD);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(735.3, refT);
  EXPECT_DOUBLE_EQ(64.79, refD);
  EXPECT_DOUBLE_EQ(9209495.778318001, refP);

  // GIVEN (INPUTS)
  height = 130.0;

  // RUN
  venusIRA.getReferenceValues(height, refT, refP, refD);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(166.8, refT);
  EXPECT_DOUBLE_EQ(1.8451151725569871e-08, refD);
  EXPECT_DOUBLE_EQ(0.00062731044945863932, refP);

  // GIVEN (INPUTS)
  height = 230.0;

  // RUN
  venusIRA.getReferenceValues(height, refT, refP, refD);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(229.9, refT);
  EXPECT_DOUBLE_EQ(1.43e-13, refD);
  EXPECT_DOUBLE_EQ(3.075e-8, refP);

  // GIVEN (INPUTS)
  height = 300.0;

  // RUN
  venusIRA.getReferenceValues(height, refT, refP, refD);

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(230.0, refT);
  EXPECT_DOUBLE_EQ(1.1109110838066177e-14, refD);
  EXPECT_DOUBLE_EQ(6.9782135088046824e-09, refP);

  // TEAR-DOWN
}
