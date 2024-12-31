#include "gtest/gtest.h"
#define GRAM_UNIT_TEST
#include "StewartModel.h"

namespace GRAM {

TEST(StewartModel, updateExosphericTemperature)
{
  // SETUP
  StewartModel stewartModel;
  MarsInputParameters params;

  // GIVEN (INPUTS)
  params.F107 = 69.9;
  params.exosphericTemperatureFactor = 1.3;
  stewartModel.setInputParameters(params);

  Position position;
  position.latitude = 66.6_deg;
  stewartModel.setPosition(position);
  stewartModel.setExosphericTemperatureOffset(11.1);

  EphemerisState ephem;
  ephem.subsolarLatitude = 22.2_deg;
  ephem.solarTime = 12.0;
  ephem.orbitalRadius = 1.5;
  stewartModel.setEphemerisState(ephem);

  // RUN
  stewartModel.updateExosphericTemperature();
  greal exoTemp = stewartModel.getExosphericTemperature();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(303.13475240523059, exoTemp);

  // TEAR-DOWN
}

TEST(StewartModel, updateThermos)
{
  // SETUP
  StewartModel stewartModel;
  MarsInputParameters params;

  // GIVEN (INPUTS)
  Position position;
  position.height = 555.5;
  position.latitude = 66.6_deg;
  position.latitudeRadius = 3380.9491876960001;
  stewartModel.setPosition(position);

  EphemerisState ephem;
  ephem.solarTime = 12.0;
  ephem.orbitalRadius = 1.5;
  stewartModel.setEphemerisState(ephem);

  params.F107 = 69.9;
  params.exosphericTemperatureFactor = 1.3;
  stewartModel.setInputParameters(params);

  stewartModel.setThermosphereBase(300.0, 300.0);
  stewartModel.exosphericTemperature = 303.13475240523059;

  // RUN
  stewartModel.updateThermos();
  AtmosphereState atmos = stewartModel.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(2.2668969813724971e-08, atmos.pressure);
  EXPECT_DOUBLE_EQ(303.13475239570329, atmos.temperature);
  EXPECT_DOUBLE_EQ(1.3629138665672936e-13, atmos.density);
  EXPECT_DOUBLE_EQ(60.448711299413283, atmos.pressureScaleHeight);
  EXPECT_DOUBLE_EQ(52.021028422111364, atmos.densityScaleHeight);
  EXPECT_DOUBLE_EQ(5416425.4903865056, atmos.totalNumberDensity);
  EXPECT_DOUBLE_EQ(15.153276207026357, atmos.averageMolecularWeight);

  // TEAR-DOWN
}

}