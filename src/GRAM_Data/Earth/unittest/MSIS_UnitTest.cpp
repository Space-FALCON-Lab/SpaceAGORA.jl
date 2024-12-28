#include "unittest.h"
#include "MSIS.h"
 
namespace GRAM {

TEST(MSIS, update)
{
  // SETUP
  MSIS msis;
  EarthInputParameters params;
  Position pos;

  // GIVEN (INPUTS)
  params.minute = 0;
  params.seconds = 0.0;
  params.hour = 0;
  params.month = 1;
  params.day = 1;
  params.year = 2019;
  params.dailyF10 = 230.0;
  params.meanF10 = 230.4;
  params.ap = 16.7;
  msis.setInputParameters(params);
  pos.height = 95.6;
  msis.setPosition(pos);
  msis.updatePosition();

  // RUN
  msis.update();
  const AtmosphereState& atmos = msis.getAtmosphereState();

  // EXPECT (OUTPUTS)
  EXPECT_DOUBLE_EQ(187.45855890918637, atmos.temperature);
  EXPECT_DOUBLE_EQ(0.080540464712296006, atmos.pressure);
  EXPECT_DOUBLE_EQ(1.4750563893845838e-06, atmos.density);
  EXPECT_DOUBLE_EQ(6.1116923880153751, atmos.ewWind);
  EXPECT_DOUBLE_EQ(14.571726052687104, atmos.nsWind);
  EXPECT_DOUBLE_EQ(0.010521386312051269, atmos.verticalWind);
  EXPECT_DOUBLE_EQ(28.545256584701857, atmos.averageMolecularWeight);
  EXPECT_DOUBLE_EQ(2.5189009841906401e+19, atmos.dinitrogen.numberDensity);

  // TEAR-DOWN
}

}