#include "unittest.h"
#define GRAM_UNIT_TEST
#include "gram.h"
#include "TesInterpolator.h"

namespace GRAM {

class TesInterpolatorTest : public TesInterpolator
{
public:
  TesInterpolatorTest() = default;
  TesInterpolatorTest(const TesInterpolatorTest& orig) = default;
  virtual ~TesInterpolatorTest() = default;

  void update() override {}
  void updateIndices(size_t heightIndex) override {}
};


TEST(TesInterpolator, updateIndices_noOffset)
{
  // SETUP
  TesInterpolatorTest tesInterp;
  Position pos;
  EphemerisState ephem;

  // GIVEN (INPUTS)
  pos.latitude = -90.0_deg;
  ephem.longitudeSun = 0.0_deg;
  tesInterp.setPosition(pos);
  tesInterp.setEphemerisState(ephem);
  tesInterp.setMapYear(1);

  // RUN  (low side)
  tesInterp.updateBaseIndices(25, 0.0_deg);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(0U, tesInterp.getBaseIndex().lat);
  EXPECT_EQ(0U, tesInterp.getBaseIndex().ls);
  EXPECT_EQ(0U, tesInterp.getBaseIndex().tyr);
  EXPECT_DOUBLE_EQ(0.0, tesInterp.getDisplacements().lat);
  EXPECT_DOUBLE_EQ(0.0, tesInterp.getDisplacements().ls);

  // GIVEN (INPUTS)
  pos.latitude = 90.0_deg;
  ephem.longitudeSun = 360.0_deg;
  tesInterp.setPosition(pos);
  tesInterp.setEphemerisState(ephem);
  tesInterp.setMapYear(2);

  // RUN  (high side)
  tesInterp.updateBaseIndices(25, 0.0_deg);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(23U, tesInterp.getBaseIndex().lat);
  EXPECT_EQ(11U, tesInterp.getBaseIndex().ls);
  EXPECT_EQ(1U, tesInterp.getBaseIndex().tyr);
  EXPECT_DOUBLE_EQ(1.0, tesInterp.getDisplacements().lat);
  EXPECT_DOUBLE_EQ(1.0, tesInterp.getDisplacements().ls);

  // GIVEN (INPUTS)
  pos.latitude = 40.0_deg;
  ephem.longitudeSun = 36.0_deg;
  tesInterp.setPosition(pos);
  tesInterp.setEphemerisState(ephem);
  tesInterp.setMapYear(1);

  // RUN (random)
  tesInterp.updateBaseIndices(25, 0.0_deg);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(17U, tesInterp.getBaseIndex().lat);
  EXPECT_EQ(1U, tesInterp.getBaseIndex().ls);
  EXPECT_EQ(0U, tesInterp.getBaseIndex().tyr);
  EXPECT_DOUBLE_EQ(2.5 / 7.5, tesInterp.getDisplacements().lat);
  EXPECT_DOUBLE_EQ(6.0 / 30.0, tesInterp.getDisplacements().ls);

  // GIVEN (INPUTS)
  pos.latitude = -55.0_deg;
  ephem.longitudeSun = 185.3_deg;
  tesInterp.setPosition(pos);
  tesInterp.setEphemerisState(ephem);
  tesInterp.setMapYear(1);

  // RUN (random)
  tesInterp.updateBaseIndices(25, 0.0);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(4U, tesInterp.getBaseIndex().lat);
  EXPECT_EQ(6U, tesInterp.getBaseIndex().ls);
  EXPECT_EQ(0U, tesInterp.getBaseIndex().tyr);
  EXPECT_DOUBLE_EQ(5.0 / 7.5, tesInterp.getDisplacements().lat);
  EXPECT_NEAR(5.3 / 30.0, tesInterp.getDisplacements().ls, 1.0e-12);

  // TEAR-DOWN
}

} // namespace
