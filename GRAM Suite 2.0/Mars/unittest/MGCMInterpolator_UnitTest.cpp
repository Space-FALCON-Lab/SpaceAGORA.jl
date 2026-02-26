#include "unittest.h"
#include "MGCMInterpolator.h"

namespace GRAM {

class MGCMInterpolatorTest : public MGCMInterpolator
{
public:
  MGCMInterpolatorTest() = default;
  MGCMInterpolatorTest(const MGCMInterpolatorTest& orig) = default;
  virtual ~MGCMInterpolatorTest() = default;

  void update() override {}
  void updateIndices(size_t heightIndex) override {}
};

TEST(MGCMInterpolator, updateIndices_noOffset)
{
  // SETUP
  MGCMInterpolatorTest marsInterp;
  Position pos;
  EphemerisState ephem;

  // GIVEN (INPUTS)
  pos.latitude = -90.0_deg;
  ephem.longitudeSun = 0.0_deg;
  marsInterp.setPosition(pos);
  marsInterp.setEphemerisState(ephem);
  marsInterp.setDustOpticalDepth(0.3);

  // RUN  (low side)
  marsInterp.updateBaseIndices(25, 0.0_deg);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(0U, marsInterp.getBaseIndex().lat);
  EXPECT_EQ(1U, marsInterp.getBaseIndex().wlat);
  EXPECT_EQ(0U, marsInterp.getBaseIndex().ls);
  EXPECT_EQ(0U, marsInterp.getBaseIndex().od);
  EXPECT_DOUBLE_EQ(0.0, marsInterp.getDisplacements().lat);
  EXPECT_DOUBLE_EQ(0.0, marsInterp.getDisplacements().wlat);
  EXPECT_DOUBLE_EQ(0.0, marsInterp.getDisplacements().ls);
  EXPECT_DOUBLE_EQ(0.0, marsInterp.getDisplacements().od);

  // GIVEN (INPUTS)
  pos.latitude = 90.0_deg;
  ephem.longitudeSun = 360.0_deg;
  marsInterp.setPosition(pos);
  marsInterp.setEphemerisState(ephem);
  marsInterp.setDustOpticalDepth(3.0);

  // RUN  (high side)
  marsInterp.updateBaseIndices(25, 0.0_deg);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(23U, marsInterp.getBaseIndex().lat);
  EXPECT_EQ(23U, marsInterp.getBaseIndex().wlat);
  EXPECT_EQ(11U, marsInterp.getBaseIndex().ls);
  EXPECT_EQ(1U, marsInterp.getBaseIndex().od);
  EXPECT_DOUBLE_EQ(1.0, marsInterp.getDisplacements().lat);
  EXPECT_DOUBLE_EQ(1.0, marsInterp.getDisplacements().wlat);
  EXPECT_DOUBLE_EQ(1.0, marsInterp.getDisplacements().ls);
  EXPECT_DOUBLE_EQ(1.0, marsInterp.getDisplacements().od);

  // GIVEN (INPUTS)
  pos.latitude = 40.0_deg;
  ephem.longitudeSun = 36.0_deg;
  marsInterp.setPosition(pos);
  marsInterp.setEphemerisState(ephem);
  marsInterp.setDustOpticalDepth(2.3);

  // RUN (random)
  marsInterp.updateBaseIndices(25, 0.0_deg);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(17U, marsInterp.getBaseIndex().lat);
  EXPECT_EQ(17U, marsInterp.getBaseIndex().wlat);
  EXPECT_EQ(1U, marsInterp.getBaseIndex().ls);
  EXPECT_EQ(1U, marsInterp.getBaseIndex().od);
  EXPECT_DOUBLE_EQ(2.5 / 7.5, marsInterp.getDisplacements().lat);
  EXPECT_DOUBLE_EQ(6.25 / 7.5, marsInterp.getDisplacements().wlat);
  EXPECT_DOUBLE_EQ(6.0 / 30.0, marsInterp.getDisplacements().ls);
  EXPECT_DOUBLE_EQ(0.75814655591088631, marsInterp.getDisplacements().od);

  // GIVEN (INPUTS)
  pos.latitude = -55.0_deg;
  ephem.longitudeSun = 185.3_deg;
  marsInterp.setPosition(pos);
  marsInterp.setEphemerisState(ephem);
  marsInterp.setDustOpticalDepth(0.4);

  // RUN (random)
  marsInterp.updateBaseIndices(25, 0.0);

  // EXPECT (OUTPUTS)
  EXPECT_EQ(4U, marsInterp.getBaseIndex().lat);
  EXPECT_EQ(5U, marsInterp.getBaseIndex().wlat);
  EXPECT_EQ(6U, marsInterp.getBaseIndex().ls);
  EXPECT_EQ(0U, marsInterp.getBaseIndex().od);
  EXPECT_DOUBLE_EQ(5.0 / 7.5, marsInterp.getDisplacements().lat);
  EXPECT_DOUBLE_EQ(1.25 / 7.5, marsInterp.getDisplacements().wlat);
  EXPECT_NEAR(5.3 / 30.0, marsInterp.getDisplacements().ls, 1.0e-12);
  EXPECT_DOUBLE_EQ(0.23894399559369162, marsInterp.getDisplacements().od);

  // TEAR-DOWN
}

} // namespace
