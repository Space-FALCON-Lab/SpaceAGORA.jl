#include "unittest.h"
#define GRAM_UNIT_TEST
#include "gram.h"
#include "MarsInterpolatorBase.h"

namespace GRAM {

class MarsInterpolatorTest : public MarsInterpolatorBase
{
public:
  MarsInterpolatorTest() = default;
  MarsInterpolatorTest(const MarsInterpolatorTest& orig) = default;
  virtual ~MarsInterpolatorTest() = default;

  void update() override {}
  void updateIndices(size_t heightIndex) override {}
};

TEST(MarsInterpolatorBase, getTideValue_UVT)
{
  // SETUP
  MarsInterpolatorTest mtm;
  MarsTideParameters tp;
  greal ltst;

  // GIVEN (INPUTS)
  tp.diurnalMean = 3.3;
  tp.amplitude1D = 1.2;
  tp.phase1D = 4.5;
  tp.amplitude2D = 6.7;
  tp.phase2D = 4.3;
  ltst = 23.4;

  // RUN
  greal tideValue = mtm.getTideValue(tp, ltst, MarsInterpolatorBase::UVT);

  // EXPECT (OUTPUTS)
  EXPECT_NEAR(-2.03895836860726, tideValue, 1.0e-13);

  // GIVEN (INPUTS)
  tp.diurnalMean = -5.3;
  tp.amplitude1D = -2.2;
  tp.phase1D = -4.9;
  tp.amplitude2D = -23.7;
  tp.phase2D = -44.3;
  ltst = 2.4;

  // RUN
  tideValue = mtm.getTideValue(tp, ltst, MarsInterpolatorBase::UVT);

  // EXPECT (OUTPUTS)
  EXPECT_NEAR(-22.9839841962159, tideValue, 1.0e-13);

  // TEAR-DOWN
}

TEST(MarsInterpolatorBase, getTideValue_DP)
{
  // SETUP
  MarsInterpolatorTest mtm;
  MarsTideParameters tp;
  greal ltst;

  // GIVEN (INPUTS)
  tp.diurnalMean = 3.3;
  tp.amplitude1D = 1.2;
  tp.phase1D = 4.5;
  tp.amplitude2D = 6.7;
  tp.phase2D = 4.3;
  ltst = 23.4;

  // RUN
  greal tideValue = mtm.getTideValue(tp, ltst, MarsInterpolatorBase::DP);

  // EXPECT (OUTPUTS)
  EXPECT_NEAR(3.12381437383596, tideValue, 1.0e-13);

  // GIVEN (INPUTS)
  tp.diurnalMean = 5.3;
  tp.amplitude1D = -2.2;
  tp.phase1D = -4.9;
  tp.amplitude2D = -23.7;
  tp.phase2D = -44.3;
  ltst = 2.4;

  // RUN
  tideValue = mtm.getTideValue(tp, ltst, MarsInterpolatorBase::DP);

  // EXPECT (OUTPUTS)
  EXPECT_NEAR(4.36274883760056, tideValue, 1.0e-13);

  // GIVEN (INPUTS)
  tp.diurnalMean = -5.3;
  tp.amplitude1D = -2.2;
  tp.phase1D = -4.9;
  tp.amplitude2D = -23.7;
  tp.phase2D = -44.3;
  ltst = 2.4;

  // RUN
  tideValue = mtm.getTideValue(tp, ltst, MarsInterpolatorBase::DP);

  // EXPECT (OUTPUTS)
  EXPECT_NEAR(-0.53, tideValue, 1.0e-13);

  // TEAR-DOWN
}

TEST(MarsInterpolatorBase, getTideValue_BIG_DP)
{
  // SETUP
  MarsInterpolatorTest mtm;
  MarsTideParameters tp;
  greal ltst;

  // GIVEN (INPUTS)
  tp.diurnalMean = 3.3;
  tp.amplitude1D = 1.2;
  tp.phase1D = 4.5;
  tp.amplitude2D = 6.7;
  tp.phase2D = 4.3;
  ltst = 23.4;

  // RUN
  greal tideValue = mtm.getTideValue(tp, ltst, MarsInterpolatorBase::BIG_DP);

  // EXPECT (OUTPUTS)
  EXPECT_NEAR(3.13402630163318, tideValue, 1.0e-13);

  // GIVEN (INPUTS)
  tp.diurnalMean = 5.3;
  tp.amplitude1D = -2.2;
  tp.phase1D = -4.9;
  tp.amplitude2D = -23.7;
  tp.phase2D = -44.3;
  ltst = 2.4;

  // RUN
  tideValue = mtm.getTideValue(tp, ltst, MarsInterpolatorBase::BIG_DP);

  // EXPECT (OUTPUTS)
  EXPECT_NEAR(4.32718283058706, tideValue, 1.0e-13);

  // GIVEN (INPUTS)
  tp.diurnalMean = -5.3;
  tp.amplitude1D = -2.2;
  tp.phase1D = -4.9;
  tp.amplitude2D = -23.7;
  tp.phase2D = -44.3;
  ltst = 2.4;

  // RUN
  tideValue = mtm.getTideValue(tp, ltst, MarsInterpolatorBase::BIG_DP);

  // EXPECT (OUTPUTS)
  EXPECT_NEAR(-4.32718283058706, tideValue, 1.0e-13);

  // TEAR-DOWN
}

TEST(MarsInterpolatorBase, getSolMinMax)
{
  // SETUP
  MarsInterpolatorTest mtm;
  MarsTideParameters tp;
  greal minValue, maxValue;

  // GIVEN (INPUTS)
  tp.diurnalMean = 3.3;
  tp.amplitude1D = 1.2;
  tp.phase1D = 4.5;
  tp.amplitude2D = 6.7;
  tp.phase2D = 4.3;

  // RUN
  mtm.getSolMinMax(tp, MarsInterpolatorBase::UVT, minValue, maxValue);

  // EXPECT (OUTPUTS)
  EXPECT_NEAR(-3.4741433126514845, minValue, 1.0e-13);
  EXPECT_NEAR(11.107245715635994, maxValue, 1.0e-13);

  // RUN
  mtm.getSolMinMax(tp, MarsInterpolatorBase::DP, minValue, maxValue);

  // EXPECT (OUTPUTS)
  EXPECT_NEAR(3.0764532706825007, minValue, 1.0e-13);
  EXPECT_NEAR(3.5576391086159878, maxValue, 1.0e-13);

  // RUN
  mtm.getSolMinMax(tp, MarsInterpolatorBase::BIG_DP, minValue, maxValue);

  // EXPECT (OUTPUTS)
  EXPECT_NEAR(3.0904383086271894, minValue, 1.0e-13);
  EXPECT_NEAR(3.5601459385437999, maxValue, 1.0e-13);

  // TEAR-DOWN
}

} // namespace
