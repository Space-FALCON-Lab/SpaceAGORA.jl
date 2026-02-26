//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Jupiter-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once
#ifdef _MSC_VER
#pragma warning( disable : 4250 )
#endif

#include "unittest_friend.h"
#include "JupiterCommon.h"
#include "PerturbedAtmosphere.h"
#include "Ephemeris.h"
#include "JupiterInputParameters.h"
#include "HeightModel.h"

//! \defgroup Cpp_Jupiter The C++ Interface for Jupiter 
//! \defgroup JupiterGRAM The Full List of JupiterGRAM Classes.

namespace GRAM {

//! \brief The primary interface for the Jupiter atmosphere model.
//!
//! The Jupiter atmosphere atmosphere model utilizes the simple height based data model.
//! \ingroup Cpp_Jupiter JupiterGRAM
class JupiterAtmosphere : public PerturbedAtmosphere, public JupiterCommon
{
//! \example Jupiter/examples/Jupiter.cpp
public:
  JupiterAtmosphere();
  JupiterAtmosphere(const JupiterAtmosphere& orig);
  virtual ~JupiterAtmosphere() override;

  void setInputParameters(const JupiterInputParameters& params);
  const InputParameters& getInputParameters() const { return inputParameters; }

  void update() override;

  static const std::string& getVersionString();

private:
  void updateReferenceValues() override;
  void updatePressureAtSurface() override;

  void getPerturbationFactors(greal& pertLow, greal& pertHigh) override;
  void getScaleParameters(greal& verticalScale, greal& horizontalScale) override;
  void getWindDeviations(greal& ewStdDev, greal& nsStdDev, greal& vertStdDev) override;

  HeightModel* jupiterModelPtr = NULL;     //!< The Jupiter height based data model.
  JupiterInputParameters inputParameters;  //!< User supplied parameters.

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(JupiterAtmosphere, getPerturbationFactors);
  FRIEND_TEST(JupiterAtmosphere, getScaleParameters);
  FRIEND_TEST(JupiterAtmosphere, getWindDeviations);
#endif // GRAM_UNIT_TEST
};

} // namespace
