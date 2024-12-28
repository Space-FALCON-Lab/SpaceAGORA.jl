//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Uranus-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once
#ifdef _MSC_VER
#pragma warning( disable : 4250 )
#endif

#include "unittest_friend.h"
#include "UranusCommon.h"
#include "PerturbedAtmosphere.h"
#include "HeightModel.h"
#include "Ephemeris.h"
#include "UranusInputParameters.h"

//! \defgroup Cpp_Uranus The C++ Interface for Uranus 
//! \defgroup UranusGRAM The Full List of UranusGRAM Classes.

namespace GRAM {

//! \brief The primary interface for the Uranus atmosphere model.
//!
//! The Uranus atmosphere atmosphere model utilizes the simple height based data model.
//! \ingroup Cpp_Uranus UranusGRAM
class UranusAtmosphere : public PerturbedAtmosphere, public UranusCommon
{
//! \example Uranus/examples/Uranus.cpp
public:
  UranusAtmosphere();
  UranusAtmosphere(const UranusAtmosphere& orig);
  virtual ~UranusAtmosphere() override;

  void setInputParameters(const UranusInputParameters& params);
  const InputParameters& getInputParameters() const { return inputParameters; }

  void update() override;

  static const std::string& getVersionString();

private:
  void updateReferenceValues() override;
  void updatePressureAtSurface() override;

  void getPerturbationFactors(greal& pertLow, greal& pertHigh) override;
  void getScaleParameters(greal& verticalScale, greal& horizontalScale) override;
  void getWindDeviations(greal& ewStdDev, greal& nsStdDev, greal& vertStdDev) override;

  HeightModel* uranusModelPtr = NULL;     //!< The Uranus height based data model.
  UranusInputParameters inputParameters;  //!< User supplied parameters.

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(UranusAtmosphere, getPerturbationFactors);
  FRIEND_TEST(UranusAtmosphere, getScaleParameters);
  FRIEND_TEST(UranusAtmosphere, getWindDeviations);
#endif // GRAM_UNIT_TEST
};

} // namespace
