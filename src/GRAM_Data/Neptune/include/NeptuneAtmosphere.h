//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Neptune-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once
#ifdef _MSC_VER
#pragma warning( disable : 4250 )
#endif

#include "gram.h"

#include "unittest_friend.h"
#include "NeptuneCommon.h"
#include "PerturbedAtmosphere.h"
#include "NeptuneMinMax.h"
#include "Ephemeris.h"
#include "NeptuneInputParameters.h"

//! \defgroup Cpp_Neptune The C++ Interface for Neptune
//! \brief Class references for the NeptuneGRAM interface.
//!
//! These are the classes needed for building an application that interfaces with NeptuneGRAM.
//! The interface for Neptune is declared in the following:
//! \code #include "NeptuneAtmosphere.h" \endcode 
//! An example using this class can be found \ref Neptune/examples/Neptune.cpp "here".

//! \defgroup NeptuneGRAM The Full List of NeptuneGRAM Classes.
//! \brief Class references for all Neptune models.
//!
//! These are the Neptune specific classes needed for building NeptuneGRAM.

namespace GRAM {

//! \brief The primary interface for the Neptune atmosphere model.
//!
//! The Neptune atmosphere atmosphere model utilizes the min-max data model.
//! \ingroup Cpp_Neptune NeptuneGRAM
class NeptuneAtmosphere : public PerturbedAtmosphere, public NeptuneCommon
{
//! \example Neptune/examples/Neptune.cpp
public:
  NeptuneAtmosphere();
  NeptuneAtmosphere(const NeptuneAtmosphere& orig);
  virtual ~NeptuneAtmosphere() override;

  void update() override;

  void setInputParameters(const NeptuneInputParameters& params);
  const InputParameters& getInputParameters() const { return inputParameters; }

  void setMinMaxFactor(greal factor, bool computeFlag);
  void setDinitrogenMoleFraction(greal n2mf);

  static const std::string& getVersionString();

private:
  void updateReferenceValues() override;
  void updatePressureAtSurface() override;
  void updateCompressibilityFactor() override;

  void getPerturbationFactors(greal& pertLow, greal& pertHigh) override;
  void getScaleParameters(greal& verticalScale, greal& horizontalScale) override;
  void getWindDeviations(greal& ewStdDev, greal& nsStdDev, greal& vertStdDev) override;

  NeptuneMinMax *neptuneMinMaxPtr = NULL; //!< The Neptune min-max data model.
  NeptuneInputParameters inputParameters; //!< User supplied parameters.

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(NeptuneAtmosphere, getPerturbationFactors);
  FRIEND_TEST(NeptuneAtmosphere, getScaleParameters);
  FRIEND_TEST(NeptuneAtmosphere, getWindDeviations);
#endif // GRAM_UNIT_TEST
};
} // namespace

