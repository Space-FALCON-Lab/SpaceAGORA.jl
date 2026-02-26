//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once
#ifdef _MSC_VER
#pragma warning( disable : 4250 )
#endif

#include "VenusCommon.h"
#include "PerturbedAtmosphere.h"
#include "VenusIRA.h"
#include "Ephemeris.h"
#include "VenusInputParameters.h"
#include "VenusTopography.h"

//! \defgroup Cpp_Venus The C++ Interface for Venus
//! \brief Class references for the VenusGRAM interface.
//!
//! These are the classes needed for building an application that interfaces with VenusGRAM.
//! The interface for Neptune is declared in the following:
//! \code #include "VenusAtmosphere.h" \endcode 
//! An example using this class can be found \ref Venus/examples/Venus.cpp "here".

//! \defgroup VenusGRAM The Full List of VenusGRAM Classes.
//! \brief Class references for all Venus models.
//!
//! These are the Venus specific classes needed for building VenusGRAM.

namespace GRAM {

//! \brief The main Venus Atmosphere class.
//!
//! The main Venus Atmosphere which chooses the appropriate model to run.
//! \ingroup Cpp_Venus VenusGRAM
class VenusAtmosphere : public PerturbedAtmosphere, public VenusCommon
{
//! \example Venus/examples/Venus.cpp
public:
  VenusAtmosphere();
  VenusAtmosphere(const VenusAtmosphere& orig);
  virtual ~VenusAtmosphere() override;

  void setInputParameters(const VenusInputParameters& params);
  const InputParameters& getInputParameters() const { return inputParameters; }

  void update() override;

  static const std::string& getVersionString();

private:
  void updateReferenceValues() override;
  void updatePressureAtSurface() override;
  void updateCompressibilityFactor() override;

  void getPerturbationFactors(greal& pertLow, greal& pertHigh) override;
  void getScaleParameters(greal& verticalScale, greal& horizontalScale) override;
  void getWindDeviations(greal& ewStdDev, greal& nsStdDev, greal& vertStdDev) override;

  VenusIRA* venusIRAPtr = NULL;          //!< The Venus IRA data model.
  VenusInputParameters inputParameters;  //!< User supplied parameters.
  VenusTopography topography;

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(VenusAtmosphere, getPerturbationFactors);
  FRIEND_TEST(VenusAtmosphere, getScaleParameters);
  FRIEND_TEST(VenusAtmosphere, getWindDeviations);
#endif // GRAM_UNIT_TEST
};

} // namespace
