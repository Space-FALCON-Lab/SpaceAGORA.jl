//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Titan-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once
#ifdef _MSC_VER
#pragma warning( disable : 4250 )
#endif

#include "unittest_friend.h"
#include "TitanCommon.h"
#include "PerturbedAtmosphere.h"
#include "Yelle.h"
#include "TitanGCM.h"
#include "Ephemeris.h"
#include "TitanInputParameters.h"

//! \defgroup Cpp_Titan The C++ Interface for Titan 
//! \brief Class references for the TitanGRAM interface.
//!
//! These are the classes needed for building an application that interfaces with TitanGRAM.
//! The interface for Titan is declared in the following:
//! \code #include "TitanAtmosphere.h" \endcode 
//! An example using this class can be found \ref Titan/examples/Titan.cpp "here".

//! \defgroup TitanGRAM The Full List of TitanGRAM Classes.
//! \brief Class references for all Titan models.
//!
//! These are the Titan specific classes needed for building TitanGRAM.

namespace GRAM {

//! \brief The main Titan Atmosphere class.
//!
//! This is the primary interface into the Titan Atmosphere models.  This class
//! allows the user to choose the appropriate model, Yelle or GCM, to run.  The 
//! Yelle model is based on the min-max data model.  The GCM model is a layered
//! model with a surface to 0.3 mbar layer, a middle layer up to 925 km, and an
//! exoshere layer above 925 km.
//! \ingroup Cpp_Titan TitanGRAM
class TitanAtmosphere : public PerturbedAtmosphere, public TitanCommon
{
//! \example Titan/examples/Titan.cpp
public:
  TitanAtmosphere();
  TitanAtmosphere(const TitanAtmosphere& orig);
  virtual ~TitanAtmosphere() override;

  void setInputParameters(const TitanInputParameters& params);
  const InputParameters& getInputParameters() const { return inputParameters; }

  void setModelType(TitanModelType type);
  void setMinMaxFactor(greal minMaxFactor, bool computeFlag);
  void setMethaneMoleFraction(greal mmf);

  void update() override;

  static const std::string& getVersionString();

private:
  void updateReferenceValues() override;
  void updatePressureAtSurface() override;

  void getPerturbationFactors(greal& pertLow, greal& pertHigh) override;
  void getScaleParameters(greal& verticalScale, greal& horizontalScale) override;
  void getWindDeviations(greal& ewStdDev, greal& nsStdDev, greal& vertStdDev) override;

  TitanModelType modelType = Yelle97;   //!< Select between Yelle97 model and GCM95.
  Atmosphere* atmosModelPtr = NULL;
  Yelle* yellePtr = NULL;               //!< The Yelle atmosphere model.
  TitanGCM* titanGCMPtr = NULL;         //!< The GCM atmosphere model.
  TitanInputParameters inputParameters; //!< User supplied parameters.
  greal userMethaneMoleFraction;        //!< A user supplied override value.
  bool userSuppliedMethaneMoleFraction = false;  //!< Set to true if the user overrides the methane mole fraction.

#ifdef GRAM_UNIT_TEST
  FRIEND_TEST(TitanAtmosphere, getPerturbationFactors);
  FRIEND_TEST(TitanAtmosphere, getScaleParameters);
  FRIEND_TEST(TitanAtmosphere, getWindDeviations);
#endif // GRAM_UNIT_TEST

};

} // namespace
