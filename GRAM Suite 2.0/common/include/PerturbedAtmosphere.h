//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"
#include "PerturbationState.h"
#include "Atmosphere.h"
#include "Ephemeris.h"
//#include "WHRandomNumberGenerator.h"  // The legacy version.
#include "RandomNumberGenerator.h"    // The new version.
#include "AuxiliaryAdapter.h"

namespace GRAM {

//! \brief Produces random perturbations on an Atmosphere.
//!
//! An atmosphere model that inherits PerturbedAtmosphere is capable of producing
//! randomized perturbations on atmospheric values.
//! \ingroup CommonGRAM
class PerturbedAtmosphere : public Atmosphere, public AuxiliaryAdapter
{
//! \example common/examples/UnitTestRunner.cpp
//! \example common/examples/AtmosphereExample.cpp
//! \example common/examples/TrajectoryExample.cpp
//! \example common/examples/MonteCarloExample.cpp
//! \example common/examples/NamelistExample.cpp
//! \example common/examples/EphemerisExample.cpp
//! \example common/examples/PerturbationExample.cpp
public:
  PerturbedAtmosphere();
  PerturbedAtmosphere(const PerturbedAtmosphere& orig);
  virtual ~PerturbedAtmosphere() override = default;

  virtual void setPosition(const Position& pos) override;
  virtual void setDelta(const Position& delta);

  virtual void setStartTime(const GramTime& gramTime) { time = gramTime; }
  const GramTime& getStartTime() const { return time; }

  virtual void setSeed(int seed);
  int getSeed() const { return randomNumberGenerator.getSeed(); }

  virtual void setInputParameters(const InputParameters& params);
  virtual const InputParameters& getInputParameters() const = 0;

  virtual void setMinRelativeStepSize(greal stepSize);
  greal getMinRelativeStepSize() const { return minRelativeStepSize; }

  virtual void setPerturbationScales(greal scaleFactor);
  greal getPerturbationScaleFactor() const { return densityPerturbationScale; }

  virtual void setDensityPerturbationScale(greal scaleFactor);
  greal getDensityPerturbationScale() const { return densityPerturbationScale; }

  virtual void setEWWindPerturbationScale(greal scaleFactor);
  greal getEWWindPerturbationScale() const { return ewWindPerturbationScale; }

  virtual void setNSWindPerturbationScale(greal scaleFactor);
  greal getNSWindPerturbationScale() const { return nsWindPerturbationScale; }

  virtual void setVerticalWindPerturbationScale(greal scaleFactor);
  greal getVerticalWindPerturbationScale() const { return verticalWindPerturbationScale; }

  virtual void setPerturbationState(PerturbationState& pState);
  const PerturbationState& getPerturbationState() const  { return pertState; }

  virtual void setPerturbationAction(PerturbationAction action) { perturbationAction = action; }

  void setSpiceDataPath(const std::string& path);
  void setEphemerisState(const EphemerisState& state) override { ephem = state; userSuppliedEphemeris = true; }
  void setEphemerisFastModeOn(bool flag) { ephemeris.setFastModeOn(flag); }
  void setSubsolarUpdateTime(greal utime) { ephemeris.setSubsolarUpdateTime(utime); }

  GRAM_BODY getBody() const { return gramBody; }
  bool hasVerticalWindsInModel() const { return hasVerticalWinds; }

  void findDates(greal targetLongitudeSun, greal targetSolarTime, GramTime gramTime[3], greal lonSun[3], greal tlst[3]);

protected:
  virtual void setGramBody(GRAM_BODY body) { gramBody = body; }
  virtual void setVerticalWinds(bool hasVertical) { hasVerticalWinds = hasVertical; }
  virtual void updateEphemeris();
  virtual void updatePerturbations(greal meanDensity = -9999.0, greal meanEWWind = -9999.0, greal meanNSWind = -9999.0);
  virtual void getPerturbationFactors(greal& pertLow, greal& pertHigh) = 0;
  virtual void getScaleParameters(greal& verticalScale, greal& horizontalScale) = 0;
  virtual void getWindDeviations(greal& ewStdDev, greal& nsStdDev, greal& vertStdDev) = 0;
  greal normalPercentagePoint(greal P);

  GRAM_BODY gramBody = NO_BODY;  //!< Identification of the planet or moon.
  bool userSuppliedEphemeris = false;      //!< Set to true if the user supplied the ephemeris data.
  Ephemeris ephemeris;           //!< The ephemeris computational engine.
  GramTime time;                 //!< Set start time prior to first call to update.
  //WHRandomNumberGenerator randomNumberGenerator;
  RandomNumberGenerator randomNumberGenerator;
  PerturbationAction perturbationAction = UPDATE_PERTS; //!< Designates update action with perturbations.
  greal densityPerturbationScale = 1.0;                 //!< Scales density perturbations.
  greal ewWindPerturbationScale = 1.0;                  //!< Scales East/West wind perturbations.
  greal nsWindPerturbationScale = 1.0;                  //!< Scales North/South wind perturbations.
  greal verticalWindPerturbationScale = 1.0;            //!< Scales vertical wind perturbations.
  greal minRelativeStepSize = 0.0;                      //!< Minimum relative step size for perturbations (0.0 - 1.0) (was corlmin)
  
  // Outputs
  greal& relativeStepSize = atmos.relativeStepSize;     //!< Ratio of step size to minimum step size for assured accuracy in perturbations (should be >= 1)  (was corlim)

  Position planetodeticPosition;
  Position savedPosition;
  Position previousPosition;                            //!< The position prior to the current position.
  Position delta;                                       //!< The differnce in the prior and current positions.

  PerturbationState pertState;                          //!< The perturbation state.
  greal& densityRandomNumber = pertState.densityRandomNumber;             //!< \copydoc PerturbationState::densityRandomNumber
  greal& ewWindRandomNumber = pertState.ewWindRandomNumber;               //!< \copydoc PerturbationState::ewWindRandomNumber
  greal& nsWindRandomNumber = pertState.nsWindRandomNumber;               //!< \copydoc PerturbationState::nsWindRandomNumber
  greal& verticalWindRandomNumber = pertState.verticalWindRandomNumber;   //!< \copydoc PerturbationState::verticalWindRandomNumber
  greal& densityRho = pertState.densityRho;                               //!< \copydoc PerturbationState::densityRho
  greal& ewWindRho = pertState.ewWindRho;                                 //!< \copydoc PerturbationState::ewWindRho
  greal& nsWindRho = pertState.nsWindRho;                                 //!< \copydoc PerturbationState::nsWindRho
  greal& verticalWindRho = pertState.verticalWindRho;                     //!< \copydoc PerturbationState::verticalWindRho

  greal perturbationStep = 0.0;                         //!< The perturbation step size.
  greal RHOd = 0.0;                                     //!< The density rho.
  greal RHOu = 0.0;                                     //!< The EW rho.
  greal RHOv = 0.0;                                     //!< The NS rho.
  greal RHOw = 0.0;                                     //!< The vertical rho.
  bool userPertState = false;                           //!< If true, user overrides the perturbation state.
  bool hasVerticalWinds = false;                        //!< If true, this model has vertical winds.
};

} // namespace
