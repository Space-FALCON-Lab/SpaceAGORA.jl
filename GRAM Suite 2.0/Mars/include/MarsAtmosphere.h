//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once
#ifdef _MSC_VER
#pragma warning( disable : 4250 )
#endif

#include <vector>
#include "MarsCommon.h"
#include "PerturbedAtmosphere.h"
#include "Ephemeris.h"
#include "MarsInputParameters.h"
#include "MarsGCM.h"
#include "MarsAtmosphereState.h"

//! \defgroup Cpp_Mars The C++ Interface for Mars
//! \brief Class references for the MarsGRAM interface.
//!
//! These are the classes needed for building an application that interfaces with MarsGRAM.
//! The interface for Mars is declared in the following:
//! \code #include "MarsAtmosphere.h" \endcode 
//! An example using this class can be found \ref Mars/examples/Mars.cpp "here".

//! \defgroup MarsGRAM The Full List of MarsGRAM Classes.
//! \brief Class references for all Mars models.
//!
//! These are the Mars specific classes needed for building MarsGRAM.

namespace GRAM {

//! \brief The primary interface for the Mars atmosphere model.
//!
//! The Mars atmosphere atmosphere model.
//!
//! \ingroup Cpp_Mars MarsGRAM
class MarsAtmosphere : public PerturbedAtmosphere, public MarsCommon
{
//! \example Mars/examples/Mars.cpp
public:
  MarsAtmosphere();
  MarsAtmosphere(const MarsAtmosphere& orig);
  virtual ~MarsAtmosphere() override;

  static const std::string& getVersionString();
  void setInputParameters(const MarsInputParameters& params);

  void setPlanetaryRadii(greal eqrad, greal prad);
  void setDataPath(const std::string& path);
  void setMapYear(int year);
  void setHeightOffsetModel(MarsOffsetModel model, greal hgtOffset = 0.0);
  void setHeightAboveSurface(greal heightAboveSurface);
  void setMGCMDustLevels(greal constantDustLevel, greal minDustLevel, greal maxDustLevel);
  void setF107(greal f107);
  void setPerturbationWaveLengthScale(greal scale);
  void setMOLAHeights(bool isMola);
  void setMinMax(bool minMax);
  void setDustStorm(greal lonSun, greal duration, greal intensity, greal maxRadius, greal lat, greal lon);
  void setDustDensity(greal nu, greal diameter, greal dens);
  void setExosphericTemperature(greal offset, greal factor);
  void setWaveDefaults(greal date, greal scale, greal mean, greal a1, greal p1, greal r1, greal a2, greal p2, greal r2, greal a3, greal p3, greal r3);
  void setWaveFile(const std::string& waveFile);
  void setWindScales(greal meanWinds, greal boundaryLayerWinds);
  void setAreoidRadiusCallback(TopoCallback callback) { getAreoidRadiusCallback = callback; }
  void setTopographicHeightCallback(TopoCallback callback) { getTopographicHeightCallback = callback; }
  void setCallbackData(void* dataPointer) { callbackDataPointer = dataPointer; }

  void update() override;

  const InputParameters& getInputParameters() const { return inputParameters; }
  const MarsAtmosphereState& getMarsAtmosphereState() const { return marsAtmos; }

private:
  void updatePosition() override;
  void updateReferenceValues() override;
  void updatePressureAtSurface() override;
  void updateMetrics() override;

  void getPerturbationFactors(greal& pertLow, greal& pertHigh) override;
  void getScaleParameters(greal& verticalScale, greal& horizontalScale) override;
  void getWindDeviations(greal& ewStdDev, greal& nsStdDev, greal& vertStdDev) override;

  void updateSlopeWinds();
  void initializeWaveData(const MarsInputParameters& params);
  void updateWavePertubation();
  void applyWavePerturbation();
  inline greal cp(greal t) { return 639.5 + t * (0.123687 + t * 0.00200225); }

  MarsGCM* marsgcmPtr = NULL;           //!< Pointer to MGCM or TesGCM model.
  MarsInputParameters inputParameters;  //!< User supplied parameters.
  Position refPosition;                 //!< Position relative to the reference ellipsoid.
  MarsAtmosphereState marsAtmos;        //!< Mars specific metrics added to the AtmosphereState.

  // Input Parameters
  bool isPlanetoCentric = true; //!< True for planeto-centric inputs, false for planeto-graphic.
  bool isMolaHeights = true;    //!< True for input heights relative to MOLA areoid, false for relative to reference ellipsoid.
  bool computeMinMax = true;    //!< True to compute daily extrema, false for none.
  greal waveScale = 20.0;       //!< Vertical scale \units{km} of longitude-dependent wave damping at altitudes below 100 km (10 < Wscale < 10, 000 km).
  greal dustNu = 0;             //!< Parameter for vertical distribution of dust density.
  greal dustDensity = 0;        //!< Dust particle density \units{ kg/m^3 }.
  greal dustDiameter = 0;       //!< Dust particle diameter (micrometers, assumed monodisperse).
  greal perturbationWaveLengthScale = 1.0; //!< Scale factor for perturbation wavelengths (0.1-10).

  struct WaveData {
    greal waveDate;        //!< Julian day for primary peak(s) of LDW traveling component.
    greal waveTime;        //!< Seconds from waveDate.
    greal waveMeanOffset;  //!< Diurnal mean value of longitude-dependent wave (LDW).
    greal waveAmplitude1;  //!< Amplitude of the wave-1 LDW component.
    greal wavePhase1;      //!< Phase (degrees) of the wave-1 LDW component.
    greal wavePhase1Rate;  //!< Rate of movement of wave-1 LDW peak.
    greal waveAmplitude2;  //!< Amplitude of the wave-1 LDW component.
    greal wavePhase2;      //!< Phase (degrees) of the wave-1 LDW component.
    greal wavePhase2Rate;  //!< Rate of movement of wave-1 LDW peak.
    greal waveAmplitude3;  //!< Amplitude of the wave-1 LDW component.
    greal wavePhase3;      //!< Phase (degrees) of the wave-1 LDW component.
    greal wavePhase3Rate;  //!< Rate of movement of wave-1 LDW peak.
  };

  WaveData waveData;                          //!< Internal variable to hold an instance of wave data.
  static bool dataLoaded;                     //!< Flags whether wave data has been loaded.
  static std::vector<WaveData> waveDataList;  //!< Data provided by a wave file.

  TopoCallback getAreoidRadiusCallback = nullptr;      //!< Override of MOLATopography::getAreoidRadius.
  TopoCallback getTopographicHeightCallback = nullptr; //!< Override of MOLATopography::getTopographicHeight.
  void* callbackDataPointer = nullptr;                 //!< Data pointer for overrides.

  // Convenience references
  greal& wavePerturbation = marsAtmos.wavePerturbation;             //!< The wave perturbation factor.
  greal& albedo = marsAtmos.albedo;                                 //!< \copydoc MarsAtmosphereState::albedo
  greal& dustColumnArealDensity = marsAtmos.dustColumnArealDensity; //!< \copydoc MarsAtmosphereState::dustColumnArealDensity
  greal& dustMixingRatio = marsAtmos.dustMixingRatio;               //!< \copydoc MarsAtmosphereState::dustMixingRatio
  greal& dustMassDensity = marsAtmos.dustMassDensity;               //!< \copydoc MarsAtmosphereState::dustMassDensity
  greal& dustNumberDensity = marsAtmos.dustNumberDensity;           //!< \copydoc MarsAtmosphereState::dustNumberDensity
  greal& ewWindDaily = marsAtmos.ewWindDaily;                       //!< \copydoc MarsAtmosphereState::ewWindDaily
  greal& nsWindDaily = marsAtmos.nsWindDaily;                       //!< \copydoc MarsAtmosphereState::nsWindDaily
};

} // namespace
