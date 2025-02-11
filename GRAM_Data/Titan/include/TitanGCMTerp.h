//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Titan-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include "gram.h"
#include "Atmosphere.h"
#include "TitanCommon.h"

namespace GRAM {

struct TGCMData;
typedef struct TGCMData TitanGCMData;
struct TN2Data;
typedef struct TN2Data TitanN2Data;

//! \brief The GCM Terp model of a Titan Atmosphere.
//!
//! This model is for heights from the surface to the 0.3 mbar level.
//! \ingroup TitanGRAM
class TitanGCMTerp : public Atmosphere, public TitanCommon
{
public:
  TitanGCMTerp();
  TitanGCMTerp(const TitanGCMTerp& orig) = default;
  virtual ~TitanGCMTerp() override = default;

  void setMethaneMoleFraction(greal mmf) { methane.moleFraction = mmf; }

  void updateGCMTables();
  void update() override;

  greal getHeightAt03mb() const { return height03mb; }
  greal getTemperatureAt03mb() const { return temperature03mb; }

private:
  void initializeData();
  greal getNitrogenCompressibiltyFactor(greal temp, greal pres);
  void LagrangeQuad(greal x1, greal x2, greal x3, greal y1, greal y2, greal y3, greal x, greal& y);
  void LagrangeLin(greal x1, greal x2, greal y1, greal y2, greal x, greal& y);


  static bool initialized;                                //!< \brief Data initialization flag.
  static const int latSize = 19;                          //!< \brief Number of latitudes.
  static const int presSize = 31;                         //!< \brief Number of pressure levels.
  static const int zPresSize = 6;                         //!< \brief Number of zeta pressures.
  static const int zTempSize = 16;                        //!< \brief Number of zeta temperatures.
  static const TitanGCMData data[4 * presSize];           //!< \brief  Input temperature data table.
  static const TitanN2Data n2data[zTempSize * zPresSize]; //!< \brief Input nitrogen data table

  // Table Data
  static greal gcmPressure[presSize];             //!< \brief GCM data pressure levels (N / m^2)
  static greal temperature000[presSize][latSize]; //!< \brief Temperature table for Ls = 0
  static greal temperature270[presSize][latSize]; //!< \brief Temperature table for Ls = 270
  static greal ewWind000[presSize][latSize];      //!< \brief East/west wind table for Ls = 0
  static greal ewWind270[presSize][latSize];      //!< \brief East/west wind table for Ls = 270
  static greal zetaPressure[zPresSize];           //!< \brief Pressure values (atmospheres) for tabular Zeta data
  static greal zeta[zTempSize][zPresSize];        //!< \brief Nitrogen compressibilty factor table

  // Table Data (generated)
  greal gcmHeight[presSize] = {};               //!< \brief GCM height table.
  greal gcmTemperature[presSize] = {};          //!< \brief GCM temperature table.
  greal gcmDensity[presSize] = {};              //!< \brief GCM density table.
  greal gcmEwWind[presSize] = {};               //!< \brief GCM east/west wind table.
  greal gcmPressureScaleHeight[presSize] = {};  //!< \brief GCM pressure scale height table.
  greal height03mb = 0.0;                       //!< \brief GCM height at 0.3 mb pressure
  greal temperature03mb = 0.0;                  //!< \brief GCM temperature at 0.3 mb pressure

  // Constants
  const greal avgMolWt = 27.808;         //!< \brief Constant molecular weight over GCM height range
  const greal mfNitrogenDefault = 0.95;  //!< \brief Default mole fraction for Nitrogen (N2).
  const greal mfMethaneDefault = 0.03;   //!< \brief Default mole fraction for Methane.
  const greal mfArgonDefault = 0.02;     //!< \brief Default mole fraction for Argon.

};

} // namespace
