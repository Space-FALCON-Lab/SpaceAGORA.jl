//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include "gram.h"
#include "VenusCommon.h"
#include "Atmosphere.h"

namespace GRAM {

struct VLowData;
typedef struct VLowData VIRALowData;

//! \brief The VIRA model of the lower Venus Atmosphere.
//!
//! The VIRA model of the lower Venus Atmosphere.
//! \ingroup VenusGRAM
class VenusIRALow : public Atmosphere, public VenusCommon
{
public:
  VenusIRALow();
  VenusIRALow(const VenusIRALow& orig) = default;
  virtual ~VenusIRALow() override = default;

  void update() override;

  void initializeData();

  greal setLowerHeightIndex();
  greal setUpperHeightIndex();
  void getInterpolationLimits(greal& a, greal& b, greal& y, greal& z);

  // This is needed so reference values can be initialized.
  friend class VenusIRA;

private:
  size_t heightIndex = 0;  //!< Height index for data table lookups.

  static const size_t heightSize = 81;  //!< Size of the height dimension.
  static const size_t latSize = 5;      //!< Size of the latitude dimension.
  static const VIRALowData data[heightSize * latSize];  //!< VIRA lower atmosphere data.

  static bool initialized;                        //!< Flags when data has been initialized.
  static std::vector<greal> zlo;                  //!< Height layers for lower data \units{km}.
  static std::vector<greal> vlat;                 //!< Latitude layers for lower data \units{degrees}.
  static std::vector<std::vector<greal>> dlo;     //!< Density values for lower data \units{kg/m^3}.
  static std::vector<std::vector<greal>> plo;     //!< Pressure values for lower data \units{Pa}.
  static std::vector<std::vector<greal>> tlo;     //!< Temperature values for lower data \units{K}.
  static std::vector<std::vector<greal>> zetalo;  //!< Compressibility factor for lower data.
  static std::vector<std::vector<greal>> Rlo;     //!< Specific gas constants for lower data \units{J/(kg K)}.
};

} // namespace

