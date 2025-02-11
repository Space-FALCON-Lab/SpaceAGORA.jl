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

struct VHighData;
typedef struct VHighData VIRAHighData;

//! \brief The VIRA model of the upper Venus Atmosphere.
//!
//! The VIRA model of the upper Venus Atmosphere.
//! \ingroup VenusGRAM
class VenusIRAHigh : public Atmosphere, public VenusCommon
{
public:
  VenusIRAHigh();
  VenusIRAHigh(const VenusIRAHigh& orig) = default;
  virtual ~VenusIRAHigh() override = default;

  void update() override;

  void initializeData();
  greal setLowerHeightIndex();
  greal setUpperHeightIndex();
  void getInterpolationLimits(greal& a, greal& b, greal& y, greal& z);

  friend class VenusIRA;

private:
  size_t heightIndex = 0;  //!< Height index for data table lookups.

  static const size_t heightSize = 21;  //!< Size of the height dimension.
  static const size_t szaSize = 7;      //!< Size of the solar zenith angle dimension.
  static const VIRAHighData data[heightSize * szaSize];  //!< VIRA upper atmosphere data.

  static bool initialized;                         //!< Flags when data has been initialized.
  static std::vector<greal> zhi;                   //!< Height layers for upper data \units{km}.
  static std::vector<greal> vsza;                  //!< Solar zenith angle layers for upper data \units{degrees}.
  static std::vector<std::vector<greal>> dhi;      //!< Density values for upper data \units{kg/m^3}.
  static std::vector<std::vector<greal>> phi;      //!< Pressure values for upper data \units{Pa}.
  static std::vector<std::vector<greal>> thi;      //!< Temperature values for upper data \units{K}.
  static std::vector<std::vector<greal>> yco2hi;   //!< Number density of CO2 for upper data.
  static std::vector<std::vector<greal>> yn2hi;    //!< Number density of N2 for upper data.
  static std::vector<std::vector<greal>> yohi;     //!< Number density of O for upper data.
  static std::vector<std::vector<greal>> ycohi;    //!< Number density of CO for upper data.
  static std::vector<std::vector<greal>> yhehi;    //!< Number density of He for upper data.
  static std::vector<std::vector<greal>> ynhi;     //!< Number density of N for upper data.
  static std::vector<std::vector<greal>> yhhi;     //!< Number density of H for upper data.


};

} // namespace

