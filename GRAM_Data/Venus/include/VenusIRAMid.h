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

struct VMidData;
typedef struct VMidData VIRAMidData;

//! \brief The VIRA model of the middle Venus Atmosphere.
//!
//! The VIRA model of the middle Venus Atmosphere.
//! \ingroup VenusGRAM
class VenusIRAMid : public Atmosphere, public VenusCommon
{
public:
  VenusIRAMid();
  VenusIRAMid(const VenusIRAMid& orig) = default;
  virtual ~VenusIRAMid() override = default;

  void update() override;

  void initializeData();
  greal setLowerHeightIndex();
  greal setUpperHeightIndex();
  void getInterpolationLimits(greal& a, greal& b, greal& y, greal& z);

  friend class VenusIRA;

private:
  greal getLinearMean(greal noonVal, greal midnightVal, greal amplitudeFactor);
  greal getLogMean(greal noonVal, greal midnightVal, greal amplitudeFactor);

  size_t heightIndex = 0;  //!< Height index for data table lookups.

  static const size_t heightSize = 11;  //!< Size of the height dimension.
  static const size_t lstSize = 2;      //!< Size of the solar time dimension.
  static const VIRAMidData data[heightSize * lstSize];  //!< VIRA middle atmosphere data.

  static bool initialized;                        //!< Flags when data has been initialized.
  static std::vector<greal> zmd;                  //!< Height layers for middle data \units{km}.
  static std::vector<greal> lstmd;                //!< Solar time layers for middle data \units{hr}.
  static std::vector<std::vector<greal>> dmd;     //!< Density values for middle data \units{kg/m^3}.
  static std::vector<std::vector<greal>> pmd;     //!< Pressure values for middle data \units{Pa}.
  static std::vector<std::vector<greal>> tmd;     //!< Temperature values for middle data \units{K}.
  static std::vector<std::vector<greal>> yco2md;  //!< Number density of CO2 for middle data.
  static std::vector<std::vector<greal>> yn2md;   //!< Number density of N2 for middle data.
  static std::vector<std::vector<greal>> yomd;    //!< Number density of O for middle data.
  static std::vector<std::vector<greal>> ycomd;   //!< Number density of CO for middle data.
  static std::vector<std::vector<greal>> yhemd;   //!< Number density of He for middle data.
  static std::vector<std::vector<greal>> ynmd;    //!< Number density of N for middle data.
};

} // namespace

