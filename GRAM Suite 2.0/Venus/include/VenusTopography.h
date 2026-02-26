//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Venus-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"
#include "Topography.h"

namespace GRAM {

//! \brief This is a VenusTopography.
//!
//! This is a VenusTopography.
//! \ingroup VenusGRAM
class VenusTopography : public Topography
{
public:
  VenusTopography();
  VenusTopography(const VenusTopography& orig) = delete;
  virtual ~VenusTopography() = default;

  double getTopographicHeight(greal lat, greal lon) override;

protected:

private:
  static const char topoFile[];
};

} // namespace
