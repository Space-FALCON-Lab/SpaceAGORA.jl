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
#include "NamelistReader.h"
#include "VenusInputParameters.h"

namespace GRAM {

//! \brief This is a VenusNamelistReader.
//!
//! This is a VenusNamelistReader.
//! \ingroup VenusGRAM
class VenusNamelistReader : public NamelistReader
{
public:
  VenusNamelistReader();
  VenusNamelistReader(const VenusNamelistReader& orig) = delete;
  virtual ~VenusNamelistReader() = default;

  void getParameters(const std::string& fileName, VenusInputParameters& inputParameters);

protected:

private:
};

} // namespace
