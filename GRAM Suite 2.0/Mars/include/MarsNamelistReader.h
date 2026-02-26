//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Mars-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"
#include "NamelistReader.h"
#include "MarsInputParameters.h"

namespace GRAM {

//! \brief This is a FORTRAN namelist reader for the Mars atmosphere.
//!
//! The legacy GRAM models were written in FORTRAN and used a language specific feature, known as
//! namelist files, to ingest parameters.  This class supports the form of namelists as given in
//! the legacy examples.  Note that new FORTRAN versions have expanded the definition of namelists.
//! This class does not attempt to support all namelist formats.
//! \ingroup MarsGRAM
class MarsNamelistReader : public NamelistReader
{
public:
  MarsNamelistReader();
  MarsNamelistReader(const MarsNamelistReader& orig) = delete;
  virtual ~MarsNamelistReader() = default;

  void getParameters(const std::string& fileName, MarsInputParameters& inputParameters);

protected:

private:
};

} // namespace
