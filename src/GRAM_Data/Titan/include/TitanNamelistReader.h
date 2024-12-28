//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Titan-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"
#include "NamelistReader.h"
#include "TitanInputParameters.h"

namespace GRAM {

//! \brief This is a FORTRAN namelist reader for the Titan atmosphere.
//!
//! The legacy GRAM models were written in FORTRAN and used a language specific feature, known as
//! namelist files, to ingest parameters.  This class supports the form of namelists as given in
//! the legacy examples.  Note that new FORTRAN versions have expanded the definition of namelists.
//! This class does not attempt to support all namelist formats.
//! \ingroup TitanGRAM
class TitanNamelistReader : public NamelistReader
{
public:
  TitanNamelistReader();
  TitanNamelistReader(const TitanNamelistReader& orig) = delete;
  virtual ~TitanNamelistReader() = default;

  void getParameters(const std::string& fileName, TitanInputParameters& inputParameters);

protected:

private:
};

} // namespace
