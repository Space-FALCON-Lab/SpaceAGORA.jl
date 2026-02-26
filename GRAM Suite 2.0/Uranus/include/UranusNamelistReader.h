//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Uranus-GRAM
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "gram.h"
#include "NamelistReader.h"
#include "UranusInputParameters.h"

namespace GRAM {

//! \brief This is a FORTRAN namelist reader for the Uranus atmosphere.
//!
//! The legacy GRAM models were written in FORTRAN and used a language specific feature, known as
//! namelist files, to ingest parameters.  This class supports the form of namelists as given in
//! the legacy examples.  Note that new FORTRAN versions have expanded the definition of namelists.
//! This class does not attempt to support all namelist formats.
//! \ingroup UranusGRAM
class UranusNamelistReader : public NamelistReader
{
public:
  UranusNamelistReader();
  UranusNamelistReader(const UranusNamelistReader& orig) = delete;
  virtual ~UranusNamelistReader() = default;

  void getParameters(const std::string& fileName, UranusInputParameters& inputParameters);

protected:

private:
};

} // namespace
