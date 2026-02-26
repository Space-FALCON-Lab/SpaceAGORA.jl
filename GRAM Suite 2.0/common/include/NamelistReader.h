//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "gram.h"
#include "InputParameters.h"

namespace GRAM {

//! \brief The NamelistReader class provides support for FORTRAN namelist files.
//!
//! The legacy GRAM models were written in FORTRAN and used a language specific feature, known as
//! namelist files, to ingest parameters.  This class supports the form of namelists as given in
//! the legacy examples.  Note that new FORTRAN versions have expanded the definition of namelists.
//! This class does not attempt to support all namelist formats.
//!
//! This class is base class that must be subclassed by each model.  Since each model may also
//! subclass the InputParameters, it is up to the subclass to provide an interface for that type.
//! \ingroup CommonGRAM
class NamelistReader
{
public:
  NamelistReader();
  NamelistReader(const NamelistReader& orig) = delete;
  virtual ~NamelistReader() = default;

  void tryGetSpicePath(std::string& spicePath);
  void tryGetSpicePath(std::string& spicePath, std::string& dataPath);
  void tryGetSpicePath(InputParameters& params);
  std::string getWarnings() { return warnings.str(); }

protected:
  void readFile(const std::string& fileName);
  void parseLine(const std::string& lineBuffer);
  void getParameters(InputParameters& inputParameters);
  void checkForErrors();

  bool getItem(const std::string& itemName, int& value);
  bool getItem(const std::string& itemName, double& value);
  bool getItem(const std::string& itemName, float& value);
  bool getItem(const std::string& itemName, bool& value);
  bool getItem(const std::string& itemName, std::string& value);
  const std::string& findItem(const std::string& itemName);

  void processPath(std::string& path);

private:
  bool dateHasError(const InputParameters& inputParameters, std::string msg);
  std::map<std::string, std::string> itemMap;  //!< A map of key(name)-value strings.
  std::vector<std::string> validKeyNames;      //!< A list of valid keys sought by the reader.
  std::ostringstream warnings;                 //!< A stream for warnings issued by the reader.
};

} // namespace
