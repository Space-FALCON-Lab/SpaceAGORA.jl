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
#include <vector>
#include "gram.h"
#include "AuxiliaryAtmosphere.h"
#include "Position.h"

namespace GRAM {

//! \brief A container interface for multiple auxiliary atmospheres.
//!
//! Inheriting the AuxiliaryAdapter class will give an atmosphere class the ability
//! to add and utilize mulitple auxiliary atmospheres.  The order in which AuxiliaryAtmospheres
//! are added, with addAuxiliaryAtmosphere(), is important.  If A is added before B, then the
//! current AtmosphereState is faired with A. The result is then faired with B.  The fairings are
//! performed by calling updateAuxiliaryAtmosphere().
//! \ingroup CommonGRAM
class AuxiliaryAdapter
{
public:
  AuxiliaryAdapter();
  AuxiliaryAdapter(const AuxiliaryAdapter& orig);
  virtual ~AuxiliaryAdapter();

  void setInputParameters(const InputParameters& params);
  void setSDFlag(bool dataHasStandardDeviations) { hasStandardDeviations = dataHasStandardDeviations; }

  void addAuxiliaryAtmosphere(const std::string& dataFileName, greal innerRadius, greal outerRadius, bool eastLonPositive);
  void addAuxiliaryAtmosphere(AuxiliaryAtmosphere& auxAtmos);

  AuxiliaryAtmosphere& getAuxiliaryAtmosphere(size_t index);

  bool hasAuxiliaryAtmopheres() { return (!auxList.empty()); }

protected:
  void updateAuxiliaryAtmospheres(const Position& position, AtmosphereState& atmos);
  void updateAuxiliaryStandardDeviations(const Position& position, AtmosphereState& atmos);

private:
  std::vector<AuxiliaryAtmosphere> auxList;  //!< A list of auxiliary atmospheres.
  bool hasStandardDeviations = false;
};

} // namespace
