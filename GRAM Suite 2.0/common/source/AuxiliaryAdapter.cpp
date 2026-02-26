//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include "AuxiliaryAdapter.h"
#include "Position.h"
#include "AuxiliaryAtmosphere.h"

namespace GRAM {

//! \copydoc Atmosphere::Atmosphere()
AuxiliaryAdapter::AuxiliaryAdapter()
{
}

//! \copydoc Atmosphere::Atmosphere(const Atmosphere& orig)
AuxiliaryAdapter::AuxiliaryAdapter(const AuxiliaryAdapter& orig)
{
  auxList = orig.auxList;
  hasStandardDeviations = orig.hasStandardDeviations;
}

//! \copydoc Atmosphere::~Atmosphere()
AuxiliaryAdapter::~AuxiliaryAdapter()
{
  auxList.clear();
}

//! \copydoc PerturbedAtmosphere::setInputParameters()
void AuxiliaryAdapter::setInputParameters(const InputParameters& params)
{
  size_t i = 0;
  while (i < params.auxiliaryAtmosphereFileName.size() 
    && !params.auxiliaryAtmosphereFileName[i].empty()) {
    addAuxiliaryAtmosphere(params.auxiliaryAtmosphereFileName[i],
      params.innerRadius[i], params.outerRadius[i], params.isEastLongitudePositiveOnInput);
    ++i;
  }
}

//! \brief Adds an auxiliary atmosphere to the list.
//!
//! Use this function to append to the list of auxiliary atmospheres.
//! Fairing between the current atmospheric state and auxiliary atmospheres
//! is performed in the order in which the atmospheres are added to the list.
//! \param dataFileName The name of the file containing the auxiliary atmosphere data.
//! \param innerRadius The inner radius of the fairing region.
//! \param outerRadius The outer radius of the fairing region.
//! \param eastLonPositive True if data file uses east longitude positive convention.
void AuxiliaryAdapter::addAuxiliaryAtmosphere(const std::string& dataFileName, greal innerRadius, greal outerRadius, bool eastLonPositive)
{
  AuxiliaryAtmosphere auxAtmos;
  auxAtmos.setDataFile(dataFileName);
  auxAtmos.setBounds(innerRadius, outerRadius);
  auxAtmos.setEastLongitudePositive(eastLonPositive);
  auxAtmos.setSDFlag(hasStandardDeviations);
  addAuxiliaryAtmosphere(auxAtmos);
}

//! \brief Adds an auxiliary atmosphere to the list.
//!
//! Use this function to append to the list of auxiliary atmospheres.
//! Fairing between the current atmospheric state and auxiliary atmospheres
//! is performed in the order in which the atmospheres are added to the list.
//! \param auxAtmos An auxiliary atmosphere.
void AuxiliaryAdapter::addAuxiliaryAtmosphere(AuxiliaryAtmosphere& auxAtmos)
{
  auxAtmos.setAuxId(auxList.size());
  auxList.push_back(auxAtmos);
  auxList.back().loadData();
}

AuxiliaryAtmosphere& AuxiliaryAdapter::getAuxiliaryAtmosphere(size_t index)
{
  if (index > auxList.size() - 1) {
    // Error
  }
  return auxList[index];
}

//! \brief Performs an update() on each auxiliary atmosphere.
//!
//! In the order in which atmospheres have been added to the list, each atmosphere
//! gets the position and the current state (as recieved from the previous atmosphere).
//! The atmosphere is updated and the current state is retrieved.
//! \param position The current position.
//! \param atmos The current atmospheric state.
//! \retval The updated atmospheric state.
void AuxiliaryAdapter::updateAuxiliaryAtmospheres(const Position& position, AtmosphereState& atmos)
{
  for (auto& aux : auxList) {
    aux.setPosition(position);
    aux.setInputState(atmos);
    aux.update();
    atmos = aux.getAtmosphereState();
  }
}

void AuxiliaryAdapter::updateAuxiliaryStandardDeviations(const Position& position, AtmosphereState& atmos)
{
  for (auto& aux : auxList) {
    aux.setPosition(position);
    aux.setInputState(atmos);
    aux.updateStandardDeviations();
    atmos = aux.getAtmosphereState();
  }
}

} // namespace
