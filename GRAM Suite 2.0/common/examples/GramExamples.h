//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: GRAM common
//////////////////////////////////////////////////////////////////////////

#include <string>
#include "PerturbedAtmosphere.h"

extern int unitTesting(int argc, char** argv); 
extern void atmosphereExample(GRAM::PerturbedAtmosphere& atmosphere, const std::string& name);
extern void trajectoryExample(GRAM::PerturbedAtmosphere& atmosphere, const std::string& name);
extern void monteCarloExample(GRAM::PerturbedAtmosphere& atmosphere, const std::string& name);
extern void namelistExample(GRAM::PerturbedAtmosphere& atmosphere, const std::string& name);
extern void ephemerisExample(GRAM::PerturbedAtmosphere& atmosphere, const std::string& name);
extern void perturbationExample(GRAM::PerturbedAtmosphere& atmosphere, const std::string& name);
