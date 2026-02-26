%module gram

%{
#include "Position.h"
#include "EphemerisState.h"
#include "AtmosphereState.h"
#include "ConstituentGas.h"
#include "gram.h"
#include "Many.h"
#include "InputParameters.h"
#include "Atmosphere.h"
#include "AuxiliaryAdapter.h"
#include "PerturbedAtmosphere.h"
#include "NamelistReader.h"
#include "GramTime.h"
#include "PerturbationState.h"
#include "AuxiliaryAtmosphere.h"
#include "SpiceLoader.h"
#include "Ephemeris.h"

%}

%include <std_string.i>
%include "gram.h"
%import "Many.h"
%include "ConstituentGas.h"
%include "Position.h"
%include "EphemerisState.h"
%include "AtmosphereState.h"
%include "InputParameters.h"
%include "Atmosphere.h"
%import "AuxiliaryAtmosphere.h"
%include "AuxiliaryAdapter.h"
%include "GramTime.h"
%include "PerturbationState.h"
%include "PerturbedAtmosphere.h"
%include "NamelistReader.h"
%include "SpiceLoader.h"
%include "Ephemeris.h"


%{
#include "EarthAtmosphere.h"
#include "EarthAtmosphereState.h"
#include "EarthCommon.h"
#include "EarthInputParameters.h"
#include "EarthNamelistReader.h"
%}

%include "EarthCommon.h"
%include "EarthInputParameters.h"
%include "EarthNamelistReader.h"
%include "EarthAtmosphereState.h"
%include "EarthAtmosphere.h"

%{
#include "JupiterAtmosphere.h"
#include "JupiterInputParameters.h"
#include "JupiterCommon.h"
#include "JupiterNamelistReader.h"
%}

%include "JupiterCommon.h"
%include "JupiterInputParameters.h"
%include "JupiterAtmosphere.h"
%include "JupiterNamelistReader.h"

%{
#include "MarsAtmosphere.h"    
#include "MarsAtmosphereState.h"
#include "MarsCommon.h"
#include "MarsNamelistReader.h"
#include "MarsInputParameters.h"
%}

%include "MarsCommon.h"
%include "MarsInputParameters.h"
%include "MarsAtmosphereState.h"
%include "MarsAtmosphere.h"
%include "MarsNamelistReader.h"

%{
#include "NeptuneAtmosphere.h"
#include "NeptuneInputParameters.h"
#include "NeptuneNamelistReader.h"
#include "NeptuneCommon.h"
%}

%include "NeptuneCommon.h"
%include "NeptuneInputParameters.h"
%include "NeptuneAtmosphere.h"
%include "NeptuneNamelistReader.h"

%{
#include "TitanAtmosphere.h"
#include "TitanCommon.h"
#include "TitanInputParameters.h"
#include "TitanNamelistReader.h"
%}

%include "TitanCommon.h"
%include "TitanInputParameters.h"
%include "TitanNamelistReader.h"
%include "TitanAtmosphere.h"


%{
#include "UranusAtmosphere.h"
#include "UranusCommon.h"
#include "UranusInputParameters.h"
#include "UranusNamelistReader.h"
%}

%include "UranusCommon.h"
%include "UranusInputParameters.h"
%include "UranusNamelistReader.h"
%include "UranusAtmosphere.h"

%{
#include "VenusAtmosphere.h"    
#include "VenusCommon.h"
#include "VenusInputParameters.h"
#include "VenusNamelistReader.h"
%}

%include "VenusCommon.h"
%include "VenusInputParameters.h"
%include "VenusNamelistReader.h"
%include "VenusAtmosphere.h"