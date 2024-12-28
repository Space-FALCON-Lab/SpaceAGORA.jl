//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: mars-GRAM
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <string>
#include "MarsAtmosphere.h"
#include "MarsNamelistReader.h" 

using namespace std;
using namespace GRAM;

// Helper function for printing
void print(const string& label, greal value1, greal value2)
{
  cout << left << setw(30) << label << right << setw(13) << value1 << "   " << setw(13) << value2 << endl;
}

int main(int argc, char** argv)
{
  cout << endl;
  cout << "=============" << endl;
  cout << "A C++ Example" << endl;
  cout << "=============" << endl;
  cout << endl;

  try {
    // Set the Spice data path (this is critical).
    MarsInputParameters inputParameters;
    MarsNamelistReader reader;
    inputParameters.dataPath = "./data/";
    reader.tryGetSpicePath(inputParameters);

    // Create a mars model and set default parameters.
    MarsAtmosphere mars;
    mars.setInputParameters(inputParameters);

    // More than one MarsAtmosphere can be created.
    MarsAtmosphere mars2;
    mars2.setInputParameters(inputParameters);

    // Get and print the version string.
    string version = mars.getVersionString();
    cout << version << endl;
    cout << endl;

    // Set a mars specific parameter
    mars.setMapYear(0);
    mars.setMGCMDustLevels(2.0, 0.0, 0.0);
    mars2.setMapYear(2);

    // Set the perturbation scale factor, relative step size, and the initial seed.
    mars.setPerturbationScales(1.5);
    mars.setMinRelativeStepSize(0.5);
    mars.setSeed(1001);
    mars2.setPerturbationScales(1.5);
    mars2.setMinRelativeStepSize(0.5);
    mars2.setSeed(1001);

    // Set the start time of the trajectory
    GramTime ttime;
    ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT);
    mars.setStartTime(ttime);
    mars2.setStartTime(ttime);

    // Set the position
    Position position;
    position.height = 2.0_km;
    position.latitude = 22.0_deg;
    position.longitude = 48.0_deg;
    position.elapsedTime = 100.0_sec;
    position.isPlanetoCentric = true;
    mars.setPosition(position);
    mars2.setPosition(position);

    // Update the atmosphere data
    mars.update();
    mars2.update();

    // The position and other output is returned in the Position class.
    const Position& pos = mars.getPosition();
    const Position& pos2 = mars2.getPosition();
    cout << "                                 Mars #1         Mars #2" << scientific << endl;
    print("Total Radius: ", pos.totalRadius, pos2.totalRadius);
    print("Gravity: ", pos.gravity, pos2.gravity);
    cout << endl;

    // The ephemeris values are returned in the EphemerisState class.
    const EphemerisState& estate = mars.getEphemerisState();
    const EphemerisState& estate2 = mars2.getEphemerisState();
    print("Solar Time:", estate.solarTime, estate2.solarTime);
    print("Longitude of the Sun:", estate.longitudeSun, estate2.longitudeSun);
    cout << endl;

    // The atmosphere values are returned in the AtmosphereState class.
    const AtmosphereState& atmos = mars.getAtmosphereState();
    const AtmosphereState& atmos2 = mars2.getAtmosphereState();
    print("Temperature:", atmos.temperature, atmos2.temperature);
    print("Pressure:", atmos.pressure, atmos2.pressure);
    print("Density:", atmos.density, atmos2.density);
    print("Pressure Scale Height:", atmos.pressureScaleHeight, atmos2.pressureScaleHeight);
    print("Density Scale Height:", atmos.densityScaleHeight, atmos2.densityScaleHeight);
    cout << endl;

    // Print perturbed density
    print("Mean Density:", atmos.density, atmos2.density);
    print("Perturbed Density:", atmos.perturbedDensity, atmos2.perturbedDensity);
    print("Perturbation Percent:", atmos.densityPerturbation * 100.0, atmos2.densityPerturbation * 100.0);
    cout << endl;

    // Get and print gases.
    print("Average Molecular Weight:", atmos.averageMolecularWeight, atmos2.averageMolecularWeight);
    print("Carbon Dioxide Mole Fraction:", atmos.carbonDioxide.moleFraction * 100.0, atmos2.carbonDioxide.moleFraction * 100.0);
    print("Dinitrogen Mole Fraction:", atmos.dinitrogen.moleFraction * 100.0, atmos2.dinitrogen.moleFraction * 100.0);
    print("Argon Mole Fraction:", atmos.argon.moleFraction * 100.0, atmos2.argon.moleFraction * 100.0);
    print("Carbon Monoxide Mole Fraction:", atmos.carbonMonoxide.moleFraction * 100.0, atmos2.carbonMonoxide.moleFraction * 100.0);
    print("Dioxygen Mole Fraction:", atmos.dioxygen.moleFraction * 100.0, atmos2.dioxygen.moleFraction * 100.0);
    print("Dihydrogen Mole Fraction:", atmos.dihydrogen.moleFraction * 100.0, atmos2.dihydrogen.moleFraction * 100.0);
    print("Hydrogen Mole Fraction:", atmos.hydrogen.moleFraction * 100.0, atmos2.hydrogen.moleFraction * 100.0);
    print("Oxygen Mole Fraction:", atmos.oxygen.moleFraction * 100.0, atmos2.oxygen.moleFraction * 100.0);
    print("Helium Mole Fraction:", atmos.helium.moleFraction * 100.0, atmos2.helium.moleFraction * 100.0);
    print("Water Vapor Mole Fraction:", atmos.water.moleFraction * 100.0, atmos2.water.moleFraction * 100.0);
    cout << endl;

    // Print wind data
    print("EW Wind:", atmos.ewWind, atmos2.ewWind);
    print("NS Wind:", atmos.nsWind, atmos2.nsWind);
    print("Vertical Wind:", atmos.verticalWind, atmos2.verticalWind);
    print("Perturbed EW Wind:", atmos.perturbedEWWind, atmos2.perturbedEWWind);
    print("Perturbed NS Wind:", atmos.perturbedNSWind, atmos2.perturbedNSWind);
    cout << endl;

    // Print Mars specific metrics.
    const MarsAtmosphereState& marsAtmos = mars.getMarsAtmosphereState();
    const MarsAtmosphereState& marsAtmos2 = mars2.getMarsAtmosphereState();
    print("Dust Optical Depth:", marsAtmos.dustOpticalDepth, marsAtmos2.dustOpticalDepth);
    print("Daily Mean Temperature:", marsAtmos.temperatureDaily, marsAtmos2.temperatureDaily);
    print("Daily Max Temperature:", marsAtmos.temperatureMax, marsAtmos2.temperatureMax);
    print("Daily Min Temperature:", marsAtmos.temperatureMin, marsAtmos2.temperatureMin);
  }
  catch (const string& msg) {
    cerr << msg << endl;
  }
  catch (const std::runtime_error& err) {
    cerr << "An unanticipated error occurred." << endl;
    cerr << err.what() << endl;
  }


  return 0;
}

