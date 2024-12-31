//////////////////////////////////////////////////////////////////////////
// The Global Reference Atmospheric Model (GRAM) Framework
//
// No recipient of this code should forward copies outside of the United 
// States without explicit approval by NASA Marshall Space Flight Center.
// 
// Module: Earth-GRAM
//////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>
#include "EarthAtmosphere.h" 
#include "EarthNamelistReader.h" 

using namespace std;
using namespace GRAM;

// Helper function for printing
void print(const string& label, greal value1, greal value2)
{
  cout << left << setw(27) << label << right << setw(13) << value1 << "   " << setw(13) << value2 << endl;
}

int main(int argc, char** argv)
{
  typedef std::chrono::high_resolution_clock Clock;
  auto t1 = Clock::now();

  cout << endl;
  cout << "=============" << endl;
  cout << "A C++ Example" << endl;
  cout << "=============" << endl;
  cout << endl;

  try {
    // Set the Spice data path (this is critical).
    EarthInputParameters inputParameters;
    EarthNamelistReader reader;
    reader.tryGetSpicePath(inputParameters);

    // Create a Earth model and set default parameters.
    EarthAtmosphere earth;
    earth.setInputParameters(inputParameters);

    // More than one EarthAtmosphere can be created.
    EarthAtmosphere earth2;
    earth2.setInputParameters(inputParameters);

    // Get and print the version string.
    string version = earth.getVersionString();
    cout << version << endl;
    cout << endl;

    // Set the MERRA-2 hour and data extents
    earth.setMERRA2Parameters(1, 0.0, 40.0, 20.0, 60.0);
    earth2.setMERRA2Parameters(9, 0.0, 40.0, 20.0, 60.0);

    // For NCEP model, uncomment lines below
    //earth.setUseNCEP(true);
    //earth2.setUseNCEP(true);
    //earth.setNCEPParameters(9715, 1);
    //earth2.setNCEPParameters(9715, 5);

    // Set the perturbation scale factor, relative step size, and the initial seed.
    earth.setRandomPerturbationScale(1.5);
    earth.setHorizontalWindPerturbationScale(1.5);
    earth.setVerticalWindPerturbationScale(1.5);
    earth.setSeed(1001);
    earth2.setRandomPerturbationScale(1.5);
    earth2.setHorizontalWindPerturbationScale(1.5);
    earth2.setVerticalWindPerturbationScale(1.5);
    earth2.setSeed(1001);

    // Set the start time of the trajectory
    GramTime ttime;
    ttime.setStartTime(2020, 12, 15, 0, 0, 0.0, UTC, PET);
    earth.setStartTime(ttime);
    earth2.setStartTime(ttime);

    // Set the position
    Position position;
    position.height = 5.0_km;
    position.latitude = 22.0_deg;
    position.longitude = 48.0_deg;
    position.elapsedTime = 100.0_sec;
    earth.setPosition(position);
    earth2.setPosition(position);

    // Update the atmosphere data
    earth.update();
    earth2.update();

    // The position and other output is returned in the Position class.
    const Position& pos = earth.getPosition();
    const Position& pos2 = earth2.getPosition();
    cout << "                              Earth #1        Earth #2" << scientific << endl;
    print("Total Radius: ", pos.totalRadius, pos2.totalRadius);
    print("Gravity: ", pos.gravity, pos2.gravity);
    cout << endl;

    // The ephemeris values are returned in the EphemerisState class.
    const EphemerisState& estate = earth.getEphemerisState();
    const EphemerisState& estate2 = earth2.getEphemerisState();
    print("Solar Time:", estate.solarTime, estate2.solarTime);
    print("Longitude of the Sun:", estate.longitudeSun, estate2.longitudeSun);
    cout << endl;

    // The atmosphere values are returned in the AtmosphereState class.
    const AtmosphereState& atmos = earth.getAtmosphereState();
    const AtmosphereState& atmos2 = earth2.getAtmosphereState();
    print("Temperature:", atmos.temperature, atmos2.temperature);
    print("Pressure:", atmos.pressure, atmos2.pressure);
    print("Density:", atmos.density, atmos2.density);
    print("Pressure Scale Height:", atmos.pressureScaleHeight, atmos2.pressureScaleHeight);
    print("Density Scale Height:", atmos.densityScaleHeight, atmos2.densityScaleHeight);
    cout << endl;

    // Print perturbed density
    const EarthAtmosphereState& eatmos = earth.getEarthAtmosphereState();
    const EarthAtmosphereState& eatmos2 = earth2.getEarthAtmosphereState();
    print("Perturbed Density:", atmos.perturbedDensity, atmos2.perturbedDensity);
    print("Perturbed Pressure:", eatmos.perturbedPressure, eatmos2.perturbedPressure);
    print("Perturbed Temperature:", eatmos.perturbedTemperature, eatmos2.perturbedTemperature);
    print("Density Perturbation:", atmos.densityPerturbation * 100.0, atmos2.densityPerturbation * 100.0);
    print("Pressure Perturbation:", eatmos.pressurePerturbation * 100.0, eatmos2.pressurePerturbation * 100.0);
    print("Temperature Perturbation:", eatmos.temperaturePerturbation * 100.0, eatmos2.temperaturePerturbation * 100.0);
    cout << endl;

    // Get and print gases.
    print("Average Molecular Weight:", atmos.averageMolecularWeight, atmos2.averageMolecularWeight);
    print("O2 Mole Fraction:", atmos.dioxygen.moleFraction * 100.0, atmos2.dioxygen.moleFraction * 100.0);
    print("N2 Mole Fraction:", atmos.dinitrogen.moleFraction * 100.0, atmos2.dinitrogen.moleFraction * 100.0);
    print("CO2 Mole Fraction:", atmos.carbonDioxide.moleFraction * 100.0, atmos2.carbonDioxide.moleFraction * 100.0);
    print("Helium Mole Fraction:", atmos.helium.moleFraction * 100.0, atmos2.helium.moleFraction * 100.0);
    print("Argon Mole Fraction:", atmos.argon.moleFraction * 100.0, atmos2.argon.moleFraction * 100.0);
    cout << endl;

    // Print wind data
    print("EW Wind:", atmos.ewWind, atmos2.ewWind);
    print("NS Wind:", atmos.nsWind, atmos2.nsWind);
    print("Vertical Wind:", atmos.verticalWind, atmos2.verticalWind);
    print("Perturbed EW Wind:", atmos.perturbedEWWind, atmos2.perturbedEWWind);
    print("Perturbed NS Wind:", atmos.perturbedNSWind, atmos2.perturbedNSWind);
    print("Perturbed Vertical Wind:", atmos.perturbedVerticalWind, atmos2.perturbedVerticalWind);
    cout << endl;

    // Print Earth specific data
    print("Vapor Pressure:", eatmos.vaporPressure, eatmos2.vaporPressure);
    print("Vapor Density:", eatmos.vaporDensity, eatmos2.vaporDensity);
    print("Dew Point:", eatmos.dewPoint, eatmos2.dewPoint);
    print("Relative Humidity:", eatmos.relativeHumidity, eatmos2.relativeHumidity);
  }
  catch (const string& msg) {
    cerr << msg << endl;
  }
  catch (const std::runtime_error& err) {
    cerr << "An unanticipated error occurred." << endl;
    cerr << err.what() << endl;
  }

  auto t2 = Clock::now();
  cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0 << " seconds." << endl;
  return 0;
}

