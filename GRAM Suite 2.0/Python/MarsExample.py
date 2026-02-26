from GRAMpy import gram

MARS_DATA_PATH = '../Mars/data/'

def mars_example():
    
    print("=============")
    print("Mars example (Python)")
    print("=============")

    # Set the Spice data path (this is critical)
    inputParameters = gram.MarsInputParameters()
    inputParameters.dataPath = MARS_DATA_PATH
    reader = gram.MarsNamelistReader()
    reader.tryGetSpicePath(inputParameters)

    # Create a Mars model and set default parameters
    mars = gram.MarsAtmosphere()
    mars.setInputParameters(inputParameters)
    print(mars)
    # More than one MarsAtmosphere can be created
    mars2 = gram.MarsAtmosphere()
    mars2.setInputParameters(inputParameters)

    print(mars.getVersionString())

    # Set a Mars specific parameter
    mars.setMapYear(0)
    mars.setMGCMDustLevels(2, 0, 0)
    mars2.setMapYear(2)

    # Set the perturbation scale factors
    mars.setPerturbationScales(1.5)
    mars.setMinRelativeStepSize(0.5)
    mars.setSeed(1001)
    mars2.setPerturbationScales(1.5)
    mars2.setMinRelativeStepSize(0.5)
    mars2.setSeed(1001)

    # Set the start time of the trajectory
    ttime = gram.GramTime()
    ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, gram.UTC, gram.ERT)
    mars.setStartTime(ttime)
    mars2.setStartTime(ttime)

    # Set the position
    position = gram.Position()
    position.height = 2.2
    position.latitude = 22.345
    position.longitude = 48.5675
    position.elapsedTime = 100.13245
    mars.setPosition(position)
    mars2.setPosition(position)

    # Update the atmosphere data
    mars.update()
    mars2.update()

    pos = mars.getPosition()
    pos2 = mars2.getPosition()

    print("\n                                   Mars #1            Mars #2")
    print("Total Radius: " + "{0:>30.6e}".format(pos.totalRadius) + "  " + "{0:>17.6e}".format(pos2.totalRadius))
    print("Gravity: " + "{0:>35.6e}".format(pos.gravity) + "   " + "{0:>16.6e}".format(pos2.gravity) + '\n')

    estate = mars.getEphemerisState()
    estate2 = mars2.getEphemerisState()
    print("Solar Time: " + "{0:32.6e}".format(estate.solarTime) + "  " + "{0:>17.6e}".format(estate2.solarTime)) 
    print("Longitude of the Sun: " + "{0:>22.6e}".format(estate.longitudeSun) + "  " + "{0:>17.6e}".format(estate2.longitudeSun) + '\n')

    atmos = mars.getAtmosphereState()
    atmos2 = mars2.getAtmosphereState()
    print("Temperature: " + "{0:>31.6e}".format(atmos.temperature) + "    " + "{0:>15.6e}".format(atmos2.temperature))
    print("Pressure:    " + "{0:>31.6e}".format(atmos.pressure) + "   " + "{0:>16.6e}".format(atmos2.pressure))
    print("Density: " + "{0:>35.6e}".format(atmos.density) + "    " + "{0:>15.6e}".format(atmos2.density))
    print("Pressure Scale Height:" + "{0:>22.6e}".format(atmos.pressureScaleHeight) + "   " + "{0:>16.6e}".format(atmos2.pressureScaleHeight))
    print("Density Scale Height:" + "{0:>23.6e}".format(atmos.densityScaleHeight) + " " + "{0:>18.6e}".format(atmos2.densityScaleHeight) + '\n')

    # Print perturbed density
    print("Mean Density: " + "{0:>30.6e}".format(atmos.density) + "    " + "{0:>15.6e}".format(atmos2.density))
    print("Perturbed Density: " + "{0:>25.6e}".format(atmos.perturbedDensity) + "    " + "{0:>15.6e}".format(atmos2.perturbedDensity))
    print("Perturbation Percent:  " + "{0:>21.6e}".format(atmos.densityPerturbation * 100) + "    " + "{0:>15.6e}".format(atmos2.densityPerturbation * 100) + '\n')


    # Get and print gases
    print("Average Molecular Weight:  " + "{0:>17.6e}".format(atmos.averageMolecularWeight) + "    " + "{0:>15.6e}".format(atmos2.averageMolecularWeight))
    print("Carbon Dioxide Mole Fraction:  " + "{0:>13.6e}".format(atmos.carbonDioxide.moleFraction * 100) + "    " + "{0:>15.6e}".format(atmos2.carbonDioxide.moleFraction * 100))
    print("Dinitrogen Mole Fraction:  " + "{0:>17.6e}".format(atmos.dinitrogen.moleFraction * 100) + "    " + "{0:>15.6e}".format(atmos2.dinitrogen.moleFraction * 100))
    print("Argon Mole Fraction:  " + "{0:>22.6e}".format(atmos.argon.moleFraction * 100) + "    " + "{0:>15.6e}".format(atmos2.argon.moleFraction * 100))
    print("Carbon Monoxide Mole Fraction:  " + "{0:>12.6e}".format(atmos.carbonMonoxide.moleFraction * 100) + "    " + "{0:>15.6e}".format(atmos2.carbonMonoxide.moleFraction * 100))
    print("Dioxygen Mole Fraction:  " + "{0:>19.6e}".format(atmos.dioxygen.moleFraction * 100) + "    " + "{0:>15.6e}".format(atmos2.dioxygen.moleFraction * 100))
    print("Dihydrogen Mole Fraction:  " + "{0:>17.6e}".format(atmos.dihydrogen.moleFraction * 100) + "    " + "{0:>15.6e}".format(atmos2.dihydrogen.moleFraction * 100))
    print("Hydrogen Mole Fraction:  " + "{0:>19.6e}".format(atmos.hydrogen.moleFraction * 100) + "    " + "{0:>15.6e}".format(atmos2.hydrogen.moleFraction * 100))
    print("Oxygen Mole Fraction:  " + "{0:>21.6e}".format(atmos.oxygen.moleFraction * 100) + "    " + "{0:>15.6e}".format(atmos2.oxygen.moleFraction * 100))
    print("Helium Mole Fraction:  " + "{0:>21.6e}".format(atmos.helium.moleFraction * 100) + "    " + "{0:>15.6e}".format(atmos2.helium.moleFraction * 100))
    print("Water Vapor Mole Fraction:  " + "{0:>16.6e}".format(atmos.water.moleFraction * 100) + "    " + "{0:>15.6e}".format(atmos2.water.moleFraction * 100) + '\n')

    # Print wind data
    print("EW Wind: " + "{0:>35.6e}".format(atmos.ewWind) + "   " + "{0:>16.6e}".format(atmos2.ewWind))
    print("NS Wind: " + "{0:>35.6e}".format(atmos.nsWind) + "   " + "{0:>16.6e}".format(atmos2.nsWind))
    print("Vertical Wind: " + "{0:>29.6e}".format(atmos.verticalWind) + "   " + "{0:>16.6e}".format(atmos2.verticalWind))
    print("Perturbed EW Wind: " + "{0:>25.6e}".format(atmos.perturbedEWWind) + "   " + "{0:>16.6e}".format(atmos2.perturbedEWWind))
    print("Perturbed NS Wind: " + "{0:>25.6e}".format(atmos.perturbedNSWind) + "   " + "{0:>16.6e}".format(atmos2.perturbedNSWind) + '\n')

    # Print Mars specific metrics
    marsAtmos = mars.getMarsAtmosphereState()
    marsAtmos2 = mars2.getMarsAtmosphereState()
    print("Dust Optical Depth:" + "{0:>25.6e}".format(marsAtmos.dustOpticalDepth) + "{0:>19.6e}".format(marsAtmos2.dustOpticalDepth))
    print("Daily Mean Temperature:" + "{0:>21.6e}".format(marsAtmos.temperatureDaily) + "{0:>19.6e}".format(marsAtmos2.temperatureDaily))
    print("Daily Min Temperature:" + "{0:>22.6e}".format(marsAtmos.temperatureMin) + "{0:>19.6e}".format(marsAtmos2.temperatureMin))
    print("Daily Max Temperature:" + "{0:>22.6e}".format(marsAtmos.temperatureMax) + "{0:>19.6e}".format(marsAtmos2.temperatureMax))

    print("=============")
    print("End example")
    print("=============")

if __name__ == '__main__':
    mars_example()