from GRAMpy import gram

def jupiter_example():
    print("=============")
    print("Jupiter Example (Python)")
    print("=============")

    # Set the Spice data path (this is critical)
    inputParameters = gram.JupiterInputParameters()
    reader = gram.JupiterNamelistReader()
    reader.tryGetSpicePath(inputParameters)

    # Create a Jupiter model and set default parameters
    jupiter = gram.JupiterAtmosphere()
    jupiter.setInputParameters(inputParameters)

    # More than one JupiterAtmosphere can be created
    jupiter2 = gram.JupiterAtmosphere()
    jupiter2.setInputParameters(inputParameters)

    print(jupiter.getVersionString())

    # Set the perturbation scale factors
    jupiter.setPerturbationScales(1.5)
    jupiter.setMinRelativeStepSize(0.5)
    jupiter.setSeed(1001)
    jupiter2.setPerturbationScales(1.5)
    jupiter2.setMinRelativeStepSize(0.5)
    jupiter2.setSeed(1001)

    # Set the start time of the trajectory
    ttime = gram.GramTime()
    ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, gram.UTC, gram.ERT)
    jupiter.setStartTime(ttime)
    jupiter2.setStartTime(ttime)

    # Set the position
    position = gram.Position()
    position.height = 50.0 
    position.latitude = 22
    position.longitude = 48 
    position.elapsedTime = 100
    jupiter.setPosition(position)
    position.height = 1000
    jupiter2.setPosition(position)

    # Update the atmosphere data
    jupiter.update()
    jupiter2.update()

    pos = jupiter.getPosition()
    pos2 = jupiter2.getPosition()

    print("\n                           Jupiter #1      Jupiter #2")
    print("Total Radius: " + "{0:>22.6e}".format(pos.totalRadius) + "  " + "{0:>17.6e}".format(pos2.totalRadius))
    print("Gravity: " + "{0:>27.6e}".format(pos.gravity) + "   " + "{0:>16.6e}".format(pos2.gravity) + '\n')

    estate = jupiter.getEphemerisState()
    estate2 = jupiter2.getEphemerisState()
    print("Solar Time: " + "{0:24.6e}".format(estate.solarTime) + "  " + "{0:>17.6e}".format(estate2.solarTime)) 
    print("Longitude of the Sun: " + "{0:>14.6e}".format(estate.longitudeSun) + "  " + "{0:>17.6e}".format(estate2.longitudeSun) + '\n')

    atmos = jupiter.getAtmosphereState()
    atmos2 = jupiter2.getAtmosphereState()
    print("Temperature: " + "{0:>23.6e}".format(atmos.temperature) + "    " + "{0:>15.6e}".format(atmos2.temperature))
    print("Pressure:    " + "{0:>23.6e}".format(atmos.pressure) + "   " + "{0:>16.6e}".format(atmos2.pressure))
    print("Density: " + "{0:>27.6e}".format(atmos.density) + "    " + "{0:>15.6e}".format(atmos2.density))
    print("Pressure Scale Height:" + "{0:>14.6e}".format(atmos.pressureScaleHeight) + "   " + "{0:>16.6e}".format(atmos2.pressureScaleHeight))
    print("Density Scale Height:" + "{0:>15.6e}".format(atmos.densityScaleHeight) + " " + "{0:>18.6e}".format(atmos2.densityScaleHeight))

    print("=============")
    print("End example")
    print("=============")

if __name__ == '__main__':
    jupiter_example()
   