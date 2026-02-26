import unittest
from gram import *

EARTH_DATA_PATH = '../../Earth/data/'

class EarthAtmosphereAtmosphereTestCase(unittest.TestCase):

    # def setUp(self):
    #     # load path to SPICE data (note: make sure spice.txt copied to appropriate locations)
    #     inputParameters = EarthInputParameters()
    #     EarthNamelistReader().tryGetSpicePath(inputParameters)
    #     SpiceLoader().setInputParameters(inputParameters)

    def test_input_parameters(self):

        earth = EarthAtmosphere() 
        inputParameters = EarthInputParameters()
        inputParameters.dataPath = EARTH_DATA_PATH
        reader = EarthNamelistReader()
        reader.tryGetSpicePath(inputParameters)
        earth.setInputParameters(inputParameters)

        self.assertEqual(str(earth.getInputParameters().spicePath), str(inputParameters.spicePath))
        inputParameters.minute = 33
        self.assertNotEqual(earth.getInputParameters().minute, inputParameters.minute)
        earth.setInputParameters(inputParameters)
        self.assertEqual(earth.getInputParameters().minute, 33)

    
class EarthSystemTestCase(unittest.TestCase):
    # system test comparing inputs from Python interface to outputs from C++ interface
    # reference Earth/examples/Earth.cpp

    e1_expected = {
        'totalRadius': 6.380123e+03, 
        'gravity': 9.772222e+00, 
        'solarTime': 3.304106e+00, 
        'longitudeSun': 2.634718e+02, 
        'temperature': 2.694625e+02, 
        'pressure': 5.570779e+04, 
        'density': 7.198355e-01, 
        'pressureScaleHeight': 7.919346e+00, 
        'densityScaleHeight': 9.641089e+00, 
        'perturbedDensity': 7.134327e-01, 
        'perturbedPressure': 5.548469e+04, 
        'perturbedTemperature': 2.707920e+02, 
        'densityPerturbation': -8.894850e-01, 
        'pressurePerturbation': -4.004782e-01, 
        'temperaturePerturbation': 4.933954e-01, 
        'averageMolecularWeight': 2.894918e+01, 
        'dioxygenMoleFraction': 2.091346e+01, 
        'dinitrogenMoleFraction': 7.795798e+01, 
        'carbonDioxideMoleFraction': 4.103220e-02, 
        'heliumMoleFraction': 5.191674e-04, 
        'argonMoleFraction': 9.325046e-01, 
        'EWWind': 1.027716e+01, 
        'NSWind': 2.355525e-01, 
        'verticalWind': 3.299024e-02, 
        'perturbedEWWind': 7.652076e+00, 
        'perturbedNSWind': -1.364762e-01, 
        'perturbedVerticalWind':-1.008657e+00, 
        'vaporPressure': 8.576418e+01, 
        'vaporDensity': 6.896226e-04, 
        'dewPoint': 2.457944e+02, 
        'relativeHumidity': 1.969167e+01
    }
    
    e2_expected = {
        'totalRadius': 6.380123e+03,
        'gravity': 9.772222e+00,
        'solarTime': 3.304106e+00,
        'longitudeSun': 2.634718e+02,
        'temperature': 2.695653e+02,
        'pressure': 5.574820e+04,
        'density': 7.200775e-01,
        'pressureScaleHeight': 7.922428e+00,
        'densityScaleHeight': 9.665737e+00,
        'perturbedDensity': 7.137008e-01,
        'perturbedPressure': 5.553302e+04,
        'perturbedTemperature': 2.709240e+02,
        'densityPerturbation': -8.855612e-01,
        'pressurePerturbation': -3.859920e-01,
        'temperaturePerturbation': 5.040327e-01,
        'averageMolecularWeight': 2.895028e+01,
        'dioxygenMoleFraction': 2.091557e+01,
        'dinitrogenMoleFraction': 7.796584e+01,
        'carbonDioxideMoleFraction': 4.103633e-02,
        'heliumMoleFraction': 5.192198e-04,
        'argonMoleFraction': 9.325986e-01,
        'EWWind': 1.029080e+01,
        'NSWind': 4.632215e-01,
        'verticalWind': 1.381595e-02,
        'perturbedEWWind': 7.711405e+00,
        'perturbedNSWind': 1.060719e-01,
        'perturbedVerticalWind':-1.027831e+00,
        'vaporPressure': 8.022889e+01,
        'vaporDensity': 6.448678e-04,
        'dewPoint': 2.451434e+02,
        'relativeHumidity': 1.811133e+01
    }

    def setUp(self):
        inputParameters = EarthInputParameters()
        inputParameters.dataPath = EARTH_DATA_PATH
        reader = EarthNamelistReader()
        reader.tryGetSpicePath(inputParameters)

        self.earth = EarthAtmosphere()
        self.earth.setInputParameters(inputParameters)

        self.earth2 = EarthAtmosphere()
        self.earth2.setInputParameters(inputParameters)

        self.earth.setNCEPParameters(9715, 1)
        self.earth2.setNCEPParameters(9715, 5)

        self.earth.setRandomPerturbationScale(1.5)
        self.earth.setSeed(1001)

        self.earth2.setRandomPerturbationScale(1.5)
        self.earth2.setSeed(1001)

        ttime =GramTime()
        ttime.setStartTime(2020, 12, 15, 0, 0, 0.0, UTC, PET)
        self.earth.setStartTime(ttime)
        self.earth2.setStartTime(ttime)

        position = Position()
        position.height = 5
        position.latitude = 22
        position.longitude = 48
        position.elapsedTime = 100
        self.earth.setPosition(position)
        self.earth2.setPosition(position)

        self.earth.update()
        self.earth2.update()

        self.pos = self.earth.getPosition()
        self.pos2 = self.earth2.getPosition()

        self.estate = self.earth.getEphemerisState()
        self.estate2 = self.earth2.getEphemerisState()

        self.atmos = self.earth.getAtmosphereState()
        self.atmos2 = self.earth2.getAtmosphereState()

        self.eatmos = self.earth.getEarthAtmosphereState()
        self.eatmos2 = self.earth2.getEarthAtmosphereState()

    # testing equality to penultimate digit of expected output. 
    # This will lead to us testing some figures that are not physically meaningful, but here
    # we are testing that the calculations done via the Python interface match with those done 
    # done directly in C++, and nothing else
    def test_position(self):
        self.assertAlmostEqual(self.pos.totalRadius, self.e1_expected['totalRadius'], places=2)
        self.assertAlmostEqual(self.pos2.totalRadius, self.e2_expected['totalRadius'], places=2)
        self.assertAlmostEqual(self.pos.gravity, self.e1_expected['gravity'], places=5)
        self.assertAlmostEqual(self.pos2.gravity, self.e2_expected['gravity'], places=5)

    def test_ephemeris_state(self):
        self.assertAlmostEqual(self.estate.solarTime, self.e1_expected['solarTime'], places=5)
        self.assertAlmostEqual(self.estate2.solarTime, self.e2_expected['solarTime'], places=5)
        self.assertAlmostEqual(self.estate.longitudeSun, self.e1_expected['longitudeSun'], places=3)
        self.assertAlmostEqual(self.estate2.longitudeSun, self.e2_expected['longitudeSun'], places=3)
        
    def test_atmosphere_state(self):
        self.assertAlmostEqual(self.atmos.temperature, self.e1_expected['temperature'], places=3)
        self.assertAlmostEqual(self.atmos2.temperature, self.e2_expected['temperature'], places=3)
        self.assertAlmostEqual(self.atmos.pressure, self.e1_expected['pressure'], places=1)
        self.assertAlmostEqual(self.atmos2.pressure, self.e2_expected['pressure'], places=1)
        self.assertAlmostEqual(self.atmos.density, self.e1_expected['density'], places=6)
        self.assertAlmostEqual(self.atmos2.density, self.e2_expected['density'], places=6)
        self.assertAlmostEqual(self.atmos.pressureScaleHeight, self.e1_expected['pressureScaleHeight'], places=5)
        self.assertAlmostEqual(self.atmos2.pressureScaleHeight, self.e2_expected['pressureScaleHeight'], places=5)
        self.assertAlmostEqual(self.atmos.densityScaleHeight, self.e1_expected['densityScaleHeight'], places=5)
        self.assertAlmostEqual(self.atmos2.densityScaleHeight, self.e2_expected['densityScaleHeight'], places=5)
        self.assertAlmostEqual(self.atmos.perturbedDensity, self.e1_expected['perturbedDensity'], places=6)
        self.assertAlmostEqual(self.atmos2.perturbedDensity, self.e2_expected['perturbedDensity'], places=6)
        self.assertAlmostEqual(self.atmos.densityPerturbation * 100, self.e1_expected['densityPerturbation'], places=6)
        self.assertAlmostEqual(self.atmos2.densityPerturbation * 100, self.e2_expected['densityPerturbation'], places = 6)


    def test_earth_atmosphere_state(self):
        self.assertAlmostEqual(self.eatmos.perturbedPressure, self.e1_expected['perturbedPressure'], places=1)
        self.assertAlmostEqual(self.eatmos2.perturbedPressure, self.e2_expected['perturbedPressure'], places=1)
        self.assertAlmostEqual(self.eatmos.perturbedTemperature, self.e1_expected['perturbedTemperature'], places=3)
        self.assertAlmostEqual(self.eatmos2.perturbedTemperature, self.e2_expected['perturbedTemperature'], places=3)
        self.assertAlmostEqual(self.atmos.densityPerturbation * 100, self.e1_expected['densityPerturbation'], places=6)
        self.assertAlmostEqual(self.atmos2.densityPerturbation * 100, self.e2_expected['densityPerturbation'], places=6)
        self.assertAlmostEqual(self.eatmos.pressurePerturbation * 100, self.e1_expected['pressurePerturbation'], places=6)
        self.assertAlmostEqual(self.eatmos2.pressurePerturbation * 100, self.e2_expected['pressurePerturbation'], places=6)
        self.assertAlmostEqual(self.eatmos.temperaturePerturbation * 100, self.e1_expected['temperaturePerturbation'], places=6)
        self.assertAlmostEqual(self.eatmos2.temperaturePerturbation * 100, self.e2_expected['temperaturePerturbation'], places=6)
        self.assertAlmostEqual(self.eatmos.vaporPressure, self.e1_expected['vaporPressure'], places=4)
        self.assertAlmostEqual(self.eatmos2.vaporPressure, self.e2_expected['vaporPressure'], places=4)
        self.assertAlmostEqual(self.eatmos.vaporDensity, self.e1_expected['vaporDensity'], places=9)
        self.assertAlmostEqual(self.eatmos2.vaporDensity, self.e2_expected['vaporDensity'], places=9)
        self.assertAlmostEqual(self.eatmos.dewPoint, self.e1_expected['dewPoint'], places=3)
        self.assertAlmostEqual(self.eatmos2.dewPoint, self.e2_expected['dewPoint'], places=3)
        self.assertAlmostEqual(self.eatmos.relativeHumidity, self.e1_expected['relativeHumidity'], places=4)
        self.assertAlmostEqual(self.eatmos2.relativeHumidity, self.e2_expected['relativeHumidity'], places=4)


    def test_gasses(self):
        self.assertAlmostEqual(self.atmos.averageMolecularWeight, self.e1_expected['averageMolecularWeight'], places=5)
        self.assertAlmostEqual(self.atmos2.averageMolecularWeight, self.e2_expected['averageMolecularWeight'], places=5)
        self.assertAlmostEqual(self.atmos.helium.moleFraction * 100, self.e1_expected['heliumMoleFraction'], places=9)
        self.assertAlmostEqual(self.atmos2.helium.moleFraction * 100, self.e2_expected['heliumMoleFraction'], places=9)
        self.assertAlmostEqual(self.atmos.dioxygen.moleFraction * 100, self.e1_expected['dioxygenMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos2.dioxygen.moleFraction * 100, self.e2_expected['dioxygenMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos.dinitrogen.moleFraction * 100, self.e1_expected['dinitrogenMoleFraction'], places=4)
        self.assertAlmostEqual(self.atmos2.dinitrogen.moleFraction * 100, self.e2_expected['dinitrogenMoleFraction'], places=4)
        self.assertAlmostEqual(self.atmos.carbonDioxide.moleFraction * 100, self.e1_expected['carbonDioxideMoleFraction'], places=7)
        self.assertAlmostEqual(self.atmos2.carbonDioxide.moleFraction * 100, self.e2_expected['carbonDioxideMoleFraction'], places=7)
        self.assertAlmostEqual(self.atmos.argon.moleFraction * 100.0, self.e1_expected['argonMoleFraction'], places=6)
        self.assertAlmostEqual(self.atmos2.argon.moleFraction * 100.0, self.e2_expected['argonMoleFraction'], places=6)
        
    
    def test_winds(self):
        self.assertAlmostEqual(self.atmos.ewWind, self.e1_expected['EWWind'], places=4)
        self.assertAlmostEqual(self.atmos2.ewWind, self.e2_expected['EWWind'], places=4)
        self.assertAlmostEqual(self.atmos.nsWind, self.e1_expected['NSWind'], places=6)
        self.assertAlmostEqual(self.atmos2.nsWind, self.e2_expected['NSWind'], places=6)
        self.assertAlmostEqual(self.atmos.verticalWind, self.e1_expected['verticalWind'], places=7)
        self.assertAlmostEqual(self.atmos2.verticalWind, self.e2_expected['verticalWind'], places=7)
        self.assertAlmostEqual(self.atmos.perturbedEWWind, self.e1_expected['perturbedEWWind'], places=5)
        self.assertAlmostEqual(self.atmos2.perturbedEWWind, self.e2_expected['perturbedEWWind'], places=5)
        self.assertAlmostEqual(self.atmos.perturbedNSWind, self.e1_expected['perturbedNSWind'], places=6)
        self.assertAlmostEqual(self.atmos2.perturbedNSWind, self.e2_expected['perturbedNSWind'], places=6)
        self.assertAlmostEqual(self.atmos.perturbedVerticalWind, self.e1_expected['perturbedVerticalWind'], places=5)
        self.assertAlmostEqual(self.atmos2.perturbedVerticalWind, self.e2_expected['perturbedVerticalWind'], places=5)

if __name__ == '__main__':
    unittest.main()    
