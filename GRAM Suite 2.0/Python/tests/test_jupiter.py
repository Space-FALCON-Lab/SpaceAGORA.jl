import unittest
from gram import *

class JupiterAtmosphereTestCase(unittest.TestCase):

    def setUp(self):
        # load path to SPICE data (note: make sure spice.txt copied to appropriate locations)
        inputParameters = InputParameters()
        NamelistReader().tryGetSpicePath(inputParameters)
        SpiceLoader().setInputParameters(inputParameters)

    def test_input_parameters(self):

        jupiter = JupiterAtmosphere() 
        inputParameters = JupiterInputParameters()
        reader = JupiterNamelistReader()
        reader.tryGetSpicePath(inputParameters)
        jupiter.setInputParameters(inputParameters)
        self.assertEqual(str(jupiter.getInputParameters().spicePath), str(inputParameters.spicePath))
        inputParameters.minute = 33
        self.assertNotEqual(jupiter.getInputParameters().minute, inputParameters.minute)
        jupiter.setInputParameters(inputParameters)
        self.assertEqual(jupiter.getInputParameters().minute, 33)

class JupiterSystemTestCase(unittest.TestCase):
    # system test comparing outputs from Python interface to outputs from C++ interface
    # Reference Jupiter/examples/Jupiter.cpp

    j1_expected = {
            'totalRadius': 7.083256e+04, 
            'gravity': 2.369636e+01,
            'solarTime': 1.018768e+01, 
            'longitudeSun': 3.246603e+02, 
            'temperature': 1.138488e+02,
            'pressure': 7.695541e+03, 
            'density': 1.878207e-02, 
            'pressureScaleHeight': 1.763863e+01,
            'densityScaleHeight': 7.697741e+01
        }

    j2_expected = {
            'totalRadius': 7.178256e+04, 
            'gravity': 2.299010e+01,
            'solarTime': 1.018768e+01, 
            'longitudeSun': 3.246603e+02, 
            'temperature': 8.969180e+02,
            'pressure': 1.147822e-04, 
            'density': 3.081254e-11, 
            'pressureScaleHeight': 1.615685e+02,
            'densityScaleHeight': 1.485278e+02
        }

    def setUp(self):

        inputParameters = JupiterInputParameters()
        reader = JupiterNamelistReader()
        reader.tryGetSpicePath(inputParameters)

        self.jupiter = JupiterAtmosphere()
        self.jupiter.setInputParameters(inputParameters)

        self.jupiter2 = JupiterAtmosphere()
        self.jupiter2.setInputParameters(inputParameters)

        self.jupiter.setPerturbationScales(1.5)
        self.jupiter.setMinRelativeStepSize(0.5)
        self.jupiter.setSeed(1001)

        self.jupiter2.setPerturbationScales(1.5)
        self.jupiter2.setMinRelativeStepSize(0.5)
        self.jupiter2.setSeed(1001)

        ttime = GramTime()
        ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT)
        self.jupiter.setStartTime(ttime)
        self.jupiter2.setStartTime(ttime)

        position = Position()
        position.height = 50.0
        position.latitude = 22.0
        position.longitude = 48.0
        position.elapsedTime = 100.0
        self.jupiter.setPosition(position)
        position.height = 1000.0
        self.jupiter2.setPosition(position)

        self.jupiter.update()
        self.jupiter2.update()

        self.pos = self.jupiter.getPosition()
        self.pos2 = self.jupiter2.getPosition()

        self.estate = self.jupiter.getEphemerisState()
        self.estate2 = self.jupiter2.getEphemerisState()

        self.atmos = self.jupiter.getAtmosphereState()
        self.atmos2 = self.jupiter2.getAtmosphereState()

    def test_position(self):
        self.assertAlmostEqual(self.pos.totalRadius, self.j1_expected['totalRadius'], places=1)
        self.assertAlmostEqual(self.pos2.totalRadius, self.j2_expected['totalRadius'], places=1)

        self.assertAlmostEqual(self.pos.gravity, self.j1_expected['gravity'], places=5)
        self.assertAlmostEqual(self.pos2.gravity, self.j2_expected['gravity'], places=5)

    def test_ephemeris_state(self):
        self.assertAlmostEqual(self.estate.solarTime, self.j1_expected['solarTime'], places=5)
        self.assertAlmostEqual(self.estate2.solarTime, self.j2_expected['solarTime'], places=5)
        self.assertAlmostEqual(self.estate.longitudeSun, self.j1_expected['longitudeSun'], places=4)
        self.assertAlmostEqual(self.estate2.longitudeSun, self.j2_expected['longitudeSun'], places=4)
        
    def test_atmosphere_state(self):
        self.assertAlmostEqual(self.atmos.temperature, self.j1_expected['temperature'], places=4)
        self.assertAlmostEqual(self.atmos2.temperature, self.j2_expected['temperature'], places=4)
        self.assertAlmostEqual(self.atmos.pressure, self.j1_expected['pressure'], places=2)
        self.assertAlmostEqual(self.atmos2.pressure, self.j2_expected['pressure'], places=9)
        self.assertAlmostEqual(self.atmos.density, self.j1_expected['density'], places=7)
        self.assertAlmostEqual(self.atmos2.density, self.j2_expected['density'], places=16)
        self.assertAlmostEqual(self.atmos.pressureScaleHeight, self.j1_expected['pressureScaleHeight'], places=5)
        self.assertAlmostEqual(self.atmos2.pressureScaleHeight, self.j2_expected['pressureScaleHeight'], places=4)
        self.assertAlmostEqual(self.atmos.densityScaleHeight, self.j1_expected['densityScaleHeight'], places=5)
        self.assertAlmostEqual(self.atmos2.densityScaleHeight, self.j2_expected['densityScaleHeight'], places=4)


if __name__ == '__main__':
    unittest.main()