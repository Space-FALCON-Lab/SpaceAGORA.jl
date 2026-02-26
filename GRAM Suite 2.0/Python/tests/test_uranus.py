import unittest
from gram import *

class UranusAtmosphereTestCase(unittest.TestCase):

    def setUp(self):
        # load path to SPICE data (note: make sure spice.txt copied to appropriate locations)
        inputParameters = InputParameters()
        NamelistReader().tryGetSpicePath(inputParameters)
        SpiceLoader().setInputParameters(inputParameters)

    def test_input_parameters(self):

        uranus = UranusAtmosphere() 
        inputParameters = UranusInputParameters()
        reader = UranusNamelistReader()
        reader.tryGetSpicePath(inputParameters)
        uranus.setInputParameters(inputParameters)
        self.assertEqual(str(uranus.getInputParameters().spicePath), str(inputParameters.spicePath))
        inputParameters.minute = 33
        self.assertNotEqual(uranus.getInputParameters().minute, inputParameters.minute)
        uranus.setInputParameters(inputParameters)
        self.assertEqual(uranus.getInputParameters().minute, 33)

class UranusSystemTestCase(unittest.TestCase):
    # system test comparing outputs from Python interface to outputs from C++ interface. 
    # Reference Uranus/examples/Uranus.cpp
    u1_expected = {
        'totalRadius': 2.552427e+04, 
        'gravity': 8.694414e+00, 
        'solarTime': 2.014505e+01, 
        'longitudeSun': 4.830880e+01, 
        'temperature': 5.310000e+01, 
        'pressure': 1.260000e+04, 
        'density': 6.641600e-02, 
        'pressureScaleHeight': 2.149407e+01, 
        'densityScaleHeight': 2.225584e+01, 
        'perturbedDensity': 7.041548e-02, 
        'perturbationPercent': 6.021858e+00, 
        'averageMolecularWeight': 2.330286e+00, 
        'dihydrogenMoleFraction': 8.438978e+01, 
        'methaneMoleFraction': 1.000046e-02, 
        'heliumMoleFraction': 1.560022e+01
    }

    u2_expected = {
        'totalRadius': 2.647427e+04,
        'gravity': 8.055715e+00,
        'solarTime': 2.014505e+01,
        'longitudeSun': 4.830880e+01,
        'temperature': 5.214000e+02,
        'pressure': 1.887100e-02,
        'density': 8.864400e-09,
        'pressureScaleHeight': 2.226212e+02,
        'densityScaleHeight': 2.121697e+02,
        'perturbedDensity': 9.576135e-09,
        'perturbationPercent': 8.029144e+00,
        'averageMolecularWeight': 2.040410e+00,
        'dihydrogenMoleFraction': 9.903001e+01,
        'methaneMoleFraction': 1.000009e-02,
        'heliumMoleFraction': 9.599899e-01
    }

    def setUp(self):

        inputParameters = UranusInputParameters()
        reader = UranusNamelistReader()
        reader.tryGetSpicePath(inputParameters)

        self.uranus = UranusAtmosphere()
        self.uranus.setInputParameters(inputParameters)

        self.uranus2 = UranusAtmosphere()
        self.uranus2.setInputParameters(inputParameters)

        self.uranus.setPerturbationScales(1.5)
        self.uranus.setMinRelativeStepSize(0.5)
        self.uranus.setSeed(1001)
        self.uranus2.setPerturbationScales(1.5)
        self.uranus2.setMinRelativeStepSize(0.5)
        self.uranus2.setSeed(1001)

        ttime = GramTime()
        ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT)
        self.uranus.setStartTime(ttime)
        self.uranus2.setStartTime(ttime)

        position = Position()
        position.height = 50
        position.latitude = 22
        position.longitude = 48
        position.elapsedTime = 100
        self.uranus.setPosition(position)
        position.height = 1000
        self.uranus2.setPosition(position)

        self.uranus.update()
        self.uranus2.update()

        self.pos = self.uranus.getPosition()
        self.pos2 = self.uranus2.getPosition()

        self.estate = self.uranus.getEphemerisState()
        self.estate2 = self.uranus2.getEphemerisState()

        self.atmos = self.uranus.getAtmosphereState()
        self.atmos2 = self.uranus2.getAtmosphereState()


    # testing equality to penultimate digit of expected output. 
    # This will lead to us testing some figures that are not physically meaningful, but here
    # we are testing that the calculations done via the Python interface match with those done 
    # done directly in C++, and nothing else
    def test_position(self):
        self.assertAlmostEqual(self.pos.totalRadius, self.u1_expected['totalRadius'], places=1)
        self.assertAlmostEqual(self.pos2.totalRadius, self.u2_expected['totalRadius'], places=1)

        self.assertAlmostEqual(self.pos.gravity, self.u1_expected['gravity'], places=5)
        self.assertAlmostEqual(self.pos2.gravity, self.u2_expected['gravity'], places=5)

    def test_ephemeris_state(self):
        self.assertAlmostEqual(self.estate.solarTime, self.u1_expected['solarTime'], places=4)
        self.assertAlmostEqual(self.estate2.solarTime, self.u2_expected['solarTime'], places=4)
        self.assertAlmostEqual(self.estate.longitudeSun, self.u1_expected['longitudeSun'], places=4)
        self.assertAlmostEqual(self.estate2.longitudeSun, self.u2_expected['longitudeSun'], places=4)
        
    def test_atmosphere_state(self):
        self.assertAlmostEqual(self.atmos.temperature, self.u1_expected['temperature'], places=4)
        self.assertAlmostEqual(self.atmos2.temperature, self.u2_expected['temperature'], places=3)
        self.assertAlmostEqual(self.atmos.pressure, self.u1_expected['pressure'], places=1)
        self.assertAlmostEqual(self.atmos2.pressure, self.u2_expected['pressure'], places=7)
        self.assertAlmostEqual(self.atmos.density, self.u1_expected['density'], places=7)
        self.assertAlmostEqual(self.atmos2.density, self.u2_expected['density'], places=14)
        self.assertAlmostEqual(self.atmos.pressureScaleHeight, self.u1_expected['pressureScaleHeight'], places=4)
        self.assertAlmostEqual(self.atmos2.pressureScaleHeight, self.u2_expected['pressureScaleHeight'], places=3)
        self.assertAlmostEqual(self.atmos.densityScaleHeight, self.u1_expected['densityScaleHeight'], places=4)
        self.assertAlmostEqual(self.atmos2.densityScaleHeight, self.u2_expected['densityScaleHeight'], places=3)
        self.assertAlmostEqual(self.atmos.perturbedDensity, self.u1_expected['perturbedDensity'], places=7)
        self.assertAlmostEqual(self.atmos2.perturbedDensity, self.u2_expected['perturbedDensity'], places=14)
        self.assertAlmostEqual(self.atmos.densityPerturbation * 100, self.u1_expected['perturbationPercent'], places=5)
        self.assertAlmostEqual(self.atmos2.densityPerturbation * 100, self.u2_expected['perturbationPercent'], places = 5)
    
    def test_gasses(self):
        self.assertAlmostEqual(self.atmos.averageMolecularWeight, self.u1_expected['averageMolecularWeight'], places = 5)
        self.assertAlmostEqual(self.atmos2.averageMolecularWeight, self.u2_expected['averageMolecularWeight'], places = 5)
        self.assertAlmostEqual(self.atmos.dihydrogen.moleFraction * 100, self.u1_expected['dihydrogenMoleFraction'], places=4)
        self.assertAlmostEqual(self.atmos2.dihydrogen.moleFraction * 100, self.u2_expected['dihydrogenMoleFraction'], places=4)
        self.assertAlmostEqual(self.atmos.methane.moleFraction * 100, self.u1_expected['methaneMoleFraction'], places=7)
        self.assertAlmostEqual(self.atmos2.methane.moleFraction * 100, self.u2_expected['methaneMoleFraction'], places=7)
        self.assertAlmostEqual(self.atmos.helium.moleFraction * 100, self.u1_expected['heliumMoleFraction'], places=4)
        self.assertAlmostEqual(self.atmos2.helium.moleFraction * 100, self.u2_expected['heliumMoleFraction'], places=6)


if __name__ == '__main__':
    unittest.main()