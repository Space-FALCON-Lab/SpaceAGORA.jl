import unittest
from gram import *

class VenusAtmosphereTestCase(unittest.TestCase):

    def setUp(self):
        # load path to SPICE data (note: make sure spice.txt copied to appropriate locations)
        inputParameters = InputParameters()
        NamelistReader().tryGetSpicePath(inputParameters)
        SpiceLoader().setInputParameters(inputParameters)

    def test_input_parameters(self):

        venus = VenusAtmosphere() 
        inputParameters = VenusInputParameters()
        reader = VenusNamelistReader()
        reader.tryGetSpicePath(inputParameters)
        venus.setInputParameters(inputParameters)
        self.assertEqual(str(venus.getInputParameters().spicePath), str(inputParameters.spicePath))
        inputParameters.minute = 33
        self.assertNotEqual(venus.getInputParameters().minute, inputParameters.minute)
        venus.setInputParameters(inputParameters)
        self.assertEqual(venus.getInputParameters().minute, 33)

class VenusSystemTestCase(unittest.TestCase):
    # system test comparing outputs from Python interface to outputs from C++ interface
    # Reference Venus/examples/Venus.cpp
    v1_expected = {
        'totalRadius': 6.101800e+03, 
        'gravity': 8.725294e+00,
        'solarTime': 2.848482e+00, 
        'longitudeSun': 2.454489e+02, 
        'temperature': 3.505000e+02, 
        'pressure': 1.066138e+05, 
        'density': 1.594000e+00, 
        'pressureScaleHeight': 7.747573e+00, 
        'densityScaleHeight': 9.396161e+00, 
        'meanDensity': 1.594000e+00, 
        'perturbedDensity': 1.612961e+00, 
        'perturbationPercent': 1.189503e+00, 
        'averageMolecularWeight': 4.344024e+01, 
        'hydrogenMoleFraction': 0.000000e+00, 
        'heliumMoleFraction': 0.000000e+00, 
        'oxygenMoleFraction': 0.000000e+00, 
        'nitrogenMoleFraction': 0.000000e+00, 
        'dinitrogenMoleFraction': 3.500000e+00, 
        'carbonMonoxideMoleFraction': 0.000000e+00, 
        'carbonDioxideMoleFraction': 9.650000e+01, 
        'EWWind': -6.265051e+01, 
        'NSWind': 1.671433e+00, 
        'perturbedEWWind': -7.510286e+01, 
        'perturbedNSWind': 5.188443e-01
    }


    v2_expected = { 
        'totalRadius': 6.251800e+03,
        'gravity': 8.311621e+00,
        'solarTime': 2.848482e+00,
        'longitudeSun': 2.454489e+02,
        'temperature': 1.311531e+02,
        'pressure': 9.390470e-08,
        'density': 2.363795e-13,
        'pressureScaleHeight': 4.433449e+01,
        'densityScaleHeight': 1.842438e+01,
        'meanDensity': 2.363795e-13,
        'perturbedDensity': 2.574676e-13,
        'perturbationPercent': 8.921271e+00,
        'averageMolecularWeight': 2.744955e+00,
        'hydrogenMoleFraction': 5.962104e+01,
        'heliumMoleFraction': 3.595756e+01,
        'oxygenMoleFraction': 4.352081e+00,
        'nitrogenMoleFraction': 5.957803e-02,
        'dinitrogenMoleFraction': 5.489343e-03,
        'carbonMonoxideMoleFraction': 4.252110e-03,
        'carbonDioxideMoleFraction': 4.698771e-06,
        'EWWind': 1.602325e+02,
        'NSWind': -6.521839e+01,
        'perturbedEWWind': 1.301058e+02,
        'perturbedNSWind': -6.800691e+01
    }


    def setUp(self):

        inputParameters = VenusInputParameters()
        reader = VenusNamelistReader()
        reader.tryGetSpicePath(inputParameters)

        self.venus = VenusAtmosphere()
        self.venus.setInputParameters(inputParameters)

        self.venus2 = VenusAtmosphere()
        self.venus2.setInputParameters(inputParameters)

        self.venus.setPerturbationScales(1.5)
        self.venus.setMinRelativeStepSize(0.5)
        self.venus.setSeed(1001)
        self.venus2.setPerturbationScales(1.5)
        self.venus2.setMinRelativeStepSize(0.5)
        self.venus2.setSeed(1001)

        ttime = GramTime()
        ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT)
        self.venus.setStartTime(ttime)
        self.venus2.setStartTime(ttime)

        position = Position()
        position.height = 50
        position.latitude = 22
        position.longitude = 48
        position.elapsedTime = 100
        self.venus.setPosition(position)
        position.height = 200
        self.venus2.setPosition(position)

        self.venus.update()
        self.venus2.update()

        self.pos = self.venus.getPosition()
        self.pos2 = self.venus2.getPosition()

        self.estate = self.venus.getEphemerisState()
        self.estate2 = self.venus2.getEphemerisState()

        self.atmos = self.venus.getAtmosphereState()
        self.atmos2 = self.venus2.getAtmosphereState()

    # testing equality to penultimate digit of expected output. 
    # This will lead to us testing some figures that are not physically meaningful, but here
    # we are testing that the calculations done via the Python interface match with those done 
    # done directly in C++, and nothing else
    def test_position(self):
        self.assertAlmostEqual(self.pos.totalRadius, self.v1_expected['totalRadius'], places=2)
        self.assertAlmostEqual(self.pos2.totalRadius, self.v2_expected['totalRadius'], places=2)
        self.assertAlmostEqual(self.pos.gravity, self.v1_expected['gravity'], places=5)
        self.assertAlmostEqual(self.pos2.gravity, self.v2_expected['gravity'], places=5)

    def test_ephemeris_state(self):
        self.assertAlmostEqual(self.estate.solarTime, self.v1_expected['solarTime'], places=5)
        self.assertAlmostEqual(self.estate2.solarTime, self.v2_expected['solarTime'], places=5)
        self.assertAlmostEqual(self.estate.longitudeSun, self.v1_expected['longitudeSun'], places=3)
        self.assertAlmostEqual(self.estate2.longitudeSun, self.v2_expected['longitudeSun'], places=3)
        
    def test_atmosphere_state(self):
        self.assertAlmostEqual(self.atmos.temperature, self.v1_expected['temperature'], places=3)
        self.assertAlmostEqual(self.atmos2.temperature, self.v2_expected['temperature'], places=3)
        self.assertAlmostEqual(self.atmos.pressure, self.v1_expected['pressure'], places=0)
        self.assertAlmostEqual(self.atmos2.pressure, self.v2_expected['pressure'], places=13)
        self.assertAlmostEqual(self.atmos.density, self.v1_expected['density'], places=5)
        self.assertAlmostEqual(self.atmos2.density, self.v2_expected['density'], places=18)
        self.assertAlmostEqual(self.atmos.pressureScaleHeight, self.v1_expected['pressureScaleHeight'], places=5)
        self.assertAlmostEqual(self.atmos2.pressureScaleHeight, self.v2_expected['pressureScaleHeight'], places=4)
        self.assertAlmostEqual(self.atmos.densityScaleHeight, self.v1_expected['densityScaleHeight'], places=5)
        self.assertAlmostEqual(self.atmos2.densityScaleHeight, self.v2_expected['densityScaleHeight'], places=4)
        self.assertAlmostEqual(self.atmos.perturbedDensity, self.v1_expected['perturbedDensity'], places=5)
        self.assertAlmostEqual(self.atmos2.perturbedDensity, self.v2_expected['perturbedDensity'], places=18)
        self.assertAlmostEqual(self.atmos.densityPerturbation * 100, self.v1_expected['perturbationPercent'], places=5)
        self.assertAlmostEqual(self.atmos2.densityPerturbation * 100, self.v2_expected['perturbationPercent'], places = 5)
    
    def test_gasses(self):
        self.assertAlmostEqual(self.atmos.averageMolecularWeight, self.v1_expected['averageMolecularWeight'], places=4)
        self.assertAlmostEqual(self.atmos2.averageMolecularWeight, self.v2_expected['averageMolecularWeight'], places = 5)
        self.assertAlmostEqual(self.atmos.hydrogen.moleFraction * 100, self.v1_expected['hydrogenMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos2.hydrogen.moleFraction * 100, self.v2_expected['hydrogenMoleFraction'], places=4)
        self.assertAlmostEqual(self.atmos.helium.moleFraction * 100, self.v1_expected['heliumMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos2.helium.moleFraction * 100, self.v2_expected['heliumMoleFraction'], places=4)
        self.assertAlmostEqual(self.atmos.oxygen.moleFraction * 100, self.v1_expected['oxygenMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos2.oxygen.moleFraction * 100, self.v2_expected['oxygenMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos.nitrogen.moleFraction * 100, self.v1_expected['nitrogenMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos2.nitrogen.moleFraction * 100, self.v2_expected['nitrogenMoleFraction'], places=7)
        self.assertAlmostEqual(self.atmos.dinitrogen.moleFraction * 100, self.v1_expected['dinitrogenMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos2.dinitrogen.moleFraction * 100, self.v2_expected['dinitrogenMoleFraction'], places=8)
        self.assertAlmostEqual(self.atmos.carbonMonoxide.moleFraction * 100, self.v1_expected['carbonMonoxideMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos2.carbonMonoxide.moleFraction * 100, self.v2_expected['carbonMonoxideMoleFraction'], places=8)
        self.assertAlmostEqual(self.atmos.carbonDioxide.moleFraction * 100, self.v1_expected['carbonDioxideMoleFraction'], places=4)
        self.assertAlmostEqual(self.atmos2.carbonDioxide.moleFraction * 100, self.v2_expected['carbonDioxideMoleFraction'], places=11)

    
    def test_winds(self):
        self.assertAlmostEqual(self.atmos.ewWind, self.v1_expected['EWWind'], places=4)
        self.assertAlmostEqual(self.atmos2.ewWind, self.v2_expected['EWWind'], places=3)
        self.assertAlmostEqual(self.atmos.nsWind, self.v1_expected['NSWind'], places=5)
        self.assertAlmostEqual(self.atmos2.nsWind, self.v2_expected['NSWind'], places=4)
        self.assertAlmostEqual(self.atmos.perturbedEWWind, self.v1_expected['perturbedEWWind'], places=4)
        self.assertAlmostEqual(self.atmos2.perturbedEWWind, self.v2_expected['perturbedEWWind'], places=3)
        self.assertAlmostEqual(self.atmos.perturbedNSWind, self.v1_expected['perturbedNSWind'], places=6)
        self.assertAlmostEqual(self.atmos2.perturbedNSWind, self.v2_expected['perturbedNSWind'], places=4)
         

if __name__ == '__main__':
    unittest.main()