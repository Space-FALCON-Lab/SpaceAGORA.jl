import unittest
from gram import *

MARS_DATA_PATH = "../../Mars/data/"

class MarsAtmosphereTestCase(unittest.TestCase):

    def setUp(self):

        inputParameters = MarsInputParameters()
        inputParameters.dataPath = MARS_DATA_PATH
        # load path to SPICE data (note: make sure spice.txt copied to appropriate locations)
        NamelistReader().tryGetSpicePath(inputParameters)
        SpiceLoader().setInputParameters(inputParameters)

    def test_input_parameters(self):

        mars = MarsAtmosphere() 
        inputParameters = MarsInputParameters()
        inputParameters.dataPath = MARS_DATA_PATH
        reader = MarsNamelistReader()
        reader.tryGetSpicePath(inputParameters)
        mars.setInputParameters(inputParameters)
        self.assertEqual(str(mars.getInputParameters().spicePath), str(inputParameters.spicePath))
        inputParameters.minute = 33
        self.assertNotEqual(mars.getInputParameters().minute, inputParameters.minute)
        mars.setInputParameters(inputParameters)
        self.assertEqual(mars.getInputParameters().minute, 33)

class MarsSystemTestCase(unittest.TestCase):
    m1_expected = {
        'totalRadius': 3.395766e+03,
        'gravity': 3.705779e+00, 
        'solarTime': 2.051137e+01, 
        'longitudeSun': 1.664260e+02, 
        'temperature': 2.105041e+02, 
        'pressure': 4.312749e+02, 
        'density': 1.084362e-02, 
        'pressureScaleHeight': 1.074235e+01, 
        'densityScaleHeight': 8.754779e+00, 
        'perturbedDensity': 1.097261e-02, 
        'perturbationPercent': 1.189503e+00, 
        'averageMolecularWeight': 4.338581e+01, 
        'carbonDioxideMoleFraction': 9.468056e+01, 
        'dinitrogenMoleFraction': 3.139150e+00, 
        'argonMoleFraction': 1.883490e+00, 
        'carbonMonoxideMoleFraction': 1.141509e-01, 
        'dioxygenMoleFraction': 1.483962e-01, 
        'dihydrogenMoleFraction': 0.000000e+00, 
        'hydrogenMoleFraction': 0.000000e+00, 
        'oxygenMoleFraction': 0.000000e+00, 
        'heliumMoleFraction': 0.000000e+00, 
        'waterMoleFraction': 3.426654e-02, 
        'EWWind': 2.255073e-01, 
        'NSWind': -3.719453e+00, 
        'verticalWind': 1.807776e-03, 
        'perturbedEWWind': -3.589672e+00, 
        'perturbedNSWind': -4.425719e+00, 
        'dustOpticalDepth': 2.000000e+00, 
        'temperatureDaily': 2.160810e+02, 
        'temperatureMax': 2.419776e+02, 
        'temperatureMin': 1.977304e+02
    }

    m2_expected = {
        'totalRadius': 3.395766e+03,
        'gravity': 3.705779e+00,
        'solarTime': 2.051137e+01,
        'longitudeSun': 1.664260e+02,
        'temperature': 2.238048e+02,
        'pressure': 4.359997e+02,
        'density': 1.026755e-02,
        'pressureScaleHeight': 1.140048e+01,
        'densityScaleHeight': 1.310461e+01,
        'perturbedDensity': 1.038969e-02,
        'perturbationPercent': 1.189503e+00,
        'averageMolecularWeight': 4.334565e+01,
        'carbonDioxideMoleFraction': 9.453043e+01,
        'dinitrogenMoleFraction': 3.134173e+00,
        'argonMoleFraction': 1.880504e+00,
        'carbonMonoxideMoleFraction': 1.139699e-01,
        'dioxygenMoleFraction': 1.481609e-01,
        'dihydrogenMoleFraction': 0.000000e+00,
        'hydrogenMoleFraction': 0.000000e+00,
        'oxygenMoleFraction': 0.000000e+00,
        'heliumMoleFraction': 0.000000e+00,
        'waterMoleFraction': 1.931316e-01,
        'EWWind': -2.060489e-01,
        'NSWind': -6.664234e+00,
        'verticalWind': 2.299431e-03,
        'perturbedEWWind': -4.021228e+00,
        'perturbedNSWind': -7.370499e+00,
        'dustOpticalDepth': 1.956481e-01,
        'temperatureDaily': 2.199744e+02,
        'temperatureMax': 2.282888e+02,
        'temperatureMin': 2.132719e+02
    }

    

    def setUp(self):
        inputParameters = MarsInputParameters()
        reader = MarsNamelistReader()
        inputParameters.dataPath = MARS_DATA_PATH
        reader.tryGetSpicePath(inputParameters)

        self.mars = MarsAtmosphere()
        self.mars.setInputParameters(inputParameters)

        self.mars2 = MarsAtmosphere()
        self.mars2.setInputParameters(inputParameters)

        self.mars.setMapYear(0)
        self.mars.setMGCMDustLevels(2.0, 0.0, 0.0)
        self.mars2.setMapYear(2)

        self.mars.setPerturbationScales(1.5)
        self.mars.setMinRelativeStepSize(0.5)
        self.mars.setSeed(1001)
        self.mars2.setPerturbationScales(1.5)
        self.mars2.setMinRelativeStepSize(0.5)
        self.mars2.setSeed(1001)

        ttime = GramTime()
        ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT)
        self.mars.setStartTime(ttime)
        self.mars2.setStartTime(ttime)

        position = Position()
        position.height = 2
        position.latitude = 22
        position.longitude = 48
        position.elapsedTime = 100
        position.isPlanetoCentric = True
        self.mars.setPosition(position)
        self.mars2.setPosition(position)

        self.mars.update()
        self.mars2.update()

        self.pos = self.mars.getPosition()
        self.pos2 = self.mars2.getPosition()

        self.estate = self.mars.getEphemerisState()
        self.estate2 = self.mars2.getEphemerisState()

        self.atmos = self.mars.getAtmosphereState()
        self.atmos2 = self.mars2.getAtmosphereState()

        self.marsAtmos = self.mars.getMarsAtmosphereState()
        self.marsAtmos2 = self.mars2.getMarsAtmosphereState()


    def test_position(self):
        self.assertAlmostEqual(self.pos.totalRadius, self.m1_expected['totalRadius'], places=2)
        self.assertAlmostEqual(self.pos2.totalRadius, self.m2_expected['totalRadius'], places=2)
        self.assertAlmostEqual(self.pos.gravity, self.m1_expected['gravity'], places=5)
        self.assertAlmostEqual(self.pos2.gravity, self.m2_expected['gravity'], places=5)

    def test_ephemeris_state(self):
        self.assertAlmostEqual(self.estate.solarTime, self.m1_expected['solarTime'], places=4)
        self.assertAlmostEqual(self.estate2.solarTime, self.m2_expected['solarTime'], places=4)
        self.assertAlmostEqual(self.estate.longitudeSun, self.m1_expected['longitudeSun'], places=3)
        self.assertAlmostEqual(self.estate2.longitudeSun, self.m2_expected['longitudeSun'], places=3)
        
    def test_atmosphere_state(self):
        self.assertAlmostEqual(self.atmos.temperature, self.m1_expected['temperature'], places=3)
        self.assertAlmostEqual(self.atmos2.temperature, self.m2_expected['temperature'], places=3)
        self.assertAlmostEqual(self.atmos.pressure, self.m1_expected['pressure'], places=3)
        self.assertAlmostEqual(self.atmos2.pressure, self.m2_expected['pressure'], places=3)
        self.assertAlmostEqual(self.atmos.density, self.m1_expected['density'], places=7)
        self.assertAlmostEqual(self.atmos2.density, self.m2_expected['density'], places=7)
        self.assertAlmostEqual(self.atmos.pressureScaleHeight, self.m1_expected['pressureScaleHeight'], places=4)
        self.assertAlmostEqual(self.atmos2.pressureScaleHeight, self.m2_expected['pressureScaleHeight'], places=4)
        self.assertAlmostEqual(self.atmos.densityScaleHeight, self.m1_expected['densityScaleHeight'], places=5)
        self.assertAlmostEqual(self.atmos2.densityScaleHeight, self.m2_expected['densityScaleHeight'], places=4)
        self.assertAlmostEqual(self.atmos.perturbedDensity, self.m1_expected['perturbedDensity'], places=7)
        self.assertAlmostEqual(self.atmos2.perturbedDensity, self.m2_expected['perturbedDensity'], places=7)
        self.assertAlmostEqual(self.atmos.densityPerturbation * 100, self.m1_expected['perturbationPercent'], places=5)
        self.assertAlmostEqual(self.atmos2.densityPerturbation * 100, self.m2_expected['perturbationPercent'], places = 5)
    
    def test_gasses(self):
        self.assertAlmostEqual(self.atmos.averageMolecularWeight, self.m1_expected['averageMolecularWeight'], places = 4)
        self.assertAlmostEqual(self.atmos2.averageMolecularWeight, self.m2_expected['averageMolecularWeight'], places = 4)
        self.assertAlmostEqual(self.atmos.carbonDioxide.moleFraction * 100, self.m1_expected['carbonDioxideMoleFraction'], places=4)
        self.assertAlmostEqual(self.atmos2.carbonDioxide.moleFraction * 100, self.m2_expected['carbonDioxideMoleFraction'], places=4)
        self.assertAlmostEqual(self.atmos.dinitrogen.moleFraction * 100, self.m1_expected['dinitrogenMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos2.dinitrogen.moleFraction * 100, self.m2_expected['dinitrogenMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos.argon.moleFraction * 100, self.m1_expected['argonMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos2.argon.moleFraction * 100, self.m2_expected['argonMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos.carbonMonoxide.moleFraction * 100, self.m1_expected['carbonMonoxideMoleFraction'], places=6)
        self.assertAlmostEqual(self.atmos2.carbonMonoxide.moleFraction * 100, self.m2_expected['carbonMonoxideMoleFraction'], places=6)
        self.assertAlmostEqual(self.atmos.dioxygen.moleFraction * 100, self.m1_expected['dioxygenMoleFraction'], places=6)
        self.assertAlmostEqual(self.atmos2.dioxygen.moleFraction * 100, self.m2_expected['dioxygenMoleFraction'], places=6)
        self.assertAlmostEqual(self.atmos.dihydrogen.moleFraction * 100, self.m1_expected['dihydrogenMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos2.dihydrogen.moleFraction * 100, self.m2_expected['dihydrogenMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos.hydrogen.moleFraction * 100, self.m1_expected['hydrogenMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos2.hydrogen.moleFraction * 100, self.m2_expected['hydrogenMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos.oxygen.moleFraction * 100, self.m1_expected['oxygenMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos2.oxygen.moleFraction * 100, self.m2_expected['oxygenMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos.helium.moleFraction * 100, self.m1_expected['heliumMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos2.helium.moleFraction * 100, self.m2_expected['heliumMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos.water.moleFraction * 100, self.m1_expected['waterMoleFraction'], places=7)
        self.assertAlmostEqual(self.atmos2.water.moleFraction * 100, self.m2_expected['waterMoleFraction'], places=6)

    
    def test_winds(self):
        self.assertAlmostEqual(self.atmos.ewWind, self.m1_expected['EWWind'], places=6)
        self.assertAlmostEqual(self.atmos2.ewWind, self.m2_expected['EWWind'], places=6)
        self.assertAlmostEqual(self.atmos.nsWind, self.m1_expected['NSWind'], places=5)
        self.assertAlmostEqual(self.atmos2.nsWind, self.m2_expected['NSWind'], places=5)
        self.assertAlmostEqual(self.atmos.verticalWind, self.m1_expected['verticalWind'], places=8)
        self.assertAlmostEqual(self.atmos2.verticalWind, self.m2_expected['verticalWind'], places=8)
        self.assertAlmostEqual(self.atmos.perturbedEWWind, self.m1_expected['perturbedEWWind'], places=5)
        self.assertAlmostEqual(self.atmos2.perturbedEWWind, self.m2_expected['perturbedEWWind'], places=5)
        self.assertAlmostEqual(self.atmos.perturbedNSWind, self.m1_expected['perturbedNSWind'], places=5)
        self.assertAlmostEqual(self.atmos2.perturbedNSWind, self.m2_expected['perturbedNSWind'], places=5)

    
    def test_mars_specific(self):
        self.assertAlmostEqual(self.marsAtmos.dustOpticalDepth, self.m1_expected['dustOpticalDepth'], places=5)
        self.assertAlmostEqual(self.marsAtmos2.dustOpticalDepth, self.m2_expected['dustOpticalDepth'], places=6)
        self.assertAlmostEqual(self.marsAtmos.temperatureDaily, self.m1_expected['temperatureDaily'], places=3)
        self.assertAlmostEqual(self.marsAtmos2.temperatureDaily, self.m2_expected['temperatureDaily'], places=3)
        self.assertAlmostEqual(self.marsAtmos.temperatureMax, self.m1_expected['temperatureMax'], places=3)
        self.assertAlmostEqual(self.marsAtmos2.temperatureMax, self.m2_expected['temperatureMax'], places=3)
        self.assertAlmostEqual(self.marsAtmos.temperatureMin, self.m1_expected['temperatureMin'], places=3)
        self.assertAlmostEqual(self.marsAtmos2.temperatureMin, self.m2_expected['temperatureMin'], places=3)
        
if __name__ == '__main__':
    unittest.main()
        
