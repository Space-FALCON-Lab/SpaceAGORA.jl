import unittest
from gram import * 

class NeptuneAtmosphereTestCase(unittest.TestCase):

    def setUp(self):
        # load path to SPICE data (note: make sure spice.txt copied to appropriate locations)
        inputParameters = InputParameters()
        NamelistReader().tryGetSpicePath(inputParameters)
        SpiceLoader().setInputParameters(inputParameters)

    def test_input_parameters(self):

        neptune = NeptuneAtmosphere() 
        inputParameters = NeptuneInputParameters()
        reader = NeptuneNamelistReader()
        reader.tryGetSpicePath(inputParameters)
        neptune.setInputParameters(inputParameters)
        self.assertEqual(str(neptune.getInputParameters().spicePath), str(inputParameters.spicePath))
        inputParameters.minute = 33
        self.assertNotEqual(neptune.getInputParameters().minute, inputParameters.minute)
        neptune.setInputParameters(inputParameters)
        self.assertEqual(neptune.getInputParameters().minute, 33)

class NeptuneSystemTestCase(unittest.TestCase):

    n1_expected = {
        'totalRadius': 2.475331e+04, 
        'gravity': 1.093854e+01, 
        'solarTime': 1.838570e+01, 
        'longitudeSun': 3.024877e+02, 
        'temperature': 5.337300e+01, 
        'pressure': 4.698319e+03,
        'density': 2.762573e-02, 
        'pressureScaleHeight': 1.537637e+01, 
        'densityScaleHeight': 1.437031e+01, 
        'meanDensity': 2.762573e-02, 
        'perturbedDensity': 2.928931e-02, 
        'perturbationPercent': 6.021858e+00, 
        'averageMolecularWeight': 2.608515e+00, 
        'dihydrogenMoleFraction': 7.973899e+01, 
        'methaneMoleFraction': 1.556037e+00, 
        'heliumMoleFraction': 1.870498e+01,
        'dinitrogenMoleFraction': 0.000000e+00, 
        'EWWind': -2.774529e+02, 
        'perturbedEWWind': -3.160150e+02
    }

    n2_expected = {
        'totalRadius': 2.475331e+04, 
        'gravity': 1.093854e+01, 
        'solarTime': 1.838570e+01, 
        'longitudeSun': 3.024877e+02, 
        'temperature': 5.098401e+01, 
        'pressure': 4.304924e+03,
        'density': 2.649969e-02, 
        'pressureScaleHeight': 1.472063e+01, 
        'densityScaleHeight': 1.398864e+01, 
        'meanDensity': 2.649969e-02, 
        'perturbedDensity': 2.656020e-02, 
        'perturbationPercent': 2.283518e-01, 
        'averageMolecularWeight': 2.608515e+00, 
        'dihydrogenMoleFraction': 8.215301e+01, 
        'methaneMoleFraction': 1.556039e+00, 
        'heliumMoleFraction': 1.609180e+01,
        'dinitrogenMoleFraction': 1.991578e-01, 
        'EWWind': -2.774529e+02, 
        'perturbedEWWind': -2.671564e+02
    }

    def setUp(self):
        inputParameters = NeptuneInputParameters()
        reader = NeptuneNamelistReader()
        reader.tryGetSpicePath(inputParameters)

        self.neptune = NeptuneAtmosphere()
        self.neptune.setInputParameters(inputParameters)

        self.neptune2 = NeptuneAtmosphere()
        self.neptune2.setInputParameters(inputParameters)

        self.neptune.setMinMaxFactor(0.5, True)
        self.neptune2.setMinMaxFactor(-0.5, True)

        self.neptune2.setDinitrogenMoleFraction(0.002)

        self.neptune.setPerturbationScales(1.5)
        self.neptune.setMinRelativeStepSize(0.5)
        self.neptune.setSeed(1001)
        self.neptune2.setPerturbationScales(0.5)
        self.neptune2.setMinRelativeStepSize(0.5)
        self.neptune2.setSeed(3333)

        ttime = GramTime()
        ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT)
        self.neptune.setStartTime(ttime)
        self.neptune2.setStartTime(ttime)

        position = Position()
        position.height = 50.0
        position.latitude = 22.0
        position.longitude = 48.0
        position.elapsedTime = 100.0
        self.neptune.setPosition(position)
        position.latitude = -22.0
        self.neptune2.setPosition(position)

        self.neptune.update()
        self.neptune2.update()

        self.pos = self.neptune.getPosition()
        self.pos2 = self.neptune2.getPosition()

        self.estate = self.neptune.getEphemerisState()
        self.estate2 = self.neptune2.getEphemerisState()

        self.atmos = self.neptune.getAtmosphereState()
        self.atmos2 = self.neptune2.getAtmosphereState()

    def test_position(self):
        self.assertAlmostEqual(self.pos.totalRadius, self.n1_expected['totalRadius'], places=2)
        self.assertAlmostEqual(self.pos2.totalRadius, self.n2_expected['totalRadius'], places=2)

        self.assertAlmostEqual(self.pos.gravity, self.n1_expected['gravity'], places=5)
        self.assertAlmostEqual(self.pos2.gravity, self.n2_expected['gravity'], places=5)

    def test_ephemeris_state(self):
        self.assertAlmostEqual(self.estate.solarTime, self.n1_expected['solarTime'], places=5)
        self.assertAlmostEqual(self.estate2.solarTime, self.n2_expected['solarTime'], places=5)
        self.assertAlmostEqual(self.estate.longitudeSun, self.n1_expected['longitudeSun'], places=4)
        self.assertAlmostEqual(self.estate2.longitudeSun, self.n2_expected['longitudeSun'], places=4)
        
    def test_atmosphere_state(self):
        self.assertAlmostEqual(self.atmos.temperature, self.n1_expected['temperature'], places=5)
        self.assertAlmostEqual(self.atmos2.temperature, self.n2_expected['temperature'], places=5)
        self.assertAlmostEqual(self.atmos.pressure, self.n1_expected['pressure'], places=3)
        self.assertAlmostEqual(self.atmos2.pressure, self.n2_expected['pressure'], places=3)
        self.assertAlmostEqual(self.atmos.density, self.n1_expected['density'], places=8)
        self.assertAlmostEqual(self.atmos2.density, self.n2_expected['density'], places=8)
        self.assertAlmostEqual(self.atmos.pressureScaleHeight, self.n1_expected['pressureScaleHeight'], places=5)
        self.assertAlmostEqual(self.atmos2.pressureScaleHeight, self.n2_expected['pressureScaleHeight'], places=5)
        self.assertAlmostEqual(self.atmos.densityScaleHeight, self.n1_expected['densityScaleHeight'], places=5)
        self.assertAlmostEqual(self.atmos2.densityScaleHeight, self.n2_expected['densityScaleHeight'], places=5)
        self.assertAlmostEqual(self.atmos.perturbedDensity, self.n1_expected['perturbedDensity'], places=8)
        self.assertAlmostEqual(self.atmos2.perturbedDensity, self.n2_expected['perturbedDensity'], places=8)
        self.assertAlmostEqual(self.atmos.densityPerturbation * 100, self.n1_expected['perturbationPercent'], places=6)
        self.assertAlmostEqual(self.atmos2.densityPerturbation * 100, self.n2_expected['perturbationPercent'], places = 7)
    
    def test_gasses(self):
        self.assertAlmostEqual(self.atmos.averageMolecularWeight, self.n1_expected['averageMolecularWeight'], places = 6)
        self.assertAlmostEqual(self.atmos2.averageMolecularWeight, self.n2_expected['averageMolecularWeight'], places = 6)
        self.assertAlmostEqual(self.atmos.dihydrogen.moleFraction * 100, self.n1_expected['dihydrogenMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos2.dihydrogen.moleFraction * 100, self.n2_expected['dihydrogenMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos.methane.moleFraction * 100, self.n1_expected['methaneMoleFraction'], places=6)
        self.assertAlmostEqual(self.atmos2.methane.moleFraction * 100, self.n2_expected['methaneMoleFraction'], places=6)
        self.assertAlmostEqual(self.atmos.helium.moleFraction * 100, self.n1_expected['heliumMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos2.helium.moleFraction * 100, self.n2_expected['heliumMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos.dinitrogen.moleFraction * 100, self.n1_expected['dinitrogenMoleFraction'], places=6)
        self.assertAlmostEqual(self.atmos2.dinitrogen.moleFraction * 100, self.n2_expected['dinitrogenMoleFraction'], places=5)
    
    def test_winds(self):
        self.assertAlmostEqual(self.atmos.ewWind, self.n1_expected['EWWind'], places=4)
        self.assertAlmostEqual(self.atmos2.ewWind, self.n2_expected['EWWind'], places=4)
        self.assertAlmostEqual(self.atmos.perturbedEWWind, self.n1_expected['perturbedEWWind'], places=4)
        self.assertAlmostEqual(self.atmos2.perturbedEWWind, self.n2_expected['perturbedEWWind'], places=4)









        



if __name__ == '__main__':
    unittest.main()
