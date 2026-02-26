import unittest
from gram import * 

class TitanAtmosphereTestCase(unittest.TestCase):

    def setUp(self):
        # load path to SPICE data (note: make sure spice.txt copied to appropriate locations)
        inputParameters = InputParameters()
        NamelistReader().tryGetSpicePath(inputParameters)
        SpiceLoader().setInputParameters(inputParameters)

    def test_input_parameters(self):

        titan = TitanAtmosphere() 
        inputParameters = TitanInputParameters()
        reader = TitanNamelistReader()
        reader.tryGetSpicePath(inputParameters)
        titan.setInputParameters(inputParameters)
        self.assertEqual(str(titan.getInputParameters().spicePath), str(inputParameters.spicePath))
        inputParameters.minute = 33
        self.assertNotEqual(titan.getInputParameters().minute, inputParameters.minute)
        titan.setInputParameters(inputParameters)
        self.assertEqual(titan.getInputParameters().minute, 33)

class TitanSystemTestCase(unittest.TestCase):
    # system test comparing outputs from Python interface to outputs from C++ interface
    # Reference Titan/examples/Titan.cpp

    t1_expected = {
        'totalRadius': 2.625500e+03, 
        'gravity': 1.302425e+00, 
        'solarTime': 2.186178e+01, 
        'longitudeSun': 1.208846e+02, 
        'temperature': 7.118267e+01, 
        'pressure': 7.564649e+03, 
        'density': 3.565945e-01, 
        'pressureScaleHeight': 1.624122e+01, 
        'densityScaleHeight': 1.579936e+01, 
        'meanDensity': 3.565945e-01, 
        'perturbedDensity': 3.841655e-01, 
        'perturbationPercent': 7.731769e+00, 
        'averageMolecularWeight': 2.778150e+01, 
        'argonMoleFraction': 1.742032e+00, 
        'methaneMoleFraction': 3.040221e+00, 
        'dinitrogenMoleFraction': 9.521775e+01, 
        'EWWind': 4.068481e+01, 
        'perturbedEWWind': 2.381389e+01
    }

    t2_expected = {
        'totalRadius': 2.625500e+03, 
        'gravity': 1.302425e+00, 
        'solarTime':2.186178e+01, 
        'longitudeSun':1.208846e+02, 
        'temperature': 7.385170e+01,
        'pressure': 7.748606e+03, 
        'density': 3.513546e-01, 
        'pressureScaleHeight': 1.719681e+01, 
        'densityScaleHeight': 1.533102e+01, 
        'meanDensity': 3.513546e-01, 
        'perturbedDensity': 3.785205e-01, 
        'perturbationPercent': 7.731769e+00, 
        'averageMolecularWeight': 2.780800e+01, 
        'argonMoleFraction': 2.000000e+00, 
        'methaneMoleFraction': 3.000000e+00, 
        'dinitrogenMoleFraction': 9.500000e+01,
        'EWWind': 2.423240e+01, 
        'perturbedEWWind': 7.361473e+00
    }

    def setUp(self):

        inputParameters = TitanInputParameters()
        reader = TitanNamelistReader()
        reader.tryGetSpicePath(inputParameters)

        self.titan = TitanAtmosphere()
        self.titan.setInputParameters(inputParameters)

        self.titan2 = TitanAtmosphere()
        self.titan2.setInputParameters(inputParameters)

        self.titan.setModelType(Yelle97)
        self.titan.setMinMaxFactor(0.0,True)
        self.titan2.setModelType(GCM95)

        self.titan.setPerturbationScales(1.5)
        self.titan.setMinRelativeStepSize(0.5)
        self.titan.setSeed(1001)
        self.titan2.setPerturbationScales(1.5)
        self.titan2.setMinRelativeStepSize(0.5)
        self.titan2.setSeed(1001)

        ttime = GramTime()
        ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT)
        self.titan.setStartTime(ttime)
        self.titan2.setStartTime(ttime)

        position = Position()
        position.height = 50
        position.latitude = 22
        position.longitude = 48
        position.elapsedTime = 100
        self.titan.setPosition(position)
        self.titan2.setPosition(position)

        self.titan.update()
        self.titan2.update()

        self.pos = self.titan.getPosition()
        self.pos2 = self.titan2.getPosition()

        self.estate = self.titan.getEphemerisState()
        self.estate2 = self.titan2.getEphemerisState()

        self.atmos = self.titan.getAtmosphereState()
        self.atmos2 = self.titan2.getAtmosphereState()



    # testing equality to penultimate digit of expected output. 
    # This will lead to us testing some figures that are not physically meaningful, but here
    # we are testing that the calculations done via the Python interface match with those done 
    # done directly in C++, and nothing else
    def test_position(self):
        self.assertAlmostEqual(self.pos.totalRadius, self.t1_expected['totalRadius'], places=1)
        self.assertAlmostEqual(self.pos2.totalRadius, self.t2_expected['totalRadius'], places=1)

        self.assertAlmostEqual(self.pos.gravity, self.t1_expected['gravity'], places=5)
        self.assertAlmostEqual(self.pos2.gravity, self.t2_expected['gravity'], places=5)

    def test_ephemeris_state(self):
        self.assertAlmostEqual(self.estate.solarTime, self.t1_expected['solarTime'], places=4)
        self.assertAlmostEqual(self.estate2.solarTime, self.t2_expected['solarTime'], places=4)
        self.assertAlmostEqual(self.estate.longitudeSun, self.t1_expected['longitudeSun'], places=3)
        self.assertAlmostEqual(self.estate2.longitudeSun, self.t2_expected['longitudeSun'], places=3)
        
    def test_atmosphere_state(self):
        self.assertAlmostEqual(self.atmos.temperature, self.t1_expected['temperature'], places=4)
        self.assertAlmostEqual(self.atmos2.temperature, self.t2_expected['temperature'], places=4)
        self.assertAlmostEqual(self.atmos.pressure, self.t1_expected['pressure'], places=2)
        self.assertAlmostEqual(self.atmos2.pressure, self.t2_expected['pressure'], places=2)
        self.assertAlmostEqual(self.atmos.density, self.t1_expected['density'], places=6)
        self.assertAlmostEqual(self.atmos2.density, self.t2_expected['density'], places=6)
        self.assertAlmostEqual(self.atmos.pressureScaleHeight, self.t1_expected['pressureScaleHeight'], places=4)
        self.assertAlmostEqual(self.atmos2.pressureScaleHeight, self.t2_expected['pressureScaleHeight'], places=4)
        self.assertAlmostEqual(self.atmos.densityScaleHeight, self.t1_expected['densityScaleHeight'], places=4)
        self.assertAlmostEqual(self.atmos2.densityScaleHeight, self.t2_expected['densityScaleHeight'], places=4)
        self.assertAlmostEqual(self.atmos.perturbedDensity, self.t1_expected['perturbedDensity'], places=6)
        self.assertAlmostEqual(self.atmos2.perturbedDensity, self.t2_expected['perturbedDensity'], places=6)
        self.assertAlmostEqual(self.atmos.densityPerturbation * 100, self.t1_expected['perturbationPercent'], places=5)
        self.assertAlmostEqual(self.atmos2.densityPerturbation * 100, self.t2_expected['perturbationPercent'], places = 5)
    
    def test_gasses(self):
        self.assertAlmostEqual(self.atmos.averageMolecularWeight, self.t1_expected['averageMolecularWeight'], places = 4)
        self.assertAlmostEqual(self.atmos2.averageMolecularWeight, self.t2_expected['averageMolecularWeight'], places = 4)
        self.assertAlmostEqual(self.atmos.methane.moleFraction * 100, self.t1_expected['methaneMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos2.methane.moleFraction * 100, self.t2_expected['methaneMoleFraction'], places=5)
        self.assertAlmostEqual(self.atmos.dinitrogen.moleFraction * 100, self.t1_expected['dinitrogenMoleFraction'], places=4)
        self.assertAlmostEqual(self.atmos2.dinitrogen.moleFraction * 100, self.t2_expected['dinitrogenMoleFraction'], places=4)
    
    def test_winds(self):
        self.assertAlmostEqual(self.atmos.ewWind, self.t1_expected['EWWind'], places=4)
        self.assertAlmostEqual(self.atmos2.ewWind, self.t2_expected['EWWind'], places=4)
        self.assertAlmostEqual(self.atmos.perturbedEWWind, self.t1_expected['perturbedEWWind'], places=4)
        self.assertAlmostEqual(self.atmos2.perturbedEWWind, self.t2_expected['perturbedEWWind'], places=5)
        

if __name__ == '__main__':
    unittest.main()
