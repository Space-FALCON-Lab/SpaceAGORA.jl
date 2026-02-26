import unittest
import ctypes
from gram import *

# Using JupiterAtmosphere because Atmosphere is an abstract class 
class AtmosphereTestCase(unittest.TestCase):

    def test_get_ephemeris_state(self):

        estate = EphemerisState()
        atmos = JupiterAtmosphere()
        atmos.setEphemerisState(estate)

        self.assertEqual(atmos.getEphemerisState().solarTime, estate.solarTime)

        estate.solarTime = 33
        self.assertNotEqual(atmos.getEphemerisState().solarTime, estate.solarTime)

        atmos.setEphemerisState(estate)
        self.assertEqual(atmos.getEphemerisState().solarTime, 33)

    def test_get_constituent_gas(self):
        atmos = JupiterAtmosphere()
        hydrogen = atmos.getConstituentGas(HYDROGEN)
        self.assertIsInstance(hydrogen, ConstituentGas)

    def test_set_position(self):
        atmos = JupiterAtmosphere()
        pos = Position()
        pos.height = 33
        atmos.setPosition(pos)
        self.assertEqual(atmos.getPosition().height, 33)

class PerturbedAtmosphereTestCase(unittest.TestCase):

    def setUp(self):
        inputParameters = InputParameters()
        NamelistReader().tryGetSpicePath(inputParameters)
        SpiceLoader().setInputParameters(inputParameters)

    def test_set_start_time(self):

        atmos = JupiterAtmosphere()
        ttime = GramTime()
        ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT)
        atmos.setStartTime(ttime)

        self.assertIsInstance(atmos.getStartTime(), GramTime)
        self.assertEqual(atmos.getStartTime().getDayOfYear(), ttime.getDayOfYear())


class GramTimeTestCase(unittest.TestCase):

    def setUp(self):
        inputParameters = InputParameters()
        NamelistReader().tryGetSpicePath(inputParameters)
        SpiceLoader().setInputParameters(inputParameters)

    def test_set_start_time(self):
        ttime =GramTime()
        ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT)
        # This is just the output we get by doing the same directly in C++
        self.assertAlmostEqual(ttime.getSpiceTime(), 6.37502e+08, delta=500)

    def test_set_start_time_overload(self):
        ttime = GramTime()
        ttime.setStartTime(333, UTC, ERT)
        self.assertAlmostEqual(ttime.getSpiceTime(), -2.11785e+11, delta=300000)

    def test_get_time_scale_frame(self):
        ttime = GramTime()
        ttime.setStartTime(2020, 3, 15, 0, 0, 0.0, UTC, ERT)
        self.assertEqual(ttime.getTimeScale(), UTC)
        self.assertEqual(ttime.getTimeFrame(), ERT)

    # This test passes when using ctypes as below, but this may not be expected behavior
    @unittest.skip
    def test_get_time(self):
        ttime = GramTime()
        ttime.setSpiceTime(1132971254.14)
        year = ctypes.c_int(0)
        month = ctypes.c_int(0)
        day = ctypes.c_int(0)
        hour = ctypes.c_int(0)
        minute = ctypes.c_int(0)
        seconds = ctypes.c_double(0)
        # use ctypes.c_int for pass-by-ref of ints
        ttime.getTime(TDB, PET, year, month, day, hour, minute, seconds)
        self.assertEqual(year.value, 2035)
        self.assertEqual(month.value, 11)
        self.assertEqual(day.value, 26)
        self.assertEqual(hour.value, 14)
        self.assertEqual(minute.value, 14)
        self.assertAlmostEqual(seconds.value, 14.14, places=2)

    # This test passes when using ctypes as below, but this may not be expected behavior
    @unittest.skip
    def test_get_time_julian_date(self):
        # same issue with pass by reference ints as above
        ttime= GramTime()
        julianDate = 0
        ttime.setSpiceTime(1132971254.14)
        ttime.getTime(TDB, PET, julianDate)
        self.assertAlmostEqual(julianDate, 2464658.0932192, places=6)

        julianDate = 0
        ttime.setSpiceTime(1132971254.1389618)
        ttime.getTime(TDT, PET, julianDate)
        self.assertAlmostEqual(julianDate, 2464658.0932192, places=6)

        julianDate - 0
        ttime.setSpiceTime(132971323.3229617)
        ttime.getTime(UTC, PET, julianDate)
        self.assertAlmostEqual(julianDate, 2464658.0932192)


class PositionTestCase(unittest.TestCase):

    def test_get_set_longitude(self):
        pos = Position()
        pos.setLongitude(100)
        self.assertEqual(pos.getLongitude(), 100)


class ConstituentGasTestCase(unittest.TestCase):

    def test_update_specific_heat_capacity_argon(self):
        # Just testing one gas type here - if interface handles inputs/ouputs correctly, 
        # that will hold for other gas types because they have the same signature
        gas = ConstituentGas()

        temperature = [30.0,  150.0,  450.0, 1111.0, 432.1]
        pressure = [0.04,  1.0e4,  3.3e7,  2.2e9, 6.7e5]
        expected_cp = [0.5203, 0.5210, 0.5637, 0.5280, 0.52372722]
        expected_cv = [0.3122, 0.3124, 0.3187, 0.3140, 0.3127253466666667]

        for i in range(len(temperature)):
            gas.updateSpecificHeatCapacity(ARGON, temperature[i], pressure[i])

            self.assertAlmostEqual(gas.specificHeatPressure, expected_cp[i])
            self.assertAlmostEqual(gas.specificHeatVolume, expected_cv[i])

if __name__ == '__main__':
    unittest.main()