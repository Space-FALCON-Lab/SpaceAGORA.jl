#include "unittest.h"
//Not the best way to get to private member auxDataTable
#define private protected
#include "AuxiliaryAtmosphere.h"
#undef private

namespace GRAM {

	TEST(AuxiliaryAtmosphere, HASINPUT(loadData))
	{
		class AuxiliaryAtmosphereExposed : public AuxiliaryAtmosphere
		{
		public:
			using AuxiliaryAtmosphere::auxDataTable;
		};
		// SETUP
		AuxiliaryAtmosphereExposed atmo;

		// GIVEN (INPUTS)
		// need to put file in $(OutDir) defined in compiler
		atmo.setDataFile("auxatm_test.txt");

		// RUN
		atmo.loadData();

		// EXPECTED (OUTPUTS)
		EXPECT_DOUBLE_EQ(1500.0, atmo.auxDataTable[150].height);
		EXPECT_DOUBLE_EQ(67.0, atmo.auxDataTable[150].latitude);
		EXPECT_DOUBLE_EQ(123.0, atmo.auxDataTable[150].longitude);
		EXPECT_DOUBLE_EQ(523.0, atmo.auxDataTable[150].temperature);
		EXPECT_DOUBLE_EQ(2.91E-12, atmo.auxDataTable[150].pressure);
		EXPECT_DOUBLE_EQ(6.18E-06, atmo.auxDataTable[150].density);
		EXPECT_DOUBLE_EQ(435.9, atmo.auxDataTable[150].ewWind);
		EXPECT_DOUBLE_EQ(0.0, atmo.auxDataTable[150].nsWind);
		// TEAR DOWN
		

	}
	

} //namespace