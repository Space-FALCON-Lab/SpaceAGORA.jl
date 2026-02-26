#include <vector>
#include "unittest.h"
#include "Atmosphere.h"
#include "PerturbedAtmosphere.h"


using namespace GRAM;

// PerturbedAtmosphere is abstract, so make a test class.
//! \cond Hide_this_from_doxygen
class PerturbedAtmosphereTest : public PerturbedAtmosphere
{
public:
	using PerturbedAtmosphere::normalPercentagePoint;
	//using PerturbedAtmosphere::uniformRandomNumber;
  using PerturbedAtmosphere::randomNumberGenerator;
 // using PerturbedAtmosphere::ix;
 // using PerturbedAtmosphere::iy;
	//using PerturbedAtmosphere::iz;
	using PerturbedAtmosphere::RHOd;
	using PerturbedAtmosphere::RHOu;
	using PerturbedAtmosphere::RHOv;

	PerturbedAtmosphereTest() { }
	virtual ~PerturbedAtmosphereTest() override = default;

	virtual void update() override {}
	
	InputParameters ip;
	virtual const InputParameters& getInputParameters() const override { return ip;};
	virtual void getPerturbationFactors(greal& pertLow, greal& pertHigh) override {};
	virtual void getScaleParameters(greal& verticalScale, greal& horizontalScale) override {};
	virtual void getWindDeviations(greal& ewStdDev, greal& nsStdDev, greal& vertStdDev) override {};
	
};
//! \endcond

TEST(PerturbedAtmosphere, normalPercentagePoint)
{
	// SETUP
	PerturbedAtmosphereTest a;
	greal res = 0;
	// GIVEN (INPUTS)
	greal p[10] = {0.4376, 0.0159, 0.8772, 0.0021, 0.0359,
		0.9192, 0.811, 0.0075, 0.0625, 0.2876 };

	// RUN & EXPECT (OUTPUTS)
	greal val[10] = {-0.1570569050554871, -2.146915617426713,
		1.161103038740104, -2.862736224507694,
		-1.800384161311344, 1.399710593591071,
		0.8815873475968675, -2.432379062093879,
		-1.534120543709685, -0.5604097234915217 };

	for (int i = 0; i < 10; i++)
	{
		res = a.normalPercentagePoint(p[i]);
		EXPECT_DOUBLE_EQ(val[i], res);
	}
}

TEST(PerturbedAtmosphere, setSeed)
{
	// SETUP
	PerturbedAtmosphereTest a;

	// GIVEN (INPUTS)
	int s[5] = { 2511, 43, 7811, 167, 29011 };

	// RUN & EXPECT (OUTPUTS)
	greal ansd[5] = { 0.60550707547150651, -0.71481041302557413,
   -0.36437824971788441, 0.45842156855761607, 0.58886020757009938 };
	greal ansu[5] = { 1.3645917463845345, -0.44475853957367045,
   -0.34451299660974899, 1.3592658206044466, 0.05437436017925816 };
	greal ansv[5] = { 0.53772294225549777, -1.2063756109420474,
    1.0731106270573088, 1.0134837303802571, 0.84817847781978317 };
	for (int i = 0; i < 5; i++)
	{
		a.setSeed(s[i]);
		EXPECT_NEAR(ansd[i], a.RHOd, 0.0000000001);
		EXPECT_NEAR(ansu[i], a.RHOu, 0.0000000001);
		EXPECT_NEAR(ansv[i], a.RHOv, 0.0000000001);
	}
}


//TEST(PerturbedAtmosphere, uniformRandomNumber)
//{
//	// SETUP
//
//	PerturbedAtmosphereTest a;
//	
//	// GIVEN (INPUTS)
//	int s[10] = {2511, 43, 7811, 167, 29011, 551, 1211, 110, 4511, 10020};
//
//	// RUN & EXPECT (OUTPUTS)
//	greal ans[10] = {0.7867453311953283, 0.9035563716612502,
//		0.8762516057215088, 9.055614110299492E-02, 0.234276703826231,
//		0.5083618787290426, 0.6327154902738124, 0.5207256019241284,
//		0.2544835479976606, 0.4333684661796848 };
//	
//  std::vector<greal> rvec(1);
//	for (int i = 0; i < 10; i++)
//	{
//		a.setSeed(s[i]);
//		a.randomNumberGenerator.getRandomNumbers(rvec);
//		EXPECT_DOUBLE_EQ(ans[i], rvec[0]);
//	}
//}
