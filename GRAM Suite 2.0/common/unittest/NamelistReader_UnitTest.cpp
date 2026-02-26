#include <vector>
#include <string>
#include <algorithm>
#include "unittest.h"
#include "gram.h"
// needed to get access to private member - itemMap
#define private protected
#include "InputParameters.h"
#include "NamelistReader.h"
#undef private

using namespace std;
using namespace GRAM;

//! \cond Hide_this_from_doxygen
class NamelistReaderExposed : public NamelistReader
{
public:
  using NamelistReader::parseLine;
  using NamelistReader::findItem;
  using NamelistReader::getItem;
	using NamelistReader::getParameters;
	using NamelistReader::readFile;
	using NamelistReader::itemMap;
};
//! \endcond

TEST(NamelistReader, parseLine)
{
  // SETUP
  NamelistReaderExposed r;

  // GIVEN (INPUTS)
  vector<string> line = { 
    "MYVALUE = '/mypath/tothe/spicedata'",
    "MYVALUE = '/mypath/to the/spice data'",
    " myvalue='/mypath/to the/spice data'    ",
    " MyValue\t=\t'/mypath/to the/spice data'  //This is a comment",
    "\t myvalue=''",
    " myvalue = 'null'",
    "MyValue = 34.12   ",
    "\tMyValue\t\t\t= \t34.12\t ",
    "MyValue = 34.12  degrees ",
    " MyValue\t =\t 34   ",
    " MyValue = 34.   ",
  };

  vector<string> values = { 
    "/mypath/tothe/spicedata",
    "/mypath/to the/spice data",
    "/mypath/to the/spice data",
    "/mypath/to the/spice data",
    "null",
    "null",
    "34.12",
    "34.12",
    "34.12",
    "34",
    "34.",
  };

  for (size_t i = 0; i < line.size(); ++i) {
    // RUN
    r.parseLine(line[i]);
    std::string item = r.itemMap["MYVALUE"];

    // EXPECT (OUTPUTS)
    EXPECT_EQ(values[i], item);
  }

  // TEAR-DOWN
}

TEST(NamelistReader, findItem)
{
  // SETUP
  NamelistReaderExposed r;

  // GIVEN (INPUTS)
  std::string itemNames[33] = { "SPICEDIR", "LSTFL", "OUTFL", "TRAJFL", "PROFILE",
    "DATADIR", "IERT", "IUTC", "MONTH", "MDAY", "MYEAR", "IHOUR", "IMIN", "SEC",
    "NPOS", "LONEAST", "NR1", "NVARX", "NVARY", "LOGSCALE", "FLAT", "FLON",
    "FHGT", "DELHGT", "DELLAT", "DELLON", "DELTIME", "PROFNEAR", "PROFFAR",
    "RPSCALE", "NMONTE", "IUP", "CORLMIN" };
  std::string itemValues[33] = { DEFAULT_SPICE_PATH, "LIST.txt", "OUTPUT.txt", "TRAJDATA.txt",
    "null", "D:\\neptune\\data\\", "1", "1", "8", "25", "89", "0", "0", "0.0",
    "201", "0", "1001", "1", "0", "0", "22.0", "48.0", "0.0", "20.0", "0.3",
    "0.5", "500.0", "0.0", "0.0", "1.0", "1", "13", "0.0" };
  for (int i = 0; i < 33; i++)
    r.itemMap[itemNames[i]] = itemValues[i];

  // RUN & EXPECT (OUTPUTS)
  std::string item;
  for (int i = 0; i < 33; i++)
  {
    item = r.findItem(itemNames[i]);
    EXPECT_EQ(itemValues[i], item);
  }
  // TEST FOR ITEM NOT IN MAP
  item = r.findItem("NOINMAP");
  EXPECT_EQ("", item);
  EXPECT_NE("NOTINMAP", item);
  // TEAR-DOWN
}

TEST(NamelistReader, getItem_string)
{
	// SETUP
	NamelistReaderExposed r;

	// GIVEN (INPUTS)
	std::string itemNames[33] = { "SPICEDIR", "LSTFL", "OUTFL", "TRAJFL", "PROFILE",
		"DATADIR", "IERT", "IUTC", "MONTH", "MDAY", "MYEAR", "IHOUR", "IMIN", "SEC",
		"NPOS", "LONEAST", "NR1", "NVARX", "NVARY", "LOGSCALE", "FLAT", "FLON",
		"FHGT", "DELHGT", "DELLAT", "DELLON", "DELTIME", "PROFNEAR", "PROFFAR",
		"RPSCALE", "NMONTE", "IUP", "CORLMIN" };
	std::string itemValues[33] = { DEFAULT_SPICE_PATH, "LIST.txt", "OUTPUT.txt", "TRAJDATA.txt",
		"null", "D:\\neptune\\data\\", "1", "1", "8", "25", "89", "0", "0", "0.0",
		"201", "0", "1001", "1", "0", "0", "22.0", "48.0", "0.0", "20.0", "0.3",
		"0.5", "500.0", "0.0", "0.0", "1.0", "1", "13", "0.0" };
	for (int i = 0; i < 33; i++)
		r.itemMap[itemNames[i]] = itemValues[i];

	// RUN & EXPECT (OUTPUTS)
	std::string item;
	for (int i = 0; i < 33; i++)
	{
    item.clear();
		r.getItem(itemNames[i], item);
    if (i != 4)
		  EXPECT_EQ(itemValues[i], item);
    else
      EXPECT_EQ("", item);
  }
	// TEST FOR ITEM NOT IN MAP
	std::string emptyItem;
	r.getItem("NOINMAP", emptyItem);
	EXPECT_EQ("", emptyItem);
	EXPECT_NE("NOTINMAP", emptyItem);

	// TEAR-DOWN
}
TEST(NamelistReader, getItem_int)
{
	// SETUP
	NamelistReaderExposed r;

	// GIVEN (INPUTS)
	std::string itemNames[15] = { "IERT", "IUTC", "MONTH", "MDAY", "MYEAR", "IHOUR", "IMIN",
		"NPOS", "LONEAST", "NR1", "NVARX", "NVARY", "LOGSCALE", "NMONTE", "IUP" };
	std::string itemValues[15] = { "1", "1", "8", "25", "89", "0", "0",
		"201", "0", "1001", "1", "0", "0","1", "13" };
	for (int i = 0; i < 15; i++)
		r.itemMap[itemNames[i]] = itemValues[i];

	// RUN & EXPECT (OUTPUTS)
	int item;
	for (int i = 0; i < 15; i++)
	{
		r.getItem(itemNames[i], item);
		EXPECT_EQ(std::stoi(itemValues[i]), item);
	}
	// TEST FOR ITEM NOT IN MAP
	int value = -1;
	r.getItem("NOINMAP", value);
	EXPECT_EQ(-1, value);

	// TEAR-DOWN
}

TEST(NamelistReader, getItem_double)
{
	// SETUP
	NamelistReaderExposed r;

	// GIVEN (INPUTS)
	std::string itemNames[12] = { "SEC", "FLAT", "FLON", "FHGT", "DELHGT",
		"DELLAT", "DELLON", "DELTIME", "PROFNEAR", "PROFFAR",
		"RPSCALE", "CORLMIN" };
	std::string itemValues[12] = { "0.0", "22.0", "48.0", "0.0", "20.0", "0.3",
		"0.5", "500.0", "0.0", "0.0", "1.0", "0.0" };
	for (int i = 0; i < 12; i++)
		r.itemMap[itemNames[i]] = itemValues[i];

	// RUN & EXPECT (OUTPUTS)
	double item;
	for (int i = 0; i < 12; i++)
	{
		r.getItem(itemNames[i], item);
		EXPECT_EQ(std::stod(itemValues[i]), item);
	}
	// TEST FOR ITEM NOT IN MAP
	double value = -1.0;
	r.getItem("NOINMAP", value);
	EXPECT_EQ(-1.0, value);

	// TEAR-DOWN
}

TEST(NamelistReader, getItem_float)
{
	// SETUP
	NamelistReaderExposed r;

	// GIVEN (INPUTS)
	std::string itemNames[12] = { "SEC", "FLAT", "FLON", "FHGT", "DELHGT",
		"DELLAT", "DELLON", "DELTIME", "PROFNEAR", "PROFFAR",
		"RPSCALE", "CORLMIN" };
	std::string itemValues[12] = { "0.0", "22.0", "48.0", "0.0", "20.0", "0.3",
		"0.5", "500.0", "0.0", "0.0", "1.0", "0.0" };
	for (int i = 0; i < 12; i++)
		r.itemMap[itemNames[i]] = itemValues[i];

	// RUN & EXPECT (OUTPUTS)
	float item;
	for (int i = 0; i < 12; i++)
	{
		r.getItem(itemNames[i], item);
		EXPECT_EQ(std::stof(itemValues[i]), item);
	}
	// TEST FOR ITEM NOT IN MAP
	float value = -1.0;
	r.getItem("NOINMAP", value);
	EXPECT_EQ(-1.0, value);

	// TEAR-DOWN
}

TEST(NamelistReader, getItem_bool)
{
	// SETUP
	NamelistReaderExposed r;

	// GIVEN (INPUTS)
	r.itemMap["LOGSCALE"] = "0";

	// RUN & EXPECT (OUTPUTS)
	bool value;
	r.getItem("LOGSCALE", value);
	EXPECT_EQ(false, value);
	EXPECT_NE(true, value);

	// TEST FOR ITEM NOT IN MAP
	value = true;
	r.getItem("NOINMAP", value);
	EXPECT_EQ(true, value);

	// TEAR-DOWN
}

TEST(NamelistReader, getParameters)
{
	// SETUP
	NamelistReaderExposed r;
	InputParameters p;

	// GIVEN (INPUTS)
	std::string itemNames[33] = { "SPICEDIR", "LSTFL", "OUTFL", "TRAJFL", "PROFILE",
		"DATADIR", "IERT", "IUTC", "MONTH", "MDAY", "MYEAR", "IHOUR", "IMIN", "SEC",
		"NPOS", "LONEAST", "NR1", "NVARX", "NVARY", "LOGSCALE", "FLAT", "FLON",
		"FHGT", "DELHGT", "DELLAT", "DELLON", "DELTIME", "PROFNEAR", "PROFFAR",
		"RPSCALE", "NMONTE", "IUP", "CORLMIN" };
	std::string itemValues[33] = { DEFAULT_SPICE_PATH, "LIST.txt", "OUTPUT.txt", "TRAJDATA.txt",
		"null", "D:\\neptune\\data\\", "1", "1", "8", "25", "89", "0", "0", "0.0",
		"201", "0", "1001", "1", "0", "0", "22.0", "48.0", "0.0", "20.0", "0.3",
		"0.5", "500.0", "0.0", "0.0", "1.0", "1", "13", "0.0" };
	for (int i = 0; i < 33; i++)
		r.itemMap[itemNames[i]] = itemValues[i];

	// RUN
	r.getParameters(p);
	
	// EXPECT (OUTPUTS)
	EXPECT_EQ(DEFAULT_SPICE_PATH, p.spicePath);
	EXPECT_EQ("LIST.txt", p.listFileName);
	EXPECT_EQ("OUTPUT.txt", p.columnFileName);
	EXPECT_EQ("TRAJDATA.txt", p.trajectoryFileName);
	EXPECT_EQ("", p.auxiliaryAtmosphereFileName[0]);
	EXPECT_EQ("D:/neptune/data/", p.dataPath);
	EXPECT_EQ(ERT, p.timeFrame);
	EXPECT_EQ(UTC, p.timeScale);
	EXPECT_EQ(8, p.month);
	EXPECT_EQ(25, p.day);
	EXPECT_EQ(1989, p.year);
	EXPECT_EQ(0, p.hour);
	EXPECT_EQ(0, p.minute);
	EXPECT_EQ(0.0, p.seconds);
	EXPECT_EQ(201, p.numberOfPositions);
  EXPECT_EQ(false, p.isEastLongitudePositiveOnInput);
  EXPECT_EQ(false, p.isEastLongitudePositiveOnOutput);
  EXPECT_EQ(1001, p.initialRandomSeed);
	EXPECT_EQ(1, p.nvarx);
	EXPECT_EQ(0, p.nvary);
	EXPECT_EQ(0, p.densityPrintScale);
	EXPECT_EQ(22.0, p.initialLatitude);
	EXPECT_EQ(48.0, p.initialLongitude);
	EXPECT_EQ(0.0, p.initialHeight);
	EXPECT_EQ(20.0, p.deltaHeight);
	EXPECT_EQ(0.3, p.deltaLatitude);
	EXPECT_EQ(0.5, p.deltaLongitude);
	EXPECT_EQ(500.0, p.deltaTime);
	EXPECT_EQ(0.0, p.innerRadius[0]);
	EXPECT_EQ(0.0, p.outerRadius[0]);
	EXPECT_EQ(1.0, p.densityPerturbationScale);
	EXPECT_EQ(1, p.numberOfMonteCarloRuns);
	EXPECT_EQ(13, p.iup);
	EXPECT_EQ(0.0, p.minRelativeStepSize);
	
	// TEAR-DOWN
}

TEST(NamelistReader, HASINPUT(readFile))
{
	// SETUP
	NamelistReaderExposed r;
	InputParameters p;
	// GIVEN (INPUTS)
	// need to put file in $(OutDir) defined in compiler
	std::string file = "namelist_test.txt";

	// RUN
	r.readFile(file);
	r.getParameters(p);

	// EXPECT (OUTPUTS)
	EXPECT_EQ(DEFAULT_SPICE_PATH, p.spicePath);
	EXPECT_EQ("LIST.txt", p.listFileName);
	EXPECT_EQ("OUTPUT.txt", p.columnFileName);
	EXPECT_EQ("TRAJDATA.txt", p.trajectoryFileName);
	EXPECT_EQ("null", p.auxiliaryAtmosphereFileName[0]);
	EXPECT_EQ("D:\\neptune\\data\\", p.dataPath);
	EXPECT_EQ(ERT, p.timeFrame);
	EXPECT_EQ(UTC, p.timeScale);
	EXPECT_EQ(8, p.month);
	EXPECT_EQ(25, p.day);
	EXPECT_EQ(1989, p.year);
	EXPECT_EQ(2, p.hour);
	EXPECT_EQ(45, p.minute);
	EXPECT_EQ(30.0, p.seconds);
	EXPECT_EQ(201, p.numberOfPositions);
  EXPECT_EQ(false, p.isEastLongitudePositiveOnInput);
  EXPECT_EQ(false, p.isEastLongitudePositiveOnOutput);
  EXPECT_EQ(1001, p.initialRandomSeed);
	EXPECT_EQ(1, p.nvarx);
	EXPECT_EQ(0, p.nvary);
	EXPECT_EQ(0, p.densityPrintScale);
	EXPECT_EQ(22.0, p.initialLatitude);
	EXPECT_EQ(48.0, p.initialLongitude);
	EXPECT_EQ(0.0, p.initialHeight);
	EXPECT_EQ(20.0, p.deltaHeight);
	EXPECT_EQ(0.3, p.deltaLatitude);
	EXPECT_EQ(0.5, p.deltaLongitude);
	EXPECT_EQ(500.0, p.deltaTime);
	EXPECT_EQ(0.0, p.innerRadius[0]);
	EXPECT_EQ(0.0, p.outerRadius[0]);
	EXPECT_EQ(1.0, p.densityPerturbationScale);
	EXPECT_EQ(1, p.numberOfMonteCarloRuns);
	EXPECT_EQ(13, p.iup);
	EXPECT_EQ(0.0, p.minRelativeStepSize);

	// TEAR-DOWN
}