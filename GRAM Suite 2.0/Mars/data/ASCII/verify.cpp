// Program to convert ASCII files to binary values, write out
// data for map plots and/or profile plots, and write binary
// version NCEP data arrays

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char** argv)
{
  string fileName;
  cout << " Enter the file name: ";
  cin >> fileName;

  cout << " Reading binary file = " << fileName << endl;

  ifstream file(fileName, ios::binary);
  ifstream file2(fileName+"2", ios::binary);

  ofstream oFile(fileName + ".csv");
  oFile << scientific << setprecision(15);

  // Read the front block size tag
  size_t bytesFront = 0;
  size_t bytesFront2 = 0;
  double maxDiff = 0;
  int block = 1;
  while (!file.eof()) {
    file.read((char *)(&bytesFront), sizeof(size_t));
    file2.read((char *)(&bytesFront2), sizeof(size_t));
    
    if (!file.good()) break;

    size_t size = bytesFront / sizeof(double);
    
    cout << "Block " << block << " size: " << size << endl;

    double* blockOfDoubles = new double[size];
    double* blockOfDoubles2 = new double[size];

    // Read the data block of doubles.
    file.read((char *)(blockOfDoubles), bytesFront);
    file2.read((char *)(blockOfDoubles2), bytesFront2);

//    // Read the back block size tag
//    size_t bytesBack = 0;
//    file.read((char *)(&bytesBack), sizeof(int));
//    file2.read((char *)(&bytesBack), sizeof(int));

    for (int i = 0; i < size; ++i) {
      double relDiff = (blockOfDoubles[i] - blockOfDoubles2[i]) / (blockOfDoubles[i] + blockOfDoubles2[i]);
      if (relDiff > 1.0e-14) {
        oFile << block << ", " << i << ", " << blockOfDoubles[i] << ", " <<  blockOfDoubles2[i] << ", " << relDiff << '\n';
      }
      if (relDiff > maxDiff) {
        maxDiff = relDiff;
      }
    }
    ++block;
  }
  oFile.close();
  cout << "Max diff: " << maxDiff << endl;
  return 0;
}
