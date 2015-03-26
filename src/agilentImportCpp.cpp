#include <Rcpp.h>
#include <fstream>
using namespace Rcpp;
using namespace std;

#define HEADERSIZE 5758
#define MZSIZE 16
#define INTSIZE 16
#define SCANSIZE 16


int getFileSize(string filename){
  streampos begin,end;
  ifstream myfile (filename.c_str(), ios::binary);
  begin = myfile.tellg();
  myfile.seekg (0, ios::end);
  end = myfile.tellg();
  myfile.close();
  return(end - begin);
}

int16_t read16_be(istream& stream)
{
    uint8_t b[2];
    stream.read((char *)b,2);
    return static_cast<uint16_t>(
        (b[1])      |
        (b[0] << 8) );
}

unsigned long long readCertainBits(unsigned long long raw, int numOfBits, int readFrom, int readTo){
  unsigned long long result;
  int shift_number = numOfBits - readTo;
  result = (unsigned long long)pow(2, (readTo - readFrom + 1)) - 1;
	result = result << shift_number;
	result = raw & result;
	result = result >> shift_number;
	return result;
}

float encodeMZ(int16_t mz){
  return((float)((float)mz / 20));
}

int encodeIntensity(int16_t intensity){
  int base = readCertainBits(intensity, 16, 1, 2);
  int data = readCertainBits(intensity, 16, 3, 16);
  int encodedIntensity = data << (base * 3);
  return(encodedIntensity);
}

// [[Rcpp::export]]
List agilentImportCpp(string filename){
  List agilent;
  vector<int> fixedPattern;
  vector<int> counts;
  vector<float> mz;
  vector<int> intensity;
  
  ifstream fbin (filename.c_str(), ios::binary | ios::in);
  if (!fbin) {
    cout << "Could not open file " << filename;
    return agilent;
  }
  
  // Skip the header
  fbin.seekg (0, ios::beg);
  fbin.seekg(HEADERSIZE * sizeof(char));
  
  int count = -1;
  vector<float> scanMZ;
  vector<int> scanIntensity;
  while(!fbin.eof()){
    for(int i = 0; i < 11; i++)
      fixedPattern.push_back(read16_be(fbin));
    count = read16_be(fbin);
    if(count == 0)
      break;
      
    counts.push_back(count);
    // Highest peak in scan
    read16_be(fbin);
    read16_be(fbin);
    
    scanMZ.clear();
    scanIntensity.clear();
    // Get spectrometry data
    for(int j = 0; j < count; j++){
      scanMZ.push_back(encodeMZ(read16_be(fbin)));
      scanIntensity.push_back(encodeIntensity(read16_be(fbin)));
    }
    reverse(scanMZ.begin(), scanMZ.end());
    reverse(scanIntensity.begin(), scanIntensity.end());
    mz.insert(mz.end(), scanMZ.begin(), scanMZ.end());
    intensity.insert(intensity.end(), scanIntensity.begin(), scanIntensity.end());
  }
  cout << "\n";
  agilent["fixed_pattern"] = fixedPattern;
  agilent["counts"] = counts;
  agilent["mz"] = mz;
  agilent["intensity"] = intensity;
  
  fbin.close();
  return agilent;
}





