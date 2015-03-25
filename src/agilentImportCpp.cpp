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

float encodeMZ(int16_t mz){
  return((float)(mz / 20));
}

int encodeIntensity(int16_t intensity){
  int shiftNumber = (intensity >> 13) * 3;
  return(intensity << shiftNumber);
}

// [[Rcpp::export]]
List agilentImportCpp(string filename){
  List agilent;
  vector<int> fixedPattern;
  vector<int> counts;
  vector<int> mz;
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
  int scanIndex = 0;
  vector<int> scanMZ;
  vector<int> scanIntensity;
  while(!fbin.eof()){
    for(int i = 0; i < 11; i++)
      fixedPattern.push_back(read16_be(fbin));
    count = read16_be(fbin);
    counts.push_back(count);
    // Highest peak in scan
    read16_be(fbin);
    read16_be(fbin);
    
    // Get spectrometry data
    for(int j = 0; j < count; j++){
      scanMZ.push_back(encodeMZ(read16_be(fbin)));
      scanIntensity.push_back(encodeIntensity(read16_be(fbin)));
    }
    reverse(scanMZ.begin(), scanMZ.end());
    mz.push_back(scanMZ);
    reverse(scanIntensity.begin(), scanIntensity.end());
    intensity.push_back(scanIntensity);
    scanIndex++;
    if(scanIndex > 8000)
      break;
  }
  cout << "\n";
  agilent["fixed_pattern"] = fixedPattern;
  agilent["counts"] = counts;
  agilent["mz"] = mz;
  agilent["intensity"] = intensity;
  
  fbin.close();
  return agilent;
}







