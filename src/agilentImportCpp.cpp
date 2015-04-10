#include <Rcpp.h>
#include <fstream>
using namespace Rcpp;
using namespace std;

#define HEADERSIZE 5758
#define MZSIZE 16
#define INTSIZE 16
#define SCANSIZE 16
#define MZBASE 20
#define INTENSITYHEADER 2
#define SHIFTBASE 3
#define MILISECOND 1000
#define MINMZ 100000
#define MAXMZ -1

double round(float d) {
  return floor(d + 0.5);
}

IntegerMatrix assembleMassSpectrometry(List importedData){
  vector<float> scanTime = importedData["scanTime"];
  int scanNum = scanTime.size();
  vector<double> uniqueMz = importedData["uniqueMz"];
  int mzNum = uniqueMz.size();
  IntegerMatrix massSpectrometry(scanNum, mzNum);
  int count;
  int prevCount = 0;
  vector<int> counts = importedData["counts"];
  vector<float> mzs = importedData["mz"];
  vector<float> intensity = importedData["intensity"];
  vector<double>::iterator it;
  double mz;
  int mzPos;
  
  for(int i = 0; i < scanNum; i++){
    count = counts[i];
    if(i != 0)
      prevCount = prevCount + counts[i-1];
    for(int j = 0; j < count; j++){
      mz = round(mzs[prevCount + j]);
      it = find(uniqueMz.begin(), uniqueMz.end(), mz);
      mzPos = distance(uniqueMz.begin(), it);
      massSpectrometry(i, mzPos) = massSpectrometry(i, mzPos) + intensity[prevCount + j];
    }
  }
  
  return(massSpectrometry);
}

int getFileSize(string filename){
  streampos begin,end;
  ifstream myfile (filename.c_str(), ios::binary);
  begin = myfile.tellg();
  myfile.seekg (0, ios::end);
  end = myfile.tellg();
  myfile.close();
  return(end - begin);
}

uint16_t read16_be(istream& stream)
{
    uint8_t b[2];
    stream.read((char *)b,2);
    return static_cast<uint16_t>(
        (b[1])      |
        (b[0] << 8) );
}

unsigned long readCertainBits(unsigned long raw, int numOfBits, int readFrom, int readTo){
  unsigned long result;
  int shift_number = numOfBits - readTo;
  result = (unsigned long)pow(2, (readTo - readFrom + 1)) - 1;
	result = result << shift_number;
	result = raw & result;
	result = result >> shift_number;
	return result;
}

float encodeMZ(int16_t mz){
  return((float)((float)mz / MZBASE));
}

int encodeIntensity(int16_t intensity){
  int intensityHeader = readCertainBits(intensity, 16, 1, INTENSITYHEADER);
  int intensityData = readCertainBits(intensity, 16, INTENSITYHEADER + 1, 16);
  int encodedIntensity = intensityData << (intensityHeader * SHIFTBASE);
  return(encodedIntensity);
}

//' @importFrom Rcpp evalCpp
//' @useDynLib 'GCMS'
//' @export
// [[Rcpp::export]]
List agilentImportCpp(std::string file, bool assembleMS = true){
  List agilent;
  vector<int> fixedPatterns;
  vector<int> counts;
  vector<float> mz;
  vector<double> uniqueMz;
  vector<int> intensity;
  vector<float> scanTime;
  float maxMz = MAXMZ;
  float minMz = MINMZ;
  float currMz;
//  std::string file = Rcpp::as<std::string>(filename);

  int found = file.find_first_of("~");
  if(found != -1)
    file = file.replace(0,1,getenv("HOME"));
  
  ifstream fbin (file.c_str(), ios::binary | ios::in);
  if (!fbin) {
    cout << "Could not open file " << file;
    return agilent;
  }
  
  // Skip the header
  fbin.seekg (0, ios::beg);
  fbin.seekg(HEADERSIZE * sizeof(char));
  
  unsigned int count, timePrefix, timeSuffix;
  vector<float> scanMZ;
  vector<int> scanIntensity;
  vector<unsigned int> blockPattern;
//  bool breakFlag = false;
  while(!fbin.eof()){
    blockPattern.clear();
    for(int i = 0; i < 12; i++){
      blockPattern.push_back(read16_be(fbin));
      if(i == 6)
        timePrefix = blockPattern.back();
      else if(i == 7)
        timeSuffix = blockPattern.back();
      else if(i == 11)
        count = blockPattern.back();
    }
    if(blockPattern.at(9) != 1)
      break;
    else
      fixedPatterns.insert(fixedPatterns.end(), blockPattern.begin(), blockPattern.end());
    
    scanTime.push_back(((float)(timePrefix << 16) + timeSuffix) / MILISECOND);
    counts.push_back(count);
    // Highest peak in scan
    read16_be(fbin);
    read16_be(fbin);
    
    scanMZ.clear();
    scanIntensity.clear();
    // Get spectrometry data
    for(int j = 0; j < count; j++){
      currMz = encodeMZ(read16_be(fbin));
      
      if(currMz > maxMz)
        maxMz = currMz;
      if(currMz < minMz)
        minMz = currMz;
        
      if(find(uniqueMz.begin(), uniqueMz.end(), round(currMz)) == uniqueMz.end()) {
          /* v does not contain x */
          uniqueMz.push_back(round(currMz));
      } else {
          /* v contains x */
      }
        
      scanMZ.push_back(currMz);
      scanIntensity.push_back(encodeIntensity(read16_be(fbin)));
    }
    reverse(scanMZ.begin(), scanMZ.end());
    reverse(scanIntensity.begin(), scanIntensity.end());
    mz.insert(mz.end(), scanMZ.begin(), scanMZ.end());
    intensity.insert(intensity.end(), scanIntensity.begin(), scanIntensity.end());
  }
  fbin.close();
  
  agilent["fixed_patterns"] = fixedPatterns;
  agilent["counts"] = counts;
  agilent["mz"] = mz;
  agilent["intensity"] = intensity;
  agilent["scanTime"] = scanTime;
  agilent["minMz"] = minMz;
  agilent["maxMz"] = maxMz;
  if(assembleMS){
    sort(uniqueMz.begin(), uniqueMz.end(), std::less<float>());
    agilent["uniqueMz"] = uniqueMz;
    agilent["massSpectrometry"] = assembleMassSpectrometry(agilent);
  }
    
  return agilent;
}





