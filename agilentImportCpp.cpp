#include <Rcpp.h>
#include <fstream>
using namespace Rcpp;
using namespace std;

#define HEADERSIZE 1440 * 4
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

// [[Rcpp::export]]
List agilentImportCpp(string filename){
  List agilent;
  vector<int> fixedPattern;
  
  ifstream fbin (filename.c_str(), ios::binary | ios::in);
  if (!fbin) {
    cout << "Could not open file " << filename;
    return agilent;
  }
  
  // Skip the header
  fbin.seekg (0, ios::beg);
  fbin.seekg(HEADERSIZE * sizeof(char));
  
  int i = 0;
  while(!fbin.eof()){
    if(i > 10)
      break;
    fixedPattern.push_back(i);
    cout << fixedPattern.size() << " " << read16_be(fbin) << "\n";
    i++;
  }
  cout << "\n";
  agilent["fixed_pattern"] = fixedPattern;
  
  fbin.close();
  return agilent;
}







