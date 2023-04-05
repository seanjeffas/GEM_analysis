#include "GEM_map_decode.h"



void GEM_map_decode(TString filename){
 
  ifstream infile(filename.Data());

  string line;
  
  while(getline(infile,line)){
  
    if( line[0] == '#' || line.length() == 0)
      continue;
    
    vector<int> values;
    std::stringstream ss(line);
    std::string token;
    int itoken = 0;
    while (std::getline(ss, token, ',')) {
      if(itoken == 0 || itoken ==  10){
	itoken++;
	continue;
      }
      values.push_back(std::stoi(token));
      itoken++;
    }
    
    APVAddress apv_addr;
    apv_addr.crate = values[0];
    apv_addr.mpd = values[2];
    apv_addr.adc_id = values[5];

    APV_info temp_APV;
    temp_APV.module = values[1] + values[10]; //Layer + gem position
    temp_APV.flip = values[8];
   
    APV.insert({apv_addr,temp_APV});
  }

}
