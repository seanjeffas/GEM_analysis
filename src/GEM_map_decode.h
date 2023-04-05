#ifndef GEM_MAP_DECODE_H
#define GEM_MAP_DECODE_H

struct APVAddress{

  int crate;
  int mpd;
  int adc_id;
  
  // default constructor
APVAddress():
  crate(0), mpd(0), adc_id(0)
  {}

  // constructor
APVAddress(int cid, int mid, int aid) :
  crate(cid), mpd(mid), adc_id(aid)
  {}    

  bool operator==(const APVAddress &a) const {
        return (a.crate == crate) && 
            (a.mpd == mpd) && (a.adc_id == adc_id);
    }

    bool operator<(const APVAddress &a) const 
    {
        if(crate < a.crate) 
            return true;
        else if(crate == a.crate)
        {
            if(mpd < a.mpd) 
                return true;
            else if(mpd == a.mpd) 
            {
                if(adc_id < a.adc_id) 
                    return true;
                else 
                    return false;
            }
            else 
                return false;
        }
        else
            return false;
    }

    bool operator>(const APVAddress &a) const 
    {
        if(crate > a.crate) 
            return true;
        else if(crate == a.crate)
        {
            if(mpd > a.mpd) 
                return true;
            else if(mpd == a.mpd)
            {
                if(adc_id > a.adc_id) 
                    return true;
                else 
                    return false;
            }
            else 
                return false;
        }
        else
            return false;
    }
};

struct APV_info{

  int mpd;
  int adc_id;
  int flip;
  int module;
  TH1F *hAPV;

};


map<APVAddress, APV_info> APV;


void GEM_map_decode(TString filename);


#endif
