#ifndef GET_EVENT_DATA_H
#define GET_EVENT_DATA_H


const int nmodules = 8;
const int nfiber = 24;
const int nadc = 16;
const int nchan = 128;
const int nsamples = 6;

//index 0 is x/u and 1 is y/v
const int nAPV[nmodules][2] = {
  {30,30},
  {30,30},
  {30,30},
  {30,30},
  {10,12},
  {10,12},
  {10,12},
  {10,12},
};


const  int maxstrips = 12000;
const  int maxdata = maxstrips*nsamples;

int nstrips[nmodules];
int ndata[nmodules];
double strip_num[nmodules][maxstrips], IsU[nmodules][maxstrips], ADC[nmodules][maxdata];
double mpd[nmodules][maxstrips], adc_id[nmodules][maxstrips];

TChain *C;
TH1F *hAPV[nmodules][2][30];
APV_info APV_data[nmodules][2][30]; //[nmodules][2 axes][30 APVs max on one axis]

void Get_event_data();
void tree_init(TString rootfile);

#endif
