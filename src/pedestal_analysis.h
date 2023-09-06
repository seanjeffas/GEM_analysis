#ifndef PEDESTAL_ANALYSIS_H
#define PEDESTAL_ANALYSIS_H

#include "GEM_map_decode.h"
#include "Get_event_data.h"

double PedMean[nfiber][nadc][nchan] = {0};
double PedRMS[nfiber][nadc][nchan] = {0};
double CMmean[nfiber][nadc] = {0};
double CMrms[nfiber][nadc] = {0};
int APVMAP[nmodules][128];

void LoadPedestals( TString pedfilename );
void LoadCM( TString CMfilename );

void InitAPVMAP();


int sorting_strip_low = 54;
int sorting_strip_high = 54;
double fZeroSuppressRMS = 3;
int fCommonModeNumIterations = 3;
double fCommonModeRange_nsigma = 5.0;
double fCommonModeDanningMethod_NsigmaCut = 3.0;
double fRMS_ConversionFactor = sqrt(nsamples);
double fCommonModeBinWidth_Nsigma = 2.0; //Bin width +/- 2 sigma
double fCommonModeScanRange_Nsigma = 4.0; //Scan window +/- 4 sigma
double fCommonModeStepSize_Nsigma = 0.2; //sigma/5 for step size:
int fCommonModeMinStripsInRange = 10;

int enable_rolling_avg_danning = 0; //1/0 for yes/no
int enable_rolling_avg_histo = 1; //1/0 for yes/no
int fNeventsRollingAverage_by_APV[nfiber][nadc];
deque<double> fCommonModeResultContainer_by_APV[nfiber][nadc];
double fCommonModeRollingAverage_by_APV[nfiber][nadc];
double fCommonModeRollingRMS_by_APV[nfiber][nadc];
const int fNeventsCommonModeLookBack = 100;

double Sorting_CM(TH1F *hAPV, int isamp);
double Danning_CM_offline(APV_info APV_data, int isamp, int nsigma);
double Histogramming_CM(APV_info APV_data, int isamp, int use_rolling_avg, int nsigma, double setstepsize, double setbinwidth);
void UpdateRollingCommonModeAverage( APV_info APV_data, double CM_sample);

#endif
