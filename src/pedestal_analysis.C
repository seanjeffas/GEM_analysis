#include "pedestal_analysis.h"
#include "GEM_map_decode.h"



void LoadPedestals( TString pedfilename ){

  
  std::ifstream pedfile( pedfilename.Data() );

  if( !pedfile.good() ){
    pedfile.close();

    std::cout << "Warning: could not find ped file " << pedfilename << " in working directory, pedestals not loaded" << std::endl;

    return;
  } else {
    std::cout << "Found pedestal file " << pedfilename << endl;
  }

  std::string currentline;

  int crate=0, slot=0, mpd=0, adc_ch=0;
  APVAddress apv_addr;
  APV_info APV_data;

  //Read each line of the pedestal file
  while( std::getline(pedfile, currentline) ){

    if( pedfile.eof() ) break;

    if( currentline[0] != '#' ){    
      
      std::istringstream is(currentline);

      string dummy;
      //Read the line that tells us which APV this is
      if ( currentline.find("APV") == 0 ){
	is >> dummy >> crate >> slot >> mpd >> adc_ch;
	apv_addr.crate = crate;
	apv_addr.mpd = mpd;
	apv_addr.adc_id = adc_ch;
	APV_data = APV.find(apv_addr)->second;
	
      } else {
	//Otherwise the lines should tell us the APV channel mean and rms
	int apvchan;
	double mean, rms;
	is >> apvchan >> mean >> rms;

	//Convert apvchan to strip chan
	//APVMAP must be initialized for this to work
	if(APV_data.flip == 1){
	  PedMean[mpd][adc_ch][127 - APVMAP[APV_data.module][apvchan]] = mean;
	  PedRMS[mpd][adc_ch][127 - APVMAP[APV_data.module][apvchan]] = rms;
	}
	else{
	  PedMean[mpd][adc_ch][APVMAP[APV_data.module][apvchan]] = mean;
	  PedRMS[mpd][adc_ch][APVMAP[APV_data.module][apvchan]] = rms;
	}

      }
    }
  }

}



void LoadCM( TString CMfilename ){

  
  std::ifstream CMfile( CMfilename.Data() );

  if( !CMfile.good() ){
    CMfile.close();

    std::cout << "Warning: could not find CM file " << CMfilename << " in working directory, pedestals not loaded" << std::endl;

    return;
  } else {
    std::cout << "Found CM file " << CMfilename << endl;
  }


  
  std::string currentline;

  int crate=0, slot=0, mpd=0, adc_ch=0;
  
  //Read each line of the file
  while( std::getline(CMfile, currentline) ){
    //TString currentline;
    if( CMfile.eof() ) break;

    if( currentline[0] != '#' ){    
      
      std::istringstream is(currentline);

      string dummy;
      
      double mean, rms;
   
      //Read the APV infor and the CM mean and rms
      is >> crate >> slot >> mpd >> adc_ch >> mean >> rms;
       
      CMmean[mpd][adc_ch] = mean;
      CMrms[mpd][adc_ch] = rms;
    }
  }

}

//This is how the APV internal channels to strips
void InitAPVMAP(){
  
  for(int imod = 0; imod < nmodules; imod++){
    for( UInt_t i=0; i<128; i++ ){
      Int_t strip1 = 32*(i%4) + 8*(i/4) - 31*(i/16);
      Int_t strip2 = strip1 + 1 + strip1 % 4 - 5 * ( ( strip1/4 ) % 2 );
      Int_t strip3 = ( strip2 % 2 == 0 ) ? strip2/2 + 32 : ( (strip2<64) ? (63 - strip2)/2 : 127 + (65-strip2)/2 );
      if(imod > 3) //This should not be hard coded. Fix this later
	APVMAP[imod][i] = strip2;
      else 
	APVMAP[imod][i] = strip3;
    }
  }
}

// Calculate the common mode using the sorting method, here is the outline of the method:
// 1. Sort all strips in one time sample from lowest to highest ADC
// 2. Remove some number (user defined) of low and high strips
// 3. Take the average of the remaining strips and this is the CM
//In this particular code the APV data is saved in a histogram. So the inputs are the APV histogram and the time sample in question.
double Sorting_CM(TH1F *hAPV, int isamp){

  vector<double> adc;

  //Read all strips and ADCs for this time sample
  for(int istrip = 0; istrip < 128; istrip++){
    adc.push_back(hAPV->GetBinContent(istrip + 129*isamp));
  }
  //Sort all the strips by ADC value
  std::sort( adc.begin(), adc.end() );

  double CM = 0;
  int n_keep = 0;

  //Remove all the low and high strips
  for(int ihit = sorting_strip_low; ihit < 128 - sorting_strip_high; ihit++){
    CM += adc[ihit];
    n_keep++;
  }

  // Return the average
  return CM/n_keep;

}


//Calculate the common mode using the Danning method, here is the method outline:
// 1. Read in the CM mean and RMS from the pedestal run
// 2. Cut out all strips outside of CM_Mean +/- n*CM_sigma (n is user defined)
// 3. Take the average of the remaining strips
// 4. Take the new average and cut all strips outside of new_average +/- n*ped_rms. Here ped_rms is the pedestal rms for that partcular strip in question
// 5. Repeat strip 4 again with the new average. Overall the averaging is done 3 times, but the number of iterations is user defined.
// The histogram input has all the APV signal data. The mpd and adc_ch input is used to find which CM_mean and CM_sigma to start with
double Danning_CM_offline(APV_info APV_data, int isamp, int nsigma = 0){

  TH1F *hAPV = APV_data.hAPV;
  int mpd = APV_data.mpd;
  int adc_ch = APV_data.adc_id;

  double cm_mean = CMmean[mpd][adc_ch];
  double cm_rms = CMrms[mpd][adc_ch];

  int fNsigma = fCommonModeRange_nsigma;
  if(nsigma != 0) fNsigma = nsigma;

  if( fNeventsRollingAverage_by_APV[mpd][adc_ch] >= std::min(100, nsamples*fNeventsCommonModeLookBack ) && enable_rolling_avg_danning){
    cm_mean = fCommonModeRollingAverage_by_APV[mpd][adc_ch];
    cm_rms = fCommonModeRollingRMS_by_APV[mpd][adc_ch];
  }
  
  double cm_temp = 0.0;
  int n_keep_final = 0;

  for( int iter=0; iter<3; iter++ ){

    double cm_min = cm_mean - fNsigma*cm_rms;
    double cm_max = cm_mean + fNsigma*cm_rms;
    double sumADCinrange = 0.0;
    int n_keep = 0;
      
    for( int ihit=0; ihit<128; ihit++ ){
	
      double ADCtemp = hAPV->GetBinContent(ihit + 129*isamp);
      double rmstemp = PedRMS[mpd][adc_ch][ihit];
	
      if(iter != 0){
	cm_min = cm_temp - fCommonModeDanningMethod_NsigmaCut*2.5*rmstemp;
	cm_max = cm_temp + fCommonModeDanningMethod_NsigmaCut*2.5*rmstemp;
      }

      if( ADCtemp >= cm_min && ADCtemp <= cm_max ){
	n_keep++;
	sumADCinrange += ADCtemp;

      }
    }
   
    cm_temp = sumADCinrange / n_keep;
    n_keep_final = n_keep;
  }

  if(n_keep_final == 0)
    return Sorting_CM( hAPV, isamp);
  else
    return cm_temp;
}



//Add description for this method later. It is complicated
double Histogramming_CM(APV_info APV_data, int isamp, int use_rolling_avg = 1, int nsigma = 0, double setstepsize = 0, double setbinwidth = 0){
 
  TH1F *hAPV = APV_data.hAPV;
  int mpd = APV_data.mpd;
  int adc_ch = APV_data.adc_id;

  double cm_mean = CMmean[mpd][adc_ch];
  double cm_rms = CMrms[mpd][adc_ch];

  int fNsigma = fCommonModeScanRange_Nsigma;
  if(nsigma != 0) fNsigma = nsigma;
  
  if( fNeventsRollingAverage_by_APV[mpd][adc_ch] >= std::min(100, nsamples*fNeventsCommonModeLookBack ) && enable_rolling_avg_histo && use_rolling_avg){
    cm_mean = fCommonModeRollingAverage_by_APV[mpd][adc_ch];
    cm_rms = fCommonModeRollingRMS_by_APV[mpd][adc_ch];
  }

  //bin width/stepsize = 8 with these settings:
  double stepsize = cm_rms*fCommonModeStepSize_Nsigma; //Default is 0.2 = rms/5
  double binwidth = cm_rms*fCommonModeBinWidth_Nsigma; //Default is +/- 2 sigma, bin width / step size = 20 with these settings

  //this will actually include all ADCs within +/- (ScanRange + BinWidth) sigma of the mean since range is bin center +/- 1*RMS.
  double scan_min = cm_mean - fNsigma*cm_rms; 
  double scan_max = cm_mean + fNsigma*cm_rms;
  
  if(setstepsize != 0)
    stepsize = setstepsize;

  if(setbinwidth != 0)
    binwidth = setbinwidth;

  
  /*
  if(nsigma == 0 && use_rolling_avg == 1 && APV_data.iAPV == 2 && isamp == 0)
    cout<<"correct "<<scan_min<<" "<<scan_max<<" "<<stepsize<<" "<<binwidth<<endl;

  if(nsigma != 0 && APV_data.iAPV == 2 && isamp == 0)
    cout<<"online "<<scan_min<<" "<<scan_max<<" "<<stepsize<<" "<<binwidth<<endl;
  */
  int nbins= int( (scan_max - scan_min)/stepsize ); //Default = 8 * RMS / (rms/5) = 40 bins.
  
  //NOTE: The largest number of bins that could contain any given sample is binwidth/stepsize = 20 with default settings:
  if(isnan(stepsize) || stepsize == 0){// APV data looks weird, default to sorting
    return Sorting_CM( hAPV, isamp);
  }

  //if(stepsize == 0) return GetCommonMode( isamp, 0, apvinfo );
  //Construct std::vectors and explicitly zero-initialize them:
  std::vector<double> bincounts(nbins,0);
  std::vector<double> binADCsum(nbins,0.0);
  std::vector<double> binADCsum2(nbins,0.0);
 
  int ibinmax=-1;
  int maxcounts=0;
  //Now loop on all the strips and fill the histogram: this is the version assuming full readout:
  //for( int istrip=0; istrip<fN_APV25_CHAN; istrip++ ){
  for( int istrip=0; istrip<128; istrip++ ){
      
    double ADC =  hAPV->GetBinContent(istrip + 129*isamp);
    //double ADC = strip_ADC[istrip][isamp];
      
    //calculate the lowest bin containing this ADC value. 
    int nearestbin = std::max(0,std::min(nbins-1, int(round( (ADC - scan_min)/stepsize ) ) ) );

    int binlow = nearestbin;
    int binhigh = nearestbin+1;
      
    while( binlow >= 0 && fabs( ADC - (scan_min + binlow*stepsize) ) <= binwidth ){
      bincounts[binlow]++;
      binADCsum[binlow] += ADC;
      binADCsum2[binlow] += pow(ADC,2);

      if( ibinmax < 0 || bincounts[binlow] > maxcounts ){
	ibinmax = binlow;
	maxcounts = bincounts[binlow];
      }
      binlow--;
    }
      
    while( binhigh < nbins && fabs( ADC - (scan_min + binhigh*stepsize) ) <= binwidth ){
      bincounts[binhigh]++;
      binADCsum[binhigh] += ADC;
      binADCsum2[binhigh] += pow(ADC,2);
      if( ibinmax < 0 || bincounts[binhigh] > maxcounts ){
	ibinmax = binhigh;
	maxcounts = bincounts[binhigh];
      }
      binhigh++;
    }
      
  }
  /*
  if(nsigma != 0 && APV_data.iAPV == 1 && isamp == 0) {
    for(int ibin=0; ibin < nbins;ibin++){
      cout<<"Online: "<<binADCsum[ibin]/double(bincounts[ibin])<<" "<<bincounts[ibin]<<endl;
    }      
    cout<<"Online result: "<<binADCsum[ibinmax]/double(bincounts[ibinmax])<<endl;
  }
  if(nsigma == 0 && APV_data.iAPV == 1 && isamp == 0 && use_rolling_avg == 1) {
    for(int ibin=0; ibin < nbins;ibin++){
      cout<<"Offline: "<<binADCsum[ibin]/double(bincounts[ibin])<<" "<<bincounts[ibin]<<endl;
    }  
    cout<<"Offline result: "<<binADCsum[ibinmax]/double(bincounts[ibinmax])<<endl; 
  }
  */
  if( ibinmax >= 0 && maxcounts >= fCommonModeMinStripsInRange ){
    return binADCsum[ibinmax]/double(bincounts[ibinmax]);
  } else { //Fall back on sorting method:
    return Sorting_CM( hAPV, isamp);
  }
}


//This gets called for each time sample for each APV or each full readout event or whenever BUILD_ALL_SAMPLES is true and CM_ENABLED is false:
//There are two cases to handle:
// 1) before the container is full, meaning less than nsamples*fNeventsCommonModeLookBack have been added. In this case we increment everything.
// 2) after the container is full, meaning the earliest sample needs to roll off the average and one new sample has to be added at the end.
void UpdateRollingCommonModeAverage( APV_info APV_data, double CM_sample ){
  
 
  int mpd = APV_data.mpd;
  int iapv = APV_data.adc_id;
  
  UInt_t N = fNeventsRollingAverage_by_APV[mpd][iapv];
  UInt_t Nmax = nsamples*fNeventsCommonModeLookBack;
  
  double sum, sum2;
  
  if( N < Nmax ){
    //before reaching the size of the look back window, we just add the common-mode samples onto the end of the array

    fCommonModeResultContainer_by_APV[mpd][iapv].push_back(CM_sample);
    
    if( N == 0 ){ //First sample, initialize all sums/averages:
      fCommonModeRollingAverage_by_APV[mpd][iapv] = CM_sample;
      fCommonModeRollingRMS_by_APV[mpd][iapv] = 0.0;
      sum = CM_sample;
      sum2 = pow(CM_sample,2);
    } else { //Second and subsequent samples: increment sums, recalculate
      double oldavg = fCommonModeRollingAverage_by_APV[mpd][iapv];
      double oldrms = fCommonModeRollingRMS_by_APV[mpd][iapv];
      
      sum = N*oldavg + CM_sample;
      sum2 = N * (pow(oldrms,2) + pow(oldavg,2)) + pow(CM_sample,2);

      double newavg = sum/double(N+1);
      double newrms = sqrt( sum2/double(N+1) - pow(newavg,2) );

      fCommonModeRollingAverage_by_APV[mpd][iapv] = newavg;
      fCommonModeRollingRMS_by_APV[mpd][iapv] = newrms;

    }

    fNeventsRollingAverage_by_APV[mpd][iapv] = N+1;
  
  } else {
      
    //grab the earliest sample in the rolling average:
    double oldfirstsample = fCommonModeResultContainer_by_APV[mpd][iapv].front();
    
    //The net result of the following two operations should be to keep the container size the same:
    fCommonModeResultContainer_by_APV[mpd][iapv].pop_front(); //remove oldest sample
    fCommonModeResultContainer_by_APV[mpd][iapv].push_back( CM_sample ); //Insert newest sample at the end
    //we only need to update the calculation for the fact that the
    //earliest sample rolled off and a new sample was added: 
    double oldavg = fCommonModeRollingAverage_by_APV[mpd][iapv];
    double oldsum = oldavg * Nmax;
    
    double oldrms = fCommonModeRollingRMS_by_APV[mpd][iapv];
    // RMS^2 = sum^2/N - avg^2 --> sum^2 = N * (RMS^2 + avg^2)
    double oldsum2 = Nmax * ( pow(oldrms,2) + pow(oldavg,2) );

    //double lastsample = fCommonModeResultContainer_by_APV[mpd][iapv].back();
    double lastsample = CM_sample;
    
    double newsum = oldsum - oldfirstsample + lastsample;
    double newsum2 = oldsum2 - pow(oldfirstsample,2) + pow(lastsample,2);

    double newavg = newsum/double( Nmax );
    double newrms = sqrt( newsum2/double( Nmax ) - pow(newavg,2) );
    
    
    fCommonModeRollingAverage_by_APV[mpd][iapv] = newavg;
    fCommonModeRollingRMS_by_APV[mpd][iapv] = newrms;
    
  }
  
}
