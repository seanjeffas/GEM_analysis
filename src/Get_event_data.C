#include "Get_event_data.h"
#include "pedestal_analysis.h"



void Get_event_data(){


  // Loop over all modules
  for(int imod = 0; imod < nmodules; imod++){
    // Loop over all strip information for each module
    for(int idata=0; idata < ndata[imod]; idata++){
      
      //idata = isamp + nsamples * strip_index
      int isamp = idata % nsamples;    
      int index = (idata - isamp) / nsamples;
      int axis = IsU[imod][index] == 1 ? 0 : 1;
      int iAPV = strip_num[imod][index] / 128;
      int istrip = (int)strip_num[imod][index] % 128;

      // Subtract the pedestal offset from the signal
      //double adc_signal = ADC[imod][idata] - PedMean[(int)mpd[imod][index]][(int)adc_id[imod][index]][istrip];
      double adc_signal = ADC[imod][idata] - PedMean[(int)mpd[imod][index]][(int)adc_id[imod][index]][istrip];

	
      //Set the histogram bin content to the ADC value
      //We need to use one bin to divide the time samples so they are divided by 129 instead of 128 as seen below
      hAPV[imod][axis][iAPV]->SetBinContent(istrip + 129*isamp,adc_signal);
      
      //Get the APV mpd and adc_ch. This is needed to calculate the common mode later
      APV_data[imod][axis][iAPV].mpd = (int)mpd[imod][index];
      APV_data[imod][axis][iAPV].adc_id = (int)adc_id[imod][index];
    } // end loop on strip data
  } // end loop on modules


  for(int imod = 0; imod < nmodules; imod++){
    for(int iaxis=0; iaxis < 2; iaxis++){
      for(int iAPV=0; iAPV < nAPV[imod][iaxis]; iAPV++){
	
	//Set divisions between time samples
	for(int isamp = 0; isamp < nsamples; isamp++)
	  hAPV[imod][iaxis][iAPV]->SetBinContent(128*isamp,0);
	
	APV_data[imod][iaxis][iAPV].hAPV = hAPV[imod][iaxis][iAPV];
      }
    }
  }


}

void tree_init(TString rootfile){

  if(C != NULL) C->Reset();

  C = new TChain("T");
  C->Add(rootfile);

  C->SetBranchStatus("*",0);
  
  for(int imod=0; imod < nmodules; imod++){
    
    C->SetBranchStatus(Form("bb.gem.m%i.strip.istrip",imod),1);
    C->SetBranchStatus(Form("bb.gem.m%i.strip.IsU",imod),1);
    C->SetBranchStatus(Form("bb.gem.m%i.strip.rawADCsamples",imod),1);
    C->SetBranchStatus(Form("bb.gem.m%i.strip.mpd",imod),1);
    C->SetBranchStatus(Form("bb.gem.m%i.strip.adc_id",imod),1);
    C->SetBranchStatus(Form("Ndata.bb.gem.m%i.strip.istrip",imod),1);
    C->SetBranchStatus(Form("Ndata.bb.gem.m%i.strip.rawADCsamples",imod),1);

    C->SetBranchAddress(Form("bb.gem.m%i.strip.istrip",imod),strip_num[imod]);
    C->SetBranchAddress(Form("bb.gem.m%i.strip.IsU",imod),IsU[imod]);
    C->SetBranchAddress(Form("bb.gem.m%i.strip.rawADCsamples",imod),ADC[imod]);
    C->SetBranchAddress(Form("bb.gem.m%i.strip.mpd",imod),mpd[imod]);
    C->SetBranchAddress(Form("bb.gem.m%i.strip.adc_id",imod),adc_id[imod]);
    C->SetBranchAddress(Form("Ndata.bb.gem.m%i.strip.istrip",imod),&nstrips[imod]);
    C->SetBranchAddress(Form("Ndata.bb.gem.m%i.strip.rawADCsamples",imod),&ndata[imod]);

    //Initialize the APV histograms
    for(int iaxis=0; iaxis < 2; iaxis++){
      for(int iAPV=0; iAPV < nAPV[imod][iaxis]; iAPV++){
	hAPV[imod][iaxis][iAPV] = new TH1F(Form("hAPV_%i_%i_%i",imod,iaxis,iAPV),Form("Module %i Axis %i APV %i;Strip;ADC",imod,iaxis,iAPV),768,0,768);
      }
    }

  }


}
