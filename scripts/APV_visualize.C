////////////////////////////////////////////////////////////////
//  Written by Sean Jeffas
//  sj9ry@virginia.edu
//  Last updated March 23, 2023
//
//  The purpose of this script is to take full readout GEM data
//  generated from SBS-offline and visually display APV data.
//  It also shows the common mode results visually for each APV.
//  The script only needs a run number and a pedestal run number
//  And it will prompt the user to select the APVs they want to
//  look at.
//
//  Usage: root APV_visualize(#runnumber, #pedestal_run_number)
//
////////////////////////////////////////////////////////////////
#include "../src/pedestal_analysis.C"
#include "../src/GEM_map_decode.C"
#include "../src/Get_event_data.C"





void APV_visualize(int runnumber = 2820,int ped_run_number = 2360){

  gStyle->SetOptStat(0);
  
  int lookmod = 0;
  int lookaxis = 0;
  bool first15 = true;
  bool printpdf = false;

  // Prompt the user to find out how which APVs to look at
  TString reply = "no";
  cout << "Follow the prompts for which range of APVs you want to investigate:" << endl;
  cout<<"Module number: ";
  cin >> lookmod;
  cout<<"Axis (0/1 with 0 as X/U and 1 and Y/V): ";
  cin >> lookaxis;
  cout<<"Do you want to look at the first 15 APVs here (y/n)? If answered no then the second 15 will be plotted"<<endl;
  cin >> reply;
  first15 = reply.BeginsWith("y") ? true : false;
  cout<<"Print PDF of 100 events (y/n)?"<<endl;
  cin >> reply;  
  printpdf = reply.BeginsWith("y") ? true : false;

  if(printpdf) gROOT->SetBatch(kTRUE);
 
  
  TString pdffile = Form("../plots/Event_display_run%i.pdf",runnumber);
  int npages = 100;
  int ipage = 0;

  // Varius root files and database files
  TString rootfile = Form("/volatile/halla/sbs/jeffas/GEN_root/Rootfiles/GEM_luminosity/full_readout/e1209016_replayed_%i_stream0_2_*_nevent5000*.root",runnumber);
  TString pedfile = Form("/w/halla-scshelf2102/sbs/jeffas/SBS_OFFLINE/SBS-replay/DB/gemped/daq_ped_bb_gem_run%i.dat",ped_run_number);
  TString cmfile = Form("/w/halla-scshelf2102/sbs/jeffas/SBS_OFFLINE/SBS-replay/DB/gemped/db_cmr_bb_gem_run%i.dat",ped_run_number);

  //Get our GEM map and APV map and load pedestals
  GEM_map_decode("../src/gem_map_BigBite.txt");
  InitAPVMAP();
  LoadPedestals(pedfile);
  LoadCM(cmfile);
  
  tree_init(rootfile);
 
  
  TLine *l_Sorting[nsamples];
  TLine *l_Danning[nsamples];
  TLine *l_Danning_nsigma[nsamples];
  TLine *l_Histogramming[nsamples];
  TLine *l_Histogramming_online[nsamples];
  TLine *l_Histogramming_online_nsigma[nsamples];
  TCanvas *c = new TCanvas("c","",1600,1000);
  c->Divide(4,4);

  int nevent = 0;
  
  //Start looping over all events
  while(C->GetEntry(nevent++)){
   
    Get_event_data(); //Function puts all the APV data into histograms
    
    for(int imod = 0; imod < nmodules; imod++){
      for(int iaxis = 0; iaxis < 2; iaxis++){
	for(int iAPV = 0; iAPV < nAPV[imod][iaxis]; iAPV++){
	  for(int isamp = 0; isamp < nsamples; isamp++){
	    UpdateRollingCommonModeAverage(APV_data[imod][iaxis][iAPV],Histogramming_CM(APV_data[imod][iaxis][iAPV], isamp));
	    //UpdateRollingCommonModeAverage(APV_data[imod][iaxis][iAPV],Danning_CM_offline(APV_data[imod][iaxis][iAPV], isamp));
	  }
	}
      }
    }
    

    if(nevent > 100){
      int APV_offset = 0;
      
      //We can only reasonbly display 15 plots at once
      if(!first15 && nAPV[lookmod][lookaxis] > 15) APV_offset = 15;

      int nAPVtot = nAPV[lookmod][lookaxis];
      if(nAPV[lookmod][lookaxis] > 15) nAPVtot = 15;

      for(int iAPV=0; iAPV < nAPVtot; iAPV++){

	int lookAPV = iAPV + APV_offset;

	//Get APV data
	APV_info thisAPV = APV_data[lookmod][lookaxis][lookAPV];

	//Loop over all time samples
	
	for(int isamp = 0; isamp < nsamples; isamp++){

	  //Calculat the sorting CM and plot it
	  double CM_sorting = Sorting_CM(thisAPV.hAPV,isamp);
	  l_Sorting[isamp] = new TLine(128*isamp, CM_sorting, 128*(isamp+1), CM_sorting);
	  l_Sorting[isamp]->SetLineColor(kRed);

	  //Calculate the Danning CM and plot it
	  double CM_danning = Danning_CM_offline(thisAPV, isamp);
	  l_Danning[isamp] = new TLine(128*isamp, CM_danning, 128*(isamp+1), CM_danning);
	  l_Danning[isamp]->SetLineColor(kGreen);

	  
	  //Calculate the Danning CM with more nsigmas and plot it
	  double CM_danning_nsigma = Danning_CM_offline(thisAPV, isamp, 30);
	  l_Danning_nsigma[isamp] = new TLine(128*isamp, CM_danning_nsigma, 128*(isamp+1), CM_danning_nsigma);
	  l_Danning_nsigma[isamp]->SetLineColor(6);
	  	  

	  //Calculate the Histo CM and plot it
	  double CM_histogramming = Histogramming_CM(thisAPV, isamp, 1);
	  l_Histogramming[isamp] = new TLine(128*isamp, CM_histogramming, 128*(isamp+1), CM_histogramming);
	  l_Histogramming[isamp]->SetLineColor(kBlue);

	  //Calculate the Histo CM online and plot it
	  double CM_histogramming_online = Histogramming_CM(thisAPV, isamp, 0);
	  l_Histogramming_online[isamp] = new TLine(128*isamp, CM_histogramming_online, 128*(isamp+1), CM_histogramming_online);
	  l_Histogramming_online[isamp]->SetLineColor(7);

	  //Calculate the Histo CM online with more sigmas and plot it
	  double CM_histogramming_online_nsigma = Histogramming_CM(thisAPV, isamp, 0, 30, 10, 100);
	  l_Histogramming_online_nsigma[isamp] = new TLine(128*isamp, CM_histogramming_online_nsigma, 128*(isamp+1), CM_histogramming_online_nsigma);
	  l_Histogramming_online_nsigma[isamp]->SetLineColor(8);
      
	}
	
	TLegend *leg = new TLegend(0.65,0.7,0.9,0.9);
      
	c->cd(iAPV + 1);
	thisAPV.hAPV->Draw();

	for(int isamp = 0; isamp < nsamples; isamp++){
	  //l_Sorting[isamp]->Draw("same");
	  //l_Danning[isamp]->Draw("same");
	  //l_Danning_nsigma[isamp]->Draw("same");
	  l_Histogramming[isamp]->Draw("same");
	  l_Histogramming_online[isamp]->Draw("same");
	  l_Histogramming_online_nsigma[isamp]->Draw("same");
	}
      
	//Only draw the legend for the first plot in the pdf
	//leg->AddEntry(l_Sorting[0],"Sorting CM","l");
	//leg->AddEntry(l_Danning[0],"Online Danning CM","l");
	//leg->AddEntry(l_Danning_nsigma[0],"Danning Online CM","l");
	leg->AddEntry(l_Histogramming_online[0],"Histo CM","l");
	leg->AddEntry(l_Histogramming_online_nsigma[0],"Histo CM 30#sigma","l");
	leg->AddEntry(l_Histogramming[0],"Correct CM","l");
	if(ipage == 0) leg->Draw("same");
      
      }
   
      c->Update();

      //Either print the pdf or look at one plot at a time
      if(!printpdf){
	reply = "no";
	cout<<"Showing event "<<nevent<<endl;
	cout << "press any key to continue (q to quit):" << endl;
	reply.ReadLine(cin,kFALSE);
	if( reply.BeginsWith("q") ) break;
      } else {
	if(ipage == 0) c->Print(pdffile + "(");
	else if(ipage == npages - 1){
	  c->Print(pdffile + ")");
	  break;
	}
	else c->Print(pdffile);
      
	ipage++;
      }
    }
    
  }

}
