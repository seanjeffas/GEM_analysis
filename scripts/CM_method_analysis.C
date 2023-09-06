////////////////////////////////////////////////////////////////
//  Written by Sean Jeffas
//  sj9ry@virginia.edu
//  Last updated March 30, 2023
//
//  The purpose of this script is to take full readout GEM data
//  generated from SBS-offline and perform common mode analysis.
//  It is set up to loop over some runs with different currents 
//  And plot the CM performance.
//
//  Usage: root CM_method_analysis.C
//
////////////////////////////////////////////////////////////////
#include "../src/pedestal_analysis.C"
#include "../src/GEM_map_decode.C"
#include "../src/Get_event_data.C"

const int nruns = 4;
int runs[nruns] = {2813,2815,2817,2820};
int currents[nruns] = {5,15,30,45};

//int runs[nruns] = {2817,2820};
//int currents[nruns] = {30,45};


void process_run(int irun,vector<vector<TH2F*>> hCM_Danning,vector<vector<TH2F*>> hCM_Sorting,vector<vector<TH2F*>> hCM_Danning_Online,vector<vector<TH2F*>> hCM_Histo_Online){

  gErrorIgnoreLevel = kError;

  cout<<"Process Run "<<runs[irun]<<endl;
  TString rootfile = Form("/volatile/halla/sbs/jeffas/GEN_root/Rootfiles/GEM_luminosity/full_readout/e1209016_replayed_%i_stream0_2_*_nevent5000*.root",runs[irun]);

  tree_init(rootfile);

  int nevent = 0;
  
  //Start looping over all events
  while(C->GetEntry(nevent++) && nevent < 5000){

    if(nevent % 100 == 0) cout<<nevent<<endl;
    
    Get_event_data(); //Function puts all the APV data into histograms
    
    for(int imod=0; imod < nmodules; imod++){
      for(int iaxis=0; iaxis < 2; iaxis++){
	for(int iAPV=0; iAPV < nAPV[imod][iaxis]; iAPV++){
	  
	  APV_info thisAPV = APV_data[imod][iaxis][iAPV];

	  //Loop over all time samples
	  for(int isamp = 0; isamp < nsamples; isamp++){
	    
	    //Calculate the sorting CM
	    double CM_sorting = Sorting_CM(APV_data[imod][iaxis][iAPV].hAPV,isamp);
	    
	    //Calculate the Danning CM
	    double CM_danning = Danning_CM_offline(APV_data[imod][iaxis][iAPV], isamp);

	    //Calculate the Danning online CM
	    double CM_danning_online = Danning_CM_offline(APV_data[imod][iaxis][iAPV], isamp, 30);
	    
	    //Calculate the Histo CM 
	    double CM_histogramming = Histogramming_CM(APV_data[imod][iaxis][iAPV], isamp);

	    //Calculate the Histo Online CM 
	    double CM_histogramming_online = Histogramming_CM(APV_data[imod][iaxis][iAPV], isamp, 0, 30,10,100);
	    
	    UpdateRollingCommonModeAverage(thisAPV,CM_histogramming);


	    if(nevent > 100){ //Skip first 100 event for rolling average
	      hCM_Danning[imod][iaxis]->Fill(iAPV,CM_danning - CM_histogramming);
	      hCM_Danning_Online[imod][iaxis]->Fill(iAPV,CM_danning_online - CM_histogramming);
	      hCM_Histo_Online[imod][iaxis]->Fill(iAPV,CM_histogramming_online - CM_histogramming);
	      hCM_Sorting[imod][iaxis]->Fill(iAPV,CM_sorting - CM_histogramming);
	      
	      //if(imod == 0 && iaxis == 0 && iAPV == 14)
	      //	cout<<nevent<<" "<<isamp<<" "<<CM_histogramming_online - CM_histogramming<<endl;
	    }
	  }
	}
      }
    }
  }
}

void make_graph(vector<vector<vector<TH2F*>>> hCM,vector<vector<vector<TGraphErrors*>>> &gCM, TString title){
  
  for(int irun = 0; irun < nruns; irun++){
    for(int imod=0; imod < nmodules; imod++){
      for(int iaxis=0; iaxis < 2; iaxis++){

	double x_APV[nAPV[imod][iaxis]], y_CM[nAPV[imod][iaxis]], y_CM_err[nAPV[imod][iaxis]];

	for(int iapv=0; iapv < nAPV[imod][iaxis]; iapv++)
	  x_APV[iapv] = iapv + 1;
	
	for(int iapv=0; iapv < nAPV[imod][iaxis]; iapv++){
	  
	  TF1 *fit = new TF1("fit","gaus");

	  hCM[irun][imod][iaxis]->ProjectionY("",x_APV[iapv],x_APV[iapv])->Fit("fit","q0");
	  y_CM[iapv] = fit->GetParameter(1);
	  y_CM_err[iapv] = fit->GetParameter(2);
	}
	
	gCM[irun][imod][iaxis] = new TGraphErrors(nAPV[imod][iaxis],x_APV,y_CM,0,y_CM_err);

	gCM[irun][imod][iaxis]->SetMarkerStyle(8);

	TString axis_str = "U/X";
	if(iaxis == 1) axis_str = "V/Y";
	
	gCM[irun][imod][iaxis]->SetTitle(Form(title + " CM Module %i %s Axis;APV;" + title + " - Histo Method CM (ADC)",imod,axis_str.Data()));
      }
    }
  }
  
}


void CM_method_analysis(){

  gROOT->SetBatch(kTRUE);
  
  int ped_run_number = 2360;

  // Varius root files and database files
  TString pedfile = Form("/w/halla-scshelf2102/sbs/jeffas/SBS_OFFLINE/SBS-replay/DB/gemped/daq_ped_bb_gem_run%i.dat",ped_run_number);
  TString cmfile = Form("/w/halla-scshelf2102/sbs/jeffas/SBS_OFFLINE/SBS-replay/DB/gemped/db_cmr_bb_gem_run%i.dat",ped_run_number);

  //Get our GEM map and APV map and load pedestals
  GEM_map_decode("../src/gem_map_BigBite.txt");
  InitAPVMAP();
  LoadPedestals(pedfile);
  LoadCM(cmfile);


  vector<vector<vector<TH2F*>>> hCM_Danning;
  vector<vector<vector<TH2F*>>> hCM_Danning_Online;
  vector<vector<vector<TH2F*>>> hCM_Histo_Online;
  vector<vector<vector<TH2F*>>> hCM_Sorting;

  hCM_Danning.resize(nruns,vector<vector<TH2F*> >(nmodules,vector<TH2F*>(2)));
  hCM_Danning_Online.resize(nruns,vector<vector<TH2F*> >(nmodules,vector<TH2F*>(2)));
  hCM_Histo_Online.resize(nruns,vector<vector<TH2F*> >(nmodules,vector<TH2F*>(2)));
  hCM_Sorting.resize(nruns,vector<vector<TH2F*> >(nmodules,vector<TH2F*>(2)));

  
  for( int irun = 0; irun < nruns; irun++){
    for(int imod=0; imod < nmodules; imod++){
      for(int iaxis=0; iaxis < 2; iaxis++){
	TString axis_str = "U/X";
	if(iaxis == 1) axis_str = "V/Y";

	hCM_Danning[irun][imod][iaxis] = new TH2F(Form("hCM_Danning_%i_%i_%i",irun,imod,iaxis),Form("Danning CM Module %i %s Axis;APV;Danning - Histo Method CM (ADC)",imod,axis_str.Data()),nAPV[imod][iaxis],0,nAPV[imod][iaxis],100,-100,100);
	hCM_Sorting[irun][imod][iaxis] = new TH2F(Form("hCM_Sorting_%i_%i_%i",irun,imod,iaxis),Form("Sorting CM Module %i %s Axis;APV;Sorting - Histo Method CM (ADC)",imod,axis_str.Data()),nAPV[imod][iaxis],0,nAPV[imod][iaxis],100,-100,100);
	hCM_Danning_Online[irun][imod][iaxis] = new TH2F(Form("hCM_Danning_Online_%i_%i_%i",irun,imod,iaxis),Form("Danning Online CM Module %i %s Axis;APV;Danning Online - Histo Method CM (ADC)",imod,axis_str.Data()),nAPV[imod][iaxis],0,nAPV[imod][iaxis],100,-100,100);
	hCM_Histo_Online[irun][imod][iaxis] = new TH2F(Form("hCM_Histo_Online_%i_%i_%i",irun,imod,iaxis),Form("Histo Online CM Module %i %s Axis;APV;Histo Online - Histo Method CM (ADC)",imod,axis_str.Data()),nAPV[imod][iaxis],0,nAPV[imod][iaxis],100,-100,100);

	
      }
    }
    
    process_run(irun,hCM_Danning[irun],hCM_Sorting[irun],hCM_Danning_Online[irun],hCM_Histo_Online[irun]);
  }


  vector<vector<vector<TGraphErrors*>>> gCM_Danning;
  vector<vector<vector<TGraphErrors*>>> gCM_Danning_Online;
  vector<vector<vector<TGraphErrors*>>> gCM_Histo_Online;
  vector<vector<vector<TGraphErrors*>>> gCM_Sorting;
  
  gCM_Danning.resize(nruns,vector<vector<TGraphErrors*> >(nmodules,vector<TGraphErrors*>(2)));
  gCM_Danning_Online.resize(nruns,vector<vector<TGraphErrors*> >(nmodules,vector<TGraphErrors*>(2)));
  gCM_Histo_Online.resize(nruns,vector<vector<TGraphErrors*> >(nmodules,vector<TGraphErrors*>(2)));
  gCM_Sorting.resize(nruns,vector<vector<TGraphErrors*> >(nmodules,vector<TGraphErrors*>(2)));


  for(int irun = 0; irun < nruns; irun++){

    TString pdfname = Form("../plots/CM_comparison_run%i.pdf",runs[irun]);      
    
    for(int imod=0; imod < nmodules; imod++){
      TCanvas *c = new TCanvas("c","",1200,1000);
      c->Divide(2,3);
      c->cd(1);
      hCM_Danning[irun][imod][0]->Draw("colz");
      c->cd(2);
      hCM_Danning[irun][imod][1]->Draw("colz");
      c->cd(3);
      hCM_Danning_Online[irun][imod][0]->Draw("colz");
      c->cd(4);
      hCM_Danning_Online[irun][imod][1]->Draw("colz");
      c->cd(5);
      hCM_Histo_Online[irun][imod][0]->Draw("colz");
      c->cd(6);
      hCM_Histo_Online[irun][imod][1]->Draw("colz");
    
      if(imod == 0) c->Print(pdfname + "(");
      else if(imod == nmodules - 1) c->Print(pdfname + ")");
      else c->Print(pdfname);
    }
  }
   
  make_graph(hCM_Danning,gCM_Danning,"Danning");
  make_graph(hCM_Sorting,gCM_Sorting,"Sorting");
  make_graph(hCM_Danning_Online,gCM_Danning_Online,"Danning Online");
  make_graph(hCM_Histo_Online,gCM_Histo_Online,"Histo Online");

    
  TString pdfname_all = "../plots/CM_comparison_study.pdf";      

  for(int imod=0; imod < nmodules; imod++){
    
    TCanvas *c = new TCanvas("c","",1200,1000);
    c->Divide(2,3);
	
    TLegend *legend = new TLegend(0.60,0.1,0.9,0.40);  
    int icolor = 0;
    
    for(int irun = 0; irun < nruns; irun++){
      
      icolor++;
      if(icolor == 5 || icolor == 10) icolor ++;

      for(int iaxis=0; iaxis < 2; iaxis++){
	gCM_Danning[irun][imod][iaxis]->SetLineColor(icolor);
	gCM_Danning_Online[irun][imod][iaxis]->SetLineColor(icolor);
	gCM_Histo_Online[irun][imod][iaxis]->SetLineColor(icolor);
	gCM_Sorting[irun][imod][iaxis]->SetLineColor(icolor);
	gCM_Danning[irun][imod][iaxis]->SetMarkerColor(icolor);
	gCM_Sorting[irun][imod][iaxis]->SetMarkerColor(icolor);
	gCM_Danning_Online[irun][imod][iaxis]->SetMarkerColor(icolor);
	gCM_Histo_Online[irun][imod][iaxis]->SetMarkerColor(icolor);
      }      


      legend->AddEntry(gCM_Danning[irun][imod][0],Form("%i uA",currents[irun]),"l");
		       
      c->cd(1);
      if(irun == 0) gCM_Danning[irun][imod][0]->Draw("AP");
      else gCM_Danning[irun][imod][0]->Draw("P");
      c->cd(2);
      if(irun == 0) gCM_Danning[irun][imod][1]->Draw("AP");
      else gCM_Danning[irun][imod][1]->Draw("P");
      c->cd(3);
      if(irun == 0) gCM_Danning_Online[irun][imod][0]->Draw("AP");
      else gCM_Sorting[irun][imod][0]->Draw("P");
      c->cd(4);
      if(irun == 0) gCM_Danning_Online[irun][imod][1]->Draw("AP");
      else gCM_Sorting[irun][imod][1]->Draw("P");
      c->cd(5);
      if(irun == 0) gCM_Histo_Online[irun][imod][0]->Draw("AP");
      else gCM_Histo_Online[irun][imod][0]->Draw("P");
      c->cd(6);
      if(irun == 0) gCM_Histo_Online[irun][imod][1]->Draw("AP");
      else gCM_Histo_Online[irun][imod][1]->Draw("P");

      gCM_Danning[irun][imod][0]->GetYaxis()->SetRangeUser(-100,100);
      gCM_Danning[irun][imod][1]->GetYaxis()->SetRangeUser(-100,100);
      gCM_Danning_Online[irun][imod][0]->GetYaxis()->SetRangeUser(-100,100);
      gCM_Danning_Online[irun][imod][1]->GetYaxis()->SetRangeUser(-100,100);
      gCM_Histo_Online[irun][imod][0]->GetYaxis()->SetRangeUser(-100,100);
      gCM_Histo_Online[irun][imod][1]->GetYaxis()->SetRangeUser(-100,100);
    }

    c->cd(1);
    legend->Draw("same");
    
    if(imod == 0) c->Print(pdfname_all + "(");
    else if(imod == nmodules - 1) c->Print(pdfname_all + ")");
    else c->Print(pdfname_all);

  }
     
       
}
