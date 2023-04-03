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



void process_run(int irun,TH2F *hCM_Danning[][2],TH2F *hCM_Histo[][2]){

  TString rootfile = Form("/volatile/halla/sbs/jeffas/GEN_root/Rootfiles/e1209016_replayed_%i_*.root",runs[irun]);

  tree_init(rootfile);

  int nevent = 0;
  
  //Start looping over all events
  while(C->GetEntry(nevent++) && nevent < 1000){

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
	    
	    //Calculate the Histo CM 
	    double CM_histogramming = Histogramming_CM(APV_data[imod][iaxis][iAPV], isamp);
	    
	    UpdateRollingCommonModeAverage(thisAPV,CM_histogramming);

	    hCM_Danning[imod][iaxis]->Fill(iAPV,CM_danning - CM_sorting);
	    hCM_Histo[imod][iaxis]->Fill(iAPV,CM_histogramming - CM_sorting);
	  }
	}
      }
    }
  }


}

void CM_method_analysis(){

  gROOT->SetBatch(kTRUE);
  
  int ped_run_number = 2360;

  // Varius root files and database files
  TString pedfile = Form("/w/halla-scshelf2102/sbs/jeffas/SBS_OFFLINE/SBS-replay/DB/gemped/daq_ped_bb_gem_run%i.dat",ped_run_number);
  TString cmfile = Form("/w/halla-scshelf2102/sbs/jeffas/SBS_OFFLINE/SBS-replay/DB/gemped/daq_cmr_bb_gem_run%i.dat",ped_run_number);

  //Get our GEM map and APV map and load pedestals
  GEM_map_decode("gem_map_BigBite.txt");
  InitAPVMAP();
  LoadPedestals(pedfile);
  LoadCM(cmfile);



  TH2F *hCM_Danning[nruns][nmodules][2];
  TH2F *hCM_Histo[nruns][nmodules][2];

  for( int irun = 0; irun < nruns; irun++){
    for(int imod=0; imod < nmodules; imod++){
      for(int iaxis=0; iaxis < 2; iaxis++){
	TString axis_str = "U/X";
	if(iaxis == 1) axis_str = "V/Y";

	hCM_Danning[irun][imod][iaxis] = new TH2F(Form("hCM_Danning_%i_%i_%i",irun,imod,iaxis),Form("Danning CM Module %i %s Axis;APV;Histo - Sorting Method CM (ADC)",imod,axis_str.Data()),nAPV[imod][iaxis],0,nAPV[imod][iaxis],100,-100,100);
	hCM_Histo[irun][imod][iaxis] = new TH2F(Form("hCM_Histo_%i_%i_%i",irun,imod,iaxis),Form("Histo CM Module %i %s Axis;APV;Histo - Sorting Method CM (ADC)",imod,axis_str.Data()),nAPV[imod][iaxis],0,nAPV[imod][iaxis],100,-100,100);

	
      }
    }
  
    process_run(irun,hCM_Danning[irun],hCM_Histo[irun]);
  }


  TGraphErrors *gCM_Danning[nruns][nmodules][2];
  TGraphErrors *gCM_Histo[nruns][nmodules][2];
  

  for(int irun = 0; irun < nruns; irun++){

    TString pdfname = Form("../plots/CM_comparison_run%i.pdf",runs[irun]);      
    
    for(int imod=0; imod < nmodules; imod++){
      TCanvas *c = new TCanvas("c","",1200,800);
      c->Divide(2,2);
      c->cd(1);
      hCM_Danning[irun][imod][0]->Draw("colz");
      c->cd(2);
      hCM_Histo[irun][imod][0]->Draw("colz");
      c->cd(3);
      hCM_Danning[irun][imod][1]->Draw("colz");
      c->cd(4);
      hCM_Histo[irun][imod][1]->Draw("colz");
    
      if(imod == 0) c->Print(pdfname + "(");
      else if(imod == nmodules - 1) c->Print(pdfname + ")");
      else c->Print(pdfname);

      
      for(int iaxis=0; iaxis < 2; iaxis++){

	double x_APV[nAPV[imod][iaxis]], y_CM_Danning[nAPV[imod][iaxis]], y_CM_Danning_err[nAPV[imod][iaxis]];
	double y_CM_Histo[nAPV[imod][iaxis]], y_CM_Histo_err[nAPV[imod][iaxis]];

	for(int iapv=0; iapv < nAPV[imod][iaxis]; iapv++)
	  x_APV[iapv] = iapv + 1;
	
	for(int iapv=0; iapv < nAPV[imod][iaxis]; iapv++){
	  
	  TF1 *fit = new TF1("fit","gaus");

	  hCM_Danning[irun][imod][iaxis]->ProjectionY("",x_APV[iapv],x_APV[iapv])->Fit("fit","q0");
	  y_CM_Danning[iapv] = fit->GetParameter(1);
	  y_CM_Danning_err[iapv] = fit->GetParameter(2);

	  hCM_Histo[irun][imod][iaxis]->ProjectionY("",x_APV[iapv],x_APV[iapv])->Fit("fit","q0");
	  y_CM_Histo[iapv] = fit->GetParameter(1);
	  y_CM_Histo_err[iapv] = fit->GetParameter(2);
	  

	}
	
	gCM_Danning[irun][imod][iaxis] = new TGraphErrors(nAPV[imod][iaxis],x_APV,y_CM_Danning,0,y_CM_Danning_err);
	gCM_Histo[irun][imod][iaxis] = new TGraphErrors(nAPV[imod][iaxis],x_APV,y_CM_Histo,0,y_CM_Histo_err);

	TString axis_str = "U/X";
	if(iaxis == 1) axis_str = "V/Y";
	
	gCM_Danning[irun][imod][iaxis]->SetTitle(Form("Danning CM Module %i %s Axis;APV;Danning - Sorting Method CM (ADC)",imod,axis_str.Data()));
	gCM_Histo[irun][imod][iaxis]->SetTitle(Form("Histogramming CM Module %i %s Axis;APV;Histo - Sorting Method CM (ADC)",imod,axis_str.Data()));

	gCM_Danning[irun][imod][iaxis]->SetMarkerStyle(8);
	gCM_Histo[irun][imod][iaxis]->SetMarkerStyle(8);
	
      }
    }
  }
  
  TString pdfname_all = "../plots/CM_comparison_study.pdf";      

  for(int imod=0; imod < nmodules; imod++){
    
    TCanvas *c = new TCanvas("c","",1200,800);
    c->Divide(2,2);
	
    TLegend *legend = new TLegend(0.60,0.1,0.9,0.40);  
    int icolor = 0;
    
    for(int irun = 0; irun < nruns; irun++){
      
      icolor++;
      if(icolor == 5 || icolor == 10) icolor ++;

      for(int iaxis=0; iaxis < 2; iaxis++){
	gCM_Danning[irun][imod][iaxis]->SetLineColor(icolor);
	gCM_Histo[irun][imod][iaxis]->SetLineColor(icolor);
	gCM_Danning[irun][imod][iaxis]->SetMarkerColor(icolor);
	gCM_Histo[irun][imod][iaxis]->SetMarkerColor(icolor);
      }      


      legend->AddEntry(gCM_Danning[irun][imod][0],Form("%i uA",currents[irun]),"l");
		       
      c->cd(1);
      if(irun == 0) gCM_Danning[irun][imod][0]->Draw("AP");
      else gCM_Danning[irun][imod][0]->Draw("P");
      c->cd(2);
      if(irun == 0) gCM_Histo[irun][imod][0]->Draw("AP");
      else gCM_Histo[irun][imod][0]->Draw("P");
      c->cd(3);
      if(irun == 0) gCM_Danning[irun][imod][1]->Draw("AP");
      else gCM_Danning[irun][imod][1]->Draw("P");
      c->cd(4);
      if(irun == 0) gCM_Histo[irun][imod][1]->Draw("AP");
      else gCM_Histo[irun][imod][1]->Draw("P");

      gCM_Danning[irun][imod][0]->GetYaxis()->SetRangeUser(-100,100);
      gCM_Danning[irun][imod][1]->GetYaxis()->SetRangeUser(-100,100);
      gCM_Histo[irun][imod][0]->GetYaxis()->SetRangeUser(-100,100);
      gCM_Histo[irun][imod][1]->GetYaxis()->SetRangeUser(-100,100);
    }

    c->cd(1);
    legend->Draw("same");
    
    if(imod == 0) c->Print(pdfname_all + "(");
    else if(imod == nmodules - 1) c->Print(pdfname_all + ")");
    else c->Print(pdfname_all);

  }
     
       
}
