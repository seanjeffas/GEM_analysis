////////////////////////////////////////////////////////////////
//  Written by Sean Jeffas
//  sj9ry@virginia.edu
//  Last updated March 23, 2023
//
//  The purpose of this script is to take specific runs from the 
//  GEM luminosity study and print plots showing how the GEM 
//  efficiency and gain performed
//
//  Usage: root GEM_luminosity
//
////////////////////////////////////////////////////////////////


const int nruns = 4;
const int nlayers = 5;
const int nmodules = 8;

double currents[nruns] = {5,15,30,45};
int runs_nocorr[nruns] = {2813,2815,2817,2820};  //Layer 0 has no HV corrections
int runs_withcorr[nruns] = {2814,2816,2818,2821};//Layer 0 has HV corrections 
int nstrips[nlayers] = {3840,3840,3840,3840,5120};


// This gets the efficiency from the root tree
double GEM_layer_efficiency(TFile *file, int ilayer){

  TH1D *hshould = (TH1D*)file->Get(Form("hshouldhit_xy_bb_gem_layer%i",ilayer));
  TH1D *hdid = (TH1D*)file->Get(Form("hdidhit_xy_bb_gem_layer%i",ilayer));

  return hdid->GetEntries() / hshould->GetEntries();

}

// This gets the MPV from a landau fit
double landau_fit(TH1D *ADC){

  TF1 *landau = new TF1("land","landau",0,500);

  ADC->Fit("land","qR");
  
  return landau->GetParameter(1);
}

//This function plots the ADC spectra for all modules
void make_ADC_plots(int imod, TCanvas *c,TH1D *h_adc[nmodules][nruns]){

  TLegend *legend2 = new TLegend(0.40,0.60,0.80,0.83);
  legend2->SetLineColor(0);
  
  int icolor = 0;
 
  //loop over all runs
  for(int irun = 0; irun < nruns; irun++){

    icolor++;
    if(icolor == 5) icolor++; //no yellow
    
    //Normalize the scale
    h_adc[imod][irun]->Scale(1/h_adc[imod][irun]->GetEntries());
    h_adc[imod][irun]->GetYaxis()->SetRangeUser(0.0001,0.1);

    c->SetLogy();
    h_adc[imod][irun]->SetLineColor(icolor);
    if(irun == 0) h_adc[imod][irun]->Draw("hist");
    else h_adc[imod][irun]->Draw("same hist");

    //Write the MPV on the plot
    double MPV = landau_fit(h_adc[imod][irun]); 
    legend2->AddEntry(h_adc[imod][irun],Form("%g #muA, MPV = %g",currents[irun],MPV),"l");

  }
  
  legend2->Draw("same");


}

//This is the main function
void GEM_luminosity(){

  gStyle->SetOptStat(0);

  TString Rootfiles = "/volatile/halla/sbs/jeffas/GEN_root/Rootfiles/GEM_luminosity/";

  TH2D *h2_CM[nruns];
  TH1D *h_adc[nmodules][nruns];
  double GEM_eff[nlayers+1][nruns];  //Add the HV corrected runs as an "extra layer"
  double GEM_occu[nlayers][nruns];
  double GEM_occu_err[nlayers][nruns];

  //Loop over all runs
  for(int irun = 0; irun < nruns; irun++){
        
    //Read in the root files
    TFile *file = new TFile(Rootfiles + Form("gen_replayed_%i_all.root",runs_nocorr[irun]),"read");
    TFile *file_corr = new TFile(Rootfiles + Form("gen_replayed_%i_all.root",runs_withcorr[irun]),"read");
   
    h2_CM[irun] = (TH2D*)file->Get("hcommonmodeU_diff_bb_gem_m0");
    TH2D *h2_occu = (TH2D*)file->Get("hbb_gem_NstripsU_layer");
    
    //Loop over all modules 
    for(int imod = 0; imod < nmodules; imod++){
      h_adc[imod][irun] = (TH1D*)file->Get(Form("hbb_gem_m%i_ADCmaxU_good",imod));
    }
    
    //Loop over all layers ad get efficiency and occupancy numbers
    for(int ilayer = 0; ilayer < nlayers; ilayer++){
      GEM_eff[ilayer][irun] = GEM_layer_efficiency(file,ilayer);
      GEM_occu[ilayer][irun] = h2_occu->ProjectionY("",ilayer+1,ilayer+1)->GetMean()/nstrips[ilayer];
      GEM_occu_err[ilayer][irun] = h2_occu->ProjectionY("",ilayer+1,ilayer+1)->GetRMS()/nstrips[ilayer];
    }
    
    //Add the HV corrections as another layer to the plot
    GEM_eff[5][irun] = GEM_layer_efficiency(file_corr,0);


  }
  
  
  
  TCanvas *c1 = new TCanvas("c1","",800,600);
  TCanvas *c2 = new TCanvas("c2","",800,600);
  TGraph *geff[nlayers+1];
  TGraphErrors *goccu[nlayers];
  TLegend *legend = new TLegend(0.12,0.12,0.4,0.4);
  legend->SetLineColor(0);

  int icolor = 0;

  for(int ilayer=0; ilayer < nlayers; ilayer++){
    
    icolor++;
    if(icolor == 5) icolor++; //no yellow

    c1->cd();
    geff[ilayer] = new TGraph(nruns,currents,GEM_eff[ilayer]);
    geff[ilayer]->SetMarkerStyle(8);
    geff[ilayer]->SetMarkerColor(icolor);
    if(ilayer == 0) geff[ilayer]->Draw("AP");
    else geff[ilayer]->Draw("P same");
    
    legend->AddEntry(geff[ilayer],Form("Layer %i",ilayer),"p");

    c2->cd();
    goccu[ilayer] = new TGraphErrors(nruns,currents,GEM_occu[ilayer],0,GEM_occu_err[ilayer]);
    goccu[ilayer]->SetMarkerStyle(8);
    goccu[ilayer]->SetMarkerColor(icolor);
    goccu[ilayer]->SetLineColor(icolor);
    if(ilayer == 0) goccu[ilayer]->Draw("AP");
    else goccu[ilayer]->Draw("P same");

  }

 
  c1->cd();
  geff[5] = new TGraph(nruns,currents,GEM_eff[5]);
  geff[5]->SetMarkerStyle(8);
  geff[5]->SetMarkerColor(++icolor);
  geff[5]->Draw("P same");
  legend->AddEntry(geff[5],"Layer 0 corrected","p");
  legend->Draw("same");
  geff[0]->SetTitle("BB GEM Luminosity Study; Beam Current (#muA); GEM Efficiency");
  geff[0]->GetYaxis()->SetRangeUser(0.5,1);

  goccu[0]->SetTitle("BB GEM Occupancies; Beam Current (#muA); GEM Occupancy (U axis)");
  goccu[0]->GetYaxis()->SetRangeUser(0,0.3);



  TString outfile = "../plots/GEM_luminosity_study.pdf";

  c1->Print(outfile + "(");
  c2->Print(outfile);

  TCanvas *c3[nmodules];
  
  for(int imod = 0; imod < nmodules; imod++){
    c3[imod] = new TCanvas(Form("c3_%i",imod),"",800,600);
    make_ADC_plots(imod,c3[imod],h_adc);
  
    if(imod == nmodules -1) c3[imod]->Print(outfile + ")");
    else c3[imod]->Print(outfile);
  }


  
  
  
}
