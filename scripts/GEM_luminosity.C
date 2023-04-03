
const int nruns = 4;
const int nlayers = 5;
const int nmodules = 8;

double currents[nruns] = {5,15,30,45};
int runs_nocorr[nruns] = {2813,2815,2817,2820};
int runs_withcorr[nruns] = {2814,2816,2818,2821};
int nstrips[nlayers] = {3840,3840,3840,3840,5120};


double GEM_layer_efficiency(TFile *file, int ilayer){

  TH1D *hshould = (TH1D*)file->Get(Form("hshouldhit_xy_bb_gem_layer%i",ilayer));
  TH1D *hdid = (TH1D*)file->Get(Form("hdidhit_xy_bb_gem_layer%i",ilayer));

  return hdid->GetEntries() / hshould->GetEntries();

}


double landau_fit(TH1D *ADC){

  TF1 *landau = new TF1("land","landau",0,500);

  ADC->Fit("land","qR");
  
  return landau->GetParameter(1);
}

void make_ADC_plots(int imod, TCanvas *c,TH1D *h_adc[nmodules][nruns]){

  TLegend *legend2 = new TLegend(0.40,0.60,0.80,0.83);
  legend2->SetLineColor(0);
  
  int icolor = 0;
 
  for(int irun = 0; irun < nruns; irun++){

    icolor++;
    if(icolor == 5) icolor++; //no yellow
    
    h_adc[imod][irun]->Scale(1/h_adc[imod][irun]->GetEntries());
    h_adc[imod][irun]->GetYaxis()->SetRangeUser(0.0001,0.1);

    c->SetLogy();
    h_adc[imod][irun]->SetLineColor(icolor);
    if(irun == 0) h_adc[imod][irun]->Draw("hist");
    else h_adc[imod][irun]->Draw("same hist");

    double MPV = landau_fit(h_adc[imod][irun]);
 
    legend2->AddEntry(h_adc[imod][irun],Form("%g #muA, MPV = %g",currents[irun],MPV),"l");

  }
  
  legend2->Draw("same");


}

void GEM_luminosity(){

  gStyle->SetOptStat(0);

  TString Rootfiles = "/volatile/halla/sbs/jeffas/GEN_root/Rootfiles/GEM_luminosity/";

  TH2D *h2_CM[nruns];
  TH1D *h_adc[nmodules][nruns];
  double GEM_eff[nlayers+1][nruns];
  double GEM_occu[nlayers][nruns];
  double GEM_occu_err[nlayers][nruns];

  for(int irun = 0; irun < nruns; irun++){
    
    //cout<<runs_withcorr[irun]<<endl;
    //cout<<runs_nocorr[irun]<<endl;
    //TChain *TScal = new TChain("TSsbs");
    //TChain *T = new TChain("T");
    
    TFile *file = new TFile(Rootfiles + Form("gen_replayed_%i_all.root",runs_nocorr[irun]),"read");
    
    TFile *file_corr = new TFile(Rootfiles + Form("gen_replayed_%i_all.root",runs_withcorr[irun]),"read");
   
    h2_CM[irun] = (TH2D*)file->Get("hcommonmodeU_diff_bb_gem_m0");
    TH2D *h2_occu = (TH2D*)file->Get("hbb_gem_NstripsU_layer");
    
    for(int imod = 0; imod < nmodules; imod++){
      h_adc[imod][irun] = (TH1D*)file->Get(Form("hbb_gem_m%i_ADCmaxU_good",imod));
      //h_adc[imod][irun]->SetName(Form("run_%i_mod_%i",runs_nocorr[irun],imod));
    }
    

    for(int ilayer = 0; ilayer < nlayers; ilayer++){
      GEM_eff[ilayer][irun] = GEM_layer_efficiency(file,ilayer);
      GEM_occu[ilayer][irun] = h2_occu->ProjectionY("",ilayer+1,ilayer+1)->GetMean()/nstrips[ilayer];
      GEM_occu_err[ilayer][irun] = h2_occu->ProjectionY("",ilayer+1,ilayer+1)->GetRMS()/nstrips[ilayer];
    }
    
    GEM_eff[5][irun] = GEM_layer_efficiency(file_corr,0);


    /*
    TTree *TScal = (TTree*)file->Get("TSsbs");
    TTree *T = (TTree*)file->Get("T");
    //TScal->Add(Rootfiles + Form("*%i*",runs_nocorr[irun]));
    //T->Add(Rootfiles + Form("*%i*",runs_nocorr[irun]));

    //TScal->Add(Rootfiles + Form("*%i*",runs_withcorr[irun]));
    //T->Add(Rootfiles + Form("*%i*",runs_withcorr[irun]));

    T->SetBranchStatus("*",0);    
    TScal->SetBranchStatus("*",0);    

    double evnum; setrootvar::setbranch(T,"g","evnum",&evnum);
    double evnum_scaler; setrootvar::setbranch(TScal,"evNumber","",&evnum_scaler);
    double current_scaler; setrootvar::setbranch(TScal,"sbs.bcm.unew","current",&current_scaler);
   
    int nevent = 0;
    int currenttreenum = 0;
    int ntotal = 0;
  
    int nevent_scaler = 0;
    int event_switch = 0;
    double beam_current = 0;
    bool beam_ramp = true;
    
    while(T->GetEntry(nevent++)){
      
      //if(nevent%100000 == 0) cout<<nevent<<endl;
    
      
      if(nevent_scaler == 0){
	TScal->GetEntry(nevent_scaler++);
	beam_current = current_scaler;
	TScal->GetEntry(nevent_scaler++);
	event_switch = evnum_scaler;
      }

      if(evnum == event_switch){
	beam_current = current_scaler;
	TScal->GetEntry(nevent_scaler++);
	event_switch = evnum_scaler;
      }

      if(nevent%10000==0) cout<<evnum<<" "<<beam_current<<endl;

      
      if(beam_current > currents[irun] - 1){
	if(beam_ramp) cout<<"beam good "<<beam_current<<" "<<evnum<<endl;
	beam_ramp = false;
      }
      else if(beam_current < currents[irun] - 1){
	if(!beam_ramp) cout<<"beam bad "<<beam_current<<" "<<evnum<<endl;
	beam_ramp = true;
      }
      

      if(beam_current > currents[irun] - 1) ntotal++;
  
    }
    
    cout<<runs_withcorr[irun]<<" "<<ntotal<<" "<<T->GetEntries()<<endl;
    */
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
