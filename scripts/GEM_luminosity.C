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

double modx[nmodules] = {0, -0.001931935, -0.007540193, -0.009137839, -0.7647015, -0.2523657, 0.2586995, 0.7706357};
double mody[nmodules] = {0, -0.01142927, 0.005879457, 0.004917223, 0.0020337, -0.002284228, -0.008135348, -0.01182199};
double modz[nmodules] = {0, 0.165205, 0.2790333, 0.2790333, 1.578944, 1.535895, 1.578437, 1.537177};
double modax[nmodules] = {0, 0.07756231, 0.1309815, 0.2160082, -0.2102473, -0.2292616, -0.1861139, -0.204406};
double moday[nmodules] = {0, 0.04990183, 0.07427616, 0.07427616, -0.2975331, -0.04278965, -0.183984, -0.0703861};
double modaz[nmodules] = {0, 0.06161188, -0.2973769, -0.4551869, -0.3461524, -0.2910578, -0.250579, -0.2781818};

double mod_sizex[nmodules] = {1.5, 1.5, 1.5, 1.5, 0.512, 0.512, 0.512, 0.512};
double mod_sizey[nmodules] = {0.4, 0.4, 0.4, 0.4, 0.6144, 0.6144, 0.6144, 0.6144};

double optics_pos[7] = {-0.223403, -0.148302, -0.0737851, 0.000634964, 0.0754245, 0.149334, 0.224007};
double optics_sigma[7] = {0.00813423, 0.00693141, 0.00662755, 0.00675542, 0.00718186, 0.00659294, 0.00654504};

TVector3 fXax, fYax, fZax;



TVector3 TrackIntersect( int module, TVector3 track_origin, TVector3 track_direction){
  TVector3 modpos(modx[module],mody[module],modz[module]);
  TVector3 modzaxis = fZax;
  
  double sintersect = modzaxis.Dot( modpos - track_origin )/modzaxis.Dot( track_direction );
  
  return track_origin + sintersect * track_direction;
}

bool IsInActiveArea(int module, TVector3 point){

  TVector3 fOrigin(modx[module], mody[module], modz[module]);
  TVector3 v = point - fOrigin;
  TVector3 det(v.Dot(fXax), v.Dot(fYax), v.Dot(fZax));
  

  return (abs(det.X()) <= 0.5*mod_sizex[module] && abs(det.Y()) <= 0.5*mod_sizey[module]);

}

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


void CalculateEfficiency(int irun, TFile *file, double GEM_eff[][nruns]){

  
  TTree *C = (TTree*)file->Get("T");


  C->SetBranchStatus("*",0);

  int maxtr = 1000;
  int maxhit = 10000;
  double tr_n, hit_n;
  double tr_x[maxtr], tr_y[maxtr], tr_xp[maxtr], tr_yp[maxtr], tr_vz[maxtr], px[maxtr], py[maxtr], pz[maxtr], p[maxtr];
  double hit_tr_i[maxhit], hit_mod[maxhit];
  double ps_e;

  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.x",1);
  C->SetBranchStatus("bb.tr.y",1);
  C->SetBranchStatus("bb.tr.th",1);
  C->SetBranchStatus("bb.tr.ph",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.p",1);
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.gem.hit.ngoodhits",1);
  C->SetBranchStatus("bb.gem.hit.trackindex",1);
  C->SetBranchStatus("bb.gem.hit.module",1);
  C->SetBranchStatus("bb.ps.e",1);

  C->SetBranchAddress("bb.tr.n",&tr_n);
  C->SetBranchAddress("bb.tr.x",tr_x);
  C->SetBranchAddress("bb.tr.y",tr_y);
  C->SetBranchAddress("bb.tr.th",tr_xp);
  C->SetBranchAddress("bb.tr.ph",tr_yp);
  C->SetBranchAddress("bb.tr.px",px);
  C->SetBranchAddress("bb.tr.py",py);
  C->SetBranchAddress("bb.tr.pz",pz);
  C->SetBranchAddress("bb.tr.p",p);
  C->SetBranchAddress("bb.tr.vz",tr_vz);
  C->SetBranchAddress("bb.gem.hit.ngoodhits",&hit_n);
  C->SetBranchAddress("bb.gem.hit.trackindex",hit_tr_i);
  C->SetBranchAddress("bb.gem.hit.module",hit_mod);
  C->SetBranchAddress("bb.ps.e",&ps_e);

  int nevent = 0;
  int mod_didhit[nlayers] = {0};
  int mod_shouldhit[nlayers] = {0};

  while(C->GetEntry(nevent++)){

    if(nevent % 10000 == 0) cout<<irun<<" "<<nevent<<endl;
    
    bool didhit[nmodules] = {false};
    bool shouldhit[nmodules] = {false};

    if(tr_n == 1 && ps_e > 0.15){

      bool optics_cut = false;

      for(int ioptics = 0; ioptics < 7; ioptics++)
	if(tr_vz[0] > optics_pos[ioptics] - 2*optics_sigma[ioptics] && tr_vz[0] < optics_pos[ioptics] + 2*optics_sigma[ioptics])
	  optics_cut = true;

      if(!optics_cut) continue;

      TLorentzVector Peprime(px[0],py[0],pz[0],p[0]);
      
      TVector3 track_origin(tr_x[0],tr_y[0],0.0);
      TVector3 track_dir(tr_xp[0],tr_yp[0],1.0);

      
      //double Q2recon = 2.0*Pe.E()*Peprime.E()*(1.0-acos(Peprime.Pz() / Peprime.E()));
      //double W2recon = pow(0.938,2.0) + 2.0*0.938*(Pe.E()-Peprime.E()) - Q2recon;
      
      for(int ihit = 0; ihit < hit_n; ihit++)
   	didhit[(int)hit_mod[ihit]] = true;
      
      
      for(int imod = 0; imod < nmodules; imod++){
	TVector3 Intersect = TrackIntersect(imod,track_origin,track_dir);
	bool isinactivearea = IsInActiveArea(imod,Intersect);
	
	if(isinactivearea) shouldhit[imod] = true;

	int layer;
	if(imod >= 4) layer = 4;
	else layer = imod;

	if(didhit[imod]) mod_didhit[layer]++;
	if(shouldhit[imod]) mod_shouldhit[layer]++;
      } 
      
    }
  }

  for(int ilayer = 0; ilayer < nlayers; ilayer++){
    GEM_eff[ilayer][irun] = mod_didhit[ilayer]*1.0/mod_shouldhit[ilayer];
  }
}



//This is the main function
void GEM_luminosity(){

  gStyle->SetOptStat(0);

  TString Rootfiles = "/volatile/halla/sbs/jeffas/GEN_root/Rootfiles/GEM_luminosity/zero_suppressed/";

  TH2D *h2_CM[nruns];
  TH1D *h_adc[nmodules][nruns];
  double GEM_eff[nlayers+1][nruns];  //Add the HV corrected runs as an "extra layer"
  double GEM_occu[nlayers][nruns];
  double GEM_occu_err[nlayers][nruns];

  //Initialize modules coordinates (used for efficiency calculation)
  for(int imod=0; imod < nmodules; imod++){
    
    TRotation RotTemp;

    RotTemp.RotateX( modax[imod] * TMath::DegToRad() );
    RotTemp.RotateY( moday[imod] * TMath::DegToRad() );
    RotTemp.RotateZ( modaz[imod] * TMath::DegToRad() );

    fXax.SetXYZ( RotTemp.XX(), RotTemp.YX(), RotTemp.ZX() );
    fYax.SetXYZ( RotTemp.XY(), RotTemp.YY(), RotTemp.ZY() );
    fZax.SetXYZ( RotTemp.XZ(), RotTemp.YZ(), RotTemp.ZZ() );
    
  }


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

    CalculateEfficiency(irun, file, GEM_eff);
    
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
